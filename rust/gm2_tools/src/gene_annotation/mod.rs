//! Protein-guided gene annotation for TipSeek's ordinary gene workflow.
//!
//! Coordinates are represented internally as zero-based, half-open intervals.
//! DNA readers in this module preserve IUPAC ambiguity characters so external
//! annotation coordinates always refer to the sequence that was actually read.

use std::collections::BTreeMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

mod model;
mod select;

pub use model::{
    build_gene_model, ExonSequence, GeneModel, IntronSequence, ModelConfig, ModelState, SpliceClass,
};
pub use select::select_competing_models;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DnaRecord {
    pub id: String,
    pub sequence: String,
}

pub fn read_dna_fasta(path: &Path) -> io::Result<Vec<DnaRecord>> {
    let reader = BufReader::new(File::open(path)?);
    let mut records = Vec::new();
    let mut id: Option<String> = None;
    let mut sequence = String::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = id.take() {
                records.push(DnaRecord {
                    id: previous,
                    sequence: normalize_dna(&sequence)?,
                });
            }
            id = Some(header.split_whitespace().next().unwrap_or("").to_owned());
            sequence.clear();
        } else {
            sequence.push_str(line.trim());
        }
    }
    if let Some(previous) = id {
        records.push(DnaRecord {
            id: previous,
            sequence: normalize_dna(&sequence)?,
        });
    }
    Ok(records)
}

fn normalize_dna(sequence: &str) -> io::Result<String> {
    let mut normalized = String::with_capacity(sequence.len());
    for base in sequence.bytes() {
        let base = base.to_ascii_uppercase();
        normalized.push(match base {
            b'A' | b'C' | b'G' | b'T' | b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B'
            | b'D' | b'H' | b'V' | b'-' => base as char,
            b'U' => 'T',
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid DNA character: {}", base as char),
                ));
            }
        });
    }
    Ok(normalized)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    pub fn parse(value: &str) -> Option<Self> {
        match value {
            "+" => Some(Self::Forward),
            "-" => Some(Self::Reverse),
            _ => None,
        }
    }

    pub fn as_char(self) -> char {
        match self {
            Self::Forward => '+',
            Self::Reverse => '-',
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct Interval {
    pub start: usize,
    pub end: usize,
}

impl Interval {
    pub fn new(start: usize, end: usize) -> Option<Self> {
        (start < end).then_some(Self { start, end })
    }

    pub fn len(self) -> usize {
        self.end - self.start
    }

    pub fn is_empty(self) -> bool {
        false
    }

    pub fn overlap(self, other: Self) -> usize {
        self.end
            .min(other.end)
            .saturating_sub(self.start.max(other.start))
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CdsFeature {
    pub interval: Interval,
    pub strand: Strand,
    pub phase: Option<u8>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PafRecord {
    pub query_id: String,
    pub query_len: usize,
    pub query_interval: Interval,
    pub strand: Strand,
    pub target_id: String,
    pub target_len: usize,
    pub target_interval: Interval,
    pub matches_nt: usize,
    pub block_len_nt: usize,
    pub mapq: u8,
    pub alignment_score: Option<i64>,
    pub frameshifts: usize,
    pub inframe_stops: usize,
    pub cigar: Option<String>,
    pub cs: Option<String>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RawGeneModel {
    pub gff_id: String,
    pub candidate_id: String,
    pub protein_id: String,
    pub target_interval_aa: Interval,
    pub strand: Strand,
    pub genomic_interval: Interval,
    pub gff_score: Option<i64>,
    pub identity_millionths: Option<u32>,
    pub positive_millionths: Option<u32>,
    pub frameshifts: usize,
    pub inframe_stops: usize,
    pub rank: Option<usize>,
    pub cdss: Vec<CdsFeature>,
    pub paf: Option<PafRecord>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseReport {
    pub models: Vec<RawGeneModel>,
    pub unassigned_paf: usize,
    pub warnings: Vec<String>,
}

#[derive(Clone, Debug)]
struct MrnaFeature {
    id: String,
    candidate_id: String,
    protein_id: String,
    target_interval_aa: Interval,
    strand: Strand,
    genomic_interval: Interval,
    score: Option<i64>,
    identity_millionths: Option<u32>,
    positive_millionths: Option<u32>,
    frameshifts: usize,
    inframe_stops: usize,
    rank: Option<usize>,
}

fn parse_attributes(value: &str) -> BTreeMap<&str, &str> {
    value
        .split(';')
        .filter_map(|field| field.split_once('='))
        .map(|(key, value)| (key.trim(), value.trim()))
        .collect()
}

fn decimal_millionths(value: Option<&str>) -> Option<u32> {
    let value = value?.parse::<f64>().ok()?;
    (value.is_finite() && value >= 0.0).then(|| (value.min(1.0) * 1_000_000.0).round() as u32)
}

fn parse_gff_interval(start: &str, end: &str) -> Option<Interval> {
    Interval::new(
        start.parse::<usize>().ok()?.checked_sub(1)?,
        end.parse().ok()?,
    )
}

fn parse_target(value: Option<&str>) -> Option<(String, Interval)> {
    let fields: Vec<_> = value?.split_whitespace().collect();
    Some((
        fields.first()?.to_string(),
        Interval::new(
            fields.get(1)?.parse::<usize>().ok()?.checked_sub(1)?,
            fields.get(2)?.parse().ok()?,
        )?,
    ))
}

fn paf_tag<'a>(tags: impl Iterator<Item = &'a str>, key: &str) -> Option<&'a str> {
    tags.filter_map(|tag| tag.split_once(':'))
        .find_map(|(tag_key, value)| {
            if tag_key == key {
                value.split_once(':').map(|(_, value)| value)
            } else {
                None
            }
        })
}

pub fn parse_paf_line(line: &str) -> Result<PafRecord, String> {
    let fields: Vec<_> = line.trim_start_matches("##PAF\t").split('\t').collect();
    if fields.len() < 12 {
        return Err("PAF record has fewer than 12 fields".to_owned());
    }
    let query_interval = Interval::new(
        fields[2].parse().map_err(|_| "invalid PAF query start")?,
        fields[3].parse().map_err(|_| "invalid PAF query end")?,
    )
    .ok_or_else(|| "invalid PAF query interval".to_owned())?;
    let target_interval = Interval::new(
        fields[7].parse().map_err(|_| "invalid PAF target start")?,
        fields[8].parse().map_err(|_| "invalid PAF target end")?,
    )
    .ok_or_else(|| "invalid PAF target interval".to_owned())?;
    let tags = &fields[12..];
    Ok(PafRecord {
        query_id: fields[0].to_owned(),
        query_len: fields[1].parse().map_err(|_| "invalid PAF query length")?,
        query_interval,
        strand: Strand::parse(fields[4]).ok_or_else(|| "invalid PAF strand".to_owned())?,
        target_id: fields[5].to_owned(),
        target_len: fields[6].parse().map_err(|_| "invalid PAF target length")?,
        target_interval,
        matches_nt: fields[9].parse().map_err(|_| "invalid PAF matches")?,
        block_len_nt: fields[10].parse().map_err(|_| "invalid PAF block length")?,
        mapq: fields[11].parse().map_err(|_| "invalid PAF mapq")?,
        alignment_score: paf_tag(tags.iter().copied(), "AS").and_then(|value| value.parse().ok()),
        frameshifts: paf_tag(tags.iter().copied(), "fs")
            .and_then(|value| value.parse().ok())
            .unwrap_or(0),
        inframe_stops: paf_tag(tags.iter().copied(), "st")
            .and_then(|value| value.parse().ok())
            .unwrap_or(0),
        cigar: paf_tag(tags.iter().copied(), "cg").map(str::to_owned),
        cs: paf_tag(tags.iter().copied(), "cs").map(str::to_owned),
    })
}

fn paf_compatible(mrna: &MrnaFeature, paf: &PafRecord) -> bool {
    paf.query_id == mrna.protein_id
        && paf.target_id == mrna.candidate_id
        && paf.strand == mrna.strand
        && paf.query_interval == mrna.target_interval_aa
        && endpoint_distance(paf.target_interval, mrna.genomic_interval) <= 3
}

fn endpoint_distance(left: Interval, right: Interval) -> usize {
    left.start.abs_diff(right.start) + left.end.abs_diff(right.end)
}

/// Parse miniprot `--gff` output, including embedded `##PAF` records.
pub fn parse_miniprot_output(output: &str) -> ParseReport {
    let mut pafs = Vec::new();
    let mut mrnas = Vec::new();
    let mut cds_by_parent: BTreeMap<String, Vec<CdsFeature>> = BTreeMap::new();
    let mut warnings = Vec::new();
    for (line_number, line) in output.lines().enumerate() {
        if line.starts_with("##PAF\t") {
            match parse_paf_line(line) {
                Ok(record) => pafs.push(record),
                Err(error) => warnings.push(format!("line {}: {error}", line_number + 1)),
            }
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 9 {
            warnings.push(format!(
                "line {}: invalid GFF3 field count",
                line_number + 1
            ));
            continue;
        }
        let Some(strand) = Strand::parse(fields[6]) else {
            warnings.push(format!("line {}: invalid GFF3 strand", line_number + 1));
            continue;
        };
        let Some(interval) = parse_gff_interval(fields[3], fields[4]) else {
            warnings.push(format!("line {}: invalid GFF3 interval", line_number + 1));
            continue;
        };
        let attributes = parse_attributes(fields[8]);
        match fields[2] {
            "mRNA" => {
                let Some(id) = attributes.get("ID").copied() else {
                    warnings.push(format!("line {}: mRNA lacks ID", line_number + 1));
                    continue;
                };
                let Some((protein_id, target_interval_aa)) =
                    parse_target(attributes.get("Target").copied())
                else {
                    warnings.push(format!("line {}: mRNA lacks valid Target", line_number + 1));
                    continue;
                };
                mrnas.push(MrnaFeature {
                    id: id.to_owned(),
                    candidate_id: fields[0].to_owned(),
                    protein_id,
                    target_interval_aa,
                    strand,
                    genomic_interval: interval,
                    score: (fields[5] != ".").then(|| fields[5].parse().ok()).flatten(),
                    identity_millionths: decimal_millionths(attributes.get("Identity").copied()),
                    positive_millionths: decimal_millionths(attributes.get("Positive").copied()),
                    frameshifts: attributes
                        .get("Frameshift")
                        .and_then(|value| value.parse().ok())
                        .unwrap_or(0),
                    inframe_stops: attributes
                        .get("StopCodon")
                        .and_then(|value| value.parse().ok())
                        .unwrap_or(0),
                    rank: attributes.get("Rank").and_then(|value| value.parse().ok()),
                });
            }
            "CDS" => {
                let Some(parent) = attributes.get("Parent").copied() else {
                    warnings.push(format!("line {}: CDS lacks Parent", line_number + 1));
                    continue;
                };
                let phase = match fields[7] {
                    "." => None,
                    "0" => Some(0),
                    "1" => Some(1),
                    "2" => Some(2),
                    _ => {
                        warnings.push(format!("line {}: invalid CDS phase", line_number + 1));
                        None
                    }
                };
                cds_by_parent
                    .entry(parent.to_owned())
                    .or_default()
                    .push(CdsFeature {
                        interval,
                        strand,
                        phase,
                    });
            }
            _ => {}
        }
    }
    let mut used_pafs = vec![false; pafs.len()];
    let mut models = Vec::new();
    for mrna in mrnas {
        let best_paf = pafs
            .iter()
            .enumerate()
            .filter(|(index, paf)| !used_pafs[*index] && paf_compatible(&mrna, paf))
            .min_by_key(|(_, paf)| endpoint_distance(paf.target_interval, mrna.genomic_interval))
            .map(|(index, _)| index);
        if let Some(index) = best_paf {
            used_pafs[index] = true;
        }
        models.push(RawGeneModel {
            gff_id: mrna.id.clone(),
            candidate_id: mrna.candidate_id,
            protein_id: mrna.protein_id,
            target_interval_aa: mrna.target_interval_aa,
            strand: mrna.strand,
            genomic_interval: mrna.genomic_interval,
            gff_score: mrna.score,
            identity_millionths: mrna.identity_millionths,
            positive_millionths: mrna.positive_millionths,
            frameshifts: mrna.frameshifts,
            inframe_stops: mrna.inframe_stops,
            rank: mrna.rank,
            cdss: cds_by_parent.remove(&mrna.id).unwrap_or_default(),
            paf: best_paf.map(|index| pafs[index].clone()),
        });
    }
    ParseReport {
        models,
        unassigned_paf: used_pafs.into_iter().filter(|used| !used).count(),
        warnings,
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CigarOperation {
    pub length: usize,
    pub op: char,
}

pub fn parse_protein_cigar(value: &str) -> Result<Vec<CigarOperation>, String> {
    let mut operations = Vec::new();
    let mut length = 0usize;
    for byte in value.bytes() {
        if byte.is_ascii_digit() {
            length = length
                .checked_mul(10)
                .and_then(|value| value.checked_add((byte - b'0') as usize))
                .ok_or_else(|| "protein CIGAR length overflow".to_owned())?;
            continue;
        }
        let op = byte as char;
        if length == 0 || !matches!(op, 'M' | 'I' | 'D' | 'F' | 'G' | 'N' | 'U' | 'V') {
            return Err(format!("invalid protein CIGAR operation: {length}{op}"));
        }
        operations.push(CigarOperation { length, op });
        length = 0;
    }
    if length != 0 {
        return Err("protein CIGAR ends without an operation".to_owned());
    }
    Ok(operations)
}

pub fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' | b'U' => 'A',
            b'R' => 'Y',
            b'Y' => 'R',
            b'S' => 'S',
            b'W' => 'W',
            b'K' => 'M',
            b'M' => 'K',
            b'B' => 'V',
            b'D' => 'H',
            b'H' => 'D',
            b'V' => 'B',
            b'N' => 'N',
            b'-' => '-',
            _ => 'N',
        })
        .collect()
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TranslationError {
    NonTripletLength(usize),
}

impl fmt::Display for TranslationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonTripletLength(length) => {
                write!(formatter, "CDS length {length} is not divisible by 3")
            }
        }
    }
}

pub fn translate_cds(sequence: &str) -> Result<String, TranslationError> {
    if !sequence.len().is_multiple_of(3) {
        return Err(TranslationError::NonTripletLength(sequence.len()));
    }
    Ok(sequence.as_bytes().chunks_exact(3).map(codon_aa).collect())
}

fn codon_aa(codon: &[u8]) -> char {
    fn base(value: u8) -> Option<usize> {
        match value.to_ascii_uppercase() {
            b'T' | b'U' => Some(0),
            b'C' => Some(1),
            b'A' => Some(2),
            b'G' => Some(3),
            _ => None,
        }
    }
    let (Some(a), Some(b), Some(c)) = (base(codon[0]), base(codon[1]), base(codon[2])) else {
        return 'X';
    };
    const TABLE: &str = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    TABLE.as_bytes()[a * 16 + b * 4 + c] as char
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preserves_ambiguity_and_length() {
        assert_eq!(normalize_dna("aUnR-y").unwrap(), "ATNR-Y");
    }

    #[test]
    fn parses_miniprot_paf_tags() {
        let paf = parse_paf_line("p\t100\t5\t95\t+\tc\t300\t10\t280\t240\t270\t60\tAS:i:400\tfs:i:1\tst:i:2\tcg:Z:30M1F20M").unwrap();
        assert_eq!(paf.query_interval, Interval { start: 5, end: 95 });
        assert_eq!(paf.alignment_score, Some(400));
        assert_eq!(paf.frameshifts, 1);
        assert_eq!(paf.inframe_stops, 2);
    }

    #[test]
    fn joins_gff_to_compatible_paf() {
        let input = "##PAF\tp\t10\t0\t10\t+\tc\t30\t0\t30\t30\t30\t60\tAS:i:42\tfs:i:0\tst:i:0\tcg:Z:10M\n##gff-version 3\nc\tminiprot\tmRNA\t1\t30\t42\t+\t.\tID=MP1;Target=p 1 10;Identity=1.0;Positive=1.0;Rank=0\nc\tminiprot\tCDS\t1\t30\t42\t+\t0\tParent=MP1\n";
        let report = parse_miniprot_output(input);
        assert!(report.warnings.is_empty());
        assert_eq!(report.models.len(), 1);
        assert_eq!(
            report.models[0].paf.as_ref().unwrap().alignment_score,
            Some(42)
        );
        assert_eq!(report.models[0].cdss.len(), 1);
    }

    #[test]
    fn does_not_reuse_one_paf_for_two_gff_models() {
        let input = "##PAF\tp\t10\t0\t10\t+\tc\t30\t0\t30\t30\t30\t60\tAS:i:42\tfs:i:0\tst:i:0\tcg:Z:10M\n##gff-version 3\nc\tminiprot\tmRNA\t1\t30\t42\t+\t.\tID=MP1;Target=p 1 10\nc\tminiprot\tCDS\t1\t30\t42\t+\t0\tParent=MP1\nc\tminiprot\tmRNA\t1\t30\t41\t+\t.\tID=MP2;Target=p 1 10\nc\tminiprot\tCDS\t1\t30\t41\t+\t0\tParent=MP2\n";
        let report = parse_miniprot_output(input);
        assert_eq!(report.models.len(), 2);
        assert_eq!(
            report
                .models
                .iter()
                .filter(|model| model.paf.is_some())
                .count(),
            1
        );
        assert_eq!(report.unassigned_paf, 0);
    }

    #[test]
    fn pairs_nested_gff_models_with_matching_paf_endpoints() {
        let input = "##PAF\tp\t20\t0\t10\t+\tc\t200\t100\t130\t30\t30\t60\tAS:i:10\n##PAF\tp\t20\t0\t20\t+\tc\t200\t90\t150\t60\t60\t60\tAS:i:20\n##gff-version 3\nc\tminiprot\tmRNA\t101\t130\t10\t+\t.\tID=MP1;Target=p 1 10\nc\tminiprot\tCDS\t101\t130\t10\t+\t0\tParent=MP1\nc\tminiprot\tmRNA\t91\t150\t20\t+\t.\tID=MP2;Target=p 1 20\nc\tminiprot\tCDS\t91\t150\t20\t+\t0\tParent=MP2\n";
        let report = parse_miniprot_output(input);
        assert_eq!(report.models.len(), 2);
        assert_eq!(
            report.models[0].paf.as_ref().unwrap().alignment_score,
            Some(10)
        );
        assert_eq!(
            report.models[1].paf.as_ref().unwrap().alignment_score,
            Some(20)
        );
        assert_eq!(report.unassigned_paf, 0);
    }

    #[test]
    fn parses_all_miniprot_protein_cigar_operations() {
        let operations = parse_protein_cigar("3M2I4D1F1G5N2U2V").unwrap();
        assert_eq!(operations.len(), 8);
        assert_eq!(operations[3], CigarOperation { length: 1, op: 'F' });
        assert_eq!(operations[7], CigarOperation { length: 2, op: 'V' });
    }

    #[test]
    fn strict_translation_refuses_partial_codon() {
        assert_eq!(translate_cds("ATGAAA").unwrap(), "MK");
        assert_eq!(
            translate_cds("ATGA"),
            Err(TranslationError::NonTripletLength(4))
        );
    }
}
