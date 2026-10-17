use super::{reverse_complement, translate_cds, Interval, RawGeneModel, Strand};

use super::evidence::AlignmentEvidence;

#[derive(Clone, Copy, Debug)]
pub struct ModelConfig {
    pub minimum_coverage: f64,
    pub complete_coverage: f64,
}

impl Default for ModelConfig {
    fn default() -> Self {
        Self {
            minimum_coverage: 0.20,
            complete_coverage: 0.80,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum ModelState {
    Complete,
    TerminalPartial,
    InternalPartial,
    LowCoverage,
    Frameshifted,
    InternalStop,
    PhaseInconsistent,
    AmbiguousCds,
    AmbiguousModel,
}

impl ModelState {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Complete => "complete",
            Self::TerminalPartial => "terminal_partial",
            Self::InternalPartial => "internal_partial",
            Self::LowCoverage => "low_coverage",
            Self::Frameshifted => "frameshifted",
            Self::InternalStop => "internal_stop",
            Self::PhaseInconsistent => "phase_inconsistent",
            Self::AmbiguousCds => "ambiguous_cds",
            Self::AmbiguousModel => "ambiguous_model",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum SpliceClass {
    Major,
    Minor,
    Noncanonical,
    Unavailable,
}

impl SpliceClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Major => "major",
            Self::Minor => "minor",
            Self::Noncanonical => "noncanonical",
            Self::Unavailable => "unavailable",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ExonSequence {
    pub index: usize,
    pub interval: Interval,
    pub phase: Option<u8>,
    pub sequence: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IntronSequence {
    pub index: usize,
    pub interval: Interval,
    pub sequence: String,
    pub donor: String,
    pub acceptor: String,
    pub splice_class: SpliceClass,
}

#[derive(Clone, Debug)]
pub struct GeneModel {
    pub engine: Option<super::evidence::AlignmentEngine>,
    pub gff_id: String,
    pub candidate_id: String,
    pub protein_id: String,
    pub protein_len_aa: usize,
    pub query_interval_aa: Interval,
    pub genomic_interval: Interval,
    pub strand: Strand,
    pub alignment_score: i64,
    pub identity_millionths: u32,
    pub positive_millionths: u32,
    pub coverage_millionths: u32,
    pub rank: usize,
    pub frameshifts: usize,
    pub inframe_stops: usize,
    pub phase_valid: bool,
    pub five_prime_trim_nt: usize,
    pub three_prime_trim_nt: usize,
    pub large_query_gaps: usize,
    pub exons: Vec<ExonSequence>,
    pub introns: Vec<IntronSequence>,
    pub cds: String,
    pub query_cds_offsets: Option<Vec<Option<usize>>>,
    pub protein: Option<String>,
    pub gene: String,
    pub state: ModelState,
    pub qc_flags: Vec<String>,
    pub selected: bool,
    pub eligible_for_resolve: bool,
    pub competition_group: usize,
}

impl GeneModel {
    /// A reference beginning with M does not justify inventing a genomic start
    /// codon. Report the mismatch and keep an otherwise clean model partial.
    pub fn check_reference_start(&mut self, reference_start: Option<u8>) {
        if self.query_interval_aa.start == 0
            && reference_start.is_some_and(|r| r.eq_ignore_ascii_case(&b'M'))
            && self.protein.as_ref().is_some_and(|p| !p.starts_with('M'))
        {
            if !self
                .qc_flags
                .iter()
                .any(|f| f == "reference_start_mismatch")
            {
                self.qc_flags.push("reference_start_mismatch".into());
            }
            if self.state == ModelState::Complete {
                self.state = ModelState::TerminalPartial;
            }
        }
    }

    pub fn coverage(&self) -> f64 {
        self.coverage_millionths as f64 / 1_000_000.0
    }

    pub fn identity(&self) -> f64 {
        self.identity_millionths as f64 / 1_000_000.0
    }

    pub fn positive_fraction(&self) -> f64 {
        self.positive_millionths as f64 / 1_000_000.0
    }

    pub fn noncanonical_introns(&self) -> usize {
        self.introns
            .iter()
            .filter(|intron| intron.splice_class == SpliceClass::Noncanonical)
            .count()
    }

    pub fn structural_defect_count(&self) -> usize {
        usize::from(self.frameshifts > 0)
            + usize::from(self.inframe_stops > 0 || self.state == ModelState::InternalStop)
            + usize::from(!self.phase_valid)
            + usize::from(
                self.protein
                    .as_ref()
                    .is_some_and(|protein| protein.contains('X')),
            )
            + self.noncanonical_introns()
            + self
                .introns
                .iter()
                .filter(|i| i.splice_class == SpliceClass::Unavailable)
                .count()
    }
}

fn oriented_slice(sequence: &str, interval: Interval, strand: Strand) -> Result<String, String> {
    let value = sequence.get(interval.start..interval.end).ok_or_else(|| {
        format!(
            "interval {}..{} is outside candidate",
            interval.start, interval.end
        )
    })?;
    Ok(if strand == Strand::Reverse {
        reverse_complement(value)
    } else {
        value.to_owned()
    })
}

fn transcript_ordered_cdss(
    raw: &RawGeneModel,
    sequence_len: usize,
) -> Result<Vec<super::CdsFeature>, String> {
    let mut cdss = raw.cdss.clone();
    if cdss.is_empty() {
        return Err("model contains no CDS features".to_owned());
    }
    if cdss.iter().any(|cds| {
        cds.strand != raw.strand
            || cds.interval.end > sequence_len
            || cds.interval.start >= cds.interval.end
    }) {
        return Err("CDS strand or interval is inconsistent with its mRNA".to_owned());
    }
    cdss.sort_by_key(|cds| cds.interval.start);
    if cdss
        .windows(2)
        .any(|pair| pair[0].interval.end > pair[1].interval.start)
    {
        return Err("model contains overlapping CDS features".into());
    }
    if raw.strand == Strand::Reverse {
        cdss.reverse();
    }
    Ok(cdss)
}

fn phase_is_valid(cdss: &[super::CdsFeature], query_starts_at_origin: bool) -> bool {
    if cdss.iter().any(|cds| cds.phase.is_none()) {
        return false;
    }
    let first_phase = cdss[0].phase.unwrap() as usize;
    if query_starts_at_origin && first_phase != 0 {
        return false;
    }
    let mut coding_length = cdss[0].interval.len().saturating_sub(first_phase);
    for cds in &cdss[1..] {
        let expected = (3 - coding_length % 3) % 3;
        if cds.phase.map(usize::from) != Some(expected) {
            return false;
        }
        coding_length += cds.interval.len();
    }
    true
}

fn intron_interval(current: Interval, next: Interval, strand: Strand) -> Option<Interval> {
    match strand {
        Strand::Forward => Interval::new(current.end, next.start),
        Strand::Reverse => Interval::new(next.end, current.start),
    }
}

fn splice_metadata(sequence: &str) -> (String, String, SpliceClass) {
    if sequence.len() < 4 {
        return (String::new(), String::new(), SpliceClass::Unavailable);
    }
    let donor = sequence[..2].to_owned();
    let acceptor = sequence[sequence.len() - 2..].to_owned();
    let class = match (donor.as_str(), acceptor.as_str()) {
        ("GT", "AG") => SpliceClass::Major,
        ("GC", "AG") | ("AT", "AC") => SpliceClass::Minor,
        _ => SpliceClass::Noncanonical,
    };
    (donor, acceptor, class)
}

pub fn build_gene_model(
    raw: &RawGeneModel,
    candidate_sequence: &str,
    fallback_protein_len: Option<usize>,
    config: ModelConfig,
) -> Result<GeneModel, String> {
    let evidence = match &raw.evidence {
        Some(value) => Some(value.clone()),
        None => AlignmentEvidence::from_miniprot(raw)?,
    };
    build_gene_model_with_evidence(
        raw,
        evidence.as_ref(),
        candidate_sequence,
        fallback_protein_len,
        config,
    )
}

/// Both adapters use the same sequence extraction and structural acceptance rules.
pub fn build_gene_model_with_evidence(
    raw: &RawGeneModel,
    evidence: Option<&AlignmentEvidence>,
    candidate_sequence: &str,
    fallback_protein_len: Option<usize>,
    config: ModelConfig,
) -> Result<GeneModel, String> {
    let cdss = transcript_ordered_cdss(raw, candidate_sequence.len())?;
    let protein_len = evidence
        .map(|value| value.query_len)
        .or(fallback_protein_len)
        .ok_or_else(|| "protein length is unavailable".to_owned())?;
    if protein_len == 0 {
        return Err("protein length is zero".to_owned());
    }
    let query_interval = evidence
        .map(|value| value.query_interval)
        .unwrap_or(raw.target_interval_aa);
    if query_interval.end > protein_len {
        return Err("protein interval exceeds protein length".to_owned());
    }
    let frameshifts = raw
        .frameshifts
        .max(evidence.map(|v| v.frameshifts).unwrap_or(0));
    let inframe_stops = raw
        .inframe_stops
        .max(evidence.map(|v| v.inframe_stops).unwrap_or(0));
    let large_query_gaps = evidence.map(|v| v.large_query_gaps).unwrap_or(0);

    let mut exons = Vec::new();
    let mut observed_cds = String::new();
    for (index, cds) in cdss.iter().enumerate() {
        let sequence = oriented_slice(candidate_sequence, cds.interval, raw.strand)?;
        observed_cds.push_str(&sequence);
        exons.push(ExonSequence {
            index: index + 1,
            interval: cds.interval,
            phase: cds.phase,
            sequence,
        });
    }
    let phase_valid = phase_is_valid(&cdss, query_interval.start == 0);
    let five_prime_trim_nt = if query_interval.start > 0 {
        cdss[0].phase.unwrap_or(0) as usize
    } else {
        0
    };
    if five_prime_trim_nt > observed_cds.len() {
        return Err("5-prime phase trim exceeds CDS length".to_owned());
    }
    let mut cds = observed_cds[five_prime_trim_nt..].to_owned();
    let remainder = cds.len() % 3;
    let three_prime_trim_nt = if remainder > 0 && query_interval.end < protein_len {
        remainder
    } else {
        0
    };
    if three_prime_trim_nt > 0 {
        cds.truncate(cds.len() - three_prime_trim_nt);
    }
    let phase_valid = phase_valid && cds.len().is_multiple_of(3);
    let protein = phase_valid.then(|| translate_cds(&cds).ok()).flatten();

    let mut introns = Vec::new();
    for index in 0..cdss.len().saturating_sub(1) {
        let Some(interval) =
            intron_interval(cdss[index].interval, cdss[index + 1].interval, raw.strand)
        else {
            continue;
        };
        let sequence = oriented_slice(candidate_sequence, interval, raw.strand)?;
        let (donor, acceptor, splice_class) = splice_metadata(&sequence);
        introns.push(IntronSequence {
            index: index + 1,
            interval,
            sequence,
            donor,
            acceptor,
            splice_class,
        });
    }
    let genomic_start = cdss.iter().map(|cds| cds.interval.start).min().unwrap();
    let genomic_end = cdss.iter().map(|cds| cds.interval.end).max().unwrap();
    let genomic_interval = Interval::new(genomic_start, genomic_end).unwrap();
    let gene = oriented_slice(candidate_sequence, genomic_interval, raw.strand)?;
    let query_gap_aa = evidence.map(|e| e.query_gap_aa).unwrap_or(0);
    let covered_aa = query_interval
        .len()
        .checked_sub(query_gap_aa)
        .ok_or("query_gap_exceeds_span")?;
    let coverage = covered_aa as f64 / protein_len as f64;
    let coverage_millionths = (coverage.min(1.0) * 1_000_000.0).round() as u32;
    let reaches_ends = query_interval.start == 0 && query_interval.end == protein_len;
    let has_ambiguous_codon = protein.as_ref().is_some_and(|value| value.contains('X'));
    let translated_internal_stop = protein.as_ref().is_some_and(|value| {
        value
            .strip_suffix('*')
            .unwrap_or(value.as_str())
            .contains('*')
    });
    let query_cds_offsets =
        evidence
            .and_then(|e| e.coding_map.as_ref())
            .and_then(|(offsets, length)| {
                let terminal_stop = protein.as_ref().is_some_and(|p| p.ends_with('*'));
                // Partial phase trimming must not silently shift alignment coordinates.
                // Only retain maps whose coding consumption matches the extracted CDS.
                (*length == cds.len()
                    || (terminal_stop && length.checked_add(3) == Some(cds.len())))
                .then(|| offsets.clone())
            });
    let state = if frameshifts > 0 {
        ModelState::Frameshifted
    } else if inframe_stops > 0 || translated_internal_stop {
        ModelState::InternalStop
    } else if !phase_valid {
        ModelState::PhaseInconsistent
    } else if has_ambiguous_codon {
        ModelState::AmbiguousCds
    } else if coverage < config.minimum_coverage {
        ModelState::LowCoverage
    } else if query_gap_aa > 0 {
        ModelState::InternalPartial
    } else if coverage >= config.complete_coverage && reaches_ends {
        ModelState::Complete
    } else {
        ModelState::TerminalPartial
    };
    let mut qc_flags = Vec::new();
    if query_interval.start > 0 {
        qc_flags.push("query_start_truncated".into());
    }
    if query_interval.end < protein_len {
        qc_flags.push("query_end_truncated".into());
    }
    if evidence.is_none() {
        qc_flags.push("missing_paf".to_owned());
    }
    if introns
        .iter()
        .any(|intron| intron.splice_class == SpliceClass::Noncanonical)
    {
        qc_flags.push("noncanonical_splice".to_owned());
    }
    if introns
        .iter()
        .any(|intron| intron.splice_class == SpliceClass::Unavailable)
    {
        qc_flags.push("invalid_splice_length".to_owned());
    }
    if query_gap_aa > 0 {
        qc_flags.push(format!("query_gap_aa={query_gap_aa}"));
    }
    if large_query_gaps > 0 {
        qc_flags.push("large_query_gap".to_owned());
    }
    let structurally_eligible = matches!(
        state,
        ModelState::Complete | ModelState::TerminalPartial | ModelState::InternalPartial
    ) && evidence.is_some()
        && !qc_flags
            .iter()
            .any(|flag| flag == "noncanonical_splice" || flag == "invalid_splice_length");
    Ok(GeneModel {
        engine: evidence.map(|e| e.engine),
        gff_id: raw.gff_id.clone(),
        candidate_id: raw.candidate_id.clone(),
        protein_id: raw.protein_id.clone(),
        protein_len_aa: protein_len,
        query_interval_aa: query_interval,
        genomic_interval,
        strand: raw.strand,
        alignment_score: evidence
            .and_then(|value| value.alignment_score)
            .or(raw.gff_score)
            .unwrap_or(0),
        identity_millionths: raw.identity_millionths.unwrap_or(0),
        positive_millionths: raw.positive_millionths.unwrap_or(0),
        coverage_millionths,
        rank: raw.rank.unwrap_or(usize::MAX),
        frameshifts,
        inframe_stops,
        phase_valid,
        five_prime_trim_nt,
        three_prime_trim_nt,
        large_query_gaps,
        exons,
        introns,
        cds,
        query_cds_offsets,
        protein,
        gene,
        state,
        qc_flags,
        selected: false,
        eligible_for_resolve: structurally_eligible,
        competition_group: 0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gene_annotation::{CdsFeature, PafRecord, RawGeneModel};

    fn raw_model(strand: Strand, cdss: Vec<CdsFeature>) -> RawGeneModel {
        RawGeneModel {
            gff_id: "MP1".into(),
            candidate_id: "candidate_1".into(),
            protein_id: "protein_1".into(),
            target_interval_aa: Interval::new(0, 4).unwrap(),
            strand,
            genomic_interval: Interval::new(0, 16).unwrap(),
            gff_score: Some(100),
            identity_millionths: Some(900_000),
            positive_millionths: Some(950_000),
            frameshifts: 0,
            inframe_stops: 0,
            rank: Some(0),
            cdss,
            evidence: None,
            paf: Some(PafRecord {
                query_id: "protein_1".into(),
                query_len: 4,
                query_interval: Interval::new(0, 4).unwrap(),
                strand,
                target_id: "candidate_1".into(),
                target_len: 16,
                target_interval: Interval::new(0, 16).unwrap(),
                matches_nt: 12,
                block_len_nt: 12,
                mapq: 60,
                alignment_score: Some(100),
                frameshifts: 0,
                inframe_stops: 0,
                cigar: Some("2M4N2M".into()),
                cs: None,
            }),
        }
    }

    #[test]
    fn missing_alignment_evidence_cannot_pass_structural_acceptance() {
        let mut raw = raw_model(
            Strand::Forward,
            vec![CdsFeature {
                interval: Interval::new(0, 12).unwrap(),
                strand: Strand::Forward,
                phase: Some(0),
            }],
        );
        raw.paf = None;
        let model =
            build_gene_model(&raw, "ATGAAATTTCCC", Some(4), ModelConfig::default()).unwrap();
        assert_eq!(model.state, ModelState::Complete);
        assert!(!model.eligible_for_resolve);
        assert!(model.qc_flags.iter().any(|f| f == "missing_paf"));
    }

    #[test]
    fn extracts_forward_exons_and_intron_separately() {
        let raw = raw_model(
            Strand::Forward,
            vec![
                CdsFeature {
                    interval: Interval::new(0, 6).unwrap(),
                    strand: Strand::Forward,
                    phase: Some(0),
                },
                CdsFeature {
                    interval: Interval::new(10, 16).unwrap(),
                    strand: Strand::Forward,
                    phase: Some(0),
                },
            ],
        );
        let model =
            build_gene_model(&raw, "ATGAAAGTAGTTTCCC", None, ModelConfig::default()).unwrap();
        assert_eq!(model.cds, "ATGAAATTTCCC");
        assert_eq!(model.introns[0].sequence, "GTAG");
        assert_eq!(model.introns[0].splice_class, SpliceClass::Major);
        assert_eq!(model.state, ModelState::Complete);
    }

    #[test]
    fn extracts_reverse_model_in_transcript_orientation() {
        let raw = raw_model(
            Strand::Reverse,
            vec![
                CdsFeature {
                    interval: Interval::new(0, 6).unwrap(),
                    strand: Strand::Reverse,
                    phase: Some(0),
                },
                CdsFeature {
                    interval: Interval::new(10, 16).unwrap(),
                    strand: Strand::Reverse,
                    phase: Some(0),
                },
            ],
        );
        let model =
            build_gene_model(&raw, "GGGAAACTACATTTTC", None, ModelConfig::default()).unwrap();
        assert_eq!(model.cds, "GAAAATTTTCCC");
        assert_eq!(model.introns[0].sequence, "GTAG");
        assert_eq!(model.strand, Strand::Reverse);
    }

    #[test]
    fn extracts_three_reverse_strand_exons_with_split_codon_phase() {
        let transcript_gene = "ATGAAAGTAAAGTTTCGTAGCCGGG";
        let candidate = reverse_complement(transcript_gene);
        let mut raw = raw_model(
            Strand::Reverse,
            vec![
                CdsFeature {
                    interval: Interval::new(0, 5).unwrap(),
                    strand: Strand::Reverse,
                    phase: Some(2),
                },
                CdsFeature {
                    interval: Interval::new(9, 13).unwrap(),
                    strand: Strand::Reverse,
                    phase: Some(0),
                },
                CdsFeature {
                    interval: Interval::new(19, 25).unwrap(),
                    strand: Strand::Reverse,
                    phase: Some(0),
                },
            ],
        );
        raw.genomic_interval = Interval::new(0, 25).unwrap();
        raw.target_interval_aa = Interval::new(0, 5).unwrap();
        let paf = raw.paf.as_mut().unwrap();
        paf.query_len = 5;
        paf.query_interval = Interval::new(0, 5).unwrap();
        paf.target_len = 25;
        paf.target_interval = Interval::new(0, 25).unwrap();
        paf.cigar = Some("2M6N1M7V1M".into());
        let model = build_gene_model(&raw, &candidate, None, ModelConfig::default()).unwrap();
        assert_eq!(model.cds, "ATGAAATTTCCCGGG");
        assert_eq!(
            model
                .introns
                .iter()
                .map(|intron| intron.sequence.as_str())
                .collect::<Vec<_>>(),
            vec!["GTAAAG", "GTAG"]
        );
        assert!(model.phase_valid);
        assert_eq!(model.state, ModelState::Complete);
    }

    #[test]
    fn frameshift_is_never_resolve_eligible() {
        let mut raw = raw_model(
            Strand::Forward,
            vec![CdsFeature {
                interval: Interval::new(0, 12).unwrap(),
                strand: Strand::Forward,
                phase: Some(0),
            }],
        );
        raw.paf.as_mut().unwrap().cigar = Some("2M1F2M".into());
        let model = build_gene_model(&raw, "ATGAAATTTCCC", None, ModelConfig::default()).unwrap();
        assert_eq!(model.state, ModelState::Frameshifted);
        assert!(!model.eligible_for_resolve);
    }

    #[test]
    fn trims_five_prime_partial_cds_from_annotated_phase() {
        let mut raw = raw_model(
            Strand::Forward,
            vec![CdsFeature {
                interval: Interval::new(0, 10).unwrap(),
                strand: Strand::Forward,
                phase: Some(1),
            }],
        );
        raw.genomic_interval = Interval::new(0, 10).unwrap();
        raw.target_interval_aa = Interval::new(1, 4).unwrap();
        let paf = raw.paf.as_mut().unwrap();
        paf.query_len = 5;
        paf.query_interval = Interval::new(1, 4).unwrap();
        paf.target_len = 10;
        paf.target_interval = Interval::new(0, 10).unwrap();
        paf.cigar = Some("3M".into());
        let model = build_gene_model(&raw, "AATGAAATTT", None, ModelConfig::default()).unwrap();
        assert_eq!(model.cds, "ATGAAATTT");
        assert_eq!(model.five_prime_trim_nt, 1);
        assert_eq!(model.three_prime_trim_nt, 0);
        assert_eq!(model.state, ModelState::TerminalPartial);
        assert!(model.eligible_for_resolve);
    }

    #[test]
    fn trims_three_prime_partial_cds_to_complete_codons() {
        let mut raw = raw_model(
            Strand::Forward,
            vec![CdsFeature {
                interval: Interval::new(0, 10).unwrap(),
                strand: Strand::Forward,
                phase: Some(0),
            }],
        );
        raw.genomic_interval = Interval::new(0, 10).unwrap();
        raw.target_interval_aa = Interval::new(0, 3).unwrap();
        let paf = raw.paf.as_mut().unwrap();
        paf.query_len = 5;
        paf.query_interval = Interval::new(0, 3).unwrap();
        paf.target_len = 10;
        paf.target_interval = Interval::new(0, 10).unwrap();
        paf.cigar = Some("3M".into());
        let model = build_gene_model(&raw, "ATGAAATTTC", None, ModelConfig::default()).unwrap();
        assert_eq!(model.cds, "ATGAAATTT");
        assert_eq!(model.five_prime_trim_nt, 0);
        assert_eq!(model.three_prime_trim_nt, 1);
        assert_eq!(model.state, ModelState::TerminalPartial);
    }

    #[test]
    fn validates_a_codon_split_across_two_exons() {
        let mut raw = raw_model(
            Strand::Forward,
            vec![
                CdsFeature {
                    interval: Interval::new(0, 4).unwrap(),
                    strand: Strand::Forward,
                    phase: Some(0),
                },
                CdsFeature {
                    interval: Interval::new(8, 13).unwrap(),
                    strand: Strand::Forward,
                    phase: Some(2),
                },
            ],
        );
        raw.genomic_interval = Interval::new(0, 13).unwrap();
        let paf = raw.paf.as_mut().unwrap();
        paf.query_len = 3;
        paf.query_interval = Interval::new(0, 3).unwrap();
        paf.target_len = 13;
        paf.target_interval = Interval::new(0, 13).unwrap();
        let model = build_gene_model(&raw, "ATGAGTAGAATTT", None, ModelConfig::default()).unwrap();
        assert_eq!(model.cds, "ATGAAATTT");
        assert!(model.phase_valid);

        raw.cdss[1].phase = Some(1);
        let inconsistent =
            build_gene_model(&raw, "ATGAGTAGAATTT", None, ModelConfig::default()).unwrap();
        assert_eq!(inconsistent.state, ModelState::PhaseInconsistent);
        assert!(!inconsistent.eligible_for_resolve);
    }

    #[test]
    fn cigar_g_and_translated_stop_are_structural_failures() {
        let mut frameshift_raw = raw_model(
            Strand::Forward,
            vec![CdsFeature {
                interval: Interval::new(0, 12).unwrap(),
                strand: Strand::Forward,
                phase: Some(0),
            }],
        );
        frameshift_raw.paf.as_mut().unwrap().cigar = Some("2M1G2M".into());
        let frameshift = build_gene_model(
            &frameshift_raw,
            "ATGAAATTTCCC",
            None,
            ModelConfig::default(),
        )
        .unwrap();
        assert_eq!(frameshift.state, ModelState::Frameshifted);

        let stop_raw = raw_model(
            Strand::Forward,
            vec![CdsFeature {
                interval: Interval::new(0, 12).unwrap(),
                strand: Strand::Forward,
                phase: Some(0),
            }],
        );
        let stop =
            build_gene_model(&stop_raw, "ATGTAATTTCCC", None, ModelConfig::default()).unwrap();
        assert_eq!(stop.state, ModelState::InternalStop);
        assert!(!stop.eligible_for_resolve);
    }

    #[test]
    fn permits_a_single_terminal_stop_codon() {
        let raw = raw_model(
            Strand::Forward,
            vec![CdsFeature {
                interval: Interval::new(0, 12).unwrap(),
                strand: Strand::Forward,
                phase: Some(0),
            }],
        );
        let model = build_gene_model(&raw, "ATGAAATTTTAA", None, ModelConfig::default()).unwrap();
        assert_eq!(model.protein.as_deref(), Some("MKF*"));
        assert_eq!(model.state, ModelState::Complete);
        assert!(model.eligible_for_resolve);
    }

    #[test]
    fn missing_cds_phase_is_not_assumed_to_be_zero() {
        let raw = raw_model(
            Strand::Forward,
            vec![CdsFeature {
                interval: Interval::new(0, 12).unwrap(),
                strand: Strand::Forward,
                phase: None,
            }],
        );
        let model = build_gene_model(&raw, "ATGAAATTTCCC", None, ModelConfig::default()).unwrap();
        assert!(!model.phase_valid);
        assert_eq!(model.state, ModelState::PhaseInconsistent);
        assert!(!model.eligible_for_resolve);
    }

    #[test]
    fn miniprot_query_gap_reduces_coverage_and_prevents_complete_state() {
        let mut raw = raw_model(
            Strand::Forward,
            vec![CdsFeature {
                interval: Interval::new(0, 12).unwrap(),
                strand: Strand::Forward,
                phase: Some(0),
            }],
        );
        let paf = raw.paf.as_mut().unwrap();
        paf.query_len = 5;
        paf.query_interval = Interval::new(0, 5).unwrap();
        paf.cigar = Some("2M1I2M".into());
        let model = build_gene_model(&raw, "ATGAAATTTCCC", None, ModelConfig::default()).unwrap();
        assert_eq!(model.state, ModelState::InternalPartial);
        assert_eq!(model.coverage_millionths, 800_000);
        assert_eq!(
            model.query_cds_offsets,
            Some(vec![Some(0), Some(3), None, Some(6), Some(9)])
        );
    }

    #[test]
    fn incomplete_splice_signals_cannot_pass_qc_on_either_strand() {
        for length in 1..4 {
            for strand in [Strand::Forward, Strand::Reverse] {
                let forward = format!("ATGAAA{}TTTCCC", "G".repeat(length));
                let target = if strand == Strand::Reverse {
                    super::super::reverse_complement(&forward)
                } else {
                    forward
                };
                let raw = raw_model(
                    strand,
                    vec![
                        CdsFeature {
                            interval: Interval::new(0, 6).unwrap(),
                            strand,
                            phase: Some(0),
                        },
                        CdsFeature {
                            interval: Interval::new(6 + length, 12 + length).unwrap(),
                            strand,
                            phase: Some(0),
                        },
                    ],
                );
                let model = build_gene_model(&raw, &target, None, ModelConfig::default()).unwrap();
                assert!(!model.eligible_for_resolve);
                assert!(model.qc_flags.iter().any(|f| f == "invalid_splice_length"));
                assert_eq!(model.structural_defect_count(), 1);
            }
        }
    }

    #[test]
    fn classifies_minor_and_noncanonical_splice_sites() {
        assert_eq!(splice_metadata("GCAG").2, SpliceClass::Minor);
        assert_eq!(splice_metadata("ATAC").2, SpliceClass::Minor);
        assert_eq!(splice_metadata("CTAG").2, SpliceClass::Noncanonical);

        let raw = raw_model(
            Strand::Forward,
            vec![
                CdsFeature {
                    interval: Interval::new(0, 6).unwrap(),
                    strand: Strand::Forward,
                    phase: Some(0),
                },
                CdsFeature {
                    interval: Interval::new(10, 16).unwrap(),
                    strand: Strand::Forward,
                    phase: Some(0),
                },
            ],
        );
        let model =
            build_gene_model(&raw, "ATGAAACTAGTTTCCC", None, ModelConfig::default()).unwrap();
        assert_eq!(model.introns[0].splice_class, SpliceClass::Noncanonical);
        assert!(model.qc_flags.contains(&"noncanonical_splice".to_owned()));
        assert!(!model.eligible_for_resolve);
    }
    #[test]
    fn a_single_missing_terminal_residue_is_partial_even_at_high_coverage() {
        for (start, end, flag) in [
            (1, 100, "query_start_truncated"),
            (0, 99, "query_end_truncated"),
        ] {
            let mut raw = raw_model(
                Strand::Forward,
                vec![CdsFeature {
                    interval: Interval::new(0, 297).unwrap(),
                    strand: Strand::Forward,
                    phase: Some(0),
                }],
            );
            let paf = raw.paf.as_mut().unwrap();
            paf.query_len = 100;
            paf.query_interval = Interval::new(start, end).unwrap();
            paf.cigar = Some("99M".into());
            let model =
                build_gene_model(&raw, &"ATG".repeat(99), None, ModelConfig::default()).unwrap();
            assert_eq!(model.state, ModelState::TerminalPartial);
            assert!(model.eligible_for_resolve);
            assert!(model.qc_flags.iter().any(|f| f == flag));
        }
    }

    #[test]
    fn reference_start_mismatch_is_audited_without_inventing_sequence() {
        let raw = raw_model(
            Strand::Forward,
            vec![CdsFeature {
                interval: Interval::new(0, 12).unwrap(),
                strand: Strand::Forward,
                phase: Some(0),
            }],
        );
        let mut model =
            build_gene_model(&raw, "CCAAAATTTCCC", None, ModelConfig::default()).unwrap();
        model.check_reference_start(Some(b'P'));
        assert_eq!(model.state, ModelState::Complete);
        model.check_reference_start(Some(b'M'));
        model.check_reference_start(Some(b'M'));
        assert_eq!(model.state, ModelState::TerminalPartial);
        assert_eq!(model.cds, "CCAAAATTTCCC");
        assert_eq!(
            model
                .qc_flags
                .iter()
                .filter(|f| *f == "reference_start_mismatch")
                .count(),
            1
        );
        model.state = ModelState::Frameshifted;
        model.eligible_for_resolve = false;
        model.check_reference_start(Some(b'M'));
        assert_eq!(model.state, ModelState::Frameshifted);
        assert!(!model.eligible_for_resolve);
    }
    #[test]
    fn overlapping_cds_features_are_rejected_on_both_strands() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let raw = raw_model(
                strand,
                vec![
                    CdsFeature {
                        interval: Interval::new(0, 6).unwrap(),
                        strand,
                        phase: Some(0),
                    },
                    CdsFeature {
                        interval: Interval::new(3, 9).unwrap(),
                        strand,
                        phase: Some(0),
                    },
                ],
            );
            assert!(build_gene_model(&raw, "ATGAAATTT", None, ModelConfig::default()).is_err());
        }
    }
}
