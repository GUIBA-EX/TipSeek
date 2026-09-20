use super::{
    parse_protein_cigar, reverse_complement, translate_cds, Interval, RawGeneModel, Strand,
};

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
    pub protein: Option<String>,
    pub gene: String,
    pub state: ModelState,
    pub qc_flags: Vec<String>,
    pub selected: bool,
    pub eligible_for_resolve: bool,
    pub competition_group: usize,
}

impl GeneModel {
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
            + usize::from(self.inframe_stops > 0)
            + usize::from(!self.phase_valid)
            + usize::from(
                self.protein
                    .as_ref()
                    .is_some_and(|protein| protein.contains('X')),
            )
            + self.noncanonical_introns()
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
    let cdss = transcript_ordered_cdss(raw, candidate_sequence.len())?;
    let paf = raw.paf.as_ref();
    let protein_len = paf
        .map(|value| value.query_len)
        .or(fallback_protein_len)
        .ok_or_else(|| "protein length is unavailable".to_owned())?;
    if protein_len == 0 {
        return Err("protein length is zero".to_owned());
    }
    let query_interval = paf
        .map(|value| value.query_interval)
        .unwrap_or(raw.target_interval_aa);
    if query_interval.end > protein_len {
        return Err("protein interval exceeds protein length".to_owned());
    }
    let reported_frameshifts = raw
        .frameshifts
        .max(paf.map(|value| value.frameshifts).unwrap_or(0));
    let inframe_stops = raw
        .inframe_stops
        .max(paf.map(|value| value.inframe_stops).unwrap_or(0));
    let mut large_query_gaps = 0usize;
    let mut cigar_frameshifts = 0usize;
    if let Some(cigar) = paf.and_then(|value| value.cigar.as_deref()) {
        for operation in parse_protein_cigar(cigar)? {
            if matches!(operation.op, 'F' | 'G') {
                cigar_frameshifts += 1;
            }
            if operation.op == 'I' && operation.length >= 10 {
                large_query_gaps += 1;
            }
        }
    }
    let frameshifts = reported_frameshifts.max(cigar_frameshifts);

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
    let coverage = query_interval.len() as f64 / protein_len as f64;
    let coverage_millionths = (coverage.min(1.0) * 1_000_000.0).round() as u32;
    let terminal_tolerance = 3usize.max((protein_len as f64 * 0.02).ceil() as usize);
    let reaches_ends = query_interval.start <= terminal_tolerance
        && protein_len.saturating_sub(query_interval.end) <= terminal_tolerance;
    let has_ambiguous_codon = protein.as_ref().is_some_and(|value| value.contains('X'));
    let translated_internal_stop = protein.as_ref().is_some_and(|value| {
        value
            .strip_suffix('*')
            .unwrap_or(value.as_str())
            .contains('*')
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
    } else if coverage >= config.complete_coverage && reaches_ends {
        ModelState::Complete
    } else {
        ModelState::TerminalPartial
    };
    let mut qc_flags = Vec::new();
    if paf.is_none() {
        qc_flags.push("missing_paf".to_owned());
    }
    if introns
        .iter()
        .any(|intron| intron.splice_class == SpliceClass::Noncanonical)
    {
        qc_flags.push("noncanonical_splice".to_owned());
    }
    if large_query_gaps > 0 {
        qc_flags.push("large_query_gap".to_owned());
    }
    let structurally_eligible = matches!(state, ModelState::Complete | ModelState::TerminalPartial)
        && paf.is_some()
        && !qc_flags.iter().any(|flag| flag == "noncanonical_splice");
    Ok(GeneModel {
        gff_id: raw.gff_id.clone(),
        candidate_id: raw.candidate_id.clone(),
        protein_id: raw.protein_id.clone(),
        protein_len_aa: protein_len,
        query_interval_aa: query_interval,
        genomic_interval,
        strand: raw.strand,
        alignment_score: paf
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
}
