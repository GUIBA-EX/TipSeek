//! Protein-to-genome vulgar adapter. Coordinates are zero-based and half-open.
use super::{Interval, Strand};
use std::collections::BTreeMap;

#[derive(Debug, PartialEq)]
pub struct Hit {
    pub protein: String,
    pub query: Interval,
    pub genomic: Interval,
    pub strand: Strand,
    pub score: i64,
    pub cdss: Vec<Interval>,
    pub frameshift_ops: usize,
    operations: Vec<(String, usize, usize)>,
    target_start: usize,
}

// Vulgar endpoints are oriented, zero-based interbase coordinates. Splice
// signals belong to the intron, including split-codon phases 1 and 2.
pub fn parse_vulgar(
    text: &str,
    candidate: &str,
    length: usize,
    proteins: &BTreeMap<String, usize>,
) -> Result<Vec<Hit>, String> {
    let mut hits = Vec::new();
    for line in text.lines().filter(|line| !line.trim().is_empty()) {
        let f: Vec<_> = line.split_whitespace().collect();
        if f.len() < 13 || f[0] != "vulgar:" || (f.len() - 10) % 3 != 0 {
            return Err("invalid_vulgar_record".into());
        }
        let number = |i: usize| {
            f[i].parse::<usize>()
                .map_err(|_| "invalid_coordinate".to_owned())
        };
        let (qs, qe, ts, te) = (number(2)?, number(3)?, number(6)?, number(7)?);
        let query_len = proteins.get(f[1]).ok_or("unknown_protein")?;
        if f[5] != candidate
            || !matches!(f[4], "." | "+")
            || qs >= qe
            || qe > *query_len
            || ts > length
            || te > length
            || ts == te
        {
            return Err("invalid_alignment_bounds".into());
        }
        let strand = match f[8] {
            "+" if ts < te => Strand::Forward,
            "-" if ts > te => Strand::Reverse,
            _ => return Err("invalid_strand".into()),
        };
        let mut hit = Hit {
            protein: f[1].into(),
            query: Interval::new(qs, qe).unwrap(),
            genomic: Interval::new(ts.min(te), ts.max(te)).unwrap(),
            strand,
            score: f[9].parse().map_err(|_| "invalid_score")?,
            cdss: Vec::new(),
            frameshift_ops: 0,
            operations: Vec::new(),
            target_start: ts,
        };
        let (mut q, mut t) = (qs, ts);
        let mut exon_start = t;
        for op in f[10..].chunks_exact(3) {
            let qa: usize = op[1].parse().map_err(|_| "invalid_advance")?;
            let ta: usize = op[2].parse().map_err(|_| "invalid_advance")?;
            if !matches!(op[0], "M" | "G" | "F" | "S" | "5" | "3" | "I")
                || (matches!(op[0], "5" | "3" | "I") && qa != 0)
            {
                return Err("unsupported_operation".into());
            }
            hit.operations.push((op[0].into(), qa, ta));
            let next = if strand == Strand::Forward {
                t.checked_add(ta)
            } else {
                t.checked_sub(ta)
            }
            .ok_or("coordinate_overflow")?;
            q = q.checked_add(qa).ok_or("coordinate_overflow")?;
            if q > qe || next > length {
                return Err("operation_out_of_bounds".into());
            }
            if matches!(op[0], "5" | "3" | "I") {
                if exon_start != t {
                    hit.cdss
                        .push(Interval::new(exon_start.min(t), exon_start.max(t)).unwrap());
                }
                exon_start = next;
            }
            if op[0] == "F" {
                hit.frameshift_ops += 1;
            }
            t = next;
        }
        if q != qe || t != te {
            return Err("operation_span_mismatch".into());
        }
        if exon_start != t {
            hit.cdss
                .push(Interval::new(exon_start.min(t), exon_start.max(t)).unwrap());
        }
        if hit.cdss.is_empty() {
            return Err("missing_coding_intervals".into());
        }
        hit.cdss.sort_by_key(|x| x.start);
        hits.push(hit);
    }
    Ok(hits)
}

use super::evidence::{AlignmentEngine, AlignmentEvidence};
use super::{
    build_gene_model_with_evidence, reverse_complement, translate_cds, CdsFeature, GeneModel,
    ModelConfig, RawGeneModel,
};

/// Normalized evidence retains engine provenance; scores must not be compared
/// between engines. Identity counts aligned M/S amino-acid pairs (gaps excluded).
pub struct ExonerateModel {
    pub raw: RawGeneModel,
    pub evidence: AlignmentEvidence,
    pub aligned_amino_acids: usize,
    pub identical_amino_acids: usize,
    reference_start: Option<u8>,
}

impl ExonerateModel {
    pub fn build(&self, target: &str, config: ModelConfig) -> Result<GeneModel, String> {
        let mut model =
            build_gene_model_with_evidence(&self.raw, Some(&self.evidence), target, None, config)?;
        model.check_reference_start(self.reference_start);
        Ok(model)
    }
}

impl Hit {
    pub fn normalize(
        &self,
        candidate: &str,
        target: &str,
        protein: &str,
    ) -> Result<ExonerateModel, String> {
        if !target.is_ascii() || !protein.is_ascii() || self.query.end > protein.len() {
            return Err("invalid_input_sequence".into());
        }
        let mut q = self.query.start;
        let mut t = self.target_start;
        let mut split = String::new();
        let mut aligned = 0usize;
        let mut identical = 0usize;
        let mut query_gaps = 0usize;
        let mut query_gap_aa = 0usize;
        let mut inframe_stops = 0usize;
        let mut intron = false;
        for (op, qa, ta) in &self.operations {
            let next = if self.strand == Strand::Forward {
                t.checked_add(*ta)
            } else {
                t.checked_sub(*ta)
            }
            .ok_or("coordinate_overflow")?;
            let bases = target
                .get(t.min(next)..t.max(next))
                .ok_or("target_out_of_bounds")?;
            let bases = if !matches!(op.as_str(), "M" | "S" | "G") {
                String::new()
            } else if self.strand == Strand::Reverse {
                reverse_complement(bases)
            } else {
                bases.to_ascii_uppercase()
            };
            let mut paired = None;
            match op.as_str() {
                "M" if !intron && split.is_empty() && *qa > 0 && qa.checked_mul(3) == Some(*ta) => {
                    paired = Some(translate_cds(&bases).map_err(|e| e.to_string())?);
                }
                "G" if !intron
                    && split.is_empty()
                    && ((*qa > 0 && *ta == 0) || (*qa == 0 && *ta > 0 && ta.is_multiple_of(3))) =>
                {
                    query_gaps += usize::from(*qa >= 10);
                    query_gap_aa += qa;
                    if *ta > 0 {
                        inframe_stops += translate_cds(&bases)
                            .map_err(|e| e.to_string())?
                            .bytes()
                            .filter(|b| *b == b'*')
                            .count();
                    }
                }
                "F" if !intron && split.is_empty() && *qa == 0 && *ta > 0 => {}
                "S" if !intron && *qa == 0 && matches!(*ta, 1 | 2) && split.is_empty() => {
                    split.push_str(&bases);
                }
                "S" if !intron && *qa == 1 && !split.is_empty() && split.len() + ta == 3 => {
                    split.push_str(&bases);
                    paired = Some(translate_cds(&split).map_err(|e| e.to_string())?);
                    split.clear();
                }
                "5" if !intron && *qa == 0 && *ta == 2 => {
                    intron = true;
                }
                "I" if intron && *qa == 0 && *ta > 0 => {}
                "3" if intron && *qa == 0 && *ta == 2 => {
                    intron = false;
                }
                _ => return Err("invalid_protein2genome_operation".into()),
            }
            if let Some(translated) = paired {
                let query = protein.get(q..q + qa).ok_or("query_out_of_bounds")?;
                for (a, b) in query.bytes().zip(translated.bytes()) {
                    aligned += 1;
                    identical += usize::from(a.to_ascii_uppercase() == b && b != b'X');
                    // An aligned stop replacing an amino acid is a defect even
                    // at the path endpoint; vulgar has no miniprot st tag.
                    inframe_stops += usize::from(b == b'*' && a != b'*');
                }
            }
            q += qa;
            t = next;
        }
        if intron || !split.is_empty() || aligned == 0 {
            return Err("incomplete_alignment_evidence".into());
        }
        let mut intervals = self.cdss.clone();
        if self.strand == Strand::Reverse {
            intervals.reverse();
        }
        let mut coding_length = 0;
        let cdss = intervals
            .into_iter()
            .map(|interval| {
                let phase = ((3 - coding_length % 3) % 3) as u8;
                coding_length += interval.len();
                CdsFeature {
                    interval,
                    strand: self.strand,
                    phase: Some(phase),
                }
            })
            .collect();
        let mut offsets = Vec::new();
        let mut coding_nt = 0;
        for (op, qa, ta) in &self.operations {
            match op.as_str() {
                "M" => {
                    offsets.extend((0..*qa).map(|i| Some(coding_nt + 3 * i)));
                    coding_nt += ta;
                }
                "G" => {
                    offsets.extend(std::iter::repeat_n(None, *qa));
                    coding_nt += ta;
                }
                "S" => {
                    coding_nt += ta;
                    if *qa == 1 {
                        offsets.push(coding_nt.checked_sub(3));
                    }
                }
                "F" => coding_nt += ta,
                _ => {}
            }
        }
        let coding_map = (self.frameshift_ops == 0 && offsets.len() == self.query.len())
            .then_some((offsets, coding_nt));
        Ok(ExonerateModel {
            raw: RawGeneModel {
                gff_id: format!(
                    "EXO_{}_{}_{}_{}",
                    self.protein,
                    self.genomic.start,
                    self.genomic.end,
                    self.strand.as_char()
                ),
                candidate_id: candidate.into(),
                protein_id: self.protein.clone(),
                target_interval_aa: self.query,
                strand: self.strand,
                genomic_interval: self.genomic,
                gff_score: Some(self.score),
                identity_millionths: Some(
                    (identical as f64 / aligned as f64 * 1_000_000.0).round() as u32,
                ),
                // Vulgar does not report substitution-matrix positives.
                positive_millionths: None,
                frameshifts: self.frameshift_ops,
                inframe_stops,
                rank: None,
                cdss,
                paf: None,
                evidence: None,
            },
            evidence: AlignmentEvidence {
                engine: AlignmentEngine::Exonerate,
                query_len: protein.len(),
                query_interval: self.query,
                alignment_score: Some(self.score),
                frameshifts: self.frameshift_ops,
                inframe_stops,
                large_query_gaps: query_gaps,
                query_gap_aa,
                coding_map,
            },
            aligned_amino_acids: aligned,
            identical_amino_acids: identical,
            reference_start: protein.bytes().next(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gene_annotation::ModelState;

    fn adapt(line: &str, target: &str, protein: &str) -> Result<ExonerateModel, String> {
        parse_vulgar(
            line,
            "c",
            target.len(),
            &BTreeMap::from([("p".into(), protein.len())]),
        )?
        .remove(0)
        .normalize("c", target, protein)
    }

    #[test]
    fn split_codons_on_both_strands_preserve_phase_and_translation() {
        for phase in [1usize, 2] {
            let cds = "ATGAAATTT";
            let target = format!("{}GTCCCCAG{}", &cds[..3 + phase], &cds[3 + phase..]);
            for reverse in [false, true] {
                let sequence = if reverse {
                    reverse_complement(&target)
                } else {
                    target.clone()
                };
                let bounds = if reverse { "17 0 -" } else { "0 17 +" };
                let line = format!("vulgar: p 0 3 . c {bounds} 30 M 1 3 S 0 {phase} 5 0 2 I 0 4 3 0 2 S 1 {} M 1 3", 3 - phase);
                let adapted = adapt(&line, &sequence, "MKF").unwrap();
                assert!(adapted.raw.paf.is_none());
                assert_eq!(adapted.evidence.engine, AlignmentEngine::Exonerate);
                assert_eq!(adapted.identical_amino_acids, 3);
                let model = adapted.build(&sequence, ModelConfig::default()).unwrap();
                assert_eq!(model.cds, cds);
                assert_eq!(model.protein.as_deref(), Some("MKF"));
                assert_eq!(model.exons[1].phase, Some((3 - phase) as u8));
                assert_eq!(model.introns[0].sequence, "GTCCCCAG");
                assert!(model.eligible_for_resolve);
                assert!(model.qc_flags.is_empty());
            }
        }
    }

    #[test]
    fn frameshifts_internal_stops_and_noncanonical_splices_fail_common_qc() {
        for (line, target, protein, state) in [
            (
                "vulgar: p 0 2 . c 0 7 + 20 M 1 3 F 0 1 M 1 3",
                "ATGCAAA",
                "MK",
                ModelState::Frameshifted,
            ),
            (
                "vulgar: p 0 3 . c 0 9 + 20 M 3 9",
                "ATGTAAAAA",
                "MKK",
                ModelState::InternalStop,
            ),
        ] {
            let model = adapt(line, target, protein)
                .unwrap()
                .build(target, ModelConfig::default())
                .unwrap();
            assert_eq!(model.state, state);
            assert!(!model.eligible_for_resolve);
        }
        let target = "ATGCCCCAAA";
        let model = adapt(
            "vulgar: p 0 2 . c 0 10 + 20 M 1 3 5 0 2 3 0 2 M 1 3",
            target,
            "MK",
        )
        .unwrap()
        .build(target, ModelConfig::default())
        .unwrap();
        assert!(model.qc_flags.iter().any(|f| f == "noncanonical_splice"));
        assert!(!model.eligible_for_resolve);
    }

    #[test]
    fn stop_aligned_to_residue_or_inserted_at_path_end_is_not_a_terminal_stop() {
        for (line, target, protein) in [
            ("vulgar: p 0 3 . c 0 9 + 20 M 3 9", "ATGAAATAA", "MKF"),
            ("vulgar: p 0 3 . c 9 0 - 20 M 3 9", "TTATTTCAT", "MKF"),
            ("vulgar: p 0 2 . c 0 9 + 20 M 2 6 G 0 3", "ATGAAATAA", "MK"),
            (
                "vulgar: p 0 3 . c 0 17 + 20 M 2 6 S 0 1 5 0 2 I 0 4 3 0 2 S 1 2",
                "ATGAAATGTCCCCAGAA",
                "MKF",
            ),
        ] {
            let model = adapt(line, target, protein)
                .unwrap()
                .build(target, ModelConfig::default())
                .unwrap();
            assert_eq!(model.state, ModelState::InternalStop);
            assert_eq!(model.inframe_stops, 1);
            assert!(!model.eligible_for_resolve);
        }
    }

    #[test]
    fn reaching_reference_ends_with_a_coding_gap_is_not_complete() {
        let target = "ATGAAATTTCCC";
        let model = adapt(
            "vulgar: p 0 5 . c 0 12 + 20 M 2 6 G 1 0 M 2 6",
            target,
            "MKAFP",
        )
        .unwrap()
        .build(target, ModelConfig::default())
        .unwrap();
        assert_eq!(model.state, ModelState::InternalPartial);
        assert_eq!(model.coverage_millionths, 800_000);
        assert!(model.qc_flags.iter().any(|f| f == "query_gap_aa=1"));
    }

    #[test]
    fn gaps_and_terminal_stops_keep_explicit_semantics() {
        let target = "ATGAAA";
        let adapted = adapt(
            "vulgar: p 0 12 . c 0 6 + 20 M 1 3 G 10 0 M 1 3",
            target,
            "MAAAAAAAAAAK",
        )
        .unwrap();
        assert_eq!(adapted.evidence.large_query_gaps, 1);
        assert_eq!(adapted.raw.identity_millionths, Some(1_000_000));
        assert_eq!(adapted.raw.positive_millionths, None);
        let target = "ATGCCCAAA";
        let model = adapt("vulgar: p 0 2 . c 0 9 + 20 M 1 3 G 0 3 M 1 3", target, "MK")
            .unwrap()
            .build(target, ModelConfig::default())
            .unwrap();
        assert_eq!(model.protein.as_deref(), Some("MPK"));
        for (protein, endpoint, expected) in [("MK", 6, "MK"), ("MK*", 9, "MK*")] {
            let target = "ATGAAATAA";
            let line = format!(
                "vulgar: p 0 {} . c 0 {endpoint} + 20 M {} {endpoint}",
                protein.len(),
                protein.len()
            );
            let model = adapt(&line, target, protein)
                .unwrap()
                .build(target, ModelConfig::default())
                .unwrap();
            assert_eq!(model.protein.as_deref(), Some(expected));
            assert!(model.eligible_for_resolve);
        }
    }

    #[test]
    fn rejects_span_consistent_but_invalid_operations() {
        for (line, target) in [
            ("vulgar: p 0 2 . c 0 6 + 20 M 2 5 F 0 1", "ATGAAA"),
            ("vulgar: p 0 2 . c 0 6 + 20 S 1 3 M 1 3", "ATGAAA"),
            ("vulgar: p 0 2 . c 0 6 + 20 G 1 3 M 1 3", "ATGAAA"),
        ] {
            assert!(adapt(line, target, "MK").is_err());
        }
    }
}

#[cfg(test)]
mod repeated_hit_tests {
    use super::*;

    #[test]
    fn repeated_loci_and_query_offsets_remain_separate() {
        let target = "ATGAAAAAAATGAAA";
        let protein = "AAAAAAAAAAMKAAAAAAAAAA";
        let text = "vulgar: p 10 12 . c 0 6 + 20 M 2 6\nvulgar: p 10 12 . c 9 15 + 20 M 2 6";
        let hits = parse_vulgar(
            text,
            "c",
            target.len(),
            &BTreeMap::from([("p".into(), protein.len())]),
        )
        .unwrap();
        assert_eq!(hits.len(), 2);
        for (i, hit) in hits.iter().enumerate() {
            let adapted = hit.normalize("c", target, protein).unwrap();
            assert_eq!(adapted.identical_amino_acids, 2);
            let model = adapted
                .build(
                    target,
                    ModelConfig {
                        minimum_coverage: 0.05,
                        complete_coverage: 0.8,
                    },
                )
                .unwrap();
            assert_eq!(model.genomic_interval.start, i * 9);
            assert_eq!(model.query_interval_aa.start, 10);
            assert_eq!(model.state, super::super::ModelState::TerminalPartial);
            assert!(model.eligible_for_resolve);
        }
        assert_ne!(
            hits[0].normalize("c", target, protein).unwrap().raw.gff_id,
            hits[1].normalize("c", target, protein).unwrap().raw.gff_id
        );
    }
}
