//! Observed exon candidates survive whole-gene model competition. A candidate
//! is not a claim of complete boundaries, orthology, or physical linkage.
use super::{GeneModel, SpliceClass};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug)]
pub struct MaskInterval {
    pub start: usize,
    pub end: usize,
    pub reason: &'static str,
}

fn mask_reason(base: u8) -> Option<&'static str> {
    if base == 0 {
        Some("missing_evidence")
    } else if base & 16 != 0 {
        Some("coordinate_conflict")
    } else if base.count_ones() > 1 {
        Some("base_conflict")
    } else {
        None
    }
}

fn record_mask(intervals: &mut Vec<MaskInterval>, position: usize, base: u8) {
    if let Some(reason) = mask_reason(base) {
        if let Some(last) = intervals
            .last_mut()
            .filter(|last| last.end == position && last.reason == reason)
        {
            last.end += 1;
        } else {
            intervals.push(MaskInterval {
                start: position,
                end: position + 1,
                reason,
            });
        }
    }
}

fn masks(bases: &[u8]) -> Vec<MaskInterval> {
    let mut result = Vec::new();
    for (position, &base) in bases.iter().enumerate() {
        record_mask(&mut result, position, base);
    }
    result
}

#[derive(Debug)]
pub struct ScoSequence {
    pub reference: String,
    pub sequence: String,
    pub missing_nt: usize,
    pub conflict_nt: usize,
    pub projected_models: usize,
    pub unmappable_models: usize,
    pub cross_reference_models: usize,
    pub coordinate_conflict_codons: usize,
    pub dominated_models: usize,
    pub reference_selection: &'static str,
    pub maps: BTreeMap<String, super::protein_projection::ResidueMap>,
    pub mask_intervals: Vec<MaskInterval>,
    ignored_terminal_sources: BTreeSet<(String, String, String, super::Interval)>,
    terminal_stop: [u8; 3],
}

impl ScoSequence {
    /// Preserve a terminal-stop coordinate present in the supplied nucleotide
    /// reference, using only an observed sample stop (otherwise NNN).
    pub fn append_terminal_stop(&mut self, models: &[GeneModel]) {
        let end = self.sequence.len() / 3;
        let mut bases = self.terminal_stop;
        for model in models.iter().filter(|m| {
            m.protein_id == self.reference
                && m.query_interval_aa.end == end
                && eligible(m)
                && !self.ignored_terminal_sources.contains(&(
                    m.candidate_id.clone(),
                    m.protein_id.clone(),
                    m.gff_id.clone(),
                    m.query_interval_aa,
                ))
                && m.protein.as_ref().is_some_and(|p| p.ends_with('*'))
        }) {
            for (slot, base) in bases
                .iter_mut()
                .zip(model.cds.as_bytes().iter().rev().take(3).rev())
            {
                *slot |= match base {
                    b'A' => 1,
                    b'C' => 2,
                    b'G' => 4,
                    b'T' => 8,
                    _ => 15,
                };
            }
        }
        for base in bases {
            record_mask(&mut self.mask_intervals, self.sequence.len(), base);
            self.missing_nt += usize::from(base == 0);
            self.conflict_nt += usize::from(base.count_ones() > 1);
            self.sequence.push(match base {
                1 => 'A',
                2 => 'C',
                4 => 'G',
                8 => 'T',
                _ => 'N',
            });
        }
    }
}

/// A reference-coordinate evidence composite, NOT a phased CDS or an accepted
/// gene model. Never copies reference bases. Both missing and disagreeing
/// sample evidence become N. Choosing the backbone is independent of sample.
pub fn reconstruct(
    models: &[GeneModel],
    lengths: &BTreeMap<String, usize>,
    requested: Option<&str>,
) -> Result<ScoSequence, String> {
    let references: Vec<_> = lengths
        .iter()
        .filter(|(id, _)| {
            requested.is_none_or(|name| {
                id.as_str() == name
                    || id
                        .strip_prefix(name)
                        .is_some_and(|suffix| suffix.starts_with("__auto"))
            })
        })
        .collect();
    if references.is_empty() || (requested.is_some() && references.len() != 1) {
        return Err("SCO reference ID must match exactly one protein in every family".into());
    }
    let (reference, &length) = references[0];
    let mut bases = vec![0u8; length * 3];
    let mut projected_models = 0;
    let mut unmappable_models = 0;
    for model in models
        .iter()
        .filter(|m| &m.protein_id == reference && eligible(m))
    {
        let Some(offsets) = &model.query_cds_offsets else {
            unmappable_models += 1;
            continue;
        };
        if offsets.len() != model.query_interval_aa.len()
            || model.query_interval_aa.end > length
            || offsets
                .iter()
                .flatten()
                .any(|&offset| model.cds.get(offset..offset + 3).is_none())
        {
            unmappable_models += 1;
            continue;
        }
        projected_models += 1;
        for (residue, offset) in offsets.iter().enumerate() {
            let Some(offset) = offset else { continue };
            let target = (model.query_interval_aa.start + residue) * 3;
            for (slot, base) in bases[target..target + 3]
                .iter_mut()
                .zip(model.cds.as_bytes()[*offset..*offset + 3].iter())
            {
                *slot |= match base {
                    b'A' => 1,
                    b'C' => 2,
                    b'G' => 4,
                    b'T' => 8,
                    _ => 15,
                };
            }
        }
    }

    let missing_nt = bases.iter().filter(|&&b| b == 0).count();
    let conflict_nt = bases.iter().filter(|&&b| b.count_ones() > 1).count();
    let mask_intervals = masks(&bases);
    let sequence = bases
        .into_iter()
        .map(|b| match b {
            1 => 'A',
            2 => 'C',
            4 => 'G',
            8 => 'T',
            _ => 'N',
        })
        .collect();
    Ok(ScoSequence {
        reference: reference.clone(),
        sequence,
        missing_nt,
        conflict_nt,
        projected_models,
        unmappable_models,
        cross_reference_models: 0,
        coordinate_conflict_codons: 0,
        dominated_models: 0,
        reference_selection: "legacy_explicit_reference",
        maps: BTreeMap::new(),
        mask_intervals,
        ignored_terminal_sources: BTreeSet::new(),
        terminal_stop: [0; 3],
    })
}

/// HybPiper-like best-supported target selection, using Captus-like identity
/// weighting. Each query residue contributes once, so contained/duplicate
/// models cannot inflate support. The caller pools samples to freeze one
/// backbone per SCO for the cohort.
pub fn choose_reference<'a>(
    models: impl IntoIterator<Item = &'a GeneModel>,
    proteins: &BTreeMap<String, String>,
    requested: Option<&str>,
) -> Result<(String, &'static str), String> {
    if proteins.is_empty() {
        return Err("no protein reference for SCO".into());
    }
    if let Some(name) = requested {
        let matches: Vec<_> = proteins
            .keys()
            .filter(|id| {
                id.as_str() == name
                    || id
                        .strip_prefix(name)
                        .is_some_and(|s| s.starts_with("__auto"))
            })
            .collect();
        return if matches.len() == 1 {
            Ok((matches[0].clone(), "explicit_reference"))
        } else {
            Err("SCO reference ID must match exactly one protein in every family".into())
        };
    }
    let mut lengths: Vec<_> = proteins.values().map(String::len).collect();
    lengths.sort_unstable();
    let median = lengths[lengths.len() / 2].max(1);
    let mut evidence: BTreeMap<&str, Vec<u32>> = BTreeMap::new();
    for model in models.into_iter().filter(|m| eligible(m)) {
        let Some(protein) = proteins.get(&model.protein_id) else {
            continue;
        };
        let Some(offsets) = &model.query_cds_offsets else {
            continue;
        };
        if offsets.len() != model.query_interval_aa.len()
            || model.query_interval_aa.end > protein.len()
        {
            continue;
        }
        let support = evidence
            .entry(&model.protein_id)
            .or_insert_with(|| vec![0; protein.len()]);
        for (i, offset) in offsets.iter().enumerate() {
            if offset.is_some_and(|p| model.cds.get(p..p + 3).is_some()) {
                support[model.query_interval_aa.start + i] =
                    support[model.query_interval_aa.start + i].max(model.identity_millionths);
            }
        }
    }
    let quality = |id: &String| {
        let p = &proteins[id];
        (
            std::cmp::Reverse(p.bytes().filter(|&b| b == b'X').count()),
            std::cmp::Reverse(p.len().abs_diff(median)),
            p.starts_with('M'),
        )
    };
    let score = |id: &String| {
        evidence.get(id.as_str()).map_or(0u128, |v| {
            v.iter()
                .map(|&x| u128::from(x) * u128::from(x))
                .sum::<u128>()
                * proteins[id].len().min(median) as u128
                / proteins[id].len().max(median) as u128
        })
    };
    let best = proteins
        .keys()
        .max_by(|a, b| {
            score(a)
                .cmp(&score(b))
                .then_with(|| quality(a).cmp(&quality(b)))
                .then_with(|| b.cmp(a))
        })
        .unwrap();
    Ok((
        best.clone(),
        if score(best) > 0 {
            "identity_weighted_unique_residue_support"
        } else {
            "reference_length_quality_fallback"
        },
    ))
}

fn codon_coordinates(model: &GeneModel, offset: usize) -> Option<[usize; 3]> {
    let mut result = [0; 3];
    for (base, slot) in result.iter_mut().enumerate() {
        let mut position = offset
            .checked_add(model.five_prime_trim_nt)?
            .checked_add(base)?;
        let exon = model.exons.iter().find(|e| {
            if position < e.interval.len() {
                true
            } else {
                position -= e.interval.len();
                false
            }
        })?;
        *slot = if model.strand == super::Strand::Forward {
            exon.interval.start + position
        } else {
            exon.interval.end - 1 - position
        };
    }
    Some(result)
}

/// Resolve only clearly inferior alternative placements of the same observed
/// codons, including unique exact contig containment. Other contigs, complementary
/// bases and equally strong placements still preserve disagreement.
fn dominated_placements(
    models: &[GeneModel],
    maps: &BTreeMap<String, super::protein_projection::ResidueMap>,
    equivalence: &super::contig_equivalence::Equivalence<'_>,
) -> BTreeSet<usize> {
    type Physical = BTreeMap<[usize; 3], u8>;
    let physical = |model: &GeneModel| -> Physical {
        let Some(offsets) = &model.query_cds_offsets else {
            return BTreeMap::new();
        };
        let Some(map) = maps.get(&model.protein_id) else {
            return BTreeMap::new();
        };
        offsets
            .iter()
            .enumerate()
            .filter_map(|(i, offset)| {
                let offset = (*offset)?;
                map.positions
                    .get(model.query_interval_aa.start + i)
                    .copied()
                    .flatten()?;
                let codon = model.cds.as_bytes().get(offset..offset + 3)?;
                if !codon.iter().all(|b| b"ACGT".contains(b)) {
                    return None;
                }
                let aa = *model.protein.as_ref()?.as_bytes().get(offset / 3)?;
                Some((
                    equivalence
                        .map(
                            &model.candidate_id,
                            model.strand,
                            codon_coordinates(model, offset)?,
                        )
                        .2,
                    aa,
                ))
            })
            .collect()
    };
    let mut groups: BTreeMap<(&str, super::Strand), Vec<usize>> = BTreeMap::new();
    for (i, m) in models.iter().enumerate().filter(|(_, m)| eligible(m)) {
        let (id, strand, _) = equivalence.map(&m.candidate_id, m.strand, [0; 3]);
        groups.entry((id, strand)).or_default().push(i);
    }
    let lengths: Vec<_> = models
        .iter()
        .map(|m| {
            m.query_cds_offsets
                .as_ref()
                .map_or(0, |v| v.iter().flatten().count())
        })
        .collect();
    let mut dominated = BTreeSet::new();
    for indices in groups.values() {
        let mut cache: BTreeMap<usize, Physical> = BTreeMap::new();
        for &weak in indices {
            let w = &models[weak];
            let weak_length = lengths[weak];
            if weak_length < 10 {
                continue;
            }
            for &strong in indices {
                let s = &models[strong];
                let strong_length = lengths[strong];
                if s.identity_millionths < 900_000
                    || s.identity_millionths < w.identity_millionths.saturating_add(100_000)
                    || strong_length < weak_length + 10
                {
                    continue;
                }
                cache.entry(weak).or_insert_with(|| physical(w));
                cache.entry(strong).or_insert_with(|| physical(s));
                let (a, b) = (&cache[&weak], &cache[&strong]);
                if a.len() != weak_length || b.len() < a.len() + 10 {
                    continue;
                }
                if !a.keys().all(|key| b.contains_key(key)) {
                    continue;
                }
                let context: BTreeSet<_> = b
                    .iter()
                    .filter(|(key, _)| !a.contains_key(*key))
                    .map(|(_, aa)| *aa)
                    .collect();
                // Longer low-complexity repeats cannot resolve placement.
                if context.len() < 4 {
                    continue;
                }
                dominated.insert(weak);
                break;
            }
        }
    }
    dominated
}

/// Protein alignments normally exclude the terminal stop. Recover it only
/// from the immediately adjacent, in-frame sample DNA, with an exact unique
/// ten-residue terminal anchor. Never search downstream or cross contigs.
fn adjacent_terminal_stop(
    model: &GeneModel,
    map: &super::protein_projection::ResidueMap,
    backbone: &str,
    candidates: &BTreeMap<String, String>,
) -> Option<[u8; 3]> {
    const ANCHOR: usize = 10;
    let offsets = model.query_cds_offsets.as_ref()?;
    let start = offsets.len().checked_sub(ANCHOR)?;
    let target_start = backbone.len().checked_sub(ANCHOR)?;
    let peptide = model.protein.as_ref()?;
    if model.three_prime_trim_nt != 0 || peptide.ends_with('*') {
        return None;
    }
    let last = (*offsets.last()?)?;
    if last + 3 != model.cds.len() {
        return None;
    }
    let first = last.checked_sub(3 * (ANCHOR - 1))?;
    for (k, offset) in offsets[start..].iter().enumerate() {
        if *offset != Some(first + 3 * k)
            || map.positions.get(model.query_interval_aa.start + start + k)
                != Some(&Some(target_start + k))
        {
            return None;
        }
    }
    let anchor = backbone.get(target_start..)?;
    if peptide.get(first / 3..last / 3 + 1)? != anchor
        || anchor.bytes().collect::<BTreeSet<_>>().len() < 5
        || backbone
            .as_bytes()
            .windows(ANCHOR)
            .filter(|w| *w == anchor.as_bytes())
            .count()
            != 1
    {
        return None;
    }
    let dna = candidates.get(&model.candidate_id)?;
    let coordinate = codon_coordinates(model, last)?[2];
    let triplet = if model.strand == super::Strand::Forward {
        dna.get(coordinate.checked_add(1)?..coordinate.checked_add(4)?)?
            .to_owned()
    } else {
        super::reverse_complement(dna.get(coordinate.checked_sub(3)?..coordinate)?)
    };
    if !matches!(triplet.as_str(), "TAA" | "TAG" | "TGA") {
        return None;
    }
    triplet.as_bytes().try_into().ok()
}

pub fn reconstruct_projected(
    models: &[GeneModel],
    proteins: &BTreeMap<String, String>,
    requested: Option<&str>,
) -> Result<ScoSequence, String> {
    reconstruct_with_candidates(models, proteins, requested, &BTreeMap::new())
}

/// Resolve redundant contained-contig models using observed DNA only.
pub fn reconstruct_with_candidates(
    models: &[GeneModel],
    proteins: &BTreeMap<String, String>,
    requested: Option<&str>,
    candidates: &BTreeMap<String, String>,
) -> Result<ScoSequence, String> {
    let (reference, selection) = choose_reference(models, proteins, requested)?;
    let backbone = &proteins[&reference];
    let mut maps = BTreeMap::new();
    for model in models.iter().filter(|m| eligible(m)) {
        if let Some(source) = proteins.get(&model.protein_id) {
            maps.entry(model.protein_id.clone())
                .or_insert_with(|| super::protein_projection::map(source, backbone));
        }
    }
    let equivalence = super::contig_equivalence::Equivalence::new(models, candidates);
    let dominated = dominated_placements(models, &maps, &equivalence);
    let mut bases = vec![0u8; backbone.len() * 3];
    type ObservationKey<'a> = (&'a str, super::Strand, [usize; 3]);
    type Observation = (BTreeSet<usize>, [u8; 3]);
    let mut observations: BTreeMap<ObservationKey<'_>, Observation> = BTreeMap::new();
    let mut projected_models = 0;
    let mut unmappable_models = 0;
    let mut cross_reference_models = 0;
    let mut terminal_stop = [0u8; 3];
    for (_, model) in models
        .iter()
        .enumerate()
        .filter(|(i, m)| eligible(m) && !dominated.contains(i))
    {
        let Some(map) = maps.get(&model.protein_id) else {
            unmappable_models += 1;
            continue;
        };
        let Some(offsets) = &model.query_cds_offsets else {
            unmappable_models += 1;
            continue;
        };
        if offsets.len() != model.query_interval_aa.len()
            || model.query_interval_aa.end > map.positions.len()
        {
            unmappable_models += 1;
            continue;
        }
        let mut used = false;
        for (i, offset) in offsets.iter().enumerate() {
            let Some(offset) = offset else { continue };
            let Some(target) = map.positions[model.query_interval_aa.start + i] else {
                continue;
            };
            let Some(codon) = model.cds.as_bytes().get(*offset..*offset + 3) else {
                continue;
            };
            let Some(coordinates) = codon_coordinates(model, *offset) else {
                continue;
            };
            if !codon.iter().all(|b| b"ACGT".contains(b)) {
                continue;
            }
            let entry = observations
                .entry((model.candidate_id.as_str(), model.strand, coordinates))
                .or_insert_with(|| (BTreeSet::new(), [0; 3]));
            entry.0.insert(target);
            for (slot, base) in entry.1.iter_mut().zip(codon) {
                *slot |= match base {
                    b'A' => 1,
                    b'C' => 2,
                    b'G' => 4,
                    b'T' => 8,
                    _ => 15,
                };
            }
            used = true;
        }
        if used {
            projected_models += 1;
            cross_reference_models += usize::from(model.protein_id != reference);
        } else {
            unmappable_models += 1;
        }
        if model.query_interval_aa.end == map.positions.len() {
            if let Some(stop) = adjacent_terminal_stop(model, map, backbone, candidates) {
                for (slot, base) in terminal_stop.iter_mut().zip(stop) {
                    *slot |= match base {
                        b'A' => 1,
                        b'C' => 2,
                        b'G' => 4,
                        b'T' => 8,
                        _ => 0,
                    };
                }
            }
        }
        if model.query_interval_aa.end == map.positions.len()
            && map.positions.last().copied().flatten() == backbone.len().checked_sub(1)
            && model.protein.as_ref().is_some_and(|p| p.ends_with('*'))
        {
            for (slot, base) in terminal_stop
                .iter_mut()
                .zip(model.cds.as_bytes().iter().rev().take(3).rev())
            {
                *slot |= match base {
                    b'A' => 1,
                    b'C' => 2,
                    b'G' => 4,
                    b'T' => 8,
                    _ => 15,
                };
            }
        }
    }
    let mut coordinate_conflict_codons = 0;
    for (_, (targets, codon)) in observations {
        if targets.len() > 1 {
            coordinate_conflict_codons += 1;
        }
        for target in &targets {
            for (slot, base) in bases[target * 3..target * 3 + 3].iter_mut().zip(codon) {
                *slot |= if targets.len() == 1 { base } else { 16 };
            }
        }
    }
    let missing_nt = bases.iter().filter(|&&x| x == 0).count();
    let conflict_nt = bases
        .iter()
        .filter(|&&x| x & 16 != 0 || x.count_ones() > 1)
        .count();
    let mask_intervals = masks(&bases);
    let sequence = bases
        .into_iter()
        .map(|x| match x {
            1 => 'A',
            2 => 'C',
            4 => 'G',
            8 => 'T',
            _ => 'N',
        })
        .collect();
    Ok(ScoSequence {
        reference,
        sequence,
        missing_nt,
        conflict_nt,
        projected_models,
        unmappable_models,
        cross_reference_models,
        coordinate_conflict_codons,
        dominated_models: dominated.len(),
        reference_selection: selection,
        maps,
        mask_intervals,
        ignored_terminal_sources: dominated
            .iter()
            .map(|&i| {
                let m = &models[i];
                (
                    m.candidate_id.clone(),
                    m.protein_id.clone(),
                    m.gff_id.clone(),
                    m.query_interval_aa,
                )
            })
            .collect(),
        terminal_stop,
    })
}

pub fn eligible(model: &GeneModel) -> bool {
    model.identity_millionths >= 800_000
        && model.structural_defect_count() == 0
        && model
            .protein
            .as_ref()
            .is_some_and(|p| !p.strip_suffix('*').unwrap_or(p).contains(['X', '*']))
        && !model
            .qc_flags
            .iter()
            .any(|f| matches!(f.as_str(), "padded_fragment_join" | "missing_paf"))
}

pub fn boundaries(model: &GeneModel, position: usize) -> (&'static str, &'static str) {
    let splice = |i: usize| {
        model.introns.get(i).is_some_and(|intron| {
            matches!(intron.splice_class, SpliceClass::Major | SpliceClass::Minor)
                && intron.sequence.bytes().all(|b| b"ACGT".contains(&b))
        })
    };
    let left = if position > 0 && splice(position - 1) {
        "observed_splice"
    } else {
        "unverified_alignment_end"
    };
    let right = if position + 1 < model.exons.len() && splice(position) {
        "observed_splice"
    } else {
        "unverified_alignment_end"
    };
    (left, right)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gene_annotation::{
        exonerate::parse_vulgar, reverse_complement, ModelConfig, ModelState,
    };
    fn model(
        query: &str,
        dna: &str,
        operations: &str,
        start: usize,
        end: usize,
        reverse: bool,
    ) -> GeneModel {
        let (left, right, strand) = if reverse {
            (dna.len(), 0, '-')
        } else {
            (0, dna.len(), '+')
        };
        let line = format!("vulgar: p {start} {end} . c {left} {right} {strand} 100 {operations}");
        let hit = parse_vulgar(
            &line,
            "c",
            dna.len(),
            &BTreeMap::from([("p".into(), query.len())]),
        )
        .unwrap()
        .remove(0);
        hit.normalize("c", dna, query)
            .unwrap()
            .build(dna, ModelConfig::default())
            .unwrap()
    }
    fn recover(models: &[GeneModel], length: usize) -> ScoSequence {
        reconstruct(models, &BTreeMap::from([("p".into(), length)]), Some("p")).unwrap()
    }
    #[test]
    fn terminal_stop_requires_unique_exact_consecutive_anchor() {
        let query = "MKAWDVFPGH";
        let coding = "ATGAAAGCTTGGGATGTTTTTCCTGGTCAT";
        let m = model(query, coding, "M 10 30", 0, 10, false);
        let candidates = BTreeMap::from([("c".into(), format!("{coding}TAA"))]);
        let mut map = super::super::protein_projection::map(query, query);
        assert_eq!(
            adjacent_terminal_stop(&m, &map, query, &candidates),
            Some(*b"TAA")
        );
        map.positions[5] = None;
        assert_eq!(adjacent_terminal_stop(&m, &map, query, &candidates), None);
        let repeated = query.repeat(2);
        let mut map = super::super::protein_projection::map(query, query);
        for p in &mut map.positions {
            *p = p.map(|x| x + 10);
        }
        assert_eq!(
            adjacent_terminal_stop(&m, &map, &repeated, &candidates),
            None
        );
        let map = super::super::protein_projection::map(query, query);
        let mut mismatch = m.clone();
        mismatch.protein.as_mut().unwrap().replace_range(5..6, "A");
        assert_eq!(
            adjacent_terminal_stop(&mismatch, &map, query, &candidates),
            None
        );
        let mut gap = m.clone();
        gap.query_cds_offsets.as_mut().unwrap()[5] = None;
        assert_eq!(adjacent_terminal_stop(&gap, &map, query, &candidates), None);
        assert_eq!(
            adjacent_terminal_stop(&m, &map, query, &BTreeMap::new()),
            None
        );
    }

    #[test]
    fn terminal_stop_uses_adjacent_sample_dna_on_both_strands_and_preserves_conflicts() {
        let query = "MKAWDVFPGH";
        let coding = "ATGAAAGCTTGGGATGTTTTTCCTGGTCAT";
        let proteins = BTreeMap::from([("p".into(), query.into())]);
        for reverse in [false, true] {
            let build = |stop: &str, id: &str| {
                let full = format!("{coding}{stop}");
                let dna = if reverse {
                    reverse_complement(&full)
                } else {
                    full
                };
                let bounds = if reverse {
                    format!("{} {} -", dna.len(), stop.len())
                } else {
                    format!("0 {} +", coding.len())
                };
                let hit = parse_vulgar(
                    &format!("vulgar: p 0 10 . {id} {bounds} 100 M 10 30"),
                    id,
                    dna.len(),
                    &BTreeMap::from([("p".into(), 10)]),
                )
                .unwrap()
                .remove(0);
                let m = hit
                    .normalize(id, &dna, query)
                    .unwrap()
                    .build(&dna, ModelConfig::default())
                    .unwrap();
                (m, dna)
            };
            for stop in ["TAA", "TAG", "TGA", "AAA", "NAA", "", "ATGA"] {
                let (m, dna) = build(stop, "c");
                let mut result = reconstruct_with_candidates(
                    std::slice::from_ref(&m),
                    &proteins,
                    Some("p"),
                    &BTreeMap::from([("c".into(), dna)]),
                )
                .unwrap();
                // Stop coordinates are added only when requested by the nucleotide reference.
                assert_eq!(result.sequence, coding);
                result.append_terminal_stop(&[m]);
                let expected = if matches!(stop, "TAA" | "TAG" | "TGA") {
                    stop
                } else {
                    "NNN"
                };
                assert_eq!(result.sequence, format!("{coding}{expected}"));
            }
            let (a, da) = build("TAA", "a");
            let (b, db) = build("TAG", "b");
            let models = [a, b];
            let mut result = reconstruct_with_candidates(
                &models,
                &proteins,
                Some("p"),
                &BTreeMap::from([("a".into(), da), ("b".into(), db)]),
            )
            .unwrap();
            result.append_terminal_stop(&models);
            assert_eq!(result.sequence, format!("{coding}TAN"));
            assert_eq!(result.conflict_nt, 1);
        }
    }

    #[test]
    fn longer_informative_same_contig_placement_resolves_only_clearly_dominated_hits() {
        let query = "MKAWDVFPGHYSNTRECILQMKAWDVFPGH";
        let dna = "ATGAAAGCTTGGGATGTTTTTCCTGGTCATTATTCTAATACTCGTGAATGTATTCTGCAA";
        let strong = model(query, dna, "M 20 60", 0, 20, false);
        let mut weak = model(query, &dna[..30], "M 10 30", 20, 30, false);
        weak.identity_millionths = 850_000;
        let proteins = BTreeMap::from([("p".into(), query.into())]);
        let recovered =
            reconstruct_projected(&[strong.clone(), weak.clone()], &proteins, Some("p")).unwrap();
        assert_eq!(recovered.dominated_models, 1);
        assert_eq!(recovered.sequence, format!("{dna}{}", "N".repeat(30)));
        assert_eq!(recovered.coordinate_conflict_codons, 0);
        weak.identity_millionths = 950_000;
        let unresolved =
            reconstruct_projected(&[strong.clone(), weak], &proteins, Some("p")).unwrap();
        assert_eq!(unresolved.dominated_models, 0);
        assert!(unresolved.coordinate_conflict_codons > 0);
        assert!(unresolved
            .mask_intervals
            .iter()
            .any(|m| m.reason == "coordinate_conflict"));
        let mut allele_dna = dna[..30].to_owned();
        allele_dna.replace_range(8..9, "C");
        let mut allele = model(query, &allele_dna, "M 10 30", 0, 10, false);
        allele.candidate_id = "other_contig".into();
        allele.identity_millionths = 850_000;
        let conflict = reconstruct_projected(&[strong, allele], &proteins, Some("p")).unwrap();
        assert_eq!(conflict.dominated_models, 0);
        assert_eq!(&conflict.sequence[6..9], "GCN");
        assert!(conflict
            .mask_intervals
            .iter()
            .any(|m| m.reason == "base_conflict"));
    }
    #[test]
    fn exact_containment_resolves_weak_cross_contig_placement_but_not_variants_or_repeats() {
        let query = "MKAWDVFPGHYSNTRECILQMKAWDVFPGH";
        let dna = "ATGAAAGCTTGGGATGTTTTTCCTGGTCATTATTCTAATACTCGTGAATGTATTCTGCAA";
        let strong = model(query, dna, "M 20 60", 0, 20, false);
        let mut weak = model(query, &dna[..30], "M 10 30", 20, 30, false);
        weak.candidate_id = "short".into();
        weak.identity_millionths = 850_000;
        let proteins = BTreeMap::from([("p".into(), query.into())]);
        let candidates =
            BTreeMap::from([("c".into(), dna.into()), ("short".into(), dna[..30].into())]);
        let r = reconstruct_with_candidates(
            &[strong.clone(), weak.clone()],
            &proteins,
            Some("p"),
            &candidates,
        )
        .unwrap();
        assert_eq!(r.dominated_models, 1);
        assert_eq!(r.sequence, format!("{dna}{}", "N".repeat(30)));
        // Reverse-complement containment must map each physical base correctly.
        let reverse = reverse_complement(&dna[..30]);
        let mut reverse_weak = model(query, &reverse, "M 10 30", 20, 30, true);
        reverse_weak.candidate_id = "short".into();
        reverse_weak.identity_millionths = 850_000;
        let mut reversed_candidates = candidates.clone();
        reversed_candidates.insert("short".into(), reverse);
        let r = reconstruct_with_candidates(
            &[strong.clone(), reverse_weak],
            &proteins,
            Some("p"),
            &reversed_candidates,
        )
        .unwrap();
        assert_eq!(r.dominated_models, 1);
        assert_eq!(r.sequence, format!("{dna}{}", "N".repeat(30)));
        // Same-length models cannot use containment to break a biological tie.
        weak.identity_millionths = 950_000;
        assert_eq!(
            reconstruct_with_candidates(
                &[strong.clone(), weak.clone()],
                &proteins,
                Some("p"),
                &candidates
            )
            .unwrap()
            .dominated_models,
            0
        );
        weak.identity_millionths = 850_000;
        let mut variant_candidates = candidates.clone();
        variant_candidates
            .get_mut("short")
            .unwrap()
            .replace_range(8..9, "C");
        assert_eq!(
            reconstruct_with_candidates(
                &[strong.clone(), weak.clone()],
                &proteins,
                Some("p"),
                &variant_candidates
            )
            .unwrap()
            .dominated_models,
            0
        );
        let mut another_parent = strong.clone();
        another_parent.candidate_id = "repeat".into();
        let mut repeated = candidates;
        repeated.insert("repeat".into(), dna.into());
        assert_eq!(
            reconstruct_with_candidates(
                &[strong, weak, another_parent],
                &proteins,
                Some("p"),
                &repeated
            )
            .unwrap()
            .dominated_models,
            0
        );
    }

    #[test]
    fn longer_low_complexity_placement_does_not_resolve_repeat_ambiguity() {
        let query = "A".repeat(40);
        let strong = model(&query, &"GCT".repeat(20), "M 20 60", 0, 20, false);
        let mut weak = model(&query, &"GCT".repeat(10), "M 10 30", 20, 30, false);
        weak.identity_millionths = 850_000;
        let proteins = BTreeMap::from([("p".into(), query)]);
        let recovered = reconstruct_projected(&[strong, weak], &proteins, Some("p")).unwrap();
        assert_eq!(recovered.dominated_models, 0);
        assert!(recovered.coordinate_conflict_codons > 0);
        assert_eq!(
            recovered
                .mask_intervals
                .iter()
                .map(|r| r.end - r.start)
                .sum::<usize>(),
            recovered.missing_nt + recovered.conflict_nt
        );
    }

    #[test]
    fn supported_backbone_selection_does_not_count_duplicate_hits() {
        let mut a = model("AAAA", &"GCT".repeat(3), "M 3 9", 0, 3, false);
        a.protein_id = "z".into();
        let mut b = model("AAAA", &"GCT".repeat(2), "M 2 6", 0, 2, false);
        b.protein_id = "a".into();
        let proteins = BTreeMap::from([("a".into(), "AAAA".into()), ("z".into(), "AAAA".into())]);
        let mut models = vec![b; 30];
        models.push(a);
        assert_eq!(choose_reference(&models, &proteins, None).unwrap().0, "z");
        assert_eq!(
            choose_reference(&models, &proteins, Some("a")).unwrap().0,
            "a"
        );
        assert!(choose_reference(&models, &proteins, Some("missing")).is_err());
    }

    #[test]
    fn other_reference_evidence_recovers_observed_bases_and_masks_alleles() {
        let mut a = model("AAAA", &"GCT".repeat(2), "M 2 6", 0, 2, false);
        a.protein_id = "other".into();
        let proteins = BTreeMap::from([
            ("fixed".into(), "AAAA".into()),
            ("other".into(), "AAAA".into()),
        ]);
        let r = reconstruct_projected(&[a.clone()], &proteins, Some("fixed")).unwrap();
        assert_eq!(r.sequence, "GCTGCTNNNNNN");
        assert_eq!(r.cross_reference_models, 1);
        assert_eq!(
            reconstruct_projected(&[a.clone(), a.clone()], &proteins, Some("fixed"))
                .unwrap()
                .sequence,
            r.sequence
        );
        let mut b = model("AAAA", "GCC", "M 1 3", 1, 2, false);
        b.protein_id = "fixed".into();
        b.candidate_id = "different_contig".into();
        assert_eq!(
            reconstruct_projected(&[a, b], &proteins, Some("fixed"))
                .unwrap()
                .sequence,
            "GCTGCNNNNNNN"
        );
    }

    #[test]
    fn same_sample_codon_cannot_be_reused_at_conflicting_reference_positions() {
        let mut a = model("AAAA", "GCT", "M 1 3", 0, 1, false);
        a.protein_id = "other".into();
        let mut b = a.clone();
        b.protein_id = "fixed".into();
        b.query_interval_aa = crate::gene_annotation::Interval { start: 1, end: 2 };
        let proteins = BTreeMap::from([
            ("fixed".into(), "AAAA".into()),
            ("other".into(), "AAAA".into()),
        ]);
        let r = reconstruct_projected(&[a, b], &proteins, Some("fixed")).unwrap();
        assert_eq!(r.sequence, "N".repeat(12));
        assert_eq!(r.coordinate_conflict_codons, 1);
        assert_eq!(r.conflict_nt, 6);
        assert_eq!(r.missing_nt, 6);
    }

    #[test]
    fn cross_reference_indel_preserves_downstream_sample_codon_coordinates() {
        let backbone = "MKWVTFISLLFLFSSAYSR";
        let source = "MKWVTFGQISLLFLFSSAYSR";
        // A terminal fragment with three observed codons, beyond the insertion.
        let mut a = model(source, "TATAGTCGT", "M 3 9", 18, 21, false);
        a.protein_id = "source".into();
        let proteins = BTreeMap::from([
            ("fixed".into(), backbone.into()),
            ("source".into(), source.into()),
        ]);
        let r = reconstruct_projected(&[a], &proteins, Some("fixed")).unwrap();
        assert_eq!(&r.sequence[48..], "TATAGTCGT");
        assert_eq!(r.sequence[..48], "N".repeat(48));
    }

    #[test]
    fn projected_split_codons_on_both_strands_match_direct_evidence() {
        let dna = "GCTGGTAAAGCTGCT";
        for reverse in [false, true] {
            let sequence = if reverse {
                reverse_complement(dna)
            } else {
                dna.into()
            };
            let m = model(
                "AAA",
                &sequence,
                "M 1 3 S 0 1 5 0 2 I 0 2 3 0 2 S 1 2 M 1 3",
                0,
                3,
                reverse,
            );
            let proteins = BTreeMap::from([("p".into(), "AAA".into())]);
            assert_eq!(
                reconstruct_projected(&[m.clone()], &proteins, Some("p"))
                    .unwrap()
                    .sequence,
                recover(&[m], 3).sequence
            );
        }
    }
    #[test]
    fn one_sco_keeps_order_masks_conflicts_and_never_fills_missing_reference() {
        let a = model("AAAAAAAAAA", &"GCT".repeat(4), "M 4 12", 0, 4, false);
        let mut b = model("AAAAAAAAAA", "GCCGCTGCT", "M 3 9", 3, 6, false);
        b.state = ModelState::AmbiguousModel;
        b.selected = false;
        let result = recover(&[a.clone(), b.clone()], 10);
        assert_eq!(
            result.sequence,
            format!(
                "{}GCN{}{}",
                "GCT".repeat(3),
                "GCT".repeat(2),
                "N".repeat(12)
            )
        );
        assert_eq!((result.missing_nt, result.conflict_nt), (12, 1));
        assert_eq!(recover(&[b, a], 10).sequence, result.sequence);
        assert_eq!(recover(&[], 10).sequence, "N".repeat(30));
    }
    #[test]
    fn insertion_and_deletion_do_not_shift_downstream_exons() {
        let m = model(
            "AAAAA",
            "GCTGCTAAAGCTGCT",
            "M 1 3 G 1 0 M 1 3 G 0 3 M 2 6",
            0,
            5,
            false,
        );
        assert_eq!(
            m.query_cds_offsets,
            Some(vec![Some(0), None, Some(3), Some(9), Some(12)])
        );
        let result = recover(&[m], 5);
        assert_eq!(result.sequence, "GCTNNNGCTGCTGCT");
        assert_eq!(result.unmappable_models, 0);
    }
    #[test]
    fn split_codons_and_reverse_strand_preserve_reference_order() {
        let dna = "GCTGGTAAAGCTGCT";
        for reverse in [false, true] {
            let sequence = if reverse {
                reverse_complement(dna)
            } else {
                dna.into()
            };
            let m = model(
                "AAA",
                &sequence,
                "M 1 3 S 0 1 5 0 2 I 0 2 3 0 2 S 1 2 M 1 3",
                0,
                3,
                reverse,
            );
            assert_eq!(
                boundaries(&m, 0),
                ("unverified_alignment_end", "observed_splice")
            );
            assert_eq!(
                boundaries(&m, 1),
                ("observed_splice", "unverified_alignment_end")
            );
            assert_eq!(recover(&[m], 3).sequence, "GCTGCTGCT");
        }
    }
    #[test]
    fn requested_reference_must_exist_and_padding_is_not_observation() {
        assert!(reconstruct(&[], &BTreeMap::from([("p".into(), 4)]), Some("missing")).is_err());
        let mut m = model("AAAA", &"GCT".repeat(4), "M 4 12", 0, 4, false);
        m.qc_flags.push("padded_fragment_join".into());
        assert_eq!(recover(&[m], 4).sequence, "N".repeat(12));
    }

    #[test]
    fn reference_stop_is_not_copied_without_sample_evidence() {
        let mut result = recover(&[], 4);
        result.append_terminal_stop(&[]);
        assert_eq!(result.sequence, "N".repeat(15));
        assert_eq!(result.missing_nt, 15);
        let mut m = model("AAAA", &"GCT".repeat(4), "M 4 12", 0, 4, false);
        m.cds.push_str("TAA");
        m.protein.as_mut().unwrap().push('*');
        let mut result = recover(&[m.clone()], 4);
        result.append_terminal_stop(&[m]);
        assert_eq!(result.sequence, format!("{}TAA", "GCT".repeat(4)));
        assert_eq!(result.missing_nt, 0);
    }
}
