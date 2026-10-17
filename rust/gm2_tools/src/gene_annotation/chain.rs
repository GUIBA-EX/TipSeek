//! Reference-ordered fragment recovery. A chain is a hypothesis, not linkage evidence.
use super::{GeneModel, ModelState};
use std::collections::{BTreeMap, BTreeSet};

pub const MIN_FRAGMENT_AA: usize = 10;
const MIN_IDENTITY: u32 = 800_000;
const MAX_NODES: usize = 128;
const MAX_PATHS: usize = 64;

#[derive(Debug, PartialEq, Eq)]
pub struct FragmentChain {
    pub protein: String,
    pub models: Vec<usize>,
    pub covered_aa: usize,
    pub ambiguous: bool,
    pub search_limited: bool,
}

pub fn source_eligible(m: &GeneModel) -> bool {
    matches!(
        m.state,
        ModelState::TerminalPartial | ModelState::LowCoverage | ModelState::AmbiguousModel
    ) && m.query_interval_aa.len() >= MIN_FRAGMENT_AA
        && m.identity_millionths >= MIN_IDENTITY
        && m.structural_defect_count() == 0
        && m.large_query_gaps == 0
        && !m.qc_flags.iter().any(|f| f.starts_with("query_gap_aa="))
        && m.protein
            .as_ref()
            .is_some_and(|p| !p.strip_suffix('*').unwrap_or(p).contains(['X', '*']))
        && !m
            .qc_flags
            .iter()
            .any(|f| matches!(f.as_str(), "missing_paf" | "reference_start_mismatch"))
}

fn coding(m: &GeneModel) -> Option<&str> {
    let sequence = if m.protein.as_ref()?.ends_with('*') {
        m.cds.get(..m.cds.len().checked_sub(3)?)?
    } else {
        &m.cds
    };
    // Without a residue-to-CDS map, do not trim indel-bearing paths by multiplying
    // a query offset by three. Genomic reannotation may still validate these.
    (sequence.len() == m.query_interval_aa.len() * 3).then_some(sequence)
}

fn compatible(a: &GeneModel, b: &GeneModel) -> bool {
    connection_status(a, b) == "compatible"
}

/// Explain ordered graph edges without hiding conflicting coding overlaps.
pub fn connection_status(a: &GeneModel, b: &GeneModel) -> &'static str {
    let (x, y) = (a.query_interval_aa, b.query_interval_aa);
    if a.candidate_id == b.candidate_id
        || x.start >= y.start
        || x.end >= y.end
        || a.identity_millionths.abs_diff(b.identity_millionths) > 150_000
    {
        return "order_contig_or_identity";
    }
    let overlap = x.overlap(y);
    // De Bruijn fragments commonly share k-1 bases. A fixed 12-aa / one-third
    // cap rejects legitimate long exact overlaps (e.g. k=51 unitigs). Require
    // exact coordinate-compatible coding agreement instead of a length cap.
    if overlap > 0 {
        // Never discard a discordant overlap simply to obtain a longer model.
        if let (Some(left), Some(right)) = (coding(a), coding(b)) {
            let nt = overlap * 3;
            if left[left.len() - nt..] != right[..nt] {
                return "coding_overlap_conflict";
            }
        } else {
            return "coding_coordinate_map_missing";
        }
    }
    "compatible"
}

fn covered(path: &[usize], models: &[GeneModel]) -> usize {
    let mut end: usize = 0;
    let mut count = 0;
    for &i in path {
        let q = models[i].query_interval_aa;
        count += q.end.saturating_sub(q.start.max(end));
        end = end.max(q.end);
    }
    count
}

/// Bounded DAG search; competing near-equal paths are reported, never silently
/// chosen by score. A contig cannot re-enter a chain. Search overflow is unresolved.
pub fn plan_chains(models: &[GeneModel]) -> Vec<FragmentChain> {
    plan_chains_with_veto(models, |_, _| false)
}

/// Apply negative read evidence before path ranking and bounding, so an
/// unsupported winning path cannot hide a viable alternative. Missing read
/// evidence must never be used as a veto for a long intron.
pub fn plan_chains_with_veto(
    models: &[GeneModel],
    veto: impl FnMut(&GeneModel, &GeneModel) -> bool,
) -> Vec<FragmentChain> {
    plan_chains_impl::<true>(models, veto)
}

fn plan_chains_impl<const EARLY_BOUND: bool>(
    models: &[GeneModel],
    mut veto: impl FnMut(&GeneModel, &GeneModel) -> bool,
) -> Vec<FragmentChain> {
    let complete: BTreeSet<_> = models
        .iter()
        .filter(|m| m.eligible_for_resolve && m.state == ModelState::Complete)
        .map(|m| m.protein_id.as_str())
        .collect();
    let mut groups: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (i, m) in models.iter().enumerate() {
        if source_eligible(m) && !complete.contains(m.protein_id.as_str()) {
            groups.entry(&m.protein_id).or_default().push(i);
        }
    }
    let mut result = Vec::new();
    for (protein, mut nodes) in groups {
        // Compare fragments against a common protein before cross-reference
        // selection; otherwise each contig can choose a different species' ortholog.
        let mut local: Vec<_> = nodes.iter().map(|&i| models[i].clone()).collect();
        // Global competition mixes ortholog references. Recompute it within
        // this reference, without carrying its previous state/flags forward.
        // This does not accept a whole-gene model: these are chain hypotheses.
        for model in &mut local {
            model.qc_flags.retain(|flag| {
                !matches!(
                    flag.as_str(),
                    "ambiguous_competitor" | "competing_model_not_selected"
                )
            });
            if model.state == ModelState::AmbiguousModel {
                model.state = ModelState::TerminalPartial;
            }
        }
        super::select_competing_models(&mut local);
        nodes = nodes
            .into_iter()
            .zip(local)
            .filter(|(_, m)| m.selected && source_eligible(m))
            .map(|(i, _)| i)
            .collect();
        nodes.sort_by_key(|&i| (models[i].query_interval_aa, i));
        if nodes.len() > MAX_NODES {
            result.push(FragmentChain {
                protein: protein.into(),
                models: Vec::new(),
                covered_aa: 0,
                ambiguous: true,
                search_limited: true,
            });
            continue;
        }
        let mut paths: Vec<Vec<Vec<usize>>> = Vec::new();
        let mut limited = false;
        for (j, &node) in nodes.iter().enumerate() {
            let mut ending = vec![vec![node]];
            for k in 0..j {
                if !compatible(&models[nodes[k]], &models[node])
                    || veto(&models[nodes[k]], &models[node])
                {
                    continue;
                }
                for path in &paths[k] {
                    if path
                        .iter()
                        .any(|&i| models[i].candidate_id == models[node].candidate_id)
                    {
                        continue;
                    }
                    let mut next = path.clone();
                    next.push(node);
                    ending.push(next);
                }
            }
            let mut ending: Vec<_> = ending
                .into_iter()
                .map(|p| (covered(&p, models), p))
                .collect();
            ending.sort_by(|(a, p), (b, q)| {
                b.cmp(a)
                    .then_with(|| p.len().cmp(&q.len()))
                    .then_with(|| p.cmp(q))
            });
            // A worse subset reaching this same endpoint cannot improve a
            // future extension unless its omitted contig can reappear later.
            // Remove those dominated paths before bounding the search, so ten
            // non-overlapping exons do not enumerate every possible subset.
            let future_contigs: BTreeSet<_> = nodes[j + 1..]
                .iter()
                .map(|&i| models[i].candidate_id.as_str())
                .collect();
            let mut maximal: Vec<(usize, Vec<usize>)> = Vec::new();
            for (score, path) in ending {
                if maximal.iter().any(|(winner_score, winner)| {
                    *winner_score >= score
                        && (winner.iter().all(|i| path.contains(i))
                            || (path.iter().all(|i| winner.contains(i))
                                && winner.iter().filter(|i| !path.contains(i)).all(|&i| {
                                    !future_contigs.contains(models[i].candidate_id.as_str())
                                })))
                }) {
                    continue;
                }
                maximal.push((score, path));
                // Survivors are append-only in sorted order. Once one more
                // than the retained prefix survives, overflow is certain and
                // later paths cannot change that prefix or the ambiguity flag.
                if EARLY_BOUND && maximal.len() > MAX_PATHS {
                    break;
                }
            }
            let mut ending: Vec<_> = maximal.into_iter().map(|(_, p)| p).collect();
            if ending.len() > MAX_PATHS {
                limited = true;
                ending.truncate(MAX_PATHS);
            }
            paths.push(ending);
        }
        let mut candidates: Vec<_> = paths
            .into_iter()
            .flatten()
            .filter(|p| p.len() > 1)
            .collect();
        candidates
            .sort_by_cached_key(|p| (std::cmp::Reverse(covered(p, models)), p.len(), p.clone()));
        if let Some(best) = candidates.first() {
            let coverage = covered(best, models);
            let tolerance = (models[best[0]].protein_len_aa / 50).max(1);
            let ambiguous = limited
                || candidates.iter().skip(1).any(|p| {
                    covered(p, models).saturating_add(tolerance) >= coverage
                        && !p.iter().all(|i| best.contains(i))
                });
            result.push(FragmentChain {
                protein: protein.into(),
                models: best.clone(),
                covered_aa: coverage,
                ambiguous,
                search_limited: limited,
            });
        }
    }
    result.sort_by(|a, b| {
        let fraction = |c: &FragmentChain| {
            c.models
                .first()
                .map(|&i| c.covered_aa as f64 / models[i].protein_len_aa as f64)
                .unwrap_or(0.0)
        };
        fraction(b)
            .total_cmp(&fraction(a))
            .then_with(|| a.protein.cmp(&b.protein))
    });
    let mut seen = BTreeSet::new();
    result.retain(|c| {
        c.ambiguous
            || seen.insert(
                c.models
                    .iter()
                    .map(|&i| (models[i].candidate_id.clone(), models[i].strand))
                    .collect::<Vec<_>>(),
            )
    });
    result
}

/// Retain observed splice flanks at internal joins without carrying an arbitrary
/// contig tail (which can contain an unrelated competing locus) into reannotation.
/// Returns the oriented sequence, its source model anchor, and retained original coordinates.
pub fn join_source(
    sequence: &str,
    model: &GeneModel,
    first: bool,
    last: bool,
) -> (String, super::Interval, super::Interval) {
    let oriented = if model.strand == super::Strand::Forward {
        sequence.to_owned()
    } else {
        super::reverse_complement(sequence)
    };
    let anchor = if model.strand == super::Strand::Forward {
        model.genomic_interval
    } else {
        super::Interval {
            start: sequence.len() - model.genomic_interval.end,
            end: sequence.len() - model.genomic_interval.start,
        }
    };
    let start = if first {
        0
    } else {
        anchor.start.saturating_sub(120)
    };
    let end = if last {
        sequence.len()
    } else {
        (anchor.end + 120).min(sequence.len())
    };
    let retained = if model.strand == super::Strand::Forward {
        super::Interval { start, end }
    } else {
        super::Interval {
            start: sequence.len() - end,
            end: sequence.len() - start,
        }
    };
    (
        oriented[start..end].to_owned(),
        super::Interval {
            start: anchor.start - start,
            end: anchor.end - start,
        },
        retained,
    )
}

/// CDS-only hypothesis with codon-sized Ns for missing *coding* sequence.
/// This is never promoted to a validated genomic model by itself.
pub fn stitch_cds(chain: &FragmentChain, models: &[GeneModel]) -> Option<(String, usize, usize)> {
    if chain.ambiguous || chain.models.len() < 2 {
        return None;
    }
    let mut sequence = String::new();
    let mut end: usize = 0;
    let mut gaps = 0;
    let mut trimmed = 0;
    for (part, &i) in chain.models.iter().enumerate() {
        let m = &models[i];
        let q = m.query_interval_aa;
        let cds = coding(m)?;
        let overlap = if part == 0 {
            0
        } else {
            end.saturating_sub(q.start)
        };
        let gap = if part == 0 {
            0
        } else {
            q.start.saturating_sub(end)
        };
        if overlap > 0 && sequence[sequence.len() - overlap * 3..] != cds[..overlap * 3] {
            return None;
        }
        sequence.push_str(&"N".repeat(gap * 3));
        sequence.push_str(&cds[overlap * 3..]);
        gaps += gap;
        trimmed += overlap * 3;
        end = q.end;
    }
    Some((sequence, gaps, trimmed))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gene_annotation::{exonerate::parse_vulgar, ModelConfig};
    fn model(candidate: &str, start: usize, end: usize) -> GeneModel {
        let target = "GCT".repeat(end - start);
        let protein = "A".repeat(100);
        let hit = parse_vulgar(
            &format!(
                "vulgar: p {start} {end} . {candidate} 0 {} + 100 M {} {}",
                target.len(),
                end - start,
                target.len()
            ),
            candidate,
            target.len(),
            &BTreeMap::from([("p".into(), 100)]),
        )
        .unwrap()
        .remove(0);
        let mut m = hit
            .normalize(candidate, &target, &protein)
            .unwrap()
            .build(&target, ModelConfig::default())
            .unwrap();
        m.selected = true;
        m
    }
    #[test]
    fn early_path_bound_preserves_exhaustive_prefix_and_overflow_with_vetoes() {
        for choices in [2, 3, 4] {
            let models: Vec<_> = (0..5)
                .flat_map(|layer| {
                    (0..choices).map(move |copy| {
                        model(
                            &format!("contig_{layer}_{copy}"),
                            layer * 15,
                            layer * 15 + 12,
                        )
                    })
                })
                .collect();
            for vetoed in [false, true] {
                let veto = |a: &GeneModel, b: &GeneModel| {
                    vetoed && a.candidate_id.ends_with('1') && b.candidate_id.ends_with('2')
                };
                let before = plan_chains_impl::<false>(&models, veto);
                let after = plan_chains_impl::<true>(&models, veto);
                assert_eq!(before, after);
                if choices == 4 {
                    assert!(after.iter().any(|c| c.search_limited));
                }
            }
        }
    }

    #[test]
    fn short_internal_exons_survive_reference_local_competition() {
        let mut models = vec![model("a", 0, 30), model("b", 30, 60), model("c", 60, 100)];
        assert_eq!(models[1].state, ModelState::TerminalPartial);
        let chain = plan_chains(&models).remove(0);
        assert_eq!(chain.models, vec![0, 1, 2]);
        assert!(!chain.ambiguous);
        // A short internal exon marked ambiguous against a different reference
        // is reconsidered only after competition is scoped to one reference.
        models[1].state = ModelState::AmbiguousModel;
        models[1].qc_flags.push("ambiguous_competitor".into());
        let chain = plan_chains(&models).remove(0);
        assert_eq!(chain.models, vec![0, 1, 2]);
        assert!(!chain.ambiguous);
        // Genuine same-reference competitors must still block this bridge.
        let mut competitor = models[1].clone();
        competitor.cds.replace_range(2..3, "C");
        models.push(competitor);
        let chains = plan_chains(&models);
        assert!(chains
            .iter()
            .all(|c| !c.models.contains(&1) && !c.models.contains(&3)));
        assert!(models[1]
            .qc_flags
            .iter()
            .any(|f| f == "ambiguous_competitor"));
    }

    #[test]
    fn read_veto_precedes_best_path_selection() {
        let models = vec![
            model("a", 0, 45),
            model("b", 45, 55),
            model("copy", 45, 55),
            model("c", 55, 100),
        ];
        assert!(plan_chains(&models)[0].ambiguous);
        let chain = plan_chains_with_veto(&models, |a, b| {
            a.candidate_id == "copy" || b.candidate_id == "copy"
        })
        .remove(0);
        assert!(!chain.ambiguous);
        assert_eq!(chain.models, vec![0, 1, 3]);
        assert_eq!(chain.covered_aa, 100);
        // No evidence must retain the unresolved competing copy.
        assert!(plan_chains_with_veto(&models, |_, _| false)[0].ambiguous);
    }

    #[test]
    fn internal_join_trims_only_unaligned_tails_and_preserves_reverse_anchor() {
        let mut m = model("a", 0, 30);
        m.genomic_interval = super::super::Interval {
            start: 1000,
            end: 1090,
        };
        let sequence = "A".repeat(3000);
        let (seq, anchor, retained) = join_source(&sequence, &m, false, false);
        assert_eq!(seq.len(), 330);
        assert_eq!((anchor.start, anchor.end), (120, 210));
        assert_eq!((retained.start, retained.end), (880, 1210));
        m.strand = super::super::Strand::Reverse;
        let (seq, anchor, retained) = join_source(&sequence, &m, false, false);
        assert_eq!(seq, "T".repeat(330));
        assert_eq!((anchor.start, anchor.end), (120, 210));
        assert_eq!((retained.start, retained.end), (880, 1210));
        let (seq, _, _) = join_source(&sequence, &m, true, true);
        assert_eq!(seq.len(), 3000);
    }
    #[test]
    fn dominance_preserves_an_omitted_contig_needed_by_a_future_extension() {
        let mut later = model("b", 50, 100);
        later.genomic_interval = crate::gene_annotation::Interval::new(1000, 1150).unwrap();
        let models = vec![
            model("a", 0, 10),
            model("b", 10, 30),
            model("c", 30, 50),
            later,
        ];
        let chain = plan_chains(&models).remove(0);
        assert_eq!(chain.models, vec![0, 2, 3]);
        assert_eq!(chain.covered_aa, 80);
    }

    #[test]
    fn ten_disconnected_exons_do_not_exhaust_search_on_redundant_subpaths() {
        let models: Vec<_> = (0..10)
            .map(|i| model(&format!("c{i}"), i * 10, (i + 1) * 10))
            .collect();
        let chain = plan_chains(&models).remove(0);
        assert!(!chain.ambiguous);
        assert!(!chain.search_limited);
        assert_eq!(chain.models.len(), 10);
        assert_eq!(chain.covered_aa, 100);
    }

    #[test]
    fn equivalent_reference_orthologs_share_one_physical_chain() {
        let mut models = vec![model("a", 0, 45), model("b", 45, 55), model("c", 55, 100)];
        let mut other = models.clone();
        for m in &mut other {
            m.protein_id = "other_reference".into();
            m.selected = false;
        }
        models.extend(other);
        let chains = plan_chains(&models);
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].models.len(), 3);
    }

    #[test]
    fn long_exact_coding_overlap_is_trimmed_once_but_disagreement_is_rejected() {
        let mut models = vec![model("a", 0, 70), model("b", 30, 100)];
        let chains = plan_chains(&models);
        assert_eq!(chains.len(), 1);
        let (sequence, gaps, trimmed) = stitch_cds(&chains[0], &models).unwrap();
        assert_eq!(sequence, "GCT".repeat(100));
        assert_eq!(gaps, 0);
        assert_eq!(trimmed, 120);
        models[1].cds.replace_range(0..3, "GCC");
        assert!(plan_chains(&models).is_empty());
    }
    #[test]
    fn short_middle_exon_joins_three_contigs() {
        let models = vec![model("a", 0, 45), model("b", 45, 55), model("c", 55, 100)];
        assert_eq!(models[1].state, ModelState::LowCoverage);
        let chain = plan_chains(&models).remove(0);
        assert!(!chain.ambiguous);
        assert_eq!(chain.models, vec![0, 1, 2]);
        assert_eq!(chain.covered_aa, 100);
        assert_eq!(stitch_cds(&chain, &models), Some(("GCT".repeat(100), 0, 0)));
    }
    #[test]
    fn alternative_copy_is_ambiguous_and_not_stitched() {
        let models = vec![
            model("a", 0, 45),
            model("b", 45, 55),
            model("copy", 45, 55),
            model("c", 55, 100),
        ];
        let chain = plan_chains(&models).remove(0);
        assert!(chain.ambiguous);
        assert!(stitch_cds(&chain, &models).is_none());
    }
    #[test]
    fn coding_gaps_and_identical_overlaps_are_explicit() {
        let models = vec![model("a", 0, 45), model("b", 42, 65), model("c", 70, 100)];
        let chain = plan_chains(&models).remove(0);
        let (cds, gaps, trim) = stitch_cds(&chain, &models).unwrap();
        assert_eq!(gaps, 5);
        assert_eq!(trim, 9);
        assert_eq!(
            cds,
            format!("{}{}{}", "GCT".repeat(65), "N".repeat(15), "GCT".repeat(30))
        );
    }
    #[test]
    fn discordant_overlap_and_repeated_contig_are_not_merged() {
        let mut models = vec![model("a", 0, 45), model("b", 42, 100)];
        models[1].cds.replace_range(..3, "GCC");
        assert!(plan_chains(&models).is_empty());
        let models = vec![model("a", 0, 45), model("b", 45, 55), model("a", 55, 100)];
        assert!(plan_chains(&models).iter().all(|c| c.models.len() < 3));
    }
}
