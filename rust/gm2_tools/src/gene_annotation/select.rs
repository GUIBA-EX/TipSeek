use std::cmp::Ordering;

use super::{GeneModel, ModelState};

fn competing(left: &GeneModel, right: &GeneModel) -> bool {
    if left.candidate_id != right.candidate_id {
        return false;
    }
    let genomic_shorter = left
        .genomic_interval
        .len()
        .min(right.genomic_interval.len());
    left.genomic_interval.overlap(right.genomic_interval) * 2 >= genomic_shorter
}

fn model_order(left: &GeneModel, right: &GeneModel) -> Ordering {
    left.structural_defect_count()
        .cmp(&right.structural_defect_count())
        .then_with(|| right.coverage_millionths.cmp(&left.coverage_millionths))
        .then_with(|| right.alignment_score.cmp(&left.alignment_score))
        .then_with(|| right.identity_millionths.cmp(&left.identity_millionths))
        .then_with(|| {
            left.noncanonical_introns()
                .cmp(&right.noncanonical_introns())
        })
        .then_with(|| left.large_query_gaps.cmp(&right.large_query_gaps))
        .then_with(|| left.rank.cmp(&right.rank))
        .then_with(|| left.gff_id.cmp(&right.gff_id))
}

fn nearly_tied(left: &GeneModel, right: &GeneModel) -> bool {
    left.structural_defect_count() == right.structural_defect_count()
        && left.coverage_millionths.abs_diff(right.coverage_millionths) < 20_000
        && left.alignment_score.abs_diff(right.alignment_score)
            < (left.alignment_score.unsigned_abs().max(1) / 20).max(1)
        && left.cds != right.cds
}

/// Select the best model among direct competitors. A rejected bridge model
/// cannot merge two otherwise independent loci into one competition group.
/// Near-tied models with different CDS sequences are retained for audit but
/// excluded from resolve.
pub fn select_competing_models(models: &mut [GeneModel]) {
    for model in models.iter_mut() {
        model.selected = false;
        model.competition_group = 0;
    }
    let mut order: Vec<_> = (0..models.len()).collect();
    order.sort_by(|left, right| model_order(&models[*left], &models[*right]));
    let mut assigned = vec![false; models.len()];
    let mut group_number = 0usize;
    for seed in order {
        if assigned[seed] {
            continue;
        }
        group_number += 1;
        let mut group = vec![seed];
        assigned[seed] = true;
        let mut competitors: Vec<_> = (0..models.len())
            .filter(|other| !assigned[*other] && competing(&models[seed], &models[*other]))
            .collect();
        competitors.sort_by(|left, right| model_order(&models[*left], &models[*right]));
        for competitor in competitors {
            assigned[competitor] = true;
            group.push(competitor);
        }
        for index in &group {
            models[*index].competition_group = group_number;
        }
        if group.len() > 1 && nearly_tied(&models[group[0]], &models[group[1]]) {
            for index in group.iter().copied().take(2) {
                models[index].state = ModelState::AmbiguousModel;
                models[index].eligible_for_resolve = false;
                models[index]
                    .qc_flags
                    .push("ambiguous_competitor".to_owned());
            }
            for index in group.into_iter().skip(2) {
                models[index].eligible_for_resolve = false;
                models[index]
                    .qc_flags
                    .push("competing_model_not_selected".to_owned());
            }
        } else {
            let winner = group[0];
            models[winner].selected = true;
            for index in group.into_iter().skip(1) {
                models[index].eligible_for_resolve = false;
                models[index]
                    .qc_flags
                    .push("competing_model_not_selected".to_owned());
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gene_annotation::{GeneModel, Interval, ModelState, Strand};

    fn model(
        id: &str,
        genomic: Interval,
        query: Interval,
        score: i64,
        defective: bool,
    ) -> GeneModel {
        GeneModel {
            gff_id: id.into(),
            candidate_id: "c1".into(),
            protein_id: "p1".into(),
            protein_len_aa: 100,
            query_interval_aa: query,
            genomic_interval: genomic,
            strand: Strand::Forward,
            alignment_score: score,
            identity_millionths: 900_000,
            positive_millionths: 950_000,
            coverage_millionths: (query.len() as u32) * 10_000,
            rank: 0,
            frameshifts: usize::from(defective),
            inframe_stops: 0,
            phase_valid: true,
            five_prime_trim_nt: 0,
            three_prime_trim_nt: 0,
            large_query_gaps: 0,
            exons: Vec::new(),
            introns: Vec::new(),
            cds: format!("ATG{id}"),
            protein: Some("M".into()),
            gene: String::new(),
            state: if defective {
                ModelState::Frameshifted
            } else {
                ModelState::Complete
            },
            qc_flags: Vec::new(),
            selected: false,
            eligible_for_resolve: !defective,
            competition_group: 0,
        }
    }

    #[test]
    fn clean_model_beats_longer_frameshifted_model() {
        let mut models = vec![
            model(
                "bad",
                Interval::new(0, 300).unwrap(),
                Interval::new(0, 100).unwrap(),
                200,
                true,
            ),
            model(
                "good",
                Interval::new(10, 280).unwrap(),
                Interval::new(5, 95).unwrap(),
                180,
                false,
            ),
        ];
        select_competing_models(&mut models);
        assert!(!models[0].selected);
        assert!(models[1].selected);
    }

    #[test]
    fn nonoverlapping_loci_are_both_retained() {
        let mut models = vec![
            model(
                "one",
                Interval::new(0, 90).unwrap(),
                Interval::new(0, 30).unwrap(),
                100,
                false,
            ),
            model(
                "two",
                Interval::new(200, 290).unwrap(),
                Interval::new(60, 90).unwrap(),
                100,
                false,
            ),
        ];
        select_competing_models(&mut models);
        assert!(models.iter().all(|model| model.selected));
        assert_ne!(models[0].competition_group, models[1].competition_group);
    }

    #[test]
    fn tandem_loci_covering_the_same_protein_are_both_retained() {
        let mut models = vec![
            model(
                "one",
                Interval::new(0, 300).unwrap(),
                Interval::new(0, 100).unwrap(),
                100,
                false,
            ),
            model(
                "two",
                Interval::new(500, 800).unwrap(),
                Interval::new(0, 100).unwrap(),
                100,
                false,
            ),
        ];
        select_competing_models(&mut models);
        assert!(models.iter().all(|model| model.selected));
    }

    #[test]
    fn rejected_bridge_does_not_merge_independent_loci() {
        let mut models = vec![
            model(
                "left",
                Interval::new(0, 100).unwrap(),
                Interval::new(0, 34).unwrap(),
                200,
                false,
            ),
            model(
                "bridge",
                Interval::new(50, 250).unwrap(),
                Interval::new(0, 67).unwrap(),
                300,
                true,
            ),
            model(
                "right",
                Interval::new(200, 300).unwrap(),
                Interval::new(66, 100).unwrap(),
                190,
                false,
            ),
        ];
        select_competing_models(&mut models);
        assert!(models[0].selected);
        assert!(!models[1].selected);
        assert!(models[2].selected);
        assert_ne!(models[0].competition_group, models[2].competition_group);
    }
}
