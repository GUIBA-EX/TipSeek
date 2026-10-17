//! Original versus indexed selection on frozen HX50 models.
use std::{cmp::Ordering, collections::BTreeMap};

use gm2_tools::gene_annotation::{GeneModel, ModelState};

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
fn legacy_select_competing_models(models: &mut [GeneModel]) {
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
        let ambiguous: Vec<bool> = group
            .iter()
            .map(|index| *index == seed || nearly_tied(&models[seed], &models[*index]))
            .collect();
        if ambiguous.iter().skip(1).any(|value| *value) {
            for (index, is_ambiguous) in group.into_iter().zip(ambiguous) {
                models[index].eligible_for_resolve = false;
                if is_ambiguous {
                    models[index].state = ModelState::AmbiguousModel;
                    models[index]
                        .qc_flags
                        .push("ambiguous_competitor".to_owned());
                } else {
                    models[index]
                        .qc_flags
                        .push("competing_model_not_selected".to_owned());
                }
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    use gm2_tools::{
        fasta::read_fasta,
        gene_annotation::{
            build_gene_model, parse_miniprot_output, read_dna_fasta, select_competing_models,
            ModelConfig,
        },
    };
    use std::{
        collections::BTreeSet,
        fs::{self, File},
        io::{BufReader, BufWriter, Write},
        path::PathBuf,
        time::Instant,
    };
    let args: Vec<_> = std::env::args().collect();
    let annotation = PathBuf::from(&args[1]);
    let gene = PathBuf::from(&args[2]);
    let families: BTreeSet<_> = fs::read_to_string(&args[3])?
        .lines()
        .map(str::to_owned)
        .collect();
    let mut ids: BTreeMap<String, Vec<(String, String)>> = BTreeMap::new();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(annotation.join("manifest/id_map.tsv"))?;
    for row in reader.records() {
        let r = row?;
        if families.contains(&r[1]) {
            ids.entry(r[1].into())
                .or_default()
                .push((r[2].into(), r[3].into()));
        }
    }
    let mut out = BufWriter::new(File::create(&args[4])?);
    writeln!(
        out,
        "family_id\tmodels\trepeats\tlegacy_seconds\tindexed_seconds\tidentical_decisions"
    )?;
    for (family, ids) in ids {
        let dna: BTreeMap<_, _> =
            read_dna_fasta(&gene.join(format!("samples/1_Hx50/candidates/{family}.fasta")))?
                .into_iter()
                .map(|r| (r.id.split('|').next_back().unwrap().to_owned(), r.sequence))
                .collect();
        let lengths: BTreeMap<_, _> = read_fasta(BufReader::new(File::open(
            annotation.join(format!("manifest/derived_proteins/{family}.faa")),
        )?))?
        .into_iter()
        .map(|r| {
            (
                r.name.split_whitespace().next().unwrap().to_owned(),
                r.sequence.len(),
            )
        })
        .collect();
        let mut models = Vec::new();
        for (internal, candidate) in ids {
            let path = annotation.join(format!("raw_miniprot/1_Hx50/{family}/{internal}.gff3"));
            let Some(sequence) = dna.get(&candidate) else {
                continue;
            };
            if !path.is_file() {
                continue;
            }
            for mut raw in parse_miniprot_output(&fs::read_to_string(path)?).models {
                raw.candidate_id = candidate.clone();
                if let Ok(model) = build_gene_model(
                    &raw,
                    sequence,
                    lengths.get(&raw.protein_id).copied(),
                    ModelConfig::default(),
                ) {
                    models.push(model);
                }
            }
        }
        let mut legacy = 0.0;
        let mut indexed = 0.0;
        for round in 0..10 {
            let mut a = models.clone();
            let mut b = models.clone();
            if round % 2 == 0 {
                let t = Instant::now();
                legacy_select_competing_models(&mut a);
                legacy += t.elapsed().as_secs_f64();
                let t = Instant::now();
                select_competing_models(&mut b);
                indexed += t.elapsed().as_secs_f64();
            } else {
                let t = Instant::now();
                select_competing_models(&mut b);
                indexed += t.elapsed().as_secs_f64();
                let t = Instant::now();
                legacy_select_competing_models(&mut a);
                legacy += t.elapsed().as_secs_f64();
            }
            for (x, y) in a.iter().zip(&b) {
                assert_eq!(
                    (
                        x.selected,
                        x.state,
                        x.eligible_for_resolve,
                        x.competition_group,
                        &x.qc_flags
                    ),
                    (
                        y.selected,
                        y.state,
                        y.eligible_for_resolve,
                        y.competition_group,
                        &y.qc_flags
                    )
                );
            }
        }
        writeln!(
            out,
            "{family}\t{}\t10\t{legacy:.9}\t{indexed:.9}\ttrue",
            models.len()
        )?;
    }
    out.flush()?;
    Ok(())
}
