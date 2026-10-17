//! Replay saved shadow inputs through separate engine selectors. No cross-engine scores.
use gm2_tools::gene_annotation::{
    build_gene_model, exonerate::parse_vulgar, parse_miniprot_output, read_dna_fasta,
    select_competing_models, GeneModel, ModelConfig,
};
use std::{
    collections::BTreeMap,
    error::Error,
    fs,
    path::{Path, PathBuf},
};

fn files(path: &Path, name: &str, out: &mut Vec<PathBuf>) -> Result<(), Box<dyn Error>> {
    for entry in fs::read_dir(path)? {
        let path = entry?.path();
        if path.is_dir() {
            files(&path, name, out)?;
        } else if path.file_name().is_some_and(|v| v == name) {
            out.push(path);
        }
    }
    Ok(())
}
fn proteins(path: &Path) -> Result<BTreeMap<String, String>, Box<dyn Error>> {
    let mut result = BTreeMap::new();
    let mut id = String::new();
    for line in fs::read_to_string(path)?.lines() {
        if let Some(header) = line.strip_prefix('>') {
            id = header
                .split_whitespace()
                .next()
                .ok_or("empty header")?
                .into();
            result.insert(id.clone(), String::new());
        } else if !line.trim().is_empty() {
            result
                .get_mut(&id)
                .ok_or("sequence without header")?
                .push_str(line.trim());
        }
    }
    Ok(result)
}
fn accepted(models: &[GeneModel]) -> Vec<&GeneModel> {
    models
        .iter()
        .filter(|m| m.selected && m.eligible_for_resolve)
        .collect()
}
fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 3 {
        return Err("usage: validate_exonerate_selection SHADOW_OUTPUT REPORT.tsv".into());
    }
    let root = Path::new(&args[1]);
    let shadow = root.join("raw_exonerate_shadow");
    let mut inputs = Vec::new();
    files(&shadow, "target.fasta", &mut inputs)?;
    inputs.sort();
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&args[2])?;
    writer.write_record([
        "comparison",
        "miniprot_models",
        "exonerate_models",
        "miniprot_accepted",
        "exonerate_accepted",
        "exact_cds_matches",
        "terminal_stop_normalized_matches",
        "unmatched_exonerate",
        "missing_miniprot",
        "miniprot_ambiguous",
        "exonerate_ambiguous",
    ])?;
    let mut details = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(format!("{}.models.tsv", args[2]))?;
    details.write_record([
        "comparison",
        "engine",
        "protein",
        "strand",
        "start",
        "end",
        "state",
        "selected",
        "eligible",
        "qc",
        "protein_sequence",
    ])?;
    for input in inputs {
        let directory = input.parent().ok_or("missing directory")?;
        let relative = directory.strip_prefix(&shadow)?;
        let target = read_dna_fasta(&input)?.remove(0);
        let proteins = proteins(&directory.join("proteins.faa"))?;
        let lengths = proteins.iter().map(|(k, v)| (k.clone(), v.len())).collect();
        let output = fs::read_to_string(
            root.join("raw_miniprot")
                .join(relative)
                .with_extension("gff3"),
        )?;
        let parsed = parse_miniprot_output(&output);
        if !parsed.warnings.is_empty() || parsed.unassigned_paf != 0 {
            return Err(format!("miniprot parse failure: {}", relative.display()).into());
        }
        let mut mini = parsed
            .models
            .iter()
            .map(|raw| {
                build_gene_model(
                    raw,
                    &target.sequence,
                    proteins.get(&raw.protein_id).map(String::len),
                    ModelConfig::default(),
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let status = fs::read_to_string(directory.join("status.txt"))?;
        if !status.starts_with("ok\n") {
            return Err(format!("Exonerate failed: {}", relative.display()).into());
        }
        let hits = parse_vulgar(
            &fs::read_to_string(directory.join("stdout.txt"))?,
            &target.id,
            target.sequence.len(),
            &lengths,
        )?;
        let mut exo = hits
            .iter()
            .map(|h| {
                h.normalize(&target.id, &target.sequence, &proteins[&h.protein])?
                    .build(&target.sequence, ModelConfig::default())
            })
            .collect::<Result<Vec<_>, _>>()?;
        for model in &mut mini {
            model.check_reference_start(
                proteins
                    .get(&model.protein_id)
                    .and_then(|p| p.bytes().next()),
            );
        }
        select_competing_models(&mut mini);
        select_competing_models(&mut exo);
        for (engine, models) in [("miniprot", &mini), ("exonerate", &exo)] {
            for m in models {
                details.write_record([
                    relative.display().to_string(),
                    engine.into(),
                    m.protein_id.clone(),
                    m.strand.as_char().to_string(),
                    m.genomic_interval.start.to_string(),
                    m.genomic_interval.end.to_string(),
                    m.state.as_str().into(),
                    m.selected.to_string(),
                    m.eligible_for_resolve.to_string(),
                    m.qc_flags.join(","),
                    m.protein.clone().unwrap_or_default(),
                ])?;
            }
        }
        let ma = accepted(&mini);
        let ea = accepted(&exo);
        let same = |a: &&GeneModel, b: &&GeneModel| {
            a.protein_id == b.protein_id
                && a.strand == b.strand
                && a.genomic_interval.overlap(b.genomic_interval) > 0
        };
        let without_terminal_stop = |model: &&GeneModel| {
            if model.protein.as_deref().is_some_and(|p| p.ends_with('*')) {
                model.cds[..model.cds.len() - 3].to_owned()
            } else {
                model.cds.clone()
            }
        };
        let equivalent = |a: &&GeneModel, b: &&GeneModel| {
            same(a, b) && without_terminal_stop(a) == without_terminal_stop(b)
        };
        let exact = ea
            .iter()
            .filter(|e| ma.iter().any(|m| same(e, m) && e.cds == m.cds))
            .count();
        let normalized = ea
            .iter()
            .filter(|e| ma.iter().any(|m| equivalent(e, m)))
            .count();
        let missing = ma
            .iter()
            .filter(|m| !ea.iter().any(|e| equivalent(e, m)))
            .count();
        writer.write_record([
            relative.display().to_string(),
            mini.len().to_string(),
            exo.len().to_string(),
            ma.len().to_string(),
            ea.len().to_string(),
            exact.to_string(),
            normalized.to_string(),
            (ea.len() - normalized).to_string(),
            missing.to_string(),
            mini.iter()
                .filter(|m| m.qc_flags.iter().any(|f| f == "ambiguous_competitor"))
                .count()
                .to_string(),
            exo.iter()
                .filter(|m| m.qc_flags.iter().any(|f| f == "ambiguous_competitor"))
                .count()
                .to_string(),
        ])?;
    }
    writer.flush()?;
    details.flush()?;
    Ok(())
}
