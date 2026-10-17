//! Replay saved alignments, never invoke an aligner. One backbone per family,
//! chosen over the entire cohort; original candidate DNA is the only base source.
use super::*;
use gm2_tools::gene_annotation::exon_recovery;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

const HEADER: &str = "sample\tfamily_id\tprotein_reference\tlength_nt\tobserved_nt\tmissing_nt\tconflict_nt\tprojected_models\tunmappable_models\tstatus\tcross_reference_models\tcoordinate_conflict_codons\treference_selection\tdominated_models";

struct Input {
    sample: String,
    fasta: PathBuf,
    ids: Vec<(String, String, usize)>,
}
struct Job {
    family: String,
    samples: Vec<Input>,
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

// New runs also freeze external protein references here. Older runs without
// their exact protein snapshot fail rather than substitute a new reference.
fn protein_paths(annotation: &Path) -> io::Result<BTreeMap<String, PathBuf>> {
    let mut result = BTreeMap::new();
    let derived = annotation.join("manifest/derived_proteins");
    if derived.is_dir() {
        for path in sorted_directory_paths(&derived)? {
            if let Some(stem) = path.file_stem().and_then(|x| x.to_str()) {
                result.insert(stem.to_owned(), path);
            }
        }
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn recover(
    annotation: &Path,
    input: &Path,
    out: &Path,
    workers: usize,
    requested: Option<&str>,
    fixed_manifest: Option<&Path>,
    config: ModelConfig,
) -> io::Result<()> {
    recover_with_reads(
        annotation,
        input,
        out,
        workers,
        requested,
        fixed_manifest,
        config,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn recover_with_reads(
    annotation: &Path,
    input: &Path,
    out: &Path,
    workers: usize,
    requested: Option<&str>,
    fixed_manifest: Option<&Path>,
    config: ModelConfig,
    reads_root: Option<&Path>,
) -> io::Result<()> {
    let annotation = fs::canonicalize(annotation)?;
    let input = fs::canonicalize(input)?;
    if out.exists() && fs::canonicalize(out)? != annotation {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            "recovery destination already exists",
        ));
    }
    let in_place = out.exists();
    type CachedIds = BTreeMap<(String, String), Vec<(String, String, usize)>>;
    let mut ids: CachedIds = BTreeMap::new();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(annotation.join("manifest/id_map.tsv"))?;
    for record in reader.records() {
        let r = record?;
        ids.entry((r[0].to_owned(), r[1].to_owned()))
            .or_default()
            .push((
                r[2].to_owned(),
                r[3].to_owned(),
                r[4].parse()
                    .map_err(|_| invalid("invalid cached candidate length"))?,
            ));
    }
    let mut excluded = BTreeSet::new();
    let states = annotation.join("manifest/candidate_manifest.tsv");
    if states.is_file() {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(states)?;
        let header = reader.headers()?.clone();
        let state = header
            .iter()
            .position(|x| x == "state")
            .ok_or_else(|| invalid("cached manifest lacks state"))?;
        for record in reader.records() {
            let r = record?;
            if matches!(
                r.get(state),
                Some(
                    "protein_screen_rejected"
                        | "miniprot_failed"
                        | "exonerate_failed"
                        | "exonerate_rs_failed"
                )
            ) {
                excluded.insert((r[0].to_owned(), r[1].to_owned(), r[2].to_owned()));
            }
        }
    }
    let mut jobs: BTreeMap<String, Vec<Input>> = BTreeMap::new();
    for sample_path in sorted_directory_paths(&input.join("samples"))? {
        let candidates = sample_path.join("candidates");
        if !candidates.is_dir() {
            continue;
        }
        let sample = sample_path
            .file_name()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        for fasta in sorted_directory_paths(&candidates)? {
            let Some(family) = family_id(&fasta) else {
                continue;
            };
            // No usable protein reference produces no SCO in the original run.
            if !annotation
                .join("sco")
                .join(format!("{family}.fasta"))
                .exists()
            {
                continue;
            }
            jobs.entry(family.clone()).or_default().push(Input {
                sample: sample.clone(),
                fasta,
                ids: ids
                    .remove(&(sample.clone(), family))
                    .ok_or_else(|| invalid("missing cached ID map for annotated family"))?,
            });
        }
    }
    let jobs: Vec<_> = jobs
        .into_iter()
        .map(|(family, samples)| Job { family, samples })
        .collect();
    let proteins = protein_paths(&annotation)?;
    let mut stops = BTreeSet::new();
    let mut translations = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(annotation.join("manifest/reference_translation.tsv"))?;
    for record in translations.records() {
        let r = record?;
        if r.get(11) == Some("1") && r.get(13) == Some("1") {
            stops.insert((r[0].to_owned(), r[2].to_owned()));
        }
    }
    let mut fixed = BTreeMap::new();
    if let Some(path) = fixed_manifest {
        let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
        let header = reader.headers()?.clone();
        let family = header
            .iter()
            .position(|x| x == "family_id")
            .ok_or_else(|| invalid("backbone manifest lacks family_id"))?;
        let reference = header
            .iter()
            .position(|x| x == "protein_reference")
            .ok_or_else(|| invalid("backbone manifest lacks protein_reference"))?;
        for record in reader.records() {
            let r = record?;
            if let Some(old) = fixed.insert(r[family].to_owned(), r[reference].to_owned()) {
                if old != r[reference] {
                    return Err(invalid("inconsistent cohort backbone manifest"));
                }
            }
        }
        if jobs.iter().any(|job| !fixed.contains_key(&job.family)) {
            return Err(invalid(
                "backbone manifest does not cover every annotated family",
            ));
        }
    }
    // Write to a separate staging area even for in-place cohort normalization.
    let staging = out.join(".sco_rebuild");
    if staging.exists() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            "recovery staging area already exists",
        ));
    }
    fs::create_dir_all(staging.join("sco"))?;
    let next = AtomicUsize::new(0);
    let completed = AtomicUsize::new(0);
    let rows = Mutex::new(Vec::<(String, String, String, String)>::new());
    std::thread::scope(|scope| -> io::Result<()> {
        let handles: Vec<_> = (0..workers.max(1).min(jobs.len().max(1)))
            .map(|_| {
                scope.spawn(|| -> io::Result<()> {
                    loop {
                        let index = next.fetch_add(1, Ordering::Relaxed);
                        let Some(job) = jobs.get(index) else {
                            break;
                        };
                        let protein_path = proteins.get(&job.family).ok_or_else(|| {
                            invalid(format!(
                                "{}: exact cached reference proteins unavailable",
                                job.family
                            ))
                        })?;
                        let proteins: BTreeMap<String, String> = read_raw_fasta(protein_path)?
                            .into_iter()
                            .map(|(h, s)| (h.split_whitespace().next().unwrap().to_owned(), s))
                            .collect();
                        let lengths = proteins
                            .iter()
                            .map(|(id, s)| (id.clone(), s.len()))
                            .collect::<BTreeMap<_, _>>();
                        let mut sample_models = Vec::new();
                        let mut sample_dna = Vec::new();
                        let mut replay_audit = String::new();
                        for sample in &job.samples {
                            let dna: BTreeMap<_, _> = read_dna_fasta(&sample.fasta)?
                                .into_iter()
                                .map(|r| {
                                    (r.id.split('|').next_back().unwrap().to_owned(), r.sequence)
                                })
                                .collect();
                            let mut models = Vec::new();
                            let mut skipped = 0;
                            let mut raw_count = 0;
                            let mut rejected = 0;
                            for (internal, candidate, length) in &sample.ids {
                                // Derived padded contigs have no original observed DNA.
                                let Some(sequence) = dna.get(candidate) else {
                                    if internal.starts_with("TSKPAD") {
                                        continue;
                                    }
                                    return Err(invalid(
                                        "original cached candidate missing from DNA input",
                                    ));
                                };
                                if excluded.contains(&(
                                    sample.sample.clone(),
                                    job.family.clone(),
                                    candidate.clone(),
                                )) {
                                    skipped += 1;
                                    continue;
                                }
                                if sequence.len() != *length {
                                    return Err(invalid("cached DNA length mismatch"));
                                }
                                let mini = annotation
                                    .join("raw_miniprot")
                                    .join(&sample.sample)
                                    .join(&job.family)
                                    .join(format!("{internal}.gff3"));
                                let exo = annotation
                                    .join("raw_exonerate")
                                    .join(&sample.sample)
                                    .join(&job.family)
                                    .join(format!("{internal}.vulgar"));
                                let (backend, raw) = if mini.is_file() {
                                    (AnnotationBackend::Miniprot(String::new()), mini)
                                } else if exo.is_file() {
                                    (AnnotationBackend::Exonerate, exo)
                                } else {
                                    return Err(invalid(format!(
                                        "{}/{}: expected alignment cache missing for {internal}",
                                        sample.sample, job.family
                                    )));
                                };
                                raw_count += 1;
                                let report = backend.parse(
                                    &fs::read_to_string(raw)?,
                                    internal,
                                    sequence,
                                    &proteins,
                                )?;
                                for mut raw in report.models {
                                    if raw.candidate_id != *internal {
                                        return Err(invalid(
                                            "cached alignment targets unknown candidate",
                                        ));
                                    }
                                    raw.candidate_id = candidate.clone();
                                    match build_gene_model(
                                        &raw,
                                        sequence,
                                        lengths.get(&raw.protein_id).copied(),
                                        config,
                                    ) {
                                        Ok(mut model) => {
                                            model.check_reference_start(
                                                proteins
                                                    .get(&raw.protein_id)
                                                    .and_then(|p| p.bytes().next()),
                                            );
                                            models.push(model);
                                        }
                                        Err(_) => rejected += 1,
                                    }
                                }
                            }
                            replay_audit.push_str(&format!(
                                "{}\t{}\t{raw_count}\t{skipped}\t{rejected}\t{}\n",
                                sample.sample,
                                job.family,
                                models.len()
                            ));
                            if let Some(root) = reads_root {
                                audit_chains(
                                    &staging,
                                    root,
                                    &sample.sample,
                                    &job.family,
                                    &models,
                                    &dna,
                                )?;
                            }
                            sample_models.push(models);
                            sample_dna.push(dna);
                        }
                        // Candidate names are sample-local; only selection pools models.
                        let cohort = sample_models.iter().flatten();
                        let requested = fixed.get(&job.family).map(String::as_str).or(requested);
                        let (reference, selection) =
                            exon_recovery::choose_reference(cohort, &proteins, requested)
                                .map_err(invalid)?;
                        let mut manifest = String::new();
                        let mut projection = String::new();
                        let mut controls = String::new();
                        let family_path = staging.join("sco").join(format!("{}.fasta", job.family));
                        for ((sample, models), dna) in
                            job.samples.iter().zip(&sample_models).zip(&sample_dna)
                        {
                            let mut recovered = exon_recovery::reconstruct_with_candidates(
                                models,
                                &proteins,
                                Some(&reference),
                                dna,
                            )
                            .map_err(invalid)?;
                            recovered.reference_selection = selection;
                            if stops.contains(&(job.family.clone(), reference.clone())) {
                                recovered.append_terminal_stop(models);
                            }
                            append_fasta(
                                &family_path,
                                &format!(
                                    "{}|{} reference={} status=reference_coordinate_composite",
                                    sample.sample, job.family, reference
                                ),
                                &recovered.sequence,
                            )?;
                            write_masks(&staging, &job.family, &sample.sample, &recovered)?;
                            manifest.push_str(&row(&sample.sample, &job.family, &recovered));
                            for (source, map) in &recovered.maps {
                                projection.push_str(&format!(
                                    "{}\t{}\t{}\t{source}\t{}\t{}\t{}\t{}\n",
                                    sample.sample,
                                    job.family,
                                    reference,
                                    map.status,
                                    map.cost
                                        .map(|x| x.to_string())
                                        .unwrap_or_else(|| "NA".into()),
                                    map.positions.iter().flatten().count(),
                                    map.positions.len()
                                ));
                            }
                            // Separate backbone selection gains from cross-reference gains.
                            let mut direct =
                                exon_recovery::reconstruct(models, &lengths, Some(&reference))
                                    .map_err(invalid)?;
                            if stops.contains(&(job.family.clone(), reference.clone())) {
                                direct.append_terminal_stop(models);
                            }
                            controls.push_str(&format!(
                                "{}\t{}\t{}\t{}\t{}\n",
                                sample.sample,
                                job.family,
                                reference,
                                direct.sequence.len(),
                                direct.sequence.len() - direct.missing_nt - direct.conflict_nt
                            ));
                        }
                        rows.lock().unwrap().push((
                            job.family.clone(),
                            manifest,
                            projection,
                            format!("{controls}\u{0}{replay_audit}"),
                        ));
                        let count = completed.fetch_add(1, Ordering::Relaxed) + 1;
                        if count % 100 == 0 || count == jobs.len() {
                            eprintln!("SCO cache recovery completed {count}/{}", jobs.len());
                        }
                    }
                    Ok(())
                })
            })
            .collect();
        let mut failure = None;
        for h in handles {
            if let Err(e) = h
                .join()
                .map_err(|_| io::Error::other("recovery worker panicked"))
                .and_then(|r| r)
            {
                failure.get_or_insert(e);
            }
        }
        failure.map_or(Ok(()), Err)
    })?;
    let mut rows = rows.into_inner().unwrap();
    rows.sort_by(|a, b| a.0.cmp(&b.0));
    fs::create_dir_all(staging.join("manifest"))?;
    for (file,header,column) in [("sco_recovery.tsv",HEADER,1),("reference_projection.tsv","sample\tfamily_id\tprotein_reference\tsource_reference\tstatus\talignment_cost\tmapped_residues\tsource_residues",2),("sco_selection_control.tsv","sample\tfamily_id\tprotein_reference\tlength_nt\tdirect_only_observed_nt",3),("sco_replay.tsv","sample\tfamily_id\traw_candidates\tno_raw_candidates\tbuild_rejected\tmodels",4)] {
        let mut writer = BufWriter::new(File::create(staging.join("manifest").join(file))?);
        writeln!(writer,"{header}")?;
        for (_,manifest,projection,other) in &rows {
            let (control,replay) = other.split_once('\0').unwrap();
            writer.write_all(match column {1=>manifest,2=>projection,3=>control,_=>replay}.as_bytes())?;
        }
        writer.flush()?;
    }
    if in_place {
        fs::remove_dir_all(out.join("sco"))?;
    }
    fs::rename(staging.join("sco"), out.join("sco"))?;
    if out.join("sco_masks").exists() {
        fs::remove_dir_all(out.join("sco_masks"))?;
    }
    if staging.join("sco_masks").exists() {
        fs::rename(staging.join("sco_masks"), out.join("sco_masks"))?;
    }
    fs::create_dir_all(out.join("manifest"))?;
    for path in sorted_directory_paths(&staging.join("manifest"))? {
        fs::rename(&path, out.join("manifest").join(path.file_name().unwrap()))?;
    }
    if staging.join("chain_audit").exists() {
        if out.join("chain_audit").exists() {
            fs::remove_dir_all(out.join("chain_audit"))?;
        }
        fs::rename(staging.join("chain_audit"), out.join("chain_audit"))?;
    }
    fs::remove_dir_all(staging)?;
    Ok(())
}

// Recompute read-vetoed chain hypotheses from raw caches. Sequence recovery
// remains independent: a reference-ordered chain is not physical linkage.
fn audit_chains(
    out: &Path,
    root: &Path,
    sample: &str,
    family: &str,
    models: &[GeneModel],
    dna: &BTreeMap<String, String>,
) -> io::Result<()> {
    use gm2_tools::gene_annotation::{chain, reads::PairAudit};
    let contigs = dna
        .iter()
        .map(|(id, s)| (id.clone(), (id.clone(), s.clone())))
        .collect();
    let reads = PairAudit::load(Some(root), sample, family, &contigs)?;
    // Production plans chains before global competition; mirror that order.
    let selected = models.to_vec();
    // Exact previous source filter; the new planner is otherwise unchanged.
    let previous: Vec<_> = selected
        .iter()
        .filter(|m| {
            m.state == ModelState::Complete
                || (matches!(
                    m.state,
                    ModelState::TerminalPartial | ModelState::LowCoverage
                ) && !m.qc_flags.iter().any(|f| f == "ambiguous_competitor"))
        })
        .cloned()
        .collect();
    let directory = out.join("chain_audit").join(sample);
    fs::create_dir_all(&directory)?;
    let mut writer = BufWriter::new(File::create(directory.join(format!("{family}.tsv")))?);
    writeln!(writer, "mode\tfamily_id\tprotein_reference\tcontigs\tcovered_aa\tambiguous\tsearch_limited\tpaired_supported_edges\tdiscordant_edges\treads_available\ttotal_pairs\tmembers\tunique_mapped_pairs\tdistinct_cross_contig_pairs")?;
    for (mode, sources) in [("before", &previous), ("after", &selected)] {
        let chains = chain::plan_chains_with_veto(sources, |a, b| reads.edge(a, b).warning());
        if chains.is_empty() {
            writeln!(
                writer,
                "{mode}\t{family}\t.\t0\t0\tfalse\tfalse\t0\t0\t{}\t{}\t.\t{}\t{}",
                reads.available,
                reads.total_pairs,
                reads.unique_mapped_pairs,
                reads.distinct_cross_contig_pairs
            )?;
        }
        for c in chains {
            let mut supported = 0;
            let mut discordant = 0;
            for edge in c.models.windows(2) {
                let evidence = reads.edge(&sources[edge[0]], &sources[edge[1]]);
                supported += usize::from(evidence.concordant >= 2);
                discordant += usize::from(evidence.warning());
            }
            let members = c
                .models
                .iter()
                .map(|&i| sources[i].candidate_id.as_str())
                .collect::<Vec<_>>()
                .join(";");
            writeln!(writer, "{mode}\t{family}\t{}\t{}\t{}\t{}\t{}\t{supported}\t{discordant}\t{}\t{}\t{members}\t{}\t{}", c.protein, c.models.len(), c.covered_aa, c.ambiguous, c.search_limited, reads.available, reads.total_pairs, reads.unique_mapped_pairs, reads.distinct_cross_contig_pairs)?;
        }
    }
    writer.flush()
}

pub(super) fn write_masks(
    out: &Path,
    family: &str,
    sample: &str,
    recovered: &exon_recovery::ScoSequence,
) -> io::Result<()> {
    let directory = out.join("sco_masks");
    fs::create_dir_all(&directory)?;
    let path = directory.join(format!("{family}.tsv"));
    let mut writer = fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(path)?;
    if writer.metadata()?.len() == 0 {
        writeln!(writer, "sample\tfamily_id\tstart_0\tend_0\treason")?;
    }
    for interval in &recovered.mask_intervals {
        writeln!(
            writer,
            "{sample}\t{family}\t{}\t{}\t{}",
            interval.start, interval.end, interval.reason
        )?;
    }
    Ok(())
}

pub(super) fn row(sample: &str, family: &str, recovered: &exon_recovery::ScoSequence) -> String {
    format!("{sample}\t{family}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\treference_coordinate_composite\t{}\t{}\t{}\t{}\n",recovered.reference,recovered.sequence.len(),recovered.sequence.len()-recovered.missing_nt-recovered.conflict_nt,recovered.missing_nt,recovered.conflict_nt,recovered.projected_models,recovered.unmappable_models,recovered.cross_reference_models,recovered.coordinate_conflict_codons,recovered.reference_selection,recovered.dominated_models)
}

pub(super) const fn header() -> &'static str {
    HEADER
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn replay_uses_one_cohort_backbone_without_running_an_aligner() {
        let root = super::super::tests::temporary_test_directory("sco_cohort_replay");
        let input = root.join("input");
        let cached = root.join("cached");
        fs::create_dir_all(cached.join("manifest/derived_proteins")).unwrap();
        fs::create_dir_all(cached.join("sco")).unwrap();
        fs::write(cached.join("sco/fam.fasta"), ">old\nNNN\n").unwrap();
        fs::write(
            cached.join("manifest/derived_proteins/fam.faa"),
            ">p1\nAAAAAAAAAA\n>p2\nAAAAAAAAAA\n",
        )
        .unwrap();
        let mut ids =
            String::from("sample\tfamily_id\tinternal_candidate\tcandidate\tsequence_length\n");
        for (sample, protein, start, end) in [("A", "p1", 0, 4), ("B", "p2", 4, 10)] {
            let dna = "GCT".repeat(end - start);
            fs::create_dir_all(input.join(format!("samples/{sample}/candidates"))).unwrap();
            fs::write(
                input.join(format!("samples/{sample}/candidates/fam.fasta")),
                format!(">{sample}|fam|candidate_1\n{dna}\n"),
            )
            .unwrap();
            let raw = cached.join(format!("raw_exonerate/{sample}/fam"));
            fs::create_dir_all(&raw).unwrap();
            fs::write(
                raw.join("TSK000001.vulgar"),
                format!(
                    "vulgar: {protein} {start} {end} . TSK000001 0 {} + 100 M {} {}\n",
                    dna.len(),
                    end - start,
                    dna.len()
                ),
            )
            .unwrap();
            ids.push_str(&format!(
                "{sample}\tfam\tTSK000001\tcandidate_1\t{}\n",
                dna.len()
            ));
        }
        fs::write(cached.join("manifest/id_map.tsv"), ids).unwrap();
        fs::write(cached.join("manifest/reference_translation.tsv"),"family_id\treference_sequence\tprotein_id\tsource\tstrand\tframe\tnucleotide_length\tamino_acid_length\tinternal_stops\tambiguous_amino_acids\tstarts_with_methionine\tterminal_stop\tselection_score\tselected\tstatus\n").unwrap();
        for workers in [1, 3] {
            let out = root.join(format!("recovery{workers}"));
            recover(
                &cached,
                &input,
                &out,
                workers,
                None,
                None,
                ModelConfig::default(),
            )
            .unwrap();
            let fasta = read_raw_fasta(&out.join("sco/fam.fasta")).unwrap();
            assert_eq!(fasta.len(), 2);
            assert!(fasta
                .iter()
                .all(|(h, s)| h.contains("reference=p2") && s.len() == 30));
            assert_eq!(fasta[0].1, format!("{}{}", "GCT".repeat(4), "N".repeat(18)));
            assert_eq!(fasta[1].1, format!("{}{}", "N".repeat(12), "GCT".repeat(6)));
        }
        assert_eq!(
            fs::read(root.join("recovery1/manifest/sco_recovery.tsv")).unwrap(),
            fs::read(root.join("recovery3/manifest/sco_recovery.tsv")).unwrap()
        );
        let audited = root.join("read_audited");
        recover_with_reads(
            &cached,
            &input,
            &audited,
            2,
            None,
            None,
            ModelConfig::default(),
            Some(&root),
        )
        .unwrap();
        assert_eq!(
            fs::read(audited.join("sco/fam.fasta")).unwrap(),
            fs::read(root.join("recovery1/sco/fam.fasta")).unwrap()
        );
        for sample in ["A", "B"] {
            let text =
                fs::read_to_string(audited.join(format!("chain_audit/{sample}/fam.tsv"))).unwrap();
            assert_eq!(text.lines().count(), 3);
            assert!(text.contains("reads_available"));
            assert!(text.contains("false\t0\t.\t0\t0"));
        }
        let frozen = root.join("frozen.tsv");
        fs::write(&frozen, "family_id\tprotein_reference\nfam\tp1\n").unwrap();
        let out = root.join("frozen");
        recover(
            &cached,
            &input,
            &out,
            2,
            None,
            Some(&frozen),
            ModelConfig::default(),
        )
        .unwrap();
        assert!(read_raw_fasta(&out.join("sco/fam.fasta"))
            .unwrap()
            .iter()
            .all(|(h, _)| h.contains("reference=p1")));
        fs::remove_file(cached.join("raw_exonerate/A/fam/TSK000001.vulgar")).unwrap();
        assert!(recover(
            &cached,
            &input,
            &root.join("missing"),
            1,
            None,
            None,
            ModelConfig::default()
        )
        .is_err());
        fs::remove_dir_all(root).unwrap();
    }
}
