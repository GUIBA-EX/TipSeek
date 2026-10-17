//! GeneMiner2-style queue: each worker handles one complete SCO, then immediately
//! claims another. Workers own isolated files; merge in sample/family order.
use super::*;
use std::cell::RefCell;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;

thread_local! {
    static FAMILY_BUDGET: RefCell<Option<(usize, Arc<AtomicUsize>)>> = const { RefCell::new(None) };
}

pub(super) fn candidate_width(fallback: usize) -> usize {
    FAMILY_BUDGET.with(|slot| {
        slot.borrow()
            .as_ref()
            .map_or(fallback, |(limit, remaining)| {
                // Include queued jobs: borrowing starts only as the pool drains.
                // Existing batches have allocations from a larger remaining count.
                (limit / remaining.load(Ordering::Relaxed).max(1)).max(1)
            })
    })
}

struct Job {
    sample: String,
    family: String,
    input: PathBuf,
    references: PathBuf,
    output: PathBuf,
}

fn link_or_copy(source: &Path, target: &Path) -> io::Result<()> {
    fs::hard_link(source, target).or_else(|_| fs::copy(source, target).map(|_| ()))
}

fn merge_tree(
    source: &Path,
    target: &Path,
    relative: &Path,
    seen: &mut HashMap<PathBuf, HashSet<String>>,
) -> io::Result<()> {
    fs::create_dir_all(target)?;
    for path in sorted_directory_paths(source)? {
        let name = path.file_name().unwrap();
        let rel = relative.join(name);
        let destination = target.join(name);
        if path.is_dir() {
            merge_tree(&path, &destination, &rel, seen)?;
        } else if rel == Path::new("manifest/execution.tsv") {
            continue;
        } else if path.extension().is_some_and(|x| x == "tsv") {
            let mut reader = BufReader::new(File::open(&path)?);
            let mut header = String::new();
            reader.read_line(&mut header)?;
            let deduplicate = matches!(
                rel.to_str(),
                Some("manifest/reference_translation.tsv" | "manifest/backend.tsv")
            );
            let mut writer = fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(&destination)?;
            if destination.metadata()?.len() == 0 {
                writer.write_all(header.as_bytes())?;
            } else {
                let mut existing = String::new();
                BufReader::new(File::open(&destination)?).read_line(&mut existing)?;
                if existing != header {
                    return Err(io::Error::other(format!(
                        "inconsistent worker table header: {}",
                        rel.display()
                    )));
                }
            }
            for line in reader.lines() {
                let line = line?;
                if !deduplicate || seen.entry(rel.clone()).or_default().insert(line.clone()) {
                    writeln!(writer, "{line}")?;
                }
            }
        } else if !destination.exists() {
            fs::rename(&path, &destination)?;
        } else if rel.starts_with("manifest/derived_proteins")
            || path.extension().is_none_or(|x| x != "fasta")
        {
            if fs::read(&path)? != fs::read(&destination)? {
                return Err(io::Error::other(format!(
                    "conflicting worker output: {}",
                    rel.display()
                )));
            }
        } else {
            let mut file = fs::OpenOptions::new().append(true).open(&destination)?;
            io::copy(&mut File::open(&path)?, &mut file)?;
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn annotate(
    input: &Path,
    references: &Path,
    proteins: Option<&Path>,
    out: &Path,
    backend: &AnnotationBackend,
    workers: usize,
    max_intron: usize,
    minimum_coverage: f64,
    complete_coverage: f64,
    flank: usize,
    fragment_padding: usize,
    shadow: Option<&Path>,
    reads: Option<&Path>,
    sco_reference: Option<&str>,
) -> io::Result<()> {
    let reference_files: BTreeMap<_, _> = sorted_directory_paths(references)?
        .into_iter()
        .filter_map(|p| family_id(&p).map(|id| (id, p)))
        .collect();
    if out.exists() {
        fs::remove_dir_all(out)?;
    }
    let temporary = out.join(".family_workers");
    fs::create_dir_all(&temporary)?;
    let mut jobs = Vec::new();
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
        // Parse the sample-wide assembly audit once. Copying it into every
        // task would make 7,000 workers repeatedly scan the entire table.
        let assembly = sample_path.join("candidate_assembly.tsv");
        let mut assembly_header = None;
        let mut assembly_rows: BTreeMap<String, Vec<csv::StringRecord>> = BTreeMap::new();
        if assembly.is_file() {
            let mut reader = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .from_path(&assembly)?;
            let header = reader.headers()?.clone();
            if let Some(column) = header.iter().position(|name| name == "family_id") {
                for row in reader.records() {
                    let row = row?;
                    assembly_rows
                        .entry(row[column].to_owned())
                        .or_default()
                        .push(row);
                }
            }
            assembly_header = Some(header);
        }
        for path in sorted_directory_paths(&candidates)? {
            let Some(family) = family_id(&path) else {
                continue;
            };
            let task = temporary.join(format!("task_{:06}", jobs.len()));
            let task_input = task.join("input");
            let task_sample = task_input.join("samples").join(&sample);
            fs::create_dir_all(task_sample.join("candidates"))?;
            link_or_copy(
                &path,
                &task_sample
                    .join("candidates")
                    .join(path.file_name().unwrap()),
            )?;
            if let Some(header) = &assembly_header {
                let mut writer = csv::WriterBuilder::new()
                    .delimiter(b'\t')
                    .from_path(task_sample.join("candidate_assembly.tsv"))?;
                writer.write_record(header)?;
                for row in assembly_rows.remove(&family).unwrap_or_default() {
                    writer.write_record(&row)?;
                }
                writer.flush()?;
            }
            let task_references = task.join("references");
            fs::create_dir_all(&task_references)?;
            if let Some(reference) = reference_files.get(&family) {
                link_or_copy(
                    reference,
                    &task_references.join(reference.file_name().unwrap()),
                )?;
            }
            jobs.push(Job {
                sample: sample.clone(),
                family,
                input: task_input,
                references: task_references,
                output: task.join("output"),
            });
        }
    }
    let next = AtomicUsize::new(0);
    let stopped = AtomicBool::new(false);
    let completed = AtomicUsize::new(0);
    let remaining = Arc::new(AtomicUsize::new(jobs.len()));
    std::thread::scope(|scope| -> io::Result<()> {
        let handles: Vec<_> = (0..workers.min(jobs.len()))
            .map(|_| {
                scope.spawn(|| -> io::Result<()> {
                    while !stopped.load(Ordering::Relaxed) {
                        let index = next.fetch_add(1, Ordering::Relaxed);
                        let Some(job) = jobs.get(index) else { break };
                        FAMILY_BUDGET
                            .with(|slot| *slot.borrow_mut() = Some((workers, remaining.clone())));
                        let result = annotate_with_reads(
                            job.input.clone(),
                            job.references.clone(),
                            proteins.map(Path::to_path_buf),
                            job.output.clone(),
                            backend.clone(),
                            "1".into(),
                            max_intron,
                            minimum_coverage,
                            complete_coverage,
                            flank,
                            fragment_padding,
                            shadow.map(Path::to_path_buf),
                            reads.map(Path::to_path_buf),
                            sco_reference.map(str::to_owned),
                        );
                        FAMILY_BUDGET.with(|slot| slot.borrow_mut().take());
                        if let Err(error) = result {
                            stopped.store(true, Ordering::Relaxed);
                            return Err(io::Error::other(format!(
                                "SCO worker failed ({}/{}): {error}",
                                job.sample, job.family
                            )));
                        }
                        let count = completed.fetch_add(1, Ordering::Relaxed) + 1;
                        remaining.fetch_sub(1, Ordering::Relaxed);
                        if count % 100 == 0 || count == jobs.len() {
                            eprintln!("SCO annotation completed {count}/{}", jobs.len());
                        }
                    }
                    Ok(())
                })
            })
            .collect();
        let mut failure = None;
        for handle in handles {
            let result = handle
                .join()
                .map_err(|_| io::Error::other("SCO annotation worker panicked"))
                .and_then(|r| r);
            if let Err(error) = result {
                failure.get_or_insert(error);
            }
        }
        failure.map_or(Ok(()), Err)
    })?;
    let mut seen = HashMap::new();
    for job in &jobs {
        merge_tree(&job.output, out, Path::new(""), &mut seen)?;
    }
    if jobs.is_empty() {
        fs::remove_dir_all(&temporary)?;
        return annotate_with_reads(
            input.to_path_buf(),
            references.to_path_buf(),
            proteins.map(Path::to_path_buf),
            out.to_path_buf(),
            backend.clone(),
            "1".into(),
            max_intron,
            minimum_coverage,
            complete_coverage,
            flank,
            fragment_padding,
            shadow.map(Path::to_path_buf),
            reads.map(Path::to_path_buf),
            sco_reference.map(str::to_owned),
        );
    }
    fs::write(out.join("manifest/execution.tsv"), format!("candidate_workers\tthreads_per_candidate\tfamily_workers\tscheduling\n{workers}\t1\t{workers}\tgeneminer2_family_queue_shared_budget\n"))?;
    fs::remove_dir_all(&temporary)?;
    if jobs
        .iter()
        .map(|j| &j.sample)
        .collect::<BTreeSet<_>>()
        .len()
        > 1
    {
        sco_rebuild::recover(
            out,
            input,
            out,
            workers,
            sco_reference,
            None,
            ModelConfig {
                minimum_coverage,
                complete_coverage,
            },
        )?;
    }
    Ok(())
}
