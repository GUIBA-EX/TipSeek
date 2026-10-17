//! Optional, non-selecting Exonerate comparison. Never creates accepted models.
use super::exonerate::{parse_vulgar, Hit};
use super::{Interval, ModelConfig, RawGeneModel, Strand};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

fn field(value: &str) -> String {
    value.replace(['\t', '\r', '\n'], " ")
}
fn intervals(values: &[Interval]) -> String {
    values
        .iter()
        .map(|v| format!("{}-{}", v.start, v.end))
        .collect::<Vec<_>>()
        .join(",")
}
fn same_coordinates(hit: &Hit, model: &RawGeneModel) -> bool {
    let mut cdss: Vec<_> = model.cdss.iter().map(|cds| cds.interval).collect();
    cdss.sort_by_key(|v| v.start);
    model.protein_id == hit.protein
        && model.strand == hit.strand
        && model.target_interval_aa == hit.query
        && model.genomic_interval == hit.genomic
        && cdss == hit.cdss
}

pub fn execute(
    executable: &Path,
    args: &[String],
    directory: &Path,
    timeout: Duration,
) -> io::Result<String> {
    fs::create_dir_all(directory)?;
    // One argument per line, passed directly to Command (never through a shell).
    fs::write(
        directory.join("arguments.txt"),
        format!("{}\n{}\n", executable.display(), args.join("\n")),
    )?;
    let mut command = Command::new(executable);
    command
        .args(args)
        .stdin(Stdio::null())
        .stdout(File::create(directory.join("stdout.txt"))?)
        .stderr(File::create(directory.join("stderr.txt"))?);
    let mut child = match command.spawn() {
        Ok(child) => child,
        Err(error) => {
            fs::write(directory.join("spawn_error.txt"), error.to_string())?;
            fs::write(directory.join("status.txt"), "spawn_failed\n")?;
            return Ok("spawn_failed".into());
        }
    };
    // Linux can wait for process exit without timer polling. Keep the portable
    // fallback for older kernels and systems where pidfd_open is unavailable.
    #[cfg(target_os = "linux")]
    let exit_fd = {
        use std::os::fd::{FromRawFd, OwnedFd};
        // SAFETY: pidfd_open takes a live child pid and zero flags; ownership of
        // the returned descriptor is transferred exactly once to OwnedFd.
        let fd = unsafe { libc::syscall(libc::SYS_pidfd_open, child.id(), 0) };
        (fd >= 0).then(|| unsafe { OwnedFd::from_raw_fd(fd as i32) })
    };
    let start = Instant::now();
    let mut poll_delay = Duration::from_millis(1);
    let status = loop {
        if let Some(status) = child.try_wait()? {
            break if status.success() {
                "ok".into()
            } else {
                format!("exit_{status}")
            };
        }
        if start.elapsed() >= timeout {
            let _ = child.kill();
            child.wait()?;
            break "timeout".into();
        }
        #[cfg(target_os = "linux")]
        if let Some(fd) = &exit_fd {
            use std::os::fd::AsRawFd;
            let remaining = timeout.saturating_sub(start.elapsed());
            let millis = remaining
                .as_millis()
                .saturating_add(u128::from(remaining.subsec_nanos() % 1_000_000 != 0))
                .min(i32::MAX as u128) as i32;
            let mut event = libc::pollfd {
                fd: fd.as_raw_fd(),
                events: libc::POLLIN,
                revents: 0,
            };
            // SAFETY: event is a valid, initialized pollfd for this call.
            let ready = unsafe { libc::poll(&mut event, 1, millis) };
            if ready >= 0 || io::Error::last_os_error().kind() == io::ErrorKind::Interrupted {
                continue;
            }
            // An unexpected poll error still uses the bounded portable wait.
        }
        // Bound completion latency for short alignments. A 2 ms cap avoids
        // the 10–25 ms polling steps dominating these small subprocesses.
        std::thread::sleep(poll_delay.min(timeout.saturating_sub(start.elapsed())));
        poll_delay = (poll_delay * 2).min(Duration::from_millis(2));
    };
    fs::write(
        directory.join("status.txt"),
        format!("{status}\nseconds={}\n", start.elapsed().as_secs_f64()),
    )?;
    Ok(status)
}

pub struct ShadowAudit {
    executable: PathBuf,
    root: PathBuf,
    writer: BufWriter<File>,
}

impl ShadowAudit {
    pub fn new(executable: PathBuf, out: &Path) -> io::Result<Self> {
        let root = out.join("raw_exonerate_shadow");
        fs::create_dir_all(&root)?;
        let _ = execute(
            &executable,
            &["--version".into()],
            &root.join("version"),
            Duration::from_secs(10),
        )?;
        let mut writer = BufWriter::new(File::create(out.join("manifest/exonerate_shadow.tsv"))?);
        writeln!(writer, "sample\tfamily\tcandidate\tstage\tprotein\tstatus\tminiprot_status\tminiprot_model_count\texonerate_hit_count\tcoordinate_match_count\texonerate_score\tquery_start_aa\tquery_end_aa\ttarget_start_nt\ttarget_end_nt\tstrand\tcds_intervals_nt\tframeshift_operations")?;
        Ok(Self {
            executable,
            root,
            writer,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn compare(
        &mut self,
        sample: &str,
        family: &str,
        candidate: &str,
        internal: &str,
        stage: &str,
        target: &Path,
        target_sequence: &str,
        proteins: &Path,
        protein_sequences: &BTreeMap<String, String>,
        config: ModelConfig,
        max_intron: usize,
        miniprot_output: &str,
        miniprot_ok: bool,
    ) -> io::Result<()> {
        let directory = self.root.join(sample).join(family).join(internal);
        fs::create_dir_all(&directory)?;
        let saved_target = directory.join("target.fasta");
        let saved_proteins = directory.join("proteins.faa");
        fs::copy(target, &saved_target)?;
        fs::copy(proteins, &saved_proteins)?;
        let args = command_arguments(&saved_proteins, &saved_target, max_intron);
        let run_status = execute(&self.executable, &args, &directory, Duration::from_secs(60))?;
        let parsed = super::parse_miniprot_output(miniprot_output);
        let mini_status = if !miniprot_ok {
            "failed"
        } else if !parsed.warnings.is_empty() || parsed.unassigned_paf > 0 {
            "parse_warning"
        } else {
            "ok"
        };
        let protein_lengths: BTreeMap<String, usize> = protein_sequences
            .iter()
            .map(|(id, seq)| (id.clone(), seq.len()))
            .collect();
        let result = if run_status == "ok" {
            fs::read_to_string(directory.join("stdout.txt"))
                .map_err(|error| format!("invalid_output: {error}"))
                .and_then(|text| {
                    parse_vulgar(&text, internal, target_sequence.len(), &protein_lengths)
                })
        } else {
            Err(run_status)
        };
        if let Err(error) = &result {
            fs::write(directory.join("comparison_error.txt"), error)?;
        }
        // Structural audit only: never feed these models into selection or padding.
        let mut normalized = BufWriter::new(File::create(directory.join("normalized_models.tsv"))?);
        writeln!(normalized, "engine\tprotein\tstatus\tstate\tstructurally_eligible\tidentity_millionths\tidentity_basis\tpositive_millionths\tphase_valid\tcds_phases\tframeshifts\tlarge_query_gaps\tcds\tprotein_sequence\tqc_flags")?;
        if let Ok(hits) = &result {
            for hit in hits {
                let result = hit
                    .normalize(internal, target_sequence, &protein_sequences[&hit.protein])
                    .and_then(|adapted| adapted.build(target_sequence, config));
                match result {
                    Ok(model) => {
                        let phases = model
                            .exons
                            .iter()
                            .map(|e| e.phase.map_or(".".into(), |p| p.to_string()))
                            .collect::<Vec<_>>()
                            .join(",");
                        writeln!(normalized, "exonerate\t{}\tok\t{}\t{}\t{}\taligned_M_S_pairs\t.\t{}\t{}\t{}\t{}\t{}\t{}\t{}", field(&hit.protein), model.state.as_str(), model.eligible_for_resolve, model.identity_millionths, model.phase_valid, phases, model.frameshifts, model.large_query_gaps, model.cds, model.protein.as_deref().unwrap_or("."), model.qc_flags.join(","))?;
                    }
                    Err(error) => writeln!(
                        normalized,
                        "exonerate\t{}\terror\t.\tfalse\t.\t.\t.\t.\t.\t.\t.\t.\t.\t{}",
                        field(&hit.protein),
                        field(&error)
                    )?,
                }
            }
        }
        normalized.flush()?;
        for protein in protein_lengths.keys() {
            let mini: Vec<_> = parsed
                .models
                .iter()
                .filter(|m| m.candidate_id == internal && m.protein_id == *protein)
                .collect();
            let hits: Vec<_> = result
                .as_ref()
                .map(|h| h.iter().filter(|h| h.protein == *protein).collect())
                .unwrap_or_default();
            let matched = hits
                .iter()
                .filter(|h| mini.iter().any(|m| same_coordinates(h, m)))
                .count();
            let status = if result.is_err() {
                "exonerate_failed"
            } else if mini_status != "ok" {
                "miniprot_unavailable"
            } else if hits.is_empty() && mini.is_empty() {
                "both_no_hit"
            } else if hits.is_empty() {
                "miniprot_only"
            } else if mini.is_empty() {
                "exonerate_only"
            } else if matched == hits.len() {
                "reported_hits_match"
            } else {
                "coordinate_difference"
            };
            let prefix = format!(
                "{}\t{}\t{}\t{}\t{}\t{status}\t{mini_status}\t{}\t{}\t{matched}",
                field(sample),
                field(family),
                field(candidate),
                field(stage),
                field(protein),
                mini.len(),
                hits.len()
            );
            if hits.is_empty() {
                writeln!(self.writer, "{prefix}\t.\t.\t.\t.\t.\t.\t.\t.")?;
            }
            for hit in hits {
                writeln!(
                    self.writer,
                    "{prefix}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    hit.score,
                    hit.query.start,
                    hit.query.end,
                    hit.genomic.start,
                    hit.genomic.end,
                    if hit.strand == Strand::Forward {
                        "+"
                    } else {
                        "-"
                    },
                    intervals(&hit.cdss),
                    hit.frameshift_ops
                )?;
            }
        }
        self.writer.flush()
    }
}

/// Shared protein2genome invocation for primary and shadow modes.
pub fn command_arguments(proteins: &Path, target: &Path, max_intron: usize) -> Vec<String> {
    [
        "--model",
        "protein2genome",
        "--querytype",
        "protein",
        "--targettype",
        "dna",
        "--exhaustive",
        "no",
        "--subopt",
        "no",
        "--revcomp",
        "yes",
        "--score",
        "1",
        "--dpmemory",
        "64",
        "--verbose",
        "0",
        "--showalignment",
        "no",
        "--showsugar",
        "no",
        "--showcigar",
        "no",
        "--showvulgar",
        "yes",
        "--maxintron",
    ]
    .into_iter()
    .map(str::to_owned)
    .chain([
        max_intron.to_string(),
        proteins.display().to_string(),
        target.display().to_string(),
    ])
    .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn split_codons_and_reverse_interbase_coordinates() {
        let proteins = BTreeMap::from([("p".into(), 3)]);
        for (bounds, expected) in [
            ("0 17 +", vec![(0, 4), (12, 17)]),
            ("17 0 -", vec![(0, 5), (13, 17)]),
        ] {
            let hit = parse_vulgar(
                &format!("vulgar: p 0 3 . c {bounds} 30 M 1 3 S 0 1 5 0 2 I 0 4 3 0 2 S 1 2 M 1 3"),
                "c",
                17,
                &proteins,
            )
            .unwrap()
            .remove(0);
            assert_eq!(
                hit.cdss
                    .iter()
                    .map(|v| (v.start, v.end))
                    .collect::<Vec<_>>(),
                expected
            );
        }
        assert!(parse_vulgar("vulgar: p 0 3 . c 0 17 + 30 M 3 9", "c", 17, &proteins).is_err());
        assert!(parse_vulgar("vulgar: p 0 3 . wrong 0 9 + 30 M 3 9", "c", 17, &proteins).is_err());
    }
    #[cfg(unix)]
    #[test]
    fn timeouts_and_failed_processes_are_distinct_from_no_hits() {
        use std::os::unix::fs::PermissionsExt;
        let root = std::env::temp_dir().join(format!(
            "tipseek-shadow-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let executable = root.join("slow");
        fs::write(&executable, "#!/bin/sh\nexec sleep 2\n").unwrap();
        fs::set_permissions(&executable, fs::Permissions::from_mode(0o755)).unwrap();
        assert_eq!(
            execute(
                &executable,
                &[],
                &root.join("timeout"),
                Duration::from_millis(20)
            )
            .unwrap(),
            "timeout"
        );
        assert_eq!(
            execute(
                &root.join("missing"),
                &[],
                &root.join("missing-run"),
                Duration::from_secs(1)
            )
            .unwrap(),
            "spawn_failed"
        );
        fs::write(&executable, "#!/bin/sh\nexit 3\n").unwrap();
        assert!(execute(
            &executable,
            &[],
            &root.join("failed"),
            Duration::from_secs(1)
        )
        .unwrap()
        .starts_with("exit_"));
        assert!(parse_vulgar("", "c", 17, &BTreeMap::new())
            .unwrap()
            .is_empty());
        assert!(parse_vulgar("unexpected output", "c", 17, &BTreeMap::new()).is_err());
        fs::remove_dir_all(root).unwrap();
    }
}
