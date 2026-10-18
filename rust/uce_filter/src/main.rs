use std::env;
use std::path::PathBuf;
use std::process;
use uce_filter_core::{run, Config};

fn next_value(argv: &[String], i: &mut usize, key: &str) -> Result<String, String> {
    *i += 1;
    argv.get(*i)
        .cloned()
        .ok_or_else(|| format!("{key} requires a value"))
}

fn parse(argv: Vec<String>) -> Result<Config, String> {
    if argv.iter().any(|v| v == "-h" || v == "--help") {
        println!(
            "uce_filter -r REFERENCES -q1 R1 -q2 R2 -o SAMPLE_DIR [--kmer-size 31 --step 4]\n\
Fused UCE rolling-kmer recruitment, run-k verification and adaptive per-locus selection.\n\
Recruitment and verification may be tuned independently with --verification-kmer-size.\n\
Use --max-locus-count 1 to reject fragments supported by more than one locus.\n\
Use --retain-loci-file FILE to emit only named loci after the ambiguity check.\n\
Use --minimum-alignment-overlap N and --minimum-alignment-identity F to require reference alignment support.\n\
Use --selection legacy for regression, or --reference-role contig during rescue.\n\
Optional evidence-only mode: --alignment-shadow [--shadow-per-locus 64 --shadow-band 32 --terminal-window 150].\n\
Use --profile to print read/decode, recruitment, evidence and candidate-store timings."
        );
        process::exit(0);
    }
    let mut references = None;
    let mut read1 = None;
    let mut read2 = None;
    let mut output = None;
    let mut config = Config {
        references: PathBuf::new(),
        recruit_references: None,
        read1: PathBuf::new(),
        read2: PathBuf::new(),
        output: PathBuf::new(),
        kmer_size: 31,
        verify_kmer_size: None,
        step: 4,
        max_locus_count: 0,
        retain_loci_file: None,
        minimum_alignment_overlap: 0,
        minimum_alignment_identity: 0.0,
        min_depth: 50,
        max_depth: 768,
        max_size_mb: 6,
        max_fragments: 0,
        memory_limit_mib: 512,
        selection_auto: true,
        selection_all: false,
        reference_is_contig: false,
        alignment_shadow: false,
        shadow_per_locus: 64,
        shadow_band: 32,
        terminal_window: 150,
        profile: false,
    };
    let mut i = 0;
    while i < argv.len() {
        let key = &argv[i];
        match key.as_str() {
            "-r" | "--references" => {
                references = Some(PathBuf::from(next_value(&argv, &mut i, key)?))
            }
            "--recruit-references" => {
                config.recruit_references = Some(PathBuf::from(next_value(&argv, &mut i, key)?))
            }
            "-q1" => read1 = Some(PathBuf::from(next_value(&argv, &mut i, key)?)),
            "-q2" => read2 = Some(PathBuf::from(next_value(&argv, &mut i, key)?)),
            "-o" | "--output" => output = Some(PathBuf::from(next_value(&argv, &mut i, key)?)),
            "-kf" | "--kmer-size" => {
                config.kmer_size = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--verification-kmer-size" => {
                config.verify_kmer_size = Some(
                    next_value(&argv, &mut i, key)?
                        .parse()
                        .map_err(|_| format!("invalid {key}"))?,
                )
            }
            "-s" | "--step" => {
                config.step = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--max-locus-count" => {
                config.max_locus_count = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--retain-loci-file" => {
                config.retain_loci_file = Some(PathBuf::from(next_value(&argv, &mut i, key)?))
            }
            "--minimum-alignment-overlap" => {
                config.minimum_alignment_overlap = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--minimum-alignment-identity" => {
                config.minimum_alignment_identity = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--min-depth" => {
                config.min_depth = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--max-depth" => {
                config.max_depth = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--max-size" => {
                config.max_size_mb = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--max-fragments" | "-m_reads" => {
                config.max_fragments = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--memory-limit-mib" => {
                config.memory_limit_mib = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--alignment-shadow" => config.alignment_shadow = true,
            "--profile" => config.profile = true,
            "--shadow-per-locus" => {
                config.shadow_per_locus = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--shadow-band" => {
                config.shadow_band = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--terminal-window" => {
                config.terminal_window = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?
            }
            "--selection" => {
                let value = next_value(&argv, &mut i, key)?;
                config.selection_auto = match value.as_str() {
                    "auto" => true,
                    "legacy" | "all" => false,
                    _ => return Err("--selection must be auto, legacy or all".to_string()),
                };
                config.selection_all = value == "all";
            }
            "--reference-role" => {
                let value = next_value(&argv, &mut i, key)?;
                config.reference_is_contig = match value.as_str() {
                    "bait" => false,
                    "contig" => true,
                    _ => return Err("--reference-role must be bait or contig".to_string()),
                };
            }
            "--threads" => {
                let value: usize = next_value(&argv, &mut i, key)?
                    .parse()
                    .map_err(|_| format!("invalid {key}"))?;
                if value != 1 {
                    return Err(
                        "UCEFilter uses one compute thread plus two bounded background decode workers per sample".to_string(),
                    );
                }
            }
            _ => return Err(format!("unknown argument {key}")),
        }
        i += 1;
    }
    config.references = references.ok_or_else(|| "-r is required".to_string())?;
    config.read1 = read1.ok_or_else(|| "-q1 is required".to_string())?;
    config.read2 = read2.ok_or_else(|| "-q2 is required".to_string())?;
    config.output = output.ok_or_else(|| "-o is required".to_string())?;
    if config.step == 0 {
        return Err("--step must be positive".to_string());
    }
    if config
        .verify_kmer_size
        .is_some_and(|verify_k| !(1..=64).contains(&verify_k))
    {
        return Err("--verification-kmer-size must be in 1..=64".to_string());
    }
    if config.min_depth < 0 {
        return Err("--min-depth must be non-negative".to_string());
    }
    if config.max_depth <= 0 {
        return Err("--max-depth must be positive".to_string());
    }
    if config.max_size_mb <= 0 {
        return Err("--max-size must be positive".to_string());
    }
    if config.alignment_shadow && config.shadow_per_locus == 0 {
        return Err("--shadow-per-locus must be positive".to_string());
    }
    if config.alignment_shadow && config.shadow_band == 0 {
        return Err("--shadow-band must be positive".to_string());
    }
    if !config.minimum_alignment_identity.is_finite()
        || !(0.0..=1.0).contains(&config.minimum_alignment_identity)
    {
        return Err("--minimum-alignment-identity must be in 0..=1".to_string());
    }
    if (config.minimum_alignment_overlap > 0 || config.minimum_alignment_identity > 0.0)
        && config.shadow_band == 0
    {
        return Err(
            "--shadow-band must be positive when alignment filtering is enabled".to_string(),
        );
    }
    Ok(config)
}

fn main() {
    let result = parse(env::args().skip(1).collect()).and_then(|config| {
        let profile = config.profile;
        let summary = run(&config)?;
        Ok((summary, profile))
    });
    match result {
        Ok((summary, profile)) => {
            eprintln!(
                "UCEFilter finished: {} fragments read, {} retained once, {} assignments, {} loci, {:.1} MiB fragments + {:.1} MiB candidates + {:.1} MiB spill ({:.3}s)",
                summary.fragments_read, summary.fragments_retained_once, summary.assignments,
                summary.loci_written,
                summary.fragment_memory_bytes as f64 / 1_048_576.0,
                summary.candidate_memory_bytes as f64 / 1_048_576.0,
                summary.fragment_spill_bytes as f64 / 1_048_576.0,
                summary.elapsed_seconds,
            );
            if profile {
                eprintln!(
                    "UCEFilter stages: index {:.3}s, scan/evidence {:.3}s, selection/routes {:.3}s, output {:.3}s",
                    summary.index_seconds,
                    summary.scan_seconds,
                    summary.selection_seconds,
                    summary.output_seconds,
                );
                eprintln!(
                    "UCEFilter scan profile: FASTQ wait/parse {:.3}s, recruit {:.3}s, run-k/exact {:.3}s, candidate store {:.3}s",
                    summary.decode_seconds,
                    summary.recruit_seconds,
                    summary.evidence_seconds,
                    summary.store_seconds,
                );
                eprintln!(
                    "UCEFilter storage profile: evidence scratch {:.3} MiB, candidate pool {:.3} MiB",
                    summary.evidence_scratch_bytes as f64 / 1_048_576.0,
                    summary.candidate_memory_bytes as f64 / 1_048_576.0,
                );
                let profile = &summary.index_profile;
                eprintln!(
                    "UCEFilter index profile: {} recruit probes / {} Bloom negatives / {} hits, {} exact locus queries / {} per-locus index queries, {} run-k windows / {} matching, {} MEM starts / {} MEM bp",
                    profile.recruit_probes,
                    profile.recruit_bloom_rejected,
                    profile.recruit_hits,
                    profile.exact_locus_queries,
                    profile.exact_index_queries,
                    profile.exact_run_windows,
                    profile.exact_matching_windows,
                    profile.mem_starts,
                    profile.mem_bases,
                );
            }
            if summary.shadow_sampled_assignments > 0 {
                eprintln!(
                    "Alignment shadow: {} sampled locus assignments, {} aligned mates ({:.3}s alignment time)",
                    summary.shadow_sampled_assignments,
                    summary.shadow_aligned_mates,
                    summary.shadow_seconds,
                );
            }
        }
        Err(error) => {
            eprintln!("Error: {error}");
            process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::parse;
    use std::path::PathBuf;

    fn required() -> Vec<String> {
        ["-r", "refs", "-q1", "r1.fq", "-q2", "r2.fq", "-o", "out"]
            .into_iter()
            .map(str::to_string)
            .collect()
    }

    #[test]
    fn copy_recovery_can_disable_abundance_thinning() {
        let mut args = required();
        args.extend(["--selection".into(), "all".into()]);
        let config = parse(args).unwrap();
        assert!(config.selection_all);
        assert!(!config.selection_auto);
    }

    #[test]
    fn rejects_invalid_resource_and_depth_limits() {
        for extra in [
            ["--step", "0"],
            ["--min-depth", "-1"],
            ["--max-depth", "0"],
            ["--max-size", "0"],
            ["--threads", "2"],
            ["--minimum-alignment-identity", "1.1"],
            ["--alignment-shadow", "--shadow-per-locus"],
        ] {
            let mut args = required();
            args.extend(extra.into_iter().map(str::to_string));
            assert!(parse(args).is_err());
        }
    }

    #[test]
    fn accepts_default_off_and_bounded_alignment_shadow() {
        let plain = parse(required()).unwrap();
        assert!(!plain.alignment_shadow);
        assert!(!plain.profile);
        assert!(plain.selection_auto);
        assert!(!plain.reference_is_contig);
        let mut args = required();
        args.extend(
            [
                "--alignment-shadow",
                "--shadow-per-locus",
                "8",
                "--shadow-band",
                "24",
                "--terminal-window",
                "100",
                "--profile",
            ]
            .into_iter()
            .map(str::to_string),
        );
        let shadow = parse(args).unwrap();
        assert!(shadow.alignment_shadow);
        assert_eq!(shadow.shadow_per_locus, 8);
        assert_eq!(shadow.shadow_band, 24);
        assert_eq!(shadow.terminal_window, 100);
        assert!(shadow.profile);
    }

    #[test]
    fn accepts_legacy_regression_and_contig_rescue_roles() {
        let mut args = required();
        args.extend(
            ["--selection", "legacy", "--reference-role", "contig"]
                .into_iter()
                .map(str::to_string),
        );
        let config = parse(args).unwrap();
        assert!(!config.selection_auto);
        assert!(config.reference_is_contig);
    }

    #[test]
    fn accepts_independent_verification_and_unique_locus_filter() {
        let mut args = required();
        args.extend(
            [
                "--kmer-size",
                "21",
                "--verification-kmer-size",
                "19",
                "--max-locus-count",
                "1",
                "--retain-loci-file",
                "unresolved.txt",
            ]
            .into_iter()
            .map(str::to_string),
        );
        let config = parse(args).unwrap();
        assert_eq!(config.kmer_size, 21);
        assert_eq!(config.verify_kmer_size, Some(19));
        assert_eq!(config.max_locus_count, 1);
        assert_eq!(config.minimum_alignment_overlap, 0);
        assert_eq!(config.minimum_alignment_identity, 0.0);
        assert_eq!(
            config.retain_loci_file,
            Some(PathBuf::from("unresolved.txt"))
        );

        let mut args = required();
        args.extend(
            [
                "--minimum-alignment-overlap",
                "45",
                "--minimum-alignment-identity",
                "0.8",
            ]
            .into_iter()
            .map(str::to_string),
        );
        let config = parse(args).unwrap();
        assert_eq!(config.minimum_alignment_overlap, 45);
        assert_eq!(config.minimum_alignment_identity, 0.8);
    }
}
