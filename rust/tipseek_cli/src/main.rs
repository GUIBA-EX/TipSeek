//! Native Rust command dispatcher.
//!
//! The public CLI is implemented in Rust and does not require a Python runtime.

mod rescue_qc;
mod uce_recruit;

use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;
use std::io::{self, BufRead, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};
use std::sync::{
    atomic::{AtomicBool, Ordering},
    mpsc, Arc, Condvar, Mutex,
};
use std::thread;
use std::time::Instant;
use uce_recruit::RecruitPass;

const UCE_RESCUE_ASSEMBLY_KMER: &str = "21";
const DEFAULT_UCE_FILTER_KMER: &str = "23";
const DEFAULT_UCE_RESCUE_ROUNDS: &str = "1";
const UCE_TERMINAL_MIN_EXTENSION: usize = 30;
const UCE_TERMINAL_MIN_BREADTH: f64 = 0.85;
const UCE_TERMINAL_MAX_GAP: usize = 30;
const UCE_TERMINAL_MIN_FRAGMENTS: usize = 2;
const UCE_TERMINAL_MIN_BRIDGES: usize = 1;

const COMMANDS: &[&str] = &[
    "filter",
    "refilter",
    "assemble",
    "gene",
    "stats",
    "te",
    "population",
    "consensus",
    "trim",
    "combine",
    "tree",
    "gene-annotate",
    "gene-resolve",
    "gene-tree",
    "profiling",
    "mito",
    "rad",
    "rad-probe",
    "rad-validate",
];
const FLAG_OPTIONS: &[&str] = &[
    "--te-estimate-genome-fraction",
    "--rad-linked-recruitment",
    "--uce-alignment-shadow",
    "--uce-rescue-reads",
    "--no-uce-rescue-reads",
    "--stats-count-input-reads",
    "--stats-no-heatmap",
    "--population-panrefv2-include-low-confidence",
    "--population-skip-mark-duplicates",
    "--population-skip-plink",
    "--population-skip-admixture",
    "--strict-combine-errors",
    "--no-alignment",
    "--no-trimal",
    "--profile-force-rebuild",
    "--cleanup-intermediates",
    "--cleanup-dry-run",
    "--no-mito-adaptive-stop",
    "--reuse-reference-cache",
    "--legacy-uce-filter",
    "--workflow-profile",
    "--resume",
    "--rad-denovo",
];
const VALUE_OPTIONS: &[&str] = &[
    "-f",
    "-r",
    "-o",
    "-p",
    "--log-format",
    "--assembly-mode",
    "-kf",
    "-s",
    "--step-size",
    "-ka",
    "--min-ka",
    "--max-ka",
    "-e",
    "--error-threshold",
    "-i",
    "--search-depth",
    "-sb",
    "--soft-boundary",
    "--min-coverage",
    "--depth-low-water-mark",
    "--depth-limit",
    "--file-size-limit",
    "--max-reads",
    "--assembler-graph-format",
    "--uce-side-candidates",
    "--uce-max-contig-length",
    "--uce-min-read-density",
    "--uce-density-check-min-length",
    "--uce-max-depth-cv",
    "--uce-max-depth-ratio",
    "--uce-shadow-per-locus",
    "--uce-shadow-band",
    "--uce-shadow-terminal-window",
    "--uce-recruit-mode",
    "--uce-fallback-kmer-size",
    "--uce-fallback-step",
    "--uce-fallback-verify-kmer-size",
    "--uce-fallback-min-alignment-overlap",
    "--uce-fallback-min-alignment-identity",
    "--te-stage",
    "--te-kmer",
    "--te-min-kmer-count",
    "--te-catalog-pairs",
    "--te-read-ledger",
    "--te-library",
    "--te-quantify-pairs",
    "--te-bootstrap-replicates",
    "--te-annotate-min-fragment",
    "--te-annotate-max-fragment",
    "--te-annotate-min-support",
    "--te-annotate-min-identity",
    "--te-annotate-min-coverage",
    "--te-annotate-min-delta",
    "--te-assemble-min-kmer-count",
    "--te-assemble-branch-ratio",
    "--te-assemble-max-fragments",
    "--engine",
    "--population-reference-strategy",
    "--population-reference-fasta",
    "--population-min-mapq",
    "--population-min-baseq",
    "--population-min-dp",
    "--population-min-gq",
    "--population-min-qual",
    "--population-min-call-rate",
    "--population-min-mac",
    "--population-ld-window",
    "--population-ld-step",
    "--population-ld-r2",
    "--population-admixture-k-min",
    "--population-admixture-k-max",
    "--population-admixture-cv",
    "--population-start-at",
    "--population-stop-after",
    "--population-minibwa",
    "--population-samtools",
    "--population-bcftools",
    "--population-plink",
    "--population-admixture",
    "-c",
    "--consensus-threshold",
    "-ts",
    "--trim-source",
    "-tm",
    "--trim-mode",
    "-tr",
    "--trim-retention",
    "-cs",
    "--combine-source",
    "-cd",
    "--clean-difference",
    "-cn",
    "--clean-sequences",
    "--msa-program",
    "--msa-threads",
    "--alignment-filter",
    "--filter-processes",
    "--alifilter-model",
    "-m",
    "--tree-method",
    "-b",
    "--bootstrap",
    "--phylo-program",
    "--gene-protein-reference",
    "--gene-miniprot",
    "--gene-input",
    "--gene-mafft",
    "--gene-iqtree",
    "--gene-min-taxa",
    "--gene-min-aa-length",
    "--gene-min-effective-codon-sites",
    "--gene-outgroup",
    "--gene-ufboot",
    "--gene-taper",
    "--gene-julia",
    "--gene-species-mode",
    "--gene-aster",
    "--profile-kmer-size",
    "--profile-pseudoalign-threshold",
    "--profile-relevant-kmer-fraction",
    "--profile-group-map",
    "--profile-decoy",
    "--profile-index-dir",
    "--profile-index-memory-gb",
    "--profile-themisto",
    "--reference-cache-dir",
    "--mito-genbank",
    "--mito-flank",
    "--mito-tile-length",
    "--mito-tile-step",
    "--mito-min-overlap",
    "--mito-min-overlap-identity",
    "--mito-min-junction-support",
    "--mito-terminal-window",
    "--mito-link-kmer",
    "--mito-min-link-hits",
    "--mito-min-pair-support",
    "--mito-bridge-kmer",
    "--mito-bridge-min-depth",
    "--mito-max-bridge",
    "--mito-initial-reads",
    "--mito-max-reads",
    "--uce-rescue-rounds",
    "--uce-rescue-min-contig-length",
    "--uce-rescue-terminal-window",
    "--uce-rescue-min-density-ratio",
    "--uce-rescue-reverse-reuse-reference-scale",
    "--uce-rescue-inverted-repeat-min-bp",
    "--assembler-implementation",
    "--assembler-read-chunk-size",
    "--uce-path-strategy",
    "--uce-backbone-lookahead",
    "--min-depth",
    "--max-depth",
    "--ipyrad-loci",
    "--rad-min-arm-breadth",
    "--rad-probe",
    "--ipyrad-params",
    "--ipyrad-executable",
    "--ipyrad-steps",
    "--rad-overhang",
    "--rad-overhang-r2",
    "--rad-kmer",
    "--rad-min-count",
    "--rad-min-samples",
    "--rad-min-length",
    "--rad-max-arm-distance",
    "--rad-fallback-kmers",
    "--rad-link-max-fragments",
    "--rad-recovery",
    "--rad-validate-min-identity",
    "--rad-validate-min-breadth",
    "--rad-validate-min-delta",
];

#[derive(Clone, Debug)]
struct Sample {
    name: String,
    read1: String,
    read2: Option<String>,
}

#[derive(Clone, Debug)]
struct Options {
    raw: Vec<String>,
    commands: Vec<String>,
    reference: String,
    samples: String,
    output: String,
    assembly_mode: String,
    workers: usize,
    worker_source: String,
    kf: String,
    step: String,
    ka: String,
    min_ka: String,
    max_ka: String,
    error_threshold: String,
    search_depth: String,
    soft_boundary: String,
    min_coverage: String,
    low_depth: String,
    depth_limit: String,
    size_limit: String,
    max_reads: String,
    graph_format: String,
    side_candidates: String,
    max_contig_length: String,
    min_density: String,
    density_min_length: String,
    max_depth_cv: String,
    max_depth_ratio: String,
    alignment_shadow: bool,
    shadow_per_locus: String,
    shadow_band: String,
    shadow_terminal_window: String,
    uce_recruit_mode: String,
    uce_fallback_kmer_size: String,
    uce_fallback_step: String,
    uce_fallback_verify_kmer_size: String,
    uce_fallback_min_alignment_overlap: String,
    uce_fallback_min_alignment_identity: String,
    uce_memory_limit_mib: u64,
    rescue: bool,
    stats_count_input_reads: bool,
    stats_no_heatmap: bool,
    cleanup_intermediates: bool,
    cleanup_dry_run: bool,
    reuse_reference_cache: bool,
    legacy_uce_filter: bool,
    workflow_profile: bool,
    resume: bool,
    log_format: String,
}

fn value(args: &[String], names: &[&str], default: &str) -> Result<String, String> {
    for (index, arg) in args.iter().enumerate() {
        if names.contains(&arg.as_str()) {
            return args
                .get(index + 1)
                .cloned()
                .ok_or_else(|| format!("{arg} requires a value"));
        }
        for name in names {
            if let Some(value) = arg.strip_prefix(&format!("{name}=")) {
                return Ok(value.to_owned());
            }
        }
    }
    Ok(default.to_owned())
}

fn resolve_worker_budget(request: &str) -> Result<(usize, String), String> {
    if request == "auto" {
        return Ok(auto_worker_budget());
    }
    let workers = request
        .parse::<usize>()
        .map_err(|_| "-p must be 'auto' or a positive integer")?;
    if workers == 0 {
        return Err("-p must be at least 1".into());
    }
    Ok((workers, format!("explicit -p {workers}")))
}

fn flag(args: &[String], name: &str) -> Result<bool, String> {
    for arg in args {
        if arg == name {
            return Ok(true);
        }
        if arg
            .strip_prefix(name)
            .is_some_and(|suffix| suffix.starts_with('='))
        {
            return Err(format!("{name} does not take a value"));
        }
    }
    Ok(false)
}

fn commands(args: &[String]) -> Result<Vec<String>, String> {
    let mut commands = Vec::new();
    let mut index = 0;
    while index < args.len() {
        let arg = &args[index];
        if COMMANDS.contains(&arg.as_str()) {
            commands.push(arg.clone());
            index += 1;
        } else if arg.starts_with('-') {
            let option = arg.split_once('=').map_or(arg.as_str(), |(name, _)| name);
            if FLAG_OPTIONS.contains(&option) {
                if arg.contains('=') {
                    return Err(format!("{option} does not take a value"));
                }
                index += 1;
            } else if VALUE_OPTIONS.contains(&option) {
                if arg.contains('=') {
                    index += 1;
                } else {
                    if index + 1 == args.len() {
                        return Err(format!("{arg} requires a value"));
                    }
                    index += 2;
                }
            } else {
                return Err(format!("Rust CLI does not support option '{option}'"));
            }
        } else {
            return Err(format!("Rust CLI does not support command '{arg}'"));
        }
    }
    Ok(commands)
}

fn parse(args: &[String]) -> Result<Options, String> {
    let mut commands = commands(args)?;
    let assembly_mode = value(args, &["--assembly-mode"], "original")?;
    if !matches!(assembly_mode.as_str(), "original" | "uce") {
        return Err("--assembly-mode must be original or uce".into());
    }
    let log_format = value(args, &["--log-format"], "text")?;
    if !matches!(log_format.as_str(), "text" | "json") {
        return Err("--log-format must be text or json".into());
    }
    let (workers, worker_source) = resolve_worker_budget(&value(args, &["-p"], "auto")?)?;
    let legacy_uce_filter = flag(args, "--legacy-uce-filter")?;
    let default_uce_recruit_mode = if assembly_mode == "uce" && !legacy_uce_filter {
        "auto"
    } else {
        "fast"
    };
    let uce_recruit_mode = value(args, &["--uce-recruit-mode"], default_uce_recruit_mode)?;
    if !matches!(uce_recruit_mode.as_str(), "fast" | "auto") {
        return Err("--uce-recruit-mode must be fast or auto".into());
    }
    if uce_recruit_mode == "auto" && assembly_mode != "uce" {
        return Err("--uce-recruit-mode auto requires --assembly-mode uce".into());
    }
    let uce_fallback_kmer_size = value(args, &["--uce-fallback-kmer-size"], "21")?;
    let uce_fallback_step = value(args, &["--uce-fallback-step"], "1")?;
    let uce_fallback_verify_kmer_size = value(args, &["--uce-fallback-verify-kmer-size"], "19")?;
    let uce_fallback_min_alignment_overlap =
        value(args, &["--uce-fallback-min-alignment-overlap"], "45")?;
    let uce_fallback_min_alignment_identity =
        value(args, &["--uce-fallback-min-alignment-identity"], "0.80")?;
    for (name, raw) in [
        ("--uce-fallback-kmer-size", &uce_fallback_kmer_size),
        (
            "--uce-fallback-verify-kmer-size",
            &uce_fallback_verify_kmer_size,
        ),
    ] {
        let parsed = raw
            .parse::<usize>()
            .map_err(|_| format!("{name} must be an integer in 1..=64"))?;
        if !(1..=64).contains(&parsed) {
            return Err(format!("{name} must be an integer in 1..=64"));
        }
    }
    if uce_fallback_step
        .parse::<usize>()
        .ok()
        .is_none_or(|step| step == 0)
    {
        return Err("--uce-fallback-step must be a positive integer".into());
    }
    if uce_fallback_min_alignment_overlap.parse::<usize>().is_err() {
        return Err("--uce-fallback-min-alignment-overlap must be a non-negative integer".into());
    }
    let fallback_alignment_identity = uce_fallback_min_alignment_identity
        .parse::<f64>()
        .map_err(|_| "--uce-fallback-min-alignment-identity must be in 0..=1")?;
    if !fallback_alignment_identity.is_finite()
        || !(0.0..=1.0).contains(&fallback_alignment_identity)
    {
        return Err("--uce-fallback-min-alignment-identity must be in 0..=1".into());
    }
    if legacy_uce_filter && uce_recruit_mode == "auto" {
        return Err("--uce-recruit-mode auto is unavailable with --legacy-uce-filter".into());
    }
    let rescue_requested = flag(args, "--uce-rescue-reads")?;
    let rescue_disabled = flag(args, "--no-uce-rescue-reads")?;
    if rescue_requested && rescue_disabled {
        return Err("--uce-rescue-reads and --no-uce-rescue-reads cannot be used together".into());
    }
    let rescue = if rescue_disabled {
        false
    } else {
        rescue_requested || assembly_mode == "uce"
    };
    if commands == ["gene"] {
        commands = vec![
            "filter".into(),
            "refilter".into(),
            "assemble".into(),
            "gene".into(),
        ];
    } else if commands.is_empty() {
        commands = if assembly_mode == "uce" {
            vec![
                "filter".into(),
                "refilter".into(),
                "assemble".into(),
                "combine".into(),
                "tree".into(),
            ]
        } else {
            vec![
                "filter".into(),
                "refilter".into(),
                "assemble".into(),
                "trim".into(),
                "combine".into(),
                "tree".into(),
            ]
        };
    }
    if commands.iter().any(|command| command == "gene") {
        let missing = ["filter", "refilter", "assemble"]
            .into_iter()
            .filter(|stage| !commands.iter().any(|command| command == stage))
            .collect::<Vec<_>>();
        if !missing.is_empty() {
            return Err(format!(
                "gene requires filter, refilter, and assemble; missing {}",
                missing.join(", ")
            ));
        }
    }
    let default_kf = if assembly_mode == "uce" {
        DEFAULT_UCE_FILTER_KMER
    } else {
        "31"
    };
    Ok(Options {
        raw: args.to_vec(),
        commands,
        reference: value(args, &["-r"], "")?,
        samples: value(args, &["-f"], "")?,
        output: value(args, &["-o"], "")?,
        assembly_mode,
        workers,
        worker_source,
        kf: value(args, &["-kf"], default_kf)?,
        step: value(args, &["-s", "--step-size"], "4")?,
        ka: value(args, &["-ka"], "0")?,
        min_ka: value(args, &["--min-ka"], "21")?,
        max_ka: value(args, &["--max-ka"], "51")?,
        error_threshold: value(args, &["-e", "--error-threshold"], "2")?,
        search_depth: value(args, &["-i", "--search-depth"], "4096")?,
        soft_boundary: value(args, &["-sb", "--soft-boundary"], "auto")?,
        min_coverage: value(args, &["--min-coverage"], "0")?,
        low_depth: value(args, &["--depth-low-water-mark"], "50")?,
        depth_limit: value(args, &["--depth-limit"], "768")?,
        size_limit: value(args, &["--file-size-limit"], "6")?,
        max_reads: value(args, &["--max-reads"], "0")?,
        graph_format: value(args, &["--assembler-graph-format"], "none")?,
        side_candidates: value(args, &["--uce-side-candidates"], "8")?,
        max_contig_length: value(args, &["--uce-max-contig-length"], "0")?,
        min_density: value(args, &["--uce-min-read-density"], "0.003")?,
        density_min_length: value(args, &["--uce-density-check-min-length"], "1000")?,
        max_depth_cv: value(args, &["--uce-max-depth-cv"], "0")?,
        max_depth_ratio: value(args, &["--uce-max-depth-ratio"], "0")?,
        alignment_shadow: flag(args, "--uce-alignment-shadow")?,
        shadow_per_locus: value(args, &["--uce-shadow-per-locus"], "64")?,
        shadow_band: value(args, &["--uce-shadow-band"], "32")?,
        shadow_terminal_window: value(args, &["--uce-shadow-terminal-window"], "150")?,
        uce_recruit_mode,
        uce_fallback_kmer_size,
        uce_fallback_step,
        uce_fallback_verify_kmer_size,
        uce_fallback_min_alignment_overlap,
        uce_fallback_min_alignment_identity,
        uce_memory_limit_mib: 0,
        rescue,
        stats_count_input_reads: flag(args, "--stats-count-input-reads")?,
        stats_no_heatmap: flag(args, "--stats-no-heatmap")?,
        cleanup_intermediates: flag(args, "--cleanup-intermediates")?,
        cleanup_dry_run: flag(args, "--cleanup-dry-run")?,
        reuse_reference_cache: flag(args, "--reuse-reference-cache")?,
        legacy_uce_filter,
        workflow_profile: flag(args, "--workflow-profile")?,
        resume: flag(args, "--resume")?,
        log_format,
    })
}

fn sample_name(raw: &str) -> String {
    let value: String = raw
        .trim()
        .chars()
        .filter_map(|c| match c {
            ' ' | '-' => Some('_'),
            c if c.is_alphanumeric() || c == '_' || c == '.' => Some(c),
            _ => None,
        })
        .collect();
    let mut chars = value.chars();
    let Some(first) = chars.next() else {
        return value;
    };
    let mut normalized = first.to_ascii_uppercase().to_string();
    normalized.extend(chars.map(|c| c.to_ascii_lowercase()));
    normalized
}

fn read_samples(path: &str, output: &Path) -> Result<Vec<Sample>, String> {
    read_samples_with_directory_creation(path, output, true)
}

fn read_samples_with_directory_creation(
    path: &str,
    output: &Path,
    create_directories: bool,
) -> Result<Vec<Sample>, String> {
    let file =
        fs::File::open(path).map_err(|e| format!("Unable to read sample list '{path}': {e}"))?;
    let mut rows = Vec::new();
    let mut names = std::collections::BTreeSet::new();
    for (index, line) in io::BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if !matches!(fields.len(), 2 | 3) {
            return Err(format!(
                "Sample row {} must be sample<TAB>R1[<TAB>R2]",
                index + 1
            ));
        }
        let name = sample_name(fields[0]);
        if name.is_empty() {
            return Err(format!(
                "Invalid sample name on row {}: '{}'",
                index + 1,
                fields[0]
            ));
        }
        if !names.insert(name.clone()) {
            return Err(format!("Duplicate sample name after normalization: {name}"));
        }
        if fields[1].is_empty() {
            return Err(format!("Sample row {} has an empty R1 path", index + 1));
        }
        if fields.len() == 3 && fields[2].is_empty() {
            return Err(format!("Sample row {} has an empty R2 path", index + 1));
        }
        for read in fields.iter().skip(1) {
            if !Path::new(read).is_file() {
                return Err(format!(
                    "Sample row {} read file does not exist: {read}",
                    index + 1
                ));
            }
        }
        rows.push((
            index,
            name,
            fields[1].to_owned(),
            fields.get(2).map(|read| (*read).to_owned()),
        ));
    }
    if rows.is_empty() {
        return Err("Sample list is empty or invalid".into());
    }

    let mut samples = Vec::with_capacity(rows.len());
    for (index, name, read1, read2) in rows {
        let numbered = format!("{}_{}", index + 1, name);
        if create_directories {
            fs::create_dir_all(output.join(&numbered)).map_err(|e| e.to_string())?;
        }
        // Keep the legacy two-column convention: a single supplied FASTX path is
        // deliberately used for both mates, rather than silently changing mode.
        samples.push(Sample {
            name: numbered,
            read1: read1.clone(),
            read2: read2.or(Some(read1)),
        });
    }
    Ok(samples)
}

fn read_rad_samples(path: &str) -> Result<Vec<Sample>, String> {
    let file = fs::File::open(path)
        .map_err(|e| format!("Unable to read RAD sample list '{path}': {e}"))?;
    let mut samples = Vec::new();
    let mut names = std::collections::BTreeSet::new();
    for (index, line) in io::BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 3 || fields[1].is_empty() || fields[2].is_empty() {
            return Err(format!(
                "RAD sample row {} must be: sample<TAB>R1.fastq<TAB>R2.fastq",
                index + 1
            ));
        }
        let name = sample_name(fields[0]);
        if name.is_empty() {
            return Err(format!("Invalid RAD sample name '{}'", fields[0]));
        }
        if !names.insert(name.clone()) {
            return Err(format!(
                "Duplicate RAD sample name after normalization: {name}"
            ));
        }
        for read in [fields[1], fields[2]] {
            if !Path::new(read).is_file() {
                return Err(format!("RAD read file does not exist: {read}"));
            }
        }
        samples.push(Sample {
            name,
            read1: fields[1].to_owned(),
            read2: Some(fields[2].to_owned()),
        });
    }
    if samples.is_empty() {
        return Err("RAD sample list has no paired reads".into());
    }
    Ok(samples)
}

fn components() -> Result<PathBuf, String> {
    if let Some(path) = env::var_os("GM2_COMPONENT_DIR") {
        return Ok(PathBuf::from(path));
    }
    let exe = env::current_exe().map_err(|e| e.to_string())?;
    Ok(exe
        .parent()
        .ok_or("Cannot locate GeneMiner components")?
        .to_path_buf())
}

fn run(binary_dir: &Path, name: &str, args: &[String]) -> Result<(), String> {
    let program = binary_dir.join(name);
    let status = Command::new(&program)
        .args(args)
        .status()
        .map_err(|e| format!("Unable to run {}: {e}", program.display()))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{} exited with {status}", program.display()))
    }
}

fn soft_boundary(value: &str) -> Result<String, String> {
    match value {
        "auto" => Ok("-1".into()),
        "unlimited" => Ok("10000".into()),
        other => other
            .parse::<i32>()
            .map_err(|_| "invalid --soft-boundary".into())
            .map(|v| v.to_string()),
    }
}

struct UceRecruitInvocation<'a> {
    sample_dir: &'a Path,
    verify_reference: &'a Path,
    recruit_reference: &'a Path,
    role: &'a str,
    retain_loci_file: Option<&'a Path>,
}

fn uce_filter_args_for_pass(
    opt: &Options,
    sample: &Sample,
    pass: &RecruitPass,
    invocation: &UceRecruitInvocation<'_>,
) -> Vec<String> {
    let mut args = vec![
        "-r".into(),
        invocation.verify_reference.display().to_string(),
        "--recruit-references".into(),
        invocation.recruit_reference.display().to_string(),
        "-q1".into(),
        sample.read1.clone(),
    ];
    if let Some(read2) = &sample.read2 {
        args.extend(["-q2".into(), read2.clone()]);
    }
    args.extend([
        "-o".into(),
        invocation.sample_dir.display().to_string(),
        "-kf".into(),
        pass.kmer_size.clone(),
        "-s".into(),
        pass.step.clone(),
        "--selection".into(),
        "auto".into(),
        "--reference-role".into(),
        invocation.role.into(),
        "--threads".into(),
        "1".into(),
        "--memory-limit-mib".into(),
        opt.uce_memory_limit_mib.to_string(),
        "--min-depth".into(),
        opt.low_depth.clone(),
        "--max-depth".into(),
        opt.depth_limit.clone(),
        "--max-size".into(),
        opt.size_limit.clone(),
    ]);
    if let Some(verify_kmer_size) = &pass.verify_kmer_size {
        args.extend(["--verification-kmer-size".into(), verify_kmer_size.clone()]);
    }
    if let Some(minimum_overlap) = &pass.minimum_alignment_overlap {
        args.extend([
            "--minimum-alignment-overlap".into(),
            minimum_overlap.clone(),
        ]);
    }
    if let Some(minimum_identity) = &pass.minimum_alignment_identity {
        args.extend([
            "--minimum-alignment-identity".into(),
            minimum_identity.clone(),
        ]);
    }
    if pass.max_locus_count > 0 {
        args.extend(["--max-locus-count".into(), pass.max_locus_count.to_string()]);
    }
    if let Some(path) = invocation.retain_loci_file {
        args.extend(["--retain-loci-file".into(), path.display().to_string()]);
    }
    if opt.max_reads != "0" {
        args.extend(["--max-fragments".into(), opt.max_reads.clone()]);
    }
    if opt.alignment_shadow {
        args.extend([
            "--alignment-shadow".into(),
            "--shadow-per-locus".into(),
            opt.shadow_per_locus.clone(),
            "--shadow-band".into(),
            opt.shadow_band.clone(),
            "--terminal-window".into(),
            opt.shadow_terminal_window.clone(),
        ]);
    }
    args
}

fn uce_filter_args_for_recruit(
    opt: &Options,
    sample: &Sample,
    sample_dir: &Path,
    verify_reference: &Path,
    recruit_reference: &Path,
    role: &str,
) -> Vec<String> {
    uce_filter_args_for_pass(
        opt,
        sample,
        &RecruitPass::fast(&opt.kf, &opt.step),
        &UceRecruitInvocation {
            sample_dir,
            verify_reference,
            recruit_reference,
            role,
            retain_loci_file: None,
        },
    )
}

fn uce_filter_args_for(
    opt: &Options,
    sample: &Sample,
    sample_dir: &Path,
    reference: &Path,
    role: &str,
) -> Vec<String> {
    uce_filter_args_for_recruit(opt, sample, sample_dir, reference, reference, role)
}

fn uce_filter_args(
    opt: &Options,
    sample: &Sample,
    sample_dir: &Path,
    threads: usize,
) -> Vec<String> {
    let mut args = uce_filter_args_for(opt, sample, sample_dir, Path::new(&opt.reference), "bait");
    let thread_index = args
        .iter()
        .position(|argument| argument == "--threads")
        .expect("UCEFilter arguments include --threads");
    args[thread_index + 1] = threads.to_string();
    args
}

fn uce_assembler_args(
    opt: &Options,
    sample_dir: &Path,
    threads: usize,
) -> Result<Vec<String>, String> {
    Ok(vec![
        "-r".into(),
        opt.reference.clone(),
        "-o".into(),
        sample_dir.display().to_string(),
        "-ka".into(),
        opt.ka.clone(),
        "-k_min".into(),
        opt.min_ka.clone(),
        "-k_max".into(),
        opt.max_ka.clone(),
        "-limit_count".into(),
        opt.error_threshold.clone(),
        "-iteration".into(),
        opt.search_depth.clone(),
        "-sb".into(),
        soft_boundary(&opt.soft_boundary)?,
        "-cov_min".into(),
        opt.min_coverage.clone(),
        "-p".into(),
        threads.to_string(),
        "--assembly-mode".into(),
        "uce".into(),
        "--uce-side-candidates".into(),
        opt.side_candidates.clone(),
        "--uce-max-contig-length".into(),
        opt.max_contig_length.clone(),
        "--uce-min-read-density".into(),
        opt.min_density.clone(),
        "--uce-density-check-min-length".into(),
        opt.density_min_length.clone(),
        "--uce-max-depth-cv".into(),
        opt.max_depth_cv.clone(),
        "--uce-max-depth-ratio".into(),
        opt.max_depth_ratio.clone(),
        "--uce-path-strategy".into(),
        value(&opt.raw, &["--uce-path-strategy"], "backbone")?,
        "--uce-backbone-lookahead".into(),
        value(&opt.raw, &["--uce-backbone-lookahead"], "24")?,
        "--assembler-read-chunk-size".into(),
        value(&opt.raw, &["--assembler-read-chunk-size"], "8192")?,
        "--assembler-kmer-count-threads".into(),
        // Zero delegates to the assembler's per-locus calculation.  Passing
        // the sample budget here would create nested `threads × threads`
        // k-mer workers on high-core hosts.
        "0".into(),
        "--assembler-graph-format".into(),
        opt.graph_format.clone(),
    ])
}

fn uce_rescue_assembler_args(
    opt: &Options,
    sample_dir: &Path,
    reference: &Path,
    threads: usize,
) -> Result<Vec<String>, String> {
    let scale = value(
        &opt.raw,
        &["--uce-rescue-reverse-reuse-reference-scale"],
        "1.0",
    )?;
    let parsed = scale
        .parse::<f64>()
        .map_err(|_| "--uce-rescue-reverse-reuse-reference-scale must be a number".to_string())?;
    if !(0.0..=1.0).contains(&parsed) {
        return Err("--uce-rescue-reverse-reuse-reference-scale must be in [0, 1]".into());
    }
    let mut rescue_opt = opt.clone();
    rescue_opt.reference = reference.display().to_string();
    rescue_opt.ka = UCE_RESCUE_ASSEMBLY_KMER.into();
    let mut args = uce_assembler_args(&rescue_opt, sample_dir, threads)?;
    args.extend(["--uce-reverse-reuse-reference-scale".into(), scale]);
    Ok(args)
}

fn build_uce_rescue_reference(
    reference: &Path,
    sample: &Path,
    rescue: &Path,
    minimum: usize,
    active: Option<&std::collections::BTreeSet<String>>,
) -> Result<usize, String> {
    if rescue.exists() {
        fs::remove_dir_all(rescue).map_err(|e| e.to_string())?;
    }
    if active.is_some() {
        fs::create_dir_all(rescue).map_err(|e| e.to_string())?;
    } else {
        copy_tree(reference, rescue)?;
    }
    let summary = read_uce_summary(&sample.join("uce_assembly_summary.csv"))?;
    let mut added = 0;
    for (locus, source) in reference_loci(reference)? {
        if active.is_some_and(|loci| !loci.contains(&locus)) {
            continue;
        }
        let file_name = source.file_name().ok_or("invalid UCE reference filename")?;
        let rescue_path = rescue.join(file_name);
        if active.is_some() {
            fs::copy(&source, &rescue_path).map_err(|e| e.to_string())?;
        }
        if !uce_row_accepted(summary.rows.get(&locus)) {
            continue;
        }
        let contig = sample.join("results").join(format!("{locus}.fasta"));
        if !contig.is_file() {
            continue;
        }
        let mut target = fs::OpenOptions::new()
            .append(true)
            .open(rescue_path)
            .map_err(|e| e.to_string())?;
        use std::io::Write;
        for (index, (_, sequence)) in fasta_records(&contig)?.into_iter().enumerate() {
            if sequence.len() >= minimum {
                writeln!(
                    target,
                    ">{locus}_gm2_rescue_contig_{}\n{sequence}",
                    index + 1
                )
                .map_err(|e| e.to_string())?;
                added += 1;
            }
        }
    }
    Ok(added)
}

fn review_only_provisional_cores(summary: &UceSummary) -> BTreeSet<String> {
    summary
        .rows
        .iter()
        .filter(|(_, row)| {
            row.get("auto_recruit_core_anchor_status")
                .is_some_and(|status| status == "anchored_with_review")
        })
        .map(|(locus, _)| locus.clone())
        .collect()
}

fn build_uce_terminal_baits(
    sample: &Path,
    baits: &Path,
    active: &std::collections::BTreeSet<String>,
    window: usize,
    minimum: usize,
) -> Result<std::collections::BTreeSet<String>, String> {
    if baits.exists() {
        fs::remove_dir_all(baits).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(baits).map_err(|e| e.to_string())?;
    let summary = read_uce_summary(&sample.join("uce_assembly_summary.csv"))?;
    let mut written = std::collections::BTreeSet::new();
    for locus in active {
        if !uce_row_accepted(summary.rows.get(locus)) {
            continue;
        }
        let Some(sequence) =
            first_fasta_sequence(&sample.join("results").join(format!("{locus}.fasta")))?
        else {
            continue;
        };
        if sequence.len() < minimum {
            continue;
        }
        let flank = window.max(minimum).min(sequence.len());
        let left = &sequence[..flank];
        let right = &sequence[sequence.len() - flank..];
        let mut text = format!(">{locus}_gm2_left_terminal\n{left}\n");
        if left != right {
            text.push_str(&format!(">{locus}_gm2_right_terminal\n{right}\n"));
        }
        fs::write(baits.join(format!("{locus}.fasta")), text).map_err(|e| e.to_string())?;
        written.insert(locus.clone());
    }
    Ok(written)
}

fn restore_locus_file(
    sample: &Path,
    backup: &Path,
    directory: &str,
    locus: &str,
) -> Result<(), String> {
    let original = backup.join(directory).join(format!("{locus}.fasta"));
    let current = sample.join(directory).join(format!("{locus}.fasta"));
    if current.exists() {
        fs::remove_file(&current).map_err(|e| e.to_string())?;
    }
    if original.is_file() {
        if let Some(parent) = current.parent() {
            fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        }
        fs::copy(original, current).map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn locus_file_name_matches(name: &str, locus: &str, paired: bool) -> bool {
    let stem = Path::new(name)
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or_default();
    if paired {
        stem == format!("{locus}_1") || stem == format!("{locus}_2")
    } else {
        stem == locus
    }
}

fn restore_locus_directory_files(
    sample: &Path,
    backup: &Path,
    directory: &str,
    locus: &str,
) -> Result<(), String> {
    let source_dir = backup.join(directory);
    let destination_dir = sample.join(directory);
    let mut names = std::collections::BTreeSet::new();
    for root in [&source_dir, &destination_dir] {
        if !root.is_dir() {
            continue;
        }
        for entry in fs::read_dir(root).map_err(|e| e.to_string())? {
            let entry = entry.map_err(|e| e.to_string())?;
            if entry.path().is_file() {
                names.insert(entry.file_name());
            }
        }
    }
    for name in names {
        let text = name.to_string_lossy();
        if !locus_file_name_matches(&text, locus, directory == "filtered_pe") {
            continue;
        }
        let source = source_dir.join(&name);
        let destination = destination_dir.join(&name);
        if destination.exists() {
            fs::remove_file(&destination).map_err(|e| e.to_string())?;
        }
        if source.is_file() {
            fs::create_dir_all(&destination_dir).map_err(|e| e.to_string())?;
            fs::copy(source, destination).map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

fn count_rows(path: &Path) -> Result<Vec<String>, String> {
    if !path.is_file() {
        return Ok(Vec::new());
    }
    Ok(fs::read_to_string(path)
        .map_err(|e| e.to_string())?
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(str::to_owned)
        .collect())
}

fn restore_locus_read_count(sample: &Path, backup: &Path, locus: &str) -> Result<(), String> {
    let filename = "ref_reads_count_dict.txt";
    let source = backup.join(filename);
    let destination = sample.join(filename);
    let backup_rows = count_rows(&source)?;
    let current_rows = count_rows(&destination)?;
    let mut merged = current_rows
        .into_iter()
        .filter(|line| line.split(',').next() != Some(locus))
        .collect::<Vec<_>>();
    merged.extend(
        backup_rows
            .into_iter()
            .filter(|line| line.split(',').next() == Some(locus)),
    );
    if merged.is_empty() {
        if destination.exists() {
            fs::remove_file(destination).map_err(|e| e.to_string())?;
        }
    } else {
        fs::write(destination, format!("{}\n", merged.join("\n"))).map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn write_result_dict_from_uce_summary(sample: &Path, summary: &UceSummary) -> Result<(), String> {
    let mut text = String::new();
    for (locus, row) in &summary.rows {
        if row.get("status").is_some_and(|status| status == "skipped") {
            continue;
        }
        text.push_str(&format!(
            "{},{},{},\n",
            locus,
            row.get("status").map(String::as_str).unwrap_or_default(),
            row.get("read_count")
                .map(String::as_str)
                .unwrap_or_default()
        ));
    }
    fs::write(sample.join("result_dict.txt"), text).map_err(|e| e.to_string())
}

fn restore_rescue_locus(sample: &Path, backup: &Path, locus: &str) -> Result<(), String> {
    for directory in ["results", "contigs_all", "contigs_all_low"] {
        restore_locus_file(sample, backup, directory, locus)?;
    }
    for directory in ["filtered", "filtered_pe"] {
        restore_locus_directory_files(sample, backup, directory, locus)?;
    }
    restore_locus_read_count(sample, backup, locus)
}

#[derive(Clone, Debug, Default)]
struct TerminalSideEvidence {
    length: usize,
    breadth: f64,
    max_gap: usize,
    fragments: usize,
    bridges: usize,
    accepted: bool,
}

#[derive(Clone, Debug, Default)]
struct TerminalEvidence {
    left: TerminalSideEvidence,
    right: TerminalSideEvidence,
}

#[derive(Default)]
struct RescueReportContext {
    round_statuses: std::collections::BTreeMap<(usize, String), String>,
    terminal_audits: std::collections::BTreeMap<(usize, String), TerminalEvidence>,
    status_by_locus: std::collections::BTreeMap<String, String>,
    overall_status: String,
}

struct RescueRoundOutcome {
    after: UceSummary,
    statuses: std::collections::BTreeMap<String, String>,
    evidence_by_locus: std::collections::BTreeMap<String, TerminalEvidence>,
}

fn reverse_complement_text(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

fn locus_result_sequence(root: &Path, locus: &str) -> Result<Option<String>, String> {
    let path = root.join("results").join(format!("{locus}.fasta"));
    if !path.is_file() {
        return Ok(None);
    }
    first_fasta_sequence(&path)
}

fn rescue_introduces_long_inverted_repeat(
    sample: &Path,
    backup: &Path,
    locus: &str,
    minimum_span: usize,
) -> Result<bool, String> {
    if minimum_span == 0 {
        return Ok(false);
    }
    let before = locus_result_sequence(backup, locus)?;
    let after = locus_result_sequence(sample, locus)?;
    let before_has_repeat = before
        .as_deref()
        .is_some_and(|sequence| rescue_qc::has_long_inverted_repeat(sequence, minimum_span));
    let after_has_repeat = after
        .as_deref()
        .is_some_and(|sequence| rescue_qc::has_long_inverted_repeat(sequence, minimum_span));
    Ok(!before_has_repeat && after_has_repeat)
}

fn rescue_introduces_unsupported_internal_gap(
    sample: &Path,
    backup: &Path,
    locus: &str,
) -> Result<bool, String> {
    let before = locus_result_sequence(backup, locus)?;
    let after = locus_result_sequence(sample, locus)?;
    let (Some(before), Some(after)) = (before, after) else {
        return Ok(false);
    };
    let reads = read_locus_fastq(sample, locus)?;
    Ok(rescue_qc::introduces_unsupported_internal_gap(
        &before, &after, &reads,
    ))
}

fn read_locus_fastq(sample: &Path, locus: &str) -> Result<Vec<(String, String)>, String> {
    let path = sample.join("filtered").join(format!("{locus}.fq"));
    if !path.is_file() {
        return Ok(Vec::new());
    }
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let lines = text.lines().collect::<Vec<_>>();
    let mut reads = Vec::new();
    for record in lines.chunks_exact(4) {
        let title = record[0]
            .trim_start_matches('@')
            .split_whitespace()
            .next()
            .unwrap_or_default();
        let fragment = title
            .rsplit_once('/')
            .map_or(title, |(prefix, _)| prefix)
            .to_owned();
        reads.push((fragment, record[1].trim().to_ascii_uppercase()));
    }
    Ok(reads)
}

fn maximum_false_run(values: &[bool]) -> usize {
    let mut longest = 0;
    let mut current = 0;
    for value in values {
        if *value {
            current = 0;
        } else {
            current += 1;
            longest = longest.max(current);
        }
    }
    longest
}

fn terminal_support_metrics(
    sequence: &str,
    old_start: usize,
    old_end: usize,
    reads: &[(String, String)],
) -> (
    TerminalSideEvidence,
    TerminalSideEvidence,
    Vec<bool>,
    std::collections::HashMap<String, u8>,
) {
    const LEFT_EXTENSION: u8 = 1;
    const RIGHT_EXTENSION: u8 = 2;
    const LEFT_CORE: u8 = 4;
    const RIGHT_CORE: u8 = 8;
    const KMER_SIZE: usize = 21;

    let bytes = sequence.as_bytes();
    let mut covered = vec![false; bytes.len()];
    let left_core_end = old_end.min(old_start.saturating_add(150));
    let right_core_start = old_start.max(old_end.saturating_sub(150));
    let mut positions: std::collections::HashMap<Vec<u8>, Vec<usize>> =
        std::collections::HashMap::new();
    if bytes.len() >= KMER_SIZE {
        for start in 0..=bytes.len() - KMER_SIZE {
            let kmer = &bytes[start..start + KMER_SIZE];
            if kmer
                .iter()
                .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            {
                positions
                    .entry(kmer.iter().map(u8::to_ascii_uppercase).collect())
                    .or_default()
                    .push(start);
            }
        }
    }

    let mut fragment_regions = std::collections::HashMap::<String, u8>::new();
    for (fragment, read) in reads {
        fragment_regions.entry(fragment.clone()).or_default();
        let mut observed = std::collections::HashSet::<Vec<u8>>::new();
        for oriented in [read.clone(), reverse_complement_text(read)] {
            let oriented = oriented.as_bytes();
            if oriented.len() < KMER_SIZE {
                continue;
            }
            for offset in 0..=oriented.len() - KMER_SIZE {
                observed.insert(oriented[offset..offset + KMER_SIZE].to_vec());
            }
        }
        for kmer in observed {
            for start in positions.get(&kmer).into_iter().flatten().copied() {
                let end = start + KMER_SIZE;
                covered[start..end].fill(true);
                let regions = fragment_regions.entry(fragment.clone()).or_default();
                if start < old_start {
                    *regions |= LEFT_EXTENSION;
                }
                if end > old_end {
                    *regions |= RIGHT_EXTENSION;
                }
                if end > old_start && start < left_core_end {
                    *regions |= LEFT_CORE;
                }
                if end > right_core_start && start < old_end {
                    *regions |= RIGHT_CORE;
                }
            }
        }
    }

    let side = |start: usize, end: usize, extension: u8, core: u8| {
        let length = end.saturating_sub(start);
        if length == 0 {
            return TerminalSideEvidence {
                breadth: 1.0,
                ..TerminalSideEvidence::default()
            };
        }
        let side_coverage = &covered[start..end];
        let fragments = fragment_regions
            .values()
            .filter(|regions| **regions & extension != 0)
            .count();
        let bridges = fragment_regions
            .values()
            .filter(|regions| **regions & extension != 0 && **regions & core != 0)
            .count();
        let breadth = side_coverage.iter().filter(|value| **value).count() as f64 / length as f64;
        let max_gap = maximum_false_run(side_coverage);
        TerminalSideEvidence {
            length,
            breadth,
            max_gap,
            fragments,
            bridges,
            accepted: length >= UCE_TERMINAL_MIN_EXTENSION
                && breadth >= UCE_TERMINAL_MIN_BREADTH
                && max_gap <= UCE_TERMINAL_MAX_GAP
                && fragments >= UCE_TERMINAL_MIN_FRAGMENTS
                && bridges >= UCE_TERMINAL_MIN_BRIDGES,
        }
    };
    (
        side(0, old_start, LEFT_EXTENSION, LEFT_CORE),
        side(old_end, bytes.len(), RIGHT_EXTENSION, RIGHT_CORE),
        covered,
        fragment_regions,
    )
}

fn write_trimmed_locus_sequence(sample: &Path, locus: &str, sequence: &str) -> Result<(), String> {
    let path = sample.join("results").join(format!("{locus}.fasta"));
    if !path.is_file() {
        return Ok(());
    }
    let title = fasta_records(&path)?
        .into_iter()
        .next()
        .map(|record| record.0)
        .unwrap_or_else(|| locus.to_owned());
    fs::write(path, format!(">{title}\n{sequence}\n")).map_err(|e| e.to_string())
}

fn terminal_reconcile_locus(
    sample: &Path,
    backup: &Path,
    locus: &str,
    after_row: &mut std::collections::BTreeMap<String, String>,
) -> Result<(Option<TerminalEvidence>, String), String> {
    let Some(old_sequence) =
        first_fasta_sequence(&backup.join("results").join(format!("{locus}.fasta")))?
    else {
        return Ok((None, "missing_contig".into()));
    };
    let Some(mut new_sequence) =
        first_fasta_sequence(&sample.join("results").join(format!("{locus}.fasta")))?
    else {
        return Ok((None, "missing_contig".into()));
    };
    let old_start = if let Some(position) = new_sequence.find(&old_sequence) {
        position
    } else {
        let reverse = reverse_complement_text(&new_sequence);
        let Some(position) = reverse.find(&old_sequence) else {
            return Ok((None, "core_changed".into()));
        };
        new_sequence = reverse;
        position
    };
    let old_end = old_start + old_sequence.len();
    let reads = read_locus_fastq(sample, locus)?;
    let (left, right, covered, fragment_regions) =
        terminal_support_metrics(&new_sequence, old_start, old_end, &reads);
    let kept_left = if left.accepted {
        &new_sequence[..old_start]
    } else {
        ""
    };
    let kept_right = if right.accepted {
        &new_sequence[old_end..]
    } else {
        ""
    };
    let accepted_sequence = format!("{kept_left}{old_sequence}{kept_right}");
    let evidence = TerminalEvidence { left, right };
    if accepted_sequence == old_sequence {
        return Ok((Some(evidence), "no_supported_extension".into()));
    }

    write_trimmed_locus_sequence(sample, locus, &accepted_sequence)?;
    let accepted_start = if evidence.left.accepted { 0 } else { old_start };
    let accepted_end = if evidence.right.accepted {
        new_sequence.len()
    } else {
        old_end
    };
    let accepted_coverage = &covered[accepted_start..accepted_end];
    let supported = accepted_coverage.iter().filter(|value| **value).count();
    let supported_positions = accepted_coverage
        .iter()
        .enumerate()
        .filter_map(|(index, value)| value.then_some(index))
        .collect::<Vec<_>>();
    let supported_span = supported_positions
        .first()
        .zip(supported_positions.last())
        .map_or(0, |(first, last)| last - first + 1);
    let fragments = fragment_regions
        .values()
        .filter(|regions| **regions != 0)
        .count();
    let length = accepted_sequence.len();
    after_row.insert("selected_contig_length".into(), length.to_string());
    after_row.insert("read_supported_span".into(), supported_span.to_string());
    after_row.insert("slice_supported_bases".into(), supported.to_string());
    after_row.insert(
        "slice_support_breadth".into(),
        format!("{:.6}", supported as f64 / length as f64),
    );
    after_row.insert(
        "max_slice_support_gap".into(),
        maximum_false_run(accepted_coverage).to_string(),
    );
    after_row.insert("read_count".into(), fragments.to_string());
    if let Some(unique) = uce_number(Some(after_row), "unique_read_count") {
        after_row.insert(
            "read_density".into(),
            format!("{:.6}", fragments as f64 / length as f64),
        );
        after_row.insert(
            "unique_read_density".into(),
            format!("{:.6}", unique as f64 / length as f64),
        );
    }
    Ok((Some(evidence), "accepted".into()))
}

fn restore_prior_rescue_rounds(
    sample: &Path,
    backup: &Path,
    current_round: usize,
) -> Result<(), String> {
    for round in 1..current_round {
        let name = format!("uce_rescue_round_{round}");
        let source = backup.join(&name);
        let destination = sample.join(&name);
        if source.is_dir() && !destination.exists() {
            copy_tree(&source, &destination)?;
        }
    }
    Ok(())
}

fn restore_initial_uce_recruit_audits(sample: &Path, backup: &Path) -> Result<(), String> {
    for name in [
        "uce_filter_summary.fast.tsv",
        "uce_filter_summary.fallback.tsv",
        "uce_recruit_passes.tsv",
        "uce_recruit_contig_probe_gate.tsv",
    ] {
        let source = backup.join(name);
        let destination = sample.join(name);
        if source.is_file() && !destination.exists() {
            fs::copy(&source, &destination).map_err(|error| {
                format!(
                    "Unable to preserve UCE recruitment audit '{}' as '{}': {error}",
                    source.display(),
                    destination.display()
                )
            })?;
        }
    }

    let source = backup.join("fallback_probe_rejected");
    let destination = sample.join("fallback_probe_rejected");
    if source.is_dir() && !destination.exists() {
        copy_tree(&source, &destination)?;
    }
    Ok(())
}

fn row_density(row: Option<&std::collections::BTreeMap<String, String>>) -> Option<f64> {
    let length = uce_number(row, "selected_contig_length")? as f64;
    let reads =
        uce_number(row, "unique_read_count").or_else(|| uce_number(row, "read_count"))? as f64;
    (length > 0.0 && reads >= 0.0).then_some(reads / length)
}

fn write_rescue_reports(
    sample: &Sample,
    directory: &Path,
    initial: &UceSummary,
    final_rows: &UceSummary,
    rounds: &[(usize, String, UceSummary, UceSummary)],
    report: &RescueReportContext,
) -> Result<(), String> {
    let mut round_csv = String::from("sample,round,locus,round_status,before_status,after_status,before_length,after_length,length_delta,before_unique_reads,after_unique_reads,unique_read_delta,left_extension_length,left_breadth,left_max_gap,left_fragments,left_bridges,left_accepted,right_extension_length,right_breadth,right_max_gap,right_fragments,right_bridges,right_accepted\n");
    for (round, status, before, after) in rounds {
        for locus in before
            .rows
            .keys()
            .chain(after.rows.keys())
            .collect::<std::collections::BTreeSet<_>>()
        {
            let left = before.rows.get(locus);
            let right = after.rows.get(locus);
            let length_delta = uce_number(right, "selected_contig_length")
                .zip(uce_number(left, "selected_contig_length"))
                .map(|(a, b)| (a - b).to_string())
                .unwrap_or_default();
            let key = (*round, locus.to_string());
            let decision = report
                .round_statuses
                .get(&key)
                .map(String::as_str)
                .unwrap_or(status);
            if *round > 1 && decision == "stable_not_recruited" {
                continue;
            }
            let read_delta = uce_number(right, "unique_read_count")
                .zip(uce_number(left, "unique_read_count"))
                .map(|(a, b)| (a - b).to_string())
                .unwrap_or_default();
            let mut fields = vec![
                sample.name.clone(),
                round.to_string(),
                locus.to_string(),
                decision.to_owned(),
                left.and_then(|row| row.get("status"))
                    .cloned()
                    .unwrap_or_default(),
                right
                    .and_then(|row| row.get("status"))
                    .cloned()
                    .unwrap_or_default(),
                left.and_then(|row| row.get("selected_contig_length"))
                    .cloned()
                    .unwrap_or_default(),
                right
                    .and_then(|row| row.get("selected_contig_length"))
                    .cloned()
                    .unwrap_or_default(),
                length_delta,
                left.and_then(|row| row.get("unique_read_count"))
                    .cloned()
                    .unwrap_or_default(),
                right
                    .and_then(|row| row.get("unique_read_count"))
                    .cloned()
                    .unwrap_or_default(),
                read_delta,
            ];
            if let Some(evidence) = report.terminal_audits.get(&key) {
                for side in [&evidence.left, &evidence.right] {
                    fields.extend([
                        side.length.to_string(),
                        format!("{:.6}", side.breadth),
                        side.max_gap.to_string(),
                        side.fragments.to_string(),
                        side.bridges.to_string(),
                        u8::from(side.accepted).to_string(),
                    ]);
                }
            } else {
                fields.extend(std::iter::repeat_n(String::new(), 12));
            }
            round_csv.push_str(&fields.join(","));
            round_csv.push('\n');
        }
    }
    fs::write(directory.join("uce_rescue_rounds.csv"), round_csv).map_err(|e| e.to_string())?;
    let mut summary_csv = String::from("sample,locus,rescue_status,before_status,after_status,before_length,after_length,length_delta,before_read_density,after_read_density\n");
    for locus in initial
        .rows
        .keys()
        .chain(final_rows.rows.keys())
        .collect::<std::collections::BTreeSet<_>>()
    {
        let before = initial.rows.get(locus);
        let after = final_rows.rows.get(locus);
        let delta = uce_number(after, "selected_contig_length")
            .zip(uce_number(before, "selected_contig_length"))
            .map(|(a, b)| (a - b).to_string())
            .unwrap_or_default();
        let rescue_status = report
            .status_by_locus
            .get(locus)
            .map(String::as_str)
            .unwrap_or(&report.overall_status);
        summary_csv.push_str(&format!(
            "{},{},{},{},{},{},{},{},{:.6},{:.6}\n",
            sample.name,
            locus,
            rescue_status,
            before
                .and_then(|r| r.get("status"))
                .cloned()
                .unwrap_or_default(),
            after
                .and_then(|r| r.get("status"))
                .cloned()
                .unwrap_or_default(),
            before
                .and_then(|r| r.get("selected_contig_length"))
                .cloned()
                .unwrap_or_default(),
            after
                .and_then(|r| r.get("selected_contig_length"))
                .cloned()
                .unwrap_or_default(),
            delta,
            row_density(before).unwrap_or(0.0),
            row_density(after).unwrap_or(0.0)
        ));
    }
    fs::write(directory.join("uce_rescue_summary.csv"), summary_csv).map_err(|e| e.to_string())
}

fn execute_uce_rescue(
    opt: &Options,
    bins: &Path,
    sample: &Sample,
    sample_dir: &Path,
) -> Result<(), String> {
    let minimum = raw_number::<usize>(
        &opt.raw,
        &["--uce-rescue-min-contig-length"],
        "60",
        "--uce-rescue-min-contig-length",
    )?
    .max(opt.kf.parse().unwrap_or(1));
    let maximum_rounds = raw_number::<usize>(
        &opt.raw,
        &["--uce-rescue-rounds"],
        DEFAULT_UCE_RESCUE_ROUNDS,
        "--uce-rescue-rounds",
    )?
    .clamp(1, 2);
    let terminal_window = raw_number::<usize>(
        &opt.raw,
        &["--uce-rescue-terminal-window"],
        "350",
        "--uce-rescue-terminal-window",
    )?
    .max(minimum);
    let density_ratio = value(&opt.raw, &["--uce-rescue-min-density-ratio"], "0.5")?
        .parse::<f64>()
        .map_err(|_| "--uce-rescue-min-density-ratio must be a number")?;
    if !(0.0..=1.0).contains(&density_ratio) {
        return Err("--uce-rescue-min-density-ratio must be in [0, 1]".into());
    }
    let inverted_repeat_minimum = raw_number::<usize>(
        &opt.raw,
        &["--uce-rescue-inverted-repeat-min-bp"],
        "150",
        "--uce-rescue-inverted-repeat-min-bp",
    )?;
    let initial = read_uce_summary(&sample_dir.join("uce_assembly_summary.csv"))?;
    let review_only_cores = review_only_provisional_cores(&initial);
    let mut current = initial.clone();
    let mut previous = initial.clone();
    let mut records: Vec<(usize, String, UceSummary)> = Vec::new();
    let mut report = RescueReportContext {
        overall_status: "success".into(),
        ..RescueReportContext::default()
    };
    for round in 1..=maximum_rounds {
        let candidate = if round == 1 {
            None
        } else {
            Some(terminal_rescue_loci(&previous, &current))
        };
        if candidate.as_ref().is_some_and(|loci| loci.is_empty()) {
            break;
        }
        let active = if round == 1 {
            if review_only_cores.is_empty() {
                None
            } else {
                Some(
                    reference_loci(Path::new(&opt.reference))?
                        .into_iter()
                        .map(|(locus, _)| locus)
                        .filter(|locus| !review_only_cores.contains(locus))
                        .collect::<BTreeSet<_>>(),
                )
            }
        } else {
            let mut active = candidate.clone().unwrap_or_default();
            active.retain(|locus| !review_only_cores.contains(locus));
            Some(active)
        };
        if active.as_ref().is_some_and(|loci| loci.is_empty()) {
            break;
        }
        // Keep rescue-only inputs outside the sample directory. The sample
        // directory can then be renamed to its rollback backup instead of
        // copied byte-for-byte before this round rebuilds filtered reads and
        // assembly graphs. The stage is moved back under the sample only
        // after the round has settled, preserving the historical layout.
        let root = Path::new(&opt.output)
            .join(".uce_rescue_stage")
            .join(&sample.name)
            .join(format!("round_{round}"));
        let reference = root.join("assembly_refs");
        let added = build_uce_rescue_reference(
            Path::new(&opt.reference),
            sample_dir,
            &reference,
            minimum,
            active.as_ref(),
        )?;
        if added == 0 {
            if round == 1 {
                report.overall_status = "skipped".into();
            }
            break;
        }
        let recruit = if round > 1 {
            let active = active
                .as_ref()
                .ok_or("terminal UCE rescue has no active locus set")?;
            let terminal = root.join("terminal_baits");
            let baits =
                build_uce_terminal_baits(sample_dir, &terminal, active, terminal_window, minimum)?;
            if baits.is_empty() {
                break;
            }
            terminal
        } else {
            reference.clone()
        };
        let backup = Path::new(&opt.output)
            .join(".uce_rescue_backups")
            .join(format!("{}_r{round}", sample.name));
        if backup.exists() {
            fs::remove_dir_all(&backup).map_err(|e| e.to_string())?;
        }
        move_tree(sample_dir, &backup)?;
        let before = current.clone();
        let result: Result<RescueRoundOutcome, String> = (|| {
            let filtered = sample_dir.join("filtered");
            if filtered.exists() {
                fs::remove_dir_all(&filtered).map_err(|e| e.to_string())?;
            }
            run(
                bins,
                "uce_filter",
                &uce_filter_args_for_recruit(
                    opt, sample, sample_dir, &reference, &recruit, "contig",
                ),
            )?;
            run(
                bins,
                "main_assembler-rust",
                &uce_rescue_assembler_args(opt, sample_dir, &reference, 1)?,
            )?;
            let mut after = read_uce_summary(&sample_dir.join("uce_assembly_summary.csv"))?;
            let mut statuses = std::collections::BTreeMap::new();
            let mut evidence_by_locus = std::collections::BTreeMap::new();
            for (locus, before_row) in &before.rows {
                let inactive = active
                    .as_ref()
                    .is_some_and(|active| !active.contains(locus));
                if inactive {
                    restore_rescue_locus(sample_dir, &backup, locus)?;
                    after.rows.insert(locus.clone(), before_row.clone());
                    statuses.insert(
                        locus.clone(),
                        if review_only_cores.contains(locus) {
                            "stable_review_only_core"
                        } else {
                            "stable_not_recruited"
                        }
                        .into(),
                    );
                    continue;
                }

                if !uce_row_accepted(after.rows.get(locus)) {
                    let status = if uce_row_accepted(Some(before_row)) {
                        "reverted_failed_rescue"
                    } else {
                        "not_recovered"
                    };
                    restore_rescue_locus(sample_dir, &backup, locus)?;
                    after.rows.insert(locus.clone(), before_row.clone());
                    statuses.insert(locus.clone(), status.into());
                    continue;
                }

                let density_drop = uce_row_accepted(Some(before_row))
                    && row_density(Some(before_row))
                        .zip(row_density(after.rows.get(locus)))
                        .is_some_and(|(old, new)| old > 0.0 && new / old < density_ratio);
                if density_drop {
                    restore_rescue_locus(sample_dir, &backup, locus)?;
                    after.rows.insert(locus.clone(), before_row.clone());
                    statuses.insert(locus.clone(), "reverted_density_drop".into());
                    continue;
                }

                let accepted_status = if round > 1 && uce_row_accepted(Some(before_row)) {
                    let after_row = after.rows.get_mut(locus).ok_or_else(|| {
                        format!("accepted rescue locus {locus} is missing from summary")
                    })?;
                    let (evidence, terminal_status) =
                        terminal_reconcile_locus(sample_dir, &backup, locus, after_row)?;
                    if let Some(evidence) = evidence {
                        evidence_by_locus.insert(locus.clone(), evidence);
                    }
                    if terminal_status == "missing_contig" || terminal_status == "core_changed" {
                        restore_rescue_locus(sample_dir, &backup, locus)?;
                        after.rows.insert(locus.clone(), before_row.clone());
                        statuses.insert(locus.clone(), format!("reverted_{terminal_status}"));
                        continue;
                    }
                    if terminal_status == "no_supported_extension" {
                        restore_rescue_locus(sample_dir, &backup, locus)?;
                        after.rows.insert(locus.clone(), before_row.clone());
                        statuses.insert(locus.clone(), "stable_no_supported_extension".into());
                        continue;
                    }
                    let evidence = evidence_by_locus
                        .get(locus)
                        .ok_or("terminal rescue accepted without evidence")?;
                    format!(
                        "terminal_left_{}_right_{}",
                        if evidence.left.accepted {
                            "kept"
                        } else {
                            "trimmed"
                        },
                        if evidence.right.accepted {
                            "kept"
                        } else {
                            "trimmed"
                        }
                    )
                } else {
                    "accepted".into()
                };

                if rescue_introduces_long_inverted_repeat(
                    sample_dir,
                    &backup,
                    locus,
                    inverted_repeat_minimum,
                )? {
                    restore_rescue_locus(sample_dir, &backup, locus)?;
                    after.rows.insert(locus.clone(), before_row.clone());
                    evidence_by_locus.remove(locus);
                    statuses.insert(locus.clone(), "reverted_inverted_repeat".into());
                    continue;
                }
                if uce_row_accepted(Some(before_row))
                    && rescue_introduces_unsupported_internal_gap(sample_dir, &backup, locus)?
                {
                    restore_rescue_locus(sample_dir, &backup, locus)?;
                    after.rows.insert(locus.clone(), before_row.clone());
                    evidence_by_locus.remove(locus);
                    statuses.insert(locus.clone(), "reverted_unsupported_internal_gap".into());
                    continue;
                }
                statuses.insert(locus.clone(), accepted_status);
            }
            write_uce_summary(&sample_dir.join("uce_assembly_summary.csv"), &after)?;
            write_result_dict_from_uce_summary(sample_dir, &after)?;
            restore_prior_rescue_rounds(sample_dir, &backup, round)?;
            restore_initial_uce_recruit_audits(sample_dir, &backup)?;
            Ok(RescueRoundOutcome {
                after,
                statuses,
                evidence_by_locus,
            })
        })();
        match result {
            Ok(outcome) => {
                let RescueRoundOutcome {
                    after,
                    statuses,
                    evidence_by_locus,
                } = outcome;
                for (locus, status) in statuses {
                    report
                        .round_statuses
                        .insert((round, locus.clone()), status.clone());
                    if status != "stable_not_recruited" {
                        report.status_by_locus.insert(locus, status);
                    }
                }
                for (locus, evidence) in evidence_by_locus {
                    report.terminal_audits.insert((round, locus), evidence);
                }
                records.push((
                    round,
                    if round == 1 {
                        "whole-contig"
                    } else {
                        "terminal-only"
                    }
                    .into(),
                    before,
                ));
                records.push((round, "accepted".into(), after.clone()));
                previous = current;
                current = after;
            }
            Err(error) => {
                if sample_dir.exists() {
                    fs::remove_dir_all(sample_dir).map_err(|e| e.to_string())?;
                }
                move_tree(&backup, sample_dir)?;
                report.overall_status = if records.is_empty() {
                    "failed_rolled_back".into()
                } else {
                    format!(
                        "success_round_{}_round_{round}_failed_rolled_back",
                        round - 1
                    )
                };
                eprintln!(
                    "Warning: UCE rescue round {round} rolled back for {}: {error}",
                    sample.name
                );
                break;
            }
        }
        if backup.exists() {
            fs::remove_dir_all(&backup).map_err(|e| e.to_string())?;
        }
        let settled_root = sample_dir.join(format!("uce_rescue_round_{round}"));
        if settled_root.exists() {
            fs::remove_dir_all(&settled_root).map_err(|e| e.to_string())?;
        }
        if root.exists() {
            move_tree(&root, &settled_root)?;
        }
    }
    // Pair each before/after state for compact reports.
    let pairs = records
        .chunks(2)
        .filter(|pair| pair.len() == 2)
        .map(|pair| {
            (
                pair[0].0,
                pair[0].1.clone(),
                pair[0].2.clone(),
                pair[1].2.clone(),
            )
        })
        .collect::<Vec<_>>();
    write_rescue_reports(sample, sample_dir, &initial, &current, &pairs, &report)
}

#[derive(Clone, Default)]
struct WorkflowProfile {
    rows: Arc<Mutex<Vec<WorkflowProfileRow>>>,
}

#[derive(Clone)]
struct WorkflowProfileRow {
    sample: String,
    round: u32,
    stage: String,
    wall_ms: u128,
    input_bytes: u64,
    output_bytes: u64,
    status: &'static str,
}

fn profile_path_size(path: &Path) -> u64 {
    if path.is_file() {
        return fs::metadata(path)
            .map(|metadata| metadata.len())
            .unwrap_or(0);
    }
    if path.is_dir() {
        return directory_size(path).unwrap_or(0);
    }
    0
}

fn record_profile_named(
    profile: Option<&WorkflowProfile>,
    sample: &str,
    stage: &str,
    started: Instant,
    input_bytes: u64,
    output: &Path,
    result: &Result<(), String>,
) {
    let Some(profile) = profile else { return };
    profile
        .rows
        .lock()
        .expect("workflow profile poisoned")
        .push(WorkflowProfileRow {
            sample: sample.into(),
            round: 0,
            stage: stage.into(),
            wall_ms: started.elapsed().as_millis(),
            input_bytes,
            output_bytes: profile_path_size(output),
            status: if result.is_ok() { "ok" } else { "failed" },
        });
}

fn record_profile_stage(
    profile: Option<&WorkflowProfile>,
    sample: &Sample,
    stage: &str,
    started: Instant,
    input_bytes: u64,
    output: &Path,
    result: &Result<(), String>,
) {
    record_profile_named(
        profile,
        &sample.name,
        stage,
        started,
        input_bytes,
        output,
        result,
    );
}

fn run_profiled_action<F>(
    profile: Option<&WorkflowProfile>,
    sample: &str,
    stage: &str,
    input: &Path,
    output: &Path,
    action: F,
) -> Result<(), String>
where
    F: FnOnce() -> Result<(), String>,
{
    let input_bytes = profile_path_size(input);
    let started = Instant::now();
    let result = action();
    record_profile_named(
        profile,
        sample,
        stage,
        started,
        input_bytes,
        output,
        &result,
    );
    result
}

#[allow(clippy::too_many_arguments)]
fn run_profiled(
    profile: Option<&WorkflowProfile>,
    sample: &Sample,
    stage: &str,
    input: &Path,
    output: &Path,
    bins: &Path,
    binary: &str,
    args: &[String],
) -> Result<(), String> {
    run_profiled_action(profile, &sample.name, stage, input, output, || {
        run(bins, binary, args)
    })
}

fn execute_uce(
    opt: &Options,
    bins: &Path,
    sample: &Sample,
    profile: Option<&WorkflowProfile>,
    filter_threads: usize,
    assembler_threads: usize,
) -> Result<(), String> {
    let sample_dir = Path::new(&opt.output).join(&sample.name);
    if opt.commands.iter().any(|c| c == "filter") {
        let args = uce_filter_args(opt, sample, &sample_dir, filter_threads);
        run_profiled(
            profile,
            sample,
            "filter",
            Path::new(&sample.read1),
            &sample_dir,
            bins,
            "uce_filter",
            &args,
        )?;
        if opt.uce_recruit_mode == "auto" {
            execute_uce_auto_recruit(opt, bins, sample, profile, &sample_dir, filter_threads)?;
        }
    }
    if opt.commands.iter().any(|c| c == "assemble") {
        let args = uce_assembler_args(opt, &sample_dir, assembler_threads)?;
        run_profiled(
            profile,
            sample,
            "assemble",
            &sample_dir.join("filtered"),
            &sample_dir,
            bins,
            "main_assembler-rust",
            &args,
        )?;
        if opt.uce_recruit_mode == "auto" {
            let started = Instant::now();
            let input_bytes = profile_path_size(&sample_dir.join("results"));
            let result = execute_uce_fallback_probe_gate(opt, &sample_dir, assembler_threads);
            record_profile_stage(
                profile,
                sample,
                "fallback-probe-gate",
                started,
                input_bytes,
                &sample_dir,
                &result,
            );
            result?;
        }
        if opt.rescue {
            let rescue_input_bytes = profile_path_size(&sample_dir);
            let started = Instant::now();
            let result = execute_uce_rescue(opt, bins, sample, &sample_dir);
            record_profile_stage(
                profile,
                sample,
                "rescue",
                started,
                rescue_input_bytes,
                &sample_dir,
                &result,
            );
            result?;
        }
    }
    Ok(())
}

fn archive_fallback_probe_rejected(sample_dir: &Path, locus: &str) -> Result<(), String> {
    let archive = sample_dir.join("fallback_probe_rejected");
    for directory in ["results", "contigs_all", "contigs_all_low"] {
        let source = sample_dir.join(directory).join(format!("{locus}.fasta"));
        if !source.is_file() {
            continue;
        }
        let destination_dir = archive.join(directory);
        fs::create_dir_all(&destination_dir).map_err(|error| error.to_string())?;
        let destination = destination_dir.join(format!("{locus}.fasta"));
        if destination.exists() {
            return Err(format!(
                "refusing to overwrite archived UCE fallback candidate '{}'",
                destination.display()
            ));
        }
        fs::rename(&source, &destination).map_err(|error| {
            format!(
                "Unable to archive rejected UCE fallback candidate '{}' as '{}': {error}",
                source.display(),
                destination.display()
            )
        })?;
    }
    Ok(())
}

fn execute_uce_fallback_probe_gate(
    opt: &Options,
    sample_dir: &Path,
    workers: usize,
) -> Result<(), String> {
    let recruit_audit = sample_dir.join("uce_recruit_passes.tsv");
    let fallback_loci = uce_recruit::fallback_recruited_loci(&recruit_audit)?;
    let mut summary = read_uce_summary(&sample_dir.join("uce_assembly_summary.csv"))?;
    let mut contigs = BTreeMap::new();
    for locus in &fallback_loci {
        if !uce_row_accepted(summary.rows.get(locus)) {
            continue;
        }
        let path = sample_dir.join("results").join(format!("{locus}.fasta"));
        let Some(sequence) = first_fasta_sequence(&path)? else {
            continue;
        };
        contigs.insert(locus.clone(), sequence);
    }
    let inverted_repeat_minimum = raw_number::<usize>(
        &opt.raw,
        &["--uce-rescue-inverted-repeat-min-bp"],
        "150",
        "--uce-rescue-inverted-repeat-min-bp",
    )?;
    let gate_workers = if contigs.is_empty() {
        0
    } else {
        workers.max(1).min(contigs.len())
    };
    let evaluation_started = Instant::now();
    let mut evaluated = uce_recruit::evaluate_contig_probe_support_parallel(
        Path::new(&opt.reference),
        &contigs,
        workers,
    )?
    .into_iter()
    .map(|row| (row.locus.clone(), row))
    .collect::<BTreeMap<_, _>>();
    let evaluation_seconds = evaluation_started.elapsed().as_secs_f64();
    let mut evidence = Vec::with_capacity(fallback_loci.len());
    for locus in &fallback_loci {
        if !uce_row_accepted(summary.rows.get(locus)) {
            evidence.push(uce_recruit::ContigProbeEvidence::unavailable(
                locus,
                "assembler_rejected",
            ));
            continue;
        }
        let Some(sequence) = contigs.get(locus) else {
            evidence.push(uce_recruit::ContigProbeEvidence::unavailable(
                locus,
                "accepted_contig_missing",
            ));
            continue;
        };
        let mut row = evaluated
            .remove(locus)
            .ok_or_else(|| format!("missing UCE provisional core evidence for '{locus}'"))?;
        if row.accepted {
            let reads = read_locus_fastq(sample_dir, locus)?;
            row.apply_structure_checks(
                rescue_qc::has_long_inverted_repeat(sequence, inverted_repeat_minimum),
                rescue_qc::maximum_unsupported_internal_gap(sequence, &reads),
            );
        }
        evidence.push(row);
    }
    uce_recruit::write_contig_probe_audit(
        &sample_dir.join("uce_recruit_contig_probe_gate.tsv"),
        &evidence,
    )?;
    for field in [
        "auto_recruit_probe_gate",
        "auto_recruit_probe_gate_reason",
        "auto_recruit_core_anchor_status",
        "auto_recruit_core_structural_review",
        "auto_recruit_core_long_inverted_repeat",
        "auto_recruit_core_maximum_unsupported_internal_gap",
    ] {
        if !summary.headers.iter().any(|header| header == field) {
            summary.headers.push(field.to_owned());
        }
    }
    let mut anchored = 0_usize;
    let mut review = 0_usize;
    let mut probe_rejected = 0_usize;
    let mut structure_rejected = 0_usize;
    let mut assembler_rejected = 0_usize;
    for row in &evidence {
        let summary_row = summary.rows.get_mut(&row.locus).ok_or_else(|| {
            format!(
                "UCE provisional core gate has no assembly summary row for '{}'",
                row.locus
            )
        })?;
        let was_assembler_accepted = uce_row_accepted(Some(summary_row));
        summary_row.insert(
            "auto_recruit_probe_gate".into(),
            if row.accepted { "pass" } else { "reject" }.into(),
        );
        summary_row.insert("auto_recruit_probe_gate_reason".into(), row.reason.into());
        summary_row.insert(
            "auto_recruit_core_anchor_status".into(),
            row.core_anchor_status.into(),
        );
        summary_row.insert(
            "auto_recruit_core_structural_review".into(),
            row.structural_review.into(),
        );
        summary_row.insert(
            "auto_recruit_core_long_inverted_repeat".into(),
            u8::from(row.long_inverted_repeat).to_string(),
        );
        summary_row.insert(
            "auto_recruit_core_maximum_unsupported_internal_gap".into(),
            row.maximum_unsupported_internal_gap.to_string(),
        );
        if row.accepted {
            anchored += 1;
            review += usize::from(row.core_anchor_status == "anchored_with_review");
            continue;
        }
        match row.core_anchor_status {
            "structure_rejected" => structure_rejected += 1,
            "probe_rejected" => probe_rejected += 1,
            _ => assembler_rejected += 1,
        }
        if !was_assembler_accepted {
            continue;
        }
        summary_row.insert("accepted".into(), "0".into());
        summary_row.insert("low_quality".into(), "1".into());
        summary_row.insert("status".into(), "fallback_core_gate_rejected".into());
        archive_fallback_probe_rejected(sample_dir, &row.locus)?;
    }
    write_uce_summary(&sample_dir.join("uce_assembly_summary.csv"), &summary)?;
    write_result_dict_from_uce_summary(sample_dir, &summary)?;
    eprintln!(
        "UCE auto provisional core gate: {} anchored ({} internal-gap review), {} probe rejected, {} structure rejected, {} assembler rejected or missing; {} worker(s), {:.3}s probe evaluation",
        anchored,
        review,
        probe_rejected,
        structure_rejected,
        assembler_rejected,
        gate_workers,
        evaluation_seconds,
    );
    Ok(())
}

fn execute_uce_auto_recruit(
    opt: &Options,
    bins: &Path,
    sample: &Sample,
    profile: Option<&WorkflowProfile>,
    sample_dir: &Path,
    filter_threads: usize,
) -> Result<(), String> {
    let summary_path = sample_dir.join("uce_filter_summary.tsv");
    let fast_selected = uce_recruit::read_selected_fragments(&summary_path)?;
    let unresolved = uce_recruit::unresolved_loci(&fast_selected);
    let fast_pass = RecruitPass::fast(&opt.kf, &opt.step);
    let fallback_pass = RecruitPass::fallback(
        &opt.uce_fallback_kmer_size,
        &opt.uce_fallback_step,
        &opt.uce_fallback_verify_kmer_size,
        &opt.uce_fallback_min_alignment_overlap,
        &opt.uce_fallback_min_alignment_identity,
    );
    uce_recruit::preserve_summary(
        &summary_path,
        &sample_dir.join("uce_filter_summary.fast.tsv"),
    )?;
    if unresolved.is_empty() {
        uce_recruit::write_recruit_audit(
            &sample_dir.join("uce_recruit_passes.tsv"),
            &fast_pass,
            &fallback_pass,
            &fast_selected,
            &BTreeMap::new(),
            &unresolved,
            &BTreeSet::new(),
        )?;
        return Ok(());
    }

    let auto_root = sample_dir.join(".uce_recruit_auto");
    if auto_root.exists() {
        fs::remove_dir_all(&auto_root).map_err(|error| error.to_string())?;
    }
    fs::create_dir_all(&auto_root).map_err(|error| error.to_string())?;
    let retain_loci_file = auto_root.join("unresolved_loci.txt");
    uce_recruit::write_locus_allowlist(&retain_loci_file, &unresolved)?;
    let recruit_reference = auto_root.join("recruit_references");
    uce_recruit::materialize_recruit_reference_subset(
        Path::new(&opt.reference),
        &recruit_reference,
        &unresolved,
    )?;
    let fallback_dir = auto_root.join("fallback");
    let mut args = uce_filter_args_for_pass(
        opt,
        sample,
        &fallback_pass,
        &UceRecruitInvocation {
            sample_dir: &fallback_dir,
            verify_reference: Path::new(&opt.reference),
            recruit_reference: &recruit_reference,
            role: "bait",
            retain_loci_file: Some(&retain_loci_file),
        },
    );
    let thread_index = args
        .iter()
        .position(|argument| argument == "--threads")
        .expect("UCEFilter arguments include --threads");
    args[thread_index + 1] = filter_threads.to_string();
    run_profiled(
        profile,
        sample,
        "filter-fallback",
        Path::new(&sample.read1),
        &fallback_dir,
        bins,
        "uce_filter",
        &args,
    )?;
    let fallback_summary_path = fallback_dir.join("uce_filter_summary.tsv");
    let fallback_selected = uce_recruit::read_selected_fragments(&fallback_summary_path)?;
    let recovered = uce_recruit::merge_fallback_outputs(
        sample_dir,
        &fallback_dir,
        &unresolved,
        &fallback_selected,
    )?;
    uce_recruit::preserve_summary(
        &fallback_summary_path,
        &sample_dir.join("uce_filter_summary.fallback.tsv"),
    )?;
    uce_recruit::write_recruit_audit(
        &sample_dir.join("uce_recruit_passes.tsv"),
        &fast_pass,
        &fallback_pass,
        &fast_selected,
        &fallback_selected,
        &unresolved,
        &recovered,
    )?;
    fs::remove_dir_all(&auto_root).map_err(|error| error.to_string())?;
    eprintln!(
        "UCE auto recruit: {} unresolved fast-pass loci, {} recovered by fallback",
        unresolved.len(),
        recovered.len()
    );
    Ok(())
}

fn execute_uce_legacy(
    opt: &Options,
    bins: &Path,
    sample: &Sample,
    dictionary: &Path,
    profile: Option<&WorkflowProfile>,
    assembler_threads: usize,
) -> Result<(), String> {
    let sample_dir = Path::new(&opt.output).join(&sample.name);
    let candidates = sample_dir.join("filtered_pe");
    if opt.commands.iter().any(|command| command == "filter") {
        let mut args = vec![
            "-r".into(),
            opt.reference.clone(),
            "-q1".into(),
            sample.read1.clone(),
        ];
        if let Some(read2) = &sample.read2 {
            args.extend(["-q2".into(), read2.clone()]);
        }
        args.extend([
            "-o".into(),
            sample_dir.display().to_string(),
            "-kf".into(),
            opt.kf.clone(),
            "-s".into(),
            opt.step.clone(),
            "-gr".into(),
            "-subdir".into(),
            "filtered_pe".into(),
            "-m".into(),
            "5".into(),
            "-lb".into(),
            "-lkd".into(),
            dictionary.display().to_string(),
        ]);
        if opt.max_reads != "0" {
            args.extend(["-m_reads".into(), opt.max_reads.clone()]);
        }
        run_profiled(
            profile,
            sample,
            "filter",
            Path::new(&sample.read1),
            &sample_dir,
            bins,
            "MainFilterNew",
            &args,
        )?;
    }
    if opt.commands.iter().any(|command| command == "refilter") {
        if !candidates.is_dir() {
            return Err("No successful filter run, cannot re-filter".into());
        }
        let args = vec![
            "-r".into(),
            opt.reference.clone(),
            "-qd".into(),
            candidates.display().to_string(),
            "-o".into(),
            sample_dir.join("filtered").display().to_string(),
            "-kf".into(),
            opt.kf.clone(),
            "-p".into(),
            assembler_threads.to_string(),
            "--log-file".into(),
            sample_dir.join("log.txt").display().to_string(),
            "--min-depth".into(),
            opt.low_depth.clone(),
            "--max-depth".into(),
            opt.depth_limit.clone(),
            "--max-size".into(),
            opt.size_limit.clone(),
            "--use-gm2-format".into(),
            "--keep-linked-mates".into(),
        ];
        run_profiled(
            profile,
            sample,
            "refilter",
            &candidates,
            &sample_dir.join("filtered"),
            bins,
            "main_refilter_new",
            &args,
        )?;
    }
    if opt.commands.iter().any(|command| command == "assemble") {
        let args = uce_assembler_args(opt, &sample_dir, assembler_threads)?;
        run_profiled(
            profile,
            sample,
            "assemble",
            &sample_dir.join("filtered"),
            &sample_dir,
            bins,
            "main_assembler-rust",
            &args,
        )?;
        if opt.rescue {
            let rescue_input_bytes = profile_path_size(&sample_dir);
            let started = Instant::now();
            let result = execute_uce_rescue(opt, bins, sample, &sample_dir);
            record_profile_stage(
                profile,
                sample,
                "rescue",
                started,
                rescue_input_bytes,
                &sample_dir,
                &result,
            );
            result?;
        }
    }
    Ok(())
}

fn execute_gene(
    opt: &Options,
    bins: &Path,
    sample: &Sample,
    dictionary: &Path,
    profile: Option<&WorkflowProfile>,
) -> Result<(), String> {
    let sample_dir = Path::new(&opt.output).join(&sample.name);
    let filtered_reads = sample_dir.join("filtered_pe");
    if opt.commands.iter().any(|c| c == "filter") {
        let mut args = vec![
            "-r".into(),
            opt.reference.clone(),
            "-q1".into(),
            sample.read1.clone(),
        ];
        if let Some(read2) = &sample.read2 {
            args.extend(["-q2".into(), read2.clone()]);
        }
        args.extend([
            "-o".into(),
            sample_dir.display().to_string(),
            "-kf".into(),
            opt.kf.clone(),
            "-s".into(),
            opt.step.clone(),
            "-gr".into(),
            "-subdir".into(),
            "filtered_pe".into(),
            "-m".into(),
            if sample.read2.is_some() { "5" } else { "0" }.into(),
            "-lb".into(),
            "-lkd".into(),
            dictionary.display().to_string(),
        ]);
        let rad_linked = opt
            .raw
            .iter()
            .any(|argument| argument == "--rad-internal-linked-recruitment");
        if rad_linked {
            args.extend([
                "--link-rad-arms".into(),
                "--link-rad-max-fragments".into(),
                value(&opt.raw, &["--rad-link-max-fragments"], "256")?,
            ]);
        }
        if let Some(fallback) = optional_value(&opt.raw, &["--rad-fallback-kmers"])? {
            let primary = opt
                .kf
                .parse::<usize>()
                .map_err(|_| "-kf must be an integer")?;
            let mut previous = primary;
            for item in fallback.split(",").filter(|item| !item.is_empty()) {
                let k = item
                    .parse::<usize>()
                    .map_err(|_| "--rad-fallback-kmers must be comma-separated integers")?;
                if k < 16 || k >= previous {
                    return Err(
                        "--rad-fallback-kmers must be >=16 and strictly decrease from -kf".into(),
                    );
                }
                args.extend(["--fallback-kmer".into(), k.to_string()]);
                previous = k;
            }
        }
        run_profiled(
            profile,
            sample,
            "filter",
            Path::new(&sample.read1),
            &sample_dir,
            bins,
            "MainFilterNew",
            &args,
        )?;
    }
    if opt.commands.iter().any(|c| c == "refilter") {
        let input_flag = if sample.read2.is_some() { "-qd" } else { "-qs" };
        let mut args = vec![
            "-r".into(),
            opt.reference.clone(),
            input_flag.into(),
            filtered_reads.display().to_string(),
            "-o".into(),
            sample_dir.join("filtered").display().to_string(),
            "-kf".into(),
            opt.kf.clone(),
            "-p".into(),
            "1".into(),
            "--log-file".into(),
            sample_dir.join("log.txt").display().to_string(),
            "--min-depth".into(),
            opt.low_depth.clone(),
            "--max-depth".into(),
            opt.depth_limit.clone(),
            "--max-size".into(),
            opt.size_limit.clone(),
        ];
        if sample.read2.is_some() {
            args.push("--use-gm2-format".into());
            if opt
                .raw
                .iter()
                .any(|argument| argument == "--rad-internal-linked-recruitment")
            {
                args.push("--keep-linked-mates".into());
            }
        }
        run_profiled(
            profile,
            sample,
            "refilter",
            &filtered_reads,
            &sample_dir.join("filtered"),
            bins,
            "main_refilter_new",
            &args,
        )?;
    }
    if opt.commands.iter().any(|c| c == "assemble") {
        let implementation = value(&opt.raw, &["--assembler-implementation"], "auto")?;
        let binary =
            match implementation.as_str() {
                "auto" | "original" | "original-rust" => "main_assembler-original-rust",
                "uce-rust" => "main_assembler-rust",
                _ => return Err(
                    "--assembler-implementation must be auto, uce-rust, original, or original-rust"
                        .into(),
                ),
            };
        let mut args = vec![
            "-r".into(),
            opt.reference.clone(),
            "-o".into(),
            sample_dir.display().to_string(),
            "-ka".into(),
            opt.ka.clone(),
            "-k_min".into(),
            opt.min_ka.clone(),
            "-k_max".into(),
            opt.max_ka.clone(),
            "-limit_count".into(),
            opt.error_threshold.clone(),
            "-iteration".into(),
            opt.search_depth.clone(),
            "-sb".into(),
            soft_boundary(&opt.soft_boundary)?,
            "-cov_min".into(),
            opt.min_coverage.clone(),
            "-p".into(),
            "1".into(),
        ];
        if implementation == "uce-rust" {
            args.extend([
                "--assembly-mode".into(),
                "original".into(),
                "--assembler-read-chunk-size".into(),
                value(&opt.raw, &["--assembler-read-chunk-size"], "8192")?,
            ]);
        }
        if implementation != "original" {
            if let Some(cache) = assembler_cache_directory(opt)? {
                fs::create_dir_all(&cache).map_err(|e| e.to_string())?;
                args.extend([
                    "--assembler-reference-cache-dir".into(),
                    cache.display().to_string(),
                ]);
            }
        }
        run_profiled(
            profile,
            sample,
            "assemble",
            &sample_dir.join("filtered"),
            &sample_dir,
            bins,
            binary,
            &args,
        )?;
    }
    if opt.commands.iter().any(|c| c == "gene") {
        let classify_input_bytes = profile_path_size(&sample_dir.join("contigs_all"));
        let started = Instant::now();
        let result = run(
            bins,
            "gene_workflow",
            &[
                "classify".into(),
                "--reference".into(),
                opt.reference.clone(),
                "--contigs".into(),
                sample_dir.join("contigs_all").display().to_string(),
                "--sample".into(),
                sample.name.clone(),
                "--out".into(),
                Path::new(&opt.output).join("gene").display().to_string(),
            ],
        );
        record_profile_stage(
            profile,
            sample,
            "gene-classify",
            started,
            classify_input_bytes,
            &Path::new(&opt.output).join("gene"),
            &result,
        );
        result?;
    }
    Ok(())
}

fn optional_value(args: &[String], names: &[&str]) -> Result<Option<String>, String> {
    let value = value(args, names, "")?;
    Ok((!value.is_empty()).then_some(value))
}

fn reference_loci(reference: &Path) -> Result<Vec<(String, PathBuf)>, String> {
    let entries = fs::read_dir(reference).map_err(|e| {
        format!(
            "Unable to read reference directory '{}': {e}",
            reference.display()
        )
    })?;
    let mut loci = entries
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.is_file())
        .filter_map(|path| {
            let extension = path.extension()?.to_str()?.to_ascii_lowercase();
            if matches!(extension.as_str(), "fa" | "fas" | "fasta") {
                let name = path.file_stem()?.to_str()?.to_owned();
                Some((name, path))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    loci.sort_by(|a, b| a.1.cmp(&b.1));
    Ok(loci)
}

fn reference_cache_directory(opt: &Options) -> Result<Option<PathBuf>, String> {
    let configured = optional_value(&opt.raw, &["--reference-cache-dir"])?;
    if configured.is_some() && !opt.reuse_reference_cache {
        return Err("--reference-cache-dir requires --reuse-reference-cache".into());
    }
    if !opt.reuse_reference_cache {
        return Ok(None);
    }
    Ok(Some(configured.map(PathBuf::from).unwrap_or_else(|| {
        Path::new(&opt.output).join(".gm2_reference_cache")
    })))
}

fn reference_dictionary_path(opt: &Options) -> Result<PathBuf, String> {
    let Some(cache) = reference_cache_directory(opt)? else {
        return Ok(Path::new(&opt.output).join(format!("kmer_dict_k{}.dict", opt.kf)));
    };
    let mut digest = Sha256::new();
    digest.update(opt.reference.as_bytes());
    digest.update([0]);
    digest.update(opt.kf.as_bytes());
    digest.update([0]);
    digest.update(opt.step.as_bytes());
    for (_, path) in reference_loci(Path::new(&opt.reference))? {
        digest.update([0]);
        digest.update(
            path.file_name()
                .and_then(|name| name.to_str())
                .unwrap_or_default()
                .as_bytes(),
        );
        digest.update([0]);
        let metadata = fs::metadata(&path).map_err(|e| e.to_string())?;
        digest.update(metadata.len().to_le_bytes());
        digest.update(fs::read(&path).map_err(|e| e.to_string())?);
    }
    let hex = format!("{:x}", digest.finalize());
    Ok(cache.join(format!(
        "reference_kmer_k{}_s{}_{}.dict",
        opt.kf,
        opt.step,
        &hex[..16]
    )))
}

fn assembler_cache_directory(opt: &Options) -> Result<Option<PathBuf>, String> {
    Ok(reference_cache_directory(opt)?.map(|root| root.join("assembler")))
}

fn fastx_output_extension(read: &str) -> &'static str {
    let path = Path::new(read);
    let extension = path
        .extension()
        .and_then(|value| value.to_str())
        .unwrap_or_default();
    let extension = if extension.eq_ignore_ascii_case("gz") {
        path.file_stem()
            .and_then(|value| Path::new(value).extension())
            .and_then(|value| value.to_str())
            .unwrap_or_default()
    } else {
        extension
    };
    if extension.eq_ignore_ascii_case("fq") || extension.eq_ignore_ascii_case("fastq") {
        ".fq"
    } else {
        ".fasta"
    }
}

fn run_program(program: &Path, args: &[String]) -> Result<(), String> {
    let status = Command::new(program)
        .args(args)
        .status()
        .map_err(|e| format!("Unable to run {}: {e}", program.display()))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{} exited with {status}", program.display()))
    }
}

fn run_program_in(program: &Path, args: &[String], directory: &Path) -> Result<(), String> {
    let status = Command::new(program)
        .args(args)
        .current_dir(directory)
        .status()
        .map_err(|e| format!("Unable to run {}: {e}", program.display()))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{} exited with {status}", program.display()))
    }
}

fn first_fasta_sequence(path: &Path) -> Result<Option<String>, String> {
    let content = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut sequence = String::new();
    let mut active = false;
    for line in content.lines() {
        if line.starts_with('>') {
            if active {
                break;
            }
            active = true;
        } else if active {
            sequence.push_str(line.trim());
        }
    }
    Ok((!sequence.is_empty()).then_some(sequence))
}

#[derive(Clone, Debug, Default)]
struct UceSummary {
    headers: Vec<String>,
    rows: std::collections::BTreeMap<String, std::collections::BTreeMap<String, String>>,
}

fn read_uce_summary(path: &Path) -> Result<UceSummary, String> {
    if !path.is_file() {
        return Ok(UceSummary::default());
    }
    let content = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut lines = content.lines();
    let headers = lines
        .next()
        .map(|line| line.split(',').map(str::to_owned).collect::<Vec<_>>())
        .unwrap_or_default();
    let Some(locus_index) = headers.iter().position(|field| field == "locus") else {
        return Ok(UceSummary {
            headers,
            rows: Default::default(),
        });
    };
    let mut rows = std::collections::BTreeMap::new();
    for line in lines.filter(|line| !line.trim().is_empty()) {
        let values = line.split(',').collect::<Vec<_>>();
        let Some(locus) = values.get(locus_index).filter(|value| !value.is_empty()) else {
            continue;
        };
        let row = headers
            .iter()
            .enumerate()
            .map(|(index, header)| {
                (
                    header.clone(),
                    values.get(index).copied().unwrap_or_default().to_owned(),
                )
            })
            .collect();
        rows.insert((*locus).to_owned(), row);
    }
    Ok(UceSummary { headers, rows })
}

fn uce_row_accepted(row: Option<&std::collections::BTreeMap<String, String>>) -> bool {
    let Some(row) = row else {
        return false;
    };
    match row
        .get("accepted")
        .map(|value| value.trim().to_ascii_lowercase())
    {
        Some(value) if !value.is_empty() => matches!(value.as_str(), "1" | "true" | "yes"),
        _ => {
            row.get("status").is_some_and(|status| status == "success")
                && !row.get("low_quality").is_some_and(|value| {
                    matches!(
                        value.trim().to_ascii_lowercase().as_str(),
                        "1" | "true" | "yes"
                    )
                })
        }
    }
}

fn uce_accepted_loci(path: &Path) -> Result<Option<std::collections::HashSet<String>>, String> {
    let summary = read_uce_summary(&path.join("uce_assembly_summary.csv"))?;
    Ok(Some(
        summary
            .rows
            .iter()
            .filter(|(_, row)| uce_row_accepted(Some(row)))
            .map(|(locus, _)| locus.clone())
            .collect(),
    ))
}

fn uce_number(row: Option<&std::collections::BTreeMap<String, String>>, key: &str) -> Option<i64> {
    row.and_then(|row| row.get(key))
        .and_then(|value| value.parse().ok())
}

fn terminal_rescue_loci(
    before: &UceSummary,
    after: &UceSummary,
) -> std::collections::BTreeSet<String> {
    after
        .rows
        .iter()
        .filter_map(|(locus, row)| {
            if !uce_row_accepted(Some(row)) {
                return None;
            }
            let previous = before.rows.get(locus);
            let length_gain = uce_number(Some(row), "selected_contig_length")
                .zip(uce_number(previous, "selected_contig_length"))
                .is_some_and(|(next, prior)| next - prior >= 30);
            let read_gain = uce_number(Some(row), "unique_read_count")
                .zip(uce_number(previous, "unique_read_count"))
                .is_some_and(|(next, prior)| next - prior >= 2);
            (previous.is_none() || !uce_row_accepted(previous) || length_gain || read_gain)
                .then(|| locus.clone())
        })
        .collect()
}

fn write_uce_summary(path: &Path, summary: &UceSummary) -> Result<(), String> {
    if summary.headers.is_empty() {
        return Ok(());
    }
    let mut text = String::new();
    text.push_str(&summary.headers.join(","));
    text.push('\n');
    for row in summary.rows.values() {
        text.push_str(
            &summary
                .headers
                .iter()
                .map(|field| row.get(field).cloned().unwrap_or_default())
                .collect::<Vec<_>>()
                .join(","),
        );
        text.push('\n');
    }
    fs::write(path, text).map_err(|e| e.to_string())
}

fn write_combined_locus(
    locus: &str,
    input_dir: &str,
    output: &Path,
    samples: &[Sample],
    uce: bool,
) -> Result<bool, String> {
    let mut records = String::new();
    for sample in samples {
        let sample_dir = output.join(&sample.name);
        if uce && uce_accepted_loci(&sample_dir)?.is_some_and(|accepted| !accepted.contains(locus))
        {
            continue;
        }
        let source = sample_dir.join(input_dir).join(format!("{locus}.fasta"));
        if !source.is_file() {
            continue;
        }
        if let Some(sequence) = first_fasta_sequence(&source)? {
            records.push('>');
            records.push_str(&sample.name);
            records.push('\n');
            records.push_str(&sequence);
            records.push('\n');
        }
    }
    if records.is_empty() {
        return Ok(false);
    }
    fs::write(
        output
            .join("combined_results")
            .join(format!("{locus}.fasta")),
        records,
    )
    .map_err(|e| e.to_string())?;
    Ok(true)
}

fn alignment_filter(raw: &[String]) -> Result<String, String> {
    if flag(raw, "--no-trimal")? {
        Ok("none".into())
    } else {
        value(raw, &["--alignment-filter"], "trimal")
    }
}

fn phylogeny_binary(program: &str) -> Result<PathBuf, String> {
    let (env_name, default) = match program {
        "raxmlng" => ("GM2_RAXMLNG", "raxml-ng"),
        "iqtree" => ("GM2_IQTREE", "iqtree"),
        "veryfasttree" => ("GM2_VERYFASTTREE", "VeryFastTree"),
        "fasttree" => ("GM2_FASTTREE", "FastTree"),
        _ => {
            return Err("--phylo-program must be raxmlng, iqtree, fasttree, or veryfasttree".into())
        }
    };
    Ok(PathBuf::from(
        env::var(env_name).unwrap_or_else(|_| default.into()),
    ))
}

fn build_tree(
    program: &str,
    binary: &Path,
    input: &Path,
    bootstrap: usize,
    threads: usize,
    quiet: bool,
) -> Result<PathBuf, String> {
    let input_text = input.display().to_string();
    let output = match program {
        "raxmlng" => format!("{input_text}.raxml.bestTree"),
        "iqtree" => format!("{input_text}.treefile"),
        "veryfasttree" => format!("{input_text}.veryfasttree.tre"),
        _ => format!("{input_text}.fasttree.tre"),
    };
    let output = PathBuf::from(output);
    if output.exists() {
        fs::remove_file(&output).map_err(|e| e.to_string())?;
    }
    let mut command = Command::new(binary);
    match program {
        "raxmlng" => {
            command.args([
                "--msa",
                &input_text,
                "--msa-format",
                "FASTA",
                "--model",
                "GTR+G",
                "--redo",
            ]);
            if bootstrap > 0 {
                command.args(["--all", "--bs-trees", &bootstrap.to_string()]);
            } else {
                command.arg("--search");
            }
            if threads > 1 {
                command.args([
                    "--threads",
                    &format!("auto{{{threads}}}"),
                    "--workers",
                    "auto",
                ]);
            } else {
                command.args(["--threads", "1"]);
            }
        }
        "iqtree" => {
            command.args(["-s", &input_text, "-redo"]);
            if bootstrap > 0 {
                command.args(["-B", &bootstrap.to_string()]);
            }
            if threads > 1 {
                command.args(["-T", "AUTO", "-ntmax", &threads.to_string()]);
            } else {
                command.args(["-T", "1"]);
            }
        }
        "veryfasttree" => {
            command.args(["-out", &output.display().to_string(), "-gtr"]);
            if bootstrap > 0 {
                command.args(["-boot", &bootstrap.to_string()]);
            } else {
                command.arg("-nosupport");
            }
            if threads > 1 {
                command.args(["-threads", &threads.to_string()]);
            }
            command.args(["-nt", &input_text]);
        }
        _ => {
            command.args(["-out", &output.display().to_string(), "-gtr"]);
            if bootstrap > 0 {
                command.args(["-boot", &bootstrap.to_string()]);
            } else {
                command.arg("-nosupport");
            }
            command.args(["-nt", &input_text]);
        }
    }
    if quiet {
        command
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null());
    }
    let status = command
        .status()
        .map_err(|e| format!("Unable to run {}: {e}", binary.display()))?;
    if !status.success() {
        return Err(format!("{} exited with {status}", binary.display()));
    }
    if output.is_file() {
        Ok(output)
    } else {
        Err(format!(
            "{} did not create {}",
            binary.display(),
            output.display()
        ))
    }
}

fn file_sha256(path: &Path) -> Result<String, String> {
    let mut file = fs::File::open(path).map_err(|e| e.to_string())?;
    let mut digest = Sha256::new();
    let mut buffer = [0u8; 65536];
    loop {
        let count = file.read(&mut buffer).map_err(|e| e.to_string())?;
        if count == 0 {
            break;
        }
        digest.update(&buffer[..count]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

fn path_sha256(path: &Path) -> Result<String, String> {
    if path.is_file() {
        return file_sha256(path);
    }
    if !path.is_dir() {
        return Err(format!("{} is not a file or directory", path.display()));
    }
    fn visit(root: &Path, directory: &Path, digest: &mut Sha256) -> Result<(), String> {
        let mut entries = fs::read_dir(directory)
            .map_err(|error| error.to_string())?
            .collect::<Result<Vec<_>, _>>()
            .map_err(|error| error.to_string())?;
        entries.sort_by_key(|entry| entry.file_name());
        for entry in entries {
            let path = entry.path();
            let relative = path.strip_prefix(root).map_err(|error| error.to_string())?;
            digest.update(relative.to_string_lossy().as_bytes());
            if path.is_dir() {
                digest.update(b"/\n");
                visit(root, &path, digest)?;
            } else if path.is_file() {
                digest.update(b"\0");
                digest.update(file_sha256(&path)?.as_bytes());
                digest.update(b"\n");
            }
        }
        Ok(())
    }

    let mut digest = Sha256::new();
    visit(path, path, &mut digest)?;
    Ok(format!("{:x}", digest.finalize()))
}

fn input_identity(path: &Path) -> Result<String, String> {
    let metadata = fs::metadata(path).map_err(|error| error.to_string())?;
    let modified = metadata
        .modified()
        .ok()
        .and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|time| format!("{}.{}", time.as_secs(), time.subsec_nanos()))
        .unwrap_or_else(|| "unknown".into());
    Ok(format!(
        "{};bytes={};modified={modified}",
        fs::canonicalize(path)
            .unwrap_or_else(|_| path.to_path_buf())
            .display(),
        metadata.len()
    ))
}

fn raw_number<T: std::str::FromStr>(
    raw: &[String],
    names: &[&str],
    default: &str,
    label: &str,
) -> Result<T, String> {
    value(raw, names, default)?
        .parse()
        .map_err(|_| format!("{label} must be numeric"))
}

fn mito_reference(opt: &Options, bins: &Path) -> Result<PathBuf, String> {
    let raw = &opt.raw;
    let input = PathBuf::from(value(raw, &["--mito-genbank"], "")?);
    if !input.is_file() {
        return Err("mito requires a readable --mito-genbank".into());
    }
    let flank: usize = raw_number(raw, &["--mito-flank"], "150", "--mito-flank")?;
    let length: usize = raw_number(raw, &["--mito-tile-length"], "1200", "--mito-tile-length")?;
    let step: usize = raw_number(raw, &["--mito-tile-step"], "600", "--mito-tile-step")?;
    if step == 0 || length == 0 || step > length {
        return Err("mito requires 0 < --mito-tile-step <= --mito-tile-length".into());
    }
    let reference = Path::new(&opt.output).join(".gm2_mito_reference");
    if reference.exists() {
        fs::remove_dir_all(&reference).map_err(|e| e.to_string())?;
    }
    run(
        bins,
        "mito_workflow",
        &[
            "prepare-reference".into(),
            "--input".into(),
            input.display().to_string(),
            "--out-dir".into(),
            reference.display().to_string(),
            "--flank".into(),
            flank.to_string(),
            "--tile-length".into(),
            length.to_string(),
            "--tile-step".into(),
            step.to_string(),
        ],
    )?;
    Ok(reference)
}

fn mito_assembler_args(
    opt: &Options,
    reference: &Path,
    sample: &Path,
) -> Result<Vec<String>, String> {
    let raw = &opt.raw;
    let ka = value(raw, &["-ka"], "0")?;
    let ka = if ka == "0" { "31".into() } else { ka };
    let args = vec![
        "-r".into(),
        reference.display().to_string(),
        "-o".into(),
        sample.display().to_string(),
        "-ka".into(),
        ka,
        "-k_min".into(),
        value(raw, &["--min-ka"], "21")?,
        "-k_max".into(),
        value(raw, &["--max-ka"], "51")?,
        "-limit_count".into(),
        value(raw, &["-e", "--error-threshold"], "2")?,
        "-iteration".into(),
        raw_number::<usize>(raw, &["-i", "--search-depth"], "4096", "--search-depth")?
            .max(30000)
            .to_string(),
        "-sb".into(),
        "10000".into(),
        "-cov_min".into(),
        value(raw, &["--min-coverage"], "0")?,
        "-p".into(),
        "1".into(),
        "--assembly-mode".into(),
        "uce".into(),
        "--uce-side-candidates".into(),
        value(raw, &["--uce-side-candidates"], "8")?,
        "--uce-max-contig-length".into(),
        value(raw, &["--uce-max-contig-length"], "0")?,
        "--uce-min-read-density".into(),
        "0".into(),
        "--uce-density-check-min-length".into(),
        value(raw, &["--uce-density-check-min-length"], "1000")?,
        "--uce-max-depth-cv".into(),
        value(raw, &["--uce-max-depth-cv"], "0")?,
        "--uce-max-depth-ratio".into(),
        value(raw, &["--uce-max-depth-ratio"], "0")?,
        "--uce-path-strategy".into(),
        value(&opt.raw, &["--uce-path-strategy"], "backbone")?,
        "--uce-backbone-lookahead".into(),
        value(&opt.raw, &["--uce-backbone-lookahead"], "24")?,
        "--assembler-read-chunk-size".into(),
        value(&opt.raw, &["--assembler-read-chunk-size"], "8192")?,
        "--assembler-kmer-count-threads".into(),
        "1".into(),
        "--assembler-graph-format".into(),
        "gfa".into(),
    ];
    Ok(args)
}

fn mito_recruit_refilter_assemble(
    opt: &Options,
    bins: &Path,
    reference: &Path,
    sample: &Sample,
    sample_dir: &Path,
    dictionary: &Path,
    max_reads: usize,
) -> Result<(), String> {
    let raw = &opt.raw;
    let paired = sample.read2.as_ref().unwrap_or(&sample.read1);
    let candidates = sample_dir.join("filtered_pe");
    if candidates.exists() {
        fs::remove_dir_all(&candidates).map_err(|e| e.to_string())?;
    }
    run(
        bins,
        "MainFilterNew",
        &[
            "-r".into(),
            reference.display().to_string(),
            "-q1".into(),
            sample.read1.clone(),
            "-q2".into(),
            paired.clone(),
            "-o".into(),
            sample_dir.display().to_string(),
            "-kf".into(),
            value(raw, &["-kf"], "31")?,
            "-s".into(),
            value(raw, &["-s", "--step-size"], "4")?,
            "-gr".into(),
            "-subdir".into(),
            "filtered_pe".into(),
            "-m".into(),
            "4".into(),
            "-lb".into(),
            "-lkd".into(),
            dictionary.display().to_string(),
            "-m_reads".into(),
            max_reads.to_string(),
        ],
    )?;
    let collapsed = sample_dir.join("filtered_pe_collapsed");
    if collapsed.exists() {
        fs::remove_dir_all(&collapsed).map_err(|e| e.to_string())?;
    }
    run(
        bins,
        "mito_workflow",
        &[
            "collapse-baits".into(),
            "--input-dir".into(),
            candidates.display().to_string(),
            "--out-dir".into(),
            collapsed.display().to_string(),
            "--output-name".into(),
            "mitochondrion".into(),
        ],
    )?;
    fs::remove_dir_all(&candidates).map_err(|e| e.to_string())?;
    fs::rename(&collapsed, &candidates).map_err(|e| e.to_string())?;
    let filtered = sample_dir.join("filtered");
    if filtered.exists() {
        fs::remove_dir_all(&filtered).map_err(|e| e.to_string())?;
    }
    run(
        bins,
        "main_refilter_new",
        &[
            "-r".into(),
            reference.display().to_string(),
            "-qd".into(),
            candidates.display().to_string(),
            "-o".into(),
            filtered.display().to_string(),
            "-kf".into(),
            value(raw, &["-kf"], "31")?,
            "-p".into(),
            "1".into(),
            "--log-file".into(),
            sample_dir.join("log.txt").display().to_string(),
            "--min-depth".into(),
            value(raw, &["--depth-low-water-mark"], "50")?,
            "--max-depth".into(),
            value(raw, &["--depth-limit"], "768")?,
            "--max-size".into(),
            value(raw, &["--file-size-limit"], "6")?,
        ],
    )?;
    run(
        bins,
        "main_assembler-rust",
        &mito_assembler_args(opt, reference, sample_dir)?,
    )
}

fn copy_tree(source: &Path, destination: &Path) -> Result<(), String> {
    fs::create_dir_all(destination).map_err(|e| e.to_string())?;
    for entry in fs::read_dir(source).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let from = entry.path();
        let to = destination.join(entry.file_name());
        if from.is_dir() {
            copy_tree(&from, &to)?;
        } else {
            fs::copy(&from, &to).map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

fn fasta_records(path: &Path) -> Result<Vec<(String, String)>, String> {
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut records = Vec::new();
    let mut id = String::new();
    let mut sequence = String::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if !id.is_empty() && !sequence.is_empty() {
                records.push((std::mem::take(&mut id), std::mem::take(&mut sequence)));
            }
            id = header
                .split_whitespace()
                .next()
                .unwrap_or("sequence")
                .to_owned();
        } else if !id.is_empty() {
            sequence.push_str(line.trim());
        }
    }
    if !id.is_empty() && !sequence.is_empty() {
        records.push((id, sequence));
    }
    Ok(records)
}

/// Builds the rescue reference at the caller-chosen `rescue_root`, independent
/// of `sample`'s own directory. Keeping it outside `sample` lets the caller
/// back up and later reassemble into `sample` with a cheap directory rename
/// instead of a byte copy, since the rescue reference is never nested inside
/// the tree being moved.
fn build_mito_rescue_reference(
    reference: &Path,
    sample: &Path,
    rescue_root: &Path,
) -> Result<Option<PathBuf>, String> {
    let contigs = sample.join("contigs_all/mitochondrion.fasta");
    if !contigs.is_file() {
        return Ok(None);
    }
    let seeds = fasta_records(&contigs)?;
    if seeds.is_empty() {
        return Ok(None);
    }
    if rescue_root.exists() {
        fs::remove_dir_all(rescue_root).map_err(|e| e.to_string())?;
    }
    copy_tree(reference, rescue_root)?;
    let mut bait = fs::OpenOptions::new()
        .append(true)
        .open(rescue_root.join("mitochondrion.fasta"))
        .map_err(|e| e.to_string())?;
    use std::io::Write;
    for (index, (_, sequence)) in seeds.iter().enumerate() {
        if sequence.len() >= 31 {
            writeln!(bait, ">sample_seed_{index}\n{sequence}").map_err(|e| e.to_string())?;
        }
    }
    Ok(Some(rescue_root.to_path_buf()))
}

/// Moves a directory tree cheaply. `destination` must not already exist.
/// A same-filesystem rename is an O(1) metadata operation, unlike `copy_tree`
/// which reads and rewrites every file's bytes. Falls back to copy-then-delete
/// only if the rename itself fails (for example across a filesystem boundary),
/// so the move still succeeds on unusual mount layouts.
fn move_tree(source: &Path, destination: &Path) -> Result<(), String> {
    if let Some(parent) = destination.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    if fs::rename(source, destination).is_ok() {
        return Ok(());
    }
    copy_tree(source, destination)?;
    fs::remove_dir_all(source).map_err(|e| e.to_string())
}

fn build_mito_dictionary(
    opt: &Options,
    bins: &Path,
    reference: &Path,
    dictionary: &Path,
    output: &Path,
) -> Result<(), String> {
    if let Some(parent) = dictionary.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    run(
        bins,
        "MainFilterNew",
        &[
            "-r".into(),
            reference.display().to_string(),
            "-o".into(),
            output.display().to_string(),
            "-kf".into(),
            value(&opt.raw, &["-kf"], "31")?,
            "-s".into(),
            value(&opt.raw, &["-s", "--step-size"], "4")?,
            "-gr".into(),
            "-lkd".into(),
            dictionary.display().to_string(),
            "-m".into(),
            "2".into(),
        ],
    )
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'R' => 'Y',
            b'Y' => 'R',
            b'S' => 'S',
            b'W' => 'W',
            b'K' => 'M',
            b'M' => 'K',
            b'B' => 'V',
            b'V' => 'B',
            b'D' => 'H',
            b'H' => 'D',
            _ => 'N',
        })
        .collect()
}
fn minimal_rotation(sequence: &str) -> String {
    let bytes = sequence.as_bytes();
    let n = bytes.len();
    if n == 0 {
        return String::new();
    }
    let mut doubled = Vec::with_capacity(n * 2);
    doubled.extend_from_slice(bytes);
    doubled.extend_from_slice(bytes);
    let (mut left, mut right, mut offset) = (0usize, 1usize, 0usize);
    while left < n && right < n && offset < n {
        let (a, b) = (doubled[left + offset], doubled[right + offset]);
        if a == b {
            offset += 1;
            continue;
        }
        if a > b {
            left += offset + 1;
            if left == right {
                left += 1;
            }
        } else {
            right += offset + 1;
            if left == right {
                right += 1;
            }
        }
        offset = 0;
    }
    String::from_utf8(doubled[left.min(right)..left.min(right) + n].to_vec()).expect("DNA UTF-8")
}
fn canonical_circular(sequence: &str) -> String {
    let sequence = sequence.to_ascii_uppercase();
    minimal_rotation(&sequence).min(minimal_rotation(&reverse_complement(&sequence)))
}
fn mito_observation(sample: &Path) -> Result<(String, String), String> {
    let summary = fs::read_to_string(sample.join("mito/mitochondrial_assembly_summary.tsv"))
        .unwrap_or_default();
    let status = summary
        .lines()
        .find_map(|line| line.strip_prefix("status\t"))
        .unwrap_or("missing")
        .to_owned();
    let fasta = sample.join("mito/mitochondrial_assembly.fasta");
    let evidence = if status == "circular" {
        fasta_records(&fasta)?
            .first()
            .map(|(_, s)| canonical_circular(s))
            .unwrap_or_default()
    } else if fasta.is_file() {
        file_sha256(&fasta)?
    } else {
        String::new()
    };
    Ok((status, evidence))
}

fn finalize_mito_sample(
    opt: &Options,
    bins: &Path,
    reference: &Path,
    sample_dir: &Path,
    require_circular: bool,
) -> Result<(), String> {
    let raw = &opt.raw;
    let mut args = vec![
        "finalize".into(),
        "--reference-genome".into(),
        reference
            .join("metadata/mitochondrial_reference.fasta")
            .display()
            .to_string(),
        "--gene-metadata".into(),
        reference
            .join("metadata/mitochondrial_genes.tsv")
            .display()
            .to_string(),
        "--contigs".into(),
        sample_dir
            .join("contigs_all/mitochondrion.fasta")
            .display()
            .to_string(),
        "--paired-reads".into(),
        sample_dir
            .join("filtered/mitochondrion.fq")
            .display()
            .to_string(),
        "--out-dir".into(),
        sample_dir.join("mito").display().to_string(),
        "--minimum-overlap".into(),
        value(raw, &["--mito-min-overlap"], "41")?,
        "--minimum-identity".into(),
        value(raw, &["--mito-min-overlap-identity"], "0.98")?,
        "--terminal-window".into(),
        value(raw, &["--mito-terminal-window"], "500")?,
        "--link-kmer".into(),
        value(raw, &["--mito-link-kmer"], "31")?,
        "--minimum-link-hits".into(),
        value(raw, &["--mito-min-link-hits"], "2")?,
        "--minimum-pair-support".into(),
        value(raw, &["--mito-min-pair-support"], "3")?,
        "--bridge-kmer".into(),
        value(raw, &["--mito-bridge-kmer"], "31")?,
        "--bridge-minimum-depth".into(),
        value(raw, &["--mito-bridge-min-depth"], "2")?,
        "--maximum-bridge".into(),
        value(raw, &["--mito-max-bridge"], "1000")?,
        "--minimum-junction-support".into(),
        value(raw, &["--mito-min-junction-support"], "3")?,
        "--require-circular".into(),
        require_circular.to_string(),
    ];
    let graph = sample_dir.join("assembly_graphs/mitochondrion.gfa");
    if graph.is_file() {
        args.extend(["--graph".into(), graph.display().to_string()]);
    }
    run(bins, "mito_workflow", &args)
}

#[allow(clippy::too_many_arguments)]
fn execute_mito_single_stage(
    opt: &Options,
    bins: &Path,
    samples: &[Sample],
    reference: &Path,
    output: &Path,
    dictionary: &Path,
    max_reads: usize,
    require_circular: bool,
) -> Result<(), String> {
    // The recruitment dictionary depends only on the fixed reference and `-kf`,
    // so it is built once by the caller and reused across every adaptive stage.
    fs::create_dir_all(output).map_err(|e| e.to_string())?;
    let failures = Arc::new(Mutex::new(Vec::new()));
    let queued_samples = samples.to_vec();
    let next = Arc::new(Mutex::new(queued_samples.into_iter()));
    let stage_options = opt.clone();
    let mut handles = Vec::new();
    for _ in 0..opt.workers.min(samples.len()).max(1) {
        let failures = Arc::clone(&failures);
        let next = Arc::clone(&next);
        let bins = bins.to_path_buf();
        let reference = reference.to_path_buf();
        let dictionary = dictionary.to_path_buf();
        let staged = output.to_path_buf();
        let mut stage_opt = stage_options.clone();
        stage_opt.output = staged.display().to_string();
        stage_opt.workers = 1;
        handles.push(thread::spawn(move || loop {
            // Each worker keeps pulling samples until the shared queue is empty,
            // so a cohort with more samples than workers still processes every
            // sample instead of silently dropping those past the worker count.
            let Some(sample) = next.lock().expect("mito queue poisoned").next() else {
                return;
            };
            let sample_dir = staged.join(&sample.name);
            let result = mito_recruit_refilter_assemble(
                &stage_opt,
                &bins,
                &reference,
                &sample,
                &sample_dir,
                &dictionary,
                max_reads,
            )
            .and_then(|_| {
                // Reference-recruited mitochondrial pools are very high
                // coverage, so the first UCE-style pass frequently closes the
                // circle on its own. Finalize once to check: a junction-verified
                // circular assembly is already a complete single molecule and
                // cannot be improved by the seed rescue, whose whole purpose is
                // to join fragments. Skipping the rescue recruit+refilter+
                // assemble in that case removes the most expensive redundant
                // work for the common "closed on the first pass" sample.
                finalize_mito_sample(&stage_opt, &bins, &reference, &sample_dir, false)?;
                if mito_observation(&sample_dir)?.0 == "circular" {
                    return Ok(());
                }
                // The rescue reference is built outside sample_dir (a sibling
                // staging area) so backing up the pre-rescue sample directory
                // can be a rename instead of a byte-for-byte copy: nothing the
                // rescue pass is about to overwrite is read back afterward
                // except on the rare hard-failure path.
                let rescue_stage_root = staged.join(".mito_rescue_stage").join(&sample.name);
                if rescue_stage_root.exists() {
                    fs::remove_dir_all(&rescue_stage_root).map_err(|e| e.to_string())?;
                }
                let rescue_ref_dir = rescue_stage_root.join("assembly_refs");
                let Some(rescue_reference) =
                    build_mito_rescue_reference(&reference, &sample_dir, &rescue_ref_dir)?
                else {
                    // No contig means this adaptive depth simply lacks enough
                    // evidence for rescue. Keep its no_contigs observation and
                    // let the scheduler advance to the next read checkpoint.
                    // A one-pass run still requires a circle and must fail.
                    return if require_circular {
                        finalize_mito_sample(&stage_opt, &bins, &reference, &sample_dir, true)
                    } else {
                        Ok(())
                    };
                };
                let backup = staged.join(".mito_seed_backups").join(&sample.name);
                if backup.exists() {
                    fs::remove_dir_all(&backup).map_err(|e| e.to_string())?;
                }
                move_tree(&sample_dir, &backup)?;
                let rescue_result = (|| {
                    let rescue_dict = rescue_stage_root.join("filter.dict");
                    build_mito_dictionary(
                        &stage_opt,
                        &bins,
                        &rescue_reference,
                        &rescue_dict,
                        &rescue_stage_root,
                    )?;
                    mito_recruit_refilter_assemble(
                        &stage_opt,
                        &bins,
                        &rescue_reference,
                        &sample,
                        &sample_dir,
                        &rescue_dict,
                        max_reads,
                    )
                })();
                if rescue_result.is_err() {
                    if sample_dir.exists() {
                        fs::remove_dir_all(&sample_dir).map_err(|e| e.to_string())?;
                    }
                    move_tree(&backup, &sample_dir)?;
                } else if backup.exists() {
                    fs::remove_dir_all(&backup).map_err(|e| e.to_string())?;
                }
                // Fold the attempted rescue reference into the sample directory
                // so the on-disk layout matches prior behaviour for inspection,
                // whether the rescue attempt succeeded or was rolled back.
                if rescue_stage_root.exists() {
                    let destination = sample_dir.join("mito_rescue_round_1");
                    if destination.exists() {
                        fs::remove_dir_all(&destination).map_err(|e| e.to_string())?;
                    }
                    move_tree(&rescue_stage_root, &destination)?;
                }
                finalize_mito_sample(&stage_opt, &bins, &reference, &sample_dir, require_circular)
            });
            if let Err(error) = result {
                failures
                    .lock()
                    .expect("mito failures poisoned")
                    .push(format!("{}: {error}", sample.name));
            }
        }));
    }
    for handle in handles {
        handle.join().map_err(|_| "Rust mito worker panicked")?;
    }
    let failures = failures.lock().map_err(|_| "mito failures poisoned")?;
    if failures.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} mitochondrial sample(s) failed:\n{}",
            failures.len(),
            failures.join("\n")
        ))
    }
}

fn inferred_ipyrad_loci(params: &Path) -> Result<PathBuf, String> {
    let text = fs::read_to_string(params)
        .map_err(|e| format!("Unable to read ipyrad params '{}': {e}", params.display()))?;
    let values = text
        .lines()
        .filter_map(|line| line.split("##").next())
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect::<Vec<_>>();
    if values.len() < 2 {
        return Err("ipyrad params must contain assembly_name [0] and project_dir [1]".into());
    }
    let project = PathBuf::from(&values[1]);
    let project = if project.is_absolute() {
        project
    } else {
        params
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(project)
    };
    Ok(project
        .join(format!("{}_outfiles", values[0]))
        .join(format!("{}.loci", values[0])))
}

fn rad_loci_input(opt: &Options) -> Result<(PathBuf, String), String> {
    let supplied = value(&opt.raw, &["--ipyrad-loci"], "")?;
    let params = value(&opt.raw, &["--ipyrad-params"], "")?;
    if !params.is_empty() {
        let params_path = PathBuf::from(&params);
        if !params_path.is_file() {
            return Err("--ipyrad-params must name a readable params file".into());
        }
        let executable = value(&opt.raw, &["--ipyrad-executable"], "ipyrad")?;
        let steps = value(&opt.raw, &["--ipyrad-steps"], "1234567")?;
        if steps.is_empty()
            || !steps
                .bytes()
                .all(|byte| byte.is_ascii_digit() && byte != b'0')
        {
            return Err("--ipyrad-steps must be a non-empty sequence of steps 1-7".into());
        }
        if !steps.bytes().all(|byte| matches!(byte, b'1'..=b'7')) {
            return Err("--ipyrad-steps may contain only digits 1 through 7".into());
        }
        let status = Command::new(&executable)
            .arg("-p")
            .arg(&params)
            .arg("-s")
            .arg(&steps)
            .status()
            .map_err(|e| format!("Unable to start ipyrad executable '{executable}': {e}"))?;
        if !status.success() {
            return Err(format!("ipyrad assembly failed with status {status}"));
        }
        let loci = if supplied.is_empty() {
            inferred_ipyrad_loci(&params_path)?
        } else {
            PathBuf::from(supplied)
        };
        if !loci.is_file() {
            return Err(format!(
                "ipyrad completed but no .loci file was found at '{}'; pass --ipyrad-loci explicitly if its output was relocated",
                loci.display()
            ));
        }
        return Ok((loci, format!("ipyrad params={params} steps={steps}")));
    }
    let loci = PathBuf::from(supplied);
    if !loci.is_file() {
        return Err("provide a readable --ipyrad-loci FILE, or --ipyrad-params FILE to assemble raw RAD reads with ipyrad".into());
    }
    Ok((loci, "existing ipyrad .loci".into()))
}

fn build_rad_reference(opt: &Options, bins: &Path, reference: &Path) -> Result<PathBuf, String> {
    let (loci, source) = rad_loci_input(opt)?;
    run(
        bins,
        "rad_workflow",
        &[
            "reference".into(),
            "--loci".into(),
            loci.display().to_string(),
            "--out".into(),
            reference.display().to_string(),
        ],
    )?;
    fs::write(
        reference.join("PROVENANCE.txt"),
        format!("source\t{source}\nloci\t{}\n", loci.display()),
    )
    .map_err(|e| e.to_string())?;
    Ok(reference.to_path_buf())
}

fn execute_rad_probe(opt: &Options, bins: &Path) -> Result<(), String> {
    if opt.commands != ["rad-probe"] {
        return Err("rad-probe cannot be combined with other subcommands".into());
    }
    let root = Path::new(&opt.output);
    fs::create_dir_all(root).map_err(|e| e.to_string())?;
    let reference = root.join("rad_reference");
    if flag(&opt.raw, "--rad-denovo")? {
        if !value(&opt.raw, &["--ipyrad-loci", "--ipyrad-params"], "")?.is_empty() {
            return Err(
                "--rad-denovo cannot be combined with --ipyrad-loci or --ipyrad-params".into(),
            );
        }
        if opt.samples.is_empty() {
            return Err("rad-probe --rad-denovo requires -f paired_rad_samples.tsv".into());
        }
        let samples = read_rad_samples(&opt.samples)?;
        let mut args = vec![
            "denovo".into(),
            "--out".into(),
            reference.display().to_string(),
        ];
        for sample in &samples {
            args.extend([
                "--sample".into(),
                sample.name.clone(),
                "--read1".into(),
                sample.read1.clone(),
                "--read2".into(),
                sample.read2.clone().expect("paired RAD samples validated"),
            ]);
        }
        let options = [
            ("--rad-overhang", "--overhang"),
            ("--rad-overhang-r2", "--overhang-r2"),
            ("--rad-kmer", "--kmer"),
            ("--rad-min-count", "--min-count"),
            ("--rad-min-samples", "--min-samples"),
            ("--rad-min-length", "--min-length"),
            ("--rad-max-arm-distance", "--max-arm-distance"),
        ];
        for (source, target) in options {
            if let Some(value) = optional_value(&opt.raw, &[source])? {
                args.extend([target.into(), value]);
            }
        }
        run(bins, "rad_workflow", &args)?;
        fs::write(
            reference.join("PROVENANCE.txt"),
            "source\tdenovo_candidate_probe\nmode\tcanonical_solid_kmer_paired_arms\n",
        )
        .map_err(|e| e.to_string())?;
    } else {
        if !opt.samples.is_empty() {
            return Err("rad-probe uses no -f unless --rad-denovo is selected".into());
        }
        build_rad_reference(opt, bins, &reference)?;
    }
    Ok(())
}

fn execute_rad_validate(opt: &Options, bins: &Path) -> Result<(), String> {
    if opt.commands != ["rad-validate"] {
        return Err("rad-validate cannot be combined with other subcommands".into());
    }
    if !opt.samples.is_empty() {
        return Err("rad-validate discovers samples from --rad-recovery; do not pass -f".into());
    }
    let reference = value(&opt.raw, &["--rad-probe"], "")?;
    if reference.is_empty() || !Path::new(&reference).join("arms").is_dir() {
        return Err("rad-validate requires --rad-probe DIR containing arms/".into());
    }
    let recovery = value(&opt.raw, &["--rad-recovery"], "")?;
    if recovery.is_empty() || !Path::new(&recovery).is_dir() {
        return Err("rad-validate requires --rad-recovery DIR from a completed rad run".into());
    }
    let mut args = vec![
        "validate".into(),
        "--reference".into(),
        reference,
        "--recovery".into(),
        recovery,
        "--out".into(),
        Path::new(&opt.output)
            .join("rad_validated")
            .display()
            .to_string(),
    ];
    for (source, target) in [
        ("--rad-validate-min-identity", "--min-identity"),
        ("--rad-validate-min-breadth", "--min-breadth"),
        ("--rad-validate-min-delta", "--min-delta"),
    ] {
        if let Some(value) = optional_value(&opt.raw, &[source])? {
            args.extend([target.into(), value]);
        }
    }
    run(bins, "rad_workflow", &args)
}

fn execute_rad(opt: &Options, bins: &Path) -> Result<(), String> {
    if opt.commands != ["rad"] {
        return Err(
            "rad is a complete workflow and cannot be combined with other subcommands".into(),
        );
    }
    let implementation = value(&opt.raw, &["--assembler-implementation"], "auto")?;
    if !matches!(
        implementation.as_str(),
        "auto" | "original" | "original-rust"
    ) {
        return Err(
            "rad requires --assembler-implementation auto, original, or original-rust".into(),
        );
    }
    let min_arm_breadth = value(&opt.raw, &["--rad-min-arm-breadth"], "0.80")?;
    let breadth = min_arm_breadth
        .parse::<f64>()
        .map_err(|_| "--rad-min-arm-breadth must be a number")?;
    if !(0.0..=1.0).contains(&breadth) {
        return Err("--rad-min-arm-breadth must be in [0, 1]".into());
    }
    let root = Path::new(&opt.output);
    fs::create_dir_all(root).map_err(|e| e.to_string())?;
    let provided_reference = value(&opt.raw, &["--rad-probe"], "")?;
    if !provided_reference.is_empty()
        && (!value(&opt.raw, &["--ipyrad-loci"], "")?.is_empty()
            || !value(&opt.raw, &["--ipyrad-params"], "")?.is_empty())
    {
        return Err("rad accepts either --rad-probe or an ipyrad input, not both".into());
    }
    let reference = if provided_reference.is_empty() {
        build_rad_reference(opt, bins, &root.join("rad_reference"))?
    } else {
        let path = PathBuf::from(provided_reference);
        if !path.join("arms").is_dir() {
            return Err("--rad-probe must name a rad_reference directory containing arms/".into());
        }
        path
    };
    let recovery = root.join("rad_recovery");
    fs::create_dir_all(&recovery).map_err(|e| e.to_string())?;
    let samples = read_samples(&opt.samples, &recovery)?;
    let mut stage = opt.clone();
    if stage
        .raw
        .iter()
        .any(|argument| argument == "--rad-linked-recruitment")
    {
        stage.raw.push("--rad-internal-linked-recruitment".into());
    }
    stage.reference = reference.join("arms").display().to_string();
    stage.output = recovery.display().to_string();
    stage.assembly_mode = "original".into();
    stage.commands = vec!["filter".into(), "refilter".into(), "assemble".into()];
    let dictionary = recovery.join(format!("rad_kmer_dict_k{}.dict", stage.kf));
    run(
        bins,
        "MainFilterNew",
        &[
            "-r".into(),
            stage.reference.clone(),
            "-o".into(),
            stage.output.clone(),
            "-kf".into(),
            stage.kf.clone(),
            "-s".into(),
            stage.step.clone(),
            "-gr".into(),
            "-lkd".into(),
            dictionary.display().to_string(),
            "-m".into(),
            "2".into(),
        ],
    )?;
    let failures = Arc::new(Mutex::new(Vec::new()));
    let pending = Arc::new(Mutex::new(samples.clone().into_iter()));
    let mut handles = Vec::new();
    for _ in 0..stage.workers.min(samples.len()).max(1) {
        let pending = Arc::clone(&pending);
        let failures = Arc::clone(&failures);
        let stage = stage.clone();
        let bins = bins.to_path_buf();
        let dictionary = dictionary.clone();
        handles.push(thread::spawn(move || loop {
            let Some(sample) = pending.lock().expect("rad sample queue poisoned").next() else {
                break;
            };
            if let Err(error) = execute_gene(&stage, &bins, &sample, &dictionary, None) {
                failures
                    .lock()
                    .expect("rad failures poisoned")
                    .push(format!("{}: {error}", sample.name));
            }
        }));
    }
    for handle in handles {
        handle.join().map_err(|_| "Rust rad worker panicked")?;
    }
    let failures = failures.lock().map_err(|_| "rad failures poisoned")?;
    if !failures.is_empty() {
        return Err(format!(
            "{} RAD sample(s) failed:\n{}",
            failures.len(),
            failures.join("\n")
        ));
    }
    let mut finalize = vec![
        "finalize".into(),
        "--reference".into(),
        reference.display().to_string(),
        "--recovery".into(),
        recovery.display().to_string(),
        "--out".into(),
        root.join("rad_matrix").display().to_string(),
        "--min-arm-breadth".into(),
        min_arm_breadth,
    ];
    for sample in &samples {
        finalize.extend(["--sample".into(), sample.name.clone()]);
    }
    run(bins, "rad_workflow", &finalize)?;
    if stage.cleanup_intermediates {
        cleanup_native_intermediates(&stage, &samples)?;
    }
    Ok(())
}

fn status_line(observations: &std::collections::BTreeMap<String, (String, String)>) -> String {
    observations
        .iter()
        .map(|(sample, (status, _))| format!("{sample}={status}"))
        .collect::<Vec<_>>()
        .join(", ")
}

fn is_stable_circular(previous: Option<&(String, String)>, observation: &(String, String)) -> bool {
    observation.0 == "circular" && previous == Some(observation)
}

fn execute_mito(opt: &Options, bins: &Path, samples: &[Sample]) -> Result<(), String> {
    if opt.commands != ["mito"] {
        return Err(
            "mito is a complete workflow and cannot be combined with other subcommands".into(),
        );
    }
    let reference = mito_reference(opt, bins)?;
    let initial: usize = raw_number(
        &opt.raw,
        &["--mito-initial-reads"],
        "10",
        "--mito-initial-reads",
    )?;
    let maximum: usize = raw_number(&opt.raw, &["--mito-max-reads"], "320", "--mito-max-reads")?;
    if initial == 0 || maximum < initial {
        return Err("--mito-max-reads must be at least --mito-initial-reads".into());
    }
    let root = Path::new(&opt.output);
    fs::create_dir_all(root).map_err(|e| e.to_string())?;
    // Build the reference k-mer dictionary once. It depends only on the fixed
    // reference and `-kf`, so every adaptive stage reuses this single build
    // instead of recomputing an identical dictionary each time.
    let dictionary = root.join(format!(
        "mito_kmer_dict_k{}.dict",
        value(&opt.raw, &["-kf"], "31")?
    ));
    build_mito_dictionary(opt, bins, &reference, &dictionary, root)?;
    if flag(&opt.raw, "--no-mito-adaptive-stop")? {
        return execute_mito_single_stage(
            opt,
            bins,
            samples,
            &reference,
            root,
            &dictionary,
            initial,
            true,
        );
    }
    let stages = root.join(".mito_adaptive");
    // `previous` holds each sample's observation from the prior depth; `settled`
    // holds only verified circles that were identical at two consecutive
    // checkpoints. A partial or no-contig observation is never an early-stop
    // condition: it remains pending and receives the next, larger read budget.
    let mut previous: std::collections::BTreeMap<String, (String, String)> =
        std::collections::BTreeMap::new();
    let mut settled: std::collections::BTreeMap<String, (String, String)> =
        std::collections::BTreeMap::new();
    let freeze_sample = |stage: &Path,
                         settled: &mut std::collections::BTreeMap<String, (String, String)>,
                         name: &str,
                         observation: (String, String)|
     -> Result<(), String> {
        let destination = root.join(name);
        if destination.exists() {
            fs::remove_dir_all(&destination).map_err(|e| e.to_string())?;
        }
        copy_tree(&stage.join(name), &destination)?;
        settled.insert(name.to_string(), observation);
        Ok(())
    };
    let mut limit = initial;
    loop {
        let pending: Vec<Sample> = samples
            .iter()
            .filter(|sample| !settled.contains_key(&sample.name))
            .cloned()
            .collect();
        let stage = stages.join(format!("{limit}m"));
        if stage.exists() {
            fs::remove_dir_all(&stage).map_err(|e| e.to_string())?;
        }
        execute_mito_single_stage(
            opt,
            bins,
            &pending,
            &reference,
            &stage,
            &dictionary,
            limit,
            false,
        )?;
        for sample in &pending {
            let observation = mito_observation(&stage.join(&sample.name))?;
            // Only the same verified circle at two consecutive checkpoints is
            // stable enough to freeze. Identical non-circular results must keep
            // advancing through the adaptive read budgets.
            if is_stable_circular(previous.get(&sample.name), &observation) {
                freeze_sample(&stage, &mut settled, &sample.name, observation)?;
            } else {
                previous.insert(sample.name.clone(), observation);
            }
        }
        let reached_max = limit >= maximum;
        if settled.len() == samples.len() {
            // Every sample produced the same verified circle at two consecutive
            // depths, so the cohort can stop before the maximum read budget.
            return Ok(());
        }
        if reached_max {
            // The read budget is spent with at least one sample still changing:
            // keep each unsettled sample's latest partial and fail, since deeper
            // stability could not be confirmed within the budget.
            let mut report = settled.clone();
            for sample in samples {
                if settled.contains_key(&sample.name) {
                    continue;
                }
                let observation = previous
                    .get(&sample.name)
                    .cloned()
                    .unwrap_or_else(|| ("missing".into(), String::new()));
                freeze_sample(&stage, &mut settled, &sample.name, observation.clone())?;
                report.insert(sample.name.clone(), observation);
            }
            let statuses = status_line(&report);
            return Err(format!(
                "mito adaptive stop did not confirm a stable circular assembly by {limit}M reads; {statuses}"
            ));
        }
        limit = (limit.saturating_mul(2)).min(maximum);
    }
}

fn profile_cache_key(paths: &[&str], kmer: &str) -> Result<String, String> {
    let mut digest = Sha256::new();
    digest.update(kmer.as_bytes());
    for path in paths.iter().filter(|path| !path.is_empty()) {
        let resolved = fs::canonicalize(path).unwrap_or_else(|_| PathBuf::from(path));
        digest.update(b"\0");
        digest.update(resolved.as_os_str().as_encoded_bytes());
        if resolved.is_file() {
            let mut file = fs::File::open(&resolved).map_err(|e| e.to_string())?;
            let mut buffer = [0u8; 65536];
            loop {
                let n = file.read(&mut buffer).map_err(|e| e.to_string())?;
                if n == 0 {
                    break;
                }
                digest.update(&buffer[..n]);
            }
        }
    }
    Ok(format!("{:x}", digest.finalize())[..16].to_owned())
}

fn materialize_profile_reference(opt: &Options) -> Result<(PathBuf, PathBuf), String> {
    let input = PathBuf::from(&opt.reference);
    if !input.is_file() {
        return Err("profiling requires -r to be exactly one marker .fa/.fasta file".into());
    }
    let extension = input
        .extension()
        .and_then(|value| value.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    if !matches!(extension.as_str(), "fa" | "fasta") {
        return Err("profiling reference must use the .fa or .fasta extension".into());
    }
    let directory = Path::new(&opt.output).join(".marker_profile_reference");
    fs::create_dir_all(&directory).map_err(|e| e.to_string())?;
    let target = directory.join(
        input
            .file_name()
            .ok_or("invalid profiling reference filename")?,
    );
    if target.exists() {
        fs::remove_file(&target).map_err(|e| e.to_string())?;
    }
    if fs::hard_link(&input, &target).is_err() {
        fs::copy(&input, &target).map_err(|e| e.to_string())?;
    }
    Ok((target, directory))
}

fn execute_profiling(opt: &Options, bins: &Path, samples: &[Sample]) -> Result<(), String> {
    if opt.commands != ["profiling"] {
        return Err(
            "profiling is a complete marker workflow and cannot be combined with other subcommands"
                .into(),
        );
    }
    let raw = &opt.raw;
    let kmer = value(raw, &["--profile-kmer-size"], "21")?;
    let kmer_number = kmer
        .parse::<usize>()
        .map_err(|_| "--profile-kmer-size must be an odd integer from 15 to 31")?;
    let threshold = value(raw, &["--profile-pseudoalign-threshold"], "0.8")?
        .parse::<f64>()
        .map_err(|_| "--profile-pseudoalign-threshold must be a number")?;
    let relevant = value(raw, &["--profile-relevant-kmer-fraction"], "0.5")?
        .parse::<f64>()
        .map_err(|_| "--profile-relevant-kmer-fraction must be a number")?;
    let memory = value(raw, &["--profile-index-memory-gb"], "2")?
        .parse::<usize>()
        .map_err(|_| "--profile-index-memory-gb must be an integer")?;
    let valid_parameters = (15..=31).contains(&kmer_number)
        && !kmer_number.is_multiple_of(2)
        && 0.0 < threshold
        && threshold <= 1.0
        && (0.0..=1.0).contains(&relevant)
        && memory > 0;
    if !valid_parameters {
        return Err("invalid profiling parameters".into());
    }
    let group = optional_value(raw, &["--profile-group-map"])?;
    if group
        .as_ref()
        .is_some_and(|path| !Path::new(path).is_file())
    {
        return Err("--profile-group-map must be a readable TSV file".into());
    }
    let decoy = optional_value(raw, &["--profile-decoy"])?;
    let (reference, reference_dir) = materialize_profile_reference(opt)?;
    let themisto = optional_value(raw, &["--profile-themisto"])?
        .or_else(|| env::var("GM2_THEMISTO").ok())
        .unwrap_or_else(|| "themisto".into());
    let key = profile_cache_key(
        &[
            &reference.display().to_string(),
            group.as_deref().unwrap_or(""),
            decoy.as_deref().unwrap_or(""),
            &themisto,
        ],
        &kmer,
    )?;
    let cache_root = optional_value(raw, &["--profile-index-dir"])?
        .or_else(|| {
            optional_value(raw, &["--reference-cache-dir"])
                .ok()
                .flatten()
        })
        .unwrap_or_else(|| {
            Path::new(&opt.output)
                .join(".gm2_reference_cache")
                .display()
                .to_string()
        });
    let cache = Path::new(&cache_root).join(format!("profile_themisto_k{kmer}_{key}"));
    let failures = Arc::new(Mutex::new(Vec::new()));
    let queued_samples = samples.to_vec();
    let concurrent_samples = opt.workers.min(samples.len()).max(1);
    let profile_uce_memory_limit_mib = resolve_uce_memory_limit_mib(concurrent_samples);
    eprintln!(
        "Auto UCEFilter memory limit: {profile_uce_memory_limit_mib} MiB per sample ({concurrent_samples} concurrent profiling job(s))"
    );
    let next = Arc::new(Mutex::new(queued_samples.into_iter()));
    let mut handles = Vec::new();
    for _ in 0..concurrent_samples {
        let failures = Arc::clone(&failures);
        let next = Arc::clone(&next);
        let reference_dir = reference_dir.clone();
        let reference = reference.clone();
        let cache = cache.clone();
        let themisto = themisto.clone();
        let group = group.clone();
        let decoy = decoy.clone();
        let kmer = kmer.clone();
        let output_root = opt.output.clone();
        let step = opt.step.clone();
        let low_depth = opt.low_depth.clone();
        let depth_limit = opt.depth_limit.clone();
        let size_limit = opt.size_limit.clone();
        let max_reads = opt.max_reads.clone();
        let uce_memory_limit_mib = profile_uce_memory_limit_mib;
        let bins = bins.to_path_buf();
        let force = flag(raw, "--profile-force-rebuild")?;
        handles.push(thread::spawn(move || {
            let Some(sample) = next.lock().expect("profiling queue poisoned").next() else {
                return;
            };
            let sample_dir = Path::new(&output_root).join(&sample.name);
            let mut filter = vec![
                "-r".into(),
                reference_dir.display().to_string(),
                "--recruit-references".into(),
                reference_dir.display().to_string(),
                "-q1".into(),
                sample.read1.clone(),
            ];
            if let Some(read2) = &sample.read2 {
                filter.extend(["-q2".into(), read2.clone()]);
            }
            filter.extend([
                "-o".into(),
                sample_dir.display().to_string(),
                "-kf".into(),
                kmer.clone(),
                "-s".into(),
                step,
                "--selection".into(),
                "auto".into(),
                "--reference-role".into(),
                "bait".into(),
                "--threads".into(),
                "1".into(),
                "--memory-limit-mib".into(),
                uce_memory_limit_mib.to_string(),
                "--min-depth".into(),
                low_depth,
                "--max-depth".into(),
                depth_limit,
                "--max-size".into(),
                size_limit,
            ]);
            if max_reads != "0" {
                filter.extend(["--max-fragments".into(), max_reads]);
            }
            let result = run(&bins, "uce_filter", &filter).and_then(|_| {
                let filtered = sample_dir.join("filtered");
                let reads = fs::read_dir(&filtered)
                    .map_err(|e| e.to_string())?
                    .filter_map(Result::ok)
                    .map(|entry| entry.path())
                    .filter(|path| path.is_file())
                    .filter(|path| {
                        matches!(
                            path.extension()
                                .and_then(|x| x.to_str())
                                .unwrap_or_default()
                                .to_ascii_lowercase()
                                .as_str(),
                            "fq" | "fastq" | "fasta" | "fa"
                        )
                    })
                    .collect::<Vec<_>>();
                if reads.len() != 1 {
                    return Err("profiling requires exactly one merged recruited-read file".into());
                }
                let profile = sample_dir.join("marker_profile");
                if profile.exists() {
                    fs::remove_dir_all(&profile).map_err(|e| e.to_string())?;
                }
                let mut args = vec![
                    "--reference".into(),
                    reference.display().to_string(),
                    "--reads".into(),
                    reads[0].display().to_string(),
                    "--output".into(),
                    profile.display().to_string(),
                    "--cache".into(),
                    cache.display().to_string(),
                    "--themisto".into(),
                    themisto,
                    "--threads".into(),
                    "1".into(),
                    "--kmer-size".into(),
                    kmer,
                    "--threshold".into(),
                    threshold.to_string(),
                    "--relevant-kmer-fraction".into(),
                    relevant.to_string(),
                    "--index-memory-gb".into(),
                    memory.to_string(),
                ];
                if let Some(group) = group {
                    args.extend(["--groups".into(), group]);
                }
                if let Some(decoy) = decoy {
                    args.extend(["--decoy".into(), decoy]);
                }
                if force {
                    args.push("--force-rebuild".into());
                }
                run(&bins, "marker_profile", &args)?;
                if profile.join("marker_reference_support.tsv").is_file() {
                    Ok(())
                } else {
                    Err("profiling failed to produce marker_reference_support.tsv".into())
                }
            });
            if let Err(error) = result {
                failures
                    .lock()
                    .expect("profiling failures poisoned")
                    .push(format!("{}: {error}", sample.name));
            }
        }));
    }
    for handle in handles {
        handle
            .join()
            .map_err(|_| "Rust profiling worker panicked")?;
    }
    let failures = failures.lock().map_err(|_| "profiling failures poisoned")?;
    if failures.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} sample(s) failed:\n{}",
            failures.len(),
            failures.join("\n")
        ))
    }
}

fn execute_gene_annotate(opt: &Options, bins: &Path) -> Result<(), String> {
    let raw = &opt.raw;
    let input = value(raw, &["--gene-input"], "")?;
    let proteins = value(raw, &["--gene-protein-reference"], "")?;
    if !Path::new(&input).is_dir() {
        return Err("--gene-input must be a gene output directory".into());
    }
    if !Path::new(&proteins).is_dir() {
        return Err("gene-annotate requires --gene-protein-reference".into());
    }
    run(
        bins,
        "gene_workflow",
        &[
            "annotate".into(),
            "--input".into(),
            input,
            "--protein-reference".into(),
            proteins,
            "--out".into(),
            opt.output.clone(),
            "--miniprot".into(),
            value(raw, &["--gene-miniprot"], "miniprot")?,
            "--threads".into(),
            opt.workers.to_string(),
        ],
    )
}

fn execute_gene_resolve(opt: &Options, bins: &Path) -> Result<(), String> {
    let raw = &opt.raw;
    let input = value(raw, &["--gene-input"], "")?;
    if !Path::new(&input).is_dir() {
        return Err("--gene-input must be an annotation directory".into());
    }
    let mut args = vec![
        "resolve".into(),
        "--input".into(),
        input,
        "--out".into(),
        opt.output.clone(),
        "--mafft".into(),
        value(raw, &["--gene-mafft"], "mafft")?,
        "--iqtree".into(),
        value(raw, &["--gene-iqtree"], "iqtree")?,
        "--threads".into(),
        opt.workers.to_string(),
        "--min-taxa".into(),
        value(raw, &["--gene-min-taxa"], "4")?,
        "--min-aa-length".into(),
        value(raw, &["--gene-min-aa-length"], "30")?,
        "--min-effective-codon-sites".into(),
        value(raw, &["--gene-min-effective-codon-sites"], "30")?,
    ];
    if let Some(path) = optional_value(raw, &["--gene-outgroup"])? {
        if !Path::new(&path).is_file() {
            return Err("--gene-outgroup must be a readable file".into());
        }
        args.extend(["--outgroup".into(), path]);
    }
    let ufboot = value(raw, &["--gene-ufboot"], "0")?;
    if ufboot != "0" {
        args.extend(["--ufboot".into(), ufboot]);
    }
    if let Some(path) = optional_value(raw, &["--gene-taper"])? {
        if !Path::new(&path).is_file() {
            return Err("--gene-taper must be a readable correction_multi.jl script".into());
        }
        args.extend([
            "--taper-script".into(),
            path,
            "--julia".into(),
            value(raw, &["--gene-julia"], "julia")?,
        ]);
    }
    run(bins, "gene_workflow", &args)
}

fn execute_gene_tree(opt: &Options) -> Result<(), String> {
    let raw = &opt.raw;
    let input = PathBuf::from(value(raw, &["--gene-input"], "")?);
    if !input.is_dir() {
        return Err("--gene-input must be a gene-resolve output directory".into());
    }
    let mode = value(raw, &["--gene-species-mode"], "strict")?;
    if !matches!(mode.as_str(), "strict" | "multicopy") {
        return Err("--gene-species-mode must be strict or multicopy".into());
    }
    let (trees, mapping, output_name) = if mode == "strict" {
        (
            input.join("astral_input/resolved_1to1.trees"),
            None,
            "gene_strict_aster.tree",
        )
    } else {
        (
            input.join("astralpro_input/multicopy.trees"),
            Some(input.join("astralpro_input/leaf_to_species.tsv")),
            "gene_multicopy_aster.tree",
        )
    };
    if !trees.is_file() || fs::metadata(&trees).map_err(|e| e.to_string())?.len() == 0 {
        return Err(format!(
            "No usable {mode} gene trees found: {}",
            trees.display()
        ));
    }
    if let Some(mapping) = &mapping {
        if !mapping.is_file() {
            return Err(format!("Missing multicopy leaf map: {}", mapping.display()));
        }
    }
    fs::create_dir_all(&opt.output).map_err(|e| e.to_string())?;
    let aster = PathBuf::from(value(raw, &["--gene-aster"], "astral")?);
    let output = Path::new(&opt.output).join(output_name);
    let log = Path::new(&opt.output).join(format!("{output_name}.log"));
    if output.exists() {
        fs::remove_file(&output).map_err(|e| e.to_string())?;
    }
    let mut args = vec![
        "-i".into(),
        trees.display().to_string(),
        "-o".into(),
        output.display().to_string(),
        "-t".into(),
        opt.workers.to_string(),
    ];
    if let Some(mapping) = &mapping {
        args.extend(["-a".into(), mapping.display().to_string()]);
    }
    let file = fs::File::create(&log).map_err(|e| e.to_string())?;
    let status = Command::new(&aster)
        .args(&args)
        .stdout(file.try_clone().map_err(|e| e.to_string())?)
        .stderr(file)
        .status()
        .map_err(|e| format!("Cannot find ASTER2 executable: {}: {e}", aster.display()))?;
    if !status.success() {
        return Err(format!(
            "ASTER2 exited with {status}; inspect {}",
            log.display()
        ));
    }
    let tree = fs::read_to_string(&output)
        .map_err(|_| {
            format!(
                "ASTER2 completed without a species tree; inspect {}",
                log.display()
            )
        })?
        .lines()
        .map(str::trim)
        .find(|line| !line.is_empty())
        .unwrap_or("")
        .to_owned();
    if !tree.starts_with('(') || !tree.ends_with(';') {
        return Err(format!(
            "ASTER2 output is not a Newick tree; inspect {}",
            log.display()
        ));
    }
    let mut provenance = format!("field\tvalue\nmode\t{mode}\naster_executable\t{}\ncommand\t{}\ngene_trees\t{}\ngene_trees_sha256\t{}\nspecies_tree\t{}\nspecies_tree_sha256\t{}\n", aster.display(), args.join(" "), trees.display(), file_sha256(&trees)?, output.display(), file_sha256(&output)?);
    if let Some(mapping) = mapping {
        provenance.push_str(&format!(
            "leaf_to_species\t{}\nleaf_to_species_sha256\t{}\n",
            mapping.display(),
            file_sha256(&mapping)?
        ));
    }
    fs::write(
        Path::new(&opt.output).join("gene_tree_provenance.tsv"),
        provenance,
    )
    .map_err(|e| e.to_string())
}

fn execute_tree(opt: &Options) -> Result<(), String> {
    let raw = &opt.raw;
    let method = value(raw, &["-m", "--tree-method"], "coalescent")?;
    if !matches!(method.as_str(), "coalescent" | "concatenation") {
        return Err("--tree-method must be coalescent or concatenation".into());
    }
    let program = value(raw, &["--phylo-program"], "fasttree")?;
    let binary = phylogeny_binary(&program)?;
    let filter = alignment_filter(raw)?;
    if !matches!(filter.as_str(), "trimal" | "alifilter" | "none") {
        return Err("--alignment-filter must be trimal, alifilter, or none".into());
    }
    let output = Path::new(&opt.output);
    if method == "concatenation" {
        let alignment = output.join(if filter == "none" {
            "combined_results.fasta"
        } else {
            "combined_trimed.fasta"
        });
        if !alignment.is_file() {
            return Err(format!(
                "Unable to find the concatenated alignment at '{}'",
                alignment.display()
            ));
        }
        let bootstrap = value(raw, &["-b", "--bootstrap"], "1000")?
            .parse::<usize>()
            .map_err(|_| "--bootstrap must be an integer")?;
        let tree = build_tree(&program, &binary, &alignment, bootstrap, opt.workers, false)?;
        fs::copy(tree, output.join("Concatenation.tree")).map_err(|e| e.to_string())?;
        return Ok(());
    }
    let alignment_dir = output.join(if filter == "none" {
        "combined_results/aligned"
    } else {
        "combined_trimed"
    });
    let loci = reference_loci(Path::new(&opt.reference))?;
    let mut trees = Vec::new();
    let mut failures = Vec::new();
    for (locus, _) in loci {
        let alignment = alignment_dir.join(format!("{locus}.fasta"));
        if !alignment.is_file() {
            continue;
        }
        match build_tree(&program, &binary, &alignment, 0, 1, true) {
            Ok(path) => trees.push(path),
            Err(error) => failures.push((locus, alignment, error)),
        }
    }
    let failure_path = output.join("failed_gene_trees.tsv");
    if failures.is_empty() {
        if failure_path.exists() {
            fs::remove_file(&failure_path).map_err(|e| e.to_string())?;
        }
    } else {
        let mut text = "locus\talignment\terror\n".to_owned();
        for (locus, alignment, error) in failures {
            text.push_str(&format!(
                "{locus}\t{}\t{}\n",
                alignment.display(),
                error.replace(['\t', '\n'], " ")
            ));
        }
        fs::write(failure_path, text).map_err(|e| e.to_string())?;
    }
    trees.sort();
    let mut all = String::new();
    for tree in trees {
        let content = fs::read_to_string(&tree).map_err(|e| e.to_string())?;
        if let Some(line) = content.lines().map(str::trim).find(|line| !line.is_empty()) {
            all.push_str(line);
            all.push('\n');
        }
    }
    if all.is_empty() {
        return Err(
            "Unable to reconstruct coalescent trees because no gene tree is available".into(),
        );
    }
    let trees_path = output.join("combined_genes.trees");
    fs::write(&trees_path, all).map_err(|e| e.to_string())?;
    let coalescent = output.join("Coalescent.tree");
    if coalescent.exists() {
        fs::remove_file(&coalescent).map_err(|e| e.to_string())?;
    }
    let astral = PathBuf::from(env::var("GM2_ASTRAL").unwrap_or_else(|_| "astral".into()));
    run_program(
        &astral,
        &[
            "-i".into(),
            trees_path.display().to_string(),
            "-o".into(),
            coalescent.display().to_string(),
            "-t".into(),
            opt.workers.to_string(),
        ],
    )
}

#[derive(Clone)]
struct PermitPool {
    state: Arc<(Mutex<usize>, Condvar)>,
}

struct Permit {
    state: Arc<(Mutex<usize>, Condvar)>,
}

impl PermitPool {
    fn new(limit: usize) -> Self {
        Self {
            state: Arc::new((Mutex::new(limit), Condvar::new())),
        }
    }

    fn acquire(&self) -> Permit {
        let (available, changed) = &*self.state;
        let mut remaining = available.lock().expect("permit pool poisoned");
        while *remaining == 0 {
            remaining = changed.wait(remaining).expect("permit pool poisoned");
        }
        *remaining -= 1;
        Permit {
            state: Arc::clone(&self.state),
        }
    }
}

impl Drop for Permit {
    fn drop(&mut self) {
        let (available, changed) = &*self.state;
        *available.lock().expect("permit pool poisoned") += 1;
        changed.notify_one();
    }
}

fn execute_combine(
    opt: &Options,
    bins: &Path,
    samples: &[Sample],
    default_source: &str,
) -> Result<(), String> {
    let raw = &opt.raw;
    let source = value(raw, &["-cs", "--combine-source"], default_source)?;
    let input_dir = match source.as_str() {
        "assembly" => "results",
        "consensus" => "consensus",
        "trimmed" => "blast",
        _ => return Err("--combine-source must be assembly, consensus, or trimmed".into()),
    };
    let no_alignment = flag(raw, "--no-alignment")?;
    let filter = if flag(raw, "--no-trimal")? {
        "none".into()
    } else {
        value(raw, &["--alignment-filter"], "trimal")?
    };
    if !matches!(filter.as_str(), "trimal" | "alifilter" | "none") {
        return Err("--alignment-filter must be trimal, alifilter, or none".into());
    }
    let strict = flag(raw, "--strict-combine-errors")?;
    let clean_difference = value(raw, &["-cd", "--clean-difference"], "1")?
        .parse::<f64>()
        .map_err(|_| "--clean-difference must be a number")?;
    let clean_sequences = value(raw, &["-cn", "--clean-sequences"], "0")?
        .parse::<usize>()
        .map_err(|_| "--clean-sequences must be an integer")?;
    if !(0.0..=1.0).contains(&clean_difference) || clean_sequences > samples.len() {
        return Err("invalid combine cleanup thresholds".into());
    }
    let combined = Path::new(&opt.output).join("combined_results");
    if combined.exists() {
        fs::remove_dir_all(&combined).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(&combined).map_err(|e| e.to_string())?;
    let loci = reference_loci(Path::new(&opt.reference))?;
    let names = loci
        .iter()
        .map(|(name, _)| name.clone())
        .collect::<Vec<_>>();
    let uce = opt.assembly_mode == "uce";
    for locus in &names {
        write_combined_locus(locus, input_dir, Path::new(&opt.output), samples, uce)?;
    }
    if no_alignment {
        return Ok(());
    }
    let msa_program = value(raw, &["--msa-program"], "mafft")?;
    if !matches!(msa_program.as_str(), "mafft" | "clustalo") {
        return Err("--msa-program must be mafft or clustalo".into());
    }
    let msa_threads = value(raw, &["--msa-threads"], "1")?
        .parse::<usize>()
        .map_err(|_| "--msa-threads must be an integer")?;
    if msa_threads == 0 {
        return Err("--msa-threads must be at least 1".into());
    }
    if msa_threads > opt.workers {
        return Err("--msa-threads cannot be greater than -p".into());
    }
    let aligned = combined.join("aligned");
    fs::create_dir_all(&aligned).map_err(|e| e.to_string())?;
    let filtered = Path::new(&opt.output).join("combined_trimed");
    if filtered.exists() {
        fs::remove_dir_all(&filtered).map_err(|e| e.to_string())?;
    }
    if filter != "none" {
        fs::create_dir_all(&filtered).map_err(|e| e.to_string())?;
    }
    let mafft = PathBuf::from(env::var("GM2_MAFFT").unwrap_or_else(|_| "mafft".into()));
    let clustalo = PathBuf::from(env::var("GM2_CLUSTALO").unwrap_or_else(|_| "clustalo".into()));
    let trimal = PathBuf::from(env::var("GM2_TRIMAL").unwrap_or_else(|_| "trimal".into()));
    let alifilter = PathBuf::from(env::var("GM2_ALIFILTER").unwrap_or_else(|_| "AliFilter".into()));
    let model = optional_value(raw, &["--alifilter-model"])?;
    if model.is_some() && filter != "alifilter" {
        return Err("--alifilter-model requires --alignment-filter alifilter".into());
    }
    let filter_processes = value(raw, &["--filter-processes"], &opt.workers.to_string())?
        .parse::<usize>()
        .map_err(|_| "--filter-processes must be an integer")?;
    if filter_processes == 0 {
        return Err("--filter-processes must be at least 1".into());
    }
    // Preserve the original scheduler: up to -p loci can make progress, while
    // MSA and column-filter subprocesses have independent resource caps.
    let msa_pool = PermitPool::new((opt.workers / msa_threads).max(1));
    let filter_pool = (filter != "none").then(|| PermitPool::new(filter_processes));
    let pending = Arc::new(Mutex::new(names.into_iter()));
    let failures = Arc::new(Mutex::new(Vec::<(String, String)>::new()));
    thread::scope(|scope| {
        for _ in 0..opt.workers {
            let pending = Arc::clone(&pending);
            let failures = Arc::clone(&failures);
            let combined = combined.clone();
            let aligned = aligned.clone();
            let filtered = filtered.clone();
            let msa_program = msa_program.clone();
            let mafft = mafft.clone();
            let clustalo = clustalo.clone();
            let trimal = trimal.clone();
            let alifilter = alifilter.clone();
            let filter = filter.clone();
            let model = model.clone();
            let msa_pool = msa_pool.clone();
            let filter_pool = filter_pool.clone();
            scope.spawn(move || loop {
                let Some(locus) = pending.lock().expect("combine queue poisoned").next() else {
                    break;
                };
                let input = combined.join(format!("{locus}.fasta"));
                if !input.is_file() {
                    continue;
                }
                let output = aligned.join(format!("{locus}.fasta"));
                let msa_permit = msa_pool.acquire();
                let result = if msa_program == "mafft" {
                    let file = match fs::File::create(&output) {
                        Ok(file) => file,
                        Err(error) => {
                            failures
                                .lock()
                                .expect("combine failures poisoned")
                                .push((locus, error.to_string()));
                            continue;
                        }
                    };
                    let status = Command::new(&mafft)
                        .args([
                            "--auto",
                            "--quiet",
                            "--nuc",
                            "--thread",
                            &msa_threads.to_string(),
                            &input.display().to_string(),
                        ])
                        .stdout(file)
                        .status()
                        .map_err(|error| error.to_string());
                    match status {
                        Ok(status) if status.success() => Ok(()),
                        Ok(status) => Err(format!("mafft exited with {status}")),
                        Err(error) => Err(error),
                    }
                } else {
                    run_program(
                        &clustalo,
                        &[
                            "-i".into(),
                            input.display().to_string(),
                            "-o".into(),
                            output.display().to_string(),
                            "--auto".into(),
                            "--force".into(),
                            "--seqtype=DNA".into(),
                            format!("--threads={msa_threads}"),
                        ],
                    )
                };
                drop(msa_permit);
                let result = result.and_then(|_| {
                    run_program(
                        &bins.join("fix_alignment"),
                        &[
                            "-f".into(),
                            output.display().to_string(),
                            "-n".into(),
                            clean_sequences.to_string(),
                            "-p".into(),
                            clean_difference.to_string(),
                        ],
                    )
                });
                let result = result.and_then(|_| {
                    let filter_permit = filter_pool.as_ref().map(PermitPool::acquire);
                    let filtered_result = if filter == "trimal" {
                        run_program(
                            &trimal,
                            &[
                                "-in".into(),
                                output.display().to_string(),
                                "-out".into(),
                                filtered
                                    .join(format!("{locus}.fasta"))
                                    .display()
                                    .to_string(),
                                "-automated1".into(),
                            ],
                        )
                    } else if filter == "alifilter" {
                        let mut args = vec![
                            "-i".into(),
                            output.display().to_string(),
                            "-o".into(),
                            filtered
                                .join(format!("{locus}.fasta"))
                                .display()
                                .to_string(),
                        ];
                        if let Some(model) = &model {
                            args.extend(["-m".into(), model.clone()]);
                        }
                        run_program(&alifilter, &args)
                    } else {
                        Ok(())
                    };
                    drop(filter_permit);
                    filtered_result
                });
                if let Err(error) = result {
                    let _ = fs::remove_file(&output);
                    failures
                        .lock()
                        .expect("combine failures poisoned")
                        .push((locus, error));
                }
            });
        }
    });
    let mut failures = failures
        .lock()
        .map_err(|_| "combine failures poisoned")?
        .clone();
    failures.sort_by(|left, right| left.0.cmp(&right.0));
    if strict && !failures.is_empty() {
        let (locus, error) = &failures[0];
        return Err(format!("combine failed on {locus}: {error}"));
    }
    for (locus, error) in failures {
        eprintln!("Warning: combine failed on {locus}: {error}");
    }
    run_program(
        &bins.join("merge_seq"),
        &[
            "-input".into(),
            aligned.display().to_string(),
            "-exts".into(),
            ".fasta".into(),
            "-missing".into(),
            "-".into(),
            "-output".into(),
            Path::new(&opt.output)
                .join("combined_results.fasta")
                .display()
                .to_string(),
        ],
    )?;
    if filter != "none" {
        run_program(
            &bins.join("merge_seq"),
            &[
                "-input".into(),
                filtered.display().to_string(),
                "-exts".into(),
                ".fasta".into(),
                "-missing".into(),
                "-".into(),
                "-output".into(),
                Path::new(&opt.output)
                    .join("combined_trimed.fasta")
                    .display()
                    .to_string(),
            ],
        )?;
    }
    Ok(())
}

fn execute_trim(
    opt: &Options,
    bins: &Path,
    samples: &[Sample],
    default_source: &str,
) -> Result<(), String> {
    let raw = &opt.raw;
    let source = value(raw, &["-ts", "--trim-source"], default_source)?;
    if !matches!(source.as_str(), "assembly" | "consensus") {
        return Err("--trim-source must be assembly or consensus".into());
    }
    let mode_name = value(raw, &["-tm", "--trim-mode"], "terminal")?;
    let mode = match mode_name.as_str() {
        "all" => "0",
        "longest" => "1",
        "terminal" => "2",
        "isoform" => "3",
        _ => return Err("--trim-mode must be all, longest, terminal, or isoform".into()),
    };
    let retention = value(raw, &["-tr", "--trim-retention"], "0")?
        .parse::<f64>()
        .map_err(|_| "--trim-retention must be a number")?;
    if !(0.0..=1.0).contains(&retention) {
        return Err("--trim-retention must be in [0, 1]".into());
    }
    let makeblastdb =
        PathBuf::from(env::var("GM2_MAKEBLASTDB").unwrap_or_else(|_| "makeblastdb".into()));
    let blast = PathBuf::from(if mode_name == "isoform" {
        env::var("GM2_MAGICBLAST").unwrap_or_else(|_| "magicblast".into())
    } else {
        env::var("GM2_BLASTN").unwrap_or_else(|_| "blastn".into())
    });
    let database_dir = Path::new(&opt.output).join("blast_db");
    fs::create_dir_all(&database_dir).map_err(|e| e.to_string())?;
    let loci = reference_loci(Path::new(&opt.reference))?;
    for (locus, reference) in &loci {
        run_program_in(
            &makeblastdb,
            &[
                "-in".into(),
                fs::canonicalize(reference)
                    .map_err(|e| e.to_string())?
                    .display()
                    .to_string(),
                "-dbtype".into(),
                "nucl".into(),
                "-out".into(),
                locus.clone(),
            ],
            &database_dir,
        )?;
    }
    let mut tasks = Vec::new();
    for sample in samples {
        let sample_dir = Path::new(&opt.output).join(&sample.name);
        let input = sample_dir.join(if source == "consensus" {
            "consensus"
        } else {
            "results"
        });
        if !input.is_dir() {
            continue;
        }
        let output = sample_dir.join("blast");
        if output.exists() {
            fs::remove_dir_all(&output).map_err(|e| e.to_string())?;
        }
        fs::create_dir_all(&output).map_err(|e| e.to_string())?;
        for (locus, reference) in &loci {
            let query = input.join(format!("{locus}.fasta"));
            if query.is_file() {
                tasks.push((
                    locus.clone(),
                    query,
                    reference.clone(),
                    output.join(format!("{locus}.fasta")),
                ));
            }
        }
    }
    let trim = bins.join("build_trimed");
    let next = Arc::new(Mutex::new(tasks.into_iter()));
    let failures = Arc::new(Mutex::new(Vec::new()));
    let mut handles = Vec::new();
    for _ in 0..opt.workers {
        let next = Arc::clone(&next);
        let failures = Arc::clone(&failures);
        let trim = trim.clone();
        let blast = blast.clone();
        let database_dir = database_dir.clone();
        handles.push(thread::spawn(move || loop {
            let Some((locus, query, reference, output)) =
                next.lock().expect("trim queue poisoned").next()
            else {
                break;
            };
            let result = run_program(
                &trim,
                &[
                    "-i".into(),
                    query.display().to_string(),
                    "-r".into(),
                    reference.display().to_string(),
                    "-o".into(),
                    output.display().to_string(),
                    "-b".into(),
                    database_dir.join(&locus).display().to_string(),
                    "-m".into(),
                    mode.into(),
                    "-p".into(),
                    (retention * 100.0).to_string(),
                    "--executable".into(),
                    blast.display().to_string(),
                ],
            );
            if let Err(error) = result {
                failures
                    .lock()
                    .expect("trim failures poisoned")
                    .push(format!("{locus}: {error}"));
            }
        }));
    }
    for handle in handles {
        handle.join().map_err(|_| "Rust trim worker panicked")?;
    }
    let failures = failures.lock().map_err(|_| "trim failures poisoned")?;
    if failures.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} trim task(s) failed:\n{}",
            failures.len(),
            failures.join("\n")
        ))
    }
}

fn execute_consensus(opt: &Options, bins: &Path, samples: &[Sample]) -> Result<(), String> {
    let threshold = value(&opt.raw, &["-c", "--consensus-threshold"], "0.75")?
        .parse::<f64>()
        .map_err(|_| "--consensus-threshold must be a number")?;
    if !(0.0 < threshold && threshold <= 1.0) {
        return Err("--consensus-threshold must be in (0, 1]".into());
    }
    let loci = reference_loci(Path::new(&opt.reference))?;
    let minimap2 = PathBuf::from(env::var("GM2_MINIMAP2").unwrap_or_else(|_| "minimap2".into()));
    let consensus = bins.join("build_consensus");
    let mut tasks = Vec::new();
    for sample in samples {
        let sample_dir = Path::new(&opt.output).join(&sample.name);
        let results = sample_dir.join("results");
        if !results.is_dir() {
            continue;
        }
        let out = sample_dir.join("consensus");
        if out.exists() {
            fs::remove_dir_all(&out).map_err(|e| e.to_string())?;
        }
        fs::create_dir_all(&out).map_err(|e| e.to_string())?;
        let read_extension = fastx_output_extension(&sample.read1);
        for (locus, assembly) in &loci {
            let assembled = results.join(assembly.file_name().ok_or("invalid reference filename")?);
            let reads = sample_dir
                .join("filtered")
                .join(format!("{locus}{read_extension}"));
            if assembled.is_file() && reads.is_file() {
                tasks.push((assembled, reads, out.join(format!("{locus}.sam"))));
            }
        }
    }
    let failures = Arc::new(Mutex::new(Vec::new()));
    let next = Arc::new(Mutex::new(tasks.into_iter()));
    let workers = opt.workers.min(1.max(loci.len() * samples.len()));
    let mut handles = Vec::new();
    for _ in 0..workers {
        let failures = Arc::clone(&failures);
        let next = Arc::clone(&next);
        let minimap2 = minimap2.clone();
        let consensus = consensus.clone();
        handles.push(thread::spawn(move || loop {
            let Some((assembly, reads, sam)) =
                next.lock().expect("consensus queue poisoned").next()
            else {
                break;
            };
            let mapped = run_program(
                &minimap2,
                &[
                    "-ax".into(),
                    "sr".into(),
                    "-t".into(),
                    "1".into(),
                    "--sam-hit-only".into(),
                    "--secondary=no".into(),
                    "-o".into(),
                    sam.display().to_string(),
                    assembly.display().to_string(),
                    reads.display().to_string(),
                ],
            );
            let result = mapped
                .and_then(|_| {
                    run_program(
                        &consensus,
                        &[
                            "-i".into(),
                            sam.display().to_string(),
                            "-c".into(),
                            threshold.to_string(),
                            "-o".into(),
                            sam.parent()
                                .ok_or("SAM has no parent")?
                                .display()
                                .to_string(),
                            "-s".into(),
                            "0".into(),
                        ],
                    )
                })
                .and_then(|_| fs::remove_file(&sam).map_err(|e| e.to_string()));
            if let Err(error) = result {
                failures
                    .lock()
                    .expect("consensus failures poisoned")
                    .push(format!("{}: {error}", sam.display()));
            }
        }));
    }
    for handle in handles {
        handle
            .join()
            .map_err(|_| "Rust consensus worker panicked")?;
    }
    let failures = failures.lock().map_err(|_| "consensus failures poisoned")?;
    if failures.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} consensus task(s) failed:\n{}",
            failures.len(),
            failures.join("\n")
        ))
    }
}

fn execute_te(opt: &Options, bins: &Path) -> Result<(), String> {
    let raw = &opt.raw;
    let mut args = vec![
        "--samples".into(),
        opt.samples.clone(),
        "--output".into(),
        opt.output.clone(),
        "--stage".into(),
        value(raw, &["--te-stage"], "all")?,
        "--threads".into(),
        opt.workers.to_string(),
        "--kmer".into(),
        value(raw, &["--te-kmer"], "25")?,
        "--min-kmer-count".into(),
        value(raw, &["--te-min-kmer-count"], "8")?,
        "--catalog-pairs".into(),
        value(raw, &["--te-catalog-pairs"], "10000")?,
        "--mainfilter".into(),
        bins.join("MainFilterNew").display().to_string(),
        "--annotation-min-fragment".into(),
        value(raw, &["--te-annotate-min-fragment"], "80")?,
        "--annotation-max-fragment".into(),
        value(raw, &["--te-annotate-max-fragment"], "800")?,
        "--annotation-min-support".into(),
        value(raw, &["--te-annotate-min-support"], "5")?,
        "--annotation-min-identity".into(),
        value(raw, &["--te-annotate-min-identity"], "0.8")?,
        "--annotation-min-coverage".into(),
        value(raw, &["--te-annotate-min-coverage"], "0.6")?,
        "--annotation-min-delta".into(),
        value(raw, &["--te-annotate-min-delta"], "0.1")?,
        "--assemble-min-kmer-count".into(),
        value(raw, &["--te-assemble-min-kmer-count"], "3")?,
        "--assemble-branch-ratio".into(),
        value(raw, &["--te-assemble-branch-ratio"], "1.5")?,
        "--assemble-max-fragments".into(),
        value(raw, &["--te-assemble-max-fragments"], "3")?,
    ];
    if let Some(path) = optional_value(raw, &["--te-read-ledger"])? {
        args.extend(["--read-ledger".into(), path]);
    }
    if let Some(path) = optional_value(raw, &["--te-library"])? {
        args.extend(["--te-library".into(), path]);
    }
    for (public, internal) in [
        ("--te-quantify-pairs", "--quantify-pairs"),
        ("--te-bootstrap-replicates", "--bootstrap-replicates"),
    ] {
        if let Some(setting) = optional_value(raw, &[public])? {
            args.extend([internal.into(), setting]);
        }
    }
    if flag(raw, "--te-estimate-genome-fraction")? {
        args.push("--estimate-genome-fraction".into());
    }
    run(bins, "main_repeat", &args)
}

fn execute_population(opt: &Options, bins: &Path) -> Result<(), String> {
    let raw = &opt.raw;
    let engine = value(raw, &["--engine"], "pseudoref")?;
    if !matches!(engine.as_str(), "pseudoref" | "panref" | "panrefv2") {
        return Err("--engine must be pseudoref, panref, or panrefv2".into());
    }
    if matches!(engine.as_str(), "panref" | "panrefv2") && opt.reference.is_empty() {
        return Err("-r is required with --engine panref or panrefv2".into());
    }
    let mut args = vec![
        "--output".into(),
        opt.output.clone(),
        "--samples".into(),
        opt.samples.clone(),
        "--engine".into(),
        engine.clone(),
        "--reference-strategy".into(),
        value(raw, &["--population-reference-strategy"], "sqcl-longest")?,
        "--start-at".into(),
        value(raw, &["--population-start-at"], "reference")?,
        "--threads".into(),
        opt.workers.to_string(),
        "--min-mapq".into(),
        value(raw, &["--population-min-mapq"], "20")?,
        "--min-baseq".into(),
        value(raw, &["--population-min-baseq"], "20")?,
        "--min-dp".into(),
        value(raw, &["--population-min-dp"], "5")?,
        "--min-gq".into(),
        value(raw, &["--population-min-gq"], "20")?,
        "--min-qual".into(),
        value(raw, &["--population-min-qual"], "20")?,
        "--min-call-rate".into(),
        value(raw, &["--population-min-call-rate"], "0.8")?,
        "--min-mac".into(),
        value(raw, &["--population-min-mac"], "2")?,
        "--ld-window".into(),
        value(raw, &["--population-ld-window"], "50")?,
        "--ld-step".into(),
        value(raw, &["--population-ld-step"], "5")?,
        "--ld-r2".into(),
        value(raw, &["--population-ld-r2"], "0.2")?,
        "--admixture-k-min".into(),
        value(raw, &["--population-admixture-k-min"], "2")?,
        "--admixture-k-max".into(),
        value(raw, &["--population-admixture-k-max"], "6")?,
        "--admixture-cv".into(),
        value(raw, &["--population-admixture-cv"], "10")?,
        "--stop-after".into(),
        value(raw, &["--population-stop-after"], "selection")?,
        "--minibwa".into(),
        value(raw, &["--population-minibwa"], "minibwa")?,
        "--samtools".into(),
        value(raw, &["--population-samtools"], "samtools")?,
        "--bcftools".into(),
        value(raw, &["--population-bcftools"], "bcftools")?,
        "--plink".into(),
        value(raw, &["--population-plink"], "plink")?,
        "--admixture".into(),
        value(raw, &["--population-admixture"], "admixture")?,
    ];
    if matches!(engine.as_str(), "panref" | "panrefv2") {
        args.extend(["--panref-baits".into(), opt.reference.clone()]);
    }
    if flag(raw, "--population-panrefv2-include-low-confidence")? {
        args.push("--panrefv2-include-low-confidence".into());
    }
    if let Some(path) = optional_value(raw, &["--population-reference-fasta"])? {
        args.extend(["--reference-fasta".into(), path]);
    }
    if flag(raw, "--population-skip-mark-duplicates")? {
        args.push("--skip-mark-duplicates".into());
    }
    if flag(raw, "--population-skip-plink")? {
        args.push("--skip-plink".into());
    }
    if flag(raw, "--population-skip-admixture")? {
        args.push("--skip-admixture".into());
    }
    run(bins, "main_population", &args)
}

fn execute_stats(opt: &Options, bins: &Path, samples: &[Sample]) -> Result<(), String> {
    let mut args = vec![
        "--output".into(),
        opt.output.clone(),
        "--reference".into(),
        opt.reference.clone(),
    ];
    for sample in samples {
        args.extend([
            "--sample".into(),
            sample.name.clone(),
            sample.read1.clone(),
            sample.read2.clone().unwrap_or_default(),
        ]);
    }
    if opt.stats_count_input_reads {
        args.push("--count-input-reads".into());
    }
    if opt.stats_no_heatmap {
        args.push("--no-heatmap".into());
    }
    run(bins, "gm2_stats", &args)
}

fn cleanup_native_intermediates(opt: &Options, samples: &[Sample]) -> Result<(), String> {
    validate_cleanup_options(opt)?;
    if !opt.cleanup_intermediates {
        return Ok(());
    }
    let root = fs::canonicalize(&opt.output).map_err(|e| e.to_string())?;
    let mut rows = String::from("path\tbytes\taction\treason\n");
    for sample in samples {
        let sample_dir = root.join(&sample.name);
        for (name, reason) in [
            ("filtered", "reproducible filtered reads"),
            ("filtered_pe", "reproducible filter candidates"),
        ] {
            let path = sample_dir.join(name);
            if path.is_dir() && !path.is_symlink() {
                let bytes = directory_size(&path)?;
                let action = if opt.cleanup_dry_run {
                    "would_remove"
                } else {
                    fs::remove_dir_all(&path).map_err(|e| e.to_string())?;
                    "removed"
                };
                rows.push_str(&format!(
                    "{}\t{bytes}\t{action}\t{reason}\n",
                    path.display()
                ));
            }
        }
    }
    let manifest = if opt.cleanup_dry_run {
        "cleanup_preview.tsv"
    } else {
        "cleanup_manifest.tsv"
    };
    fs::write(root.join(manifest), rows).map_err(|e| e.to_string())
}

fn validate_cleanup_options(opt: &Options) -> Result<(), String> {
    if opt.cleanup_dry_run && !opt.cleanup_intermediates {
        return Err("--cleanup-dry-run requires --cleanup-intermediates".into());
    }
    if !opt.cleanup_intermediates {
        return Ok(());
    }
    if !opt.commands.iter().any(|command| command == "filter")
        || !opt.commands.iter().any(|command| command == "assemble")
    {
        return Err(
            "--cleanup-intermediates requires filter and assemble in the same invocation".into(),
        );
    }
    Ok(())
}

fn directory_size(path: &Path) -> Result<u64, String> {
    let mut bytes = 0;
    for entry in fs::read_dir(path).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let child = entry.path();
        if child.is_symlink() {
            continue;
        }
        if child.is_dir() {
            bytes += directory_size(&child)?;
        } else {
            bytes += fs::metadata(&child).map_err(|e| e.to_string())?.len();
        }
    }
    Ok(bytes)
}

fn write_native_workflow_profile(
    output: &Path,
    profile: &WorkflowProfile,
    elapsed_ms: u128,
) -> Result<(), String> {
    let mut rows = profile
        .rows
        .lock()
        .map_err(|_| "workflow profile poisoned")?
        .clone();
    rows.sort_by(|left, right| {
        (&left.sample, left.round, &left.stage).cmp(&(&right.sample, right.round, &right.stage))
    });
    let mut text =
        String::from("sample\tround\tstage\twall_ms\tinput_bytes\toutput_bytes\tstatus\n");
    for row in rows {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.sample,
            row.round,
            row.stage,
            row.wall_ms,
            row.input_bytes,
            row.output_bytes,
            row.status
        ));
    }
    text.push_str(&format!(
        "__workflow__\t0\tnative_dispatch\t{elapsed_ms}\t0\t0\tok\n"
    ));
    let path = output.join("workflow_profile.tsv");
    let temporary = output.join("workflow_profile.tsv.tmp");
    fs::write(&temporary, text).map_err(|e| e.to_string())?;
    fs::rename(&temporary, path).map_err(|e| e.to_string())
}

fn manifest_value(value: impl std::fmt::Display) -> String {
    value.to_string().replace(['\t', '\n', '\r'], " ")
}

fn manifest_arguments(raw: &[String]) -> String {
    raw.iter()
        // `--resume` changes control flow only; including it would make a
        // completed non-resume invocation impossible to resume.
        .filter(|argument| argument.as_str() != "--resume")
        .cloned()
        .collect::<Vec<_>>()
        .join(" ")
}

fn workflow_manifest_text(opt: &Options, samples: &[Sample]) -> Result<String, String> {
    let reference = Path::new(&opt.reference);
    let sample_list = Path::new(&opt.samples);
    let mut rows = vec![
        "field\tvalue".into(),
        "schema_version\t1".into(),
        format!("tool_version\t{}", env!("CARGO_PKG_VERSION")),
        format!("commands\t{}", manifest_value(opt.commands.join(","))),
        format!("assembly_mode\t{}", manifest_value(&opt.assembly_mode)),
        format!("workers\t{}", opt.workers),
        format!("worker_source\t{}", manifest_value(&opt.worker_source)),
        format!("reference_path\t{}", manifest_value(reference.display())),
        format!("reference_sha256\t{}", path_sha256(reference)?),
        format!(
            "sample_list_path\t{}",
            manifest_value(sample_list.display())
        ),
        format!("sample_list_sha256\t{}", file_sha256(sample_list)?),
        format!("sample_count\t{}", samples.len()),
        format!(
            "raw_arguments\t{}",
            manifest_value(manifest_arguments(&opt.raw))
        ),
    ];
    for sample in samples {
        rows.push(format!(
            "sample.{}.read1\t{}",
            manifest_value(&sample.name),
            manifest_value(input_identity(Path::new(&sample.read1))?)
        ));
        if let Some(read2) = &sample.read2 {
            rows.push(format!(
                "sample.{}.read2\t{}",
                manifest_value(&sample.name),
                manifest_value(input_identity(Path::new(read2))?)
            ));
        }
    }
    Ok(format!("{}\n", rows.join("\n")))
}

fn write_workflow_manifest(opt: &Options, samples: &[Sample]) -> Result<(), String> {
    let output = Path::new(&opt.output);
    let text = workflow_manifest_text(opt, samples)?;
    let path = output.join("workflow_manifest.tsv");
    let temporary = output.join("workflow_manifest.tsv.tmp");
    fs::write(&temporary, text).map_err(|error| error.to_string())?;
    fs::rename(&temporary, path).map_err(|error| error.to_string())
}

fn resume_completed_workflow(opt: &Options, samples: &[Sample]) -> Result<bool, String> {
    if !opt.resume {
        return Ok(false);
    }
    let output = Path::new(&opt.output);
    let existing_manifest = fs::read_to_string(output.join("workflow_manifest.tsv"))
        .map_err(|_| "--resume requires an existing workflow_manifest.tsv".to_string())?;
    if existing_manifest != workflow_manifest_text(opt, samples)? {
        return Err(
            "--resume refused: current inputs or options do not match workflow_manifest.tsv".into(),
        );
    }
    let status = fs::read_to_string(output.join("workflow_status.tsv"))
        .map_err(|_| "--resume requires an existing successful workflow_status.tsv".to_string())?;
    if !status.lines().any(|line| line == "state\tsucceeded") {
        return Err("--resume refused: previous workflow did not complete successfully; rerun without --resume".into());
    }
    Ok(true)
}

fn error_kind(error: &str) -> &'static str {
    if error.contains("Unable to run") || error.contains("exited with") {
        "component"
    } else if error.contains("does not exist")
        || error.contains("must be")
        || error.contains("requires")
        || error.contains("Invalid")
        || error.contains("invalid")
    {
        "input"
    } else if error.contains("permission") || error.contains("No space") || error.contains("I/O") {
        "io"
    } else {
        "workflow"
    }
}

fn json_string(value: impl std::fmt::Display) -> String {
    let mut escaped = String::new();
    for character in value.to_string().chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            character if character.is_control() => {
                escaped.push_str(&format!("\\u{:04x}", character as u32))
            }
            character => escaped.push(character),
        }
    }
    escaped
}

fn write_workflow_status(
    output: &Path,
    commands: &[String],
    result: &Result<(), String>,
) -> Result<(), String> {
    let (state, kind, error) = match result {
        Ok(()) => ("succeeded", "none", String::new()),
        Err(error) => ("failed", error_kind(error), manifest_value(error)),
    };
    let text = format!(
        "field\tvalue\nschema_version\t1\nstate\t{state}\nerror_kind\t{kind}\ncommands\t{}\nerror\t{error}\n",
        manifest_value(commands.join(",")),
    );
    let path = output.join("workflow_status.tsv");
    let temporary = output.join("workflow_status.tsv.tmp");
    fs::write(&temporary, text).map_err(|write_error| write_error.to_string())?;
    fs::rename(&temporary, path).map_err(|rename_error| rename_error.to_string())
}

fn execute_with_status(opt: Options) -> Result<(), String> {
    let output = PathBuf::from(&opt.output);
    let commands = opt.commands.clone();
    let log_format = opt.log_format.clone();
    let resume = opt.resume;
    let result = execute_native(opt);
    // Do not create a new directory merely to report an early validation
    // failure. Once a workflow has created its output root, the status file is
    // an atomic terminal record for users and batch schedulers.
    if !resume && !output.as_os_str().is_empty() && output.is_dir() {
        if let Err(status_error) = write_workflow_status(&output, &commands, &result) {
            return match result {
                Ok(()) => Err(format!("Unable to write workflow status: {status_error}")),
                Err(error) => Err(format!(
                    "{error}\nUnable to write workflow status: {status_error}"
                )),
            };
        }
    }
    if log_format == "json" {
        let (state, kind, error) = match &result {
            Ok(()) => ("succeeded", "none", String::new()),
            Err(error) => ("failed", error_kind(error), error.clone()),
        };
        eprintln!(
            "{{\"event\":\"workflow_finished\",\"state\":\"{state}\",\"error_kind\":\"{kind}\",\"commands\":\"{}\",\"output\":\"{}\",\"error\":\"{}\"}}",
            json_string(commands.join(",")),
            json_string(output.display()),
            json_string(error),
        );
    }
    result
}

fn validate_parallelism(opt: &Options) -> Result<(), String> {
    let msa_threads = value(&opt.raw, &["--msa-threads"], "1")?
        .parse::<usize>()
        .map_err(|_| "--msa-threads must be an integer")?;
    if msa_threads == 0 {
        return Err("--msa-threads must be at least 1".into());
    }
    if msa_threads > opt.workers {
        return Err("--msa-threads cannot be greater than -p".into());
    }
    if let Some(filter_processes) = optional_value(&opt.raw, &["--filter-processes"])? {
        let filter_processes = filter_processes
            .parse::<usize>()
            .map_err(|_| "--filter-processes must be an integer")?;
        if filter_processes == 0 {
            return Err("--filter-processes must be at least 1".into());
        }
    }
    Ok(())
}

fn parse_cpu_list(raw: &str) -> Option<Vec<usize>> {
    let mut cpus = BTreeSet::new();
    for item in raw.trim().split(',').filter(|item| !item.is_empty()) {
        if let Some((start, end)) = item.split_once('-') {
            let start = start.trim().parse::<usize>().ok()?;
            let end = end.trim().parse::<usize>().ok()?;
            if start > end {
                return None;
            }
            cpus.extend(start..=end);
        } else {
            cpus.insert(item.trim().parse::<usize>().ok()?);
        }
    }
    (!cpus.is_empty()).then(|| cpus.into_iter().collect())
}

fn allowed_cpu_ids() -> Option<Vec<usize>> {
    let status = fs::read_to_string("/proc/self/status").ok()?;
    status.lines().find_map(|line| {
        line.strip_prefix("Cpus_allowed_list:")
            .and_then(parse_cpu_list)
    })
}

fn physical_core_count(cpu_ids: &[usize], topology_root: &Path) -> Option<usize> {
    let mut cores = BTreeSet::new();
    for cpu in cpu_ids {
        let topology = topology_root.join(format!("cpu{cpu}/topology"));
        let package = fs::read_to_string(topology.join("physical_package_id")).ok()?;
        let core = fs::read_to_string(topology.join("core_id")).ok()?;
        cores.insert((package.trim().to_owned(), core.trim().to_owned()));
    }
    (!cores.is_empty()).then_some(cores.len())
}

fn parse_cpu_max(raw: &str) -> Option<usize> {
    let mut fields = raw.split_whitespace();
    let quota = fields.next()?;
    if quota == "max" {
        return None;
    }
    let quota = quota.parse::<u64>().ok()?;
    let period = fields.next()?.parse::<u64>().ok()?;
    (period > 0).then(|| (quota / period).max(1) as usize)
}

fn parse_cpu_cfs(quota: &str, period: &str) -> Option<usize> {
    let quota = quota.trim().parse::<i64>().ok()?;
    let period = period.trim().parse::<u64>().ok()?;
    (quota > 0 && period > 0).then(|| ((quota as u64 / period).max(1)) as usize)
}

fn cgroup_tree_cpu_quota(
    root: &Path,
    relative: &Path,
    quota_name: &str,
    period_name: Option<&str>,
) -> Option<usize> {
    let mut directory = root.join(relative);
    let mut limit = None;
    loop {
        let value = if let Some(period_name) = period_name {
            fs::read_to_string(directory.join(quota_name))
                .ok()
                .zip(fs::read_to_string(directory.join(period_name)).ok())
                .and_then(|(quota, period)| parse_cpu_cfs(&quota, &period))
        } else {
            fs::read_to_string(directory.join(quota_name))
                .ok()
                .and_then(|value| parse_cpu_max(&value))
        };
        if let Some(value) = value {
            limit = Some(limit.map_or(value, |prior: usize| prior.min(value)));
        }
        if directory == root || !directory.pop() {
            break;
        }
    }
    limit
}

fn cgroup_cpu_quota() -> Option<usize> {
    let v2 = cgroup_relative_path(None).and_then(|relative| {
        cgroup_tree_cpu_quota(Path::new("/sys/fs/cgroup"), &relative, "cpu.max", None)
    });
    let v1 = cgroup_relative_path(Some("cpu")).and_then(|relative| {
        ["/sys/fs/cgroup/cpu", "/sys/fs/cgroup/cpu,cpuacct"]
            .into_iter()
            .filter_map(|root| {
                cgroup_tree_cpu_quota(
                    Path::new(root),
                    &relative,
                    "cpu.cfs_quota_us",
                    Some("cpu.cfs_period_us"),
                )
            })
            .min()
    });
    match (v2, v1) {
        (Some(v2), Some(v1)) => Some(v2.min(v1)),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

fn scheduler_cpu_limit() -> Option<(String, usize)> {
    ["SLURM_CPUS_PER_TASK", "PBS_NP", "NSLOTS"]
        .into_iter()
        .find_map(|name| {
            let value = env::var(name).ok()?.parse::<usize>().ok()?;
            (value > 0).then(|| (name.to_owned(), value))
        })
}

fn auto_worker_budget() -> (usize, String) {
    let allowed = allowed_cpu_ids();
    let logical = allowed
        .as_ref()
        .map(Vec::len)
        .or_else(|| thread::available_parallelism().ok().map(usize::from))
        .unwrap_or(1)
        .max(1);
    let physical = allowed
        .as_deref()
        .and_then(|cpus| physical_core_count(cpus, Path::new("/sys/devices/system/cpu")))
        .unwrap_or(logical)
        .max(1);
    let cgroup = cgroup_cpu_quota();
    let scheduler = scheduler_cpu_limit();
    let mut workers = physical;
    if let Some(limit) = cgroup {
        workers = workers.min(limit);
    }
    if let Some((_, limit)) = scheduler.as_ref() {
        workers = workers.min(*limit);
    }
    let mut source =
        format!("auto: {physical} physical core(s) from {logical} allowed logical CPU(s)");
    if let Some(limit) = cgroup {
        source.push_str(&format!(", cgroup quota {limit}"));
    }
    if let Some((name, limit)) = scheduler {
        source.push_str(&format!(", {name}={limit}"));
    }
    (workers.max(1), source)
}

fn available_memory_mib() -> Option<u64> {
    let contents = fs::read_to_string("/proc/meminfo").ok()?;
    contents.lines().find_map(|line| {
        let mut fields = line.split_whitespace();
        (fields.next()? == "MemAvailable:")
            .then(|| fields.next()?.parse::<u64>().ok())??
            .checked_div(1024)
    })
}

fn read_memory_bytes(path: &Path) -> Option<u64> {
    let value = fs::read_to_string(path).ok()?;
    let bytes = value.trim().parse::<u64>().ok()?;
    // cgroup v1 uses a very large sentinel for an unlimited controller.
    (bytes < (1_u64 << 60)).then_some(bytes)
}

fn cgroup_relative_path(controller: Option<&str>) -> Option<PathBuf> {
    let contents = fs::read_to_string("/proc/self/cgroup").ok()?;
    contents.lines().find_map(|line| {
        let mut fields = line.splitn(3, ':');
        let _hierarchy = fields.next()?;
        let controllers = fields.next()?;
        let path = fields.next()?;
        let matches = match controller {
            None => controllers.is_empty(),
            Some(controller) => controllers.split(',').any(|name| name == controller),
        };
        matches.then(|| PathBuf::from(path.trim_start_matches('/')))
    })
}

fn cgroup_tree_available_memory_mib(
    root: &Path,
    relative: &Path,
    limit_name: &str,
    current_name: &str,
) -> Option<u64> {
    let mut directory = root.join(relative);
    let mut available = None;
    loop {
        if let (Some(limit), Some(current)) = (
            read_memory_bytes(&directory.join(limit_name)),
            read_memory_bytes(&directory.join(current_name)),
        ) {
            let remaining = limit.saturating_sub(current) / 1024 / 1024;
            available = Some(available.map_or(remaining, |prior: u64| prior.min(remaining)));
        }
        if directory == root || !directory.pop() {
            break;
        }
    }
    available
}

fn cgroup_available_memory_mib() -> Option<u64> {
    // `/proc/self/cgroup` identifies the process's leaf cgroup.  Inspecting
    // only the mount root misses systemd/Kubernetes parent limits; walk from
    // the leaf back to the controller root and retain the tightest remaining
    // allowance.
    let v2 = cgroup_relative_path(None).and_then(|relative| {
        cgroup_tree_available_memory_mib(
            Path::new("/sys/fs/cgroup"),
            &relative,
            "memory.max",
            "memory.current",
        )
    });
    let v1 = cgroup_relative_path(Some("memory")).and_then(|relative| {
        cgroup_tree_available_memory_mib(
            Path::new("/sys/fs/cgroup/memory"),
            &relative,
            "memory.limit_in_bytes",
            "memory.usage_in_bytes",
        )
    });
    match (v2, v1) {
        (Some(v2), Some(v1)) => Some(v2.min(v1)),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

fn effective_available_memory_mib() -> Option<u64> {
    match (available_memory_mib(), cgroup_available_memory_mib()) {
        (Some(host), Some(cgroup)) => Some(host.min(cgroup)),
        (Some(host), None) => Some(host),
        (None, Some(cgroup)) => Some(cgroup),
        (None, None) => None,
    }
}

fn resolve_uce_memory_limit_mib(concurrent_samples: usize) -> u64 {
    // Reserve half of effective availability for the OS, filesystem cache,
    // assemblers, and decode buffers.  Cap an individual bank at 4 GiB to
    // keep a small cohort on a large host predictable.  If memory accounting
    // is unavailable, retain the former conservative 512 MiB behavior.
    resolve_uce_memory_limit_from_available(effective_available_memory_mib(), concurrent_samples)
}

fn resolve_uce_memory_limit_from_available(
    available_memory_mib: Option<u64>,
    concurrent_samples: usize,
) -> u64 {
    available_memory_mib
        .map(|available| (available / 2 / concurrent_samples.max(1) as u64).clamp(1, 4096))
        .unwrap_or(512)
}

/// Match the upstream GeneMiner2 scheduler: `-p` is a shared CPU budget, not
/// simply a count of whole-sample workers.  Filter, refilter, and assembly
/// jobs advance through separate queues as soon as their predecessor finishes.
fn execute_uce_original_schedule(
    opt: Arc<Options>,
    bins: PathBuf,
    samples: &[Sample],
    dictionary: PathBuf,
    profile: Option<WorkflowProfile>,
) -> Vec<String> {
    let has_filter = opt.commands.iter().any(|command| command == "filter");
    let has_refilter = opt.commands.iter().any(|command| command == "refilter");
    let has_assemble = opt.commands.iter().any(|command| command == "assemble");
    let fused_filter = !opt.legacy_uce_filter;
    // Preserve the upstream GeneMiner2 1--2-unit recruitment budget while
    // keeping the UCEFilter component contract separate: UCEFilter currently
    // implements exactly one compute worker.
    let filter_threads = if opt.workers < 4 { 1 } else { 2 };
    let filter_compute_threads = 1;
    let assembler_threads = if opt.workers == 1 {
        1
    } else {
        (opt.workers / 2).clamp(2, 6)
    };
    let mut available = opt.workers;
    let mut filter_queue = Vec::new();
    let mut refilter_queue = Vec::new();
    let mut assemble_queue = Vec::new();
    let mut failures = Vec::new();
    let mut stop_launching = false;

    // Reverse-and-pop preserves sample-list order, as in the original Python
    // dispatcher, while still allowing completed samples to join the next
    // stage immediately.
    if has_filter {
        filter_queue.extend(samples.iter().rev().cloned());
    } else if has_refilter && opt.legacy_uce_filter {
        refilter_queue.extend(samples.iter().rev().cloned());
    } else if has_assemble {
        assemble_queue.extend(samples.iter().rev().cloned());
    }

    let (sender, receiver) = mpsc::channel::<(Sample, &'static str, usize, Result<(), String>)>();
    let mut running = 0usize;
    while !filter_queue.is_empty()
        || !refilter_queue.is_empty()
        || !assemble_queue.is_empty()
        || running > 0
    {
        let minimum_next = if !filter_queue.is_empty() {
            filter_threads.min(assembler_threads)
        } else {
            assembler_threads
        };
        let task_threads = |available: usize| {
            if available.saturating_sub(assembler_threads) < minimum_next {
                available
            } else {
                assembler_threads
            }
        };
        let launch = |sample: Sample, stage: &'static str, threads: usize| {
            let sender = sender.clone();
            let opt = Arc::clone(&opt);
            let bins = bins.clone();
            let dictionary = dictionary.clone();
            let profile = profile.clone();
            thread::spawn(move || {
                let mut stage_opt = (*opt).clone();
                stage_opt.commands = vec![stage.into()];
                let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    if stage_opt.legacy_uce_filter {
                        execute_uce_legacy(
                            &stage_opt,
                            &bins,
                            &sample,
                            &dictionary,
                            profile.as_ref(),
                            threads,
                        )
                    } else {
                        execute_uce(
                            &stage_opt,
                            &bins,
                            &sample,
                            profile.as_ref(),
                            filter_compute_threads,
                            threads,
                        )
                    }
                }))
                .unwrap_or_else(|_| Err("UCE workflow worker panicked".into()));
                let _ = sender.send((sample, stage, threads, result));
            });
        };

        // Upstream gives ready refilter jobs first choice, then assembly, then
        // recruitment.  Fused UCE filtering advances directly to assembly.
        while !stop_launching && !refilter_queue.is_empty() && available >= filter_threads {
            let sample = refilter_queue.pop().expect("queue was checked");
            let threads = task_threads(available);
            available -= threads;
            launch(sample, "refilter", threads);
            running += 1;
        }
        while !stop_launching && !assemble_queue.is_empty() && available >= assembler_threads {
            let sample = assemble_queue.pop().expect("queue was checked");
            let threads = task_threads(available);
            available -= threads;
            launch(sample, "assemble", threads);
            running += 1;
        }
        while !stop_launching && !filter_queue.is_empty() && available >= filter_threads {
            let sample = filter_queue.pop().expect("queue was checked");
            available -= filter_threads;
            launch(sample, "filter", filter_threads);
            running += 1;
        }

        if running == 0 {
            break;
        }
        let (sample, stage, threads, result) = receiver.recv().expect("UCE worker channel closed");
        running -= 1;
        available += threads;
        match result {
            Ok(()) if !stop_launching && stage == "filter" && has_assemble => {
                if fused_filter || !has_refilter {
                    assemble_queue.push(sample);
                } else {
                    refilter_queue.push(sample);
                }
            }
            Ok(()) if !stop_launching && stage == "refilter" && has_assemble => {
                assemble_queue.push(sample)
            }
            Ok(()) => {}
            Err(error) => {
                failures.push(format!("{} {stage}: {error}", sample.name));
                stop_launching = true;
                filter_queue.clear();
                refilter_queue.clear();
                assemble_queue.clear();
            }
        }
    }
    failures
}

fn execute_native(mut opt: Options) -> Result<(), String> {
    let workflow_started = Instant::now();
    if opt.output.is_empty() {
        return Err("-o is required".into());
    }
    if opt.workers == 0 {
        return Err("-p must be at least 1".into());
    }
    eprintln!("CPU budget: {} ({})", opt.workers, opt.worker_source);
    validate_cleanup_options(&opt)?;
    validate_parallelism(&opt)?;
    let bins = components()?;
    let standalone = [
        "te",
        "gene-annotate",
        "gene-resolve",
        "gene-tree",
        "profiling",
        "mito",
        "rad",
        "rad-probe",
        "rad-validate",
    ];
    if opt.commands.len() > 1
        && opt
            .commands
            .iter()
            .any(|command| standalone.contains(&command.as_str()))
    {
        return Err("this Rust migration route currently requires the selected post-processing command to run alone".into());
    }
    if opt.commands == ["gene-annotate"] {
        return execute_gene_annotate(&opt, &bins);
    }
    if opt.commands == ["gene-resolve"] {
        return execute_gene_resolve(&opt, &bins);
    }
    if opt.commands == ["gene-tree"] {
        return execute_gene_tree(&opt);
    }
    if opt.commands == ["mito"] {
        if opt.samples.is_empty() {
            return Err("-f is required for mito".into());
        }
        let samples = read_samples(&opt.samples, Path::new(&opt.output))?;
        return execute_mito(&opt, &bins, &samples);
    }
    if opt.commands == ["rad-probe"] {
        return execute_rad_probe(&opt, &bins);
    }
    if opt.commands == ["rad-validate"] {
        return execute_rad_validate(&opt, &bins);
    }
    if opt.commands == ["rad"] {
        if opt.samples.is_empty() {
            return Err("-f is required for rad".into());
        }
        return execute_rad(&opt, &bins);
    }
    if opt.commands == ["profiling"] {
        if opt.samples.is_empty() {
            return Err("-f is required for profiling".into());
        }
        let samples = read_samples(&opt.samples, Path::new(&opt.output))?;
        return execute_profiling(&opt, &bins, &samples);
    }
    if opt.samples.is_empty() {
        return Err("-f is required for this command".into());
    }
    if opt.commands == ["te"] {
        return execute_te(&opt, &bins);
    }
    if opt.commands == ["population"] {
        return execute_population(&opt, &bins);
    }
    if opt.reference.is_empty() {
        return Err("-r is required for this command".into());
    }
    if opt.resume {
        if !Path::new(&opt.output).is_dir() {
            return Err("--resume requires an existing workflow output directory".into());
        }
        let samples =
            read_samples_with_directory_creation(&opt.samples, Path::new(&opt.output), false)?;
        if resume_completed_workflow(&opt, &samples)? {
            return Ok(());
        }
    }
    fs::create_dir_all(&opt.output).map_err(|e| e.to_string())?;
    let samples = read_samples(&opt.samples, Path::new(&opt.output))?;
    write_workflow_manifest(&opt, &samples)?;
    if opt.commands == ["stats"] {
        return execute_stats(&opt, &bins, &samples);
    }
    if opt.commands == ["consensus"] {
        return execute_consensus(&opt, &bins, &samples);
    }
    if opt.commands == ["trim"] {
        return execute_trim(&opt, &bins, &samples, "assembly");
    }
    if opt.commands == ["combine"] {
        return execute_combine(&opt, &bins, &samples, "assembly");
    }
    if opt.commands == ["tree"] {
        return execute_tree(&opt);
    }
    let cohort_samples = samples
        .iter()
        .map(|sample| sample.name.clone())
        .collect::<Vec<_>>();
    let is_uce = opt.assembly_mode == "uce";
    if is_uce {
        let filter_threads = if opt.workers < 4 { 1 } else { 2 };
        let concurrent_samples = samples.len().min((opt.workers / filter_threads).max(1));
        opt.uce_memory_limit_mib = resolve_uce_memory_limit_mib(concurrent_samples);
        eprintln!(
            "Auto UCEFilter memory limit: {} MiB per sample ({} concurrent filter job(s))",
            opt.uce_memory_limit_mib, concurrent_samples
        );
    }
    if is_uce {
        let implementation = value(&opt.raw, &["--assembler-implementation"], "auto")?;
        if !matches!(implementation.as_str(), "auto" | "uce-rust") {
            return Err("UCE assembly requires --assembler-implementation auto or uce-rust".into());
        }
    }
    if is_uce && samples.iter().any(|sample| sample.read2.is_none()) {
        return Err(
            "UCE workflow requires paired input; two-column sample lists retain the legacy duplicated-mate convention"
                .into(),
        );
    }
    if is_uce
        && !opt.legacy_uce_filter
        && opt.commands.iter().any(|command| command == "refilter")
        && !opt.commands.iter().any(|command| command == "filter")
    {
        return Err("UCE refilter is fused into the filter stage; run filter first".into());
    }
    let profiler = opt.workflow_profile.then(WorkflowProfile::default);
    let dictionary = reference_dictionary_path(&opt)?;
    if let Some(parent) = dictionary.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    if (!is_uce || opt.legacy_uce_filter) && opt.commands.iter().any(|c| c == "filter") {
        let index_args = vec![
            "-r".into(),
            opt.reference.clone(),
            "-o".into(),
            opt.output.clone(),
            "-kf".into(),
            opt.kf.clone(),
            "-s".into(),
            opt.step.clone(),
            "-gr".into(),
            "-lkd".into(),
            dictionary.display().to_string(),
            "-m".into(),
            "2".into(),
        ];
        let result = run_profiled_action(
            profiler.as_ref(),
            "__reference__",
            "mainfilter_index",
            Path::new(&opt.reference),
            &dictionary,
            || run(&bins, "MainFilterNew", &index_args),
        );
        if let Err(error) = result {
            if let Some(profile) = profiler.as_ref() {
                write_native_workflow_profile(
                    Path::new(&opt.output),
                    profile,
                    workflow_started.elapsed().as_millis(),
                )?;
            }
            return Err(error);
        }
    }
    let failures = Arc::new(Mutex::new(Vec::new()));
    let shared = Arc::new(opt);
    macro_rules! profile_try {
        ($result:expr) => {
            if let Err(error) = $result {
                if let Some(profile) = profiler.as_ref() {
                    write_native_workflow_profile(
                        Path::new(&shared.output),
                        profile,
                        workflow_started.elapsed().as_millis(),
                    )?;
                }
                return Err(error);
            }
        };
    }
    if is_uce {
        failures
            .lock()
            .expect("failure list poisoned")
            .extend(execute_uce_original_schedule(
                Arc::clone(&shared),
                bins.clone(),
                &samples,
                dictionary.clone(),
                profiler.clone(),
            ));
    } else {
        let next = Arc::new(Mutex::new(samples.clone().into_iter()));
        let mut handles = Vec::new();
        let cancelled = Arc::new(AtomicBool::new(false));
        for _ in 0..shared.workers {
            let bins = bins.clone();
            let opt = Arc::clone(&shared);
            let next = Arc::clone(&next);
            let failures = Arc::clone(&failures);
            let cancelled = Arc::clone(&cancelled);
            let dictionary = dictionary.clone();
            let profiler = profiler.clone();
            handles.push(thread::spawn(move || loop {
                if cancelled.load(Ordering::Acquire) {
                    break;
                }
                let Some(sample) = next.lock().expect("sample queue poisoned").next() else {
                    break;
                };
                let result = execute_gene(&opt, &bins, &sample, &dictionary, profiler.as_ref());
                if let Err(error) = result {
                    failures
                        .lock()
                        .expect("failure list poisoned")
                        .push(format!("{}: {error}", sample.name));
                    cancelled.store(true, Ordering::Release);
                }
            }));
        }
        for handle in handles {
            handle.join().map_err(|_| "Rust workflow worker panicked")?;
        }
    }
    let failures = failures.lock().map_err(|_| "failure list poisoned")?;
    if !failures.is_empty() {
        if let Some(profile) = profiler.as_ref() {
            write_native_workflow_profile(
                Path::new(&shared.output),
                profile,
                workflow_started.elapsed().as_millis(),
            )?;
        }
        return Err(format!(
            "{} sample(s) failed:\n{}",
            failures.len(),
            failures.join("\n")
        ));
    }
    if !is_uce && shared.commands.iter().any(|c| c == "gene") {
        let mut cohort = vec![
            "cohort".into(),
            "--reference".into(),
            shared.reference.clone(),
            "--out".into(),
            Path::new(&shared.output).join("gene").display().to_string(),
        ];
        for name in cohort_samples {
            cohort.extend(["--sample".into(), name]);
        }
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__cohort__",
            "gene-cohort",
            Path::new(&shared.output),
            &Path::new(&shared.output).join("gene"),
            || run(&bins, "gene_workflow", &cohort),
        ));
    }
    if shared.commands.iter().any(|command| command == "consensus") {
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__workflow__",
            "consensus",
            Path::new(&shared.output),
            Path::new(&shared.output),
            || execute_consensus(&shared, &bins, &samples),
        ));
    }
    if shared.commands.iter().any(|command| command == "trim") {
        let source = if shared.commands.iter().any(|command| command == "consensus") {
            "consensus"
        } else {
            "assembly"
        };
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__workflow__",
            "trim",
            Path::new(&shared.output),
            Path::new(&shared.output),
            || execute_trim(&shared, &bins, &samples, source),
        ));
    }
    if shared.commands.iter().any(|command| command == "combine") {
        let source = if shared.commands.iter().any(|command| command == "trim") {
            "trimmed"
        } else if shared.commands.iter().any(|command| command == "consensus") {
            "consensus"
        } else {
            "assembly"
        };
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__workflow__",
            "combine",
            Path::new(&shared.output),
            Path::new(&shared.output),
            || execute_combine(&shared, &bins, &samples, source),
        ));
    }
    if shared.commands.iter().any(|command| command == "tree") {
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__workflow__",
            "tree",
            Path::new(&shared.output),
            Path::new(&shared.output),
            || execute_tree(&shared),
        ));
    }
    if shared.commands.iter().any(|command| command == "stats") {
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__workflow__",
            "stats",
            Path::new(&shared.output),
            Path::new(&shared.output),
            || execute_stats(&shared, &bins, &samples),
        ));
    }
    if shared
        .commands
        .iter()
        .any(|command| command == "population")
    {
        profile_try!(run_profiled_action(
            profiler.as_ref(),
            "__workflow__",
            "population",
            Path::new(&shared.output),
            Path::new(&shared.output),
            || execute_population(&shared, &bins),
        ));
    }
    profile_try!(run_profiled_action(
        profiler.as_ref(),
        "__workflow__",
        "cleanup",
        Path::new(&shared.output),
        Path::new(&shared.output),
        || cleanup_native_intermediates(&shared, &samples),
    ));
    if let Some(profile) = profiler.as_ref() {
        write_native_workflow_profile(
            Path::new(&shared.output),
            profile,
            workflow_started.elapsed().as_millis(),
        )?;
    }
    Ok(())
}

fn print_help() {
    let (workers, source) = auto_worker_budget();
    println!(
        "TipSeek CLI\n\nNative Rust command dispatcher; no Python runtime is required.\n\n\
Usage: tipseek [COMMAND ...] -f SAMPLES -r REFERENCES -o OUTPUT [-p INT|auto]\n\n\
Parallelism:\n  \
-p INT|auto  Shared CPU budget. The default is auto, which counts physical\n               \
cores allowed by affinity/cpuset and caps them by cgroup or scheduler limits.\n               \
Use an integer to override automatic detection.\n\n\
UCE recruitment:\n  \
--uce-recruit-mode fast|auto\n               \
Use the default automatic two-pass recruitment, or select fast to keep only\n               \
the initial pass. UCE defaults to k=23, step=4, and auto; original assembly\n               \
keeps k=31, step=4, and fast unless the user explicitly overrides them.\n  \
--uce-fallback-kmer-size INT  Sensitive-pass recruitment k (default: 21).\n  \
--uce-fallback-step INT       Sensitive-pass read-scan step (default: 1).\n  \
--uce-fallback-verify-kmer-size INT\n               \
Independent exact-match verification k (default: 19). Auto mode checks\n               \
ambiguity against the complete probe panel and retains unique-locus reads.\n  \
--uce-fallback-min-alignment-overlap INT\n               \
Minimum local alignment overlap for a fallback read pair (default: 45 bp).\n  \
--uce-fallback-min-alignment-identity FLOAT\n               \
Minimum local alignment identity for a fallback read pair (default: 0.80).\n               \
Fallback-only assembled contigs must be at least 200 bp, align to the target\n               \
probe at >=80% coverage and >=80% identity, have no near-tied panel locus,\n               \
and pass the provisional-core inverted-repeat guard before rescue. Internal\n               \
read-chain gaps >=40 bp are reported for review rather than rejected alone.\n\n\
UCE rescue:\n  \
--uce-rescue-reads  Explicitly enable rescue after the initial UCE assembly\n               \
               (already enabled by default in UCE mode). In auto mode, only\n               \
               anchored provisional cores seed rescue.\n  \
--no-uce-rescue-reads  Disable the default UCE rescue stage.\n  \
--uce-rescue-rounds 1|2  Number of rescue rounds (default: 1).\n  \
--uce-rescue-reverse-reuse-reference-scale FLOAT\n               \
Scale only the reference bonus when a reverse-complement node is already\n               \
present in either rescue assembly arm (default: 1.0; range: 0-1; 1 disables).\n  \
--uce-rescue-inverted-repeat-min-bp INT\n               \
Reject a provisional fallback core or roll back a rescue round when it contains\n               \
or newly introduces an exact long inverted repeat (default: 150 bp; 0 disables).\n\n\
Detected auto budget: {workers} ({source})"
    );
}

fn main() -> ExitCode {
    let args = env::args().skip(1).collect::<Vec<_>>();
    if args.iter().any(|arg| arg == "-h" || arg == "--help") {
        print_help();
        return ExitCode::SUCCESS;
    }
    match parse(&args).and_then(execute_with_status) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("Error: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpu_lists_expand_ranges_and_reject_invalid_input() {
        assert_eq!(
            parse_cpu_list("0-3,8,10-11"),
            Some(vec![0, 1, 2, 3, 8, 10, 11])
        );
        assert_eq!(parse_cpu_list("3-1"), None);
        assert_eq!(parse_cpu_list(""), None);
    }

    #[test]
    fn adaptive_mito_only_settles_identical_circular_observations() {
        let partial = ("linear_single_contig".into(), "same-sequence".into());
        let empty = ("no_contigs".into(), String::new());
        let circular = ("circular".into(), "canonical-circle".into());
        assert!(!is_stable_circular(Some(&partial), &partial));
        assert!(!is_stable_circular(Some(&empty), &empty));
        assert!(!is_stable_circular(None, &circular));
        assert!(is_stable_circular(Some(&circular), &circular));
    }

    #[test]
    fn physical_cores_are_deduplicated_across_smt_siblings() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gm2_cpu_topology_{}_{}",
            std::process::id(),
            unique
        ));
        for (cpu, package, core) in [(0, 0, 0), (1, 0, 0), (2, 0, 1), (3, 1, 0)] {
            let topology = root.join(format!("cpu{cpu}/topology"));
            fs::create_dir_all(&topology).unwrap();
            fs::write(topology.join("physical_package_id"), package.to_string()).unwrap();
            fs::write(topology.join("core_id"), core.to_string()).unwrap();
        }
        assert_eq!(physical_core_count(&[0, 1, 2, 3], &root), Some(3));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn cpu_quota_parsers_are_conservative_and_handle_unlimited_values() {
        assert_eq!(parse_cpu_max("max 100000"), None);
        assert_eq!(parse_cpu_max("250000 100000"), Some(2));
        assert_eq!(parse_cpu_max("50000 100000"), Some(1));
        assert_eq!(parse_cpu_cfs("-1", "100000"), None);
        assert_eq!(parse_cpu_cfs("1200000", "100000"), Some(12));
    }

    #[test]
    fn worker_budget_defaults_to_auto_and_explicit_values_override_it() {
        let auto = parse(&[
            "stats".into(),
            "-f".into(),
            "samples.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap();
        assert!(auto.workers >= 1);
        assert!(auto.worker_source.starts_with("auto:"));
        let explicit = resolve_worker_budget("12").unwrap();
        assert_eq!(explicit.0, 12);
        assert_eq!(explicit.1, "explicit -p 12");
        assert!(resolve_worker_budget("0").is_err());
        assert!(resolve_worker_budget("physical").is_err());
    }

    #[test]
    fn uce_defaults_to_k23_auto_one_round_rescue_and_builds_a_conservative_fallback() {
        let base = [
            "filter",
            "--assembly-mode",
            "uce",
            "-f",
            "samples.tsv",
            "-r",
            "references",
            "-o",
            "out",
        ]
        .into_iter()
        .map(str::to_string)
        .collect::<Vec<_>>();
        let defaults = parse(&base).unwrap();
        assert_eq!(defaults.kf, "23");
        assert_eq!(defaults.step, "4");
        assert_eq!(defaults.uce_recruit_mode, "auto");
        assert!(defaults.rescue);
        assert_eq!(DEFAULT_UCE_RESCUE_ROUNDS, "1");

        let mut fast_options = base.clone();
        fast_options.extend(["--uce-recruit-mode".into(), "fast".into()]);
        let fast = parse(&fast_options).unwrap();
        assert_eq!(fast.uce_recruit_mode, "fast");
        let sample = Sample {
            name: "sample".into(),
            read1: "r1.fq".into(),
            read2: Some("r2.fq".into()),
        };
        let fast_args = uce_filter_args(&fast, &sample, Path::new("out/sample"), 1);
        assert!(!fast_args
            .iter()
            .any(|argument| argument == "--verification-kmer-size"));
        assert!(!fast_args
            .iter()
            .any(|argument| argument == "--max-locus-count"));
        assert!(!fast_args
            .iter()
            .any(|argument| argument == "--retain-loci-file"));

        let mut auto_args = base;
        auto_args.extend(
            [
                "--uce-recruit-mode",
                "auto",
                "--uce-fallback-kmer-size",
                "23",
                "--uce-fallback-step",
                "2",
                "--uce-fallback-verify-kmer-size",
                "17",
            ]
            .into_iter()
            .map(str::to_string),
        );
        let auto = parse(&auto_args).unwrap();
        let pass = RecruitPass::fallback(
            &auto.uce_fallback_kmer_size,
            &auto.uce_fallback_step,
            &auto.uce_fallback_verify_kmer_size,
            &auto.uce_fallback_min_alignment_overlap,
            &auto.uce_fallback_min_alignment_identity,
        );
        let fallback_args = uce_filter_args_for_pass(
            &auto,
            &sample,
            &pass,
            &UceRecruitInvocation {
                sample_dir: Path::new("fallback"),
                verify_reference: Path::new("references"),
                recruit_reference: Path::new("fallback_references"),
                role: "bait",
                retain_loci_file: Some(Path::new("unresolved.txt")),
            },
        );
        let argument_value = |name: &str| {
            let index = fallback_args
                .iter()
                .position(|argument| argument == name)
                .expect("fallback argument is present");
            fallback_args[index + 1].as_str()
        };
        assert_eq!(argument_value("-kf"), "23");
        assert_eq!(argument_value("-s"), "2");
        assert_eq!(
            argument_value("--recruit-references"),
            "fallback_references"
        );
        assert_eq!(argument_value("--verification-kmer-size"), "17");
        assert_eq!(argument_value("--minimum-alignment-overlap"), "45");
        assert_eq!(argument_value("--minimum-alignment-identity"), "0.80");
        assert_eq!(argument_value("--max-locus-count"), "1");
        assert_eq!(argument_value("--retain-loci-file"), "unresolved.txt");
    }

    #[test]
    fn original_defaults_are_unchanged_and_uce_defaults_can_be_overridden() {
        let original = parse(&[
            "filter".into(),
            "-f".into(),
            "samples.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap();
        assert_eq!(original.kf, "31");
        assert_eq!(original.step, "4");
        assert_eq!(original.uce_recruit_mode, "fast");
        assert!(!original.rescue);

        let uce_override = parse(&[
            "filter".into(),
            "--assembly-mode".into(),
            "uce".into(),
            "-kf".into(),
            "31".into(),
            "--uce-recruit-mode".into(),
            "fast".into(),
            "--no-uce-rescue-reads".into(),
            "-f".into(),
            "samples.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap();
        assert_eq!(uce_override.kf, "31");
        assert_eq!(uce_override.uce_recruit_mode, "fast");
        assert!(!uce_override.rescue);

        let conflict = parse(&[
            "filter".into(),
            "--assembly-mode".into(),
            "uce".into(),
            "--uce-rescue-reads".into(),
            "--no-uce-rescue-reads".into(),
            "-f".into(),
            "samples.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap_err();
        assert!(conflict.contains("cannot be used together"));
    }

    #[test]
    fn invalid_uce_recruit_options_are_rejected() {
        let parse_uce = |extra: &[&str]| {
            let mut args = vec![
                "filter".into(),
                "--assembly-mode".into(),
                "uce".into(),
                "-f".into(),
                "samples.tsv".into(),
                "-r".into(),
                "references".into(),
                "-o".into(),
                "out".into(),
            ];
            args.extend(extra.iter().map(|value| (*value).to_owned()));
            parse(&args)
        };
        assert!(parse_uce(&["--uce-recruit-mode", "typo"]).is_err());
        assert!(parse_uce(&[
            "--uce-recruit-mode",
            "auto",
            "--uce-fallback-kmer-size",
            "0"
        ])
        .is_err());
        assert!(parse_uce(&[
            "--uce-recruit-mode",
            "auto",
            "--uce-fallback-verify-kmer-size",
            "65"
        ])
        .is_err());
        assert!(parse_uce(&["--uce-recruit-mode", "auto", "--uce-fallback-step", "0"]).is_err());
        assert!(parse_uce(&[
            "--uce-recruit-mode",
            "auto",
            "--uce-fallback-min-alignment-identity",
            "1.1"
        ])
        .is_err());
        assert!(parse_uce(&["--uce-recruit-mode", "auto", "--legacy-uce-filter"]).is_err());
    }

    #[test]
    fn uce_rescue_uses_fixed_k21_without_changing_normal_assembly_kmer() {
        let opt = parse(&[
            "assemble".into(),
            "--assembly-mode".into(),
            "uce".into(),
            "-p".into(),
            "1".into(),
            "-f".into(),
            "samples.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
            "-ka".into(),
            "51".into(),
        ])
        .unwrap();
        let normal = uce_assembler_args(&opt, Path::new("sample"), 1).unwrap();
        let rescue =
            uce_rescue_assembler_args(&opt, Path::new("sample"), Path::new("rescue"), 1).unwrap();
        let argument_value = |args: &[String], name: &str| {
            let index = args
                .iter()
                .position(|argument| argument == name)
                .expect("argument is present");
            args[index + 1].clone()
        };
        assert_eq!(argument_value(&normal, "-ka"), "51");
        assert_eq!(argument_value(&rescue, "-ka"), UCE_RESCUE_ASSEMBLY_KMER);
        assert_eq!(argument_value(&normal, "-r"), "references");
        assert_eq!(argument_value(&rescue, "-r"), "rescue");
        assert!(!normal
            .iter()
            .any(|argument| argument == "--uce-reverse-reuse-reference-scale"));
        assert_eq!(
            argument_value(&rescue, "--uce-reverse-reuse-reference-scale"),
            "1.0"
        );
    }

    fn test_dna(length: usize, seed: u64) -> String {
        let mut state = seed;
        (0..length)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(state >> 62) as usize] as char
            })
            .collect()
    }

    #[test]
    fn rescue_guard_only_flags_a_new_long_inverted_repeat() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gm2_inverted_repeat_guard_{}_{unique}",
            std::process::id()
        ));
        let sample = root.join("sample");
        let backup = root.join("backup");
        fs::create_dir_all(sample.join("results")).unwrap();
        fs::create_dir_all(backup.join("results")).unwrap();
        let arm = test_dna(180, 44);
        let before = test_dna(500, 45);
        let after = format!(
            "{arm}{before}{}{}",
            test_dna(60, 46),
            reverse_complement_text(&arm)
        );
        fs::write(
            backup.join("results/locus.fasta"),
            format!(">before\n{before}\n"),
        )
        .unwrap();
        fs::write(
            sample.join("results/locus.fasta"),
            format!(">after\n{after}\n"),
        )
        .unwrap();
        assert!(rescue_introduces_long_inverted_repeat(&sample, &backup, "locus", 150).unwrap());

        fs::write(
            backup.join("results/locus.fasta"),
            format!(">before\n{after}\n"),
        )
        .unwrap();
        assert!(!rescue_introduces_long_inverted_repeat(&sample, &backup, "locus", 150).unwrap());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn terminal_rescue_reference_contains_only_active_loci() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root =
            std::env::temp_dir().join(format!("gm2_terminal_refs_{}_{unique}", std::process::id()));
        let reference = root.join("reference");
        let sample = root.join("sample");
        let rescue = root.join("rescue");
        fs::create_dir_all(&reference).unwrap();
        fs::create_dir_all(sample.join("results")).unwrap();
        fs::write(reference.join("active.fasta"), ">ref\nAAAA\n").unwrap();
        fs::write(reference.join("inactive.fasta"), ">ref\nCCCC\n").unwrap();
        fs::write(sample.join("results/active.fasta"), ">contig\nAAAAGGGG\n").unwrap();
        fs::write(
            sample.join("uce_assembly_summary.csv"),
            "locus,status,accepted,low_quality\nactive,success,1,0\ninactive,success,1,0\n",
        )
        .unwrap();
        let active = ["active".to_owned()].into_iter().collect();
        assert_eq!(
            build_uce_rescue_reference(&reference, &sample, &rescue, 4, Some(&active)).unwrap(),
            1
        );
        assert!(rescue.join("active.fasta").is_file());
        assert!(!rescue.join("inactive.fasta").exists());
        assert!(fs::read_to_string(rescue.join("active.fasta"))
            .unwrap()
            .contains("AAAAGGGG"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn review_only_provisional_core_is_not_rescue_eligible() {
        let clean = BTreeMap::from([
            ("locus".into(), "clean".into()),
            ("auto_recruit_core_anchor_status".into(), "anchored".into()),
        ]);
        let review = BTreeMap::from([
            ("locus".into(), "review".into()),
            (
                "auto_recruit_core_anchor_status".into(),
                "anchored_with_review".into(),
            ),
        ]);
        let summary = UceSummary {
            headers: vec!["locus".into(), "auto_recruit_core_anchor_status".into()],
            rows: BTreeMap::from([("clean".into(), clean), ("review".into(), review)]),
        };
        assert_eq!(
            review_only_provisional_cores(&summary),
            ["review".to_owned()].into_iter().collect()
        );
    }

    #[test]
    fn terminal_reconcile_keeps_supported_side_and_preserves_candidates() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gm2_terminal_reconcile_{}_{unique}",
            std::process::id()
        ));
        let sample = root.join("sample");
        let backup = root.join("backup");
        for directory in [
            sample.join("results"),
            sample.join("filtered"),
            sample.join("contigs_all"),
            sample.join("contigs_all_low"),
            backup.join("results"),
        ] {
            fs::create_dir_all(directory).unwrap();
        }
        let left = test_dna(40, 1);
        let core = test_dna(120, 2);
        let right = test_dna(40, 3);
        let assembled = format!("{left}{core}{right}");
        fs::write(
            backup.join("results/locus.fasta"),
            format!(">old\n{core}\n"),
        )
        .unwrap();
        fs::write(
            sample.join("results/locus.fasta"),
            format!(">new\n{assembled}\n"),
        )
        .unwrap();
        let candidates = format!(">candidate_1\n{assembled}\n>candidate_2\nACGT\n");
        fs::write(sample.join("contigs_all/locus.fasta"), &candidates).unwrap();
        fs::write(sample.join("contigs_all_low/locus.fasta"), &candidates).unwrap();
        let spanning = &assembled[..100];
        let fastq = ["frag1", "frag2"]
            .into_iter()
            .map(|fragment| {
                format!(
                    "@{fragment}/1\n{spanning}\n+\n{}\n",
                    "I".repeat(spanning.len())
                )
            })
            .collect::<String>();
        fs::write(sample.join("filtered/locus.fq"), fastq).unwrap();
        let mut after = std::collections::BTreeMap::from([
            ("selected_contig_length".into(), assembled.len().to_string()),
            ("unique_read_count".into(), "2".into()),
        ]);
        let (evidence, status) =
            terminal_reconcile_locus(&sample, &backup, "locus", &mut after).unwrap();
        let evidence = evidence.unwrap();
        assert_eq!(status, "accepted");
        assert!(evidence.left.accepted);
        assert!(!evidence.right.accepted);
        assert_eq!(
            first_fasta_sequence(&sample.join("results/locus.fasta")).unwrap(),
            Some(format!("{left}{core}"))
        );
        assert_eq!(
            fs::read_to_string(sample.join("contigs_all/locus.fasta")).unwrap(),
            candidates
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn prior_rescue_round_evidence_survives_a_later_round() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gm2_rescue_evidence_{}_{unique}",
            std::process::id()
        ));
        let sample = root.join("sample");
        let backup = root.join("backup");
        fs::create_dir_all(&sample).unwrap();
        fs::create_dir_all(backup.join("uce_rescue_round_1")).unwrap();
        fs::write(backup.join("uce_rescue_round_1/evidence.txt"), "round one").unwrap();
        restore_prior_rescue_rounds(&sample, &backup, 2).unwrap();
        assert_eq!(
            fs::read_to_string(sample.join("uce_rescue_round_1/evidence.txt")).unwrap(),
            "round one"
        );
        assert!(backup.join("uce_rescue_round_1/evidence.txt").is_file());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn initial_uce_recruit_audits_survive_rescue_rebuild() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gm2_recruit_audits_{}_{unique}",
            std::process::id()
        ));
        let sample = root.join("sample");
        let backup = root.join("backup");
        fs::create_dir_all(&sample).unwrap();
        fs::create_dir_all(backup.join("fallback_probe_rejected/results")).unwrap();
        for name in [
            "uce_filter_summary.fast.tsv",
            "uce_filter_summary.fallback.tsv",
            "uce_recruit_passes.tsv",
            "uce_recruit_contig_probe_gate.tsv",
        ] {
            fs::write(backup.join(name), name).unwrap();
        }
        fs::write(
            backup.join("fallback_probe_rejected/results/uce-1.fasta"),
            ">uce-1\nACGT\n",
        )
        .unwrap();

        restore_initial_uce_recruit_audits(&sample, &backup).unwrap();

        for name in [
            "uce_filter_summary.fast.tsv",
            "uce_filter_summary.fallback.tsv",
            "uce_recruit_passes.tsv",
            "uce_recruit_contig_probe_gate.tsv",
        ] {
            assert_eq!(fs::read_to_string(sample.join(name)).unwrap(), name);
            assert!(backup.join(name).is_file());
        }
        assert_eq!(
            fs::read_to_string(sample.join("fallback_probe_rejected/results/uce-1.fasta")).unwrap(),
            ">uce-1\nACGT\n"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn rescue_summary_uses_per_locus_status() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let directory =
            std::env::temp_dir().join(format!("gm2_rescue_report_{}_{unique}", std::process::id()));
        fs::create_dir_all(&directory).unwrap();
        let row = std::collections::BTreeMap::from([
            ("locus".into(), "locus".into()),
            ("status".into(), "success".into()),
            ("accepted".into(), "1".into()),
            ("selected_contig_length".into(), "100".into()),
            ("unique_read_count".into(), "2".into()),
        ]);
        let summary = UceSummary {
            headers: row.keys().cloned().collect(),
            rows: std::collections::BTreeMap::from([("locus".into(), row)]),
        };
        let sample = Sample {
            name: "sample".into(),
            read1: String::new(),
            read2: None,
        };
        let report = RescueReportContext {
            status_by_locus: std::collections::BTreeMap::from([(
                "locus".into(),
                "reverted_density_drop".into(),
            )]),
            overall_status: "success".into(),
            ..RescueReportContext::default()
        };
        write_rescue_reports(&sample, &directory, &summary, &summary, &[], &report).unwrap();
        assert!(fs::read_to_string(directory.join("uce_rescue_summary.csv"))
            .unwrap()
            .contains("sample,locus,reverted_density_drop"));
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn move_tree_relocates_nested_content_and_clears_the_source() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root =
            std::env::temp_dir().join(format!("gm2_move_tree_{}_{unique}", std::process::id()));
        let source = root.join("source");
        let nested = source.join("nested");
        fs::create_dir_all(&nested).unwrap();
        fs::write(nested.join("leaf.txt"), b"payload").unwrap();
        let destination = root.join("destination");
        move_tree(&source, &destination).unwrap();
        assert!(!source.exists());
        assert_eq!(
            fs::read_to_string(destination.join("nested/leaf.txt")).unwrap(),
            "payload"
        );
        fs::remove_dir_all(&root).unwrap();
    }

    #[test]
    fn uce_memory_limit_reserves_capacity_and_applies_bounds() {
        assert_eq!(resolve_uce_memory_limit_from_available(None, 4), 512);
        assert_eq!(resolve_uce_memory_limit_from_available(Some(8192), 2), 2048);
        assert_eq!(resolve_uce_memory_limit_from_available(Some(1), 1), 1);
        assert_eq!(
            resolve_uce_memory_limit_from_available(Some(65536), 1),
            4096
        );
        assert_eq!(resolve_uce_memory_limit_from_available(Some(8192), 0), 4096);
    }

    #[test]
    fn sample_table_rejects_invalid_rows_before_creating_output_directories() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gm2_sample_table_validation_{}_{}",
            std::process::id(),
            unique
        ));
        fs::create_dir_all(&root).unwrap();
        let read = root.join("reads.fq");
        fs::write(&read, b"@r\nAC\n+\n!!\n").unwrap();
        for (name, contents, expected) in [
            ("one_column", "sample\n", "must be sample<TAB>R1"),
            (
                "too_many_columns",
                "sample\ta\tb\tc\n",
                "must be sample<TAB>R1",
            ),
            ("empty_r1", "sample\t\n", "empty R1 path"),
            ("empty_r2", "sample\ta\t\n", "empty R2 path"),
            (
                "missing_file",
                "sample\tmissing.fq\n",
                "read file does not exist",
            ),
            (
                "duplicate_name",
                &format!("a-b\t{}\na b\t{}\n", read.display(), read.display()),
                "Duplicate sample name after normalization",
            ),
        ] {
            let table = root.join(format!("{name}.tsv"));
            fs::write(&table, contents).unwrap();
            let output = root.join(format!("out_{name}"));
            let error = read_samples(&table.display().to_string(), &output).unwrap_err();
            assert!(error.contains(expected), "{error}");
            assert!(!output.exists(), "invalid table created {output:?}");
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn uce_two_sample_fixture_respects_budget_and_stops_after_a_failure() {
        use std::os::unix::fs::PermissionsExt;
        use std::sync::OnceLock;

        // Environment variables select workflow components, so serialize this
        // end-to-end test with any future component-directory tests.
        static COMPONENT_ENV_LOCK: OnceLock<Mutex<()>> = OnceLock::new();
        let _guard = COMPONENT_ENV_LOCK
            .get_or_init(|| Mutex::new(()))
            .lock()
            .unwrap();
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gm2_uce_schedule_fixture_{}_{}",
            std::process::id(),
            unique
        ));
        let components = root.join("components");
        fs::create_dir_all(&components).unwrap();
        let capture = root.join("component_calls.tsv");
        for component in ["uce_filter", "main_assembler-rust"] {
            let path = components.join(component);
            fs::write(
                &path,
                "#!/bin/sh\nprintf '%s\\t%s\\n' \"$(basename \"$0\")\" \"$*\" >> \"$GM2_CAPTURE\"\nif [ -n \"$GM2_FAIL_MATCH\" ]; then case \"$*\" in *\"$GM2_FAIL_MATCH\"*) exit 7;; esac; fi\nsleep 0.01\n",
            )
            .unwrap();
            let mut permissions = fs::metadata(&path).unwrap().permissions();
            permissions.set_mode(0o755);
            fs::set_permissions(&path, permissions).unwrap();
        }
        let reference = root.join("reference.fasta");
        let read1 = root.join("reads_1.fastq");
        let read2 = root.join("reads_2.fastq");
        fs::write(
            &reference,
            include_str!("../tests/fixtures/uce/reference.fasta"),
        )
        .unwrap();
        fs::write(&read1, include_str!("../tests/fixtures/uce/reads.fastq")).unwrap();
        fs::write(&read2, include_str!("../tests/fixtures/uce/reads.fastq")).unwrap();
        let samples = root.join("samples.tsv");
        fs::write(
            &samples,
            format!(
                "one\t{}\t{}\ntwo\t{}\t{}\n",
                read1.display(),
                read2.display(),
                read1.display(),
                read2.display()
            ),
        )
        .unwrap();
        let prior_components = env::var_os("GM2_COMPONENT_DIR");
        let prior_capture = env::var_os("GM2_CAPTURE");
        let prior_failure = env::var_os("GM2_FAIL_MATCH");
        env::set_var("GM2_COMPONENT_DIR", &components);
        env::set_var("GM2_CAPTURE", &capture);
        env::remove_var("GM2_FAIL_MATCH");
        for (workers, expected_filter_threads, mut expected_assembler_threads) in [
            (1, 1, vec![1, 1]),
            (2, 1, vec![2, 2]),
            (4, 1, vec![2, 2]),
            // Once both filters release their 2-unit reservations, the first
            // ready assembler may consume the six available budget units;
            // the remaining job receives the normal 4-thread share.
            (8, 1, vec![4, 6]),
        ] {
            fs::write(&capture, "").unwrap();
            let output = root.join(format!("output_p{workers}"));
            let result = parse(&[
                "filter".into(),
                "assemble".into(),
                "--assembly-mode".into(),
                "uce".into(),
                "--uce-recruit-mode".into(),
                "fast".into(),
                "--no-uce-rescue-reads".into(),
                "-p".into(),
                workers.to_string(),
                "-f".into(),
                samples.display().to_string(),
                "-r".into(),
                reference.display().to_string(),
                "-o".into(),
                output.display().to_string(),
            ])
            .and_then(execute_with_status);
            result.unwrap();
            let manifest = fs::read_to_string(output.join("workflow_manifest.tsv")).unwrap();
            assert!(manifest.contains("schema_version\t1\n"));
            assert!(manifest.contains("assembly_mode\tuce\n"));
            assert!(manifest.contains("worker_source\texplicit -p"));
            assert!(manifest.contains("sample_count\t2\n"));
            let status = fs::read_to_string(output.join("workflow_status.tsv")).unwrap();
            assert!(status.contains("state\tsucceeded\n"));
            let calls = fs::read_to_string(&capture).unwrap();
            let filters = calls
                .lines()
                .filter(|line| line.starts_with("uce_filter\t"))
                .collect::<Vec<_>>();
            let assemblers = calls
                .lines()
                .filter(|line| line.starts_with("main_assembler-rust\t"))
                .collect::<Vec<_>>();
            assert_eq!(filters.len(), 2, "-p {workers}: {calls}");
            assert!(filters
                .iter()
                .all(|line| { line.contains(&format!("--threads {expected_filter_threads}")) }));
            assert_eq!(assemblers.len(), 2, "-p {workers}: {calls}");
            let mut actual_assembler_threads = assemblers
                .iter()
                .map(|line| {
                    line.split(" -p ")
                        .nth(1)
                        .and_then(|tail| tail.split_whitespace().next())
                        .and_then(|value| value.parse::<usize>().ok())
                        .expect("assembler invocation includes a numeric -p value")
                })
                .collect::<Vec<_>>();
            actual_assembler_threads.sort_unstable();
            expected_assembler_threads.sort_unstable();
            assert_eq!(
                actual_assembler_threads, expected_assembler_threads,
                "-p {workers}: {calls}"
            );
        }
        fs::write(&capture, "").unwrap();
        let failed_output = root.join("output_failure");
        env::set_var("GM2_FAIL_MATCH", failed_output.join("1_One"));
        let failure = parse(&[
            "filter".into(),
            "--assembly-mode".into(),
            "uce".into(),
            "--uce-recruit-mode".into(),
            "fast".into(),
            "--no-uce-rescue-reads".into(),
            "-p".into(),
            "1".into(),
            "-f".into(),
            samples.display().to_string(),
            "-r".into(),
            reference.display().to_string(),
            "-o".into(),
            failed_output.display().to_string(),
        ])
        .and_then(execute_with_status)
        .unwrap_err();
        assert!(failure.contains("1_One"), "{failure}");
        let calls = fs::read_to_string(&capture).unwrap();
        assert_eq!(calls.lines().count(), 1, "{calls}");
        assert!(calls.contains("output_failure/1_One"), "{calls}");
        let status = fs::read_to_string(failed_output.join("workflow_status.tsv")).unwrap();
        assert!(status.contains("state\tfailed\n"));
        match prior_components {
            Some(value) => env::set_var("GM2_COMPONENT_DIR", value),
            None => env::remove_var("GM2_COMPONENT_DIR"),
        }
        match prior_capture {
            Some(value) => env::set_var("GM2_CAPTURE", value),
            None => env::remove_var("GM2_CAPTURE"),
        }
        match prior_failure {
            Some(value) => env::set_var("GM2_FAIL_MATCH", value),
            None => env::remove_var("GM2_FAIL_MATCH"),
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn consensus_uses_legacy_filtered_fastx_extension() {
        assert_eq!(fastx_output_extension("reads.fastq.gz"), ".fq");
        assert_eq!(fastx_output_extension("reads.fa"), ".fasta");
        let consensus = parse(&[
            "consensus".into(),
            "-f".into(),
            "samples.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
            "-c".into(),
            "0.9".into(),
        ])
        .unwrap();
        assert_eq!(consensus.commands, ["consensus"]);
    }

    #[test]
    fn te_and_population_are_native_standalone_commands() {
        let te = parse(&[
            "te".into(),
            "-f".into(),
            "taxa.tsv".into(),
            "-o".into(),
            "out".into(),
            "--te-stage".into(),
            "discover".into(),
        ])
        .unwrap();
        assert_eq!(te.commands, ["te"]);
        let population = parse(&[
            "population".into(),
            "-f".into(),
            "samples.tsv".into(),
            "-o".into(),
            "out".into(),
            "--engine".into(),
            "pseudoref".into(),
        ])
        .unwrap();
        assert_eq!(population.commands, ["population"]);
    }

    #[test]
    fn stats_is_a_native_standalone_command() {
        let parsed = parse(&[
            "stats".into(),
            "-f".into(),
            "a".into(),
            "-r".into(),
            "r".into(),
            "-o".into(),
            "o".into(),
            "--stats-count-input-reads".into(),
        ])
        .unwrap();
        assert_eq!(parsed.commands, ["stats"]);
        assert!(parsed.stats_count_input_reads);
    }

    #[test]
    fn gene_expands_to_recovery_stages() {
        let parsed = parse(&[
            "gene".into(),
            "-f".into(),
            "a".into(),
            "-r".into(),
            "r".into(),
            "-o".into(),
            "o".into(),
        ])
        .unwrap();
        assert_eq!(parsed.commands, ["filter", "refilter", "assemble", "gene"]);
    }
    #[test]
    fn uce_default_stages_are_complete() {
        let parsed = parse(&[
            "--assembly-mode".into(),
            "uce".into(),
            "-f".into(),
            "a".into(),
            "-r".into(),
            "r".into(),
            "-o".into(),
            "o".into(),
        ])
        .unwrap();
        assert_eq!(
            parsed.commands,
            ["filter", "refilter", "assemble", "combine", "tree"]
        );
    }
    #[test]
    fn sample_names_match_legacy_rule() {
        assert_eq!(sample_name("foo bar-1"), "Foo_bar_1");
    }

    #[test]
    fn commands_can_follow_options() {
        let parsed = parse(&[
            "--assembly-mode".into(),
            "original".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
            "gene".into(),
        ])
        .unwrap();
        assert_eq!(parsed.commands, ["filter", "refilter", "assemble", "gene"]);
    }

    #[test]
    fn unicode_sample_names_do_not_panic() {
        assert_eq!(sample_name("样本-A"), "样本_a");
    }

    #[test]
    fn equals_style_options_are_parsed() {
        let parsed = parse(&[
            "--assembly-mode=uce".into(),
            "--max-reads=7".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap();
        assert_eq!(
            parsed.commands,
            ["filter", "refilter", "assemble", "combine", "tree"]
        );
        assert_eq!(parsed.max_reads, "7");
    }

    #[test]
    fn invalid_assembly_mode_is_rejected() {
        let error = parse(&[
            "gene".into(),
            "--assembly-mode".into(),
            "typo".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap_err();
        assert!(error.contains("--assembly-mode must be original or uce"));
    }

    #[test]
    fn boolean_options_reject_explicit_values() {
        let error = parse(&[
            "--assembly-mode".into(),
            "uce".into(),
            "--uce-rescue-reads=true".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap_err();
        assert!(error.contains("--uce-rescue-reads does not take a value"));
    }

    #[test]
    fn unknown_options_are_rejected() {
        let error = parse(&[
            "--assembly-mode".into(),
            "uce".into(),
            "--max-read=7".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap_err();
        assert!(error.contains("does not support option '--max-read'"));
    }

    #[test]
    fn incomplete_gene_stage_set_is_rejected() {
        let error = parse(&[
            "gene".into(),
            "filter".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap_err();
        assert!(error.contains("gene requires filter, refilter, and assemble"));
    }

    #[test]
    fn python_compatibility_options_are_all_accepted() {
        let parsed = parse(&[
            "filter".into(),
            "refilter".into(),
            "assemble".into(),
            "--assembly-mode".into(),
            "uce".into(),
            "--assembler-implementation".into(),
            "uce-rust".into(),
            "--assembler-read-chunk-size".into(),
            "4096".into(),
            "--uce-path-strategy".into(),
            "search".into(),
            "--uce-backbone-lookahead".into(),
            "12".into(),
            "--min-depth".into(),
            "1".into(),
            "--max-depth".into(),
            "2".into(),
            "--reuse-reference-cache".into(),
            "--legacy-uce-filter".into(),
            "--workflow-profile".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap();
        assert!(parsed.reuse_reference_cache);
        assert!(parsed.legacy_uce_filter);
        assert!(parsed.workflow_profile);
    }

    #[test]
    fn cleanup_dry_run_is_an_explicit_opt_in() {
        let parsed = parse(&[
            "filter".into(),
            "assemble".into(),
            "--cleanup-intermediates".into(),
            "--cleanup-dry-run".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap();
        assert!(parsed.cleanup_intermediates);
        assert!(parsed.cleanup_dry_run);
        assert!(validate_cleanup_options(&parsed).is_ok());

        let invalid = parse(&[
            "filter".into(),
            "assemble".into(),
            "--cleanup-dry-run".into(),
            "-f".into(),
            "reads.tsv".into(),
            "-r".into(),
            "references".into(),
            "-o".into(),
            "out".into(),
        ])
        .unwrap();
        assert!(validate_cleanup_options(&invalid)
            .unwrap_err()
            .contains("requires --cleanup-intermediates"));
    }

    #[test]
    fn resume_requires_an_exact_successful_workflow() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gm2_resume_fixture_{}_{}",
            std::process::id(),
            unique
        ));
        fs::create_dir_all(&root).unwrap();
        let reference = root.join("reference.fasta");
        let reads = root.join("reads.fastq");
        let table = root.join("samples.tsv");
        let output = root.join("output");
        fs::write(&reference, ">reference\nACGT\n").unwrap();
        fs::write(&reads, "@read\nACGT\n+\n!!!!\n").unwrap();
        fs::write(&table, format!("sample\t{}\n", reads.display())).unwrap();
        fs::create_dir_all(&output).unwrap();
        let initial_options = parse(&[
            "filter".into(),
            "assemble".into(),
            "-f".into(),
            table.display().to_string(),
            "-r".into(),
            reference.display().to_string(),
            "-o".into(),
            output.display().to_string(),
        ])
        .unwrap();
        let samples = read_samples(&initial_options.samples, &output).unwrap();
        write_workflow_manifest(&initial_options, &samples).unwrap();
        write_workflow_status(&output, &initial_options.commands, &Ok(())).unwrap();
        let resume_options = parse(&[
            "filter".into(),
            "assemble".into(),
            "--resume".into(),
            "-f".into(),
            table.display().to_string(),
            "-r".into(),
            reference.display().to_string(),
            "-o".into(),
            output.display().to_string(),
        ])
        .unwrap();
        assert!(resume_completed_workflow(&resume_options, &samples).unwrap());
        fs::write(
            output.join("workflow_status.tsv"),
            "field\tvalue\nstate\tfailed\n",
        )
        .unwrap();
        assert!(resume_completed_workflow(&resume_options, &samples)
            .unwrap_err()
            .contains("did not complete successfully"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn unsupported_stage_is_not_silently_ignored() {
        let error = parse(&[
            "rescue".into(),
            "--assembly-mode".into(),
            "uce".into(),
            "-f".into(),
            "a".into(),
            "-r".into(),
            "r".into(),
            "-o".into(),
            "o".into(),
        ])
        .unwrap_err();
        assert!(error.contains("does not support command 'rescue'"));
    }
}
