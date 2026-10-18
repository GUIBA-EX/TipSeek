mod assembly;
mod hash;
mod io_utils;
mod model;
mod pipeline;
mod seq;
mod unitig;

use crate::hash::HashSet;
use io_utils::{discover_references, find_filtered};
use model::{Args, AssemblyMode, GraphFormat, PathStrategy};
use pipeline::{
    log_line, process_locus, read_result_dict, read_summary_lines, run_manifest, summary_line,
    write_assembly_profile, write_result_dict, write_summary, ProfileStats,
};
use std::env;
use std::io;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{mpsc, Arc, Mutex};
use std::time::Instant;

const HELP: &str = "GeneMiner2 Rust assembler

Usage:
  main_assembler -r PATH -o DIR [options]

Required:
  -r PATH                         Reference FASTA or directory
  -o DIR                          Sample output directory containing filtered/

Assembly:
  -ka INT                         Assembly k-mer; 0 selects automatically (default: 39)
  -k_min INT                      Minimum automatic k-mer (default: 21)
  -k_max INT                      Maximum automatic k-mer, at most 63 (default: 39)
  -limit_count INT                K-mer error threshold (default: 2)
  -iteration INT                  Search/extension limit (default: 8192)
  -cov_min FLOAT                  Minimum contig coverage (default: 0)
  -sb, --soft_boundary INT        Reference soft boundary; -1 uses half slice (default: 0)
  -p, --processes INT             Parallel locus workers (default: 1)

Assembly mode:
  --assembly-mode original|uce|its2
                                   Assembly mode (default: original)
  --uce-path-strategy search|backbone
                                   UCE path handling (default: backbone)
  --uce-backbone-lookahead INT    Bounded branch look-ahead (default: 24)
  --uce-reverse-reuse-reference-scale FLOAT
                                   Scale reference bonus on reverse path reuse
                                   (default: 1.0; no scaling)
  --uce-side-candidates INT       Legacy search candidates per side (default: 8)
  --uce-max-contig-length INT     UCE length guardrail; 0 disables (default: 0)
  --uce-min-read-density FLOAT    Long-contig unique-read density (default: 0.003)
  --uce-density-check-min-length INT
                                   Length where density guardrail applies (default: 1000)
  --uce-max-depth-cv FLOAT        Depth CV guardrail; 0 disables
  --uce-max-depth-ratio FLOAT     Max/median depth guardrail; 0 disables
  --assembler-reference-cache-dir DIR
                                   Versioned Rust reference k-mer cache
  --assembler-read-chunk-size INT  Reads per streaming batch (default: 8192)
  --assembler-kmer-count-threads INT
                                   Sort/count workers per locus; 0=auto
  --assembler-graph-format none|gfa|dot|both
                                   Write compact assembly graphs (default: none)
  --profile                        Write aggregate stage timings to assembly_profile.tsv

Other:
  -h, --help                      Show this help
  -V, --version                   Show version
";

// 命令行参数一个萝卜一个坑，缺值就别往下凑合。
fn next_value(arguments: &[String], index: &mut usize, flag: &str) -> Result<String, String> {
    *index += 1;
    arguments
        .get(*index)
        .cloned()
        .ok_or_else(|| format!("{flag} requires a value"))
}

// 数值开关在入口统一验，免得组装跑半道才发现参数不对。
fn parse_number<T: std::str::FromStr>(
    arguments: &[String],
    index: &mut usize,
    flag: &str,
) -> Result<T, String> {
    let value = next_value(arguments, index, flag)?;
    value
        .parse()
        .map_err(|_| format!("invalid value for {flag}: {value}"))
}

// 把 CLI 配置归到 Args，默认走 UCE backbone，用户指定才换路子。
fn parse_args() -> Result<Args, String> {
    let arguments: Vec<String> = env::args().collect();
    let mut reference = None;
    let mut output = None;
    let mut args = Args {
        reference: PathBuf::new(),
        output: PathBuf::new(),
        kmer_size: 39,
        kmer_min: 21,
        kmer_max: 39,
        error_limit: 2,
        iteration: 8192,
        min_coverage: 0.0,
        soft_boundary: 0,
        threads: 1,
        assembly_mode: AssemblyMode::Reference,
        side_candidates: 8,
        path_strategy: PathStrategy::Backbone,
        backbone_lookahead: 24,
        reverse_reuse_reference_scale: 1.0,
        max_contig_length: 0,
        min_read_density: 0.003,
        density_check_min_length: 1000,
        max_depth_cv: 0.0,
        max_depth_ratio: 0.0,
        reference_cache_dir: None,
        read_chunk_size: 8192,
        kmer_count_threads: 0,
        graph_format: GraphFormat::None,
        profile: false,
    };

    let mut index = 1;
    while index < arguments.len() {
        let flag = arguments[index].as_str();
        match flag {
            "-h" | "--help" => {
                print!("{HELP}");
                std::process::exit(0);
            }
            "-V" | "--version" => {
                println!("main_assembler 0.7.0");
                std::process::exit(0);
            }
            "-r" => reference = Some(PathBuf::from(next_value(&arguments, &mut index, flag)?)),
            "-o" => output = Some(PathBuf::from(next_value(&arguments, &mut index, flag)?)),
            "-ka" => args.kmer_size = parse_number(&arguments, &mut index, flag)?,
            "-k_min" => args.kmer_min = parse_number(&arguments, &mut index, flag)?,
            "-k_max" => args.kmer_max = parse_number(&arguments, &mut index, flag)?,
            "-limit_count" => args.error_limit = parse_number(&arguments, &mut index, flag)?,
            "-iteration" => args.iteration = parse_number(&arguments, &mut index, flag)?,
            "-cov_min" => args.min_coverage = parse_number(&arguments, &mut index, flag)?,
            "-sb" | "--soft_boundary" => {
                args.soft_boundary = parse_number(&arguments, &mut index, flag)?
            }
            "-p" | "--processes" => args.threads = parse_number(&arguments, &mut index, flag)?,
            "--assembly-mode" => {
                args.assembly_mode = match next_value(&arguments, &mut index, flag)?.as_str() {
                    "original" => AssemblyMode::Reference,
                    "uce" => AssemblyMode::Uce,
                    "its2" => AssemblyMode::Its2,
                    value => return Err(format!("invalid --assembly-mode: {value}")),
                }
            }
            "--uce-path-strategy" => {
                args.path_strategy = match next_value(&arguments, &mut index, flag)?.as_str() {
                    "search" => PathStrategy::Search,
                    "backbone" => PathStrategy::Backbone,
                    value => return Err(format!("invalid --uce-path-strategy: {value}")),
                }
            }
            "--uce-backbone-lookahead" => {
                args.backbone_lookahead = parse_number(&arguments, &mut index, flag)?
            }
            "--uce-reverse-reuse-reference-scale" => {
                args.reverse_reuse_reference_scale = parse_number(&arguments, &mut index, flag)?
            }
            "--uce-side-candidates" => {
                args.side_candidates = parse_number(&arguments, &mut index, flag)?
            }
            "--uce-max-contig-length" => {
                args.max_contig_length = parse_number(&arguments, &mut index, flag)?
            }
            "--uce-min-read-density" => {
                args.min_read_density = parse_number(&arguments, &mut index, flag)?
            }
            "--uce-density-check-min-length" => {
                args.density_check_min_length = parse_number(&arguments, &mut index, flag)?
            }
            "--uce-max-depth-cv" => args.max_depth_cv = parse_number(&arguments, &mut index, flag)?,
            "--uce-max-depth-ratio" => {
                args.max_depth_ratio = parse_number(&arguments, &mut index, flag)?
            }
            "--assembler-reference-cache-dir" => {
                args.reference_cache_dir =
                    Some(PathBuf::from(next_value(&arguments, &mut index, flag)?))
            }
            "--assembler-read-chunk-size" => {
                args.read_chunk_size = parse_number(&arguments, &mut index, flag)?
            }
            "--assembler-kmer-count-threads" => {
                args.kmer_count_threads = parse_number(&arguments, &mut index, flag)?
            }
            "--profile" => args.profile = true,
            "--assembler-graph-format" => {
                args.graph_format = match next_value(&arguments, &mut index, flag)?.as_str() {
                    "none" => GraphFormat::None,
                    "gfa" => GraphFormat::Gfa,
                    "dot" => GraphFormat::Dot,
                    "both" => GraphFormat::Both,
                    value => return Err(format!("invalid --assembler-graph-format: {value}")),
                }
            }
            unknown => return Err(format!("unknown argument: {unknown}")),
        }
        index += 1;
    }

    args.reference = reference.ok_or_else(|| "-r is required".to_string())?;
    args.output = output.ok_or_else(|| "-o is required".to_string())?;
    args.threads = args.threads.max(1);
    args.side_candidates = args.side_candidates.max(3);
    args.backbone_lookahead = args.backbone_lookahead.max(1);
    args.density_check_min_length = args.density_check_min_length.max(1);
    args.read_chunk_size = args.read_chunk_size.max(1);

    if args.assembly_mode == AssemblyMode::Its2 {
        if !arguments.iter().any(|arg| arg == "-ka") {
            args.kmer_size = 21;
        }
        if !arguments.iter().any(|arg| arg == "-k_min") {
            args.kmer_min = 21;
        }
        if !arguments.iter().any(|arg| arg == "-k_max") {
            args.kmer_max = 21;
        }
    }

    if args.kmer_size > 63 || args.kmer_min > 63 || args.kmer_max > 63 {
        return Err("Rust u128 assembler supports k-mer sizes up to 63".to_string());
    }
    if args.kmer_min == 0 || args.kmer_max < args.kmer_min {
        return Err("invalid automatic k-mer range".to_string());
    }
    if args.min_read_density < 0.0 || args.max_depth_cv < 0.0 || args.max_depth_ratio < 0.0 {
        return Err("UCE guardrail values must be non-negative".to_string());
    }
    if !(0.0..=1.0).contains(&args.reverse_reuse_reference_scale) {
        return Err("--uce-reverse-reuse-reference-scale must be in [0, 1]".to_string());
    }
    Ok(args)
}

// 发现全部 locus 后分派工人；每个 locus 独立干活儿，结果再统一汇总。
fn run(mut args: Args) -> io::Result<()> {
    std::fs::create_dir_all(args.output.join("results"))?;
    std::fs::create_dir_all(args.output.join("contigs_all"))?;
    std::fs::create_dir_all(args.output.join("contigs_all_low"))?;

    let log_lock = Arc::new(Mutex::new(()));
    log_line(
        &args.output,
        &log_lock,
        "======================== Assemble =========================",
    );
    let started = Instant::now();
    let tasks = discover_references(&args.reference)?;
    let valid_keys: HashSet<String> = tasks.iter().map(|task| task.key.clone()).collect();

    let result_path = args.output.join("result_dict.txt");
    let manifest_path = args.output.join("assembly_run_manifest.txt");
    let manifest = run_manifest(&args, &tasks)?;
    let resume_safe =
        std::fs::read_to_string(&manifest_path).is_ok_and(|previous| previous == manifest);
    std::fs::write(&manifest_path, &manifest)?;
    if !resume_safe && result_path.exists() {
        log_line(
            &args.output,
            &log_lock,
            "Assembly inputs or parameters changed; ignoring prior completion state.",
        );
    }
    let summary_path = args
        .output
        .join(if args.assembly_mode == AssemblyMode::Its2 {
            "its2_assembly_summary.csv"
        } else {
            "uce_assembly_summary.csv"
        });
    let mut result_dict = if resume_safe {
        read_result_dict(&result_path)?
    } else {
        Default::default()
    };
    result_dict.retain(|key, _| valid_keys.contains(key));
    let mut summary_rows =
        if resume_safe && matches!(args.assembly_mode, AssemblyMode::Uce | AssemblyMode::Its2) {
            read_summary_lines(&summary_path)?
        } else {
            Default::default()
        };
    summary_rows.retain(|key, _| valid_keys.contains(key));

    let mut completed: HashSet<String> = result_dict.keys().cloned().collect();
    if matches!(args.assembly_mode, AssemblyMode::Uce | AssemblyMode::Its2) {
        completed.retain(|key| summary_rows.contains_key(key));
    }

    // UCE loci are independent and highly uneven in recruited-read volume.
    // Run the largest ones first to reduce parallel tail time; ordinal remains
    // the reference order for deterministic user-facing labels.
    let mut scheduled_tasks: Vec<_> = tasks
        .into_iter()
        .map(|task| {
            let bytes = find_filtered(&args.output, &task.key)
                .and_then(|(path, _)| std::fs::metadata(path).ok())
                .map_or(0, |metadata| metadata.len());
            (task, bytes)
        })
        .collect();
    scheduled_tasks.sort_by(|(left, left_bytes), (right, right_bytes)| {
        right_bytes
            .cmp(left_bytes)
            .then_with(|| left.ordinal.cmp(&right.ordinal))
    });
    let tasks = Arc::new(
        scheduled_tasks
            .into_iter()
            .map(|(task, _)| task)
            .collect::<Vec<_>>(),
    );
    let completed = Arc::new(completed);
    let next_task = Arc::new(AtomicUsize::new(0));
    let worker_count = args.threads.min(tasks.len().max(1));
    let default_kmer_workers = (args.threads / worker_count).max(1);
    if args.kmer_count_threads == 0 {
        args.kmer_count_threads = default_kmer_workers;
    }
    let machine_threads = std::thread::available_parallelism()
        .map(|count| count.get())
        .unwrap_or(args.threads.max(1));
    let per_locus_cap = (machine_threads / worker_count).max(1);
    if args.kmer_count_threads > per_locus_cap {
        log_line(
            &args.output,
            &log_lock,
            &format!(
                "Cap k-mer count workers per locus at {} so {} locus workers stay within {} logical CPUs.",
                per_locus_cap, worker_count, machine_threads
            ),
        );
        args.kmer_count_threads = per_locus_cap;
    }
    let profile = args.profile.then(|| Arc::new(ProfileStats::default()));
    let (sender, receiver) = mpsc::channel();

    std::thread::scope(|scope| {
        for _ in 0..worker_count {
            let tasks = Arc::clone(&tasks);
            let completed = Arc::clone(&completed);
            let next_task = Arc::clone(&next_task);
            let sender = sender.clone();
            let log_lock = Arc::clone(&log_lock);
            let args = &args;
            let profile = profile.as_ref().map(Arc::clone);
            scope.spawn(move || {
                gm2_tools::platform::prefer_performance_cores();
                loop {
                    let index = next_task.fetch_add(1, Ordering::Relaxed);
                    let Some(task) = tasks.get(index) else {
                        break;
                    };
                    let result =
                        process_locus(args, task, &completed, &log_lock, profile.as_deref());
                    if sender.send(result).is_err() {
                        break;
                    }
                }
            });
        }
    });
    drop(sender);

    for result in receiver {
        if result.skipped {
            continue;
        }
        result_dict.insert(result.key.clone(), (result.status.clone(), result.value));
        if matches!(args.assembly_mode, AssemblyMode::Uce | AssemblyMode::Its2) {
            summary_rows.insert(result.key.clone(), summary_line(&result));
        }
    }

    write_result_dict(&result_path, &result_dict)?;
    if let Some(profile) = profile.as_deref() {
        write_assembly_profile(&args.output, profile)?;
    }
    if matches!(args.assembly_mode, AssemblyMode::Uce | AssemblyMode::Its2) {
        write_summary(&summary_path, &summary_rows)?;
    }
    log_line(
        &args.output,
        &log_lock,
        &format!("\nTime cost: {:.6}\n", started.elapsed().as_secs_f64()),
    );
    Ok(())
}

// 主程序只管报错和退出码，真正的活儿交给 run。
fn main() {
    let args = match parse_args() {
        Ok(args) => args,
        Err(error) => {
            eprintln!("error: {error}\n\n{HELP}");
            std::process::exit(2);
        }
    };
    if let Err(error) = run(args) {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}
