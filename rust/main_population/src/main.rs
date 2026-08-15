mod panref {
    pub(crate) mod backbone;
    pub(crate) mod bait;
    pub(crate) mod bait_index;
    pub(crate) mod bubble;
    pub(crate) mod dbg;
    pub(crate) mod recruit;
    pub(crate) mod v2;
}

use crate::panref::backbone::assemble_backbone;
use crate::panref::bait::BaitCatalog;
use crate::panref::bait_index::BaitIndex;
use crate::panref::recruit::recruit_pairs_to_fastq;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::env;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{self, Command, Stdio};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;

type AppResult<T> = Result<T, String>;
type PanrefBackbone = (Vec<u8>, u64);

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
enum Stage {
    Reference,
    Mapping,
    Calling,
    Selection,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ReferenceStrategy {
    SqclLongest,
    Supported,
}

impl ReferenceStrategy {
    fn label(self) -> &'static str {
        match self {
            Self::SqclLongest => "sqcl-longest",
            Self::Supported => "supported",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Engine {
    Pseudoref,
    Panref,
    PanrefV2,
}

#[derive(Debug)]
struct Args {
    output: PathBuf,
    samples_tsv: PathBuf,
    engine: Engine,
    reference_strategy: ReferenceStrategy,
    panref_baits: Option<PathBuf>,
    panrefv2_include_low_confidence: bool,
    threads: usize,
    min_mapq: u32,
    min_baseq: u32,
    min_dp: u32,
    min_gq: u32,
    min_qual: f64,
    min_call_rate: f64,
    min_mac: u32,
    ld_window: usize,
    ld_step: usize,
    ld_r2: f64,
    mark_duplicates: bool,
    reference_fasta: Option<PathBuf>,
    start_at: Stage,
    stop_after: Stage,
    skip_plink: bool,
    skip_relatedness_qc: bool,
    skip_admixture: bool,
    admixture_k_min: usize,
    admixture_k_max: usize,
    admixture_cv: usize,
    minibwa: String,
    samtools: String,
    bcftools: String,
    plink: String,
    admixture: String,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            output: PathBuf::new(),
            samples_tsv: PathBuf::new(),
            engine: Engine::Pseudoref,
            reference_strategy: ReferenceStrategy::SqclLongest,
            panref_baits: None,
            panrefv2_include_low_confidence: false,
            threads: 1,
            min_mapq: 20,
            min_baseq: 20,
            min_dp: 5,
            min_gq: 20,
            min_qual: 20.0,
            min_call_rate: 0.8,
            min_mac: 2,
            ld_window: 50,
            ld_step: 5,
            ld_r2: 0.2,
            mark_duplicates: true,
            reference_fasta: None,
            start_at: Stage::Reference,
            stop_after: Stage::Selection,
            skip_plink: false,
            skip_relatedness_qc: false,
            skip_admixture: false,
            admixture_k_min: 2,
            admixture_k_max: 6,
            admixture_cv: 10,
            minibwa: "minibwa".into(),
            samtools: "samtools".into(),
            bcftools: "bcftools".into(),
            plink: "plink".into(),
            admixture: "admixture".into(),
        }
    }
}

#[derive(Clone, Debug)]
struct Sample {
    original: String,
    internal: String,
    vcf_id: String,
    read1: PathBuf,
    read2: PathBuf,
    population: String,
    batch: String,
}

#[derive(Clone, Debug)]
struct Candidate {
    locus: String,
    safe_locus: String,
    sample: String,
    sequence: String,
    supported_bases: u64,
    support_breadth: f64,
    max_gap: u64,
    unique_reads: u64,
    unique_density: f64,
    length: usize,
}

#[derive(Clone, Debug)]
struct Site {
    chrom: String,
    pos: u64,
    reference: String,
    alternate: String,
    qual: f64,
    call_rate: f64,
    median_dp: f64,
    median_gq: f64,
    record: String,
}

fn print_help() {
    println!(
        "main_population (TipSeek Rust population pipeline)\n\
         Usage: main_population --output DIR --samples FILE [options]\n\n\
         --output DIR              Existing GeneMiner2 output directory\n\
         --samples FILE            TSV: sample, R1, [R2, population, batch]\n\
         --engine STR              pseudoref, panref, or panrefv2 (default: pseudoref)\n\
         --panref-baits DIR       Per-locus bait FASTA directory required by panref or panrefv2\n\
         --panrefv2-include-low-confidence  Include short or low-support PanRefV2 loci in FASTA\n\
         --reference-strategy STR  sqcl-longest or supported (default: sqcl-longest)\n\
         --reference-fasta FILE    Use a fixed external cohort reference\n\
         --threads INT             Threads for external tools and PanRef graph building (default: 1)\n\
         --min-mapq INT            Minimum mapping quality (default: 20)\n\
         --min-baseq INT           Minimum base quality (default: 20)\n\
         --min-dp INT              Minimum genotype depth (default: 5)\n\
         --min-gq INT              Minimum genotype quality (default: 20)\n\
         --min-qual FLOAT          Minimum site QUAL (default: 20)\n\
         --min-call-rate FLOAT     Minimum non-missing fraction (default: 0.8)\n\
         --min-mac INT             Minimum minor allele count (default: 2)\n\
         --ld-window INT           SNPs per LD window (default: 50)\n\
         --ld-step INT             SNPs shifted per LD window (default: 5)\n\
         --ld-r2 FLOAT             LD pruning r^2 threshold (default: 0.2)\n\
         --skip-mark-duplicates    Skip samtools fixmate/markdup\n\
         --skip-plink              Do not create PLINK/PCA/LD-pruned panels\n\
         --skip-relatedness-qc     Skip PLINK pairwise relatedness calculation\n\
         --skip-admixture          Do not run ADMIXTURE on the primary panel\n\
         --admixture-k-min INT     Minimum ADMIXTURE K (default: 2)\n\
         --admixture-k-max INT     Maximum ADMIXTURE K (default: 6)\n\
         --admixture-cv INT        Cross-validation folds (default: 10)\n\
         --start-at STAGE          reference, mapping, calling, or selection\n\
         --stop-after STAGE        reference, mapping, calling, or selection\n\
         --minibwa PATH            minibwa executable\n\
         --samtools PATH           samtools executable\n\
         --bcftools PATH           bcftools executable\n\
         --plink PATH              PLINK 1.9 executable\n\
         --admixture PATH          ADMIXTURE executable\n\
         --version                 Print version"
    );
}

fn take_value(argv: &[String], index: &mut usize, option: &str) -> AppResult<String> {
    *index += 1;
    argv.get(*index)
        .cloned()
        .ok_or_else(|| format!("option {option} requires an argument"))
}

fn parse_num<T: std::str::FromStr>(value: String, option: &str) -> AppResult<T> {
    value
        .parse::<T>()
        .map_err(|_| format!("invalid value {value:?} for {option}"))
}

fn parse_stage(value: &str) -> AppResult<Stage> {
    match value {
        "reference" => Ok(Stage::Reference),
        "mapping" => Ok(Stage::Mapping),
        "calling" => Ok(Stage::Calling),
        "selection" => Ok(Stage::Selection),
        _ => Err(format!(
            "invalid stage {value:?}; expected reference, mapping, calling, or selection"
        )),
    }
}

fn parse_reference_strategy(value: &str) -> AppResult<ReferenceStrategy> {
    match value {
        "sqcl-longest" | "sqcl" => Ok(ReferenceStrategy::SqclLongest),
        "supported" => Ok(ReferenceStrategy::Supported),
        _ => Err(format!(
            "invalid reference strategy {value:?}; expected sqcl-longest or supported"
        )),
    }
}

fn parse_engine(value: &str) -> AppResult<Engine> {
    match value {
        "pseudoref" => Ok(Engine::Pseudoref),
        "panref" => Ok(Engine::Panref),
        "panrefv2" => Ok(Engine::PanrefV2),
        _ => Err("invalid engine; expected pseudoref, panref, or panrefv2".to_string()),
    }
}

fn parse_args(argv: Vec<String>) -> AppResult<Args> {
    if argv.iter().any(|arg| arg == "-h" || arg == "--help") {
        print_help();
        process::exit(0);
    }
    if argv.iter().any(|arg| arg == "--version") {
        println!("main_population 0.5.0");
        process::exit(0);
    }

    let mut args = Args::default();
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "--output" => args.output = PathBuf::from(take_value(&argv, &mut i, "--output")?),
            "--samples" => {
                args.samples_tsv = PathBuf::from(take_value(&argv, &mut i, "--samples")?)
            }
            "--engine" => args.engine = parse_engine(&take_value(&argv, &mut i, "--engine")?)?,
            "--panref-baits" => {
                args.panref_baits =
                    Some(PathBuf::from(take_value(&argv, &mut i, "--panref-baits")?))
            }
            "--panrefv2-include-low-confidence" => args.panrefv2_include_low_confidence = true,
            "--reference-strategy" => {
                args.reference_strategy =
                    parse_reference_strategy(&take_value(&argv, &mut i, "--reference-strategy")?)?
            }
            "--reference-fasta" => {
                args.reference_fasta = Some(PathBuf::from(take_value(
                    &argv,
                    &mut i,
                    "--reference-fasta",
                )?))
            }
            "--threads" => {
                args.threads = parse_num(take_value(&argv, &mut i, "--threads")?, "--threads")?
            }
            "--min-mapq" => {
                args.min_mapq = parse_num(take_value(&argv, &mut i, "--min-mapq")?, "--min-mapq")?
            }
            "--min-baseq" => {
                args.min_baseq =
                    parse_num(take_value(&argv, &mut i, "--min-baseq")?, "--min-baseq")?
            }
            "--min-dp" => {
                args.min_dp = parse_num(take_value(&argv, &mut i, "--min-dp")?, "--min-dp")?
            }
            "--min-gq" => {
                args.min_gq = parse_num(take_value(&argv, &mut i, "--min-gq")?, "--min-gq")?
            }
            "--min-qual" => {
                args.min_qual = parse_num(take_value(&argv, &mut i, "--min-qual")?, "--min-qual")?
            }
            "--min-call-rate" => {
                args.min_call_rate = parse_num(
                    take_value(&argv, &mut i, "--min-call-rate")?,
                    "--min-call-rate",
                )?
            }
            "--min-mac" => {
                args.min_mac = parse_num(take_value(&argv, &mut i, "--min-mac")?, "--min-mac")?
            }
            "--ld-window" => {
                args.ld_window =
                    parse_num(take_value(&argv, &mut i, "--ld-window")?, "--ld-window")?
            }
            "--ld-step" => {
                args.ld_step = parse_num(take_value(&argv, &mut i, "--ld-step")?, "--ld-step")?
            }
            "--ld-r2" => args.ld_r2 = parse_num(take_value(&argv, &mut i, "--ld-r2")?, "--ld-r2")?,
            "--skip-mark-duplicates" => args.mark_duplicates = false,
            "--skip-plink" => args.skip_plink = true,
            "--skip-relatedness-qc" => args.skip_relatedness_qc = true,
            "--skip-admixture" => args.skip_admixture = true,
            "--admixture-k-min" => {
                args.admixture_k_min = parse_num(
                    take_value(&argv, &mut i, "--admixture-k-min")?,
                    "--admixture-k-min",
                )?
            }
            "--admixture-k-max" => {
                args.admixture_k_max = parse_num(
                    take_value(&argv, &mut i, "--admixture-k-max")?,
                    "--admixture-k-max",
                )?
            }
            "--admixture-cv" => {
                args.admixture_cv = parse_num(
                    take_value(&argv, &mut i, "--admixture-cv")?,
                    "--admixture-cv",
                )?
            }
            "--stop-after" => {
                args.stop_after = parse_stage(&take_value(&argv, &mut i, "--stop-after")?)?
            }
            "--start-at" => args.start_at = parse_stage(&take_value(&argv, &mut i, "--start-at")?)?,
            "--minibwa" => args.minibwa = take_value(&argv, &mut i, "--minibwa")?,
            "--samtools" => args.samtools = take_value(&argv, &mut i, "--samtools")?,
            "--bcftools" => args.bcftools = take_value(&argv, &mut i, "--bcftools")?,
            "--plink" => args.plink = take_value(&argv, &mut i, "--plink")?,
            "--admixture" => args.admixture = take_value(&argv, &mut i, "--admixture")?,
            unknown => return Err(format!("unknown option {unknown}")),
        }
        i += 1;
    }

    if args.output.as_os_str().is_empty() || args.samples_tsv.as_os_str().is_empty() {
        return Err("--output and --samples are required".into());
    }
    if args.threads == 0 {
        return Err("--threads must be at least 1".into());
    }
    if !(0.0 < args.min_call_rate && args.min_call_rate <= 1.0) {
        return Err("--min-call-rate must be in (0, 1]".into());
    }
    if args.ld_window == 0 || args.ld_step == 0 || !(0.0 < args.ld_r2 && args.ld_r2 < 1.0) {
        return Err("LD window/step must be positive and --ld-r2 must be in (0, 1)".into());
    }
    if args.admixture_k_min < 2
        || args.admixture_k_max < args.admixture_k_min
        || args.admixture_cv < 2
    {
        return Err("ADMIXTURE requires 2 <= K-min <= K-max and at least 2 CV folds".into());
    }
    if args.start_at > args.stop_after {
        return Err("--start-at must not be later than --stop-after".into());
    }
    if matches!(args.engine, Engine::Panref | Engine::PanrefV2) {
        if args.reference_fasta.is_some() {
            return Err(
                "--engine panref/panrefv2 cannot be combined with --reference-fasta".into(),
            );
        }
        let Some(baits) = args.panref_baits.as_ref() else {
            return Err("--engine panref/panrefv2 requires --panref-baits".into());
        };
        if !baits.is_dir() {
            return Err(format!(
                "--panref-baits is not a directory: {}",
                baits.display()
            ));
        }
    }
    Ok(args)
}

fn parse_delimited_line(line: &str, delimiter: char) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let mut chars = line.trim_end_matches(['\n', '\r']).chars().peekable();
    while let Some(ch) = chars.next() {
        if ch == '"' {
            if quoted && chars.peek() == Some(&'"') {
                field.push('"');
                chars.next();
            } else {
                quoted = !quoted;
            }
        } else if ch == delimiter && !quoted {
            fields.push(field);
            field = String::new();
        } else {
            field.push(ch);
        }
    }
    fields.push(field);
    fields
}

fn python_capitalize(value: &str) -> String {
    let lowered = value.to_lowercase();
    let mut chars = lowered.chars();
    match chars.next() {
        None => String::new(),
        Some(first) => first.to_uppercase().collect::<String>() + chars.as_str(),
    }
}

fn sanitize_sample_name(value: &str) -> String {
    let filtered: String = value
        .trim()
        .chars()
        .filter_map(|ch| {
            if ch.is_alphanumeric() || "-_. ".contains(ch) {
                Some(if ch == ' ' || ch == '-' { '_' } else { ch })
            } else {
                None
            }
        })
        .collect();
    python_capitalize(&filtered)
}

fn sanitize_vcf_id(value: &str) -> String {
    let cleaned: String = value
        .trim()
        .chars()
        .map(|ch| {
            if ch.is_alphanumeric() || ch == '_' || ch == '.' || ch == '-' {
                ch
            } else {
                '_'
            }
        })
        .collect();
    if cleaned.is_empty() {
        "sample".into()
    } else {
        cleaned
    }
}

fn read_samples(path: &Path) -> AppResult<Vec<Sample>> {
    let file = File::open(path)
        .map_err(|e| format!("unable to read sample list {}: {e}", path.display()))?;
    let mut samples = Vec::new();
    let mut used_vcf_ids = HashSet::new();
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|e| format!("unable to read sample list: {e}"))?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = parse_delimited_line(&line, '\t');
        if fields.len() < 2 {
            return Err(format!(
                "sample row requires a name and at least one read file: {line}"
            ));
        }
        let number = samples.len() + 1;
        let safe = sanitize_sample_name(&fields[0]);
        if safe.is_empty() {
            return Err(format!("invalid sample name {:?}", fields[0]));
        }
        let mut vcf_id = sanitize_vcf_id(&fields[0]);
        if used_vcf_ids.contains(&vcf_id) {
            vcf_id = format!("{vcf_id}_{number}");
        }
        used_vcf_ids.insert(vcf_id.clone());
        samples.push(Sample {
            original: fields[0].clone(),
            internal: format!("{number}_{safe}"),
            vcf_id,
            read1: PathBuf::from(&fields[1]),
            read2: PathBuf::from(fields.get(2).unwrap_or(&fields[1])),
            population: fields
                .get(3)
                .map(|value| value.trim())
                .unwrap_or("")
                .to_string(),
            batch: fields
                .get(4)
                .map(|value| value.trim())
                .unwrap_or("")
                .to_string(),
        });
    }
    if samples.is_empty() {
        return Err("sample list is empty".into());
    }
    Ok(samples)
}

fn parse_csv_records(path: &Path) -> AppResult<Vec<HashMap<String, String>>> {
    let file = File::open(path).map_err(|e| format!("unable to read {}: {e}", path.display()))?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("empty CSV {}", path.display()))?
        .map_err(|e| format!("unable to read {}: {e}", path.display()))?;
    let headers = parse_delimited_line(&header, ',');
    let mut records = Vec::new();
    for line in lines {
        let line = line.map_err(|e| format!("unable to read {}: {e}", path.display()))?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = parse_delimited_line(&line, ',');
        let mut row = HashMap::new();
        for (index, name) in headers.iter().enumerate() {
            row.insert(name.clone(), fields.get(index).cloned().unwrap_or_default());
        }
        records.push(row);
    }
    Ok(records)
}

fn truthy(value: Option<&String>) -> bool {
    matches!(
        value.map(|v| v.trim().to_ascii_lowercase()).as_deref(),
        Some("1" | "true" | "yes")
    )
}

fn accepted(row: &HashMap<String, String>) -> bool {
    if row
        .get("accepted")
        .is_some_and(|value| !value.trim().is_empty())
    {
        return truthy(row.get("accepted"));
    }
    row.get("status").is_some_and(|v| v == "success") && !truthy(row.get("low_quality"))
}

fn number<T: std::str::FromStr + Default>(row: &HashMap<String, String>, key: &str) -> T {
    row.get(key)
        .and_then(|v| v.parse::<T>().ok())
        .unwrap_or_default()
}

fn read_first_fasta(path: &Path) -> AppResult<String> {
    let file = File::open(path).map_err(|e| format!("unable to read {}: {e}", path.display()))?;
    let mut sequence = String::new();
    let mut seen_header = false;
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|e| format!("unable to read {}: {e}", path.display()))?;
        if line.starts_with('>') {
            if seen_header && !sequence.is_empty() {
                break;
            }
            seen_header = true;
        } else if seen_header {
            sequence.push_str(line.trim());
        }
    }
    let sequence = sequence.to_ascii_uppercase();
    if sequence.is_empty() {
        return Err(format!("no FASTA sequence in {}", path.display()));
    }
    if !sequence.chars().all(|ch| "ACGTRYSWKMBDHVN".contains(ch)) {
        return Err(format!("invalid nucleotide in {}", path.display()));
    }
    Ok(sequence)
}

fn safe_locus_name(locus: &str) -> String {
    let value: String = locus
        .chars()
        .map(|ch| {
            if ch.is_alphanumeric() || ch == '_' || ch == '.' || ch == '-' {
                ch
            } else {
                '_'
            }
        })
        .collect();
    if value.is_empty() {
        "locus".into()
    } else {
        value
    }
}

fn support_cmp(left: &Candidate, right: &Candidate) -> Ordering {
    left.supported_bases
        .cmp(&right.supported_bases)
        .then_with(|| {
            left.support_breadth
                .partial_cmp(&right.support_breadth)
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| right.max_gap.cmp(&left.max_gap))
        .then_with(|| {
            left.unique_density
                .partial_cmp(&right.unique_density)
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| left.unique_reads.cmp(&right.unique_reads))
}

fn candidate_cmp(strategy: ReferenceStrategy, left: &Candidate, right: &Candidate) -> Ordering {
    match strategy {
        // SqCL's make_PRG.py chooses the longest contig among acceptable
        // target matches. GeneMiner2's accepted flag and UCE guardrails are
        // the equivalent eligibility filter; support metrics break ties.
        ReferenceStrategy::SqclLongest => left
            .length
            .cmp(&right.length)
            .then_with(|| support_cmp(left, right)),
        ReferenceStrategy::Supported => {
            support_cmp(left, right).then_with(|| left.length.cmp(&right.length))
        }
    }
    .then_with(|| right.sample.cmp(&left.sample))
}

fn build_reference(args: &Args, samples: &[Sample]) -> AppResult<PathBuf> {
    match args.engine {
        Engine::Pseudoref => build_pseudoref_reference(args, samples),
        Engine::Panref => build_panref_reference(args, samples),
        Engine::PanrefV2 => panref::v2::build_reference(args, samples),
    }
}

fn assemble_panref_backbones(
    recruited_dir: &Path,
    catalog: &BaitCatalog,
    threads: usize,
) -> AppResult<Vec<Option<PanrefBackbone>>> {
    let tasks = catalog
        .loci
        .iter()
        .enumerate()
        .filter_map(|(id, locus)| {
            let path = recruited_dir.join(format!("locus_{id:05}.interleaved.fq"));
            path.is_file().then_some((id, path, locus.records.clone()))
        })
        .collect::<Vec<_>>();
    let mut results = vec![None; catalog.loci.len()];
    if tasks.is_empty() {
        return Ok(results);
    }
    let workers = threads.max(1).min(tasks.len());
    let (sender, receiver) = mpsc::channel();
    let queue = Arc::new(Mutex::new(VecDeque::from(tasks)));
    let mut handles = Vec::with_capacity(workers);
    for _ in 0..workers {
        let sender = sender.clone();
        let queue = Arc::clone(&queue);
        handles.push(thread::spawn(move || loop {
            let task = queue
                .lock()
                .expect("PanRef task queue poisoned")
                .pop_front();
            let Some((id, path, baits)) = task else { break };
            let _ = sender.send((id, assemble_backbone(&path, &baits)));
        }));
    }
    drop(sender);
    for (id, backbone) in receiver {
        results[id] = backbone?;
    }
    for handle in handles {
        handle
            .join()
            .map_err(|_| "PanRef backbone worker panicked".to_string())?;
    }
    Ok(results)
}

fn build_panref_reference(args: &Args, samples: &[Sample]) -> AppResult<PathBuf> {
    let baits = args
        .panref_baits
        .as_ref()
        .ok_or("--engine panref/panrefv2 requires --panref-baits")?;
    let reference_dir = args.output.join("population").join("reference");
    let panref_dir = reference_dir.join("panref");
    fs::create_dir_all(&panref_dir)
        .map_err(|e| format!("unable to create {}: {e}", panref_dir.display()))?;
    let catalog = BaitCatalog::read(baits)?;
    let index = BaitIndex::build_catalog(&catalog)?;
    let loci = catalog
        .loci
        .iter()
        .map(|entry| entry.name.clone())
        .collect::<Vec<_>>();
    index.write_metadata(&panref_dir.join("index_metadata.tsv"), &loci)?;
    write_sample_manifest(
        &args.output.join("population").join("sample_manifest.tsv"),
        samples,
    )?;
    let recruited_dir = panref_dir.join("recruited");
    if recruited_dir.exists() {
        fs::remove_dir_all(&recruited_dir)
            .map_err(|e| format!("unable to clear {}: {e}", recruited_dir.display()))?;
    }
    let report_path = panref_dir.join("recruitment_summary.tsv");
    let mut report = BufWriter::new(File::create(&report_path).map_err(|e| e.to_string())?);
    writeln!(report, "sample\ttotal_locus_assignments\tloci_with_pairs\tstrong_pairs\trescued_pairs\tambiguous_pairs")
        .map_err(|e| e.to_string())?;
    report.flush().map_err(|e| e.to_string())?;
    for sample in samples {
        let stats = recruit_pairs_to_fastq(
            &index,
            &catalog,
            &sample.read1,
            &sample.read2,
            &recruited_dir,
            &sample.internal,
            args.threads,
        )?;
        let total: u64 = stats.per_locus.iter().sum();
        let covered = stats.per_locus.iter().filter(|&count| *count > 0).count();
        writeln!(
            report,
            "{}\t{}\t{}\t{}\t{}\t{}",
            sample.internal,
            total,
            covered,
            stats.strong_pairs,
            stats.rescued_pairs,
            stats.ambiguous_pairs
        )
        .map_err(|e| e.to_string())?;
    }
    report.flush().map_err(|e| e.to_string())?;
    let fasta_path = reference_dir.join("population_reference.fasta");
    let mut fasta = BufWriter::new(File::create(&fasta_path).map_err(|e| e.to_string())?);
    let mut backbone_report = BufWriter::new(
        File::create(panref_dir.join("backbone_summary.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(backbone_report, "locus\tstatus\tsequence_length\tpe_links")
        .map_err(|e| e.to_string())?;
    let backbones = assemble_panref_backbones(&recruited_dir, &catalog, args.threads)?;
    let mut written = 0_usize;
    let mut used_names = HashSet::new();
    for (id, locus) in catalog.loci.iter().enumerate() {
        match backbones[id].as_ref() {
            Some((sequence, pe_links)) => {
                let name = safe_locus_name(&locus.name);
                if !used_names.insert(name.clone()) {
                    return Err(format!(
                        "PanRef locus names collide after FASTA sanitization: {}",
                        locus.name
                    ));
                }
                writeln!(fasta, ">{name}").map_err(|e| e.to_string())?;
                writeln!(fasta, "{}", String::from_utf8_lossy(sequence))
                    .map_err(|e| e.to_string())?;
                writeln!(
                    backbone_report,
                    "{}\tassembled\t{}\t{}",
                    locus.name,
                    sequence.len(),
                    pe_links
                )
                .map_err(|e| e.to_string())?;
                written += 1;
            }
            None => writeln!(backbone_report, "{}\tno_backbone\t0\t0", locus.name)
                .map_err(|e| e.to_string())?,
        }
    }
    backbone_report.flush().map_err(|e| e.to_string())?;
    if written == 0 {
        return Err(
            "PanRef found no bait-anchored local backbone; see population/reference/panref".into(),
        );
    }
    write_reference_manifest(&reference_dir, "panref", &fasta_path)?;
    Ok(fasta_path)
}

fn build_pseudoref_reference(args: &Args, samples: &[Sample]) -> AppResult<PathBuf> {
    let root = args.output.join("population");
    let reference_dir = root.join("reference");
    fs::create_dir_all(&reference_dir)
        .map_err(|e| format!("unable to create {}: {e}", reference_dir.display()))?;
    let mut candidates: BTreeMap<String, Vec<Candidate>> = BTreeMap::new();
    for sample in samples {
        let sample_dir = args.output.join(&sample.internal);
        let summary_path = sample_dir.join("uce_assembly_summary.csv");
        if !summary_path.is_file() {
            return Err(format!(
                "missing UCE summary for {}: {}",
                sample.internal,
                summary_path.display()
            ));
        }
        for row in parse_csv_records(&summary_path)? {
            if !accepted(&row) {
                continue;
            }
            let locus = row.get("locus").cloned().unwrap_or_default();
            if locus.is_empty() {
                continue;
            }
            let fasta = sample_dir.join("results").join(format!("{locus}.fasta"));
            if !fasta.is_file() {
                continue;
            }
            let sequence = read_first_fasta(&fasta)?;
            let reported_length: usize = number(&row, "selected_contig_length");
            if reported_length > 0 && reported_length != sequence.len() {
                return Err(format!(
                    "length mismatch for {locus} in {}",
                    sample.internal
                ));
            }
            let unique_reads: u64 = number(&row, "unique_read_count");
            if unique_reads == 0 {
                continue;
            }
            candidates
                .entry(locus.clone())
                .or_default()
                .push(Candidate {
                    safe_locus: safe_locus_name(&locus),
                    locus,
                    sample: sample.internal.clone(),
                    supported_bases: number(&row, "slice_supported_bases"),
                    support_breadth: number(&row, "slice_support_breadth"),
                    max_gap: number(&row, "max_slice_support_gap"),
                    unique_reads,
                    unique_density: number(&row, "unique_read_density"),
                    length: sequence.len(),
                    sequence,
                });
        }
    }
    if candidates.is_empty() {
        return Err("no accepted UCE contigs were found".into());
    }

    let fasta_path = reference_dir.join("population_reference.fasta");
    let provenance_path = reference_dir.join("population_reference_provenance.tsv");
    let name_map_path = reference_dir.join("locus_name_map.tsv");
    let contribution_path = reference_dir.join("reference_contribution.tsv");
    let mut fasta = BufWriter::new(File::create(&fasta_path).map_err(|e| e.to_string())?);
    let mut provenance = BufWriter::new(File::create(&provenance_path).map_err(|e| e.to_string())?);
    let mut name_map = BufWriter::new(File::create(&name_map_path).map_err(|e| e.to_string())?);
    writeln!(provenance, "locus\treference_locus\treference_sample\tselection_strategy\tcandidate_count\tsequence_length\tsupported_bases\tsupport_breadth\tmax_support_gap\tunique_read_count\tunique_read_density").map_err(|e| e.to_string())?;
    writeln!(name_map, "locus\treference_locus").map_err(|e| e.to_string())?;
    let mut used_names = HashSet::new();
    let mut contribution_counts: BTreeMap<String, usize> = BTreeMap::new();
    for (locus, locus_candidates) in &candidates {
        let chosen = locus_candidates
            .iter()
            .max_by(|a, b| candidate_cmp(args.reference_strategy, a, b))
            .unwrap();
        if !used_names.insert(chosen.safe_locus.clone()) {
            return Err(format!(
                "locus names collide after VCF sanitization: {locus}"
            ));
        }
        writeln!(fasta, ">{}", chosen.safe_locus).map_err(|e| e.to_string())?;
        for chunk in chosen.sequence.as_bytes().chunks(80) {
            writeln!(fasta, "{}", String::from_utf8_lossy(chunk)).map_err(|e| e.to_string())?;
        }
        writeln!(
            provenance,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{:.6}",
            chosen.locus,
            chosen.safe_locus,
            chosen.sample,
            args.reference_strategy.label(),
            locus_candidates.len(),
            chosen.length,
            chosen.supported_bases,
            chosen.support_breadth,
            chosen.max_gap,
            chosen.unique_reads,
            chosen.unique_density
        )
        .map_err(|e| e.to_string())?;
        writeln!(name_map, "{locus}\t{}", chosen.safe_locus).map_err(|e| e.to_string())?;
        *contribution_counts
            .entry(chosen.sample.clone())
            .or_default() += 1;
    }
    let total_loci: usize = contribution_counts.values().sum();
    let mut contribution =
        BufWriter::new(File::create(&contribution_path).map_err(|e| e.to_string())?);
    writeln!(
        contribution,
        "internal_sample_id\tvcf_sample_id\treference_loci\treference_fraction"
    )
    .map_err(|e| e.to_string())?;
    for sample in samples {
        let count = contribution_counts
            .get(&sample.internal)
            .copied()
            .unwrap_or_default();
        let fraction = if total_loci == 0 {
            0.0
        } else {
            count as f64 / total_loci as f64
        };
        writeln!(
            contribution,
            "{}\t{}\t{}\t{:.6}",
            sample.internal, sample.vcf_id, count, fraction
        )
        .map_err(|e| e.to_string())?;
    }
    fasta.flush().map_err(|e| e.to_string())?;
    provenance.flush().map_err(|e| e.to_string())?;
    name_map.flush().map_err(|e| e.to_string())?;
    contribution.flush().map_err(|e| e.to_string())?;
    write_sample_manifest(&root.join("sample_manifest.tsv"), samples)?;
    write_reference_manifest(&reference_dir, "pseudoref", &fasta_path)?;
    Ok(fasta_path)
}

fn write_reference_manifest(reference_dir: &Path, engine: &str, fasta: &Path) -> AppResult<()> {
    let bytes = fs::metadata(fasta).map_err(|e| e.to_string())?.len();
    fs::write(
        reference_dir.join("reference_manifest.tsv"),
        format!(
            "engine\treference_fasta\tbytes\n{}\t{}\t{}\n",
            engine,
            fasta
                .file_name()
                .and_then(|value| value.to_str())
                .unwrap_or("population_reference.fasta"),
            bytes
        ),
    )
    .map_err(|e| e.to_string())
}

fn materialize_external_reference(args: &Args, samples: &[Sample]) -> AppResult<PathBuf> {
    let source = args
        .reference_fasta
        .as_ref()
        .ok_or_else(|| "missing --reference-fasta".to_string())?;
    if !source.is_file() {
        return Err(format!(
            "external reference is not a readable file: {}",
            source.display()
        ));
    }
    let reference_dir = args.output.join("population").join("reference");
    fs::create_dir_all(&reference_dir)
        .map_err(|e| format!("unable to create {}: {e}", reference_dir.display()))?;
    let target = reference_dir.join("population_reference.fasta");
    if source != &target {
        fs::copy(source, &target).map_err(|e| {
            format!(
                "unable to copy external reference {} to {}: {e}",
                source.display(),
                target.display()
            )
        })?;
    }
    let _ = read_first_fasta(&target)?;
    let source_path = source.canonicalize().unwrap_or_else(|_| source.clone());
    fs::write(
        reference_dir.join("reference_source.tsv"),
        format!(
            "reference_mode\tsource_path\tmaterialized_path\nexternal\t{}\t{}\n",
            source_path.display(),
            target.display()
        ),
    )
    .map_err(|e| e.to_string())?;
    write_sample_manifest(
        &args.output.join("population").join("sample_manifest.tsv"),
        samples,
    )?;
    write_reference_manifest(&reference_dir, "external", &target)?;
    Ok(target)
}

fn existing_population_reference(args: &Args) -> AppResult<PathBuf> {
    let reference = args
        .output
        .join("population")
        .join("reference")
        .join("population_reference.fasta");
    if !reference.is_file() {
        return Err(format!(
            "missing existing population reference: {}; rerun with --population-start-at reference",
            reference.display()
        ));
    }
    Ok(reference)
}

fn prepare_reference(args: &Args, samples: &[Sample]) -> AppResult<PathBuf> {
    if args.reference_fasta.is_some() {
        materialize_external_reference(args, samples)
    } else if args.start_at == Stage::Reference {
        build_reference(args, samples)
    } else {
        existing_population_reference(args)
    }
}

fn existing_bams(args: &Args, sample_count: usize) -> AppResult<Vec<PathBuf>> {
    let list_path = args
        .output
        .join("population")
        .join("mapping")
        .join("bam.list");
    let list = File::open(&list_path).map_err(|_| {
        format!(
            "missing existing BAM list: {}; rerun with --population-start-at mapping",
            list_path.display()
        )
    })?;
    let bams: Vec<PathBuf> = BufReader::new(list)
        .lines()
        .map(|line| line.map(PathBuf::from).map_err(|e| e.to_string()))
        .collect::<AppResult<Vec<_>>>()?
        .into_iter()
        .filter(|path| !path.as_os_str().is_empty())
        .collect();
    if bams.len() != sample_count {
        return Err(format!(
            "existing BAM list has {} entries for {sample_count} samples; rerun with --population-start-at mapping",
            bams.len()
        ));
    }
    for bam in &bams {
        let index = PathBuf::from(format!("{}.bai", bam.display()));
        if !bam.is_file() || !index.is_file() {
            return Err(format!(
                "missing existing BAM or index for {}; rerun with --population-start-at mapping",
                bam.display()
            ));
        }
    }
    Ok(bams)
}

fn existing_filtered_vcf(args: &Args) -> AppResult<PathBuf> {
    let filtered = args
        .output
        .join("population")
        .join("variants")
        .join("cohort.filtered.vcf.gz");
    let csi = PathBuf::from(format!("{}.csi", filtered.display()));
    let tbi = PathBuf::from(format!("{}.tbi", filtered.display()));
    if !filtered.is_file() || (!csi.is_file() && !tbi.is_file()) {
        return Err(format!(
            "missing indexed filtered VCF: {}; rerun with --population-start-at calling",
            filtered.display()
        ));
    }
    Ok(filtered)
}

fn write_sample_manifest(path: &Path, samples: &[Sample]) -> AppResult<()> {
    let mut out = BufWriter::new(
        File::create(path).map_err(|e| format!("unable to write {}: {e}", path.display()))?,
    );
    writeln!(
        out,
        "original_sample_id\tinternal_sample_id\tvcf_sample_id\tread1\tread2\tlayout\tpopulation\tbatch"
    )
    .map_err(|e| e.to_string())?;
    for sample in samples {
        let layout = if sample.read1 == sample.read2 {
            "SE"
        } else {
            "PE"
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            sample.original,
            sample.internal,
            sample.vcf_id,
            sample.read1.display(),
            sample.read2.display(),
            layout,
            sample.population,
            sample.batch
        )
        .map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn run_command(program: &str, argv: &[String]) -> AppResult<()> {
    let status = Command::new(program)
        .args(argv)
        .status()
        .map_err(|e| format!("unable to run {program}: {e}"))?;
    if !status.success() {
        return Err(format!("{program} failed with status {status}"));
    }
    Ok(())
}

fn run_pipe(
    left_program: &str,
    left_args: &[String],
    right_program: &str,
    right_args: &[String],
) -> AppResult<()> {
    let mut left = Command::new(left_program)
        .args(left_args)
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|e| format!("unable to run {left_program}: {e}"))?;
    let stdout = left
        .stdout
        .take()
        .ok_or_else(|| format!("unable to capture {left_program} output"))?;
    let right_status = Command::new(right_program)
        .args(right_args)
        .stdin(Stdio::from(stdout))
        .status()
        .map_err(|e| format!("unable to run {right_program}: {e}"))?;
    let left_status = left
        .wait()
        .map_err(|e| format!("unable to wait for {left_program}: {e}"))?;
    if !left_status.success() {
        return Err(format!("{left_program} failed with status {left_status}"));
    }
    if !right_status.success() {
        return Err(format!("{right_program} failed with status {right_status}"));
    }
    Ok(())
}

fn minibwa_map_args(threads: usize, sample: &Sample, reference: &Path) -> Vec<String> {
    let rg = format!(
        "@RG\\tID:{}\\tSM:{}\\tPL:ILLUMINA",
        sample.vcf_id, sample.vcf_id
    );
    let mut mapper_args = vec![
        "map".into(),
        "-t".into(),
        threads.to_string(),
        "-R".into(),
        rg,
        reference.display().to_string(),
        sample.read1.display().to_string(),
    ];
    if sample.read1 != sample.read2 {
        mapper_args.push(sample.read2.display().to_string());
    }
    mapper_args
}

#[derive(Debug, Default, PartialEq)]
struct MappingMetrics {
    total_reads: u64,
    mapped_reads: u64,
    properly_paired_reads: u64,
    reference_bases: u64,
    covered_bases: u64,
    mean_depth: f64,
}

fn parse_flagstat(text: &str, metrics: &mut MappingMetrics) {
    for line in text.lines() {
        let count = line
            .split_whitespace()
            .next()
            .and_then(|value| value.parse::<u64>().ok())
            .unwrap_or_default();
        if line.contains(" in total ") || line.ends_with(" in total") {
            metrics.total_reads = count;
        } else if line.contains(" mapped (") {
            metrics.mapped_reads = count;
        } else if line.contains(" properly paired (") {
            metrics.properly_paired_reads = count;
        }
    }
}

fn parse_coverage(text: &str, metrics: &mut MappingMetrics) {
    let mut depth_sum = 0.0;
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 7 {
            continue;
        }
        let start = fields[1].parse::<u64>().unwrap_or(1);
        let end = fields[2].parse::<u64>().unwrap_or_default();
        let length = end.saturating_sub(start).saturating_add(1);
        let covered = fields[4].parse::<u64>().unwrap_or_default();
        let mean_depth = fields[6].parse::<f64>().unwrap_or_default();
        metrics.reference_bases += length;
        metrics.covered_bases += covered.min(length);
        depth_sum += mean_depth * length as f64;
    }
    if metrics.reference_bases > 0 {
        metrics.mean_depth = depth_sum / metrics.reference_bases as f64;
    }
}

fn rate(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

fn map_samples(args: &Args, samples: &[Sample], reference: &Path) -> AppResult<Vec<PathBuf>> {
    let mapping_dir = args.output.join("population").join("mapping");
    fs::create_dir_all(&mapping_dir)
        .map_err(|e| format!("unable to create {}: {e}", mapping_dir.display()))?;
    run_command(
        &args.minibwa,
        &[
            "index".into(),
            "-t".into(),
            args.threads.to_string(),
            reference.display().to_string(),
        ],
    )?;
    run_command(
        &args.samtools,
        &["faidx".into(), reference.display().to_string()],
    )?;

    // Alignment and sort/markdup show diminishing returns well before the
    // full core count on typical UCE-sized references, so several samples
    // are mapped concurrently at a reduced per-sample thread count instead
    // of giving every thread to one sample at a time. The reduce below
    // writes `qc`/`bams` in the original sample order afterward, so output
    // is byte-identical to the prior strictly serial version.
    let worker_count = args.threads.max(1).min(samples.len()).max(1);
    let per_sample_threads = (args.threads / worker_count).max(1);
    let queue: Mutex<VecDeque<(usize, &Sample)>> = Mutex::new(samples.iter().enumerate().collect());
    let (sender, receiver) = mpsc::channel();
    thread::scope(|scope| {
        for _ in 0..worker_count {
            let queue = &queue;
            let sender = sender.clone();
            let mapping_dir = &mapping_dir;
            scope.spawn(move || loop {
                let next = queue.lock().expect("mapping queue poisoned").pop_front();
                let Some((index, sample)) = next else {
                    break;
                };
                let outcome =
                    map_one_sample(args, sample, reference, mapping_dir, per_sample_threads);
                if sender.send((index, outcome)).is_err() {
                    break;
                }
            });
        }
    });
    drop(sender);
    let mut ordered: Vec<Option<AppResult<(PathBuf, String)>>> =
        (0..samples.len()).map(|_| None).collect();
    for (index, outcome) in receiver {
        ordered[index] = Some(outcome);
    }

    let mut bams = Vec::with_capacity(samples.len());
    let mut qc = BufWriter::new(
        File::create(mapping_dir.join("mapping_qc.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(qc, "sample\tstatus\ttotal_reads\tmapped_reads\tmapping_rate\tproperly_paired_reads\tproperly_paired_rate\treference_bases\tcovered_bases\tcoverage_breadth\tmean_depth\tbam").map_err(|e| e.to_string())?;
    for slot in ordered {
        let (bam, qc_line) = slot.expect("every queued sample reports exactly once")?;
        writeln!(qc, "{qc_line}").map_err(|e| e.to_string())?;
        bams.push(bam);
    }
    let mut list =
        BufWriter::new(File::create(mapping_dir.join("bam.list")).map_err(|e| e.to_string())?);
    for bam in &bams {
        writeln!(list, "{}", bam.display()).map_err(|e| e.to_string())?;
    }
    Ok(bams)
}

fn map_one_sample(
    args: &Args,
    sample: &Sample,
    reference: &Path,
    mapping_dir: &Path,
    threads: usize,
) -> AppResult<(PathBuf, String)> {
    let bam = mapping_dir.join(format!("{}.bam", sample.vcf_id));
    let mapper_args = minibwa_map_args(threads, sample, reference);
    if args.mark_duplicates {
        let name_bam = mapping_dir.join(format!("{}.name.bam", sample.vcf_id));
        let fixmate_bam = mapping_dir.join(format!("{}.fixmate.bam", sample.vcf_id));
        let position_bam = mapping_dir.join(format!("{}.positions.bam", sample.vcf_id));
        run_pipe(
            &args.minibwa,
            &mapper_args,
            &args.samtools,
            &[
                "sort".into(),
                "-n".into(),
                "-@".into(),
                threads.to_string(),
                "-o".into(),
                name_bam.display().to_string(),
                "-".into(),
            ],
        )?;
        run_command(
            &args.samtools,
            &[
                "fixmate".into(),
                "-m".into(),
                name_bam.display().to_string(),
                fixmate_bam.display().to_string(),
            ],
        )?;
        run_command(
            &args.samtools,
            &[
                "sort".into(),
                "-@".into(),
                threads.to_string(),
                "-o".into(),
                position_bam.display().to_string(),
                fixmate_bam.display().to_string(),
            ],
        )?;
        run_command(
            &args.samtools,
            &[
                "markdup".into(),
                position_bam.display().to_string(),
                bam.display().to_string(),
            ],
        )?;
        for temporary in [&name_bam, &fixmate_bam, &position_bam] {
            let _ = fs::remove_file(temporary);
        }
    } else {
        run_pipe(
            &args.minibwa,
            &mapper_args,
            &args.samtools,
            &[
                "sort".into(),
                "-@".into(),
                threads.to_string(),
                "-o".into(),
                bam.display().to_string(),
                "-".into(),
            ],
        )?;
    }
    run_command(
        &args.samtools,
        &["quickcheck".into(), bam.display().to_string()],
    )?;
    run_command(
        &args.samtools,
        &[
            "index".into(),
            "-@".into(),
            threads.to_string(),
            bam.display().to_string(),
        ],
    )?;
    let flagstat_path = mapping_dir.join(format!("{}.flagstat.txt", sample.vcf_id));
    let output = Command::new(&args.samtools)
        .args(["flagstat", &bam.display().to_string()])
        .output()
        .map_err(|e| format!("unable to run samtools flagstat: {e}"))?;
    if !output.status.success() {
        return Err(format!("samtools flagstat failed for {}", sample.vcf_id));
    }
    fs::write(&flagstat_path, &output.stdout).map_err(|e| e.to_string())?;
    let coverage_output = Command::new(&args.samtools)
        .args(["coverage", &bam.display().to_string()])
        .output()
        .map_err(|e| format!("unable to run samtools coverage: {e}"))?;
    if !coverage_output.status.success() {
        return Err(format!("samtools coverage failed for {}", sample.vcf_id));
    }
    let mut metrics = MappingMetrics::default();
    parse_flagstat(&String::from_utf8_lossy(&output.stdout), &mut metrics);
    parse_coverage(
        &String::from_utf8_lossy(&coverage_output.stdout),
        &mut metrics,
    );
    let qc_line = format!(
        "{}\tok\t{}\t{}\t{:.6}\t{}\t{:.6}\t{}\t{}\t{:.6}\t{:.4}\t{}",
        sample.vcf_id,
        metrics.total_reads,
        metrics.mapped_reads,
        rate(metrics.mapped_reads, metrics.total_reads),
        metrics.properly_paired_reads,
        rate(metrics.properly_paired_reads, metrics.total_reads),
        metrics.reference_bases,
        metrics.covered_bases,
        rate(metrics.covered_bases, metrics.reference_bases),
        metrics.mean_depth,
        bam.display()
    );
    Ok((bam, qc_line))
}

fn call_variants(args: &Args, reference: &Path, bams: &[PathBuf]) -> AppResult<PathBuf> {
    if bams.is_empty() {
        return Err("no BAM files for joint calling".into());
    }
    let variants = args.output.join("population").join("variants");
    fs::create_dir_all(&variants)
        .map_err(|e| format!("unable to create {}: {e}", variants.display()))?;
    let bam_list = args
        .output
        .join("population")
        .join("mapping")
        .join("bam.list");
    let raw = variants.join("cohort.raw.bcf");
    run_pipe(
        &args.bcftools,
        &[
            "mpileup".into(),
            "-f".into(),
            reference.display().to_string(),
            "-b".into(),
            bam_list.display().to_string(),
            "-q".into(),
            args.min_mapq.to_string(),
            "-Q".into(),
            args.min_baseq.to_string(),
            "-a".into(),
            "FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR".into(),
            "-Ou".into(),
        ],
        &args.bcftools,
        &[
            "call".into(),
            "-m".into(),
            "-v".into(),
            "-G".into(),
            "-".into(),
            "-f".into(),
            "GQ".into(),
            "-Ob".into(),
            "-o".into(),
            raw.display().to_string(),
        ],
    )?;
    run_command(
        &args.bcftools,
        &["index".into(), "-f".into(), raw.display().to_string()],
    )?;

    let biallelic = variants.join("cohort.biallelic.snps.vcf.gz");
    run_pipe(
        &args.bcftools,
        &[
            "norm".into(),
            "-f".into(),
            reference.display().to_string(),
            "-m".into(),
            "-any".into(),
            "-Ou".into(),
            raw.display().to_string(),
        ],
        &args.bcftools,
        &[
            "view".into(),
            "-v".into(),
            "snps".into(),
            "-m2".into(),
            "-M2".into(),
            "-Oz".into(),
            "-o".into(),
            biallelic.display().to_string(),
        ],
    )?;
    run_command(
        &args.bcftools,
        &["index".into(), "-f".into(), biallelic.display().to_string()],
    )?;

    let genotype_filtered = variants.join("cohort.genotype_filtered.vcf.gz");
    let gt_expression = genotype_filter_expression(args);
    run_command(
        &args.bcftools,
        &[
            "+setGT".into(),
            biallelic.display().to_string(),
            "-Oz".into(),
            "-o".into(),
            genotype_filtered.display().to_string(),
            "--".into(),
            "-t".into(),
            "q".into(),
            "-n".into(),
            ".".into(),
            "-i".into(),
            gt_expression,
        ],
    )?;
    run_command(
        &args.bcftools,
        &[
            "index".into(),
            "-f".into(),
            genotype_filtered.display().to_string(),
        ],
    )?;

    let tagged = variants.join("cohort.tagged.vcf.gz");
    run_command(
        &args.bcftools,
        &[
            "+fill-tags".into(),
            genotype_filtered.display().to_string(),
            "-Oz".into(),
            "-o".into(),
            tagged.display().to_string(),
            "--".into(),
            "-t".into(),
            "F_MISSING,AC,AN,MAF".into(),
        ],
    )?;
    run_command(
        &args.bcftools,
        &["index".into(), "-f".into(), tagged.display().to_string()],
    )?;

    let filtered = variants.join("cohort.filtered.vcf.gz");
    let expression = format!(
        "QUAL>={} && F_MISSING<={} && AC>={} && (AN-AC)>={}",
        args.min_qual,
        1.0 - args.min_call_rate,
        args.min_mac,
        args.min_mac
    );
    run_command(
        &args.bcftools,
        &[
            "view".into(),
            "-i".into(),
            expression,
            "-Oz".into(),
            "-o".into(),
            filtered.display().to_string(),
            tagged.display().to_string(),
        ],
    )?;
    run_command(
        &args.bcftools,
        &["index".into(), "-f".into(), filtered.display().to_string()],
    )?;
    write_variant_qc(
        args,
        &[
            ("raw_called_variants", &raw),
            ("biallelic_snps", &biallelic),
            ("genotype_filtered", &genotype_filtered),
            ("tagged", &tagged),
            ("site_filtered", &filtered),
        ],
    )?;
    Ok(filtered)
}

fn indexed_variant_count(args: &Args, path: &Path) -> AppResult<u64> {
    let output = Command::new(&args.bcftools)
        .args(["index", "-n", &path.display().to_string()])
        .output()
        .map_err(|e| format!("unable to run bcftools index -n: {e}"))?;
    if !output.status.success() {
        return Err(format!(
            "bcftools index -n failed for {} with status {}",
            path.display(),
            output.status
        ));
    }
    String::from_utf8_lossy(&output.stdout)
        .trim()
        .parse::<u64>()
        .map_err(|_| format!("invalid variant count reported for {}", path.display()))
}

fn write_variant_qc(args: &Args, stages: &[(&str, &Path)]) -> AppResult<()> {
    let path = args
        .output
        .join("population")
        .join("variants")
        .join("variant_qc.tsv");
    let mut out = BufWriter::new(File::create(&path).map_err(|e| e.to_string())?);
    writeln!(out, "stage\tvariant_records\tpath").map_err(|e| e.to_string())?;
    for (stage, variant_path) in stages {
        writeln!(
            out,
            "{}\t{}\t{}",
            stage,
            indexed_variant_count(args, variant_path)?,
            variant_path.display()
        )
        .map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn genotype_filter_expression(args: &Args) -> String {
    format!("FMT/DP<{} | FMT/GQ<{}", args.min_dp, args.min_gq)
}

fn median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let middle = values.len() / 2;
    if values.len().is_multiple_of(2) {
        (values[middle - 1] + values[middle]) / 2.0
    } else {
        values[middle]
    }
}

fn parse_site(line: &str) -> AppResult<Site> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        return Err("malformed VCF record".into());
    }
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let gt_index = format_keys.iter().position(|v| *v == "GT");
    let dp_index = format_keys.iter().position(|v| *v == "DP");
    let gq_index = format_keys.iter().position(|v| *v == "GQ");
    let mut called = 0_usize;
    let total = fields.len().saturating_sub(9);
    let mut depths = Vec::new();
    let mut qualities = Vec::new();
    for sample in &fields[9..] {
        let values: Vec<&str> = sample.split(':').collect();
        let gt = gt_index.and_then(|i| values.get(i)).copied().unwrap_or(".");
        if genotype_is_called(gt) {
            called += 1;
            if let Some(value) = dp_index
                .and_then(|i| values.get(i))
                .and_then(|v| v.parse::<f64>().ok())
            {
                depths.push(value);
            }
            if let Some(value) = gq_index
                .and_then(|i| values.get(i))
                .and_then(|v| v.parse::<f64>().ok())
            {
                qualities.push(value);
            }
        }
    }
    Ok(Site {
        chrom: fields[0].into(),
        pos: fields[1]
            .parse()
            .map_err(|_| "invalid VCF position".to_string())?,
        reference: fields[3].into(),
        alternate: fields[4].into(),
        qual: fields[5].parse().unwrap_or(0.0),
        call_rate: if total == 0 {
            0.0
        } else {
            called as f64 / total as f64
        },
        median_dp: median(&mut depths),
        median_gq: median(&mut qualities),
        record: line.to_string(),
    })
}

fn genotype_is_called(genotype: &str) -> bool {
    !genotype.is_empty()
        && genotype
            .split(['/', '|'])
            .all(|allele| !allele.is_empty() && allele != ".")
}

fn site_cmp(left: &Site, right: &Site, locus_depth: f64) -> Ordering {
    left.call_rate
        .partial_cmp(&right.call_rate)
        .unwrap_or(Ordering::Equal)
        .then_with(|| {
            left.median_gq
                .partial_cmp(&right.median_gq)
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| {
            let left_distance = (left.median_dp - locus_depth).abs();
            let right_distance = (right.median_dp - locus_depth).abs();
            right_distance
                .partial_cmp(&left_distance)
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| {
            left.qual
                .partial_cmp(&right.qual)
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| right.pos.cmp(&left.pos))
}

fn stream_vcf_sites<F>(
    args: &Args,
    filtered: &Path,
    capture_headers: bool,
    mut visit: F,
) -> AppResult<Vec<String>>
where
    F: FnMut(Site) -> AppResult<()>,
{
    let mut child = Command::new(&args.bcftools)
        .args(["view", "-Ov", &filtered.display().to_string()])
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|e| format!("unable to run bcftools view: {e}"))?;
    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| "unable to capture bcftools view output".to_string())?;
    let mut headers = Vec::new();
    let processing_result = (|| -> AppResult<()> {
        for line in BufReader::new(stdout).lines() {
            let line = line.map_err(|e| format!("unable to read VCF stream: {e}"))?;
            if line.starts_with('#') {
                if capture_headers {
                    headers.push(line);
                }
            } else if !line.is_empty() {
                visit(parse_site(&line)?)?;
            }
        }
        Ok(())
    })();

    if processing_result.is_err() {
        let _ = child.kill();
    }
    let status = child
        .wait()
        .map_err(|e| format!("unable to wait for bcftools view: {e}"))?;
    processing_result?;
    if !status.success() {
        return Err(format!("bcftools view failed with status {status}"));
    }
    Ok(headers)
}

fn select_one_snp(args: &Args, filtered: &Path) -> AppResult<PathBuf> {
    let mut locus_depths: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    stream_vcf_sites(args, filtered, false, |site| {
        locus_depths
            .entry(site.chrom)
            .or_default()
            .push(site.median_dp);
        Ok(())
    })?;
    if locus_depths.is_empty() {
        return Err("no SNP passed population filters".into());
    }

    let mut locus_stats: BTreeMap<String, (f64, usize)> = BTreeMap::new();
    for (locus, mut depths) in locus_depths {
        let candidate_count = depths.len();
        locus_stats.insert(locus, (median(&mut depths), candidate_count));
    }

    let mut selected: BTreeMap<String, Site> = BTreeMap::new();
    let mut locus_order = Vec::new();
    let mut seen_loci = HashSet::new();
    let headers = stream_vcf_sites(args, filtered, true, |site| {
        let (locus_depth, _) = locus_stats
            .get(&site.chrom)
            .ok_or_else(|| format!("missing first-pass statistics for {}", site.chrom))?;
        if seen_loci.insert(site.chrom.clone()) {
            locus_order.push(site.chrom.clone());
        }
        let should_replace = selected
            .get(&site.chrom)
            .map(|current| site_cmp(&site, current, *locus_depth) == Ordering::Greater)
            .unwrap_or(true);
        if should_replace {
            selected.insert(site.chrom.clone(), site);
        }
        Ok(())
    })?;
    if selected.len() != locus_stats.len() {
        return Err(format!(
            "selected {} SNPs for {} loci; expected exactly one per locus",
            selected.len(),
            locus_stats.len()
        ));
    }

    let structure = args.output.join("population").join("structure");
    fs::create_dir_all(&structure)
        .map_err(|e| format!("unable to create {}: {e}", structure.display()))?;
    let table_path = structure.join("selected_snps.tsv");
    let mut table = BufWriter::new(File::create(&table_path).map_err(|e| e.to_string())?);
    writeln!(
        table,
        "locus\tposition\tref\talt\tqual\tcall_rate\tmedian_dp\tmedian_gq\tcandidate_snp_count"
    )
    .map_err(|e| e.to_string())?;
    let raw_vcf_path = structure.join("selected_snps.raw.vcf");
    let mut raw_vcf = BufWriter::new(File::create(&raw_vcf_path).map_err(|e| e.to_string())?);
    for header in &headers {
        writeln!(raw_vcf, "{header}").map_err(|e| e.to_string())?;
    }
    for locus in &locus_order {
        let chosen = selected
            .get(locus)
            .ok_or_else(|| format!("missing selected SNP for {locus}"))?;
        let (_, candidate_count) = locus_stats
            .get(locus)
            .ok_or_else(|| format!("missing SNP count for {locus}"))?;
        writeln!(
            table,
            "{}\t{}\t{}\t{}\t{:.2}\t{:.6}\t{:.2}\t{:.2}\t{}",
            locus,
            chosen.pos,
            chosen.reference,
            chosen.alternate,
            chosen.qual,
            chosen.call_rate,
            chosen.median_dp,
            chosen.median_gq,
            candidate_count
        )
        .map_err(|e| e.to_string())?;
        writeln!(raw_vcf, "{}", chosen.record).map_err(|e| e.to_string())?;
    }
    drop(raw_vcf);
    let final_vcf = structure.join("one_snp_per_uce.vcf.gz");
    run_command(
        &args.bcftools,
        &[
            "annotate".into(),
            "--set-id".into(),
            "%CHROM:%POS:%REF:%FIRST_ALT".into(),
            "-Oz".into(),
            "-o".into(),
            final_vcf.display().to_string(),
            raw_vcf_path.display().to_string(),
        ],
    )?;
    run_command(
        &args.bcftools,
        &["index".into(), "-f".into(), final_vcf.display().to_string()],
    )?;
    Ok(final_vcf)
}

fn annotate_panel(args: &Args, source: &Path, target: &Path) -> AppResult<()> {
    run_command(
        &args.bcftools,
        &[
            "annotate".into(),
            "--set-id".into(),
            "%CHROM:%POS:%REF:%FIRST_ALT".into(),
            "-Oz".into(),
            "-o".into(),
            target.display().to_string(),
            source.display().to_string(),
        ],
    )?;
    run_command(
        &args.bcftools,
        &["index".into(), "-f".into(), target.display().to_string()],
    )
}

fn prepare_all_snp_panel(args: &Args, filtered: &Path) -> AppResult<PathBuf> {
    let structure = args.output.join("population").join("structure");
    fs::create_dir_all(&structure)
        .map_err(|e| format!("unable to create {}: {e}", structure.display()))?;
    let target = structure.join("all_snps.vcf.gz");
    annotate_panel(args, filtered, &target)?;
    Ok(target)
}

fn make_plink_panel(args: &Args, vcf: &Path, prefix: &Path) -> AppResult<()> {
    run_command(
        &args.plink,
        &[
            "--vcf".into(),
            vcf.display().to_string(),
            "--double-id".into(),
            "--allow-extra-chr".into(),
            "--make-bed".into(),
            "--out".into(),
            prefix.display().to_string(),
        ],
    )
}

fn run_pca(args: &Args, prefix: &Path, output: &Path) -> AppResult<()> {
    run_command(
        &args.plink,
        &[
            "--bfile".into(),
            prefix.display().to_string(),
            "--allow-extra-chr".into(),
            "--pca".into(),
            "10".into(),
            "--out".into(),
            output.display().to_string(),
        ],
    )
}

fn site_id(site: &Site) -> String {
    format!(
        "{}:{}:{}:{}",
        site.chrom, site.pos, site.reference, site.alternate
    )
}

fn write_vcf_subset(
    args: &Args,
    source: &Path,
    ids_path: &Path,
    target: &Path,
) -> AppResult<usize> {
    let ids_file =
        File::open(ids_path).map_err(|e| format!("unable to read {}: {e}", ids_path.display()))?;
    let mut wanted = HashSet::new();
    for line in BufReader::new(ids_file).lines() {
        let id = line.map_err(|e| format!("unable to read {}: {e}", ids_path.display()))?;
        if !id.trim().is_empty() {
            wanted.insert(id.trim().to_string());
        }
    }
    if wanted.is_empty() {
        return Err("PLINK LD pruning retained no SNPs".into());
    }
    let mut records = Vec::new();
    let headers = stream_vcf_sites(args, source, true, |site| {
        if wanted.contains(&site_id(&site)) {
            records.push(site.record);
        }
        Ok(())
    })?;
    if records.len() != wanted.len() {
        return Err(format!(
            "LD prune list contains {} IDs but {} exact VCF records were found",
            wanted.len(),
            records.len()
        ));
    }
    let raw = target.with_extension("raw.vcf");
    let mut writer = BufWriter::new(
        File::create(&raw).map_err(|e| format!("unable to create {}: {e}", raw.display()))?,
    );
    for header in headers {
        writeln!(writer, "{header}").map_err(|e| e.to_string())?;
    }
    for record in records {
        writeln!(writer, "{record}").map_err(|e| e.to_string())?;
    }
    drop(writer);
    annotate_panel(args, &raw, target)?;
    let _ = fs::remove_file(raw);
    Ok(wanted.len())
}

fn count_bim(path: &Path) -> AppResult<usize> {
    let file = File::open(path).map_err(|e| format!("unable to read {}: {e}", path.display()))?;
    let mut count = 0;
    for line in BufReader::new(file).lines() {
        line.map_err(|e| format!("unable to read {}: {e}", path.display()))?;
        count += 1;
    }
    Ok(count)
}

struct PanelOutputs {
    primary_bed: PathBuf,
}

fn build_structure_panels(args: &Args, all_vcf: &Path, one_vcf: &Path) -> AppResult<PanelOutputs> {
    let structure = args.output.join("population").join("structure");
    let all_prefix = structure.join("all_snps");
    make_plink_panel(args, all_vcf, &all_prefix)?;
    run_pca(args, &all_prefix, &structure.join("all_snps_pca"))?;

    // Keep the historical population.* and population_pca.* names for compatibility.
    let primary_prefix = structure.join("population");
    make_plink_panel(args, one_vcf, &primary_prefix)?;
    run_pca(args, &primary_prefix, &structure.join("population_pca"))?;

    let prune_selection = structure.join("ld_pruned_selection");
    run_command(
        &args.plink,
        &[
            "--bfile".into(),
            all_prefix.display().to_string(),
            "--allow-extra-chr".into(),
            "--indep-pairwise".into(),
            args.ld_window.to_string(),
            args.ld_step.to_string(),
            args.ld_r2.to_string(),
            "--out".into(),
            prune_selection.display().to_string(),
        ],
    )?;
    let prune_in = prune_selection.with_extension("prune.in");
    let ld_vcf = structure.join("ld_pruned.vcf.gz");
    let ld_count = write_vcf_subset(args, all_vcf, &prune_in, &ld_vcf)?;
    let ld_prefix = structure.join("ld_pruned");
    make_plink_panel(args, &ld_vcf, &ld_prefix)?;
    run_pca(args, &ld_prefix, &structure.join("ld_pruned_pca"))?;

    let mut summary = BufWriter::new(
        File::create(structure.join("panel_summary.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(summary, "panel\tvariants\tvcf\tplink_prefix\tprimary_use")
        .map_err(|e| e.to_string())?;
    writeln!(
        summary,
        "all_snps\t{}\t{}\t{}\tLi-style sensitivity, diversity and FST",
        count_bim(&all_prefix.with_extension("bim"))?,
        all_vcf.display(),
        all_prefix.display()
    )
    .map_err(|e| e.to_string())?;
    writeln!(
        summary,
        "one_snp_per_uce\t{}\t{}\t{}\tprimary PCA and ADMIXTURE",
        count_bim(&primary_prefix.with_extension("bim"))?,
        one_vcf.display(),
        primary_prefix.display()
    )
    .map_err(|e| e.to_string())?;
    writeln!(
        summary,
        "ld_pruned\t{}\t{}\t{}\tLD-pruned sensitivity analysis",
        ld_count,
        ld_vcf.display(),
        ld_prefix.display()
    )
    .map_err(|e| e.to_string())?;
    Ok(PanelOutputs {
        primary_bed: primary_prefix.with_extension("bed"),
    })
}

fn executable_available(program: &str) -> bool {
    let path = Path::new(program);
    if path.components().count() > 1 {
        return path.is_file();
    }
    env::var_os("PATH")
        .map(|paths| env::split_paths(&paths).any(|directory| directory.join(program).is_file()))
        .unwrap_or(false)
}

fn parse_cv_error(text: &str) -> Option<f64> {
    text.lines().find_map(|line| {
        if !line.contains("CV error") {
            return None;
        }
        line.rsplit([':', '='])
            .next()
            .and_then(|value| value.trim().parse::<f64>().ok())
    })
}

fn run_population_qc(args: &Args, primary_bed: &Path) -> AppResult<()> {
    let directory = args.output.join("population").join("structure").join("qc");
    fs::create_dir_all(&directory).map_err(|e| e.to_string())?;
    let prefix = directory.join("individuals");
    for mode in ["--missing", "--het"] {
        run_command(
            &args.plink,
            &[
                "--bfile".into(),
                primary_bed.with_extension("").display().to_string(),
                "--allow-extra-chr".into(),
                mode.into(),
                "--out".into(),
                prefix.display().to_string(),
            ],
        )?;
    }
    if !args.skip_relatedness_qc {
        run_command(
            &args.plink,
            &[
                "--bfile".into(),
                primary_bed.with_extension("").display().to_string(),
                "--allow-extra-chr".into(),
                "--genome".into(),
                "--out".into(),
                prefix.display().to_string(),
            ],
        )?;
    }
    fs::write(
        directory.join("README.txt"),
        "individuals.imiss: per-sample genotype missingness\nindividuals.het: observed/expected heterozygosity\nindividuals.genome: pairwise PI_HAT relatedness (absent with --skip-relatedness-qc)\n",
    ).map_err(|e| e.to_string())?;
    Ok(())
}

fn run_admixture(args: &Args, primary_bed: &Path, sample_count: usize) -> AppResult<()> {
    let directory = args
        .output
        .join("population")
        .join("structure")
        .join("admixture");
    fs::create_dir_all(&directory)
        .map_err(|e| format!("unable to create {}: {e}", directory.display()))?;
    let status_path = directory.join("status.tsv");
    if args.skip_admixture {
        fs::write(&status_path, "status\tdetail\nskipped\t--skip-admixture\n")
            .map_err(|e| e.to_string())?;
        return Ok(());
    }
    if !executable_available(&args.admixture) {
        fs::write(
            &status_path,
            format!(
                "status\tdetail\nunavailable\texecutable {} was not found\n",
                args.admixture
            ),
        )
        .map_err(|e| e.to_string())?;
        eprintln!(
            "Warning: ADMIXTURE executable {:?} was not found; PLINK panels and PCA are complete",
            args.admixture
        );
        return Ok(());
    }
    let k_max = args.admixture_k_max.min(sample_count);
    if k_max < args.admixture_k_min {
        fs::write(
            &status_path,
            "status\tdetail\nskipped\tnot enough samples for requested K range\n",
        )
        .map_err(|e| e.to_string())?;
        return Ok(());
    }
    let folds = args.admixture_cv.min(sample_count).max(2);
    let absolute_bed = primary_bed
        .canonicalize()
        .map_err(|e| format!("unable to resolve {}: {e}", primary_bed.display()))?;
    let mut cv_rows = Vec::new();
    for k in args.admixture_k_min..=k_max {
        let output = Command::new(&args.admixture)
            .current_dir(&directory)
            .arg(format!("--cv={folds}"))
            .arg(format!("-j{}", args.threads))
            .arg(&absolute_bed)
            .arg(k.to_string())
            .output()
            .map_err(|e| format!("unable to run {}: {e}", args.admixture))?;
        let mut log = output.stdout;
        log.extend_from_slice(&output.stderr);
        fs::write(directory.join(format!("K{k}.log")), &log).map_err(|e| e.to_string())?;
        if !output.status.success() {
            let log_text = String::from_utf8_lossy(&log);
            let detail_lines: Vec<&str> = log_text
                .lines()
                .rev()
                .filter(|line| !line.trim().is_empty())
                .take(8)
                .collect();
            let detail = detail_lines
                .into_iter()
                .rev()
                .collect::<Vec<&str>>()
                .join(" | ");
            let detail = if detail.is_empty() {
                "no ADMIXTURE diagnostic output".to_string()
            } else {
                detail
            };
            fs::write(
                &status_path,
                format!(
                    "status\tdetail\nfailed\tADMIXTURE failed for K={k} ({}): {detail}\n",
                    output.status
                ),
            )
            .map_err(|e| e.to_string())?;
            return Err(format!(
                "ADMIXTURE failed for K={k} ({}): {detail}; see K{k}.log",
                output.status
            ));
        }
        let text = String::from_utf8_lossy(&log);
        cv_rows.push((k, parse_cv_error(&text)));
    }
    let best = cv_rows
        .iter()
        .filter_map(|(k, error)| error.map(|value| (*k, value)))
        .min_by(|left, right| left.1.partial_cmp(&right.1).unwrap_or(Ordering::Equal));
    let mut cv =
        BufWriter::new(File::create(directory.join("cv_errors.tsv")).map_err(|e| e.to_string())?);
    writeln!(cv, "K\tcv_error\tbest").map_err(|e| e.to_string())?;
    for (k, error) in cv_rows {
        let is_best = best.is_some_and(|(best_k, _)| best_k == k);
        match error {
            Some(value) => writeln!(cv, "{k}\t{value:.10}\t{is_best}"),
            None => writeln!(cv, "{k}\tNA\t{is_best}"),
        }
        .map_err(|e| e.to_string())?;
    }
    fs::write(
        status_path,
        format!(
            "status\tdetail\ncomplete\tK={}-{}, CV folds={}\n",
            args.admixture_k_min, k_max, folds
        ),
    )
    .map_err(|e| e.to_string())?;
    Ok(())
}

fn run(args: &Args) -> AppResult<()> {
    let samples = read_samples(&args.samples_tsv)?;
    let reference = prepare_reference(args, &samples)?;
    println!("Population reference: {}", reference.display());
    if args.stop_after == Stage::Reference {
        return Ok(());
    }

    let bams = if args.start_at <= Stage::Mapping {
        map_samples(args, &samples, &reference)?
    } else {
        existing_bams(args, samples.len())?
    };
    if args.stop_after == Stage::Mapping {
        return Ok(());
    }

    let filtered = if args.start_at <= Stage::Calling {
        call_variants(args, &reference, &bams)?
    } else {
        existing_filtered_vcf(args)?
    };
    if args.stop_after == Stage::Calling {
        return Ok(());
    }

    let all_vcf = prepare_all_snp_panel(args, &filtered)?;
    let final_vcf = select_one_snp(args, &all_vcf)?;
    println!("All-SNP VCF: {}", all_vcf.display());
    println!("One-SNP-per-UCE VCF: {}", final_vcf.display());
    if !args.skip_plink {
        let panels = build_structure_panels(args, &all_vcf, &final_vcf)?;
        run_population_qc(args, &panels.primary_bed)?;
        run_admixture(args, &panels.primary_bed, samples.len())?;
    }
    Ok(())
}

fn main() {
    let result = parse_args(env::args().skip(1).collect()).and_then(|args| run(&args));
    if let Err(error) = result {
        eprintln!("Error: {error}");
        process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(unix)]
    fn write_executable(path: &Path, contents: &str) {
        use std::os::unix::fs::PermissionsExt;
        fs::write(path, contents).unwrap();
        let mut permissions = fs::metadata(path).unwrap().permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(path, permissions).unwrap();
    }

    #[test]
    fn parses_quoted_delimited_fields() {
        assert_eq!(
            parse_delimited_line("a,\"b,c\",d", ','),
            vec!["a", "b,c", "d"]
        );
    }

    #[test]
    fn sample_names_match_python_driver() {
        assert_eq!(sanitize_sample_name("A fancy-SAMPLE"), "A_fancy_sample");
    }

    #[test]
    fn best_supported_candidate_wins() {
        let base = Candidate {
            locus: "u1".into(),
            safe_locus: "u1".into(),
            sample: "a".into(),
            sequence: "ACGT".into(),
            supported_bases: 100,
            support_breadth: 0.9,
            max_gap: 3,
            unique_reads: 10,
            unique_density: 0.1,
            length: 100,
        };
        let mut better = base.clone();
        better.supported_bases = 101;
        assert_eq!(
            candidate_cmp(ReferenceStrategy::Supported, &better, &base),
            Ordering::Greater
        );
    }

    #[test]
    fn sqcl_strategy_prefers_longest_accepted_candidate() {
        let shorter = Candidate {
            locus: "u1".into(),
            safe_locus: "u1".into(),
            sample: "a".into(),
            sequence: "ACGT".into(),
            supported_bases: 100,
            support_breadth: 1.0,
            max_gap: 0,
            unique_reads: 100,
            unique_density: 1.0,
            length: 100,
        };
        let mut longer = shorter.clone();
        longer.sample = "b".into();
        longer.length = 150;
        longer.supported_bases = 90;
        longer.support_breadth = 0.9;
        assert_eq!(
            candidate_cmp(ReferenceStrategy::SqclLongest, &longer, &shorter),
            Ordering::Greater
        );
    }

    #[test]
    fn site_selection_prefers_call_rate_then_quality() {
        let low_call = Site {
            chrom: "u1".into(),
            pos: 1,
            reference: "A".into(),
            alternate: "G".into(),
            qual: 100.0,
            call_rate: 0.75,
            median_dp: 10.0,
            median_gq: 99.0,
            record: String::new(),
        };
        let high_call = Site {
            chrom: "u1".into(),
            pos: 2,
            reference: "A".into(),
            alternate: "T".into(),
            qual: 30.0,
            call_rate: 1.0,
            median_dp: 10.0,
            median_gq: 30.0,
            record: String::new(),
        };
        assert_eq!(site_cmp(&high_call, &low_call, 10.0), Ordering::Greater);
    }

    #[test]
    fn vcf_parser_computes_non_missing_fraction() {
        let site =
            parse_site("u1\t10\t.\tA\tG\t30\tPASS\t.\tGT:DP:GQ\t0/1:10:40\t./.:.:.\t1/1:8:30")
                .unwrap();
        assert!((site.call_rate - 2.0 / 3.0).abs() < 1e-9);
        assert_eq!(site.median_dp, 9.0);
        assert_eq!(site.median_gq, 35.0);
    }

    #[test]
    fn genotype_filter_is_applied_within_each_sample() {
        let args = Args::default();
        assert_eq!(genotype_filter_expression(&args), "FMT/DP<5 | FMT/GQ<20");
    }

    #[test]
    fn mapping_qc_parsers_compute_rates_and_weighted_depth() {
        let mut metrics = MappingMetrics::default();
        parse_flagstat(
            "100 + 0 in total (QC-passed reads + QC-failed reads)\n80 + 0 mapped (80.00% : N/A)\n60 + 0 properly paired (60.00% : N/A)\n",
            &mut metrics,
        );
        parse_coverage(
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\nu1\t1\t100\t20\t80\t80\t10\t40\t60\nu2\t1\t50\t10\t25\t50\t4\t40\t60\n",
            &mut metrics,
        );
        assert_eq!(metrics.total_reads, 100);
        assert_eq!(metrics.mapped_reads, 80);
        assert_eq!(metrics.properly_paired_reads, 60);
        assert_eq!(metrics.reference_bases, 150);
        assert_eq!(metrics.covered_bases, 105);
        assert!((metrics.mean_depth - 8.0).abs() < 1e-9);
        assert!((rate(metrics.covered_bases, metrics.reference_bases) - 0.7).abs() < 1e-9);
    }

    #[test]
    fn parses_admixture_cv_error_formats() {
        assert_eq!(parse_cv_error("CV error (K=2): 0.12345\n"), Some(0.12345));
        assert_eq!(parse_cv_error("CV error = 0.4\n"), Some(0.4));
        assert_eq!(parse_cv_error("Loglikelihood: -10\n"), None);
    }

    #[cfg(unix)]
    #[test]
    fn admixture_runner_records_cv_and_best_k() {
        use std::time::{SystemTime, UNIX_EPOCH};

        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = env::temp_dir().join(format!("gm2-admixture-test-{nonce}"));
        let structure = root.join("out/population/structure");
        fs::create_dir_all(&structure).unwrap();
        let bed = structure.join("population.bed");
        fs::write(&bed, b"bed").unwrap();
        let executable = root.join("admixture");
        write_executable(
            &executable,
            r##"#!/bin/sh
last=""
for value in "$@"; do last="$value"; done
printf 'CV error (K=%s): 0.%s\n' "$last" "$last"
: > "population.$last.Q"
: > "population.$last.P"
"##,
        );
        let args = Args {
            output: root.join("out"),
            admixture: executable.display().to_string(),
            admixture_k_min: 2,
            admixture_k_max: 3,
            ..Args::default()
        };
        run_admixture(&args, &bed, 4).unwrap();
        let cv = fs::read_to_string(structure.join("admixture/cv_errors.tsv")).unwrap();
        assert!(cv.contains("2\t0.2000000000\ttrue"));
        assert!(cv.contains("3\t0.3000000000\tfalse"));
        assert!(structure.join("admixture/population.2.Q").is_file());
        assert!(structure.join("admixture/population.3.P").is_file());
        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn partial_missing_genotypes_are_not_called() {
        for genotype in [".", "./.", ".|.", "0/.", "./1", "0|.", ".|1"] {
            assert!(!genotype_is_called(genotype), "{genotype}");
        }
        for genotype in ["0", "1", "0/0", "0/1", "1|0"] {
            assert!(genotype_is_called(genotype), "{genotype}");
        }
        let site =
            parse_site("u1\t10\t.\tA\tG\t30\tPASS\t.\tGT:DP:GQ\t0/.:10:40\t.|1:8:30\t1/1:12:50")
                .unwrap();
        assert!((site.call_rate - 1.0 / 3.0).abs() < 1e-9);
        assert_eq!(site.median_dp, 12.0);
    }

    #[test]
    fn minibwa_uses_native_map_cli_and_preserves_read_group() {
        let args = Args {
            threads: 4,
            ..Args::default()
        };
        let sample = Sample {
            original: "Coral 1".into(),
            internal: "1_Coral_1".into(),
            vcf_id: "Coral_1".into(),
            read1: "r1.fq.gz".into(),
            read2: "r2.fq.gz".into(),
            population: String::new(),
            batch: String::new(),
        };
        let argv = minibwa_map_args(args.threads, &sample, Path::new("ref.fa"));
        assert_eq!(argv[0], "map");
        assert_eq!(argv[1..4], ["-t", "4", "-R"]);
        assert!(argv[4].contains("ID:Coral_1"));
        assert!(argv[4].contains("\\t"));
        assert!(!argv[4].contains('\t'));
        assert_eq!(&argv[5..], ["ref.fa", "r1.fq.gz", "r2.fq.gz"]);
    }

    #[test]
    fn external_reference_is_materialized_with_a_source_manifest() {
        let root = env::temp_dir().join(format!("gm2-population-reference-{}", process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let source = root.join("external.fasta");
        fs::write(&source, ">uce_1\nACGT\n").unwrap();
        let args = Args {
            output: root.join("out"),
            reference_fasta: Some(source.clone()),
            ..Args::default()
        };
        let samples = vec![Sample {
            original: "A".into(),
            internal: "1_A".into(),
            vcf_id: "A".into(),
            read1: "a.fq".into(),
            read2: "a.fq".into(),
            population: String::new(),
            batch: String::new(),
        }];
        let reference = prepare_reference(&args, &samples).unwrap();
        assert_eq!(fs::read_to_string(&reference).unwrap(), ">uce_1\nACGT\n");
        let source_manifest =
            fs::read_to_string(root.join("out/population/reference/reference_source.tsv")).unwrap();
        assert!(source_manifest.contains("external"));
        assert!(source_manifest.contains(&source.display().to_string()));
        assert!(root.join("out/population/sample_manifest.tsv").is_file());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn start_stage_must_not_follow_stop_stage() {
        let error = parse_args(vec![
            "--output".into(),
            "out".into(),
            "--samples".into(),
            "samples.tsv".into(),
            "--start-at".into(),
            "calling".into(),
            "--stop-after".into(),
            "mapping".into(),
        ])
        .unwrap_err();
        assert!(error.contains("--start-at must not be later"));
    }

    #[test]
    fn reference_builder_selects_sqcl_longest_sequence_per_locus() {
        let root = env::temp_dir().join(format!("gm2-population-test-{}", process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(root.join("1_A/results")).unwrap();
        fs::create_dir_all(root.join("2_B/results")).unwrap();
        let header = "locus,status,accepted,selected_contig_length,slice_supported_bases,slice_support_breadth,max_slice_support_gap,unique_read_count,unique_read_density,low_quality\n";
        fs::write(
            root.join("1_A/uce_assembly_summary.csv"),
            format!("{header}uce1,success,true,4,4,1.0,0,10,2.5,false\n"),
        )
        .unwrap();
        fs::write(
            root.join("2_B/uce_assembly_summary.csv"),
            format!("{header}uce1,success,true,6,6,1.0,0,12,2.0,false\n"),
        )
        .unwrap();
        fs::write(root.join("1_A/results/uce1.fasta"), ">a\nACGT\n").unwrap();
        fs::write(root.join("2_B/results/uce1.fasta"), ">b\nACGTAA\n").unwrap();
        let args = Args {
            output: root.clone(),
            ..Args::default()
        };
        let samples = vec![
            Sample {
                original: "A".into(),
                internal: "1_A".into(),
                vcf_id: "A".into(),
                read1: "a.fq".into(),
                read2: "a.fq".into(),
                population: String::new(),
                batch: String::new(),
            },
            Sample {
                original: "B".into(),
                internal: "2_B".into(),
                vcf_id: "B".into(),
                read1: "b.fq".into(),
                read2: "b.fq".into(),
                population: String::new(),
                batch: String::new(),
            },
        ];
        let reference = build_reference(&args, &samples).unwrap();
        assert_eq!(fs::read_to_string(reference).unwrap(), ">uce1\nACGTAA\n");
        let provenance = fs::read_to_string(
            root.join("population/reference/population_reference_provenance.tsv"),
        )
        .unwrap();
        assert!(provenance.contains("uce1\tuce1\t2_B\tsqcl-longest\t2\t6"));
        fs::remove_dir_all(root).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn four_stage_pipeline_runs_with_tool_contract_fixtures() {
        use std::time::{SystemTime, UNIX_EPOCH};

        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = env::temp_dir().join(format!("gm2-population-e2e-{nonce}"));
        let tools = root.join("tools");
        fs::create_dir_all(root.join("out/1_A/results")).unwrap();
        fs::create_dir_all(root.join("out/2_B/results")).unwrap();
        fs::create_dir_all(&tools).unwrap();
        let samples_tsv = root.join("samples.tsv");
        fs::write(&samples_tsv, "A\ta_R1.fq\ta_R2.fq\nB\tb_R1.fq\tb_R2.fq\n").unwrap();
        let header = "locus,status,accepted,selected_contig_length,slice_supported_bases,slice_support_breadth,max_slice_support_gap,unique_read_count,unique_read_density,low_quality\n";
        let rows = format!(
            "{header}uce1,success,true,8,8,1.0,0,10,1.25,false\nuce2,success,true,8,8,1.0,0,10,1.25,false\n"
        );
        for sample in ["1_A", "2_B"] {
            fs::write(
                root.join(format!("out/{sample}/uce_assembly_summary.csv")),
                &rows,
            )
            .unwrap();
            fs::write(
                root.join(format!("out/{sample}/results/uce1.fasta")),
                ">x\nACGTACGT\n",
            )
            .unwrap();
            fs::write(
                root.join(format!("out/{sample}/results/uce2.fasta")),
                ">x\nTTTTCCCC\n",
            )
            .unwrap();
        }

        let minibwa = tools.join("minibwa");
        write_executable(
            &minibwa,
            r##"#!/bin/sh
if [ "$1" = "index" ]; then
  last=""
  for value in "$@"; do last="$value"; done
  : > "$last.l2b"
  : > "$last.mbw"
elif [ "$1" = "map" ]; then
  printf '@HD\tVN:1.6\tSO:unsorted\n'
fi
"##,
        );
        let samtools = tools.join("samtools");
        write_executable(
            &samtools,
            r##"#!/bin/sh
cmd="$1"; shift
case "$cmd" in
  faidx) : > "$1.fai" ;;
  sort)
    out=""; prev=""
    for value in "$@"; do
      if [ "$prev" = "-o" ]; then out="$value"; fi
      prev="$value"
    done
    cat > "$out"
    ;;
  fixmate) cp "$2" "$3" ;;
  markdup) cp "$1" "$2" ;;
  quickcheck) exit 0 ;;
  index)
    last=""
    for value in "$@"; do last="$value"; done
    : > "$last.bai"
    ;;
  flagstat) printf '2 + 0 in total\n2 + 0 mapped (100.00%% : N/A)\n2 + 0 properly paired (100.00%% : N/A)\n' ;;
  coverage) printf '#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\nuce1\t1\t8\t2\t8\t100\t2\t40\t60\n' ;;
esac
"##,
        );
        let bcftools = tools.join("bcftools");
        write_executable(
            &bcftools,
            r##"#!/bin/sh
cmd="$1"; shift
out=""; prev=""; plain=0; regions=0
for value in "$@"; do
  if [ "$prev" = "-o" ]; then out="$value"; fi
  if [ "$value" = "-Ov" ]; then plain=1; fi
  if [ "$value" = "-R" ]; then regions=1; fi
  prev="$value"
done
vcf() {
  printf '##fileformat=VCFv4.3\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n'
  printf 'uce1\t2\t.\tC\tT\t30\tPASS\t.\tGT:DP:GQ\t0/1:10:40\t0/0:9:35\n'
  printf 'uce1\t5\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t0/1:8:20\t1/1:8:20\n'
  printf 'uce2\t3\t.\tT\tC\t40\tPASS\t.\tGT:DP:GQ\t0/1:12:45\t0/0:10:40\n'
}
case "$cmd" in
  mpileup) printf 'pileup\n' ;;
  call|norm|+setGT|+fill-tags)
    if [ -n "$out" ]; then cat > "$out"; else cat; fi
    ;;
  view)
    if [ "$plain" = 1 ] || [ "$regions" = 1 ]; then
      vcf
    elif [ -n "$out" ]; then
      cat > "$out"
    else
      cat
    fi
    ;;
  annotate)
    last=""
    for value in "$@"; do last="$value"; done
    cp "$last" "$out"
    ;;
  index)
    if [ "$1" = "-n" ]; then
      printf '3\n'
    else
      last=""
      for value in "$@"; do last="$value"; done
      : > "$last.csi"
    fi
    ;;
esac
"##,
        );
        let plink = tools.join("plink");
        write_executable(
            &plink,
            r##"#!/bin/sh
out=""; prev=""; pca=0; prune=0
for value in "$@"; do
  if [ "$prev" = "--out" ]; then out="$value"; fi
  if [ "$value" = "--pca" ]; then pca=1; fi
  if [ "$value" = "--indep-pairwise" ]; then prune=1; fi
  prev="$value"
done
if [ "$pca" = 1 ]; then
  : > "$out.eigenvec"; : > "$out.eigenval"
elif [ "$prune" = 1 ]; then
  printf 'uce1:2:C:T\nuce1:5:A:G\nuce2:3:T:C\n' > "$out.prune.in"
  : > "$out.prune.out"
else
  : > "$out.bed"; : > "$out.bim"; : > "$out.fam"
fi
"##,
        );

        let args = Args {
            output: root.join("out"),
            samples_tsv,
            minibwa: minibwa.display().to_string(),
            samtools: samtools.display().to_string(),
            bcftools: bcftools.display().to_string(),
            plink: plink.display().to_string(),
            skip_admixture: true,
            ..Args::default()
        };
        run(&args).unwrap();
        let selected =
            fs::read_to_string(root.join("out/population/structure/selected_snps.tsv")).unwrap();
        assert!(selected.contains("uce1\t2\tC\tT"));
        assert!(selected.contains("uce2\t3\tT\tC"));
        assert!(root
            .join("out/population/structure/one_snp_per_uce.vcf.gz")
            .is_file());
        assert!(root
            .join("out/population/structure/population.bed")
            .is_file());
        assert!(root
            .join("out/population/structure/population_pca.eigenvec")
            .is_file());
        assert!(root
            .join("out/population/structure/all_snps.vcf.gz")
            .is_file());
        assert!(root
            .join("out/population/structure/ld_pruned.vcf.gz")
            .is_file());
        assert!(root
            .join("out/population/structure/panel_summary.tsv")
            .is_file());
        let mapping_qc =
            fs::read_to_string(root.join("out/population/mapping/mapping_qc.tsv")).unwrap();
        assert!(mapping_qc.contains("mapping_rate"));
        assert!(mapping_qc.contains("1.000000"));
        let variant_qc =
            fs::read_to_string(root.join("out/population/variants/variant_qc.tsv")).unwrap();
        assert!(variant_qc.contains("raw_called_variants\t3"));
        assert!(variant_qc.contains("site_filtered\t3"));
        fs::remove_dir_all(root).unwrap();
    }
}
