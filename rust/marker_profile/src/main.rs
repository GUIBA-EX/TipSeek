use std::collections::{BTreeMap, BTreeSet};

use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread;
use std::time::Duration;

#[derive(Debug)]
struct Args {
    reference: PathBuf,
    reads: PathBuf,
    output: PathBuf,
    cache: PathBuf,
    themisto: PathBuf,
    groups: Option<PathBuf>,
    decoy: Option<PathBuf>,
    threads: usize,
    kmer_size: usize,
    threshold: f64,
    relevant_kmer_fraction: f64,
    index_memory_gb: usize,
    force_rebuild: bool,
}

struct CacheLock(PathBuf);

impl CacheLock {
    fn acquire(path: PathBuf) -> Result<Self, String> {
        if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
            fs::create_dir_all(parent).map_err(|e| format!("cannot create cache parent: {e}"))?;
        }
        for _ in 0..600 {
            match fs::create_dir(&path) {
                Ok(()) => return Ok(Self(path)),
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => {
                    thread::sleep(Duration::from_millis(100))
                }
                Err(error) => {
                    return Err(format!(
                        "cannot acquire marker reference-cache lock: {error}"
                    ))
                }
            }
        }
        Err("timed out waiting for marker reference-cache construction".to_string())
    }
}

impl Drop for CacheLock {
    fn drop(&mut self) {
        let _ = fs::remove_dir(&self.0);
    }
}

fn usage() -> &'static str {
    "Usage: marker_profile --reference REF.fa --reads recruited.fq --output DIR --cache DIR \\
  --themisto PATH [options]\n\n\
Options:\n\
  --groups FILE                    Optional TSV: reference ID, then reporting group\n\
  --decoy FILE                     Optional non-target decoy FASTA\n\
  --threads INT                    Worker threads (default: 1)\n\
  --kmer-size INT                  Themisto k-mer size (default: 21)\n\
  --threshold FLOAT                Themisto pseudoalignment threshold (default: 0.80)\n\
  --relevant-kmer-fraction FLOAT   Minimum fraction of query k-mers found in any target (default: 0.50)\n\
  --index-memory-gb INT            Themisto build memory limit (default: 2)\n\
  --force-rebuild                  Rebuild the cached reference index\n"
}

fn next_value(args: &[String], pos: &mut usize, option: &str) -> Result<String, String> {
    *pos += 1;
    args.get(*pos)
        .cloned()
        .ok_or_else(|| format!("{option} requires a value"))
}

fn parse_args() -> Result<Args, String> {
    let raw: Vec<String> = env::args().skip(1).collect();
    if raw.iter().any(|arg| arg == "--help" || arg == "-h") {
        return Err(usage().to_string());
    }
    let mut reference = None;
    let mut reads = None;
    let mut output = None;
    let mut cache = None;
    let mut themisto = None;
    let mut groups = None;
    let mut decoy = None;
    let mut threads = 1usize;
    let mut kmer_size = 21usize;
    let mut threshold = 0.80f64;
    let mut relevant_kmer_fraction = 0.50f64;
    let mut index_memory_gb = 2usize;
    let mut force_rebuild = false;
    let mut pos = 0;
    while pos < raw.len() {
        match raw[pos].as_str() {
            "--reference" => {
                reference = Some(PathBuf::from(next_value(&raw, &mut pos, "--reference")?))
            }
            "--reads" => reads = Some(PathBuf::from(next_value(&raw, &mut pos, "--reads")?)),
            "--output" => output = Some(PathBuf::from(next_value(&raw, &mut pos, "--output")?)),
            "--cache" => cache = Some(PathBuf::from(next_value(&raw, &mut pos, "--cache")?)),
            "--themisto" => {
                themisto = Some(PathBuf::from(next_value(&raw, &mut pos, "--themisto")?))
            }
            "--groups" => groups = Some(PathBuf::from(next_value(&raw, &mut pos, "--groups")?)),
            "--decoy" => decoy = Some(PathBuf::from(next_value(&raw, &mut pos, "--decoy")?)),
            "--threads" => {
                threads = next_value(&raw, &mut pos, "--threads")?
                    .parse()
                    .map_err(|_| "invalid --threads".to_string())?
            }
            "--kmer-size" => {
                kmer_size = next_value(&raw, &mut pos, "--kmer-size")?
                    .parse()
                    .map_err(|_| "invalid --kmer-size".to_string())?
            }
            "--threshold" => {
                threshold = next_value(&raw, &mut pos, "--threshold")?
                    .parse()
                    .map_err(|_| "invalid --threshold".to_string())?
            }
            "--relevant-kmer-fraction" => {
                relevant_kmer_fraction = next_value(&raw, &mut pos, "--relevant-kmer-fraction")?
                    .parse()
                    .map_err(|_| "invalid --relevant-kmer-fraction".to_string())?
            }
            "--index-memory-gb" => {
                index_memory_gb = next_value(&raw, &mut pos, "--index-memory-gb")?
                    .parse()
                    .map_err(|_| "invalid --index-memory-gb".to_string())?
            }
            "--force-rebuild" => force_rebuild = true,
            other => return Err(format!("unknown option: {other}\n\n{}", usage())),
        }
        pos += 1;
    }
    let args = Args {
        reference: reference.ok_or_else(|| "--reference is required".to_string())?,
        reads: reads.ok_or_else(|| "--reads is required".to_string())?,
        output: output.ok_or_else(|| "--output is required".to_string())?,
        cache: cache.ok_or_else(|| "--cache is required".to_string())?,
        themisto: themisto.ok_or_else(|| "--themisto is required".to_string())?,
        groups,
        decoy,
        threads,
        kmer_size,
        threshold,
        relevant_kmer_fraction,
        index_memory_gb,
        force_rebuild,
    };
    if args.threads == 0
        || args.kmer_size < 15
        || args.kmer_size > 31
        || args.kmer_size.is_multiple_of(2)
        || !(0.0..=1.0).contains(&args.threshold)
        || !(0.0..=1.0).contains(&args.relevant_kmer_fraction)
    {
        return Err("invalid marker quantification parameter".to_string());
    }
    Ok(args)
}

fn load_group_map(path: &Path) -> Result<BTreeMap<String, String>, String> {
    let input =
        File::open(path).map_err(|e| format!("cannot open group map {}: {e}", path.display()))?;
    let mut map = BTreeMap::new();
    for (line_no, line) in BufReader::new(input).lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        if line.is_empty() || line.starts_with("#") {
            continue;
        }
        let mut fields = line.split("\t");
        let id = fields.next().unwrap_or("").trim();
        let group = fields.next().unwrap_or("").trim();
        if id.is_empty() || group.is_empty() || fields.next().is_some() {
            return Err(format!(
                "invalid group-map line {}: expected reference_id<TAB>group",
                line_no + 1
            ));
        }
        if let Some(existing) = map.get(id) {
            if existing != group {
                return Err(format!("conflicting groups for reference ID: {id}"));
            }
        } else {
            map.insert(id.to_string(), group.to_string());
        }
    }
    if map.is_empty() {
        return Err("group map is empty".to_string());
    }
    Ok(map)
}

fn validate_group_map_coverage(
    reference: &Path,
    map: &BTreeMap<String, String>,
) -> Result<(), String> {
    let input =
        File::open(reference).map_err(|e| format!("cannot open {}: {e}", reference.display()))?;
    let mut reference_ids = BTreeSet::new();
    for line in BufReader::new(input).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = header.split_whitespace().next() {
                reference_ids.insert(id.to_string());
            }
        }
    }
    let mapped_ids: BTreeSet<String> = map.keys().cloned().collect();
    let missing: Vec<_> = reference_ids
        .difference(&mapped_ids)
        .take(10)
        .cloned()
        .collect();
    let unused: Vec<_> = mapped_ids
        .difference(&reference_ids)
        .take(10)
        .cloned()
        .collect();
    if !missing.is_empty() || !unused.is_empty() {
        return Err(format!(
            "group-map coverage mismatch; missing reference IDs: {:?}; unused mapped IDs: {:?}",
            missing, unused
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_record(
    dir: &Path,
    id: usize,
    header: &str,
    sequence: &str,
    group: &str,
    list: &mut BufWriter<File>,
    metadata: &mut BufWriter<File>,
) -> io::Result<()> {
    let record_id = format!("marker_{id:06}");
    let record_path = dir.join(format!("{record_id}.fasta"));
    let mut record = BufWriter::new(File::create(&record_path)?);
    writeln!(record, ">{record_id}")?;
    writeln!(record, "{}", sequence)?;
    writeln!(list, "{}", record_path.display())?;
    writeln!(
        metadata,
        "{id}\t{record_id}\t{group}\t{}",
        header.replace('\t', " ")
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn append_fasta(
    path: &Path,
    records_dir: &Path,
    next_id: &mut usize,
    list: &mut BufWriter<File>,
    metadata: &mut BufWriter<File>,
    group_map: Option<&BTreeMap<String, String>>,
    decoy: bool,
) -> Result<(usize, usize), String> {
    let input = File::open(path).map_err(|e| format!("cannot open {}: {e}", path.display()))?;
    let mut header = String::new();
    let mut sequence = String::new();
    let mut kept = 0usize;
    let skipped = 0usize;
    let mut flush = |header: &str, sequence: &str| -> Result<(), String> {
        if header.is_empty() || sequence.is_empty() {
            return Ok(());
        }
        let group = if decoy {
            "DECOY".to_string()
        } else {
            let id = header.split_whitespace().next().unwrap_or("");
            group_map
                .and_then(|map| map.get(id))
                .cloned()
                .unwrap_or_else(|| id.to_string())
        };
        if !decoy && group == "DECOY" {
            return Err(
                "DECOY is reserved for --decoy references; rename the target group or ID".into(),
            );
        }
        write_record(
            records_dir,
            *next_id,
            header,
            sequence,
            &group,
            list,
            metadata,
        )
        .map_err(|e| e.to_string())?;
        *next_id += 1;
        kept += 1;
        Ok(())
    };
    for line in BufReader::new(input).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if let Some(next_header) = line.strip_prefix('>') {
            flush(&header, &sequence)?;
            header = next_header.trim().to_string();
            sequence.clear();
        } else {
            sequence.push_str(line.trim());
        }
    }
    flush(&header, &sequence)?;
    Ok((kept, skipped))
}

fn run(command: &mut Command, label: &str) -> Result<(), String> {
    let status = command
        .stdout(Stdio::inherit())
        .stderr(Stdio::inherit())
        .status()
        .map_err(|e| format!("cannot launch {label}: {e}"))?;
    status
        .success()
        .then_some(())
        .ok_or_else(|| format!("{label} failed with {status}"))
}

fn prepare_reference(args: &Args) -> Result<(PathBuf, PathBuf), String> {
    let _lock = CacheLock::acquire(args.cache.with_extension("lock"))?;
    let records = args.cache.join("records");
    let list_path = args.cache.join("reference_files.txt");
    let metadata_path = args.cache.join("marker_reference_metadata.tsv");
    let index_prefix = args.cache.join("themisto_index");
    let index_ready = index_prefix.with_extension("tdbg").is_file()
        && index_prefix.with_extension("tcolors").is_file();
    if !args.force_rebuild && index_ready && list_path.is_file() && metadata_path.is_file() {
        return Ok((index_prefix, metadata_path));
    }
    if args.cache.exists() {
        fs::remove_dir_all(&args.cache).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(&records).map_err(|e| e.to_string())?;
    let mut list = BufWriter::new(File::create(&list_path).map_err(|e| e.to_string())?);
    let mut metadata = BufWriter::new(File::create(&metadata_path).map_err(|e| e.to_string())?);
    writeln!(metadata, "color\treference_id\tgroup\toriginal_header").map_err(|e| e.to_string())?;
    let mut next_id = 0usize;
    let group_map = match &args.groups {
        Some(path) => {
            let map = load_group_map(path)?;
            validate_group_map_coverage(&args.reference, &map)?;
            Some(map)
        }
        None => None,
    };
    let (kept, _) = append_fasta(
        &args.reference,
        &records,
        &mut next_id,
        &mut list,
        &mut metadata,
        group_map.as_ref(),
        false,
    )?;
    if let Some(decoy) = &args.decoy {
        append_fasta(
            decoy,
            &records,
            &mut next_id,
            &mut list,
            &mut metadata,
            None,
            true,
        )?;
    }
    list.flush().map_err(|e| e.to_string())?;
    metadata.flush().map_err(|e| e.to_string())?;
    if kept == 0 {
        return Err("no marker reference records were loaded".to_string());
    }
    let temp_dir = args.cache.join("temp");
    fs::create_dir_all(&temp_dir).map_err(|e| e.to_string())?;
    let mut command = Command::new(&args.themisto);
    command
        .args(["build", "-k", &args.kmer_size.to_string(), "-i"])
        .arg(&list_path)
        .args(["--index-prefix"])
        .arg(&index_prefix)
        .args(["--temp-dir"])
        .arg(&temp_dir)
        .args([
            "--mem-gigas",
            &args.index_memory_gb.to_string(),
            "--n-threads",
            &args.threads.to_string(),
            "--file-colors",
        ]);
    run(&mut command, "Themisto index construction")?;
    Ok((index_prefix, metadata_path))
}

fn write_reference_support(
    output: &Path,
    pseudoalignments: &Path,
    metadata: &Path,
) -> Result<(), String> {
    let mut colors = BTreeMap::<usize, (String, String)>::new();
    for (line_no, line) in BufReader::new(File::open(metadata).map_err(|e| e.to_string())?)
        .lines()
        .enumerate()
    {
        let line = line.map_err(|e| e.to_string())?;
        if line_no == 0 {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 4 {
            return Err(format!("invalid reference metadata line {}", line_no + 1));
        }
        let color = fields[0]
            .parse::<usize>()
            .map_err(|_| format!("invalid reference color: {}", fields[0]))?;
        let reference_id = fields[3]
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        if reference_id.is_empty() {
            return Err(format!(
                "empty reference ID at metadata line {}",
                line_no + 1
            ));
        }
        colors.insert(color, (reference_id, fields[2].to_string()));
    }
    let has_decoy = colors.values().any(|(_, group)| group == "DECOY");
    let mut sources = BufWriter::new(
        File::create(output.join("marker_query_sources.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(
        sources,
        "query_index\ttarget_hits\tdecoy_hits\tsource_status\tgroups"
    )
    .map_err(|e| e.to_string())?;
    let mut support = BTreeMap::<(String, String), (String, usize, f64, usize)>::new();
    for line in BufReader::new(File::open(pseudoalignments).map_err(|e| e.to_string())?).lines() {
        let line = line.map_err(|e| e.to_string())?;
        let query = line
            .split_whitespace()
            .next()
            .ok_or("empty pseudoalignment record")?
            .parse::<usize>()
            .map_err(|_| "invalid pseudoalignment query index")?;
        let candidates: BTreeSet<usize> = line
            .split_whitespace()
            .skip(1)
            .map(|field| {
                field
                    .parse::<usize>()
                    .map_err(|_| format!("invalid Themisto color: {field}"))
            })
            .collect::<Result<_, _>>()?;
        let mut target_hits = 0;
        let mut decoy_hits = 0;
        let mut groups = BTreeSet::new();
        for color in &candidates {
            let (_, group) = colors
                .get(color)
                .ok_or_else(|| format!("Themisto color {color} has no reference metadata"))?;
            if group == "DECOY" {
                decoy_hits += 1;
            } else {
                target_hits += 1;
                groups.insert(group.as_str());
            }
        }
        let status = match (target_hits > 0, decoy_hits > 0, has_decoy) {
            (true, true, _) => "ambiguous_source",
            (false, true, _) => "decoy_only",
            (true, false, true) => "target_only",
            (true, false, false) => "target_match_no_decoy",
            _ => "unclassified",
        };
        writeln!(
            sources,
            "{query}\t{target_hits}\t{decoy_hits}\t{status}\t{}",
            groups.into_iter().collect::<Vec<_>>().join(";")
        )
        .map_err(|e| e.to_string())?;
        if candidates.is_empty() {
            continue;
        }
        let weight = 1.0 / candidates.len() as f64;
        for color in candidates {
            let (reference_id, group) = colors
                .get(&color)
                .ok_or_else(|| format!("Themisto color {color} has no reference metadata"))?;
            let entry = support
                .entry((reference_id.clone(), group.clone()))
                .or_insert_with(|| (group.clone(), 0, 0.0, 0));
            entry.1 += 1;
            entry.2 += weight;
            if (1.0 - weight).abs() < f64::EPSILON {
                entry.3 += 1;
            }
        }
    }
    let mut out = BufWriter::new(
        File::create(output.join("marker_reference_support.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(
        out,
        "reference_id\tgroup\thit_queries\tfractional_queries\tsingleton_queries\tambiguity_status"
    )
    .map_err(|e| e.to_string())?;
    for ((reference_id, _), (group, hits, fractional, singleton)) in support {
        let status = if singleton > 0 {
            "has_singleton_support"
        } else {
            "shared_only"
        };
        writeln!(
            out,
            "{reference_id}\t{group}\t{hits}\t{fractional:.8}\t{singleton}\t{status}"
        )
        .map_err(|e| e.to_string())?;
    }
    out.flush().map_err(|e| e.to_string())?;
    sources.flush().map_err(|e| e.to_string())?;
    Ok(())
}

fn write_reference_qc(
    output: &Path,
    pseudoalignments: &Path,
    kmer_size: usize,
    threshold: f64,
    relevant_kmer_fraction: f64,
) -> Result<(), String> {
    let mut queries = 0usize;
    let mut positive = 0usize;
    for line in BufReader::new(File::open(pseudoalignments).map_err(|e| e.to_string())?).lines() {
        let line = line.map_err(|e| e.to_string())?;
        queries += 1;
        if line.split_whitespace().nth(1).is_some() {
            positive += 1;
        }
    }
    let mut qc =
        BufWriter::new(File::create(output.join("marker_qc.tsv")).map_err(|e| e.to_string())?);
    writeln!(qc, "metric\tvalue").map_err(|e| e.to_string())?;
    writeln!(qc, "pseudoaligned_queries\t{queries}").map_err(|e| e.to_string())?;
    writeln!(qc, "queries_with_reference_hits\t{positive}").map_err(|e| e.to_string())?;
    writeln!(qc, "kmer_size\t{kmer_size}").map_err(|e| e.to_string())?;
    writeln!(qc, "pseudoalign_threshold\t{threshold}").map_err(|e| e.to_string())?;
    writeln!(qc, "relevant_kmer_fraction\t{relevant_kmer_fraction}").map_err(|e| e.to_string())?;
    writeln!(qc, "abundance_method\treference_support").map_err(|e| e.to_string())?;
    qc.flush().map_err(|e| e.to_string())?;
    Ok(())
}

fn main() -> Result<(), String> {
    if env::args()
        .skip(1)
        .any(|arg| arg == "--help" || arg == "-h")
    {
        print!("{}", usage());
        return Ok(());
    }
    let args = parse_args()?;
    if !args.reference.is_file() || !args.reads.is_file() {
        return Err("reference and reads paths must be existing files".to_string());
    }
    let (index_prefix, metadata) = prepare_reference(&args)?;
    if args.output.exists()
        && fs::read_dir(&args.output)
            .map_err(|e| e.to_string())?
            .next()
            .is_some()
    {
        return Err(format!(
            "output directory is not empty: {}",
            args.output.display()
        ));
    }
    fs::create_dir_all(&args.output).map_err(|e| e.to_string())?;
    let pseudoalignments = args.output.join("themisto_pseudoalignments.txt");
    let temp_dir = args.output.join("themisto_temp");
    fs::create_dir_all(&temp_dir).map_err(|e| e.to_string())?;
    let mut pseudoalign = Command::new(&args.themisto);
    pseudoalign
        .args(["pseudoalign", "--query-file"])
        .arg(&args.reads)
        .args(["--index-prefix"])
        .arg(&index_prefix)
        .args(["--temp-dir"])
        .arg(&temp_dir)
        .args(["--out-file"])
        .arg(&pseudoalignments)
        .args([
            "--sort-output-lines",
            "--threshold",
            &args.threshold.to_string(),
            "--relevant-kmers-fraction",
            &args.relevant_kmer_fraction.to_string(),
            "--n-threads",
            &args.threads.to_string(),
        ]);
    run(&mut pseudoalign, "Themisto pseudoalignment")?;
    write_reference_qc(
        &args.output,
        &pseudoalignments,
        args.kmer_size,
        args.threshold,
        args.relevant_kmer_fraction,
    )?;
    write_reference_support(&args.output, &pseudoalignments, &metadata)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_file(name: &str, contents: &str) -> PathBuf {
        let path = env::temp_dir().join(format!(
            "marker_profile_{}_{}_{}",
            std::process::id(),
            name,
            contents.len()
        ));
        fs::write(&path, contents).unwrap();
        path
    }

    #[test]
    fn target_and_decoy_with_same_id_do_not_merge_support() {
        let root = env::temp_dir().join(format!("marker_role_test_{}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        let metadata = root.join("metadata.tsv");
        let alignments = root.join("alignments.txt");
        fs::write(&metadata, "color\treference_id\tgroup\toriginal_header\n0\tmarker_0\thost\tsame\n1\tmarker_1\tDECOY\tsame\n").unwrap();
        fs::write(&alignments, "0 0\n1 1\n2 0 1\n3\n").unwrap();
        write_reference_support(&root, &alignments, &metadata).unwrap();
        let support = fs::read_to_string(root.join("marker_reference_support.tsv")).unwrap();
        assert!(support.contains("same\thost\t2\t"));
        assert!(support.contains("same\tDECOY\t2\t"));
        let sources = fs::read_to_string(root.join("marker_query_sources.tsv")).unwrap();
        assert!(sources.contains("0\t1\t0\ttarget_only\thost"));
        assert!(sources.contains("1\t0\t1\tdecoy_only\t"));
        assert!(sources.contains("2\t1\t1\tambiguous_source\thost"));
        assert!(sources.contains("3\t0\t0\tunclassified\t"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn cache_lock_creates_missing_parent_and_releases() {
        let root = env::temp_dir().join(format!("marker_lock_test_{}", std::process::id()));
        let path = root.join("new/cache.lock");
        assert!(!root.exists());
        {
            let _lock = CacheLock::acquire(path.clone()).unwrap();
            assert!(path.is_dir());
        }
        assert!(!path.exists());
        let _lock = CacheLock::acquire(path).unwrap();
        drop(_lock);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn group_map_accepts_consistent_duplicates() {
        let path = temp_file("duplicates.tsv", "r1\tA\nr1\tA\nr2\tB\n");
        let map = load_group_map(&path).unwrap();
        assert_eq!(map.len(), 2);
        assert_eq!(map.get("r1").map(String::as_str), Some("A"));
        let _ = fs::remove_file(path);
    }

    #[test]
    fn group_map_rejects_conflicting_duplicates() {
        let path = temp_file("conflict.tsv", "r1\tA\nr1\tB\n");
        assert!(load_group_map(&path)
            .unwrap_err()
            .contains("conflicting groups"));
        let _ = fs::remove_file(path);
    }

    #[test]
    fn group_map_requires_exact_reference_coverage() {
        let reference = temp_file("reference.fa", ">r1 first\nACGT\n>r2\nACGT\n");
        let mut map = BTreeMap::new();
        map.insert("r1".to_string(), "A".to_string());
        map.insert("unused".to_string(), "B".to_string());
        let error = validate_group_map_coverage(&reference, &map).unwrap_err();
        assert!(error.contains("r2"));
        assert!(error.contains("unused"));
        let _ = fs::remove_file(reference);
    }
}
