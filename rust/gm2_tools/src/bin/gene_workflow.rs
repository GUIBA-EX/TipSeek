//! Conservative first-pass gene-family summaries for TipSeek.
//!
//! This program deliberately does not infer phased alleles or biological copy
//! numbers. It turns the ranked candidate contigs emitted by original-rust
//! into reproducible candidate-state calls and cohort-level FASTA/TSV products.

#[path = "../resolve.rs"]
mod resolve;

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

const FASTA_EXTENSIONS: &[&str] = &["fa", "fas", "fasta"];

#[derive(Clone, Debug)]
struct Candidate {
    sequence: String,
    key: String,
    assembly_metrics: [String; 4],
}

#[derive(Clone, Debug)]
struct Call {
    state: &'static str,
}

fn usage() -> ! {
    eprintln!(
        "Usage:\n  gene_workflow classify --reference DIR --contigs DIR --sample NAME --out DIR\n  gene_workflow cohort --reference DIR --out DIR --sample NAME [--sample NAME ...]\n  gene_workflow annotate --input DIR --protein-reference DIR --out DIR --miniprot FILE [--threads N]\n  gene_workflow resolve --input DIR --out DIR --mafft FILE --iqtree FILE --min-taxa N [--threads N] [--outgroup FILE] [--ufboot N] [--min-aa-length N] [--min-effective-codon-sites N] [--taper-script FILE --julia FILE]"
    );
    std::process::exit(2);
}

fn take_value(args: &[String], index: &mut usize, flag: &str) -> String {
    *index += 1;
    args.get(*index).cloned().unwrap_or_else(|| {
        eprintln!("Missing value for {flag}");
        usage();
    })
}

fn parse_options(args: &[String]) -> HashMap<String, Vec<String>> {
    let mut options: HashMap<String, Vec<String>> = HashMap::new();
    let mut index = 0;
    while index < args.len() {
        let flag = &args[index];
        if !flag.starts_with("--") {
            eprintln!("Unexpected argument: {flag}");
            usage();
        }
        let value = take_value(args, &mut index, flag);
        options.entry(flag.clone()).or_default().push(value);
        index += 1;
    }
    options
}

fn option_path(options: &HashMap<String, Vec<String>>, name: &str) -> PathBuf {
    options
        .get(name)
        .and_then(|values| values.first())
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            eprintln!("Missing required option {name}");
            usage();
        })
}

fn option_string(options: &HashMap<String, Vec<String>>, name: &str) -> String {
    options
        .get(name)
        .and_then(|values| values.first())
        .cloned()
        .unwrap_or_else(|| {
            eprintln!("Missing required option {name}");
            usage();
        })
}

fn option_positive_usize(
    options: &HashMap<String, Vec<String>>,
    name: &str,
    default: usize,
) -> usize {
    let value = options
        .get(name)
        .and_then(|values| values.first())
        .map(String::as_str)
        .unwrap_or("");
    if value.is_empty() {
        return default;
    }
    match value.parse::<usize>() {
        Ok(n) if n > 0 => n,
        _ => {
            eprintln!("{name} must be a positive integer");
            usage();
        }
    }
}

fn family_id(path: &Path) -> Option<String> {
    let extension = path.extension()?.to_str()?.to_ascii_lowercase();
    FASTA_EXTENSIONS
        .contains(&extension.as_str())
        .then(|| path.file_stem()?.to_str().map(str::to_owned))
        .flatten()
}

fn family_lengths(reference: &Path) -> io::Result<BTreeMap<String, usize>> {
    let mut families = BTreeMap::new();
    for entry in fs::read_dir(reference)? {
        let path = entry?.path();
        let Some(id) = family_id(&path) else {
            continue;
        };
        let longest = read_fasta(&path)?
            .into_iter()
            .map(|(_, sequence)| sequence.len())
            .max()
            .unwrap_or(0);
        if longest > 0 {
            families.insert(id, longest);
        }
    }
    if families.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "reference directory contains no non-empty FASTA families",
        ));
    }
    Ok(families)
}

fn read_fasta(path: &Path) -> io::Result<Vec<(String, String)>> {
    let mut records = Vec::new();
    let reader = BufReader::new(File::open(path)?);
    let mut header: Option<String> = None;
    let mut sequence = String::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(previous) = header.take() {
                records.push((previous, normalize_sequence(&sequence)));
            }
            header = Some(rest.trim().to_owned());
            sequence.clear();
        } else {
            sequence.push_str(line.trim());
        }
    }
    if let Some(previous) = header {
        records.push((previous, normalize_sequence(&sequence)));
    }
    Ok(records)
}

fn normalize_sequence(value: &str) -> String {
    value
        .bytes()
        .filter_map(|base| match base.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' => Some(base.to_ascii_uppercase() as char),
            b'U' => Some('T'),
            _ => None,
        })
        .collect()
}

fn read_raw_fasta(path: &Path) -> io::Result<Vec<(String, String)>> {
    let mut records = Vec::new();
    let reader = BufReader::new(File::open(path)?);
    let mut header = None;
    let mut sequence = String::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(previous) = header.take() {
                records.push((previous, sequence.clone()))
            };
            header = Some(rest.trim().to_owned());
            sequence.clear();
        } else {
            sequence.push_str(line.trim());
        }
    }
    if let Some(previous) = header {
        records.push((previous, sequence));
    }
    Ok(records)
}
fn codon_aa(codon: &[u8]) -> char {
    fn b(x: u8) -> Option<usize> {
        match x {
            b'T' => Some(0),
            b'C' => Some(1),
            b'A' => Some(2),
            b'G' => Some(3),
            _ => None,
        }
    }
    let (Some(a), Some(b), Some(c)) = (b(codon[0]), b(codon[1]), b(codon[2])) else {
        return 'X';
    };
    const TABLE: &str = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    TABLE.as_bytes()[a * 16 + b * 4 + c] as char
}
fn translate_cds(sequence: &str) -> String {
    sequence.as_bytes().chunks_exact(3).map(codon_aa).collect()
}
fn codon_backtranslate(aligned: &str, cds: &str) -> Option<String> {
    let mut offset = 0usize;
    let mut out = String::new();
    for aa in aligned.bytes() {
        if aa == b'-' {
            out.push_str("---")
        } else {
            let codon = cds.get(offset..offset + 3)?;
            out.push_str(codon);
            offset += 3
        }
    }
    Some(out)
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

fn canonical(sequence: String) -> String {
    let reverse = reverse_complement(&sequence);
    if reverse < sequence {
        reverse
    } else {
        sequence
    }
}

fn assembly_metrics(header: &str) -> [String; 4] {
    let parts: Vec<_> = header.split('_').collect();
    if parts.len() == 6 && parts[0] == "contig" {
        return [
            parts[2].to_owned(),
            parts[3].to_owned(),
            parts[4].to_owned(),
            parts[5].to_owned(),
        ];
    }
    [String::new(), String::new(), String::new(), String::new()]
}

fn is_strict_path_prefix(longer: &Candidate, shorter: &Candidate) -> bool {
    longer.key.len() > shorter.key.len() && longer.key.starts_with(&shorter.key)
}

fn unique_candidates(path: &Path, reference_length: usize) -> io::Result<(Vec<Candidate>, String)> {
    if !path.is_file() {
        return Ok((Vec::new(), "no_candidates".to_owned()));
    }
    let raw = read_fasta(path)?;
    let minimum_length = (reference_length / 3).max(40);
    let mut flags = BTreeSet::new();
    if raw.is_empty() {
        flags.insert("empty_candidates");
    }
    let mut seen = HashSet::new();
    let mut candidates = Vec::new();
    for (header, sequence) in raw {
        if sequence.len() < minimum_length {
            flags.insert("short_candidate");
            continue;
        }
        let key = canonical(sequence.clone());
        if !seen.insert(key.clone()) {
            flags.insert("duplicate_candidate");
            continue;
        }
        candidates.push(Candidate {
            sequence,
            key,
            assembly_metrics: assembly_metrics(&header),
        });
    }
    candidates.sort_by(|left, right| {
        right
            .sequence
            .len()
            .cmp(&left.sequence.len())
            .then_with(|| left.key.cmp(&right.key))
    });

    // Branch enumeration can emit strict prefixes of the same assembled path.
    // Do not treat an arbitrary internal repeat as containment: retaining that
    // ambiguity is safer than collapsing an independent candidate.
    let mut retained = Vec::new();
    for candidate in candidates {
        if retained
            .iter()
            .any(|longer: &Candidate| is_strict_path_prefix(longer, &candidate))
        {
            flags.insert("contained_candidate");
        } else {
            retained.push(candidate);
        }
    }
    Ok((retained, flags.into_iter().collect::<Vec<_>>().join(";")))
}

fn write_candidates(
    path: &Path,
    sample: &str,
    family: &str,
    candidates: &[Candidate],
    evidence: &mut BufWriter<File>,
) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut output = BufWriter::new(File::create(path)?);
    for (index, candidate) in candidates.iter().enumerate() {
        let candidate_id = format!("candidate_{}", index + 1);
        writeln!(output, ">{sample}|{family}|{candidate_id}")?;
        writeln!(output, "{}", candidate.sequence)?;
        writeln!(
            evidence,
            "{sample}\t{family}\t{candidate_id}\t{}\t{}\t{}\t{}\t{}",
            candidate.sequence.len(),
            candidate.assembly_metrics[0],
            candidate.assembly_metrics[1],
            candidate.assembly_metrics[2],
            candidate.assembly_metrics[3],
        )?;
    }
    Ok(())
}

fn classify(reference: PathBuf, contigs: PathBuf, sample: String, out: PathBuf) -> io::Result<()> {
    let families = family_lengths(&reference)?;
    let sample_root = out.join("samples").join(&sample);
    let candidate_dir = sample_root.join("candidates");
    fs::create_dir_all(&candidate_dir)?;
    let mut calls = BufWriter::new(File::create(sample_root.join("family_calls.tsv"))?);
    let mut evidence = BufWriter::new(File::create(sample_root.join("candidate_assembly.tsv"))?);
    writeln!(calls, "sample\tfamily_id\tstate\tcandidate_count\tqc_flags")?;
    writeln!(evidence, "sample\tfamily_id\tcandidate\tassembled_contig_length\tassembler_seed_count\tassembler_seed_position\tassembler_path_weight\tassembler_slice_support")?;

    for (family, reference_length) in families {
        let input = contigs.join(format!("{family}.fasta"));
        let (candidates, flags) = unique_candidates(&input, reference_length)?;
        let state = if !input.is_file() {
            "not_recovered"
        } else if candidates.is_empty() {
            "uncertain"
        } else if candidates.len() == 1 {
            "one_candidate"
        } else {
            "multiple_candidates"
        };
        if !candidates.is_empty() {
            write_candidates(
                &candidate_dir.join(format!("{family}.fasta")),
                &sample,
                &family,
                &candidates,
                &mut evidence,
            )?;
        }
        writeln!(
            calls,
            "{sample}\t{family}\t{state}\t{}\t{}",
            candidates.len(),
            if flags.is_empty() { "pass" } else { &flags }
        )?;
    }
    Ok(())
}

fn parse_calls(path: &Path) -> io::Result<BTreeMap<String, Call>> {
    let reader = BufReader::new(File::open(path)?);
    let mut calls = BTreeMap::new();
    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        if line_number == 0 {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 5 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid family_calls.tsv",
            ));
        }
        fields[3]
            .parse::<usize>()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "invalid candidate_count"))?;
        calls.insert(
            fields[1].to_owned(),
            Call {
                state: match fields[2] {
                    "one_candidate" => "one_candidate",
                    "multiple_candidates" => "multiple_candidates",
                    "uncertain" => "uncertain",
                    "not_recovered" => "not_recovered",
                    _ => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "invalid call state",
                        ))
                    }
                },
            },
        );
    }
    Ok(calls)
}

fn copy_file_contents(input: &Path, output: &mut BufWriter<File>) -> io::Result<()> {
    if !input.is_file() {
        return Ok(());
    }
    for line in BufReader::new(File::open(input)?).lines() {
        writeln!(output, "{}", line?)?;
    }
    Ok(())
}

fn cohort(reference: PathBuf, out: PathBuf, samples: Vec<String>) -> io::Result<()> {
    if samples.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "cohort needs at least one sample",
        ));
    }
    let families = family_lengths(&reference)?;
    let mut by_sample = BTreeMap::new();
    for sample in &samples {
        by_sample.insert(
            sample.clone(),
            parse_calls(&out.join("samples").join(sample).join("family_calls.tsv"))?,
        );
    }
    let mut summary = BufWriter::new(File::create(out.join("family_summary.tsv"))?);
    writeln!(summary, "family_id\tfamily_state\tone_candidate_samples\tmultiple_candidate_samples\tuncertain_samples\tnot_recovered_samples")?;
    let mut matrix = BufWriter::new(File::create(out.join("family_count_matrix.tsv"))?);
    write!(matrix, "family_id")?;
    for sample in &samples {
        write!(matrix, "\t{sample}")?;
    }
    writeln!(matrix)?;

    let pseudo_dir = out.join("pseudo_sco");
    let multi_dir = out.join("multiple_candidate_families");
    fs::create_dir_all(&pseudo_dir)?;
    fs::create_dir_all(&multi_dir)?;
    let required_occupancy = (samples.len() * 7).div_ceil(10).max(1);

    for family in families.keys() {
        let mut one = 0usize;
        let mut multi = 0usize;
        let mut uncertain = 0usize;
        let mut missing = 0usize;
        write!(matrix, "{family}")?;
        for sample in &samples {
            let call = by_sample
                .get(sample)
                .and_then(|calls| calls.get(family))
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing family call"))?;
            match call.state {
                "one_candidate" => {
                    one += 1;
                    write!(matrix, "\t1")?;
                }
                "multiple_candidates" => {
                    multi += 1;
                    write!(matrix, "\t2+")?;
                }
                "uncertain" => {
                    uncertain += 1;
                    write!(matrix, "\tNA")?;
                }
                _ => {
                    missing += 1;
                    write!(matrix, "\tNA")?;
                }
            }
        }
        writeln!(matrix)?;
        let state = if multi > 0 {
            "multiple_candidate_family"
        } else if one >= required_occupancy {
            "single_candidate_family"
        } else {
            "insufficient"
        };
        writeln!(
            summary,
            "{family}\t{state}\t{one}\t{multi}\t{uncertain}\t{missing}"
        )?;

        if state == "single_candidate_family" || state == "multiple_candidate_family" {
            let target = if state == "single_candidate_family" {
                pseudo_dir.join(format!("{family}.fasta"))
            } else {
                multi_dir.join(format!("{family}.fasta"))
            };
            let mut output = BufWriter::new(File::create(target)?);
            for sample in &samples {
                let call = by_sample[sample].get(family).expect("validated call table");
                if matches!(call.state, "one_candidate" | "multiple_candidates") {
                    copy_file_contents(
                        &out.join("samples")
                            .join(sample)
                            .join("candidates")
                            .join(format!("{family}.fasta")),
                        &mut output,
                    )?;
                }
            }
        }
    }
    Ok(())
}

fn prospective_canonical_path(path: &Path) -> io::Result<PathBuf> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        env::current_dir()?.join(path)
    };
    let mut normalized = PathBuf::new();
    for component in absolute.components() {
        match component {
            std::path::Component::CurDir => {}
            std::path::Component::ParentDir => {
                normalized.pop();
            }
            _ => normalized.push(component.as_os_str()),
        }
    }
    let mut existing = normalized.clone();
    let mut suffix = Vec::new();
    while !existing.exists() {
        let Some(name) = existing.file_name() else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "path has no existing ancestor",
            ));
        };
        suffix.push(name.to_os_string());
        let Some(parent) = existing.parent() else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "path has no existing ancestor",
            ));
        };
        existing = parent.to_path_buf();
    }
    let mut resolved = fs::canonicalize(existing)?;
    for name in suffix.iter().rev() {
        resolved.push(name);
    }
    Ok(resolved)
}

fn ensure_nonoverlapping_paths(inputs: &[&Path], out: &Path) -> io::Result<()> {
    let out = prospective_canonical_path(out)?;
    for input in inputs {
        let input = prospective_canonical_path(input)?;
        if input.starts_with(&out) || out.starts_with(&input) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--out must be separate from every input path",
            ));
        }
    }
    Ok(())
}

fn annotate(
    input: PathBuf,
    proteins: PathBuf,
    out: PathBuf,
    miniprot: String,
    threads: String,
) -> io::Result<()> {
    ensure_nonoverlapping_paths(&[input.as_path(), proteins.as_path()], &out)?;
    if out.exists() {
        fs::remove_dir_all(&out)?;
    }
    let manifest_dir = out.join("manifest");
    let gff_dir = out.join("gff");
    let cds_dir = out.join("cds");
    let intron_dir = out.join("introns");
    let super_dir = out.join("supercontigs");
    fs::create_dir_all(&manifest_dir)?;
    fs::create_dir_all(&gff_dir)?;
    fs::create_dir_all(&cds_dir)?;
    fs::create_dir_all(&intron_dir)?;
    fs::create_dir_all(&super_dir)?;
    let mut manifest = BufWriter::new(File::create(manifest_dir.join("candidate_manifest.tsv"))?);
    writeln!(
        manifest,
        "sample\tfamily_id\tcandidate\tprotein_reference\tminiprot_score\tprotein_identity\tprotein_coverage\texon_count\tcds_length\tintron_length\tsupercontig_length\tstructure_state"
    )?;
    for family_entry in fs::read_dir(input.join("samples"))? {
        let sample_path = family_entry?.path();
        if !sample_path.is_dir() {
            continue;
        };
        let sample = sample_path
            .file_name()
            .unwrap()
            .to_string_lossy()
            .to_string();
        let candidates = sample_path.join("candidates");
        if !candidates.is_dir() {
            continue;
        };
        for entry in fs::read_dir(candidates)? {
            let fasta = entry?.path();
            let Some(family) = family_id(&fasta) else {
                continue;
            };
            let protein = proteins.join(format!("{family}.faa"));
            let records = read_fasta(&fasta)?;
            if !protein.is_file() {
                for (header, _) in &records {
                    let candidate = header.split('|').next_back().unwrap_or("candidate");
                    writeln!(
                        manifest,
                        "{sample}\t{family}\t{candidate}\t\t0\t0\t0\t0\t0\t0\t0\tmissing_protein_reference"
                    )?;
                }
                continue;
            }
            let protein_lengths: BTreeMap<String, usize> = read_raw_fasta(&protein)?
                .into_iter()
                .map(|(header, sequence)| {
                    (
                        header
                            .split_whitespace()
                            .next()
                            .unwrap_or(&header)
                            .to_owned(),
                        sequence.len(),
                    )
                })
                .collect();
            for (header, seq) in records {
                let candidate = header
                    .split('|')
                    .next_back()
                    .unwrap_or("candidate")
                    .to_string();
                let work = out.join(".work");
                fs::create_dir_all(&work)?;
                let dna = work.join(format!("{sample}_{family}_{candidate}.fa"));
                {
                    let mut w = BufWriter::new(File::create(&dna)?);
                    writeln!(w, ">{candidate}\n{seq}")?;
                }
                let result = Command::new(&miniprot)
                    .args(["--gff-only", "-t", &threads])
                    .arg(&dna)
                    .arg(&protein)
                    .output()?;
                if !result.status.success() {
                    writeln!(
                        manifest,
                        "{sample}\t{family}\t{candidate}\t{family}\t0\t0\t0\t0\t0\t0\t0\tminiprot_failed"
                    )?;
                    continue;
                }
                let gff = String::from_utf8_lossy(&result.stdout);
                let mut by_parent: BTreeMap<String, Vec<(usize, usize, &str)>> = BTreeMap::new();
                let mut hit_meta: BTreeMap<String, (String, String, String, String)> =
                    BTreeMap::new();
                for line in gff.lines() {
                    let f: Vec<_> = line.split('\t').collect();
                    if f.len() == 9 && f[2] == "mRNA" {
                        let id = f[8]
                            .split(';')
                            .find_map(|x| x.strip_prefix("ID="))
                            .unwrap_or("")
                            .to_string();
                        let identity = f[8]
                            .split(';')
                            .find_map(|x| x.strip_prefix("Identity="))
                            .unwrap_or("0")
                            .to_string();
                        let target_fields = f[8]
                            .split(';')
                            .find_map(|x| x.strip_prefix("Target="))
                            .unwrap_or("")
                            .split_whitespace()
                            .collect::<Vec<_>>();
                        let target = target_fields.first().copied().unwrap_or("").to_string();
                        let coverage = if let (Some(start), Some(end), Some(length)) = (
                            target_fields.get(1).and_then(|x| x.parse::<usize>().ok()),
                            target_fields.get(2).and_then(|x| x.parse::<usize>().ok()),
                            protein_lengths.get(&target),
                        ) {
                            let denominator = if end > *length {
                                length.saturating_mul(3)
                            } else {
                                *length
                            };
                            format!(
                                "{:.4}",
                                ((end.saturating_sub(start) + 1) as f64 / denominator as f64)
                                    .min(1.0)
                            )
                        } else {
                            "0".to_string()
                        };
                        hit_meta.insert(id, (target, f[5].to_string(), identity, coverage));
                    }
                    if f.len() == 9 && f[2] == "CDS" {
                        if let (Ok(a), Ok(b)) = (f[3].parse::<usize>(), f[4].parse::<usize>()) {
                            let parent = f[8]
                                .split(';')
                                .find_map(|x| x.strip_prefix("Parent="))
                                .unwrap_or("ungrouped")
                                .to_string();
                            if a >= 1 && a <= b && b <= seq.len() {
                                by_parent.entry(parent).or_default().push((a, b, f[6]));
                            }
                        }
                    }
                }
                let (best_parent, mut exons) = by_parent
                    .into_iter()
                    .max_by_key(|(_, v)| v.iter().map(|x| x.1 - x.0 + 1).sum::<usize>())
                    .unwrap_or_else(|| (String::new(), Vec::new()));
                let (protein_id, miniprot_score, protein_identity, protein_coverage) = hit_meta
                    .remove(&best_parent)
                    .unwrap_or_else(|| (String::new(), "0".into(), "0".into(), "0".into()));
                exons.sort_by_key(|x| x.0);
                if exons
                    .iter()
                    .any(|x| x.2 != exons.first().map(|y| y.2).unwrap_or("+"))
                {
                    exons.clear();
                }
                let reverse = exons.first().map(|x| x.2 == "-").unwrap_or(false);
                if reverse {
                    exons.reverse();
                }
                let mut cds = String::new();
                let mut introns = String::new();
                for (i, (a, b, strand)) in exons.iter().enumerate() {
                    let part = &seq[a - 1..*b];
                    let oriented = if *strand == "-" {
                        reverse_complement(part)
                    } else {
                        part.to_string()
                    };
                    cds.push_str(&oriented);
                    if i + 1 < exons.len() {
                        let (na, nb, _) = exons[i + 1];
                        let (l, r) = if reverse { (nb, *a) } else { (*b, na) };
                        if l < r {
                            let part = &seq[l..r - 1];
                            let oriented = if reverse {
                                reverse_complement(part)
                            } else {
                                part.to_string()
                            };
                            introns.push_str(&oriented);
                        }
                    }
                }
                let supercontig = if exons.is_empty() {
                    String::new()
                } else {
                    let lo = exons.iter().map(|x| x.0).min().unwrap();
                    let hi = exons.iter().map(|x| x.1).max().unwrap();
                    let part = &seq[lo - 1..hi];
                    if reverse {
                        reverse_complement(part)
                    } else {
                        part.to_string()
                    }
                };
                let state = if exons.is_empty() {
                    "protein_unsupported"
                } else {
                    "protein_supported"
                };
                let gf = gff_dir.join(&sample).join(&family);
                fs::create_dir_all(&gf)?;
                fs::write(gf.join(format!("{candidate}.gff3")), gff.as_bytes())?;
                if !cds.is_empty() {
                    let mut w = BufWriter::new(
                        File::options()
                            .create(true)
                            .append(true)
                            .open(cds_dir.join(format!("{family}.fasta")))?,
                    );
                    writeln!(w, ">{sample}|{family}|{candidate}\n{cds}")?;
                    let mut w = BufWriter::new(
                        File::options()
                            .create(true)
                            .append(true)
                            .open(intron_dir.join(format!("{family}.fasta")))?,
                    );
                    if !introns.is_empty() {
                        writeln!(w, ">{sample}|{family}|{candidate}\n{introns}")?;
                    }
                    let mut w = BufWriter::new(
                        File::options()
                            .create(true)
                            .append(true)
                            .open(super_dir.join(format!("{family}.fasta")))?,
                    );
                    writeln!(w, ">{sample}|{family}|{candidate}\n{supercontig}")?;
                }
                writeln!(
                    manifest,
                    "{sample}\t{family}\t{candidate}\t{protein_id}\t{miniprot_score}\t{protein_identity}\t{protein_coverage}\t{}\t{}\t{}\t{}\t{state}",
                    exons.len(),
                    cds.len(),
                    introns.len(),
                    supercontig.len()
                )?;
            }
        }
    }
    let mut multi = BufWriter::new(File::create(
        manifest_dir.join("long_multiple_candidates.tsv"),
    )?);
    writeln!(multi, "sample\tfamily_id\tprotein_supported_candidates")?;
    for entry in fs::read_dir(&cds_dir)? {
        let fasta = entry?.path();
        let Some(family) = family_id(&fasta) else {
            continue;
        };
        let mut counts: BTreeMap<String, usize> = BTreeMap::new();
        for (header, _) in read_fasta(&fasta)? {
            if let Some(sample) = header.split('|').next() {
                *counts.entry(sample.to_owned()).or_default() += 1
            }
        }
        for (sample, count) in counts {
            if count > 1 {
                writeln!(multi, "{sample}\t{family}\t{count}")?;
            }
        }
    }
    Ok(())
}

fn alignment_qc(records: &[(String, String)]) -> (usize, usize, f64) {
    let columns = records
        .first()
        .map(|(_, sequence)| sequence.len())
        .unwrap_or(0);
    if columns == 0
        || records
            .iter()
            .any(|(_, sequence)| sequence.len() != columns)
    {
        return (columns, 0, 0.0);
    }
    let effective = (0..columns)
        .filter(|&column| {
            records
                .iter()
                .filter(|(_, sequence)| {
                    sequence.as_bytes()[column].is_ascii_alphabetic()
                        && sequence.as_bytes()[column] != b'X'
                })
                .count()
                * 2
                >= records.len()
        })
        .count();
    (columns, effective, effective as f64 / columns as f64)
}

/// Summarize a codon alignment in codon, rather than nucleotide, columns.
/// A codon is effective when at least half of its sequences contain an
/// unmasked, complete codon at that site. Keeping this unit explicit makes
/// the value in family_qc.tsv directly comparable to the AA alignment.
fn codon_alignment_qc(records: &[(String, String)]) -> (usize, usize, f64) {
    let bases = records
        .first()
        .map(|(_, sequence)| sequence.len())
        .unwrap_or(0);
    if bases == 0
        || !bases.is_multiple_of(3)
        || records.iter().any(|(_, sequence)| sequence.len() != bases)
    {
        return (bases / 3, 0, 0.0);
    }
    let codons = bases / 3;
    let effective = (0..codons)
        .filter(|&column| {
            let start = column * 3;
            records
                .iter()
                .filter(|(_, sequence)| {
                    sequence.as_bytes()[start..start + 3]
                        .iter()
                        .all(|base| base.is_ascii_alphabetic() && *base != b'X')
                })
                .count()
                * 2
                >= records.len()
        })
        .count();
    (codons, effective, effective as f64 / codons as f64)
}

fn clade_support(tree: &resolve::Tree, node: usize, parent: Option<usize>) -> String {
    let support_node = match parent {
        Some(parent) if tree.nodes[parent].children.contains(&node) => node,
        Some(parent) if tree.nodes[node].children.contains(&parent) => parent,
        _ => node,
    };
    tree.nodes[support_node]
        .name
        .as_deref()
        .and_then(|label| label.parse::<f64>().ok())
        .filter(|value| value.is_finite() && (0.0..=100.0).contains(value))
        .map(|value| format!("{value:.3}"))
        .unwrap_or_else(|| "NA".to_owned())
}

fn tree_copy_counts(tree: &resolve::Tree, outgroups: &BTreeSet<String>) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for header in resolve::all_leaf_names(tree) {
        if let Some(sample) = header.split('|').next() {
            if !outgroups.contains(sample) {
                *counts.entry(sample.to_owned()).or_default() += 1;
            }
        }
    }
    counts
}

fn distinct_sample_count(records: &[(String, String)]) -> usize {
    records
        .iter()
        .filter_map(|(header, _)| header.split('|').next())
        .collect::<BTreeSet<_>>()
        .len()
}

fn median_sequence_length(records: &[(String, String)]) -> usize {
    let mut lengths: Vec<_> = records.iter().map(|(_, sequence)| sequence.len()).collect();
    lengths.sort_unstable();
    lengths.get(lengths.len() / 2).copied().unwrap_or(0)
}

#[allow(clippy::too_many_arguments)]
fn resolve_workflow(
    input: PathBuf,
    out: PathBuf,
    mafft: String,
    iqtree: String,
    threads: String,
    min_taxa: usize,
    min_aa_length: usize,
    min_effective_codon_sites: usize,
    outgroup: Option<PathBuf>,
    ufboot: usize,
    taper_script: Option<PathBuf>,
    julia: String,
) -> io::Result<()> {
    let mut input_paths = vec![input.as_path()];
    if let Some(path) = outgroup.as_deref() {
        input_paths.push(path);
    }
    ensure_nonoverlapping_paths(&input_paths, &out)?;
    if out.exists() {
        fs::remove_dir_all(&out)?;
    }
    let outgroups: BTreeSet<String> = match outgroup {
        Some(path) => {
            let mut values = BTreeSet::new();
            for line in BufReader::new(File::open(path)?)
                .lines()
                .map_while(Result::ok)
            {
                let value = line.split(['\t', ',']).next().unwrap_or("").trim();
                if value.is_empty()
                    || value.starts_with('#')
                    || matches!(
                        value.to_ascii_lowercase().as_str(),
                        "sample" | "species" | "taxon" | "taxa"
                    )
                {
                    continue;
                }
                values.insert(value.to_owned());
            }
            if values.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--outgroup contains no sample identifiers",
                ));
            }
            values
        }
        None => BTreeSet::new(),
    };
    let cds = input.join("cds");
    let work = out.join("work");
    let strict = out.join("resolved_1to1");
    let unresolved = out.join("unresolved_multicandidate");
    fs::create_dir_all(&work)?;
    fs::create_dir_all(&strict)?;
    fs::create_dir_all(&unresolved)?;
    let mut manifest = BufWriter::new(File::create(out.join("resolve_manifest.tsv"))?);
    let astral_dir = out.join("astral_input");
    let pro_dir = out.join("astralpro_input");
    fs::create_dir_all(&astral_dir)?;
    fs::create_dir_all(&pro_dir)?;
    let mut astral = BufWriter::new(File::create(astral_dir.join("resolved_1to1.trees"))?);
    let mut leaf_map = BufWriter::new(File::create(pro_dir.join("leaf_to_species.tsv"))?);
    let mut pro_trees = BufWriter::new(File::create(pro_dir.join("multicopy.trees"))?);
    let mut family_qc = BufWriter::new(File::create(out.join("family_qc.tsv"))?);
    let mut occupancy_qc = BufWriter::new(File::create(out.join("occupancy_qc.tsv"))?);
    let mut selection_qc = BufWriter::new(File::create(out.join("tree_selection_qc.tsv"))?);
    writeln!(manifest, "family_id\tstatus\tclade\ttaxa\treason")?;
    writeln!(family_qc, "family_id\tstatus\tinput_candidates\taa_alignment_columns\tcodon_alignment_columns\teffective_codon_columns\teffective_codon_fraction\ttaper_applied")?;
    writeln!(occupancy_qc, "family_id\tstage\tinput_candidates\tretained_candidates\tdistinct_samples\tmedian_sequence_length\tminimum_required\tstatus\treason")?;
    writeln!(selection_qc, "family_id\ttree_samples\tsingle_candidate_samples\tmulti_candidate_samples\tselected_clade\tclade_taxa\tclade_occupancy\tclade_support\tselected_leaves")?;
    for entry in fs::read_dir(&cds)? {
        let fasta = entry?.path();
        let Some(family) = family_id(&fasta) else {
            continue;
        };
        let records = read_fasta(&fasta)?;
        let family_work = work.join(&family);
        fs::create_dir_all(&family_work)?;
        let aa_input = family_work.join("proteins.fasta");
        let mut aa_records = Vec::new();
        let mut invalid_translation = 0usize;
        let mut short_protein = 0usize;
        for (header, cds) in &records {
            let protein = translate_cds(cds);
            if protein.contains('X') || protein.contains('*') {
                invalid_translation += 1;
            } else if protein.len() < min_aa_length {
                short_protein += 1;
            } else {
                aa_records.push((header.clone(), protein));
            }
        }
        let pre_samples = distinct_sample_count(&aa_records);
        let pre_median_length = median_sequence_length(&aa_records);
        let pre_reason = if pre_samples < min_taxa {
            "too_few_pre_alignment_taxa"
        } else if aa_records.is_empty() {
            "no_pre_alignment_proteins"
        } else {
            "pass"
        };
        writeln!(
            occupancy_qc,
            "{family}\tpre_alignment\t{}\t{}\t{pre_samples}\t{pre_median_length}\t{min_taxa}\t{}\t{pre_reason}",
            records.len(),
            aa_records.len(),
            if pre_reason == "pass" { "pass" } else { "fail" },
        )?;
        if pre_reason != "pass" {
            let reason = if invalid_translation > 0 || short_protein > 0 {
                format!("{pre_reason};invalid_translation={invalid_translation};short_protein={short_protein}")
            } else {
                pre_reason.to_owned()
            };
            writeln!(manifest, "{family}\tunresolved\t\t\t{reason}")?;
            fs::copy(&fasta, unresolved.join(format!("{family}.fasta")))?;
            continue;
        }
        {
            let mut writer = BufWriter::new(File::create(&aa_input)?);
            for (header, protein) in &aa_records {
                writeln!(writer, ">{header}\n{protein}")?;
            }
        }
        let aa_aln = family_work.join("aligned.aa.fasta");
        let output = Command::new(&mafft)
            .args(["--auto", "--thread", &threads])
            .arg(&aa_input)
            .output()?;
        if !output.status.success() {
            writeln!(manifest, "{family}\tunresolved\t\t\tmafft_failed")?;
            fs::copy(&fasta, unresolved.join(format!("{family}.fasta")))?;
            continue;
        }
        fs::write(&aa_aln, &output.stdout)?;
        let aa_for_backtranslation = if let Some(script) = &taper_script {
            let tapered = family_work.join("tapered.aa.fasta");
            let status = Command::new(&julia)
                .arg(script)
                .arg(&aa_aln)
                .stdout(File::create(&tapered)?)
                .status();
            match status {
                Ok(status)
                    if status.success()
                        && tapered.is_file()
                        && fs::metadata(&tapered)?.len() > 0 =>
                {
                    tapered
                }
                _ => {
                    writeln!(manifest, "{family}\tunresolved\t\t\ttaper_failed")?;
                    fs::copy(&fasta, unresolved.join(format!("{family}.fasta")))?;
                    continue;
                }
            }
        } else {
            aa_aln.clone()
        };
        let aligned_records = read_raw_fasta(&aa_for_backtranslation)?;
        let expected_headers: BTreeSet<String> = aa_records
            .iter()
            .map(|(header, _)| header.clone())
            .collect();
        let observed_headers: BTreeSet<String> = aligned_records
            .iter()
            .map(|(header, _)| header.clone())
            .collect();
        let output_reason = if taper_script.is_some() {
            "taper_output_mismatch"
        } else {
            "mafft_output_mismatch"
        };
        if aligned_records.len() != expected_headers.len()
            || observed_headers.len() != aligned_records.len()
            || observed_headers != expected_headers
            || aligned_records
                .iter()
                .any(|(_, sequence)| sequence.is_empty())
        {
            writeln!(manifest, "{family}\tunresolved\t\t\t{output_reason}")?;
            fs::copy(&fasta, unresolved.join(format!("{family}.fasta")))?;
            continue;
        }
        let cds_by_header: BTreeMap<String, String> = records
            .iter()
            .map(|(h, q)| (h.clone(), q.clone()))
            .collect();
        let aln = family_work.join("aligned.codon.fasta");
        let codon_alignment: Vec<(String, String)> = aligned_records
            .iter()
            .filter_map(|(header, aligned)| {
                cds_by_header
                    .get(header)
                    .and_then(|cds| codon_backtranslate(aligned, cds))
                    .map(|codon| (header.clone(), codon))
            })
            .collect();
        {
            let mut writer = BufWriter::new(File::create(&aln)?);
            for (header, codon) in &codon_alignment {
                writeln!(writer, ">{header}\n{codon}")?;
            }
        }
        let aa_qc = alignment_qc(&aligned_records);
        let codon_qc = codon_alignment_qc(&codon_alignment);
        let post_samples = distinct_sample_count(&codon_alignment);
        let post_median_length = median_sequence_length(&codon_alignment) / 3;
        let post_reason = if post_samples < min_taxa {
            "too_few_post_alignment_taxa"
        } else if codon_qc.1 < min_effective_codon_sites {
            "too_few_effective_codon_sites"
        } else {
            "pass"
        };
        writeln!(
            occupancy_qc,
            "{family}\tpost_alignment\t{}\t{}\t{post_samples}\t{post_median_length}\t{min_effective_codon_sites}\t{}\t{post_reason}",
            aa_records.len(),
            codon_alignment.len(),
            if post_reason == "pass" { "pass" } else { "fail" },
        )?;
        if post_reason != "pass" {
            writeln!(manifest, "{family}\tunresolved\t\t\t{post_reason}")?;
            fs::copy(&fasta, unresolved.join(format!("{family}.fasta")))?;
            continue;
        }
        writeln!(
            family_qc,
            "{family}\talignment_pass\t{}\t{}\t{}\t{}\t{:.6}\t{}",
            records.len(),
            aa_qc.0,
            codon_qc.0,
            codon_qc.1,
            codon_qc.2,
            taper_script.is_some()
        )?;
        let prefix = family_work.join("tree");
        let mut iqtree_command = Command::new(&iqtree);
        iqtree_command.args([
            "-s",
            aln.to_str().unwrap(),
            "-m",
            "MFP",
            "-T",
            &threads,
            "--seed",
            "1",
            "--prefix",
            prefix.to_str().unwrap(),
            "-redo",
        ]);
        let ufboot_value = ufboot.to_string();
        if ufboot > 0 {
            iqtree_command.args(["-B", &ufboot_value]);
        }
        let status = iqtree_command.status()?;
        let tree_path = family_work.join("tree.treefile");
        if !status.success() || !tree_path.is_file() {
            writeln!(manifest, "{family}\tunresolved\t\t\tiqtree_failed")?;
            fs::copy(&fasta, unresolved.join(format!("{family}.fasta")))?;
            continue;
        }
        let tree = match resolve::parse_newick(&fs::read_to_string(&tree_path)?) {
            Ok(x) => x,
            Err(e) => {
                writeln!(manifest, "{family}\tunresolved\t\t\tnewick:{e}")?;
                fs::copy(&fasta, unresolved.join(format!("{family}.fasta")))?;
                continue;
            }
        };
        writeln!(pro_trees, "{}", fs::read_to_string(&tree_path)?.trim())?;
        for header in resolve::all_leaf_names(&tree) {
            if let Some(species) = header.split('|').next() {
                writeln!(leaf_map, "{header}\t{species}")?;
            }
        }
        let clades = match resolve::select_scogs(&tree, min_taxa, &outgroups) {
            Ok(clades) => clades,
            Err(reason) => {
                writeln!(manifest, "{family}\tunresolved\t\t\t{reason}")?;
                fs::copy(&fasta, unresolved.join(format!("{family}.fasta")))?;
                continue;
            }
        };
        let copy_counts = tree_copy_counts(&tree, &outgroups);
        let tree_samples = copy_counts.len();
        let single_candidate_samples = copy_counts.values().filter(|&&count| count == 1).count();
        let multi_candidate_samples = copy_counts.values().filter(|&&count| count > 1).count();
        if clades.is_empty() {
            writeln!(
                selection_qc,
                "{family}\t{tree_samples}\t{single_candidate_samples}\t{multi_candidate_samples}\tNA\t0\t0.000000\tNA\t0"
            )?;
            writeln!(manifest, "{family}\tunresolved\t\t\tno_one_to_one_clade")?;
            fs::copy(&fasta, unresolved.join(format!("{family}.fasta")))?;
            continue;
        }
        let mut used = HashSet::new();
        for (i, c) in clades.iter().enumerate() {
            let name = format!("{family}.og{}", i + 1);
            let mut w = BufWriter::new(File::create(strict.join(format!("{name}.fasta")))?);
            for (h, q) in &records {
                if c.leaves.iter().any(|x| x == h) {
                    writeln!(w, ">{h}\n{q}")?;
                    used.insert(h.clone());
                }
            }
            let star = format!("{};", resolve::render_clade(&tree, c.node, c.parent));
            let astral_tree = format!(
                "{};",
                resolve::render_clade_samples(&tree, c.node, c.parent)
            );
            fs::write(strict.join(format!("{name}.treefile")), &star)?;
            writeln!(astral, "{astral_tree}")?;
            writeln!(
                manifest,
                "{family}\tresolved\t{name}\t{}\tpass",
                c.samples.len()
            )?;
            writeln!(
                selection_qc,
                "{family}\t{tree_samples}\t{single_candidate_samples}\t{multi_candidate_samples}\t{name}\t{}\t{:.6}\t{}\t{}",
                c.samples.len(),
                c.samples.len() as f64 / tree_samples.max(1) as f64,
                clade_support(&tree, c.node, c.parent),
                c.leaves.len()
            )?;
        }
        if used.len() < records.len() {
            let mut w = BufWriter::new(File::create(unresolved.join(format!("{family}.fasta")))?);
            for (h, q) in records {
                if !used.contains(&h) {
                    writeln!(w, ">{h}\n{q}")?
                }
            }
            writeln!(manifest, "{family}\tunresolved\t\t\tremaining_candidates")?;
        }
    }
    Ok(())
}

fn main() -> io::Result<()> {
    let raw: Vec<String> = env::args().skip(1).collect();
    let Some((command, options)) = raw.split_first() else {
        usage()
    };
    if command == "--help" || command == "-h" {
        usage();
    }
    let options = parse_options(options);
    match command.as_str() {
        "classify" => classify(
            option_path(&options, "--reference"),
            option_path(&options, "--contigs"),
            option_string(&options, "--sample"),
            option_path(&options, "--out"),
        ),
        "cohort" => cohort(
            option_path(&options, "--reference"),
            option_path(&options, "--out"),
            options.get("--sample").cloned().unwrap_or_default(),
        ),
        "resolve" => resolve_workflow(
            option_path(&options, "--input"),
            option_path(&options, "--out"),
            option_string(&options, "--mafft"),
            option_string(&options, "--iqtree"),
            options
                .get("--threads")
                .and_then(|v| v.first())
                .cloned()
                .unwrap_or_else(|| "1".into()),
            {
                let value = option_string(&options, "--min-taxa");
                match value.parse::<usize>() {
                    Ok(n) if n >= 2 => n,
                    _ => {
                        eprintln!("--min-taxa must be an integer of at least 2");
                        usage();
                    }
                }
            },
            option_positive_usize(&options, "--min-aa-length", 30),
            option_positive_usize(&options, "--min-effective-codon-sites", 30),
            options
                .get("--outgroup")
                .and_then(|v| v.first())
                .map(PathBuf::from),
            {
                let value = options
                    .get("--ufboot")
                    .and_then(|values| values.first())
                    .map(String::as_str)
                    .unwrap_or("0");
                match value.parse::<usize>() {
                    Ok(0) | Ok(1000..) => value.parse::<usize>().unwrap(),
                    _ => {
                        eprintln!("--ufboot must be 0 or an integer of at least 1000");
                        usage();
                    }
                }
            },
            options
                .get("--taper-script")
                .and_then(|v| v.first())
                .map(PathBuf::from),
            options
                .get("--julia")
                .and_then(|v| v.first())
                .cloned()
                .unwrap_or_else(|| "julia".into()),
        ),
        "annotate" => annotate(
            option_path(&options, "--input"),
            option_path(&options, "--protein-reference"),
            option_path(&options, "--out"),
            option_string(&options, "--miniprot"),
            options
                .get("--threads")
                .and_then(|v| v.first())
                .cloned()
                .unwrap_or_else(|| "1".into()),
        ),
        _ => usage(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_collapses_reverse_complements() {
        assert_eq!(canonical("AACG".into()), canonical("CGTT".into()));
    }

    #[test]
    fn retains_original_rust_assembly_metrics() {
        assert_eq!(
            assembly_metrics("contig_120_7_42_99_31"),
            ["7", "42", "99", "31"].map(str::to_owned)
        );
        assert_eq!(
            assembly_metrics("unstructured"),
            [String::new(), String::new(), String::new(), String::new()]
        );
    }

    #[test]
    fn captures_numeric_internal_clade_support() {
        let tree = resolve::parse_newick("(A|OG|candidate_1,B|OG|candidate_1)97.5;").unwrap();
        assert_eq!(clade_support(&tree, tree.root, None), "97.500");
    }

    #[test]
    fn support_comes_from_the_original_child_on_a_complement_clade() {
        let tree =
            resolve::parse_newick("(A|OG|candidate_1,(B|OG|candidate_1,C|OG|candidate_1)88.2);")
                .unwrap();
        let child = tree.nodes[tree.root].children[1];
        assert_eq!(clade_support(&tree, tree.root, Some(child)), "88.200");
    }

    #[test]
    fn occupancy_counts_distinct_samples_not_candidate_records() {
        let records = vec![
            ("A|family|candidate_1".into(), "A".repeat(30)),
            ("A|family|candidate_2".into(), "A".repeat(60)),
            ("B|family|candidate_1".into(), "A".repeat(90)),
        ];
        assert_eq!(distinct_sample_count(&records), 2);
        assert_eq!(median_sequence_length(&records), 60);
    }

    #[test]
    fn codon_qc_counts_triplets_and_excludes_masked_codons() {
        let records = vec![
            ("A".into(), "AAACCCGGG".into()),
            ("B".into(), "AAA---GGG".into()),
            ("C".into(), "AAAXXXGGG".into()),
        ];
        assert_eq!(codon_alignment_qc(&records), (3, 2, 2.0 / 3.0));
    }

    #[test]
    fn only_strict_prefixes_are_collapsed() {
        let prefix = Candidate {
            sequence: "AAAT".into(),
            key: "AAAT".into(),
            assembly_metrics: Default::default(),
        };
        let extension = Candidate {
            sequence: "AAATGG".into(),
            key: "AAATGG".into(),
            assembly_metrics: Default::default(),
        };
        let internal_repeat = Candidate {
            sequence: "CCCAAATTT".into(),
            key: "CCCAAATTT".into(),
            assembly_metrics: Default::default(),
        };
        assert!(is_strict_path_prefix(&extension, &prefix));
        assert!(!is_strict_path_prefix(&internal_repeat, &prefix));
    }
}
