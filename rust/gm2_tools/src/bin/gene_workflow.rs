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
use std::process::{Command, Output};

use gm2_tools::gene_annotation::{
    build_gene_model, read_dna_fasta, reverse_complement as reverse_complement_iupac,
    select_competing_models, GeneModel, Interval, ModelConfig, ModelState, Strand,
};

const FASTA_EXTENSIONS: &[&str] = &["fa", "fas", "fasta"];

fn complementary_union_fraction(
    left: Interval,
    right: Interval,
    protein_len: usize,
) -> Option<f64> {
    if protein_len == 0 {
        return None;
    }
    let overlap = left.overlap(right);
    let shorter = left.len().min(right.len());
    let covered_union = left.len() + right.len() - overlap;
    let improvement = covered_union.saturating_sub(left.len().max(right.len()));
    (overlap * 10 <= shorter && improvement * 5 >= protein_len)
        .then_some(covered_union as f64 / protein_len as f64)
}

fn complete_fragment_union_fraction(
    left: Interval,
    right: Interval,
    protein_len: usize,
    complete_coverage: f64,
) -> Option<f64> {
    let fraction = complementary_union_fraction(left, right, protein_len)?;
    let terminal_tolerance = 3usize.max((protein_len as f64 * 0.02).ceil() as usize);
    let start = left.start.min(right.start);
    let end = left.end.max(right.end);
    (fraction >= complete_coverage
        && start <= terminal_tolerance
        && protein_len.saturating_sub(end) <= terminal_tolerance)
        .then_some(fraction)
}

fn oriented_candidate(sequence: &str, strand: Strand) -> String {
    match strand {
        Strand::Forward => sequence.to_owned(),
        Strand::Reverse => reverse_complement_iupac(sequence),
    }
}

fn run_miniprot(
    executable: &str,
    target: &Path,
    proteins: &Path,
    threads: &str,
    max_intron: usize,
) -> io::Result<(Output, String)> {
    let max_intron = max_intron.to_string();
    let command = format!(
        "{executable} --gff -j 1 -G {max_intron} --outc 0.10 --outs 0.20 -N 32 --outn 32 -t {threads} {} {}",
        target.display(),
        proteins.display()
    );
    let output = Command::new(executable)
        .args([
            "--gff",
            "-j",
            "1",
            "-G",
            &max_intron,
            "--outc",
            "0.10",
            "--outs",
            "0.20",
            "-N",
            "32",
            "--outn",
            "32",
            "-t",
            threads,
        ])
        .arg(target)
        .arg(proteins)
        .output()?;
    Ok((output, command))
}

#[derive(Clone, Debug)]
struct FragmentOutcome {
    left_model: usize,
    right_model: usize,
    covered_union_fraction: f64,
    padding_nt: usize,
    derived_candidate: String,
    status: String,
    reason: String,
}

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
        "Usage:\n  gene_workflow classify --reference DIR --contigs DIR --sample NAME --out DIR\n  gene_workflow cohort --reference DIR --out DIR --sample NAME [--sample NAME ...]\n  gene_workflow annotate --input DIR --protein-reference DIR --out DIR --miniprot FILE [--threads N] [--max-intron N] [--minimum-coverage F] [--complete-coverage F] [--flank N] [--fragment-padding N]\n  gene_workflow resolve --input DIR --out DIR --mafft FILE --iqtree FILE --min-taxa N [--threads N] [--outgroup FILE] [--ufboot N] [--min-aa-length N] [--min-effective-codon-sites N] [--taper-script FILE --julia FILE]"
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

fn option_usize(options: &HashMap<String, Vec<String>>, name: &str, default: usize) -> usize {
    let value = options
        .get(name)
        .and_then(|values| values.first())
        .map(String::as_str)
        .unwrap_or("");
    if value.is_empty() {
        return default;
    }
    value.parse::<usize>().unwrap_or_else(|_| {
        eprintln!("{name} must be a non-negative integer");
        usage();
    })
}

fn option_fraction(options: &HashMap<String, Vec<String>>, name: &str, default: f64) -> f64 {
    let value = options
        .get(name)
        .and_then(|values| values.first())
        .map(String::as_str)
        .unwrap_or("");
    if value.is_empty() {
        return default;
    }
    match value.parse::<f64>() {
        Ok(value) if value.is_finite() && (0.0..=1.0).contains(&value) => value,
        _ => {
            eprintln!("{name} must be a finite number between 0 and 1");
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
    (offset == cds.len()).then_some(out)
}

fn trim_terminal_stop<'a>(cds: &'a str, protein: &mut String) -> &'a str {
    if protein.ends_with('*') {
        protein.pop();
        &cds[..cds.len() - 3]
    } else {
        cds
    }
}

fn resolve_eligible_headers(input: &Path) -> io::Result<Option<BTreeSet<String>>> {
    let manifest = input.join("models/gene_models.tsv");
    if !manifest.is_file() {
        return Ok(None);
    }
    let mut lines = BufReader::new(File::open(manifest)?).lines();
    let Some(header) = lines.next().transpose()? else {
        return Ok(Some(BTreeSet::new()));
    };
    let columns: Vec<_> = header.split('\t').collect();
    let index = |name: &str| {
        columns
            .iter()
            .position(|column| *column == name)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("gene_models.tsv lacks required column {name}"),
                )
            })
    };
    let sample_index = index("sample")?;
    let family_index = index("family_id")?;
    let candidate_index = index("candidate")?;
    let model_index = index("model")?;
    let eligible_index = index("eligible_for_resolve")?;
    let mut eligible = BTreeSet::new();
    for line in lines {
        let fields: Vec<_> = line?.split('\t').map(str::to_owned).collect();
        if fields.get(eligible_index).map(String::as_str) != Some("1") {
            continue;
        }
        let Some(sample) = fields.get(sample_index) else {
            continue;
        };
        let Some(family) = fields.get(family_index) else {
            continue;
        };
        let Some(candidate) = fields.get(candidate_index) else {
            continue;
        };
        let Some(model) = fields.get(model_index) else {
            continue;
        };
        eligible.insert(format!("{sample}|{family}|{candidate}|{model}"));
    }
    Ok(Some(eligible))
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

fn sorted_directory_paths(path: &Path) -> io::Result<Vec<PathBuf>> {
    let mut paths = fs::read_dir(path)?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect::<io::Result<Vec<_>>>()?;
    paths.sort();
    Ok(paths)
}

fn append_fasta(path: &Path, header: &str, sequence: &str) -> io::Result<()> {
    let mut writer = BufWriter::new(File::options().create(true).append(true).open(path)?);
    writeln!(writer, ">{header}\n{sequence}")
}

fn safe_identifier(value: &str) -> String {
    value
        .chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || matches!(character, '.' | '_' | '-') {
                character
            } else {
                '_'
            }
        })
        .collect()
}

fn write_normalized_gff(
    path: &Path,
    sample: &str,
    family: &str,
    model_name: &str,
    model: &GeneModel,
) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    let root = safe_identifier(&format!(
        "{sample}.{family}.{}.{}",
        model.candidate_id, model_name
    ));
    let gene_id = format!("{root}.gene");
    let transcript_id = format!("{root}.mrna");
    let seqid = &model.candidate_id;
    let start = model.genomic_interval.start + 1;
    let end = model.genomic_interval.end;
    let strand = model.strand.as_char();
    writeln!(writer, "##gff-version 3")?;
    writeln!(
        writer,
        "{seqid}\tTipSeek\tgene\t{start}\t{end}\t{}\t{strand}\t.\tID={gene_id};family={}",
        model.alignment_score,
        safe_identifier(family)
    )?;
    writeln!(
        writer,
        "{seqid}\tTipSeek\tmRNA\t{start}\t{end}\t{}\t{strand}\t.\tID={transcript_id};Parent={gene_id};protein_reference={};state={}",
        model.alignment_score,
        safe_identifier(&model.protein_id),
        model.state.as_str()
    )?;
    for exon in &model.exons {
        let exon_start = exon.interval.start + 1;
        let exon_end = exon.interval.end;
        let phase = exon
            .phase
            .map(|value| value.to_string())
            .unwrap_or_else(|| ".".to_owned());
        writeln!(
            writer,
            "{seqid}\tTipSeek\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\tID={transcript_id}.exon{};Parent={transcript_id}",
            exon.index
        )?;
        writeln!(
            writer,
            "{seqid}\tTipSeek\tCDS\t{exon_start}\t{exon_end}\t{}\t{strand}\t{phase}\tID={transcript_id}.cds{};Parent={transcript_id}",
            model.alignment_score,
            exon.index
        )?;
    }
    for intron in &model.introns {
        let intron_start = intron.interval.start + 1;
        let intron_end = intron.interval.end;
        writeln!(
            writer,
            "{seqid}\tTipSeek\tintron\t{intron_start}\t{intron_end}\t.\t{strand}\t.\tID={transcript_id}.intron{};Parent={transcript_id};donor={};acceptor={};splice_class={}",
            intron.index,
            intron.donor,
            intron.acceptor,
            intron.splice_class.as_str()
        )?;
    }
    Ok(())
}

fn manifest_header() -> &'static str {
    "sample\tfamily_id\tcandidate\tmodel\tprotein_reference\tquery_length_aa\tquery_start_aa\tquery_end_aa\tprotein_coverage\talignment_score\tidentity\tpositive\tstrand\tgenomic_start\tgenomic_end\tframeshifts\tinternal_stops\texon_count\tintron_count\tcds_length\tprotein_length\tgene_length\tfive_prime_trim_nt\tthree_prime_trim_nt\tstate\tqc_flags\tcompetition_group\tselected\teligible_for_resolve"
}

fn write_status_row<W: Write>(
    writer: &mut W,
    sample: &str,
    family: &str,
    candidate: &str,
    state: &str,
    qc: &str,
) -> io::Result<()> {
    let mut fields = vec![sample.to_owned(), family.to_owned(), candidate.to_owned()];
    fields.extend((0..21).map(|_| String::new()));
    fields.push(state.to_owned());
    fields.push(qc.to_owned());
    fields.extend((0..3).map(|_| "0".to_owned()));
    writeln!(writer, "{}", fields.join("\t"))
}

fn model_manifest_fields(
    sample: &str,
    family: &str,
    candidate: &str,
    model_name: &str,
    model: &GeneModel,
) -> Vec<String> {
    vec![
        sample.to_owned(),
        family.to_owned(),
        candidate.to_owned(),
        model_name.to_owned(),
        model.protein_id.clone(),
        model.protein_len_aa.to_string(),
        model.query_interval_aa.start.to_string(),
        model.query_interval_aa.end.to_string(),
        format!("{:.6}", model.coverage()),
        model.alignment_score.to_string(),
        format!("{:.6}", model.identity()),
        format!("{:.6}", model.positive_fraction()),
        model.strand.as_char().to_string(),
        (model.genomic_interval.start + 1).to_string(),
        model.genomic_interval.end.to_string(),
        model.frameshifts.to_string(),
        model.inframe_stops.to_string(),
        model.exons.len().to_string(),
        model.introns.len().to_string(),
        model.cds.len().to_string(),
        model
            .protein
            .as_ref()
            .map(String::len)
            .unwrap_or(0)
            .to_string(),
        model.gene.len().to_string(),
        model.five_prime_trim_nt.to_string(),
        model.three_prime_trim_nt.to_string(),
        model.state.as_str().to_owned(),
        model.qc_flags.join(";"),
        model.competition_group.to_string(),
        usize::from(model.selected).to_string(),
        usize::from(model.selected && model.eligible_for_resolve).to_string(),
    ]
}

fn write_unresolved_sequence(
    unresolved: &Path,
    family: &str,
    state: &str,
    header: &str,
    sequence: &str,
) -> io::Result<()> {
    let directory = unresolved.join(state);
    fs::create_dir_all(&directory)?;
    append_fasta(&directory.join(format!("{family}.fasta")), header, sequence)
}

fn validated_miniprot_version(executable: &str) -> io::Result<String> {
    let output = Command::new(executable).arg("--version").output()?;
    if !output.status.success() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "miniprot --version failed",
        ));
    }
    let text = format!(
        "{} {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let version = text
        .split_whitespace()
        .find(|token| {
            token
                .chars()
                .next()
                .is_some_and(|value| value.is_ascii_digit())
        })
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "could not parse miniprot version",
            )
        })?;
    let numeric = version.split('-').next().unwrap_or(version);
    let mut parts = numeric.split('.');
    let major = parts.next().and_then(|value| value.parse::<usize>().ok());
    let minor = parts.next().and_then(|value| value.parse::<usize>().ok());
    if !matches!((major, minor), (Some(major), Some(minor)) if major > 0 || minor >= 18) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("miniprot >=0.18 is required; found {version}"),
        ));
    }
    Ok(version.to_owned())
}

#[allow(clippy::too_many_arguments)]
fn annotate(
    input: PathBuf,
    proteins: PathBuf,
    out: PathBuf,
    miniprot: String,
    threads: String,
    max_intron: usize,
    minimum_coverage: f64,
    complete_coverage: f64,
    flank: usize,
    fragment_padding: usize,
) -> io::Result<()> {
    ensure_nonoverlapping_paths(&[input.as_path(), proteins.as_path()], &out)?;
    if !(0.0..=1.0).contains(&minimum_coverage)
        || !(0.0..=1.0).contains(&complete_coverage)
        || minimum_coverage > complete_coverage
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "gene model coverage thresholds must satisfy 0 <= minimum <= complete <= 1",
        ));
    }
    let miniprot_version = validated_miniprot_version(&miniprot)?;
    if out.exists() {
        fs::remove_dir_all(&out)?;
    }
    let manifest_dir = out.join("manifest");
    let models_dir = out.join("models");
    let raw_dir = out.join("raw_miniprot");
    let gff_dir = out.join("gff3");
    let cds_dir = out.join("cds");
    let protein_dir = out.join("proteins");
    let exon_dir = out.join("exons");
    let intron_dir = out.join("introns");
    let gene_dir = out.join("genes");
    let flanked_dir = out.join("genes_flanked");
    let super_dir = out.join("supercontigs");
    let unresolved_dir = out.join("unresolved");
    let work_dir = out.join(".work");
    for directory in [
        &manifest_dir,
        &models_dir,
        &raw_dir,
        &gff_dir,
        &cds_dir,
        &protein_dir,
        &exon_dir,
        &intron_dir,
        &gene_dir,
        &flanked_dir,
        &super_dir,
        &unresolved_dir,
        &work_dir,
    ] {
        fs::create_dir_all(directory)?;
    }
    let mut candidate_manifest =
        BufWriter::new(File::create(manifest_dir.join("candidate_manifest.tsv"))?);
    let mut gene_models = BufWriter::new(File::create(models_dir.join("gene_models.tsv"))?);
    let mut gene_segments = BufWriter::new(File::create(models_dir.join("gene_segments.tsv"))?);
    let mut id_map = BufWriter::new(File::create(manifest_dir.join("id_map.tsv"))?);
    let mut warnings = BufWriter::new(File::create(manifest_dir.join("annotation_warnings.tsv"))?);
    let mut fragments = BufWriter::new(File::create(manifest_dir.join("fragment_groups.tsv"))?);
    let mut provenance = BufWriter::new(File::create(manifest_dir.join("command_provenance.tsv"))?);
    writeln!(candidate_manifest, "{}", manifest_header())?;
    writeln!(gene_models, "{}", manifest_header())?;
    writeln!(gene_segments, "sample\tfamily_id\tcandidate\tmodel\tsegment_index\tkind\ttarget_start\ttarget_end\tstrand\tphase\tlength\tdonor\tacceptor\tsplice_class\tobserved")?;
    writeln!(
        id_map,
        "sample\tfamily_id\tinternal_candidate\tcandidate\tsequence_length"
    )?;
    writeln!(warnings, "sample\tfamily_id\tcandidate\tstage\tdetail")?;
    writeln!(fragments, "sample\tfamily_id\tcandidate_a\tmodel_a\tcandidate_b\tmodel_b\tprotein_reference\tcovered_union_fraction\tpadding_nt\tderived_candidate\tstatus\treason")?;
    writeln!(
        provenance,
        "sample\tfamily_id\tcandidate\tminiprot_version\tcommand"
    )?;
    let config = ModelConfig {
        minimum_coverage,
        complete_coverage,
    };
    let mut eligible_counts: BTreeMap<(String, String), usize> = BTreeMap::new();

    for sample_path in sorted_directory_paths(&input.join("samples"))? {
        if !sample_path.is_dir() {
            continue;
        }
        let sample = sample_path
            .file_name()
            .unwrap()
            .to_string_lossy()
            .to_string();
        let candidates_dir = sample_path.join("candidates");
        if !candidates_dir.is_dir() {
            continue;
        }
        for fasta in sorted_directory_paths(&candidates_dir)? {
            let Some(family) = family_id(&fasta) else {
                continue;
            };
            let records = read_dna_fasta(&fasta)?;
            let protein_path = proteins.join(format!("{family}.faa"));
            if !protein_path.is_file() {
                for record in &records {
                    let candidate = record.id.split('|').next_back().unwrap_or(&record.id);
                    write_status_row(
                        &mut candidate_manifest,
                        &sample,
                        &family,
                        candidate,
                        "missing_protein_reference",
                        "missing_protein_reference",
                    )?;
                    write_unresolved_sequence(
                        &unresolved_dir,
                        &family,
                        "missing_protein_reference",
                        &format!("{sample}|{family}|{candidate}"),
                        &record.sequence,
                    )?;
                }
                continue;
            }
            let protein_lengths: BTreeMap<String, usize> = read_raw_fasta(&protein_path)?
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
            let family_work = work_dir.join(&sample).join(&family);
            fs::create_dir_all(&family_work)?;
            let mut candidate_map: BTreeMap<String, (String, String)> = BTreeMap::new();
            for (index, record) in records.iter().enumerate() {
                let internal = format!("TSK{:06}", index + 1);
                let candidate = record
                    .id
                    .split('|')
                    .next_back()
                    .unwrap_or(&record.id)
                    .to_owned();
                writeln!(
                    id_map,
                    "{sample}\t{family}\t{internal}\t{candidate}\t{}",
                    record.sequence.len()
                )?;
                candidate_map.insert(internal, (candidate, record.sequence.clone()));
            }
            let family_raw = raw_dir.join(&sample).join(&family);
            fs::create_dir_all(&family_raw)?;
            let mut built_models = Vec::new();
            let mut handled_candidates = BTreeSet::new();
            for (internal, (candidate, sequence)) in &candidate_map {
                let targets_path = family_work.join(format!("{internal}.targets.fasta"));
                {
                    let mut target_writer = BufWriter::new(File::create(&targets_path)?);
                    writeln!(target_writer, ">{internal}\n{sequence}")?;
                }
                let (result, command_text) = run_miniprot(
                    &miniprot,
                    &targets_path,
                    &protein_path,
                    &threads,
                    max_intron,
                )?;
                writeln!(
                    provenance,
                    "{sample}\t{family}\t{candidate}\t{miniprot_version}\t{command_text}"
                )?;
                fs::write(family_raw.join(format!("{internal}.gff3")), &result.stdout)?;
                fs::write(
                    family_raw.join(format!("{internal}.stderr.txt")),
                    &result.stderr,
                )?;
                if !result.status.success() {
                    handled_candidates.insert(candidate.clone());
                    write_status_row(
                        &mut candidate_manifest,
                        &sample,
                        &family,
                        candidate,
                        "miniprot_failed",
                        "miniprot_failed",
                    )?;
                    write_unresolved_sequence(
                        &unresolved_dir,
                        &family,
                        "miniprot_failed",
                        &format!("{sample}|{family}|{candidate}"),
                        sequence,
                    )?;
                    continue;
                }
                let raw_text = String::from_utf8_lossy(&result.stdout);
                let report = gm2_tools::gene_annotation::parse_miniprot_output(&raw_text);
                for warning in &report.warnings {
                    writeln!(
                        warnings,
                        "{sample}\t{family}\t{candidate}\tparse\t{warning}"
                    )?;
                }
                if report.unassigned_paf > 0 {
                    writeln!(
                        warnings,
                        "{sample}\t{family}\t{candidate}\tparse\tunassigned_paf={}",
                        report.unassigned_paf
                    )?;
                }
                if !report.models.is_empty() {
                    handled_candidates.insert(candidate.clone());
                }
                for mut raw_model in report.models {
                    if raw_model.candidate_id != *internal {
                        writeln!(
                            warnings,
                            "{sample}\t{family}\t{candidate}\tbuild\tunknown_internal_candidate={}",
                            raw_model.candidate_id
                        )?;
                        continue;
                    }
                    raw_model.candidate_id = candidate.clone();
                    match build_gene_model(
                        &raw_model,
                        sequence,
                        protein_lengths.get(&raw_model.protein_id).copied(),
                        config,
                    ) {
                        Ok(model) => built_models.push(model),
                        Err(error) => {
                            writeln!(warnings, "{sample}\t{family}\t{candidate}\tbuild\t{error}")?;
                            write_status_row(
                                &mut candidate_manifest,
                                &sample,
                                &family,
                                candidate,
                                "model_build_failed",
                                &error,
                            )?;
                            write_unresolved_sequence(
                                &unresolved_dir,
                                &family,
                                "model_build_failed",
                                &format!(
                                    "{sample}|{family}|{candidate}|{}",
                                    safe_identifier(&raw_model.gff_id)
                                ),
                                sequence,
                            )?;
                        }
                    }
                }
            }
            for (candidate, sequence) in candidate_map.values() {
                if !handled_candidates.contains(candidate) {
                    write_status_row(
                        &mut candidate_manifest,
                        &sample,
                        &family,
                        candidate,
                        "protein_unsupported",
                        "no_miniprot_model",
                    )?;
                    write_unresolved_sequence(
                        &unresolved_dir,
                        &family,
                        "protein_unsupported",
                        &format!("{sample}|{family}|{candidate}"),
                        sequence,
                    )?;
                }
            }
            built_models.sort_by(|left, right| {
                left.candidate_id
                    .cmp(&right.candidate_id)
                    .then_with(|| {
                        left.genomic_interval
                            .start
                            .cmp(&right.genomic_interval.start)
                    })
                    .then_with(|| left.protein_id.cmp(&right.protein_id))
                    .then_with(|| left.gff_id.cmp(&right.gff_id))
            });
            select_competing_models(&mut built_models);
            let complete_proteins: BTreeSet<String> = built_models
                .iter()
                .filter(|model| {
                    model.selected
                        && model.eligible_for_resolve
                        && model.state == ModelState::Complete
                })
                .map(|model| model.protein_id.clone())
                .collect();
            let mut compatible_pairs = Vec::new();
            for left in 0..built_models.len() {
                for right in left + 1..built_models.len() {
                    let a = &built_models[left];
                    let b = &built_models[right];
                    if !a.selected
                        || !b.selected
                        || !a.eligible_for_resolve
                        || !b.eligible_for_resolve
                        || a.state != ModelState::TerminalPartial
                        || b.state != ModelState::TerminalPartial
                        || a.candidate_id == b.candidate_id
                        || a.protein_id != b.protein_id
                        || a.protein_len_aa != b.protein_len_aa
                        || complete_proteins.contains(&a.protein_id)
                    {
                        continue;
                    }
                    if let Some(fraction) = complete_fragment_union_fraction(
                        a.query_interval_aa,
                        b.query_interval_aa,
                        a.protein_len_aa,
                        complete_coverage,
                    ) {
                        let pair = if a.query_interval_aa.start <= b.query_interval_aa.start {
                            (left, right, fraction)
                        } else {
                            (right, left, fraction)
                        };
                        compatible_pairs.push(pair);
                    }
                }
            }
            let mut fragment_outcomes = Vec::new();
            if compatible_pairs.len() > 1 {
                for (left, right, fraction) in compatible_pairs {
                    fragment_outcomes.push(FragmentOutcome {
                        left_model: left,
                        right_model: right,
                        covered_union_fraction: fraction,
                        padding_nt: 0,
                        derived_candidate: String::new(),
                        status: "ambiguous_fragment_set".to_owned(),
                        reason: "multiple_compatible_pairs".to_owned(),
                    });
                }
            } else if let Some((left, right, fraction)) = compatible_pairs.pop() {
                if fragment_padding == 0 {
                    fragment_outcomes.push(FragmentOutcome {
                        left_model: left,
                        right_model: right,
                        covered_union_fraction: fraction,
                        padding_nt: 0,
                        derived_candidate: String::new(),
                        status: "compatible_unstitched".to_owned(),
                        reason: "fragment_padding_disabled".to_owned(),
                    });
                } else {
                    let left_model = built_models[left].clone();
                    let right_model = built_models[right].clone();
                    let left_sequence = candidate_map
                        .values()
                        .find(|(candidate, _)| candidate == &left_model.candidate_id)
                        .map(|(_, sequence)| oriented_candidate(sequence, left_model.strand));
                    let right_sequence = candidate_map
                        .values()
                        .find(|(candidate, _)| candidate == &right_model.candidate_id)
                        .map(|(_, sequence)| oriented_candidate(sequence, right_model.strand));
                    let derived_candidate = "padded_join_1".to_owned();
                    let derived_internal = "TSKPAD001".to_owned();
                    let mut accepted_model = None;
                    let mut failure_reason = "source_candidate_missing".to_owned();
                    if let (Some(left_sequence), Some(right_sequence)) =
                        (left_sequence, right_sequence)
                    {
                        let padding_start = left_sequence.len();
                        let padding_end = padding_start + fragment_padding;
                        let mut padded_sequence = String::with_capacity(
                            left_sequence.len() + fragment_padding + right_sequence.len(),
                        );
                        padded_sequence.push_str(&left_sequence);
                        padded_sequence.push_str(&"N".repeat(fragment_padding));
                        padded_sequence.push_str(&right_sequence);
                        let target_path =
                            family_work.join(format!("{derived_internal}.targets.fasta"));
                        {
                            let mut writer = BufWriter::new(File::create(&target_path)?);
                            writeln!(writer, ">{derived_internal}\n{padded_sequence}")?;
                        }
                        let (result, command_text) = run_miniprot(
                            &miniprot,
                            &target_path,
                            &protein_path,
                            &threads,
                            max_intron,
                        )?;
                        writeln!(
                            provenance,
                            "{sample}\t{family}\t{derived_candidate}\t{miniprot_version}\t{command_text}"
                        )?;
                        fs::write(
                            family_raw.join(format!("{derived_internal}.gff3")),
                            &result.stdout,
                        )?;
                        fs::write(
                            family_raw.join(format!("{derived_internal}.stderr.txt")),
                            &result.stderr,
                        )?;
                        if result.status.success() {
                            let raw_text = String::from_utf8_lossy(&result.stdout);
                            let report =
                                gm2_tools::gene_annotation::parse_miniprot_output(&raw_text);
                            for warning in &report.warnings {
                                writeln!(
                                    warnings,
                                    "{sample}\t{family}\t{derived_candidate}\tparse\t{warning}"
                                )?;
                            }
                            let mut padded_models = Vec::new();
                            for mut raw_model in report.models.into_iter().filter(|model| {
                                model.candidate_id == derived_internal
                                    && model.protein_id == left_model.protein_id
                            }) {
                                raw_model.candidate_id = derived_candidate.clone();
                                match build_gene_model(
                                    &raw_model,
                                    &padded_sequence,
                                    protein_lengths.get(&raw_model.protein_id).copied(),
                                    config,
                                ) {
                                    Ok(model) => padded_models.push(model),
                                    Err(error) => {
                                        writeln!(warnings, "{sample}\t{family}\t{derived_candidate}\tbuild\t{error}")?;
                                    }
                                }
                            }
                            select_competing_models(&mut padded_models);
                            let mut valid = padded_models.into_iter().filter(|model| {
                                model.selected
                                    && model.eligible_for_resolve
                                    && model.state == ModelState::Complete
                                    && model.strand == Strand::Forward
                                    && model.genomic_interval.start <= padding_start
                                    && model.genomic_interval.end >= padding_end
                                    && model.introns.iter().any(|intron| {
                                        intron.interval.start <= padding_start
                                            && intron.interval.end >= padding_end
                                    })
                                    && model.exons.iter().all(|exon| {
                                        exon.interval.overlap(Interval {
                                            start: padding_start,
                                            end: padding_end,
                                        }) == 0
                                    })
                            });
                            if let Some(mut model) = valid.next() {
                                if valid.next().is_none() {
                                    model.qc_flags.push("padded_fragment_join".to_owned());
                                    accepted_model = Some((model, padded_sequence));
                                    failure_reason.clear();
                                } else {
                                    failure_reason = "multiple_validated_padded_models".to_owned();
                                }
                            } else {
                                failure_reason = "padded_model_not_validated".to_owned();
                            }
                        } else {
                            failure_reason = "padded_miniprot_failed".to_owned();
                        }
                    }
                    if let Some((mut model, padded_sequence)) = accepted_model {
                        built_models[left].selected = false;
                        built_models[left].eligible_for_resolve = false;
                        built_models[left]
                            .qc_flags
                            .push("superseded_by_padded_join".to_owned());
                        built_models[right].selected = false;
                        built_models[right].eligible_for_resolve = false;
                        built_models[right]
                            .qc_flags
                            .push("superseded_by_padded_join".to_owned());
                        model.competition_group = built_models
                            .iter()
                            .map(|model| model.competition_group)
                            .max()
                            .unwrap_or(0)
                            + 1;
                        candidate_map.insert(
                            derived_internal.clone(),
                            (derived_candidate.clone(), padded_sequence.clone()),
                        );
                        writeln!(
                            id_map,
                            "{sample}\t{family}\t{derived_internal}\t{derived_candidate}\t{}",
                            padded_sequence.len()
                        )?;
                        built_models.push(model);
                        fragment_outcomes.push(FragmentOutcome {
                            left_model: left,
                            right_model: right,
                            covered_union_fraction: fraction,
                            padding_nt: fragment_padding,
                            derived_candidate,
                            status: "padded_validated".to_owned(),
                            reason: String::new(),
                        });
                    } else {
                        fragment_outcomes.push(FragmentOutcome {
                            left_model: left,
                            right_model: right,
                            covered_union_fraction: fraction,
                            padding_nt: fragment_padding,
                            derived_candidate,
                            status: "padded_validation_failed".to_owned(),
                            reason: failure_reason,
                        });
                    }
                }
            }
            let mut model_numbers: BTreeMap<String, usize> = BTreeMap::new();
            let model_names: Vec<String> = built_models
                .iter()
                .map(|model| {
                    let number = model_numbers.entry(model.candidate_id.clone()).or_default();
                    *number += 1;
                    format!("model_{number}")
                })
                .collect();
            for outcome in fragment_outcomes {
                let left = &built_models[outcome.left_model];
                let right = &built_models[outcome.right_model];
                writeln!(
                    fragments,
                    "{sample}\t{family}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}",
                    left.candidate_id,
                    model_names[outcome.left_model],
                    right.candidate_id,
                    model_names[outcome.right_model],
                    left.protein_id,
                    outcome.covered_union_fraction,
                    outcome.padding_nt,
                    outcome.derived_candidate,
                    outcome.status,
                    outcome.reason,
                )?;
            }
            for (model, model_name) in built_models.into_iter().zip(model_names) {
                let header = format!("{sample}|{family}|{}|{model_name}", model.candidate_id);
                let fields = model_manifest_fields(
                    &sample,
                    &family,
                    &model.candidate_id,
                    &model_name,
                    &model,
                );
                writeln!(candidate_manifest, "{}", fields.join("\t"))?;
                writeln!(gene_models, "{}", fields.join("\t"))?;
                for exon in &model.exons {
                    let observed = usize::from(
                        !model
                            .qc_flags
                            .iter()
                            .any(|flag| flag == "padded_fragment_join")
                            || !exon.sequence.contains('N'),
                    );
                    writeln!(gene_segments, "{sample}\t{family}\t{}\t{model_name}\t{}\tcds_exon\t{}\t{}\t{}\t{}\t{}\t\t\t\t{observed}",
                        model.candidate_id, exon.index, exon.interval.start + 1, exon.interval.end, model.strand.as_char(), exon.phase.map(|value| value.to_string()).unwrap_or_else(|| ".".into()), exon.sequence.len())?;
                }
                for intron in &model.introns {
                    let observed = usize::from(
                        !model
                            .qc_flags
                            .iter()
                            .any(|flag| flag == "padded_fragment_join")
                            || !intron.sequence.contains('N'),
                    );
                    writeln!(gene_segments, "{sample}\t{family}\t{}\t{model_name}\t{}\tintron\t{}\t{}\t{}\t.\t{}\t{}\t{}\t{}\t{observed}",
                        model.candidate_id, intron.index, intron.interval.start + 1, intron.interval.end, model.strand.as_char(), intron.sequence.len(), intron.donor, intron.acceptor, intron.splice_class.as_str())?;
                }
                let model_gff_dir = gff_dir.join(&sample).join(&family);
                fs::create_dir_all(&model_gff_dir)?;
                write_normalized_gff(
                    &model_gff_dir.join(format!("{}.{}.gff3", model.candidate_id, model_name)),
                    &sample,
                    &family,
                    &model_name,
                    &model,
                )?;
                if model.selected {
                    for exon in &model.exons {
                        append_fasta(
                            &exon_dir.join(format!("{family}.fasta")),
                            &format!("{header}|exon_{}", exon.index),
                            &exon.sequence,
                        )?;
                    }
                    for intron in &model.introns {
                        append_fasta(
                            &intron_dir.join(format!("{family}.fasta")),
                            &format!("{header}|intron_{}", intron.index),
                            &intron.sequence,
                        )?;
                    }
                    if !model
                        .qc_flags
                        .iter()
                        .any(|flag| flag == "padded_fragment_join")
                    {
                        append_fasta(
                            &gene_dir.join(format!("{family}.fasta")),
                            &header,
                            &model.gene,
                        )?;
                    }
                    append_fasta(
                        &super_dir.join(format!("{family}.fasta")),
                        &header,
                        &model.gene,
                    )?;
                    if flank > 0 {
                        if let Some((_, sequence)) = candidate_map
                            .values()
                            .find(|(candidate, _)| candidate == &model.candidate_id)
                        {
                            let start = model.genomic_interval.start.saturating_sub(flank);
                            let end = (model.genomic_interval.end + flank).min(sequence.len());
                            let part = &sequence[start..end];
                            let flanked = if model.strand == Strand::Reverse {
                                reverse_complement_iupac(part)
                            } else {
                                part.to_owned()
                            };
                            append_fasta(
                                &flanked_dir.join(format!("{family}.fasta")),
                                &header,
                                &flanked,
                            )?;
                        }
                    }
                }
                if model.selected && model.eligible_for_resolve {
                    append_fasta(
                        &cds_dir.join(format!("{family}.fasta")),
                        &header,
                        &model.cds,
                    )?;
                    append_fasta(
                        &protein_dir.join(format!("{family}.fasta")),
                        &header,
                        model.protein.as_deref().unwrap_or(""),
                    )?;
                    *eligible_counts
                        .entry((sample.clone(), family.clone()))
                        .or_default() += 1;
                } else {
                    let unresolved_state = if model
                        .qc_flags
                        .iter()
                        .any(|flag| flag == "superseded_by_padded_join")
                    {
                        "superseded_by_padded_join"
                    } else if !model.selected
                        && model
                            .qc_flags
                            .iter()
                            .any(|flag| flag == "competing_model_not_selected")
                    {
                        "competing_model_not_selected"
                    } else {
                        model.state.as_str()
                    };
                    write_unresolved_sequence(
                        &unresolved_dir,
                        &family,
                        unresolved_state,
                        &header,
                        &model.cds,
                    )?;
                }
            }
        }
    }
    let mut multi = BufWriter::new(File::create(
        manifest_dir.join("long_multiple_candidates.tsv"),
    )?);
    writeln!(multi, "sample\tfamily_id\tresolve_eligible_models")?;
    for ((sample, family), count) in eligible_counts {
        if count > 1 {
            writeln!(multi, "{sample}\t{family}\t{count}")?;
        }
    }
    fs::remove_dir_all(&work_dir)?;
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
    let eligible_headers = resolve_eligible_headers(&input)?;
    let mut seen_eligible_headers = BTreeSet::new();
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
    for fasta in sorted_directory_paths(&cds)? {
        let Some(family) = family_id(&fasta) else {
            continue;
        };
        let records = read_fasta(&fasta)?;
        if let Some(eligible) = &eligible_headers {
            for (header, _) in &records {
                if !eligible.contains(header) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("CDS record is not marked eligible_for_resolve: {header}"),
                    ));
                }
                seen_eligible_headers.insert(header.clone());
            }
        }
        let family_work = work.join(&family);
        fs::create_dir_all(&family_work)?;
        let aa_input = family_work.join("proteins.fasta");
        let mut aa_records = Vec::new();
        let mut cds_by_header = BTreeMap::new();
        let mut invalid_translation = 0usize;
        let mut short_protein = 0usize;
        for (header, cds) in &records {
            let Ok(mut protein) = gm2_tools::gene_annotation::translate_cds(cds) else {
                invalid_translation += 1;
                continue;
            };
            let alignment_cds = trim_terminal_stop(cds, &mut protein);
            if protein.contains('X') || protein.contains('*') {
                invalid_translation += 1;
            } else if protein.len() < min_aa_length {
                short_protein += 1;
            } else {
                cds_by_header.insert(header.clone(), alignment_cds.to_owned());
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
    if let Some(eligible) = eligible_headers {
        if let Some(missing) = eligible.difference(&seen_eligible_headers).next() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("eligible gene model has no CDS record: {missing}"),
            ));
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
            option_positive_usize(&options, "--max-intron", 50_000),
            option_fraction(&options, "--minimum-coverage", 0.20),
            option_fraction(&options, "--complete-coverage", 0.80),
            option_usize(&options, "--flank", 0),
            option_usize(&options, "--fragment-padding", 100),
        ),
        _ => usage(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temporary_test_directory(name: &str) -> PathBuf {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("tipseek_{name}_{}_{nonce}", std::process::id()))
    }

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

    #[test]
    fn codon_backtranslation_requires_exact_cds_consumption() {
        assert_eq!(codon_backtranslate("M-", "ATG").as_deref(), Some("ATG---"));
        assert!(codon_backtranslate("M", "ATGAAA").is_none());
        assert!(codon_backtranslate("MM", "ATG").is_none());
    }

    #[test]
    fn terminal_stop_is_removed_from_alignment_inputs_together() {
        let cds = "ATGAAATTTTAA";
        let mut protein = gm2_tools::gene_annotation::translate_cds(cds).unwrap();
        let alignment_cds = trim_terminal_stop(cds, &mut protein);
        assert_eq!(protein, "MKF");
        assert_eq!(alignment_cds, "ATGAAATTT");
        assert_eq!(
            codon_backtranslate(&protein, alignment_cds).as_deref(),
            Some("ATGAAATTT")
        );
    }

    #[test]
    fn fragment_report_requires_complementary_low_overlap_models() {
        assert_eq!(
            complementary_union_fraction(
                Interval::new(0, 45).unwrap(),
                Interval::new(50, 100).unwrap(),
                100,
            ),
            Some(0.95)
        );
        assert_eq!(
            complementary_union_fraction(
                Interval::new(0, 60).unwrap(),
                Interval::new(55, 100).unwrap(),
                100,
            ),
            None
        );
        assert_eq!(
            complementary_union_fraction(
                Interval::new(0, 70).unwrap(),
                Interval::new(75, 90).unwrap(),
                100,
            ),
            None
        );
    }

    #[test]
    fn reads_resolve_eligibility_from_structured_manifest() {
        let root = temporary_test_directory("eligible_manifest");
        fs::create_dir_all(root.join("models")).unwrap();
        fs::write(
            root.join("models/gene_models.tsv"),
            "sample\tfamily_id\tcandidate\tmodel\teligible_for_resolve\nA\tfam\tc1\tmodel_1\t1\nA\tfam\tc2\tmodel_1\t0\n",
        )
        .unwrap();
        let eligible = resolve_eligible_headers(&root).unwrap().unwrap();
        assert_eq!(eligible, BTreeSet::from(["A|fam|c1|model_1".to_owned()]));
        fs::remove_dir_all(root).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn annotation_writes_structured_coordinate_preserving_outputs() {
        use std::os::unix::fs::PermissionsExt;

        let root = temporary_test_directory("annotation_integration");
        let input = root.join("gene");
        let proteins = root.join("proteins");
        let output = root.join("annotation");
        fs::create_dir_all(input.join("samples/sampleA/candidates")).unwrap();
        fs::create_dir_all(&proteins).unwrap();
        fs::write(
            input.join("samples/sampleA/candidates/fam1.fasta"),
            ">sampleA|fam1|candidate_1\nNATGAAAGTAGTTTCCC\n",
        )
        .unwrap();
        fs::write(proteins.join("fam1.faa"), ">p1\nMKFP\n").unwrap();
        let fake_miniprot = root.join("miniprot");
        fs::write(
            &fake_miniprot,
            "#!/bin/sh\nif [ \"$1\" = \"--version\" ]; then printf '%s\\n' '0.18-r281'; exit 0; fi\nprintf '%s\\n' '##PAF\tp1\t4\t0\t4\t+\tTSK000001\t17\t1\t17\t12\t12\t60\tAS:i:100\tfs:i:0\tst:i:0\tcg:Z:2M4N2M' '##gff-version 3' 'TSK000001\tminiprot\tmRNA\t2\t17\t100\t+\t.\tID=MP1;Target=p1 1 4;Identity=1.0;Positive=1.0;Rank=0;Frameshift=0;StopCodon=0' 'TSK000001\tminiprot\tCDS\t2\t7\t100\t+\t0\tParent=MP1' 'TSK000001\tminiprot\tCDS\t12\t17\t100\t+\t0\tParent=MP1'\n",
        )
        .unwrap();
        let mut permissions = fs::metadata(&fake_miniprot).unwrap().permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(&fake_miniprot, permissions).unwrap();

        annotate(
            input,
            proteins,
            output.clone(),
            fake_miniprot.to_string_lossy().to_string(),
            "1".into(),
            50_000,
            0.20,
            0.80,
            0,
            100,
        )
        .unwrap();

        let cds = fs::read_to_string(output.join("cds/fam1.fasta")).unwrap();
        let introns = fs::read_to_string(output.join("introns/fam1.fasta")).unwrap();
        let models = fs::read_to_string(output.join("models/gene_models.tsv")).unwrap();
        let raw =
            fs::read_to_string(output.join("raw_miniprot/sampleA/fam1/TSK000001.gff3")).unwrap();
        assert!(cds.contains("sampleA|fam1|candidate_1|model_1"));
        assert!(cds.contains("ATGAAATTTCCC"));
        assert!(introns.contains("GTAG"));
        assert!(models.contains("\tcomplete\t"));
        assert!(models.lines().nth(1).unwrap().ends_with("\t1\t1"));
        assert!(raw.contains("##PAF"));
        let id_map = fs::read_to_string(output.join("manifest/id_map.tsv")).unwrap();
        assert!(id_map.contains("sampleA\tfam1\tTSK000001\tcandidate_1\t17"));
        assert!(!output.join(".work").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn annotation_validates_a_unique_n_padded_fragment_join() {
        use std::os::unix::fs::PermissionsExt;

        let root = temporary_test_directory("padded_fragment_join");
        let input = root.join("gene");
        let proteins = root.join("proteins");
        let output = root.join("annotation");
        fs::create_dir_all(input.join("samples/sampleA/candidates")).unwrap();
        fs::create_dir_all(&proteins).unwrap();
        fs::write(
            input.join("samples/sampleA/candidates/fam1.fasta"),
            ">left\nATGAAAGT\n>right\nAGTTTCCCTAA\n",
        )
        .unwrap();
        fs::write(proteins.join("fam1.faa"), ">p1\nMKFP\n").unwrap();
        let fake_miniprot = root.join("miniprot");
        fs::write(
            &fake_miniprot,
            r##"#!/bin/sh
if [ "$1" = "--version" ]; then printf '%s\n' '0.18-r281'; exit 0; fi
target=''
for argument in "$@"; do
    case "$argument" in
        *.targets.fasta) target="$argument" ;;
    esac
done
if grep -q '^>TSK000001$' "$target"; then
    printf '%s\n' '##PAF	p1	4	0	2	+	TSK000001	8	0	6	6	6	60	AS:i:100	fs:i:0	st:i:0	cg:Z:2M' '##gff-version 3' 'TSK000001	miniprot	mRNA	1	6	100	+	.	ID=L1;Target=p1 1 2;Identity=1.0;Positive=1.0;Rank=0' 'TSK000001	miniprot	CDS	1	6	100	+	0	Parent=L1'
elif grep -q '^>TSK000002$' "$target"; then
    printf '%s\n' '##PAF	p1	4	2	4	+	TSK000002	11	2	8	6	6	60	AS:i:100	fs:i:0	st:i:0	cg:Z:2M' '##gff-version 3' 'TSK000002	miniprot	mRNA	3	11	100	+	.	ID=R1;Target=p1 3 4;Identity=1.0;Positive=1.0;Rank=0' 'TSK000002	miniprot	CDS	3	11	100	+	0	Parent=R1'
elif grep -q '^>TSKPAD001$' "$target"; then
    printf '%s\n' '##PAF	p1	4	0	4	+	TSKPAD001	119	0	116	12	12	60	AS:i:200	fs:i:0	st:i:0	cg:Z:2M104N2M' '##gff-version 3' 'TSKPAD001	miniprot	mRNA	1	119	200	+	.	ID=P1;Target=p1 1 4;Identity=1.0;Positive=1.0;Rank=0' 'TSKPAD001	miniprot	CDS	1	6	100	+	0	Parent=P1' 'TSKPAD001	miniprot	CDS	111	119	100	+	0	Parent=P1'
else
    exit 1
fi
"##,
        )
        .unwrap();
        let mut permissions = fs::metadata(&fake_miniprot).unwrap().permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(&fake_miniprot, permissions).unwrap();

        annotate(
            input,
            proteins,
            output.clone(),
            fake_miniprot.to_string_lossy().to_string(),
            "1".into(),
            50_000,
            0.20,
            0.80,
            0,
            100,
        )
        .unwrap();

        let cds = read_fasta(&output.join("cds/fam1.fasta")).unwrap();
        assert_eq!(cds.len(), 1);
        assert!(cds[0].0.contains("padded_join_1"));
        assert_eq!(cds[0].1, "ATGAAATTTCCCTAA");
        let supercontigs = fs::read_to_string(output.join("supercontigs/fam1.fasta")).unwrap();
        assert_eq!(supercontigs.matches('>').count(), 1);
        assert_eq!(supercontigs.matches('N').count(), 100);
        assert!(!output.join("genes/fam1.fasta").exists());
        let fragments = fs::read_to_string(output.join("manifest/fragment_groups.tsv")).unwrap();
        assert!(fragments.contains("\t100\tpadded_join_1\tpadded_validated\t"));
        let models = fs::read_to_string(output.join("models/gene_models.tsv")).unwrap();
        assert_eq!(models.lines().count(), 4);
        assert!(models.contains("superseded_by_padded_join"));
        assert!(models
            .lines()
            .find(|line| line.contains("\tpadded_join_1\t"))
            .unwrap()
            .ends_with("\t1\t1"));
        assert!(output
            .join("raw_miniprot/sampleA/fam1/TSK000001.gff3")
            .is_file());
        assert!(output
            .join("raw_miniprot/sampleA/fam1/TSK000002.gff3")
            .is_file());
        assert!(output
            .join("raw_miniprot/sampleA/fam1/TSKPAD001.gff3")
            .is_file());
        fs::remove_dir_all(root).unwrap();
    }
}
