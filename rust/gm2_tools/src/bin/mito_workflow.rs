//! Small mitochondrial companion for TipSeek.
//! References recruit reads; the existing GM2 UCE assembler builds contigs;
//! this binary only resolves contig overlaps and read-supported mate bridges.
use gm2_tools::fastx::{FastxFormat, FastxReader, FastxRecord};
use gm2_tools::mito_merge::{assemble_and_write, LinkConfig, MitoContig};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::env;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process;

fn die(message: impl AsRef<str>) -> ! {
    eprintln!("Error: {}", message.as_ref());
    process::exit(2)
}

fn options(args: &[String]) -> HashMap<String, String> {
    let mut result = HashMap::new();
    let mut index = 0;
    while index < args.len() {
        let flag = &args[index];
        if !flag.starts_with("--") {
            die(format!("unexpected argument {flag}"));
        }
        index += 1;
        let value = args
            .get(index)
            .cloned()
            .unwrap_or_else(|| die(format!("{flag} requires a value")));
        result.insert(flag.clone(), value);
        index += 1;
    }
    result
}

fn required<'a>(options: &'a HashMap<String, String>, key: &str) -> &'a str {
    options
        .get(key)
        .map(String::as_str)
        .unwrap_or_else(|| die(format!("{key} is required")))
}

fn number(options: &HashMap<String, String>, key: &str, default: usize) -> Result<usize, String> {
    options
        .get(key)
        .map(|value| value.parse::<usize>())
        .transpose()
        .map_err(|_| format!("invalid {key}"))
        .map(|value| value.unwrap_or(default))
}

fn decimal(options: &HashMap<String, String>, key: &str, default: f64) -> Result<f64, String> {
    options
        .get(key)
        .map(|value| value.parse::<f64>())
        .transpose()
        .map_err(|_| format!("invalid {key}"))
        .map(|value| value.unwrap_or(default))
}

fn read_fasta(path: &Path) -> Result<Vec<(String, Vec<u8>)>, String> {
    let file = File::open(path).map_err(|error| error.to_string())?;
    let mut records = Vec::new();
    let mut id = None;
    let mut sequence = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|error| error.to_string())?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = id.take() {
                records.push((previous, std::mem::take(&mut sequence)));
            }
            id = Some(
                header
                    .split_whitespace()
                    .next()
                    .unwrap_or("sequence")
                    .to_string(),
            );
        } else if id.is_some() {
            sequence.extend(
                line.bytes()
                    .filter(|base| base.is_ascii_alphabetic())
                    .map(|base| base.to_ascii_uppercase()),
            );
        }
    }
    if let Some(previous) = id {
        records.push((previous, sequence));
    }
    Ok(records)
}

fn write_fasta(path: &Path, records: &[(String, Vec<u8>)]) -> Result<(), String> {
    let mut writer = BufWriter::new(File::create(path).map_err(|error| error.to_string())?);
    for (id, sequence) in records {
        writeln!(writer, ">{id}").map_err(|error| error.to_string())?;
        for line in sequence.chunks(80) {
            writeln!(writer, "{}", String::from_utf8_lossy(line))
                .map_err(|error| error.to_string())?;
        }
    }
    writer.flush().map_err(|error| error.to_string())
}

fn rc(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

fn circular_slice(sequence: &[u8], start: isize, length: usize) -> Vec<u8> {
    let start = start.rem_euclid(sequence.len() as isize) as usize;
    (0..length)
        .map(|offset| sequence[(start + offset) % sequence.len()])
        .collect()
}

fn parse_location(raw: &str) -> Option<(Vec<(usize, usize)>, bool)> {
    let numbers: Vec<usize> = raw
        .split(|character: char| !character.is_ascii_digit())
        .filter(|item| !item.is_empty())
        .filter_map(|item| item.parse().ok())
        .collect();
    if numbers.len() < 2 || !numbers.len().is_multiple_of(2) {
        return None;
    }
    let segments = numbers
        .chunks_exact(2)
        .map(|pair| (pair[0].saturating_sub(1), pair[1]))
        .collect();
    Some((segments, raw.contains("complement")))
}

fn clean_name(value: &str) -> String {
    value
        .chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || "._-".contains(character) {
                character
            } else {
                '_'
            }
        })
        .collect()
}

fn annotated_features(text: &str) -> Vec<(String, String, String)> {
    let mut records = Vec::new();
    let mut current: Option<(String, String, String)> = None;
    let flush = |item: &mut Option<(String, String, String)>,
                 records: &mut Vec<(String, String, String)>| {
        if let Some((kind, location, name)) = item.take() {
            if matches!(kind.as_str(), "gene" | "rRNA" | "tRNA") && !name.is_empty() {
                records.push((kind, location, name));
            }
        }
    };
    for line in text.lines() {
        let key = line.get(5..21).unwrap_or("").trim();
        if !key.is_empty() {
            flush(&mut current, &mut records);
            current = Some((
                key.to_string(),
                line.get(21..).unwrap_or("").trim().to_string(),
                String::new(),
            ));
        } else if let Some((_, location, name)) = current.as_mut() {
            let value = line.trim();
            if let Some((qualifier, raw)) = value.split_once('=') {
                if matches!(qualifier, "/gene" | "/product" | "/locus_tag") && name.is_empty() {
                    *name = raw.trim_matches('"').to_string();
                }
            } else if !value.is_empty() && !location.is_empty() {
                location.push_str(value);
            }
        }
    }
    flush(&mut current, &mut records);
    records
}

fn prepare(options: &HashMap<String, String>) -> Result<(), String> {
    let genbank =
        fs::read_to_string(required(options, "--input")).map_err(|error| error.to_string())?;
    let output = PathBuf::from(required(options, "--out-dir"));
    let flank = number(options, "--flank", 150)?;
    let tile_length = number(options, "--tile-length", 1200)?;
    let tile_step = number(options, "--tile-step", 600)?;
    if tile_length == 0 || tile_step == 0 || tile_step > tile_length {
        return Err("mitochondrial tiles require 0 < step <= length".into());
    }
    let origin = genbank
        .find("ORIGIN")
        .ok_or("GenBank file has no ORIGIN section")?;
    let genome: Vec<u8> = genbank[origin + 6..]
        .lines()
        .take_while(|line| !line.starts_with("//"))
        .flat_map(|line| {
            line.bytes()
                .filter(|base| base.is_ascii_alphabetic())
                .map(|base| base.to_ascii_uppercase())
        })
        .collect();
    if genome.is_empty() {
        return Err("GenBank sequence is empty".into());
    }
    fs::create_dir_all(output.join("metadata")).map_err(|error| error.to_string())?;
    let mut baits = Vec::new();
    let mut genes = BufWriter::new(
        File::create(output.join("metadata/mitochondrial_genes.tsv"))
            .map_err(|error| error.to_string())?,
    );
    writeln!(
        genes,
        "gene\tstart_0_inclusive\tend_0_exclusive\tstrand\tsegments_0_half_open"
    )
    .map_err(|error| error.to_string())?;
    let mut seen = BTreeMap::<(Vec<(usize, usize)>, bool), String>::new();
    for (_, location, label) in annotated_features(&genbank[..origin]) {
        let Some((segments, reverse)) = parse_location(&location) else {
            continue;
        };
        if segments
            .iter()
            .any(|(start, end)| start >= end || *end > genome.len())
        {
            continue;
        }
        if seen.contains_key(&(segments.clone(), reverse)) {
            continue;
        }
        let base = clean_name(&label);
        let duplicate = seen.values().filter(|name| *name == &base).count();
        let name = if duplicate == 0 {
            base
        } else {
            format!("{base}_{}", duplicate + 1)
        };
        seen.insert((segments.clone(), reverse), name.clone());
        let (first_start, _) = segments[0];
        let (_, last_end) = *segments.last().unwrap();
        let segment_text = segments
            .iter()
            .map(|(start, end)| format!("{start}..{end}"))
            .collect::<Vec<_>>()
            .join(",");
        let mut sequence = circular_slice(&genome, first_start as isize - flank as isize, flank);
        for (start, end) in &segments {
            sequence.extend_from_slice(&genome[*start..*end]);
        }
        sequence.extend(circular_slice(&genome, last_end as isize, flank));
        if reverse {
            sequence = rc(&sequence);
        }
        baits.push((format!("gene_{name}"), sequence));
        writeln!(
            genes,
            "{name}\t{first_start}\t{last_end}\t{}\t{segment_text}",
            if reverse { -1 } else { 1 },
        )
        .map_err(|error| error.to_string())?;
    }
    let padding = tile_step.min(genome.len() / 4);
    baits.push((
        "whole_mitochondrion".into(),
        circular_slice(&genome, -(padding as isize), genome.len() + padding * 2),
    ));
    for (index, start) in (0..genome.len()).step_by(tile_step).enumerate() {
        baits.push((
            format!("circular_tile_{:04}", index + 1),
            circular_slice(&genome, start as isize, tile_length),
        ));
    }
    write_fasta(&output.join("mitochondrion.fasta"), &baits)?;
    write_fasta(
        &output.join("metadata/mitochondrial_reference.fasta"),
        &[("mitochondrial_reference".into(), genome)],
    )?;
    genes.flush().map_err(|error| error.to_string())?;
    Ok(())
}

fn canonical_kmers(sequence: &[u8], k: usize) -> HashSet<Vec<u8>> {
    if k == 0 || sequence.len() < k {
        return HashSet::new();
    }
    sequence
        .windows(k)
        .filter(|part| {
            part.iter()
                .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
        })
        .map(|part| {
            let forward = part.iter().map(u8::to_ascii_uppercase).collect::<Vec<_>>();
            let reverse = rc(&forward);
            forward.min(reverse)
        })
        .collect()
}

fn write_feature_evidence(
    metadata: &Path,
    reference: &[u8],
    assembly: &Path,
    output: &Path,
) -> Result<(), String> {
    let records = read_fasta(assembly)?;
    // Only the selected primary component is evidence for a feature. Including
    // alternative contigs would turn mutually exclusive graph paths into a
    // false "recovered" call.
    let assembled = records
        .first()
        .map(|(_, sequence)| sequence.as_slice())
        .unwrap_or_default();
    let mut report = BufWriter::new(
        File::create(output.join("mitochondrial_feature_evidence.tsv"))
            .map_err(|error| error.to_string())?,
    );
    writeln!(report, "feature\treference_bases\tanchor_k\tanchor_kmers\tmatching_anchors\tanchor_fraction\tanchor_evidence\treference_similarity_interpretation\ttranslation_status")
        .map_err(|error| error.to_string())?;
    let raw = fs::read_to_string(metadata).map_err(|error| error.to_string())?;
    for line in raw.lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        let reverse = fields[3] == "-1";
        let mut feature = Vec::new();
        for interval in fields[4].split(',') {
            let Some((start, end)) = interval.split_once("..") else {
                continue;
            };
            let (Ok(start), Ok(end)) = (start.parse::<usize>(), end.parse::<usize>()) else {
                continue;
            };
            if start < end && end <= reference.len() {
                feature.extend_from_slice(&reference[start..end]);
            }
        }
        if reverse {
            feature = rc(&feature);
        }
        let k = feature.len().min(21);
        let anchors = canonical_kmers(&feature, k);
        let assembled_anchors = canonical_kmers(assembled, k);
        let matched = anchors.intersection(&assembled_anchors).count();
        let fraction = if anchors.is_empty() {
            0.0
        } else {
            matched as f64 / anchors.len() as f64
        };
        // Exact 21-mer sharing measures reference similarity, not feature
        // presence. This remains meaningful for close references, while
        // avoiding a false gene-loss claim for distant samples.
        let evidence = if fraction >= 0.80 {
            "high_anchor_similarity"
        } else if matched > 0 {
            "partial_anchor_similarity"
        } else {
            "no_exact_reference_anchor"
        };
        writeln!(
            report,
            "{}\t{}\t{}\t{}\t{}\t{fraction:.6}\t{}\treference_similarity_only\tnot_checked",
            fields[0],
            feature.len(),
            k,
            anchors.len(),
            matched,
            evidence
        )
        .map_err(|error| error.to_string())?;
    }
    report.flush().map_err(|error| error.to_string())
}

/// One annotation row parsed from `metadata/mitochondrial_genes.tsv`.
struct GeneRow {
    name: String,
    segments: Vec<(usize, usize)>,
    reverse: bool,
    start: usize,
}

fn read_gene_rows(metadata: &Path) -> Result<Vec<GeneRow>, String> {
    let raw = fs::read_to_string(metadata).map_err(|error| error.to_string())?;
    let mut rows = Vec::new();
    for line in raw.lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        let reverse = fields[3] == "-1";
        let start = fields[1].parse::<usize>().unwrap_or(usize::MAX);
        let mut segments = Vec::new();
        for interval in fields[4].split(',') {
            if let Some((from, to)) = interval.split_once("..") {
                if let (Ok(from), Ok(to)) = (from.parse::<usize>(), to.parse::<usize>()) {
                    segments.push((from, to));
                }
            }
        }
        if !segments.is_empty() {
            rows.push(GeneRow {
                name: fields[0].to_string(),
                segments,
                reverse,
                start,
            });
        }
    }
    Ok(rows)
}

/// Reconstruct an annotated gene's nucleotide sequence on its own coding strand.
fn gene_anchor_sequence(reference: &[u8], row: &GeneRow) -> Vec<u8> {
    let mut sequence = Vec::new();
    for (from, to) in &row.segments {
        if from < to && *to <= reference.len() {
            sequence.extend_from_slice(&reference[*from..*to]);
        }
    }
    if row.reverse {
        rc(&sequence)
    } else {
        sequence
    }
}

/// Preference order for the standardization anchor. tRNA-Phe is the canonical
/// vertebrate start (as in mtGrasp); the conserved protein-coding and rRNA
/// genes are near-universal fallbacks so a reproducible start can still be
/// chosen for invertebrate references that lack tRNA-Phe.
fn anchor_rank(name: &str) -> usize {
    let lowered = name.to_ascii_lowercase();
    const ORDER: [&str; 8] = ["phe", "trnf", "cox1", "coi", "co1", "nad1", "rrns", "12s"];
    ORDER
        .iter()
        .position(|key| lowered.contains(key))
        .unwrap_or(usize::MAX)
}

struct Standardized {
    sequence: Vec<u8>,
    anchor: String,
    strand: char,
    offset: usize,
    mismatches: usize,
}

/// Rotate a verified circular assembly to a reproducible gene start and place it
/// on that gene's coding strand. Only existing assembled bases are reordered or
/// reverse-complemented; no reference base ever enters the sequence. Returns
/// `None` when no annotated anchor can be located confidently, leaving the
/// audited assembly unchanged.
fn standardize_circular(
    sequence: &[u8],
    rows: &[GeneRow],
    reference: &[u8],
) -> Option<Standardized> {
    if sequence.len() < 100 {
        return None;
    }
    let anchor = rows
        .iter()
        .min_by_key(|row| (anchor_rank(&row.name), row.start))?;
    let gene = gene_anchor_sequence(reference, anchor);
    let lead_len = gene.len().min(120);
    if lead_len < 20 {
        return None;
    }
    let lead = &gene[..lead_len];
    // A divergent sample still standardizes: the leading window may differ by up
    // to 15% before the anchor is judged absent and the assembly left as is.
    let threshold = (lead_len * 15 / 100).max(1);
    let mut best: Option<(usize, usize, bool)> = None;
    for reversed in [false, true] {
        let oriented = if reversed {
            rc(sequence)
        } else {
            sequence.to_vec()
        };
        let mut doubled = oriented.clone();
        doubled.extend_from_slice(&oriented[..lead_len]);
        for offset in 0..oriented.len() {
            let mismatches = doubled[offset..offset + lead_len]
                .iter()
                .zip(lead)
                .filter(|(a, b)| a != b)
                .count();
            if best.is_none_or(|(_, best_mismatches, _)| mismatches < best_mismatches) {
                best = Some((offset, mismatches, reversed));
                if mismatches == 0 {
                    break;
                }
            }
        }
        if best.is_some_and(|(_, mismatches, _)| mismatches == 0) {
            break;
        }
    }
    let (offset, mismatches, reversed) = best?;
    if mismatches > threshold {
        return None;
    }
    let oriented = if reversed {
        rc(sequence)
    } else {
        sequence.to_vec()
    };
    let mut rotated = oriented[offset..].to_vec();
    rotated.extend_from_slice(&oriented[..offset]);
    Some(Standardized {
        sequence: rotated,
        anchor: anchor.name.clone(),
        strand: if reversed { '-' } else { '+' },
        offset,
        mismatches,
    })
}

fn write_standardized(
    output: &Path,
    reference: &[u8],
    metadata: &Path,
    assembly: &Path,
) -> Result<(), String> {
    let rows = read_gene_rows(metadata)?;
    let primary = read_fasta(assembly)?
        .into_iter()
        .next()
        .map(|(_, sequence)| sequence)
        .unwrap_or_default();
    let Some(result) = standardize_circular(&primary, &rows, reference) else {
        return Ok(());
    };
    let header = format!(
        "mito_standardized anchor={} strand={} rotation_offset={} anchor_lead_mismatches={} length={}",
        result.anchor,
        result.strand,
        result.offset,
        result.mismatches,
        result.sequence.len(),
    );
    write_fasta(
        &output.join("mitochondrial_standardized.fasta"),
        &[(header, result.sequence)],
    )
}

fn finalize(options: &HashMap<String, String>) -> Result<(), String> {
    let reference = read_fasta(Path::new(required(options, "--reference-genome")))?
        .into_iter()
        .next()
        .ok_or("empty mitochondrial reference")?
        .1;
    let require_circular = options
        .get("--require-circular")
        .is_some_and(|value| matches!(value.as_str(), "true" | "1" | "yes"));
    let contig_path = Path::new(required(options, "--contigs"));
    let contigs: Vec<MitoContig> = if contig_path.is_file() {
        read_fasta(contig_path)?
            .into_iter()
            .map(|(id, sequence)| MitoContig { id, sequence })
            .collect()
    } else {
        Vec::new()
    };
    if contigs.is_empty() {
        if require_circular {
            return Err(
                "GM2 mitochondrial assembly is no_contigs; a read-supported circular genome is required"
                    .into(),
            );
        }
        // An adaptive depth is a checkpoint, not a complete run. At shallow
        // depths there may not yet be enough recruited reads to produce a
        // contig, so record that observation and let the caller try more reads.
        let output = Path::new(required(options, "--out-dir"));
        fs::create_dir_all(output).map_err(|error| error.to_string())?;
        write_fasta(&output.join("mitochondrial_assembly.fasta"), &[])?;
        fs::write(
            output.join("mitochondrial_assembly_summary.tsv"),
            "status\tno_contigs\nresolution_reason\tno_assembled_contigs\ncomponent_count\t0\n",
        )
        .map_err(|error| error.to_string())?;
        return Ok(());
    }
    let config = LinkConfig {
        minimum_overlap: number(options, "--minimum-overlap", 41)?,
        minimum_identity: decimal(options, "--minimum-identity", 0.98)?,
        terminal_window: number(options, "--terminal-window", 500)?,
        link_kmer: number(options, "--link-kmer", 31)?,
        minimum_link_hits: number(options, "--minimum-link-hits", 2)?,
        minimum_pair_support: number(options, "--minimum-pair-support", 3)?,
        bridge_kmer: number(options, "--bridge-kmer", 31)?,
        bridge_minimum_depth: number(options, "--bridge-minimum-depth", 2)?,
        maximum_bridge: number(options, "--maximum-bridge", 1000)?,
        minimum_junction_support: number(options, "--minimum-junction-support", 3)?,
        expected_length: reference.len(),
    };
    if config.minimum_overlap == 0
        || !(0.0..=1.0).contains(&config.minimum_identity)
        || config.link_kmer == 0
        || config.link_kmer > 63
        || config.bridge_kmer == 0
        || config.bridge_kmer > 63
        || config.minimum_pair_support == 0
        || config.minimum_junction_support == 0
    {
        return Err("invalid mitochondrial overlap, link or bridge parameters".into());
    }
    let status = assemble_and_write(
        Path::new(required(options, "--out-dir")),
        contigs,
        Path::new(required(options, "--paired-reads")),
        options.get("--graph").map(Path::new),
        &config,
    )?;
    if let Some(metadata) = options.get("--gene-metadata") {
        let out_dir = Path::new(required(options, "--out-dir"));
        let assembly = out_dir.join("mitochondrial_assembly.fasta");
        write_feature_evidence(Path::new(metadata), &reference, &assembly, out_dir)?;
        // Standardize only a verified circle: reorder its bases to a reproducible
        // gene start and coding strand so assemblies are directly comparable
        // across samples, without altering the audited raw assembly.
        if status == "circular" {
            write_standardized(out_dir, &reference, Path::new(metadata), &assembly)?;
        }
    }
    if require_circular && status != "circular" {
        return Err(format!(
            "GM2 mitochondrial assembly is {status}; a read-supported circular genome is required"
        ));
    }
    Ok(())
}

fn fastq_pair_id(header: &str) -> String {
    header
        .trim_start_matches('@')
        .split_whitespace()
        .next()
        .unwrap_or("")
        .trim_end_matches("/1")
        .trim_end_matches("/2")
        .to_string()
}

fn write_fastq_record(writer: &mut BufWriter<File>, record: &FastxRecord) -> Result<(), String> {
    for line in [
        &record.header,
        &record.sequence,
        &record.plus,
        &record.quality,
    ] {
        writer.write_all(line).map_err(|error| error.to_string())?;
        writer.write_all(b"\n").map_err(|error| error.to_string())?;
    }
    Ok(())
}

fn collapse_baits(options: &HashMap<String, String>) -> Result<(), String> {
    let input = PathBuf::from(required(options, "--input-dir"));
    let output = PathBuf::from(required(options, "--out-dir"));
    let name = options
        .get("--output-name")
        .map(String::as_str)
        .unwrap_or("mitochondrion");
    fs::create_dir_all(&output).map_err(|error| error.to_string())?;
    let mut one = BufWriter::new(
        File::create(output.join(format!("{name}_1.fq"))).map_err(|error| error.to_string())?,
    );
    let mut two = BufWriter::new(
        File::create(output.join(format!("{name}_2.fq"))).map_err(|error| error.to_string())?,
    );
    let mut seen = HashSet::new();
    let mut inputs: Vec<_> = fs::read_dir(&input)
        .map_err(|error| error.to_string())?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.ends_with("_1.fq"))
        })
        .collect();
    inputs.sort();
    for first_path in inputs {
        let second_path = PathBuf::from(first_path.to_string_lossy().replace("_1.fq", "_2.fq"));
        if !second_path.is_file() {
            return Err(format!("missing mate FASTQ for {}", first_path.display()));
        }
        let mut first = FastxReader::open(&first_path, FastxFormat::Fastq)
            .map_err(|error| error.to_string())?;
        let mut second = FastxReader::open(&second_path, FastxFormat::Fastq)
            .map_err(|error| error.to_string())?;
        loop {
            let Some(left) = first.next_record().map_err(|error| error.to_string())? else {
                if second
                    .next_record()
                    .map_err(|error| error.to_string())?
                    .is_some()
                {
                    return Err("paired FASTQ files have different numbers of records".into());
                }
                break;
            };
            let right = second
                .next_record()
                .map_err(|error| error.to_string())?
                .ok_or("paired FASTQ files have different numbers of records")?;
            let id = fastq_pair_id(&String::from_utf8_lossy(&left.header));
            if id != fastq_pair_id(&String::from_utf8_lossy(&right.header)) {
                return Err("bait mate identifiers differ".into());
            }
            if seen.insert(id) {
                for (writer, record) in [(&mut one, &left), (&mut two, &right)] {
                    write_fastq_record(writer, record)?;
                }
            }
        }
    }
    one.flush().map_err(|error| error.to_string())?;
    two.flush().map_err(|error| error.to_string())?;
    Ok(())
}
fn main() {
    let arguments: Vec<String> = env::args().skip(1).collect();
    if arguments.is_empty() || arguments[0] == "--help" {
        println!("Usage: mito_workflow <prepare-reference|collapse-baits|finalize> [options]");
        return;
    }
    let options = options(&arguments[1..]);
    let result = match arguments[0].as_str() {
        "prepare-reference" => prepare(&options),
        "finalize" => finalize(&options),
        "collapse-baits" => collapse_baits(&options),
        _ => Err("unknown subcommand".into()),
    };
    if let Err(error) = result {
        die(error);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn pseudo_dna(length: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..length)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                match (state >> 33) & 3 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                }
            })
            .collect()
    }

    fn rotate(sequence: &[u8], by: usize) -> Vec<u8> {
        let mut rotated = sequence[by..].to_vec();
        rotated.extend_from_slice(&sequence[..by]);
        rotated
    }

    #[test]
    fn standardize_rotates_a_forward_anchor_to_the_start() {
        let reference = pseudo_dna(200, 11);
        let rows = vec![GeneRow {
            name: "cox1".into(),
            segments: vec![(10, 60)],
            reverse: false,
            start: 10,
        }];
        // The circular assembly is the same molecule rotated by 5 bases.
        let assembly = rotate(&reference, 5);
        let result = standardize_circular(&assembly, &rows, &reference).unwrap();
        assert_eq!(result.strand, '+');
        assert_eq!(result.mismatches, 0);
        // Standardized sequence begins exactly at the anchor gene.
        assert_eq!(&result.sequence[..50], &reference[10..60]);
        assert_eq!(result.sequence.len(), reference.len());
    }

    #[test]
    fn standardize_orients_a_reverse_strand_assembly() {
        let reference = pseudo_dna(200, 29);
        let rows = vec![GeneRow {
            name: "cox1".into(),
            segments: vec![(10, 60)],
            reverse: false,
            start: 10,
        }];
        // Same molecule, rotated then reverse-complemented onto the other strand.
        let assembly = rc(&rotate(&reference, 5));
        let result = standardize_circular(&assembly, &rows, &reference).unwrap();
        assert_eq!(result.strand, '-');
        assert_eq!(result.mismatches, 0);
        assert_eq!(&result.sequence[..50], &reference[10..60]);
    }

    #[test]
    fn standardize_declines_when_no_anchor_matches() {
        let reference = pseudo_dna(200, 7);
        let rows = vec![GeneRow {
            name: "cox1".into(),
            segments: vec![(10, 60)],
            reverse: false,
            start: 10,
        }];
        // An unrelated circular sequence shares no anchor: leave it unchanged.
        let assembly = pseudo_dna(200, 99999);
        assert!(standardize_circular(&assembly, &rows, &reference).is_none());
    }

    #[test]
    fn empty_contigs_are_a_nonfatal_adaptive_checkpoint() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = env::temp_dir().join(format!(
            "gm2_mito_empty_checkpoint_{}_{}",
            process::id(),
            unique
        ));
        fs::create_dir_all(&root).unwrap();
        let reference = root.join("reference.fasta");
        let output = root.join("output");
        fs::write(&reference, ">mitochondrion\nACGTACGTACGT\n").unwrap();
        let mut arguments = HashMap::from([
            ("--reference-genome".into(), reference.display().to_string()),
            (
                "--contigs".into(),
                root.join("missing_contigs.fasta").display().to_string(),
            ),
            (
                "--paired-reads".into(),
                root.join("missing_reads.fq").display().to_string(),
            ),
            ("--out-dir".into(), output.display().to_string()),
            ("--require-circular".into(), "false".into()),
        ]);
        finalize(&arguments).unwrap();
        assert_eq!(
            fs::read_to_string(output.join("mitochondrial_assembly_summary.tsv")).unwrap(),
            "status\tno_contigs\nresolution_reason\tno_assembled_contigs\ncomponent_count\t0\n"
        );
        assert_eq!(
            fs::read_to_string(output.join("mitochondrial_assembly.fasta")).unwrap(),
            ""
        );

        arguments.insert("--require-circular".into(), "true".into());
        assert!(finalize(&arguments)
            .unwrap_err()
            .contains("assembly is no_contigs"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn reference_is_one_mitochondrial_locus_with_many_baits() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = env::temp_dir().join(format!("gm2_mito_test_{}_{}", process::id(), unique));
        fs::create_dir_all(&root).unwrap();
        let genbank = root.join("mito.gb");
        let output = root.join("references");
        fs::write(
            &genbank,
            format!("LOCUS TEST 100 bp DNA circular\nFEATURES             Location/Qualifiers\n     gene            1..30\n                     /gene=\"cox1\"\nORIGIN\n        1 {}\n//\n", "acgt".repeat(25)),
        )
        .unwrap();
        prepare(&HashMap::from([
            ("--input".into(), genbank.display().to_string()),
            ("--out-dir".into(), output.display().to_string()),
            ("--flank".into(), "10".into()),
            ("--tile-length".into(), "40".into()),
            ("--tile-step".into(), "20".into()),
        ]))
        .unwrap();
        let fasta_files: Vec<_> = fs::read_dir(&output)
            .unwrap()
            .filter_map(Result::ok)
            .filter(|entry| {
                entry
                    .path()
                    .extension()
                    .is_some_and(|value| value == "fasta")
            })
            .collect();
        assert_eq!(fasta_files.len(), 1);
        assert_eq!(fasta_files[0].file_name(), "mitochondrion.fasta");
        assert!(
            read_fasta(&output.join("mitochondrion.fasta"))
                .unwrap()
                .len()
                > 2
        );
        assert!(output.join("metadata/mitochondrial_genes.tsv").is_file());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn parse_location_preserves_segments_across_the_origin() {
        assert_eq!(
            parse_location("join(16400..16569,1..300)"),
            Some((vec![(16399, 16569), (0, 300)], false))
        );
        assert_eq!(
            parse_location("complement(join(16400..16569,1..300))"),
            Some((vec![(16399, 16569), (0, 300)], true))
        );
    }

    #[test]
    fn feature_evidence_uses_only_the_primary_assembly_component() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = env::temp_dir().join(format!("gm2_mito_evidence_{}_{}", process::id(), unique));
        fs::create_dir_all(&root).unwrap();
        let metadata = root.join("genes.tsv");
        let assembly = root.join("assembly.fasta");
        fs::write(
            &metadata,
            "gene\tstart_0_inclusive\tend_0_exclusive\tstrand\tsegments_0_half_open\nfeature_a\t0\t28\t1\t0..28\n",
        )
        .unwrap();
        fs::write(
            &assembly,
            ">primary\nACGTTGCAGATCCGATGCTAACGGTTAACGAA\n>alternative\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
        )
        .unwrap();
        let reference = b"ACGTTGCAGATCCGATGCTAACGGTTAACGAA";
        write_feature_evidence(&metadata, reference, &assembly, &root).unwrap();
        let report = fs::read_to_string(root.join("mitochondrial_feature_evidence.tsv")).unwrap();
        assert!(report.contains("feature_a\t28\t21\t"));
        assert!(
            report.contains("\thigh_anchor_similarity\treference_similarity_only\tnot_checked\n")
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn feature_evidence_does_not_call_a_distant_reference_feature_missing() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = env::temp_dir().join(format!("gm2_mito_distant_{}_{}", process::id(), unique));
        fs::create_dir_all(&root).unwrap();
        let metadata = root.join("genes.tsv");
        let assembly = root.join("assembly.fasta");
        fs::write(&metadata, "gene\tstart_0_inclusive\tend_0_exclusive\tstrand\tsegments_0_half_open\nfeature_a\t0\t42\t1\t0..42\n").unwrap();
        fs::write(
            &assembly,
            ">primary\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
        )
        .unwrap();
        write_feature_evidence(
            &metadata,
            b"ACGTTGCAGATCCGATGCTAACGGTTAACGAACCGTTCAGGA",
            &assembly,
            &root,
        )
        .unwrap();
        let report = fs::read_to_string(root.join("mitochondrial_feature_evidence.tsv")).unwrap();
        assert!(report
            .contains("\tno_exact_reference_anchor\treference_similarity_only\tnot_checked\n"));
        assert!(!report.contains("not_detected"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn cross_origin_feature_writes_segmented_bait_and_metadata() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = env::temp_dir().join(format!("gm2_mito_cross_{}_{}", process::id(), unique));
        fs::create_dir_all(&root).unwrap();
        let genbank = root.join("mito.gb");
        let output = root.join("references");
        fs::write(
            &genbank,
            "LOCUS TEST 20 bp DNA circular\nFEATURES             Location/Qualifiers\n     gene            join(17..20,1..3)\n                     /gene=\"cross\"\nORIGIN\n        1 acgtacgtacgtacgtacgt\n//\n",
        )
        .unwrap();
        prepare(&HashMap::from([
            ("--input".into(), genbank.display().to_string()),
            ("--out-dir".into(), output.display().to_string()),
            ("--flank".into(), "2".into()),
            ("--tile-length".into(), "10".into()),
            ("--tile-step".into(), "5".into()),
        ]))
        .unwrap();
        let records = read_fasta(&output.join("mitochondrion.fasta")).unwrap();
        let bait = records.iter().find(|(id, _)| id == "gene_cross").unwrap();
        assert_eq!(bait.1, b"GTACGTACGTA");
        let metadata = fs::read_to_string(output.join("metadata/mitochondrial_genes.tsv")).unwrap();
        assert!(metadata.starts_with(
            "gene\tstart_0_inclusive\tend_0_exclusive\tstrand\tsegments_0_half_open\n"
        ));
        assert!(metadata.contains("cross\t16\t3\t1\t16..20,0..3"));
        fs::remove_dir_all(root).unwrap();
    }
}
