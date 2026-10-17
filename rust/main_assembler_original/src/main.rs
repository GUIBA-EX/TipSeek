//! 原版 GeneMiner2 reference assembler 的单线程兼容实现。
//!
//! 这不是 UCE assembler 的回退，也不使用 backbone、救援或并行搜索。
//! 目标是保留上游 Python `main_assembler.py` 的参考辅助 DBG 行为，同时把
//! 热路径中的 DNA 编码、文件读取和容器操作换成 Rust。

use std::collections::{BTreeMap, HashMap, HashSet};
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Debug, Clone)]
struct Args {
    reference: PathBuf,
    output: PathBuf,
    ka: usize,
    k_min: usize,
    k_max: usize,
    limit: u32,
    iteration: usize,
    cov_min: f64,
    soft_boundary: isize,
    _processes: usize,
    trace_dir: Option<PathBuf>,
    reference_cache_dir: Option<PathBuf>,
    exhaust_seeds: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RefInfo {
    depth: u32,
    position: u16,
    reverse: bool,
}

#[derive(Debug, Clone)]
struct KmerInfo {
    depth: u32,
    position: u16,
    reverse: bool,
    ref_depth: u32,
    insertion_order: usize,
}

#[derive(Debug, Clone)]
struct Walk {
    weights: Vec<u32>,
    bases: Vec<u8>,
}

#[derive(Debug, Clone)]
struct Candidate {
    sequence: String,
    weight: u32,
    support: u32,
}
fn usage(exit_code: i32) -> ! {
    eprintln!("Usage: main_assembler-original-rust -r <ref> -o <out> [options]\n  --assembler-reference-cache-dir <dir>  versioned reference k-mer cache");
    std::process::exit(exit_code);
}

fn next_value(iter: &mut impl Iterator<Item = String>, flag: &str) -> String {
    iter.next().unwrap_or_else(|| {
        eprintln!("Missing value for {flag}");
        usage(2);
    })
}

fn parse_args() -> Args {
    let mut reference = None;
    let mut output = None;
    let mut ka = 39;
    let mut k_min = 21;
    let mut k_max = 39;
    let mut limit = 2;
    let mut iteration = 8192;
    let mut cov_min = 0.0;
    let mut soft_boundary = 0;
    let mut trace_dir = None;
    let mut reference_cache_dir = None;
    let mut processes = 1;
    let mut exhaust_seeds = false;
    let mut it = env::args().skip(1);
    while let Some(flag) = it.next() {
        let value = |it: &mut std::iter::Skip<std::env::Args>, flag: &str| next_value(it, flag);
        match flag.as_str() {
            "-h" | "--help" => usage(0),
            "--exhaust-seeds" => exhaust_seeds = true,
            "-r" => reference = Some(PathBuf::from(value(&mut it, "-r"))),
            "-o" => output = Some(PathBuf::from(value(&mut it, "-o"))),
            "-ka" => ka = value(&mut it, "-ka").parse().unwrap_or_else(|_| usage(2)),
            "-k_min" => {
                k_min = value(&mut it, "-k_min")
                    .parse()
                    .unwrap_or_else(|_| usage(2))
            }
            "-k_max" => {
                k_max = value(&mut it, "-k_max")
                    .parse()
                    .unwrap_or_else(|_| usage(2))
            }
            "-limit_count" => {
                limit = value(&mut it, "-limit_count")
                    .parse()
                    .unwrap_or_else(|_| usage(2))
            }
            "-iteration" => {
                iteration = value(&mut it, "-iteration")
                    .parse()
                    .unwrap_or_else(|_| usage(2))
            }
            "--trace-dir" => trace_dir = Some(PathBuf::from(value(&mut it, "--trace-dir"))),
            "--assembler-reference-cache-dir" => {
                reference_cache_dir = Some(PathBuf::from(value(
                    &mut it,
                    "--assembler-reference-cache-dir",
                )))
            }
            "-cov_min" => {
                cov_min = value(&mut it, "-cov_min")
                    .parse()
                    .unwrap_or_else(|_| usage(2))
            }
            "-sb" | "--soft_boundary" => {
                soft_boundary = value(&mut it, &flag).parse().unwrap_or_else(|_| usage(2))
            }
            "-p" | "--processes" => {
                processes = value(&mut it, &flag).parse().unwrap_or_else(|_| usage(2))
            }
            _ => {
                eprintln!("Unknown option: {flag}");
                usage(2);
            }
        }
    }
    let reference = reference.unwrap_or_else(|| usage(2));
    let output = output.unwrap_or_else(|| usage(2));
    if ka > 63 || k_min > 63 || k_max > 63 || k_min == 0 || k_max < k_min {
        eprintln!("k-mer sizes must be in 1..=63 and k_max must be >= k_min");
        usage(2);
    }
    Args {
        reference,
        output,
        ka,
        k_min,
        k_max,
        limit,
        iteration,
        cov_min,
        reference_cache_dir,
        exhaust_seeds,
        soft_boundary,
        _processes: processes,
        trace_dir,
    }
}
const fn build_base_codes() -> [i8; 256] {
    let mut table = [-1; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table[b'U' as usize] = 3;
    table[b'u' as usize] = 3;
    table
}

const BASE_CODES: [i8; 256] = build_base_codes();

#[inline]
fn base_code(base: u8) -> Option<u128> {
    let code = BASE_CODES[base as usize];
    (code >= 0).then_some(code as u128)
}

#[inline]
fn base_char(code: u128) -> u8 {
    b"ACGT"[code as usize]
}

fn normalize_dna(text: &str) -> Vec<u8> {
    let mut normalized = Vec::with_capacity(text.len());
    for base in text.bytes() {
        if let Some(code) = base_code(base) {
            normalized.push(base_char(code));
        }
    }
    normalized
}

/// 从右往左产出窗口；`min_start=1` 精确保留上游图构建漏掉最左窗口的历史边界。
fn for_each_kmer_rtl(seq: &[u8], k: usize, min_start: usize, mut visit: impl FnMut(u128, usize)) {
    if seq.len() < k || min_start > seq.len() - k {
        return;
    }
    let mut start = seq.len() - k;
    let mut value = encode(&seq[start..start + k]).expect("normalized DNA");
    let high_shift = 2 * (k - 1);
    let mut ordinal = 0usize;
    loop {
        visit(value, ordinal);
        if start == min_start {
            break;
        }
        start -= 1;
        value = (base_code(seq[start]).expect("normalized DNA") << high_shift) | (value >> 2);
        ordinal += 1;
    }
}

/// 直接从正向序列生成 RC 的右到左窗口，省掉 reverse-complement 字符串分配。
fn for_each_rc_kmer_rtl(
    seq: &[u8],
    k: usize,
    min_start: usize,
    mut visit: impl FnMut(u128, usize),
) {
    if seq.len() < k || min_start > seq.len() - k {
        return;
    }
    let last_forward_start = seq.len() - k - min_start;
    let mut value = encode(&seq[..k]).expect("normalized DNA");
    let mask = (1u128 << (2 * k)) - 1;
    for start in 0..=last_forward_start {
        visit(reverse_kmer(value, k), start);
        if start < last_forward_start {
            value = ((value << 2) & mask) | base_code(seq[start + k]).expect("normalized DNA");
        }
    }
}

fn encode(seq: &[u8]) -> Option<u128> {
    let mut value = 0;
    for &base in seq {
        value = (value << 2) | base_code(base)?;
    }
    Some(value)
}

fn reverse_kmer(value: u128, k: usize) -> u128 {
    let mut source = value;
    let mut result = 0;
    for _ in 0..k {
        result = (result << 2) | (3 - (source & 3));
        source >>= 2;
    }
    result
}

fn reverse_complement(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|base| match base_code(base) {
            Some(code) => base_char(3 - code) as char,
            None => 'N',
        })
        .collect()
}

fn decode(mut value: u128, k: usize) -> String {
    let mut out = vec![b'A'; k];
    for pos in (0..k).rev() {
        out[pos] = base_char(value & 3);
        value >>= 2;
    }
    String::from_utf8(out).expect("DNA alphabet is UTF-8")
}

fn read_fasta(path: &Path) -> io::Result<Vec<(String, String)>> {
    let reader = BufReader::with_capacity(1024 * 1024, File::open(path)?);
    let mut records = Vec::new();
    let mut title = String::new();
    let mut sequence = String::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if !title.is_empty() {
                records.push((std::mem::take(&mut title), std::mem::take(&mut sequence)));
            }
            title = rest.to_owned();
        } else {
            sequence.push_str(&line);
        }
    }
    if !title.is_empty() {
        records.push((title, sequence));
    }
    Ok(records)
}

fn known_read_segments(sequence: &str) -> Vec<String> {
    let mut result = Vec::new();
    let mut segment = String::with_capacity(sequence.len());
    for base in sequence.bytes() {
        if let Some(code) = base_code(base) {
            segment.push(base_char(code) as char);
        } else if !segment.is_empty() {
            result.push(std::mem::take(&mut segment));
        }
    }
    if !segment.is_empty() {
        result.push(segment);
    }
    result
}

fn read_sequences(path: &Path, fasta: bool, preserve_boundaries: bool) -> io::Result<Vec<String>> {
    let raw = if fasta {
        read_fasta(path)?
            .into_iter()
            .map(|(_, s)| s)
            .collect::<Vec<_>>()
    } else {
        let reader = BufReader::with_capacity(1024 * 1024, File::open(path)?);
        let mut lines = reader.lines();
        let mut sequences = Vec::new();
        while let Some(header) = lines.next() {
            let _ = header?;
            let sequence = lines
                .next()
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "truncated FASTQ"))??;
            let _plus = lines
                .next()
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "truncated FASTQ"))??;
            let _quality = lines
                .next()
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "truncated FASTQ"))??;
            sequences.push(sequence);
        }
        sequences
    };
    if !preserve_boundaries {
        return Ok(raw
            .into_iter()
            .map(|s| String::from_utf8(normalize_dna(&s)).unwrap())
            .collect());
    }
    // Never manufacture an edge across an unknown read base. In exon mode
    // each known block is a path; no sequence from opposite sides is joined.
    let first_split = raw
        .first()
        .is_some_and(|s| s.bytes().any(|b| base_code(b).is_none()));
    let mut known = Vec::with_capacity(raw.len());
    for sequence in raw {
        // Almost all cleaned reads already contain uppercase A/C/G/T. Keep
        // their allocation instead of copying and growing another string.
        if sequence
            .bytes()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
        {
            if !sequence.is_empty() {
                known.push(sequence);
            }
        } else {
            known.extend(known_read_segments(&sequence));
        }
    }
    if first_split {
        // A tiny first block cannot set the representative read slice length.
        if let Some((largest, _)) = known.iter().enumerate().max_by_key(|(_, s)| s.len()) {
            known.swap(0, largest);
        }
    }
    Ok(known)
}
fn fasta_files(reference: &Path) -> io::Result<Vec<(String, PathBuf)>> {
    let mut paths = Vec::new();
    if reference.is_file() {
        paths.push(reference.to_path_buf());
    } else if reference.is_dir() {
        for entry in fs::read_dir(reference)? {
            let path = entry?.path();
            if path.is_file() {
                paths.push(path);
            }
        }
    }
    paths.sort();
    Ok(paths
        .into_iter()
        .filter_map(|path| {
            let ext = path.extension()?.to_str()?.to_ascii_lowercase();
            if ["fa", "fas", "fasta"].contains(&ext.as_str()) {
                Some((path.file_stem()?.to_str()?.to_owned(), path))
            } else {
                None
            }
        })
        .collect())
}

fn build_ref_dict(path: &Path, k: usize) -> io::Result<(HashMap<u128, RefInfo>, u32)> {
    let records = read_fasta(path)?;
    let ref_count = records.len() as u32;
    let mut dict = HashMap::new();
    for (_, raw) in records {
        let seq = normalize_dna(&raw);
        if seq.len() < k {
            continue;
        }
        let count = seq.len() - k + 1;
        for reverse in [false, true] {
            let oriented: Vec<u8> = if reverse {
                reverse_complement(&String::from_utf8_lossy(&seq)).into_bytes()
            } else {
                seq.clone()
            };
            for j in 0..count {
                // 上游把 base-4 整数从最低位往前挪；j=0 对应序列最右端的 k-mer。
                let start = oriented.len() - k - j;
                let key = encode(&oriented[start..start + k]).expect("normalized DNA");
                let entry = dict.entry(key).or_insert(RefInfo {
                    depth: 0,
                    position: 0,
                    reverse,
                });
                entry.depth += 1;
                if entry.depth == 1 {
                    entry.position = (((j + 1) * 1000) / count) as u16;
                    entry.reverse = reverse;
                }
            }
        }
    }
    Ok((dict, ref_count))
}

const REF_CACHE_MAGIC: &[u8; 8] = b"GM2ORC01";
const REF_CACHE_FORMAT: u32 = 1;
const REF_CACHE_IMPLEMENTATION: u32 = 1;

#[derive(Debug, Clone)]
struct RefCacheHeader {
    format_version: u32,
    implementation_version: u32,
    k: u32,
    reference_path: String,
    reference_size: u64,
    reference_mtime_ns: u128,
    reference_count: u32,
    entry_count: u64,
}

fn read_u32(reader: &mut impl Read) -> io::Result<u32> {
    let mut b = [0; 4];
    reader.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}
fn read_u64(reader: &mut impl Read) -> io::Result<u64> {
    let mut b = [0; 8];
    reader.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}
fn read_u128(reader: &mut impl Read) -> io::Result<u128> {
    let mut b = [0; 16];
    reader.read_exact(&mut b)?;
    Ok(u128::from_le_bytes(b))
}

fn reference_identity(path: &Path) -> io::Result<(String, u64, u128)> {
    let canonical = path.canonicalize().unwrap_or_else(|_| path.to_path_buf());
    let metadata = fs::metadata(path)?;
    let modified_ns = metadata
        .modified()
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|t| t.as_nanos())
        .unwrap_or(0);
    Ok((
        canonical.to_string_lossy().into_owned(),
        metadata.len(),
        modified_ns,
    ))
}

fn reference_cache_key(path: &Path, k: usize) -> io::Result<String> {
    let (canonical, size, modified) = reference_identity(path)?;
    let mut hash = 14695981039346656037u64;
    for byte in canonical
        .as_bytes()
        .iter()
        .chain(size.to_le_bytes().iter())
        .chain(modified.to_le_bytes().iter())
        .chain((k as u64).to_le_bytes().iter())
    {
        hash ^= *byte as u64;
        hash = hash.wrapping_mul(1099511628211);
    }
    Ok(format!(
        "ref-v{}-k{}-{hash:016x}.gm2orc",
        REF_CACHE_FORMAT, k
    ))
}

fn write_ref_cache_header(writer: &mut impl Write, h: &RefCacheHeader) -> io::Result<()> {
    writer.write_all(REF_CACHE_MAGIC)?;
    writer.write_all(&h.format_version.to_le_bytes())?;
    writer.write_all(&h.implementation_version.to_le_bytes())?;
    writer.write_all(&h.k.to_le_bytes())?;
    writer.write_all(&(h.reference_path.len() as u32).to_le_bytes())?;
    writer.write_all(h.reference_path.as_bytes())?;
    writer.write_all(&h.reference_size.to_le_bytes())?;
    writer.write_all(&h.reference_mtime_ns.to_le_bytes())?;
    writer.write_all(&h.reference_count.to_le_bytes())?;
    writer.write_all(&h.entry_count.to_le_bytes())
}

fn read_ref_cache_header(reader: &mut impl Read) -> io::Result<RefCacheHeader> {
    let mut magic = [0; 8];
    reader.read_exact(&mut magic)?;
    if &magic != REF_CACHE_MAGIC {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "cache magic"));
    }
    let format_version = read_u32(reader)?;
    let implementation_version = read_u32(reader)?;
    let k = read_u32(reader)?;
    let path_len = read_u32(reader)? as usize;
    if path_len > 1_048_576 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "cache path length",
        ));
    }
    let mut path = vec![0; path_len];
    reader.read_exact(&mut path)?;
    Ok(RefCacheHeader {
        format_version,
        implementation_version,
        k,
        reference_path: String::from_utf8(path)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "cache path"))?,
        reference_size: read_u64(reader)?,
        reference_mtime_ns: read_u128(reader)?,
        reference_count: read_u32(reader)?,
        entry_count: read_u64(reader)?,
    })
}

fn write_ref_cache_entry(writer: &mut impl Write, key: u128, info: RefInfo) -> io::Result<()> {
    writer.write_all(&key.to_le_bytes())?;
    writer.write_all(&info.depth.to_le_bytes())?;
    writer.write_all(&info.position.to_le_bytes())?;
    writer.write_all(&[u8::from(info.reverse)])
}
fn read_ref_cache_entry(reader: &mut impl Read) -> io::Result<(u128, RefInfo)> {
    let key = read_u128(reader)?;
    let depth = read_u32(reader)?;
    let mut p = [0; 2];
    reader.read_exact(&mut p)?;
    let mut s = [0; 1];
    reader.read_exact(&mut s)?;
    if s[0] > 1 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "cache strand"));
    }
    Ok((
        key,
        RefInfo {
            depth,
            position: u16::from_le_bytes(p),
            reverse: s[0] == 1,
        },
    ))
}

fn load_ref_cache(
    cache_path: &Path,
    reference_path: &Path,
    k: usize,
) -> io::Result<Option<(HashMap<u128, RefInfo>, u32)>> {
    let file = match File::open(cache_path) {
        Ok(file) => file,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(None),
        Err(error) => return Err(error),
    };
    let file_size = file.metadata()?.len();
    let mut reader = BufReader::new(file);
    let header = read_ref_cache_header(&mut reader)?;
    let (path, size, modified_ns) = reference_identity(reference_path)?;
    if header.format_version != REF_CACHE_FORMAT
        || header.implementation_version != REF_CACHE_IMPLEMENTATION
        || header.k != k as u32
        || header.reference_path != path
        || header.reference_size != size
        || header.reference_mtime_ns != modified_ns
        || header.entry_count > usize::MAX as u64
    {
        return Ok(None);
    }
    const ENTRY_BYTES: u64 = 16 + 4 + 2 + 1;
    let header_bytes = reader.stream_position()?;
    let expected_size = header
        .entry_count
        .checked_mul(ENTRY_BYTES)
        .and_then(|bytes| header_bytes.checked_add(bytes))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "cache size overflow"))?;
    if expected_size != file_size {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "cache size mismatch",
        ));
    }
    let mut dictionary = HashMap::with_capacity(header.entry_count as usize);
    for _ in 0..header.entry_count {
        let (key, info) = read_ref_cache_entry(&mut reader)?;
        if dictionary.insert(key, info).is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "duplicate cache k-mer",
            ));
        }
    }
    Ok(Some((dictionary, header.reference_count)))
}

fn write_ref_cache(
    cache_path: &Path,
    reference_path: &Path,
    k: usize,
    dictionary: &HashMap<u128, RefInfo>,
    reference_count: u32,
) -> io::Result<()> {
    let (path, size, modified_ns) = reference_identity(reference_path)?;
    let header = RefCacheHeader {
        format_version: REF_CACHE_FORMAT,
        implementation_version: REF_CACHE_IMPLEMENTATION,
        k: k as u32,
        reference_path: path,
        reference_size: size,
        reference_mtime_ns: modified_ns,
        reference_count,
        entry_count: dictionary.len() as u64,
    };
    let temporary = cache_path.with_extension(format!("gm2orc.{}.tmp", std::process::id()));
    let result = (|| -> io::Result<()> {
        let mut writer = BufWriter::new(File::create(&temporary)?);
        write_ref_cache_header(&mut writer, &header)?;
        let mut entries: Vec<_> = dictionary.iter().map(|(&key, &info)| (key, info)).collect();
        entries.sort_unstable_by_key(|entry| entry.0);
        for (key, info) in entries {
            write_ref_cache_entry(&mut writer, key, info)?;
        }
        writer.flush()?;
        writer.get_ref().sync_all()?;
        drop(writer);
        fs::rename(&temporary, cache_path)?;
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    result
}
fn make_ref_dict(
    path: &Path,
    k: usize,
    cache_dir: Option<&Path>,
) -> io::Result<(HashMap<u128, RefInfo>, u32)> {
    let Some(cache_dir) = cache_dir else {
        return build_ref_dict(path, k);
    };
    if fs::create_dir_all(cache_dir).is_err() {
        return build_ref_dict(path, k);
    }
    let cache_path = cache_dir.join(reference_cache_key(path, k)?);
    match load_ref_cache(&cache_path, path, k) {
        Ok(Some(cached)) => return Ok(cached),
        Ok(None) => {}
        Err(_) => {
            let _ = fs::remove_file(&cache_path);
        }
    }
    let built = build_ref_dict(path, k)?;
    let _ = write_ref_cache(&cache_path, path, k, &built.0, built.1);
    Ok(built)
}
fn make_reads_dict(reads: &[String]) -> (HashMap<String, u32>, usize) {
    let mut slices = HashMap::new();
    let read_len = reads
        .iter()
        .find(|read| !read.is_empty())
        .map_or(0, |read| read.len());
    let slice_len = read_len * 9 / 10;
    if slice_len == 0 {
        return (slices, 0);
    }
    for read in reads {
        for oriented in [read.clone(), reverse_complement(read)] {
            if oriented.len() < slice_len {
                continue;
            }
            let start = (oriented.len() - slice_len) / 2;
            *slices
                .entry(oriented[start..start + slice_len].to_owned())
                .or_insert(0) += 1;
        }
    }
    (slices, slice_len)
}
fn make_graph(
    reads: &[String],
    k: usize,
    reference: &HashMap<u128, RefInfo>,
) -> HashMap<u128, KmerInfo> {
    let estimated = reads
        .iter()
        .map(|read| read.len().saturating_sub(k))
        .sum::<usize>()
        .saturating_mul(2);
    let mut insertion_order = 0usize;
    let mut graph = HashMap::with_capacity(estimated.min(reference.len().saturating_mul(2)));
    for read in reads {
        let bytes = read.as_bytes();
        if bytes.len() <= k {
            continue;
        } // 保留 Python range(0, read_len-k) 的历史边界。
          // 数组排序去重仍按升序入图，与 BTreeSet 的现有顺序一致，但不再逐节点分配。
        let mut unique = Vec::with_capacity((bytes.len() - k).saturating_mul(2));
        for_each_kmer_rtl(bytes, k, 1, |key, _| unique.push(key));
        for_each_rc_kmer_rtl(bytes, k, 1, |key, _| unique.push(key));
        unique.sort_unstable();
        unique.dedup();
        for key in unique {
            let entry = graph.entry(key).or_insert_with(|| {
                let info = match reference.get(&key) {
                    Some(info) => KmerInfo {
                        depth: 0,
                        position: if info.reverse {
                            1000 - info.position
                        } else {
                            info.position
                        },
                        reverse: info.reverse,
                        ref_depth: info.depth,
                        insertion_order,
                    },
                    None => KmerInfo {
                        depth: 0,
                        position: 1023,
                        reverse: true,
                        ref_depth: 0,
                        insertion_order,
                    },
                };
                insertion_order += 1;
                info
            });
            entry.depth += 1;
        }
    }
    graph
}

/// 按上游 `Calculate_Kmer_Size` 的统计规则选自动 k；输入正是 reads_dict 的键。
fn calculate_kmer_size(
    ref_path: &Path,
    reads: &HashMap<String, u32>,
    slice_len: usize,
    mut k_min: usize,
    k_max: usize,
    error_limit: u32,
) -> io::Result<usize> {
    if slice_len <= k_min {
        return Ok(k_min);
    }
    if k_min.is_multiple_of(2) {
        k_min += 1;
    }
    if k_min > k_max {
        return Ok(k_min);
    }
    let mut observed: HashMap<u128, u32> = HashMap::with_capacity(reads.len().saturating_mul(16));
    for read in reads.keys() {
        let bytes = read.as_bytes();
        for_each_kmer_rtl(bytes, k_min, 0, |key, _| {
            *observed.entry(key).or_insert(0) += 1
        });
        for_each_rc_kmer_rtl(bytes, k_min, 0, |key, _| {
            *observed.entry(key).or_insert(0) += 1
        });
    }
    observed.retain(|_, count| *count > error_limit);
    let width = k_max - k_min + 1;
    let mut run_stats = vec![0u32; width];
    for (_, raw) in read_fasta(ref_path)? {
        let sequence = normalize_dna(&raw);
        if sequence.len() < k_min {
            continue;
        }
        let mut runs = vec![0usize];
        for start in (0..=sequence.len() - k_min).rev() {
            let key = encode(&sequence[start..start + k_min]).expect("normalized DNA");
            if observed.contains_key(&key) {
                let last = runs.last_mut().expect("initial run");
                *last += 1;
                if *last >= width {
                    runs.push(width / 2);
                }
            } else if *runs.last().expect("initial run") != 0 {
                runs.push(0);
            }
        }
        for run in runs {
            if run == 0 {
                continue;
            }
            let mut offset = run - 1;
            offset -= offset % 2;
            if offset < width {
                run_stats[offset] += 1;
            }
            for step in (2..=offset).step_by(2) {
                run_stats[offset - step] += 1;
            }
        }
    }
    let Some(upper_offset) = run_stats.iter().rposition(|&n| n > 0) else {
        return Ok(k_min);
    };
    let upper = k_min + upper_offset;
    let lower = upper.div_ceil(2);
    let scores: Vec<(usize, f64)> = run_stats
        .iter()
        .enumerate()
        .filter_map(|(offset, &count)| {
            let k = k_min + offset;
            (lower < k && k <= upper)
                .then_some((k, count as f64 * k as f64 / (slice_len - k + 1) as f64))
        })
        .collect();
    let cutoff = scores.iter().map(|(_, score)| *score).fold(0.0, f64::max) / 2.0;
    Ok(scores
        .iter()
        .rev()
        .find(|(_, score)| *score > cutoff)
        .map(|(k, _)| *k)
        .unwrap_or(k_min))
}

fn median(values: &[u32]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let n = sorted.len();
    Some(if n % 2 == 1 {
        sorted[n / 2] as f64
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) as f64 / 2.0
    })
}

fn quartile(values: &[u32]) -> Option<(f64, f64, f64, u32)> {
    if values.is_empty() {
        return None;
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let n = sorted.len();
    let mid = n / 2;
    let (left, right) = if n % 2 == 1 {
        (&sorted[..mid], &sorted[mid + 1..])
    } else {
        (&sorted[..mid], &sorted[mid..])
    };
    Some((
        median(left).unwrap_or(sorted[0] as f64),
        median(&sorted).unwrap(),
        median(right).unwrap_or(*sorted.last().unwrap() as f64),
        sorted[n - 1] + 1,
    ))
}

fn walk_graph(
    graph: &HashMap<u128, KmerInfo>,
    seed: u128,
    k: usize,
    mut iteration: usize,
) -> (Vec<Walk>, HashSet<u128>, Vec<u16>, u32) {
    let suffix_mask = (1u128 << (2 * k - 2)) - 1;
    let mut path = vec![seed];
    let mut seen: HashMap<u128, u32> = HashMap::from([(seed, 1)]);
    type WeightedNode = (u128, u16, u32);
    type WalkStackEntry = (Vec<WeightedNode>, usize, u16);
    let mut stack: Vec<WalkStackEntry> = Vec::new();
    let mut weights = Vec::new();
    let mut bases = Vec::new();
    let mut walks = Vec::new();
    let mut positions = Vec::new();
    let mut position = 0u16;
    let mut distance = 0usize;
    let mut best = 0u32;
    let mut used: HashSet<u128> = HashSet::from([seed]);
    while iteration > 0 {
        let next_start = (path.last().unwrap() & suffix_mask) << 2;
        let mut nodes: Vec<_> = (0..4u128)
            .filter_map(|base| {
                let key = next_start + base;
                let info = graph.get(&key)?;
                if seen.get(&key).copied().unwrap_or(0) != 0 {
                    None
                } else {
                    Some((key, info.position, info.depth + info.ref_depth))
                }
            })
            .collect();
        nodes.sort_by_key(|node| std::cmp::Reverse(node.2));
        if nodes.is_empty() {
            iteration -= 1;
            let total: u32 = weights.iter().sum();
            walks.push(Walk {
                weights: weights.clone(),
                bases: bases.clone(),
            });
            best = best.max(total);
            for _ in 0..distance {
                if let Some(key) = path.pop() {
                    if let Some(count) = seen.get_mut(&key) {
                        *count -= 1;
                    }
                }
                weights.pop();
                bases.pop();
            }
            let Some((alternatives, saved_distance, saved_position)) = stack.pop() else {
                break;
            };
            nodes = alternatives;
            distance = saved_distance;
            position = saved_position;
        }
        if nodes.len() >= 2 {
            stack.push((nodes[1..].to_vec(), distance, position));
            distance = 0;
        }
        let (key, next_position, weight) = nodes[0];
        if next_position > 0 {
            position = next_position;
        }
        path.push(key);
        *seen.entry(key).or_insert(0) += 1;
        used.insert(key);
        positions.push(next_position);
        weights.push(weight);
        bases.push((key & 3) as u8);
        distance += 1;
    }
    (walks, used, positions, best)
}

fn process_walks(
    mut walks: Vec<Walk>,
    max_weight: u32,
    slice_len: usize,
    reads: &HashMap<String, u32>,
    soft: usize,
) -> Vec<Candidate> {
    for walk in &mut walks {
        if walk.weights.len() > 2 {
            let q1 = quartile(&walk.weights).map(|q| q.0 as u32).unwrap_or(0);
            let cut = walk
                .weights
                .iter()
                .rposition(|&value| value >= q1)
                .unwrap_or(usize::MAX);
            if cut != usize::MAX && cut + soft + 1 < walk.weights.len() {
                walk.weights.truncate(cut + soft + 1);
                walk.bases.truncate(cut + soft + 1);
            }
        }
    }
    let mut out: Vec<Candidate> = walks
        .into_iter()
        .filter_map(|walk| {
            let weight: u32 = walk.weights.iter().sum();
            if weight <= max_weight / 2 {
                return None;
            }
            let sequence = String::from_utf8(
                walk.bases
                    .into_iter()
                    .map(|code| base_char(code as u128))
                    .collect(),
            )
            .unwrap();
            let support = slice_support(&sequence, slice_len, reads);
            Some(Candidate {
                sequence,
                weight,
                support,
            })
        })
        .collect();
    out.sort_by_key(|candidate| std::cmp::Reverse(candidate.support));
    out
}

fn slice_support(sequence: &str, slice_len: usize, reads: &HashMap<String, u32>) -> u32 {
    if slice_len == 0 || sequence.len() <= slice_len {
        return 0;
    }
    let mut total = 0;
    // 仍故意漏掉最后一个窗口，和 Python range(contig_len-slice_len) 对齐。
    for j in 0..(sequence.len() - slice_len) {
        let end = sequence.len() - j;
        total += reads
            .get(&sequence[end - slice_len..end])
            .copied()
            .unwrap_or(0);
    }
    total
}

struct ContigSearch<'a> {
    reads: &'a HashMap<String, u32>,
    slice_len: usize,
    graph: &'a HashMap<u128, KmerInfo>,
    seed: u128,
    k: usize,
    cov_min: f64,
    iteration: usize,
    soft: usize,
}

// Unanchored recovery emits only unambiguous graph paths. It must not spend
// thousands of branch combinations inventing repeat traversals without a reference
// anchor. Branch boundaries remain separate candidates for annotation/read audit.
fn observed_edges(reads: &[String], k: usize) -> HashSet<u128> {
    let mut edges = HashSet::new();
    for read in reads {
        for window in read.as_bytes().windows(k + 1) {
            if let Some(edge) = encode(window) {
                edges.insert(edge);
                edges.insert(reverse_kmer(edge, k + 1));
            }
        }
    }
    edges
}

fn get_unanchored_contig(
    search: ContigSearch<'_>,
    observed: Option<&HashSet<u128>>,
) -> (Vec<Candidate>, HashSet<u128>, i32) {
    let ContigSearch {
        reads,
        slice_len,
        graph,
        seed,
        k,
        ..
    } = search;
    let mask = (1u128 << (2 * k - 2)) - 1;
    let mut used = HashSet::from([seed, reverse_kmer(seed, k)]);
    let mut extend = |start: u128| {
        let mut current = start;
        let mut bases = String::new();
        let mut weight = 0;
        loop {
            let prefix = (current & mask) << 2;
            let next: Vec<_> = (0..4)
                .map(|b| prefix + b)
                .filter(|n| graph.contains_key(n))
                .filter(|n| observed.is_none_or(|e| e.contains(&((current << 2) | (n & 3)))))
                .collect();
            if next.len() != 1 {
                break;
            }
            let node = next[0];
            let incoming = (0..4)
                .map(|b| (b << (2 * k - 2)) | (node >> 2))
                .filter(|n| graph.contains_key(n))
                .filter(|n| observed.is_none_or(|e| e.contains(&((n << 2) | (node & 3)))))
                .count();
            if incoming != 1 || used.contains(&node) || used.contains(&reverse_kmer(node, k)) {
                break;
            }
            used.insert(node);
            used.insert(reverse_kmer(node, k));
            bases.push(base_char(node & 3) as char);
            weight += graph[&node].depth;
            current = node;
        }
        (bases, weight)
    };
    let (right, right_weight) = extend(seed);
    let (left, left_weight) = extend(reverse_kmer(seed, k));
    let sequence = format!("{}{}{}", reverse_complement(&left), decode(seed, k), right);
    let support = slice_support(&sequence, slice_len, reads);
    let candidate = Candidate {
        sequence,
        support,
        weight: left_weight + right_weight,
    };
    let position = graph
        .get(&seed)
        .filter(|info| info.ref_depth > 0)
        .map_or(-2, |info| i32::from(info.position));
    (vec![candidate], used, position)
}

fn get_contigs(search: ContigSearch<'_>) -> (Vec<Candidate>, HashSet<u128>, i32) {
    let ContigSearch {
        reads,
        slice_len,
        graph,
        seed,
        k,
        cov_min,
        iteration,
        soft,
    } = search;
    let (forward, used_f, pos_f, weight_f) = walk_graph(graph, seed, k, iteration);
    let (reverse, used_r, pos_r, weight_r) = walk_graph(graph, reverse_kmer(seed, k), k, iteration);
    let positions: Vec<u32> = pos_f
        .into_iter()
        .chain(pos_r)
        .filter(|&p| p > 0 && p < 1000)
        .map(u32::from)
        .collect();
    let pos = if positions.len() > 1 {
        quartile(&positions).map(|q| q.1 as i32).unwrap_or(-1)
    } else {
        -1
    };
    let mut left = process_walks(reverse, weight_r, slice_len, reads, soft);
    let mut right = process_walks(forward, weight_f, slice_len, reads, soft);
    if left.is_empty() {
        left.push(Candidate {
            sequence: String::new(),
            weight: 0,
            support: 0,
        });
    }
    if right.is_empty() {
        right.push(Candidate {
            sequence: String::new(),
            weight: 0,
            support: 0,
        });
    }
    let seed_seq = decode(seed, k);
    let mut candidates = Vec::new();
    for l in left.iter().take(3) {
        for r in right.iter().take(3) {
            let sequence = format!(
                "{}{}{}",
                reverse_complement(&l.sequence),
                seed_seq,
                r.sequence
            );
            let mut support = 0;
            let mut left_coord = sequence.len();
            let mut right_coord = 0;
            if slice_len > 0 && sequence.len() > slice_len {
                for j in 0..(sequence.len() - slice_len) {
                    let end = sequence.len() - j;
                    if let Some(count) = reads.get(&sequence[end - slice_len..end]) {
                        left_coord = left_coord.min(end - slice_len);
                        right_coord = right_coord.max(end);
                        support += count;
                    }
                }
            }
            let cov_len = right_coord.saturating_sub(left_coord);
            let cov_depth = support as f64 * slice_len as f64 / 0.9;
            if cov_min > 0.0
                && (cov_len == 0
                    || cov_depth / (cov_len as f64) < cov_min
                    || cov_depth / (sequence.len() as f64) < cov_min)
            {
                continue;
            }
            candidates.push(Candidate {
                sequence,
                weight: l.weight + r.weight,
                support,
            });
        }
    }
    let mut used = used_f;
    used.extend(used_r);
    (candidates, used, pos)
}

fn write_trace(
    trace_root: &Path,
    key: &str,
    graph: &HashMap<u128, KmerInfo>,
    seeds: &[(u128, u32, u16, u32, usize)],
) -> io::Result<()> {
    let trace_dir = trace_root.join(key);
    fs::create_dir_all(&trace_dir)?;
    let mut nodes: Vec<_> = graph.iter().collect();
    nodes.sort_by_key(|(_, info)| info.insertion_order);
    let mut graph_out = File::create(trace_dir.join("graph.tsv"))?;
    writeln!(graph_out, "rank\tkmer\tdepth\tposition\treverse\tref_depth")?;
    for (rank, (kmer, info)) in nodes.into_iter().enumerate() {
        writeln!(
            graph_out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            rank,
            kmer,
            info.depth,
            info.position,
            u8::from(info.reverse),
            info.ref_depth
        )?;
    }
    let mut seed_out = File::create(trace_dir.join("seeds.tsv"))?;
    writeln!(
        seed_out,
        "seed_rank\tkmer\tdepth\tposition\tref_depth\tgraph_rank"
    )?;
    for (rank, (kmer, depth, position, ref_depth, insertion)) in seeds.iter().enumerate() {
        writeln!(
            seed_out,
            "{rank}\t{kmer}\t{depth}\t{position}\t{ref_depth}\t{insertion}"
        )?;
    }
    Ok(())
}

fn write_read_trace(trace_root: &Path, key: &str, reads: &[String], k: usize) -> io::Result<()> {
    let trace_dir = trace_root.join(key);
    fs::create_dir_all(&trace_dir)?;
    let mut out = File::create(trace_dir.join("reads.tsv"))?;
    writeln!(out, "read_index\torientation\tsequence\tkmers_sorted")?;
    for (index, read) in reads.iter().take(20).enumerate() {
        for (orientation, sequence) in [
            ("forward", read.clone()),
            ("reverse", reverse_complement(read)),
        ] {
            let bytes = sequence.as_bytes();
            let mut kmers = Vec::new();
            if bytes.len() > k {
                for offset in 0..(bytes.len() - k) {
                    let start = bytes.len() - k - offset;
                    if let Some(kmer) = encode(&bytes[start..start + k]) {
                        kmers.push(kmer);
                    }
                }
            }
            kmers.sort_unstable();
            let text = kmers
                .iter()
                .map(u128::to_string)
                .collect::<Vec<_>>()
                .join(",");
            writeln!(out, "{index}\t{orientation}\t{sequence}\t{text}")?;
        }
    }
    Ok(())
}

fn process_locus(
    args: &Args,
    key: &str,
    ref_path: &Path,
    ref_count: u32,
    _index: usize,
    _total: usize,
) -> io::Result<(String, String, u32)> {
    let sample_dir = &args.output;
    let best_path = sample_dir.join("results").join(format!("{key}.fasta"));
    let all_path = sample_dir.join("contigs_all").join(format!("{key}.fasta"));
    if best_path.exists() {
        return Ok((key.to_owned(), "skipped".to_owned(), 0));
    }
    File::create(&best_path)?;
    let fasta_path = sample_dir.join("filtered").join(format!("{key}.fasta"));
    let fastq_path = sample_dir.join("filtered").join(format!("{key}.fq"));
    let (read_path, fasta) = if fasta_path.exists() {
        (fasta_path, true)
    } else if fastq_path.exists() {
        (fastq_path, false)
    } else {
        fs::remove_file(&best_path)?;
        return Ok((key.to_owned(), "no filtered file".to_owned(), 0));
    };
    let raw_reads = read_sequences(&read_path, fasta, args.exhaust_seeds)?;
    let (reads, slice_len) = make_reads_dict(&raw_reads);
    if reads.is_empty() {
        fs::remove_file(&best_path)?;
        return Ok((key.to_owned(), "no reads".to_owned(), 0));
    }
    let mut k = args.ka;
    if k == 0 {
        k = calculate_kmer_size(
            ref_path, &reads, slice_len, args.k_min, args.k_max, args.limit,
        )?;
    }
    let soft = if args.soft_boundary == -1 {
        slice_len / 2
    } else {
        args.soft_boundary.max(0) as usize
    };
    if let Some(trace_root) = &args.trace_dir {
        write_read_trace(trace_root, key, &raw_reads, k)?;
    }
    let (reference, _) = make_ref_dict(ref_path, k, args.reference_cache_dir.as_deref())?;
    let mut graph = make_graph(&raw_reads, k, &reference);
    // Node overlap alone can invent edges never traversed by any read. Only
    // exon rescue uses observed edges; legacy reference assembly is unchanged.
    let observed = args.exhaust_seeds.then(|| observed_edges(&raw_reads, k));
    if args.limit > 0 {
        graph.retain(|_, info| info.depth > args.limit || info.ref_depth > 0);
    }
    if graph.len() < 3 {
        fs::remove_file(&best_path)?;
        return Ok((key.to_owned(), "insufficient genomic kmers".to_owned(), 0));
    }
    let depths: Vec<u32> = graph.values().map(|info| info.depth).collect();
    let q = quartile(&depths).expect("non-empty graph");
    let upper = ((q.2 - q.0) * 1.5 + q.2) as u32;
    for info in graph.values_mut() {
        if info.ref_depth != 0 {
            info.ref_depth = if info.depth > args.limit {
                (ref_count as f64
                    / ((info.ref_depth as i64 - ref_count as i64).unsigned_abs() + 1) as f64
                    * upper as f64) as u32
                    + 1
            } else {
                1
            };
        }
        info.depth = info.depth.min(upper);
    }
    let mut seeds: Vec<_> = graph
        .iter()
        .filter(|(key, v)| {
            (v.position > 1 && v.position < 1000 && !v.reverse)
                // Recruited reads can encode a conserved protein while sharing
                // no exact nucleotide seed with the reference. Exon mode must
                // also visit their supported graph components. One orientation
                // per unanchored kmer avoids duplicate reverse searches.
                || (args.exhaust_seeds && v.ref_depth == 0 && **key <= reverse_kmer(**key, k))
        })
        .map(|(&kmer, v)| (kmer, v.depth, v.position, v.ref_depth, v.insertion_order))
        .collect();
    seeds.sort_by(|a, b| (b.3, b.1, a.4).cmp(&(a.3, a.1, b.4)));
    if let Some(trace_root) = &args.trace_dir {
        write_trace(trace_root, key, &graph, &seeds)?;
    }
    if seeds.is_empty() {
        fs::remove_file(&best_path)?;
        return Ok((key.to_owned(), "no seed".to_owned(), 0));
    }
    let initial = seeds.len();
    let seed_set: HashSet<u128> = seeds.iter().map(|seed| seed.0).collect();
    let mut high = Vec::new();
    let mut low = Vec::new();
    let seed_ranks: HashMap<u128, usize> = if args.exhaust_seeds {
        seeds.iter().enumerate().map(|(i, s)| (s.0, i)).collect()
    } else {
        HashMap::new()
    };
    let mut visited_seeds = vec![false; if args.exhaust_seeds { seeds.len() } else { 0 }];
    let mut cursor = 0;
    loop {
        let starting_seed = if args.exhaust_seeds {
            while cursor < seeds.len() && visited_seeds[cursor] {
                cursor += 1;
            }
            if cursor == seeds.len() {
                break;
            }
            seeds[cursor].0
        } else {
            if seeds.is_empty() || seeds.len() as f64 <= initial as f64 * 0.5 {
                break;
            }
            seeds[0].0
        };
        let search = ContigSearch {
            reads: &reads,
            slice_len,
            graph: &graph,
            seed: starting_seed,
            k,
            cov_min: args.cov_min,
            iteration: args.iteration,
            soft,
        };
        let (candidates, used, position) = if args.exhaust_seeds {
            get_unanchored_contig(search, observed.as_ref())
        } else {
            get_contigs(search)
        };
        if args.exhaust_seeds {
            visited_seeds[cursor] = true;
            for &key in &used {
                for oriented in [key, reverse_kmer(key, k)] {
                    if let Some(&rank) = seed_ranks.get(&oriented) {
                        visited_seeds[rank] = true;
                    }
                }
            }
        } else {
            seeds.retain(|seed| {
                seed.0 != starting_seed
                    && !used.contains(&seed.0)
                    && !used.contains(&reverse_kmer(seed.0, k))
            });
        }
        let used_seed_count = seed_set.intersection(&used).count();
        for candidate in candidates {
            let position =
                if args.exhaust_seeds && graph[&starting_seed].ref_depth == 0 && position < 0 {
                    -2
                } else {
                    position
                };
            let record = (candidate, used_seed_count, position);
            if record.0.support as usize * slice_len > record.0.sequence.len() {
                high.push(record);
            } else {
                low.push(record);
            }
        }
    }
    let low_quality = high.is_empty();
    let mut selected = if args.exhaust_seeds {
        high.extend(low);
        high
    } else if low_quality {
        low
    } else {
        high
    };
    if selected.is_empty() {
        fs::remove_file(&best_path)?;
        return Ok((key.to_owned(), "no contigs".to_owned(), 0));
    }
    selected.sort_by_key(|entry| std::cmp::Reverse((entry.0.support, entry.0.weight)));
    let mut all = File::create(&all_path)?;
    for (candidate, seed_count, position) in &selected {
        writeln!(
            all,
            ">contig_{}_{}_{}_{}_{}",
            candidate.sequence.len(),
            seed_count,
            position,
            candidate.weight,
            candidate.support
        )?;
        writeln!(all, "{}", candidate.sequence)?;
    }
    let (candidate, seed_count, position) = &selected[0];
    let mut best = File::create(&best_path)?;
    writeln!(
        best,
        ">contig_{}_{}_{}_{}_{}",
        candidate.sequence.len(),
        seed_count,
        position,
        candidate.weight,
        candidate.support
    )?;
    writeln!(best, "{}", candidate.sequence)?;
    Ok((
        key.to_owned(),
        if low_quality {
            "low quality"
        } else {
            "success"
        }
        .to_owned(),
        candidate.support,
    ))
}

fn main() -> io::Result<()> {
    let args = parse_args();
    fs::create_dir_all(args.output.join("results"))?;
    fs::create_dir_all(args.output.join("contigs_all"))?;
    let started = Instant::now();
    let refs = fasta_files(&args.reference)?;
    let mut results = BTreeMap::new();
    for (index, (key, path)) in refs.iter().enumerate() {
        if let Ok((_, count)) = make_ref_dict(
            path,
            args.ka.max(args.k_min),
            args.reference_cache_dir.as_deref(),
        ) {
            if let Ok((name, status, value)) =
                process_locus(&args, key, path, count, index + 1, refs.len())
            {
                if status != "skipped" {
                    results.insert(name, (status, value));
                }
            }
        }
    }
    let mut output = File::create(args.output.join("result_dict.txt"))?;
    for (key, (status, value)) in results {
        writeln!(output, "{key},{status},{value},")?;
    }
    println!("Time cost: {:?}", started.elapsed());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ambiguous_read_bases_are_path_boundaries_in_exon_mode() {
        let sequence = "AACGTNNNTTGCA";
        let segments = known_read_segments(sequence);
        assert_eq!(segments, vec!["AACGT", "TTGCA"]);
        assert!(!observed_edges(&segments, 5).contains(&encode(b"AACGTT").unwrap()));
        assert_eq!(known_read_segments("NNNN"), Vec::<String>::new());
        assert_eq!(known_read_segments("acguNta"), vec!["ACGT", "TA"]);
        assert_eq!(
            String::from_utf8(normalize_dna(sequence)).unwrap(),
            "AACGTTTGCA"
        );
    }

    #[test]
    fn observed_edges_remove_only_unobserved_branches() {
        let k = 5;
        let raw = vec!["AACGATC".to_owned(), "TCGATTG".to_owned()];
        // AACGA -> ACGAT -> CGATC is read-supported. CGATT is
        // also a graph node, but ACGATT was never observed in a read.
        let edges = observed_edges(&raw, k);
        assert!(edges.contains(&encode(b"ACGATC").unwrap()));
        assert!(!edges.contains(&encode(b"ACGATT").unwrap()));
        assert!(edges.contains(&reverse_kmer(encode(b"ACGATC").unwrap(), k + 1)));
        let graph: HashMap<_, _> = ["AACGA", "ACGAT", "CGATC", "CGATT"]
            .into_iter()
            .map(|s| {
                (
                    encode(s.as_bytes()).unwrap(),
                    KmerInfo {
                        depth: 3,
                        position: 0,
                        reverse: false,
                        ref_depth: 0,
                        insertion_order: 0,
                    },
                )
            })
            .collect();
        let reads = HashMap::new();
        let search = || ContigSearch {
            reads: &reads,
            slice_len: 10,
            graph: &graph,
            seed: encode(b"AACGA").unwrap(),
            k,
            cov_min: 0.0,
            iteration: 100,
            soft: 0,
        };
        assert_eq!(
            get_unanchored_contig(search(), None).0[0].sequence,
            "AACGAT"
        );
        assert_eq!(
            get_unanchored_contig(search(), Some(&edges)).0[0].sequence,
            "AACGATC"
        );
        let mut actual_branch = raw;
        actual_branch.push("AACGATT".into());
        let edges = observed_edges(&actual_branch, k);
        assert_eq!(
            get_unanchored_contig(search(), Some(&edges)).0[0].sequence,
            "AACGAT"
        );
    }

    #[test]
    fn unanchored_rescue_does_not_choose_a_branch_or_loop_forever() {
        let raw = vec!["AAACCC".to_owned(), "AAAGGG".to_owned()];
        let graph = make_graph(&raw, 3, &HashMap::new());
        let reads = HashMap::new();
        let (paths, _, position) = get_unanchored_contig(
            ContigSearch {
                reads: &reads,
                slice_len: 5,
                graph: &graph,
                seed: encode(b"AAA").unwrap(),
                k: 3,
                cov_min: 0.0,
                iteration: 100,
                soft: 0,
            },
            None,
        );
        assert_eq!(position, -2);
        assert_eq!(paths[0].sequence, "AAA");
        let graph = make_graph(&["ATGATGATG".to_owned()], 3, &HashMap::new());
        let (paths, used, _) = get_unanchored_contig(
            ContigSearch {
                reads: &reads,
                slice_len: 5,
                graph: &graph,
                seed: encode(b"ATG").unwrap(),
                k: 3,
                cov_min: 0.0,
                iteration: 100,
                soft: 0,
            },
            None,
        );
        assert!(used.len() <= graph.len());
        assert!(paths[0].sequence.len() <= graph.len() + 3);
    }
    #[test]
    fn reference_cache_round_trip_and_rebuilds_invalid_files() -> io::Result<()> {
        let root =
            std::env::temp_dir().join(format!("gm2-original-cache-test-{}", std::process::id()));
        let reference = root.join("reference.fasta");
        let cache = root.join("cache");
        fs::create_dir_all(&cache)?;
        fs::write(
            &reference,
            b">one\nACGTACGTACGTACGT\n>two\nTGCATGCATGCATGCA\n",
        )?;

        let expected = build_ref_dict(&reference, 5)?;
        let cold = make_ref_dict(&reference, 5, Some(&cache))?;
        assert_eq!(cold, expected);
        let mut cache_files: Vec<_> = fs::read_dir(&cache)?
            .map(|entry| entry.unwrap().path())
            .collect();
        assert_eq!(cache_files.len(), 1);
        let cache_path = cache_files.pop().unwrap();
        let original_bytes = fs::read(&cache_path)?;

        let warm = make_ref_dict(&reference, 5, Some(&cache))?;
        assert_eq!(warm, expected);
        assert_eq!(fs::read(&cache_path)?, original_bytes);

        let path_len = u32::from_le_bytes(original_bytes[20..24].try_into().unwrap()) as usize;
        let size_offset = 24 + path_len;
        let mtime_offset = size_offset + 8;
        let mut corruptions = Vec::new();
        let mut bad_magic = original_bytes.clone();
        bad_magic[0] ^= 0xff;
        corruptions.push(bad_magic);
        for offset in [8usize, 12, 16] {
            let mut bytes = original_bytes.clone();
            bytes[offset..offset + 4].copy_from_slice(&99u32.to_le_bytes());
            corruptions.push(bytes);
        }
        let mut bad_path = original_bytes.clone();
        bad_path[24] ^= 1;
        corruptions.push(bad_path);
        let mut bad_size = original_bytes.clone();
        bad_size[size_offset..size_offset + 8].copy_from_slice(&99u64.to_le_bytes());
        corruptions.push(bad_size);
        let mut bad_mtime = original_bytes.clone();
        bad_mtime[mtime_offset..mtime_offset + 16].copy_from_slice(&99u128.to_le_bytes());
        corruptions.push(bad_mtime);
        for corrupted in corruptions {
            fs::write(&cache_path, corrupted)?;
            assert_eq!(make_ref_dict(&reference, 5, Some(&cache))?, expected);
            assert_eq!(fs::read(&cache_path)?, original_bytes);
        }

        fs::write(&cache_path, &original_bytes[..12])?;
        assert_eq!(make_ref_dict(&reference, 5, Some(&cache))?, expected);
        assert_eq!(fs::read(&cache_path)?, original_bytes);
        assert!(!fs::read_dir(&cache)?.any(|entry| entry
            .unwrap()
            .path()
            .to_string_lossy()
            .ends_with(".tmp")));

        fs::remove_dir_all(root)?;
        Ok(())
    }
}
