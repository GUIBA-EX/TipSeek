use crate::hash::HashMap;
use crate::model::{LocusTask, RefKmer};
use crate::seq::{encode_kmer, reverse_complement_kmer, valid_runs};
use std::collections::HashSet;
use std::fs::{self, File};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static CACHE_TEMP_ID: AtomicU64 = AtomicU64::new(0);
const CACHE_MAGIC: &[u8; 8] = b"GM2RK001";

// 序列统一大写并去空白，后头比较时别让格式小毛病搅和。
fn clean_sequence(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .filter(|base| base.is_ascii_alphabetic())
        .map(|base| base.to_ascii_uppercase())
        .collect()
}

// 小参考文件一次读透，保留标题和序列给后面的 locus 逻辑用。
pub fn read_fasta(path: &Path) -> io::Result<Vec<(String, Vec<u8>)>> {
    let text = fs::read_to_string(path)?;
    let mut records = Vec::new();
    let mut title: Option<String> = None;
    let mut sequence = Vec::new();

    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = title.take() {
                records.push((previous, clean_sequence(&sequence)));
                sequence.clear();
            }
            title = Some(header.split_whitespace().next().unwrap_or("").to_string());
        } else if title.is_some() {
            sequence.extend_from_slice(line.as_bytes());
        }
    }
    if let Some(previous) = title {
        records.push((previous, clean_sequence(&sequence)));
    }
    Ok(records)
}

// 大 read 文件分块过，内存别一口吞得太猛。
pub fn for_each_sequence_chunk<F>(
    path: &Path,
    fasta: bool,
    chunk_size: usize,
    mut visit: F,
) -> io::Result<()>
where
    F: FnMut(&[Vec<u8>]) -> io::Result<()>,
{
    let mut reader = BufReader::new(File::open(path)?);
    let limit = chunk_size.max(1);
    let mut chunk = Vec::with_capacity(limit);
    let mut line = String::new();
    if fasta {
        let mut sequence = Vec::new();
        let mut seen_header = false;
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            if line.starts_with('>') {
                if seen_header {
                    chunk.push(clean_sequence(&sequence));
                    sequence.clear();
                    if chunk.len() == limit {
                        visit(&chunk)?;
                        chunk.clear();
                    }
                }
                seen_header = true;
            } else if seen_header {
                sequence.extend_from_slice(line.as_bytes());
            }
        }
        if seen_header {
            chunk.push(clean_sequence(&sequence));
        }
    } else {
        let mut record = Vec::with_capacity(4);
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            record.push(line.trim_end_matches(['\r', '\n']).to_string());
            if record.len() == 4 {
                if !record[0].starts_with('@') || !record[2].starts_with('+') {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid FASTQ record in {}", path.display()),
                    ));
                }
                chunk.push(clean_sequence(record[1].as_bytes()));
                record.clear();
                if chunk.len() == limit {
                    visit(&chunk)?;
                    chunk.clear();
                }
            }
        }
        if !record.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("truncated FASTQ record in {}", path.display()),
            ));
        }
    }
    if !chunk.is_empty() {
        visit(&chunk)?;
    }
    Ok(())
}

// Only explicitly matching read identifiers establish fragment linkage.
pub type LinkedReads = Vec<Vec<Vec<u8>>>;

#[cfg(test)]
pub fn read_linked_fragments(path: &Path, fasta: bool) -> io::Result<LinkedReads> {
    read_linked_reads_with_quality(path, fasta).map(|(reads, _)| reads)
}

pub fn read_linked_reads_with_quality(
    path: &Path,
    fasta: bool,
) -> io::Result<(LinkedReads, LinkedReads)> {
    let mut fragments: Vec<Vec<Vec<u8>>> = Vec::new();
    let mut qualities: LinkedReads = Vec::new();
    let mut previous: Option<(String, u8)> = None;
    let mut append = |header: &str, sequence: Vec<u8>, quality: Vec<u8>| {
        let identity = normalized_mate_header(header);
        let paired = match (&previous, identity) {
            (Some((name, 1)), Some((current, 2))) => name == current,
            _ => false,
        };
        if paired {
            if let Some(fragment) = fragments.last_mut() {
                fragment.push(sequence);
                qualities.last_mut().unwrap().push(quality);
            }
            previous = None;
        } else {
            fragments.push(vec![sequence]);
            qualities.push(vec![quality]);
            previous = identity.map(|(name, mate)| (name.to_string(), mate));
        }
    };
    if fasta {
        for (header, sequence) in read_fasta(path)? {
            append(&header, sequence, Vec::new());
        }
    } else {
        let mut lines = BufReader::new(File::open(path)?).lines();
        while let Some(header) = lines.next() {
            let header = header?;
            let mut record = Vec::with_capacity(3);
            for _ in 0..3 {
                record.push(lines.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "truncated FASTQ fragment")
                })??);
            }
            if !header.starts_with('@')
                || !record[1].starts_with('+')
                || record[0].len() != record[2].len()
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid FASTQ fragment",
                ));
            }
            if !record[0].bytes().all(|b| b.is_ascii_alphabetic())
                || !record[2].bytes().all(|b| (33..=126).contains(&b))
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid FASTQ bases or Phred+33 qualities",
                ));
            }
            append(
                &header,
                clean_sequence(record[0].as_bytes()),
                record[2].bytes().map(|b| b - 33).collect(),
            );
        }
    }
    Ok((fragments, qualities))
}

fn normalized_mate_header(header: &str) -> Option<(&str, u8)> {
    let mut fields = header.trim_start_matches(['@', '>']).split_whitespace();
    let token = fields.next()?;
    if let Some(name) = token.strip_suffix("/1") {
        Some((name, 1))
    } else if let Some(name) = token.strip_suffix("/2") {
        Some((name, 2))
    } else {
        let mate_field = fields.next()?;
        // ENA/SRA FASTQ can put /1 or /2 on the original instrument ID in
        // the second field, while the first field holds the shared run/read ID.
        if mate_field.ends_with("/1") {
            return Some((token, 1));
        }
        if mate_field.ends_with("/2") {
            return Some((token, 2));
        }
        match mate_field.split(':').next()? {
            "1" => Some((token, 1)),
            "2" => Some((token, 2)),
            _ => None,
        }
    }
}

// Refilter writes linked mates adjacently.  Only enable PE-specific graph
// evidence when the first fragment explicitly advertises that invariant.
pub fn has_interleaved_pairs(path: &Path, fasta: bool) -> io::Result<bool> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut headers = Vec::with_capacity(2);
    let mut line = String::new();
    if fasta {
        while headers.len() < 2 {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            if line.starts_with('>') {
                headers.push(line.trim_end_matches(['\r', '\n']).to_string());
            }
        }
    } else {
        while headers.len() < 2 {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            if !line.starts_with('@') {
                return Ok(false);
            }
            headers.push(line.trim_end_matches(['\r', '\n']).to_string());
            for _ in 0..3 {
                line.clear();
                if reader.read_line(&mut line)? == 0 {
                    return Ok(false);
                }
            }
        }
    }
    let Some((first, first_mate)) = headers
        .first()
        .and_then(|header| normalized_mate_header(header))
    else {
        return Ok(false);
    };
    let Some((second, second_mate)) = headers
        .get(1)
        .and_then(|header| normalized_mate_header(header))
    else {
        return Ok(false);
    };
    Ok(first == second && first_mate == 1 && second_mate == 2)
}

pub fn minimum_sequence_length(
    path: &Path,
    fasta: bool,
    chunk_size: usize,
) -> io::Result<Option<usize>> {
    let mut minimum = None;
    for_each_sequence_chunk(path, fasta, chunk_size, |chunk| {
        for sequence in chunk.iter().filter(|s| !s.is_empty()) {
            minimum = Some(minimum.map_or(sequence.len(), |old: usize| old.min(sequence.len())));
        }
        Ok(())
    })?;
    Ok(minimum)
}

pub fn discover_references(reference: &Path) -> io::Result<Vec<LocusTask>> {
    // 把参考拆成有顺序的 locus 任务，日志和结果才能对得上号。
    let mut paths = Vec::new();
    if reference.is_dir() {
        for entry in fs::read_dir(reference)? {
            let path = entry?.path();
            if path.is_file() && is_fasta(&path) {
                paths.push(path);
            }
        }
    } else if reference.is_file() && is_fasta(reference) {
        paths.push(reference.to_path_buf());
    }
    paths.sort();

    let total = paths.len();
    let mut tasks = Vec::with_capacity(total);
    let mut keys = HashSet::with_capacity(total);
    for (index, path) in paths.into_iter().enumerate() {
        let key = path
            .file_stem()
            .and_then(|value| value.to_str())
            .unwrap_or("")
            .to_string();
        if !keys.insert(key.clone()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("duplicate locus key '{key}' from reference filename stems"),
            ));
        }
        let reference_count = read_fasta(&path)?.len();
        tasks.push(LocusTask {
            key,
            reference_path: path,
            reference_count,
            ordinal: index + 1,
            total,
        });
    }
    Ok(tasks)
}

// 后缀大小写都认，免得参考文件因为写法不同被漏掉。
pub fn is_fasta(path: &Path) -> bool {
    path.extension()
        .and_then(|value| value.to_str())
        .is_some_and(|extension| {
            matches!(
                extension.to_ascii_lowercase().as_str(),
                "fa" | "fas" | "fasta"
            )
        })
}

// 按 locus 找过滤产物，FASTA 优先；找不着就让该位点正常失败。
pub fn find_filtered(output: &Path, key: &str) -> Option<(PathBuf, bool)> {
    let fasta = output.join("filtered").join(format!("{key}.fasta"));
    if fasta.is_file() {
        return Some((fasta, true));
    }
    let fastq = output.join("filtered").join(format!("{key}.fq"));
    if fastq.is_file() {
        return Some((fastq, false));
    }
    None
}

pub fn build_reference_kmers(records: &[(String, Vec<u8>)], k: usize) -> HashMap<u128, RefKmer> {
    // 参考 k-mer 连位置和方向一块儿记住，组装走路时才知道往哪边拐。
    let mut kmers = HashMap::new();
    for (_, sequence) in records {
        let total_kmers = sequence.len().saturating_sub(k) + 1;
        if total_kmers == 0 {
            continue;
        }
        for (run_start, run) in valid_runs(sequence) {
            if run.len() < k {
                continue;
            }

            // Python visits forward suffixes first.
            for start in (0..=run.len() - k).rev() {
                let encoded = encode_kmer(&run[start..start + k]).expect("valid run");
                let j = run.len() - k - start;
                let global_j = sequence.len() - run_start - run.len() + j;
                let position = (((global_j + 1) as f64 / total_kmers as f64) * 1000.0) as i32;
                insert_reference(&mut kmers, encoded, position, false);
            }

            // Then it visits reverse-complement suffixes from the original prefix.
            for start in 0..=run.len() - k {
                let encoded = encode_kmer(&run[start..start + k]).expect("valid run");
                let reverse = reverse_complement_kmer(encoded, k);
                let global_j = run_start + start;
                let position = (((global_j + 1) as f64 / total_kmers as f64) * 1000.0) as i32;
                insert_reference(&mut kmers, reverse, position, true);
            }
        }
    }
    kmers
}

// 同一个参考 k-mer 出现多次就累积深度，位置保留第一次的坐标。
fn insert_reference(
    kmers: &mut HashMap<u128, RefKmer>,
    kmer: u128,
    position: i32,
    is_reverse: bool,
) {
    if let Some(value) = kmers.get_mut(&kmer) {
        value.depth = value.depth.saturating_add(1);
    } else {
        kmers.insert(
            kmer,
            RefKmer {
                depth: 1,
                position,
                is_reverse,
            },
        );
    }
}

// 缓存名带上参考路径、大小、时间和 k，参考一变就自动换新缓存。
fn cache_path(cache_dir: &Path, reference: &Path, k: usize) -> io::Result<PathBuf> {
    let metadata = fs::metadata(reference)?;
    let modified = metadata
        .modified()
        .ok()
        .and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|duration| duration.as_nanos())
        .unwrap_or(0);
    let mut hasher = DefaultHasher::new();
    reference
        .canonicalize()
        .unwrap_or_else(|_| reference.to_path_buf())
        .hash(&mut hasher);
    metadata.len().hash(&mut hasher);
    modified.hash(&mut hasher);
    k.hash(&mut hasher);
    let name = reference
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("reference");
    Ok(cache_dir.join(format!("{name}.k{k}.{:016x}.gm2rk", hasher.finish())))
}

pub fn load_or_build_reference_kmers(
    reference: &Path,
    records: &[(String, Vec<u8>)],
    k: usize,
    cache_dir: Option<&Path>,
) -> io::Result<HashMap<u128, RefKmer>> {
    // 缓存对上就直接用，没对上才重建，别每回都从头搓参考索引。
    let Some(cache_dir) = cache_dir else {
        return Ok(build_reference_kmers(records, k));
    };
    fs::create_dir_all(cache_dir)?;
    let cache = cache_path(cache_dir, reference, k)?;
    if let Ok(kmers) = read_reference_cache(&cache, k) {
        return Ok(kmers);
    }

    let kmers = build_reference_kmers(records, k);
    let temp_id = CACHE_TEMP_ID.fetch_add(1, Ordering::Relaxed);
    let temporary = cache.with_extension(format!("{}.{}.tmp", std::process::id(), temp_id));
    write_reference_cache(&temporary, k, &kmers)?;
    fs::rename(&temporary, &cache)?;
    Ok(kmers)
}

fn read_exact_array<const N: usize>(reader: &mut impl Read) -> io::Result<[u8; N]> {
    let mut bytes = [0_u8; N];
    reader.read_exact(&mut bytes)?;
    Ok(bytes)
}

// 读缓存时逐项核验魔数和 k，旧文件或坏文件宁可不用。
fn read_reference_cache(path: &Path, expected_k: usize) -> io::Result<HashMap<u128, RefKmer>> {
    let mut reader = File::open(path)?;
    if &read_exact_array::<8>(&mut reader)? != CACHE_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "cache magic mismatch",
        ));
    }
    let k = u32::from_le_bytes(read_exact_array::<4>(&mut reader)?) as usize;
    if k != expected_k {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "cache k mismatch",
        ));
    }
    let count = u64::from_le_bytes(read_exact_array::<8>(&mut reader)?) as usize;
    let mut kmers = HashMap::with_capacity(count);
    for _ in 0..count {
        let key = u128::from_le_bytes(read_exact_array::<16>(&mut reader)?);
        let depth = u32::from_le_bytes(read_exact_array::<4>(&mut reader)?);
        let position = i32::from_le_bytes(read_exact_array::<4>(&mut reader)?);
        let reverse = read_exact_array::<1>(&mut reader)?[0] != 0;
        kmers.insert(
            key,
            RefKmer {
                depth,
                position,
                is_reverse: reverse,
            },
        );
    }
    Ok(kmers)
}

// 缓存写成定长二进制记录，后下次启动能麻溜儿读回来。
fn write_reference_cache(path: &Path, k: usize, kmers: &HashMap<u128, RefKmer>) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writer.write_all(CACHE_MAGIC)?;
    writer.write_all(&(k as u32).to_le_bytes())?;
    writer.write_all(&(kmers.len() as u64).to_le_bytes())?;
    for (key, value) in kmers {
        writer.write_all(&key.to_le_bytes())?;
        writer.write_all(&value.depth.to_le_bytes())?;
        writer.write_all(&value.position.to_le_bytes())?;
        writer.write_all(&[u8::from(value.is_reverse)])?;
    }
    writer.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linked_fastq_keeps_phred_positions_and_rejects_corruption() {
        let path = std::env::temp_dir().join(format!("its_quality_{}.fq", std::process::id()));
        fs::write(&path, "@a/1\nacgt\n+\nI-5!\n@a/2\nTGCA\n+\n5555\n").unwrap();
        let (reads, qualities) = read_linked_reads_with_quality(&path, false).unwrap();
        assert_eq!(reads[0][0], b"ACGT");
        assert_eq!(qualities[0][0], vec![40, 12, 20, 0]);
        assert_eq!(qualities[0][1], vec![20; 4]);
        fs::write(&path, "@a/1\nAC G\n+\nIIII\n").unwrap();
        assert!(read_linked_reads_with_quality(&path, false).is_err());
        fs::write(&path, "@a/1\nACGT\n+\nIII \n").unwrap();
        assert!(read_linked_reads_with_quality(&path, false).is_err());
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn sra_second_field_mate_markers_preserve_fragment_linkage() {
        assert_eq!(
            normalized_mate_header("@SRR123.45 INSTRUMENT:1:2:3/1"),
            Some(("SRR123.45", 1))
        );
        assert_eq!(
            normalized_mate_header("@SRR123.45 INSTRUMENT:1:2:3/2"),
            Some(("SRR123.45", 2))
        );
        assert_eq!(normalized_mate_header("@SRR123.45 length=151"), None);
    }

    #[test]
    fn linked_fragments_validate_every_pair_and_keep_orphans() {
        let path = std::env::temp_dir().join(format!("its2_pairs_{}.fq", std::process::id()));
        fs::write(
            &path,
            concat!(
                "@a/1\nACGT\n+\nIIII\n@a/2\nTGCA\n+\nIIII\n",
                "@orphan/1\nAAAA\n+\nIIII\n@other/2\nCCCC\n+\nIIII\n",
                "@b 1:N:0:ACGT\nGGGG\n+\nIIII\n@b 2:N:0:ACGT\nTTTT\n+\nIIII\n",
                "@single\nACAC\n+\nIIII\n@single2\nCACA\n+\nIIII\n"
            ),
        )
        .unwrap();
        let fragments = read_linked_fragments(&path, false).unwrap();
        assert_eq!(
            fragments.iter().map(Vec::len).collect::<Vec<_>>(),
            vec![2, 1, 1, 2, 1, 1]
        );
        assert_eq!(fragments[3][1], b"TTTT");
        fs::write(&path, "@bad/1\nACGT\n+\nIII\n").unwrap();
        assert!(read_linked_fragments(&path, false).is_err());
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn reference_contains_both_orientations() {
        let records = vec![("x".to_string(), b"AACCGGTT".to_vec())];
        let kmers = build_reference_kmers(&records, 4);
        let forward = encode_kmer(b"AACC").unwrap();
        let reverse = reverse_complement_kmer(forward, 4);
        assert!(kmers.contains_key(&forward));
        assert!(kmers.contains_key(&reverse));
    }
}
