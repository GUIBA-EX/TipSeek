use crate::mem_index::{
    build_locus_indexes, read_codes_into, LocusMemIndex, MemQueryProfile, MemSeed,
};
use crate::model::{Locus, LocusId};
use ahash::AHashMap;
use gm2_tools::fastx::{FastxFormat, FastxReader};
use std::fs;
use std::path::{Path, PathBuf};

const INVALID: u8 = u8::MAX;

const fn base_table() -> [u8; 256] {
    let mut table = [INVALID; 256];
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

pub const BASE_CODE: [u8; 256] = base_table();

#[inline(always)]
pub fn code(base: u8) -> Option<u8> {
    let value = BASE_CODE[base as usize];
    (value != INVALID).then_some(value)
}

pub fn valid_dna(sequence: &[u8]) -> bool {
    sequence.iter().all(|&base| code(base).is_some())
}

#[derive(Clone, Debug, Default)]
pub struct IndexProfile {
    pub recruit_probes: u64,
    pub recruit_bloom_rejected: u64,
    pub recruit_hits: u64,
    pub exact_locus_queries: u64,
    pub exact_index_queries: u64,
    pub exact_run_windows: u64,
    pub exact_matching_windows: u64,
    pub mem_starts: u64,
    pub mem_bases: u64,
}

/// Four probes confined to one 64-bit word. It has no false negatives; a
/// positive is always confirmed by the exact reference hash table.
#[derive(Debug)]
struct BlockedBloom {
    blocks: Vec<u64>,
    mask: usize,
}

impl BlockedBloom {
    fn for_keys(keys: usize) -> Self {
        let blocks = keys
            .saturating_mul(12)
            .div_ceil(64)
            .max(1)
            .next_power_of_two();
        Self {
            blocks: vec![0; blocks],
            mask: blocks - 1,
        }
    }

    #[inline(always)]
    fn mix(mut value: u64) -> u64 {
        value ^= value >> 30;
        value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value ^= value >> 27;
        value = value.wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^ (value >> 31)
    }

    #[inline(always)]
    fn bit_mask(hash: u64) -> u64 {
        (1_u64 << ((hash >> 16) & 63))
            | (1_u64 << ((hash >> 28) & 63))
            | (1_u64 << ((hash >> 40) & 63))
            | (1_u64 << ((hash >> 52) & 63))
    }

    #[inline(always)]
    fn insert_hash(&mut self, hash: u64) {
        let block = hash as usize & self.mask;
        self.blocks[block] |= Self::bit_mask(hash);
    }

    #[inline(always)]
    fn may_contain_hash(&self, hash: u64) -> bool {
        let block = hash as usize & self.mask;
        self.blocks[block] & Self::bit_mask(hash) == Self::bit_mask(hash)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ExactSeed {
    pub sequence: u32,
    pub read_start: u16,
    pub read_end: u16,
    pub reference_start: u32,
    pub reference_end: u32,
}

impl ExactSeed {
    pub fn len(self) -> usize {
        self.read_end as usize - self.read_start as usize
    }

    pub fn is_empty(self) -> bool {
        self.read_start == self.read_end
    }
}

fn exact_seed(seed: MemSeed) -> ExactSeed {
    ExactSeed {
        sequence: seed.sequence,
        read_start: seed.read_start.min(u16::MAX as usize) as u16,
        read_end: seed.read_end.min(u16::MAX as usize) as u16,
        reference_start: seed.reference_start.min(u32::MAX as usize) as u32,
        reference_end: seed.reference_end.min(u32::MAX as usize) as u32,
    }
}

#[derive(Debug, Default)]
pub struct ReadEvidenceScratch {
    orientation_events: Vec<u8>,
    windows: usize,
    best_exact: Vec<Option<ExactSeed>>,
    occurrence_counts: Vec<u32>,
    strand_masks: Vec<u8>,
    read_codes: Vec<u8>,
}

impl ReadEvidenceScratch {
    fn reset(&mut self, candidate_count: usize, windows: usize) {
        self.windows = windows;
        self.orientation_events.resize(candidate_count * windows, 0);
        self.orientation_events.fill(0);
        self.best_exact.resize(candidate_count, None);
        self.best_exact.fill(None);
        self.occurrence_counts.resize(windows, 0);
        self.strand_masks.resize(windows, 0);
    }

    pub fn orientation(&self, candidate: usize) -> &[u8] {
        let start = candidate * self.windows;
        &self.orientation_events[start..start + self.windows]
    }

    pub fn best(&self, candidate: usize) -> Option<ExactSeed> {
        self.best_exact[candidate]
    }

    pub fn allocated_bytes(&self) -> usize {
        self.orientation_events.capacity() * std::mem::size_of::<u8>()
            + self.best_exact.capacity() * std::mem::size_of::<Option<ExactSeed>>()
            + self.occurrence_counts.capacity() * std::mem::size_of::<u32>()
            + self.strand_masks.capacity() * std::mem::size_of::<u8>()
            + self.read_codes.capacity() * std::mem::size_of::<u8>()
    }
}

/// Reusable, generation-stamped candidate set for one sample. It preserves
/// the previous recruit order while eliminating per-fragment allocation and
/// linear duplicate checks.
#[derive(Debug, Default)]
pub struct RecruitScratch {
    loci: Vec<LocusId>,
    seen_generation: Vec<u32>,
    generation: u32,
}

impl RecruitScratch {
    pub fn begin(&mut self, locus_count: usize) {
        self.loci.clear();
        self.seen_generation.resize(locus_count, 0);
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.seen_generation.fill(0);
            self.generation = 1;
        }
    }

    #[inline(always)]
    fn insert(&mut self, locus: LocusId) {
        let locus_index = locus as usize;
        if self.seen_generation[locus_index] != self.generation {
            self.seen_generation[locus_index] = self.generation;
            self.loci.push(locus);
        }
    }

    pub fn loci(&self) -> &[LocusId] {
        &self.loci
    }

    pub fn sort(&mut self) {
        self.loci.sort_unstable();
    }
}

#[derive(Debug)]
enum LocusHits {
    One(LocusId),
    Many(Vec<LocusId>),
}

impl LocusHits {
    fn insert(&mut self, locus: LocusId) {
        match self {
            Self::One(existing) if *existing != locus => {
                *self = Self::Many(vec![*existing, locus]);
            }
            Self::Many(values) if !values.contains(&locus) => values.push(locus),
            _ => {}
        }
    }

    fn values(&self) -> &[LocusId] {
        match self {
            Self::One(value) => std::slice::from_ref(value),
            Self::Many(values) => values,
        }
    }
}

#[derive(Debug)]
enum RecruitIndex {
    Short(AHashMap<u64, LocusHits>),
    Long(AHashMap<u128, LocusHits>),
}

impl RecruitIndex {
    fn new(k: usize) -> Self {
        if k <= 32 {
            Self::Short(AHashMap::new())
        } else {
            Self::Long(AHashMap::new())
        }
    }

    fn insert_sequence(&mut self, sequence: &[u8], k: usize, locus: LocusId) {
        match self {
            Self::Short(index) => scan_kmers_u64(sequence, k, 1, true, |key, _| {
                index
                    .entry(key)
                    .and_modify(|hits| hits.insert(locus))
                    .or_insert(LocusHits::One(locus));
            }),
            Self::Long(index) => scan_kmers(sequence, k, 1, true, |key, _| {
                index
                    .entry(key)
                    .and_modify(|hits| hits.insert(locus))
                    .or_insert(LocusHits::One(locus));
            }),
        }
    }

    fn bloom(&self) -> BlockedBloom {
        match self {
            Self::Short(index) => {
                let mut bloom = BlockedBloom::for_keys(index.len());
                for &key in index.keys() {
                    bloom.insert_hash(BlockedBloom::mix(key));
                }
                bloom
            }
            Self::Long(index) => {
                let mut bloom = BlockedBloom::for_keys(index.len());
                for &key in index.keys() {
                    bloom.insert_hash(BlockedBloom::mix(key as u64 ^ (key >> 64) as u64));
                }
                bloom
            }
        }
    }

    fn scan(
        &self,
        sequence: &[u8],
        k: usize,
        step: usize,
        bloom: &BlockedBloom,
        mut visit: impl FnMut(Option<&LocusHits>, bool),
    ) {
        match self {
            Self::Short(index) => scan_kmers_u64(sequence, k, step, true, |key, _| {
                let may_contain = bloom.may_contain_hash(BlockedBloom::mix(key));
                visit(may_contain.then(|| index.get(&key)).flatten(), !may_contain);
            }),
            Self::Long(index) => scan_kmers(sequence, k, step, true, |key, _| {
                let may_contain =
                    bloom.may_contain_hash(BlockedBloom::mix(key as u64 ^ (key >> 64) as u64));
                visit(may_contain.then(|| index.get(&key)).flatten(), !may_contain);
            }),
        }
    }
}

fn scan_recruit_index(
    index: &RecruitIndex,
    bloom: &BlockedBloom,
    k: usize,
    sequence: &[u8],
    step: usize,
    hits: &mut RecruitScratch,
    profile: Option<&mut IndexProfile>,
) {
    let mut profile = profile;
    index.scan(sequence, k, step, bloom, |loci, bloom_rejected| {
        if let Some(profile) = profile.as_deref_mut() {
            profile.recruit_probes += 1;
            profile.recruit_bloom_rejected += u64::from(bloom_rejected);
            profile.recruit_hits += u64::from(loci.is_some());
        }
        let Some(loci) = loci else { return };
        for &locus in loci.values() {
            hits.insert(locus);
        }
    });
}

#[derive(Clone, Debug)]
pub struct OrientedReference {
    pub locus: LocusId,
    pub strand: u8,
    pub bases: Vec<u8>,
}

pub struct UceIndex {
    pub k: usize,
    pub run_k: usize,
    pub loci: Vec<Locus>,
    pub references: Vec<OrientedReference>,
    recruit: RecruitIndex,
    recruit_bloom: BlockedBloom,
    panel_recruit: Option<RecruitIndex>,
    panel_recruit_bloom: Option<BlockedBloom>,
    exact: Vec<LocusMemIndex>,
}

fn stripped_extension(path: &Path) -> Option<String> {
    let base = if path
        .extension()
        .and_then(|v| v.to_str())
        .is_some_and(|v| v.eq_ignore_ascii_case("gz"))
    {
        PathBuf::from(path.file_stem()?)
    } else {
        path.to_path_buf()
    };
    base.extension()?.to_str().map(|v| v.to_ascii_lowercase())
}

fn reference_name(path: &Path) -> Result<String, String> {
    let base = if path
        .extension()
        .and_then(|v| v.to_str())
        .is_some_and(|v| v.eq_ignore_ascii_case("gz"))
    {
        PathBuf::from(
            path.file_stem()
                .ok_or_else(|| "invalid reference path".to_string())?,
        )
    } else {
        path.to_path_buf()
    };
    base.file_stem()
        .and_then(|v| v.to_str())
        .map(str::to_owned)
        .ok_or_else(|| format!("invalid reference name: {}", path.display()))
}

fn reference_paths(reference: &Path) -> Result<Vec<PathBuf>, String> {
    let mut paths = Vec::new();
    if reference.is_dir() {
        for entry in fs::read_dir(reference).map_err(|e| e.to_string())? {
            let path = entry.map_err(|e| e.to_string())?.path();
            if path.is_file()
                && matches!(
                    stripped_extension(&path).as_deref(),
                    Some("fa" | "fas" | "fasta")
                )
            {
                paths.push(path);
            }
        }
    } else if reference.is_file() {
        paths.push(reference.to_path_buf());
    }
    paths.sort();
    if paths.is_empty() {
        return Err("no reference FASTA file found".to_string());
    }
    Ok(paths)
}

pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&base| match code(base) {
            Some(value) => b"TGCA"[value as usize],
            None => base,
        })
        .collect()
}

pub fn scan_kmers(
    sequence: &[u8],
    k: usize,
    step: usize,
    canonical: bool,
    mut visit: impl FnMut(u128, usize),
) {
    if k == 0 || k > 64 || sequence.len() < k {
        return;
    }
    let mask = if k == 64 {
        u128::MAX
    } else {
        (1_u128 << (2 * k)) - 1
    };
    let reverse_shift = 2 * (k - 1);
    let tail = sequence.len() - k;
    let (mut forward, mut reverse, mut valid, mut next_probe) = (0_u128, 0_u128, 0_usize, 0_usize);
    for (end, &base) in sequence.iter().enumerate() {
        if let Some(value) = code(base) {
            forward = ((forward << 2) | value as u128) & mask;
            reverse = (reverse >> 2) | (((3 - value) as u128) << reverse_shift);
            valid += 1;
        } else {
            forward = 0;
            reverse = 0;
            valid = 0;
        }
        if end + 1 < k {
            continue;
        }
        let start = end + 1 - k;
        let sampled = start == next_probe;
        if sampled {
            next_probe = next_probe.saturating_add(step.max(1));
        }
        if valid >= k && (sampled || start == tail) {
            visit(
                if canonical {
                    forward.min(reverse)
                } else {
                    forward
                },
                start,
            );
        }
    }
}

pub fn scan_kmers_u64(
    sequence: &[u8],
    k: usize,
    step: usize,
    canonical: bool,
    mut visit: impl FnMut(u64, usize),
) {
    debug_assert!((1..=32).contains(&k));
    if k == 0 || k > 32 || sequence.len() < k {
        return;
    }
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1_u64 << (2 * k)) - 1
    };
    let reverse_shift = 2 * (k - 1);
    let tail = sequence.len() - k;
    let (mut forward, mut reverse, mut valid, mut next_probe) = (0_u64, 0_u64, 0_usize, 0_usize);
    for (end, &base) in sequence.iter().enumerate() {
        if let Some(value) = code(base) {
            forward = ((forward << 2) | value as u64) & mask;
            reverse = (reverse >> 2) | (((3 - value) as u64) << reverse_shift);
            valid += 1;
        } else {
            forward = 0;
            reverse = 0;
            valid = 0;
        }
        if end + 1 < k {
            continue;
        }
        let start = end + 1 - k;
        let sampled = start == next_probe;
        if sampled {
            next_probe = next_probe.saturating_add(step.max(1));
        }
        if valid >= k && (sampled || start == tail) {
            visit(
                if canonical {
                    forward.min(reverse)
                } else {
                    forward
                },
                start,
            );
        }
    }
}

pub fn default_verify_k(recruit_k: usize) -> usize {
    std::cmp::max(recruit_k / 2, recruit_k.saturating_sub(13)) | 1
}

impl UceIndex {
    pub fn build(reference: &Path, k: usize) -> Result<Self, String> {
        Self::build_split(reference, reference, k)
    }

    pub fn build_split(
        recruit_reference: &Path,
        verify_reference: &Path,
        k: usize,
    ) -> Result<Self, String> {
        Self::build_split_with_verify_k(recruit_reference, verify_reference, k, default_verify_k(k))
    }

    pub fn build_split_with_verify_k(
        recruit_reference: &Path,
        verify_reference: &Path,
        k: usize,
        verify_k: usize,
    ) -> Result<Self, String> {
        if !(1..=64).contains(&k) {
            return Err("UCEFilter currently supports k-mer sizes 1..=64".to_string());
        }
        if !(1..=64).contains(&verify_k) {
            return Err("UCEFilter currently supports verification k-mer sizes 1..=64".to_string());
        }
        let mut index = Self {
            k,
            run_k: verify_k,
            loci: Vec::new(),
            references: Vec::new(),
            recruit: RecruitIndex::new(k),
            recruit_bloom: BlockedBloom::for_keys(1),
            panel_recruit: None,
            panel_recruit_bloom: None,
            exact: Vec::new(),
        };
        for path in reference_paths(verify_reference)? {
            let locus = index.loci.len() as LocusId;
            let name = reference_name(&path)?;
            let mut reader =
                FastxReader::open(&path, FastxFormat::Fasta).map_err(|e| e.to_string())?;
            let mut originals = Vec::new();
            while let Some(record) = reader.next_record().map_err(|e| e.to_string())? {
                if record.sequence.len() >= k {
                    originals.push(record.sequence);
                }
            }
            let max_len = originals.iter().map(Vec::len).max().unwrap_or(0) as f64;
            let effective_length = if originals.is_empty() {
                0.0
            } else {
                (max_len * ((originals.len() as f64).log10() + 1.0)).trunc()
            };
            index.loci.push(Locus {
                name,
                effective_length,
            });
            for original in originals {
                index.recruit.insert_sequence(&original, k, locus);
                index.add_oriented(locus, 1, original.clone());
                index.add_oriented(locus, 2, reverse_complement(&original));
            }
        }
        if recruit_reference != verify_reference {
            let panel_recruit = std::mem::replace(&mut index.recruit, RecruitIndex::new(k));
            index.panel_recruit_bloom = Some(panel_recruit.bloom());
            index.panel_recruit = Some(panel_recruit);
            let locus_by_name: AHashMap<_, _> = index
                .loci
                .iter()
                .enumerate()
                .map(|(id, locus)| (locus.name.clone(), id as LocusId))
                .collect();
            for path in reference_paths(recruit_reference)? {
                let name = reference_name(&path)?;
                let Some(&locus) = locus_by_name.get(&name) else {
                    continue;
                };
                let mut reader =
                    FastxReader::open(&path, FastxFormat::Fasta).map_err(|e| e.to_string())?;
                while let Some(record) = reader.next_record().map_err(|e| e.to_string())? {
                    index.recruit.insert_sequence(&record.sequence, k, locus);
                }
            }
        }
        index.recruit_bloom = index.recruit.bloom();
        index.exact = build_locus_indexes(&index.references, index.loci.len(), verify_k)?;
        Ok(index)
    }

    fn add_oriented(&mut self, locus: LocusId, strand: u8, bases: Vec<u8>) {
        self.references.push(OrientedReference {
            locus,
            strand,
            bases,
        });
    }

    pub fn recruit(
        &self,
        sequence: &[u8],
        step: usize,
        hits: &mut RecruitScratch,
        profile: Option<&mut IndexProfile>,
    ) {
        scan_recruit_index(
            &self.recruit,
            &self.recruit_bloom,
            self.k,
            sequence,
            step,
            hits,
            profile,
        );
    }

    pub fn requires_panel_recruitment_expansion(&self) -> bool {
        self.panel_recruit.is_some()
    }

    pub fn expand_recruitment_to_panel(
        &self,
        sequence: &[u8],
        step: usize,
        hits: &mut RecruitScratch,
        profile: Option<&mut IndexProfile>,
    ) {
        let (Some(recruit), Some(bloom)) = (&self.panel_recruit, &self.panel_recruit_bloom) else {
            return;
        };
        scan_recruit_index(recruit, bloom, self.k, sequence, step, hits, profile);
    }

    pub fn orientation_events(&self, sequence: &[u8], candidates: &[LocusId]) -> Vec<Vec<u8>> {
        let mut scratch = ReadEvidenceScratch::default();
        self.read_evidence(sequence, candidates, &mut scratch, None);
        (0..candidates.len())
            .map(|candidate| scratch.orientation(candidate).to_vec())
            .collect()
    }

    /// Queries each recruited locus directly. A run-k event contributes a
    /// direction only when it occurs exactly once in that locus; repeated
    /// windows still contribute to the MEM used for position and exact length.
    pub fn read_evidence(
        &self,
        read: &[u8],
        candidates: &[LocusId],
        result: &mut ReadEvidenceScratch,
        profile: Option<&mut IndexProfile>,
    ) {
        let windows = read.len().saturating_sub(self.run_k).saturating_add(1);
        result.reset(candidates.len(), windows);
        if !valid_dna(read) || windows == 0 || candidates.is_empty() {
            return;
        }
        let mut profile = profile;
        read_codes_into(read, &mut result.read_codes);
        for (slot, &locus) in candidates.iter().enumerate() {
            result.occurrence_counts.fill(0);
            result.strand_masks.fill(0);
            let mut query_profile = MemQueryProfile::default();
            let best = self.exact[locus as usize].collect(
                &result.read_codes,
                self.run_k,
                &mut result.occurrence_counts,
                &mut result.strand_masks,
                &mut query_profile,
            );
            let start = slot * windows;
            for position in 0..windows {
                result.orientation_events[start + position] =
                    if result.occurrence_counts[position] == 1 {
                        result.strand_masks[position]
                    } else {
                        0
                    };
            }
            result.best_exact[slot] = best.map(exact_seed);
            if let Some(profile) = profile.as_deref_mut() {
                profile.exact_locus_queries += 1;
                profile.exact_index_queries += query_profile.index_queries;
                profile.exact_run_windows += query_profile.run_windows;
                profile.exact_matching_windows += query_profile.matching_windows;
                profile.mem_starts += query_profile.mem_starts;
                profile.mem_bases += query_profile.mem_bases;
            }
        }
    }

    pub fn best_exact(&self, read: &[u8], locus: LocusId) -> Option<ExactSeed> {
        if !valid_dna(read) || read.len() < self.run_k {
            return None;
        }
        let windows = read.len() - self.run_k + 1;
        let mut occurrence_counts = vec![0_u32; windows];
        let mut strand_masks = vec![0_u8; windows];
        let mut read_codes = Vec::with_capacity(read.len());
        read_codes_into(read, &mut read_codes);
        self.exact[locus as usize]
            .collect(
                &read_codes,
                self.run_k,
                &mut occurrence_counts,
                &mut strand_masks,
                &mut MemQueryProfile::default(),
            )
            .map(exact_seed)
    }

    pub fn max_exact(&self, read: &[u8], locus: LocusId) -> usize {
        self.best_exact(read, locus).map_or(0, ExactSeed::len)
    }

    pub fn exact_index_symbols(&self) -> usize {
        self.exact.iter().map(LocusMemIndex::symbols).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_TEST_DIR: AtomicU64 = AtomicU64::new(0);

    fn test_index(k: usize, reference: &[u8]) -> UceIndex {
        let root = std::env::temp_dir().join(format!(
            "uce-filter-index-{}-{k}-{}",
            std::process::id(),
            NEXT_TEST_DIR.fetch_add(1, Ordering::Relaxed),
        ));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let mut out = fs::File::create(root.join("locus.fa")).unwrap();
        out.write_all(b">ref\n").unwrap();
        out.write_all(reference).unwrap();
        out.write_all(b"\n").unwrap();
        drop(out);
        let index = UceIndex::build(&root, k).unwrap();
        fs::remove_dir_all(root).unwrap();
        index
    }

    fn brute_hit(read: &[u8], reference: &[u8], k: usize) -> bool {
        if !valid_dna(read) || read.len() < k || reference.len() < k {
            return false;
        }
        let reverse = reverse_complement(reference);
        read.windows(k).any(|window| {
            reference.windows(k).any(|candidate| candidate == window)
                || reverse.windows(k).any(|candidate| candidate == window)
        })
    }

    fn brute_max_exact(read: &[u8], reference: &[u8]) -> usize {
        let reverse = reverse_complement(reference);
        [reference, reverse.as_slice()]
            .into_iter()
            .flat_map(|oriented| {
                (0..read.len()).flat_map(move |read_start| {
                    (0..oriented.len()).map(move |reference_start| {
                        read[read_start..]
                            .iter()
                            .zip(&oriented[reference_start..])
                            .take_while(|(left, right)| code(**left) == code(**right))
                            .count()
                    })
                })
            })
            .max()
            .unwrap_or(0)
    }

    #[test]
    fn recruitment_and_verification_kmers_are_independent() {
        let root = std::env::temp_dir().join(format!(
            "uce-filter-independent-k-{}-{}",
            std::process::id(),
            NEXT_TEST_DIR.fetch_add(1, Ordering::Relaxed),
        ));
        fs::create_dir_all(&root).unwrap();
        let reference = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCG";
        let mut out = fs::File::create(root.join("locus.fa")).unwrap();
        writeln!(out, ">ref").unwrap();
        out.write_all(reference).unwrap();
        writeln!(out).unwrap();
        drop(out);
        let sensitive = UceIndex::build_split_with_verify_k(&root, &root, 21, 19).unwrap();
        assert_eq!(sensitive.k, 21);
        assert_eq!(sensitive.run_k, 19);
        assert_eq!(default_verify_k(21), 11);
        assert_eq!(default_verify_k(31), 19);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn fm_maximum_exact_matches_brute_force_across_mutated_reads() {
        let mut state = 0x9e37_79b9_u32;
        let mut reference = Vec::with_capacity(180);
        for _ in 0..180 {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            reference.push(b"ACGT"[((state >> 29) & 3) as usize]);
        }
        let index = test_index(16, &reference);
        for iteration in 0..64 {
            let start = iteration * 2 % 110;
            let mut read = reference[start..start + 60].to_vec();
            let mutation = (iteration * 7 + 11) % read.len();
            read[mutation] = b"ACGT"[(iteration + 1) % 4];
            if iteration % 2 == 1 {
                read = reverse_complement(&read);
            }
            let brute = brute_max_exact(&read, &reference);
            let expected = if brute >= index.run_k { brute } else { 0 };
            assert_eq!(index.max_exact(&read, 0), expected, "iteration={iteration}");
        }
    }

    #[test]
    fn fm_exact_matching_is_case_insensitive() {
        let reference = b"acgttgcaacgattcggtaccatgcaagttcg";
        let read = b"ACGTTGCAACGATTCGGTACC";
        let index = test_index(16, reference);
        assert_eq!(index.max_exact(read, 0), read.len());
    }

    #[test]
    fn maximum_exact_reproduces_all_thresholds() {
        let reference = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCGATCGGATCCGTAACCGGTT";
        let read = b"TTTTGCAACGATTCGGTACCAGGG";
        let index = test_index(16, reference);
        let maximum = index.max_exact(read, 0);
        for k in [1, 5, 16, 17, 20, 24] {
            assert_eq!(
                brute_hit(read, reference, k),
                maximum >= k,
                "k={k}, M={maximum}"
            );
        }
    }

    #[test]
    fn ambiguous_read_has_no_legacy_exact_match() {
        let reference = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCG";
        let index = test_index(16, reference);
        assert_eq!(index.max_exact(b"ACGTTGCAACGATTCNGTACC", 0), 0);
    }

    #[test]
    fn maximum_exact_handles_long_and_reverse_strand_thresholds() {
        let reference =
            b"ACGTTGCAACGATTCGGTACCATGCAAGTTCGATCGGATCCGTAACCGGTTAGCTACGATGCTAGGCTTACCGATGGCATTCG";
        let read = reverse_complement(&reference[8..80]);
        let index = test_index(16, reference);
        let maximum = index.max_exact(&read, 0);
        for k in [16, 31, 32, 33, 63, 64, 67] {
            assert_eq!(
                brute_hit(&read, reference, k),
                maximum >= k,
                "k={k}, M={maximum}"
            );
        }
    }

    #[test]
    fn unique_run_windows_cast_orientation_votes() {
        let reference = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCG";
        let read = &reference[3..28];
        let index = test_index(16, reference);
        let events = index.orientation_events(read, &[0]);
        assert_eq!(events[0].len(), read.len() - index.run_k + 1);
        assert!(events[0].iter().all(|&event| event == 1));
    }

    #[test]
    fn repeated_run_windows_do_not_vote_but_still_form_a_mem() {
        let motif = b"ACGTTGCAACGATTCG";
        let reference = b"ACGTTGCAACGATTCGTTAAACGTTGCAACGATTCG";
        let index = test_index(16, reference);
        let events = index.orientation_events(motif, &[0]);
        assert!(events[0].iter().all(|&event| event == 0));
        assert_eq!(index.max_exact(motif, 0), motif.len());
    }

    #[test]
    fn maximum_exact_restarts_after_a_mismatch_on_the_same_diagonal() {
        let reference = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCGATCGGATCCGTAACCGGTT";
        let mut read = reference[4..49].to_vec();
        read[23] = if read[23] == b'A' { b'C' } else { b'A' };
        let index = test_index(16, reference);
        let maximum = index.max_exact(&read, 0);
        for k in [16, 20, 23, 24, 31] {
            assert_eq!(
                brute_hit(&read, reference, k),
                maximum >= k,
                "k={k}, M={maximum}"
            );
        }
    }
}
