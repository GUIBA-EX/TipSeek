use crate::hash::{HashMap, HashSet};
use crate::model::{
    Args, AssemblyMode, ContigRecord, KmerInfo, Node, PathContig, PathStrategy, ReadSupport,
    RefKmer, SideContig,
};
use crate::seq::{
    bits_base, decode_kmer, for_each_kmer, kmer_mask, median, quartiles, reverse_complement,
    reverse_complement_kmer, valid_runs,
};
use std::cmp::Ordering;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

pub fn add_read_slices(reads: &mut HashMap<Vec<u8>, u64>, sequences: &[Vec<u8>], slice_len: usize) {
    // 从 read 中间掐一段当身份证，正反两面都备着，后面核验 contig 用。
    for sequence in sequences {
        if sequence.len() < slice_len {
            continue;
        }
        let start = (sequence.len() - slice_len) / 2;
        let forward = sequence[start..start + slice_len].to_vec();
        let reverse_sequence = reverse_complement(sequence);
        let reverse = reverse_sequence[start..start + slice_len].to_vec();
        *reads.entry(forward.clone()).or_default() += 1;
        if reverse != forward {
            *reads.entry(reverse).or_default() += 1;
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct KmerCount {
    kmer: u128,
    depth: i64,
    fragment_support: u32,
}

pub type BranchSupport = HashMap<(u128, u128), u32>;

#[derive(Default)]
pub struct SortedKmerCounts {
    // Binary-carry runs keep memory proportional to unique k-mers while
    // avoiding a random HashMap update for every input chunk.
    levels: Vec<Option<Vec<KmerCount>>>,
}

fn merge_kmer_counts(left: Vec<KmerCount>, right: Vec<KmerCount>) -> Vec<KmerCount> {
    let mut merged = Vec::with_capacity(left.len() + right.len());
    let (mut left_index, mut right_index) = (0, 0);
    while left_index < left.len() && right_index < right.len() {
        match left[left_index].kmer.cmp(&right[right_index].kmer) {
            Ordering::Less => {
                merged.push(left[left_index]);
                left_index += 1;
            }
            Ordering::Greater => {
                merged.push(right[right_index]);
                right_index += 1;
            }
            Ordering::Equal => {
                merged.push(KmerCount {
                    kmer: left[left_index].kmer,
                    depth: left[left_index].depth + right[right_index].depth,
                    fragment_support: left[left_index]
                        .fragment_support
                        .saturating_add(right[right_index].fragment_support),
                });
                left_index += 1;
                right_index += 1;
            }
        }
    }
    merged.extend_from_slice(&left[left_index..]);
    merged.extend_from_slice(&right[right_index..]);
    merged
}

impl SortedKmerCounts {
    pub fn push(&mut self, mut run: Vec<KmerCount>) {
        if run.is_empty() {
            return;
        }
        let mut level = 0;
        loop {
            if self.levels.len() == level {
                self.levels.push(None);
            }
            if let Some(existing) = self.levels[level].take() {
                run = merge_kmer_counts(existing, run);
                level += 1;
            } else {
                self.levels[level] = Some(run);
                return;
            }
        }
    }

    pub fn into_counts(mut self) -> Vec<KmerCount> {
        let mut result = Vec::new();
        for run in self.levels.drain(..).flatten() {
            result = merge_kmer_counts(result, run);
        }
        result
    }
}

/// Sort key for k-mer counting. K-mers of k <= 32 fit in `u64`; sorting them
/// at that width gives the same order as `u128` with half the memory traffic.
trait CountKey: Copy + Ord + Send {
    fn from_kmer(kmer: u128) -> Self;
    fn into_kmer(self) -> u128;
}

impl CountKey for u64 {
    fn from_kmer(kmer: u128) -> Self {
        debug_assert!(kmer <= u128::from(u64::MAX));
        kmer as u64
    }

    fn into_kmer(self) -> u128 {
        u128::from(self)
    }
}

impl CountKey for u128 {
    fn from_kmer(kmer: u128) -> Self {
        kmer
    }

    fn into_kmer(self) -> u128 {
        self
    }
}

fn append_read_kmers<K: CountKey>(sequence: &[u8], k: usize, physical: &mut Vec<K>) {
    physical.clear();
    let expected = sequence.len().saturating_sub(k - 1).saturating_mul(2);
    if physical.capacity() < expected {
        physical.reserve(expected - physical.capacity());
    }
    for (_, run) in valid_runs(sequence) {
        for_each_kmer(run, k, |_, forward, reverse| {
            physical.push(K::from_kmer(forward));
            physical.push(K::from_kmer(reverse));
        });
    }
    physical.sort_unstable();
    physical.dedup();
}

fn counts_from_observations<K: CountKey>(
    mut observations: Vec<K>,
    mut fragments: Vec<K>,
) -> Vec<KmerCount> {
    observations.sort_unstable();
    fragments.sort_unstable();

    let mut counts = Vec::with_capacity(observations.len() / 2);
    let (mut observation_index, mut fragment_index) = (0, 0);
    while observation_index < observations.len() {
        let kmer = observations[observation_index];
        let mut depth = 0_i64;
        while observation_index < observations.len() && observations[observation_index] == kmer {
            depth += 1;
            observation_index += 1;
        }
        while fragment_index < fragments.len() && fragments[fragment_index] < kmer {
            fragment_index += 1;
        }
        let mut fragment_support = 0_u32;
        while fragment_index < fragments.len() && fragments[fragment_index] == kmer {
            fragment_support = fragment_support.saturating_add(1);
            fragment_index += 1;
        }
        counts.push(KmerCount {
            kmer: kmer.into_kmer(),
            depth,
            fragment_support,
        });
    }
    counts
}

/// Fragments (or unpaired reads) a k-mer worker claims at a time. Workers
/// claim blocks until the chunk is exhausted, so faster cores take a larger
/// share on hybrid processors instead of waiting for a fixed slice on the
/// slowest core. Counts are merged by k-mer, so the split never changes them.
const COUNT_BLOCK_UNITS: usize = 64;

/// Per-worker k-mer observations, reused across claimed blocks.
struct CountBuffers<K> {
    observations: Vec<K>,
    fragments: Vec<K>,
    physical: Vec<K>,
    first_mate: Vec<K>,
    second_mate: Vec<K>,
}

impl<K: CountKey> CountBuffers<K> {
    fn with_capacity(observations: usize) -> Self {
        Self {
            observations: Vec::with_capacity(observations),
            fragments: Vec::with_capacity(observations / 4),
            physical: Vec::new(),
            first_mate: Vec::new(),
            second_mate: Vec::new(),
        }
    }

    fn add(
        &mut self,
        part: &[Vec<u8>],
        k: usize,
        paired_fragments: bool,
        reference: Option<&HashMap<u128, RefKmer>>,
    ) {
        if paired_fragments {
            for mates in part.chunks(2) {
                if mates.len() != 2 {
                    append_read_kmers(&mates[0], k, &mut self.physical);
                    self.observations.extend_from_slice(&self.physical);
                    continue;
                }
                append_read_kmers(&mates[0], k, &mut self.first_mate);
                self.observations.extend_from_slice(&self.first_mate);
                append_read_kmers(&mates[1], k, &mut self.second_mate);
                self.observations.extend_from_slice(&self.second_mate);

                // A retained mate is useful flank evidence only
                // when its *opposite* mate anchors in reference
                // sequence.  Counting every k-mer in every pair
                // would otherwise preserve recurrent errors.
                if let Some(reference) = reference {
                    let first_anchor = self
                        .first_mate
                        .iter()
                        .any(|kmer| reference.contains_key(&kmer.into_kmer()));
                    let second_anchor = self
                        .second_mate
                        .iter()
                        .any(|kmer| reference.contains_key(&kmer.into_kmer()));
                    if first_anchor && !second_anchor {
                        self.fragments.extend_from_slice(&self.second_mate);
                    } else if second_anchor && !first_anchor {
                        self.fragments.extend_from_slice(&self.first_mate);
                    }
                }
            }
        } else {
            for sequence in part {
                append_read_kmers(sequence, k, &mut self.physical);
                self.observations.extend_from_slice(&self.physical);
            }
        }
    }

    fn into_counts(self) -> Vec<KmerCount> {
        counts_from_observations(self.observations, self.fragments)
    }
}

/// Counts the blocks this worker claims from `next_block`.
fn count_claimed_blocks<K: CountKey>(
    sequences: &[Vec<u8>],
    k: usize,
    paired_fragments: bool,
    reference: Option<&HashMap<u128, RefKmer>>,
    next_block: &AtomicUsize,
    expected_observations: usize,
) -> Vec<KmerCount> {
    let slots = if paired_fragments { 2 } else { 1 };
    let block_slots = COUNT_BLOCK_UNITS * slots;
    let mut buffers = CountBuffers::<K>::with_capacity(expected_observations);
    loop {
        let start = next_block
            .fetch_add(1, AtomicOrdering::Relaxed)
            .saturating_mul(block_slots);
        if start >= sequences.len() {
            break;
        }
        let end = (start + block_slots).min(sequences.len());
        buffers.add(&sequences[start..end], k, paired_fragments, reference);
    }
    buffers.into_counts()
}

pub fn count_assemble_chunk_parallel(
    sequences: &[Vec<u8>],
    k: usize,
    threads: usize,
    paired_fragments: bool,
    reference: Option<&HashMap<u128, RefKmer>>,
) -> Vec<KmerCount> {
    // Each read contributes a k-mer once to depth.  A PE fragment contributes
    // it once to fragment support only when both mate slots are present.
    if sequences.is_empty() {
        return Vec::new();
    }
    if 2 * k <= u64::BITS as usize {
        count_chunk::<u64>(sequences, k, threads, paired_fragments, reference)
    } else {
        count_chunk::<u128>(sequences, k, threads, paired_fragments, reference)
    }
}

fn count_chunk<K: CountKey>(
    sequences: &[Vec<u8>],
    k: usize,
    threads: usize,
    paired_fragments: bool,
    reference: Option<&HashMap<u128, RefKmer>>,
) -> Vec<KmerCount> {
    let units = if paired_fragments {
        sequences.len().div_ceil(2)
    } else {
        sequences.len()
    };
    let blocks = units.div_ceil(COUNT_BLOCK_UNITS);
    let workers = threads.max(1).min(blocks);
    let expected_observations = sequences
        .iter()
        .map(|sequence| sequence.len().saturating_sub(k - 1).saturating_mul(2))
        .sum::<usize>()
        / workers;
    let next_block = AtomicUsize::new(0);
    if workers == 1 {
        return count_claimed_blocks::<K>(
            sequences,
            k,
            paired_fragments,
            reference,
            &next_block,
            expected_observations,
        );
    }
    let mut partials: Vec<Vec<KmerCount>> = Vec::with_capacity(workers);
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                let next_block = &next_block;
                scope.spawn(move || {
                    gm2_tools::platform::prefer_performance_cores();
                    count_claimed_blocks::<K>(
                        sequences,
                        k,
                        paired_fragments,
                        reference,
                        next_block,
                        expected_observations,
                    )
                })
            })
            .collect();
        for handle in handles {
            partials.push(handle.join().expect("k-mer worker panicked"));
        }
    });

    partials.into_iter().fold(Vec::new(), merge_kmer_counts)
}

pub fn build_graph_from_counts(
    counts: Vec<KmerCount>,
    reference: &HashMap<u128, RefKmer>,
    error_limit: u32,
    min_fragment_support: u32,
) -> HashMap<u128, KmerInfo> {
    let mut graph = HashMap::with_capacity(counts.len());
    for count in counts {
        let reference_info = reference.get(&count.kmer);
        let pe_flank_rescue = error_limit > 0
            && count.depth <= error_limit as i64
            && reference_info.is_none()
            && min_fragment_support > 0
            && count.fragment_support >= min_fragment_support;
        if error_limit > 0
            && count.depth <= error_limit as i64
            && reference_info.is_none()
            && !pe_flank_rescue
        {
            continue;
        }
        let value = reference_info.map_or(
            KmerInfo {
                depth: count.depth,
                position: 1023,
                is_reverse: true,
                reference_weight: 0,
                fragment_support: count.fragment_support,
                pe_flank_rescue,
            },
            |r| KmerInfo {
                depth: count.depth,
                position: if r.is_reverse {
                    1000 - r.position
                } else {
                    r.position
                },
                is_reverse: r.is_reverse,
                reference_weight: r.depth as i64,
                fragment_support: count.fragment_support,
                pe_flank_rescue,
            },
        );
        graph.insert(count.kmer, value);
    }
    graph
}

#[cfg(test)]
pub fn build_assemble_dictionary(
    sequences: &[Vec<u8>],
    k: usize,
    reference: &HashMap<u128, RefKmer>,
) -> HashMap<u128, KmerInfo> {
    let mut graph: HashMap<u128, KmerInfo> = HashMap::new();
    for sequence in sequences {
        let mut physical_read_kmers = HashSet::new();
        for (_, run) in valid_runs(sequence) {
            for_each_kmer(run, k, |_, forward, reverse| {
                physical_read_kmers.insert(forward);
                physical_read_kmers.insert(reverse);
            });
        }
        for kmer in physical_read_kmers {
            if let Some(info) = graph.get_mut(&kmer) {
                info.depth += 1;
                continue;
            }
            if let Some(reference_info) = reference.get(&kmer) {
                let position = if reference_info.is_reverse {
                    1000 - reference_info.position
                } else {
                    reference_info.position
                };
                graph.insert(
                    kmer,
                    KmerInfo {
                        depth: 1,
                        position,
                        is_reverse: reference_info.is_reverse,
                        reference_weight: reference_info.depth as i64,
                        fragment_support: 0,
                        pe_flank_rescue: false,
                    },
                );
            } else {
                graph.insert(
                    kmer,
                    KmerInfo {
                        depth: 1,
                        position: 1023,
                        is_reverse: true,
                        reference_weight: 0,
                        fragment_support: 0,
                        pe_flank_rescue: false,
                    },
                );
            }
        }
    }
    graph
}

pub fn filter_and_weight_graph(
    graph: &mut HashMap<u128, KmerInfo>,
    error_limit: u32,
    reference_count: usize,
    min_fragment_support: u32,
) {
    // 低频毛刺先踢掉；参考上的节点留个权重，走岔路时优先认熟路。
    if error_limit > 0 {
        graph.retain(|_, value| {
            value.depth > error_limit as i64
                || value.reference_weight > 0
                || (min_fragment_support > 0 && value.fragment_support >= min_fragment_support)
        });
    }
    if graph.is_empty() {
        return;
    }

    let mut depths: Vec<i64> = graph
        .values()
        .filter(|value| !value.pe_flank_rescue)
        .map(|value| value.depth)
        .collect();
    if depths.is_empty() {
        depths.extend(graph.values().map(|value| value.depth));
    }
    let (q1, _, q3, _) = quartiles(&depths);
    let depth_upper = (((q3 - q1) * 1.5) + q3) as i64;
    for value in graph.values_mut() {
        if value.reference_weight != 0 {
            value.reference_weight = if value.depth > error_limit as i64 {
                (reference_count as i64 * depth_upper
                    / ((value.reference_weight - reference_count as i64).abs() + 1))
                    + 1
            } else {
                1
            };
        }
        value.depth = value.depth.min(depth_upper);
    }
}

// 只为真实分支保存 PE 证据；线性边无需额外索引，避免把整张图再复制一遍。
pub fn branch_edges(graph: &HashMap<u128, KmerInfo>, k: usize) -> BranchSupport {
    let suffix_mask = kmer_mask(k - 1);
    let mut support = BranchSupport::new();
    for current in graph.keys().copied() {
        let prefix = (current & suffix_mask) << 2;
        let candidates: Vec<u128> = (0..4_u128)
            .map(|base| prefix | base)
            .filter(|candidate| graph.contains_key(candidate))
            .collect();
        if candidates.len() > 1 {
            for candidate in candidates {
                support.insert((current, candidate), 0);
            }
        }
    }
    support
}

pub fn add_pe_branch_support(
    sequences: &[Vec<u8>],
    k: usize,
    reference: &HashMap<u128, RefKmer>,
    support: &mut BranchSupport,
) {
    if support.is_empty() {
        return;
    }
    for mates in sequences.chunks(2) {
        if mates.len() != 2 {
            continue;
        }
        let read_kmers = |sequence: &[u8]| {
            let mut kmers = Vec::new();
            for (run_start, run) in valid_runs(sequence) {
                for_each_kmer(run, k, |offset, forward, reverse| {
                    kmers.push((run_start + offset, forward, reverse));
                });
            }
            kmers
        };
        let first = read_kmers(&mates[0]);
        let second = read_kmers(&mates[1]);
        let first_anchor = first.iter().any(|(_, forward, reverse)| {
            reference.contains_key(forward) || reference.contains_key(reverse)
        });
        let second_anchor = second.iter().any(|(_, forward, reverse)| {
            reference.contains_key(forward) || reference.contains_key(reverse)
        });
        let mut fragment_edges = HashSet::new();
        for (kmers, mate_is_pe_anchored) in [(&first, second_anchor), (&second, first_anchor)] {
            if !mate_is_pe_anchored {
                continue;
            }
            for pair in kmers.windows(2) {
                if pair[1].0 != pair[0].0 + 1 {
                    continue;
                }
                for edge in [(pair[0].1, pair[1].1), (pair[1].2, pair[0].2)] {
                    if support.contains_key(&edge) {
                        fragment_edges.insert(edge);
                    }
                }
            }
        }
        for edge in fragment_edges {
            *support.get_mut(&edge).expect("branch edge exists") += 1;
        }
    }
}

// 从当前 k-mer 找能接上的下一跳，再按证据排个先后。
#[allow(clippy::too_many_arguments)]
fn outgoing(
    graph: &HashMap<u128, KmerInfo>,
    current: u128,
    k: usize,
    blocked: &HashSet<u128>,
    reverse_reuse_guard: &HashSet<u128>,
    discarded: Option<&HashSet<u128>>,
    reference_tie_break: bool,
    branch_support: Option<&BranchSupport>,
    reverse_reuse_reference_scale: f64,
) -> Vec<Node> {
    let suffix_mask = kmer_mask(k - 1);
    let prefix = (current & suffix_mask) << 2;
    let mut nodes = Vec::with_capacity(4);
    for base in 0..4_u128 {
        let candidate = prefix | base;
        if blocked.contains(&candidate) || discarded.is_some_and(|items| items.contains(&candidate))
        {
            continue;
        }
        if let Some(value) = graph.get(&candidate) {
            let reverse = reverse_complement_kmer(candidate, k);
            let reverse_reuse =
                blocked.contains(&reverse) || reverse_reuse_guard.contains(&reverse);
            let reference_weight = if reverse_reuse {
                (value.reference_weight as f64 * reverse_reuse_reference_scale).round() as i64
            } else {
                value.reference_weight
            };
            nodes.push(Node {
                kmer: candidate,
                position: value.position,
                weight: value.depth + reference_weight,
                pe_support: branch_support
                    .and_then(|support| support.get(&(current, candidate)).copied())
                    .unwrap_or(0),
            });
        }
    }
    if nodes.len() > 1 {
        let current_info = graph.get(&current).expect("current k-mer is in graph");
        for node in &mut nodes {
            let candidate = &graph[&node.kmer];
            // PE score is counted per independent fragment on this directed
            // branch edge; reference-adjacent transitions are preferred, while
            // abrupt coverage peaks are treated as repeat-like alternatives.
            if node.pe_support >= 2 {
                node.weight += i64::from(node.pe_support.min(4)) * 6;
            } else if candidate.pe_flank_rescue {
                // A node-level anchor is not sufficient at a branch: require
                // two fragments that directly traverse this particular edge.
                node.weight -= 16;
            }
            if current_info.position > 0
                && current_info.position < 1000
                && candidate.position > 0
                && candidate.position < 1000
            {
                let delta = (candidate.position - current_info.position).abs();
                if delta <= 2 {
                    node.weight += 8;
                } else if current_info.reference_weight > 0 && candidate.reference_weight > 0 {
                    node.weight -= 8;
                }
            }
            let surge_limit = (current_info.depth.max(2) * 4).max(8);
            if candidate.reference_weight == 0 && candidate.depth > surge_limit {
                node.weight -= (candidate.depth - surge_limit) / 2;
            }
        }
    }
    nodes.sort_by(|left, right| {
        right
            .weight
            .cmp(&left.weight)
            .then_with(|| {
                if reference_tie_break {
                    (right.position > 0).cmp(&(left.position > 0))
                } else {
                    Ordering::Equal
                }
            })
            .then_with(|| left.kmer.cmp(&right.kmer))
    });
    nodes
}

#[allow(clippy::too_many_arguments)]
pub fn walk_backbone(
    graph: &HashMap<u128, KmerInfo>,
    seed: u128,
    k: usize,
    iteration: usize,
    lookahead: usize,
    branch_support: Option<&BranchSupport>,
    reverse_reuse_guard: &HashSet<u128>,
    reverse_reuse_reference_scale: f64,
) -> (Vec<PathContig>, HashSet<u128>, Vec<i32>, i64) {
    // UCE 默认就认一条最靠谱的主干，岔路看几步再定，别在气泡里来回磨叽。
    let lookahead = lookahead.max(1);
    let mut path = vec![seed];
    let mut visited = HashSet::from([seed]);
    let mut discarded = HashSet::new();
    let mut positions = Vec::new();
    let mut result = PathContig::default();

    for _ in 0..iteration {
        let nodes = outgoing(
            graph,
            *path.last().expect("seed"),
            k,
            &visited,
            reverse_reuse_guard,
            Some(&discarded),
            true,
            branch_support,
            reverse_reuse_reference_scale,
        );
        if nodes.is_empty() {
            break;
        }
        let chosen = if nodes.len() == 1 {
            nodes[0]
        } else {
            let mut best_trace: Option<Vec<Node>> = None;
            for first in &nodes {
                let mut trace = Vec::new();
                let mut trace_seen = visited.clone();
                let mut node = *first;
                for _ in 0..lookahead {
                    if !trace_seen.insert(node.kmer) {
                        break;
                    }
                    trace.push(node);
                    let following = outgoing(
                        graph,
                        node.kmer,
                        k,
                        &trace_seen,
                        reverse_reuse_guard,
                        Some(&discarded),
                        true,
                        branch_support,
                        reverse_reuse_reference_scale,
                    );
                    let Some(next) = following.first() else {
                        break;
                    };
                    node = *next;
                }

                let replace = best_trace.as_ref().is_none_or(|best| {
                    let trace_sum: i64 = trace.iter().map(|node| node.weight).sum();
                    let best_sum: i64 = best.iter().map(|node| node.weight).sum();
                    (trace.len(), trace_sum, trace[0].weight)
                        > (best.len(), best_sum, best[0].weight)
                });
                if replace {
                    best_trace = Some(trace);
                }
            }
            let winner = best_trace.expect("outgoing branch has a trace")[0];
            for node in &nodes {
                if node.kmer != winner.kmer {
                    discarded.insert(node.kmer);
                }
            }
            winner
        };

        path.push(chosen.kmer);
        visited.insert(chosen.kmer);
        positions.push(chosen.position);
        result.weights.push(chosen.weight);
        result.bases.push((chosen.kmer & 3) as u8);

        // 后头没岔路就一口气走完，直道还一步一停那可太磨叽了。
        loop {
            let linear = outgoing(
                graph,
                *path.last().expect("path"),
                k,
                &visited,
                reverse_reuse_guard,
                Some(&discarded),
                true,
                branch_support,
                reverse_reuse_reference_scale,
            );
            if linear.len() != 1 || path.len() > iteration {
                break;
            }
            let node = linear[0];
            path.push(node.kmer);
            visited.insert(node.kmer);
            positions.push(node.position);
            result.weights.push(node.weight);
            result.bases.push((node.kmer & 3) as u8);
        }
    }

    let weight = result.weights.iter().sum();
    (vec![result], visited, positions, weight)
}

pub fn walk_search(
    graph: &HashMap<u128, KmerInfo>,
    seed: u128,
    k: usize,
    mut iteration: usize,
    reverse_reuse_guard: &HashSet<u128>,
    reverse_reuse_reference_scale: f64,
) -> (Vec<PathContig>, HashSet<u128>, Vec<i32>, i64) {
    // 兼容旧策略：遇岔路就压栈回溯，多找几条候选，不是 UCE 默认路子。
    let mut path = vec![seed];
    let mut path_set = HashSet::from([seed]);
    let mut all_visited = HashSet::from([seed]);
    let mut positions = Vec::new();
    let mut current = PathContig::default();
    let mut contigs = Vec::new();
    let mut stack: Vec<(Vec<Node>, usize)> = Vec::new();
    let mut node_distance = 0_usize;
    let mut best_weight = 0_i64;

    while iteration > 0 {
        let mut nodes = outgoing(
            graph,
            *path.last().expect("seed"),
            k,
            &path_set,
            reverse_reuse_guard,
            None,
            false,
            None,
            reverse_reuse_reference_scale,
        );
        if nodes.is_empty() {
            iteration -= 1;
            let weight: i64 = current.weights.iter().sum();
            best_weight = best_weight.max(weight);
            contigs.push(current.clone());

            for _ in 0..node_distance {
                if let Some(removed) = path.pop() {
                    path_set.remove(&removed);
                }
                current.weights.pop();
                current.bases.pop();
            }
            let Some((alternatives, distance)) = stack.pop() else {
                break;
            };
            nodes = alternatives;
            node_distance = distance;
        }

        if nodes.len() >= 2 {
            stack.push((nodes[1..].to_vec(), node_distance));
            node_distance = 0;
        }
        let chosen = nodes[0];
        path.push(chosen.kmer);
        path_set.insert(chosen.kmer);
        all_visited.insert(chosen.kmer);
        positions.push(chosen.position);
        current.weights.push(chosen.weight);
        current.bases.push((chosen.kmer & 3) as u8);
        node_distance += 1;
    }

    (contigs, all_visited, positions, best_weight)
}

// 在 contig 里找 read 身份片段，命中位置给支持度计算当凭据。
fn locate_read_slices<'a>(
    sequence: &'a [u8],
    slice_len: usize,
    reads: &HashMap<Vec<u8>, u64>,
) -> HashMap<&'a [u8], Option<usize>> {
    let mut matches = HashMap::new();
    if slice_len == 0 || sequence.len() < slice_len {
        return matches;
    }
    for position in 0..=sequence.len() - slice_len {
        let slice = &sequence[position..position + slice_len];
        if !reads.contains_key(slice) {
            continue;
        }
        if let Some(existing) = matches.get_mut(slice) {
            *existing = None;
        } else {
            matches.insert(slice, Some(position));
        }
    }
    matches
}

// 左右延伸分别核验，别一头靠谱另一头跑飞还硬拼成整条。
fn process_sides(
    mut contigs: Vec<PathContig>,
    max_weight: i64,
    slice_len: usize,
    reads: &HashMap<Vec<u8>, u64>,
    soft_boundary: usize,
    mode: AssemblyMode,
) -> Vec<SideContig> {
    if mode == AssemblyMode::Reference {
        for contig in &mut contigs {
            if contig.weights.len() > 2 {
                let (q1, _, _, _) = quartiles(&contig.weights);
                let cut_position = contig
                    .weights
                    .iter()
                    .rposition(|weight| *weight as f64 >= q1);
                if let Some(cut) = cut_position {
                    let retained = cut.saturating_add(soft_boundary).saturating_add(1);
                    if retained < contig.weights.len() {
                        contig.weights.truncate(retained);
                        contig.bases.truncate(retained);
                    }
                }
            }
        }
    }

    let minimum = match mode {
        AssemblyMode::Reference => max_weight >> 1,
        AssemblyMode::Uce => max_weight >> 2,
        AssemblyMode::Its2 => 1,
    };
    let mut processed = Vec::new();
    for contig in contigs {
        let weight: i64 = contig.weights.iter().sum();
        if weight <= minimum {
            continue;
        }
        let sequence: Vec<u8> = contig.bases.iter().map(|base| bits_base(*base)).collect();
        let matches = locate_read_slices(&sequence, slice_len, reads);
        let read_count = matches.keys().filter_map(|slice| reads.get(*slice)).sum();
        processed.push(SideContig {
            sequence,
            weight,
            read_count,
        });
    }

    if mode == AssemblyMode::Uce {
        processed.sort_by(|left, right| {
            right
                .sequence
                .len()
                .cmp(&left.sequence.len())
                .then_with(|| right.read_count.cmp(&left.read_count))
                .then_with(|| right.weight.cmp(&left.weight))
                .then_with(|| left.sequence.cmp(&right.sequence))
        });
    } else {
        processed.sort_by_key(|candidate| std::cmp::Reverse(candidate.read_count));
    }
    processed
}

// 读段支持、唯一支持和覆盖比例一块儿算，候选高低不能只看长度。
pub fn calculate_read_support(
    sequence: &[u8],
    slice_len: usize,
    reads: &HashMap<Vec<u8>, u64>,
) -> ReadSupport {
    let contig_len = sequence.len();
    let matches = locate_read_slices(sequence, slice_len, reads);
    let total_read_count = matches.keys().filter_map(|slice| reads.get(*slice)).sum();
    let mut unique_read_count = 0_u64;
    let mut multi_mapping_read_count = 0_u64;
    let mut intervals = Vec::new();

    for (slice, position) in matches {
        let count = reads.get(slice).copied().unwrap_or(0);
        if let Some(start) = position {
            unique_read_count += count;
            intervals.push((start, (start + slice_len).min(contig_len)));
        } else {
            multi_mapping_read_count += count;
        }
    }

    if intervals.is_empty() {
        return ReadSupport {
            total_read_count,
            multi_mapping_read_count,
            max_gap: contig_len,
            left_coord: contig_len,
            ..ReadSupport::default()
        };
    }

    intervals.sort_unstable();
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (start, end) in intervals {
        if merged
            .last()
            .is_none_or(|(_, previous_end)| start > *previous_end)
        {
            merged.push((start, end));
        } else if let Some((_, previous_end)) = merged.last_mut() {
            *previous_end = (*previous_end).max(end);
        }
    }

    let left_coord = merged[0].0;
    let right_coord = merged.last().expect("nonempty").1;
    let supported_extent = right_coord - left_coord;
    let supported_bases = merged.iter().map(|(start, end)| end - start).sum();
    let breadth = if contig_len > 0 {
        supported_bases as f64 / contig_len as f64
    } else {
        0.0
    };
    let mut max_gap = left_coord.max(contig_len - right_coord);
    for pair in merged.windows(2) {
        max_gap = max_gap.max(pair[1].0 - pair[0].1);
    }
    let left_extension = left_coord;
    let right_extension = contig_len - right_coord;
    let flank_balance = if supported_extent == 0 {
        0.0
    } else if left_extension == 0 && right_extension == 0 {
        1.0
    } else {
        left_extension.min(right_extension) as f64 / left_extension.max(right_extension) as f64
    };

    ReadSupport {
        total_read_count,
        unique_read_count,
        multi_mapping_read_count,
        supported_extent,
        supported_bases,
        breadth,
        max_gap,
        flank_balance,
        left_coord,
        right_coord,
    }
}

// 深度的中位数和离散度用来揪重复区，别让高拷贝把结果忽悠了。
fn depth_stats(sequence: &[u8], k: usize, graph: &HashMap<u128, KmerInfo>) -> (f64, f64, f64) {
    if sequence.len() < k {
        return (0.0, 0.0, 0.0);
    }
    let mut counts = Vec::with_capacity(sequence.len() - k + 1);
    for_each_kmer(sequence, k, |_, forward, _| {
        counts.push(graph.get(&forward).map_or(0, |value| value.depth));
    });
    let median_depth = median(&counts);
    if counts.is_empty() || median_depth <= 0.0 {
        return (median_depth, 0.0, 0.0);
    }
    let mean = counts.iter().sum::<i64>() as f64 / counts.len() as f64;
    let variance = counts
        .iter()
        .map(|count| (*count as f64 - mean).powi(2))
        .sum::<f64>()
        / counts.len() as f64;
    let cv = if mean > 0.0 {
        variance.sqrt() / mean
    } else {
        0.0
    };
    let maximum = counts.iter().copied().max().unwrap_or(0) as f64;
    (median_depth, cv, maximum / median_depth)
}

// UCE 护栏集中在这儿：太长、太稀、深度乱，都得给出退货理由。
fn rejection_reasons(
    args: &Args,
    length: usize,
    unique_read_count: u64,
    supported_bases: usize,
    unique_read_density: f64,
    depth_cv: f64,
    max_depth_ratio: f64,
) -> Vec<&'static str> {
    let mut reasons = Vec::new();
    if unique_read_count == 0 {
        reasons.push("no_unique_read_support");
    }
    if supported_bases == 0 {
        reasons.push("no_positional_support");
    }
    if args.max_contig_length > 0 && length > args.max_contig_length {
        reasons.push("contig_too_long");
    }
    if length >= args.density_check_min_length && unique_read_density < args.min_read_density {
        reasons.push("low_unique_read_density");
    }
    if args.max_depth_cv > 0.0 && depth_cv > args.max_depth_cv {
        reasons.push("high_depth_cv");
    }
    if args.max_depth_ratio > 0.0 && max_depth_ratio > args.max_depth_ratio {
        reasons.push("repeat_depth_peak");
    }
    reasons
}

#[allow(clippy::too_many_arguments)]
pub fn assemble_seed(
    args: &Args,
    reads: &HashMap<Vec<u8>, u64>,
    slice_len: usize,
    graph: &HashMap<u128, KmerInfo>,
    seed: u128,
    k: usize,
    soft_boundary: usize,
    branch_support: Option<&BranchSupport>,
) -> (Vec<ContigRecord>, HashSet<u128>, i32) {
    // 从一个 seed 往两边接，组装成候选 contig；成不成得靠后面的证据说话。
    let reverse_seed = reverse_complement_kmer(seed, k);
    let empty_reverse_reuse_guard = HashSet::new();
    let (right_paths, right_kmers, right_positions, right_weight) = if args.assembly_mode
        == AssemblyMode::Uce
        && args.path_strategy == PathStrategy::Backbone
    {
        walk_backbone(
            graph,
            seed,
            k,
            args.iteration,
            args.backbone_lookahead,
            branch_support,
            &empty_reverse_reuse_guard,
            args.reverse_reuse_reference_scale,
        )
    } else {
        walk_search(
            graph,
            seed,
            k,
            args.iteration,
            &empty_reverse_reuse_guard,
            args.reverse_reuse_reference_scale,
        )
    };
    let (left_paths, left_kmers, left_positions, left_weight) = if args.assembly_mode
        == AssemblyMode::Uce
        && args.path_strategy == PathStrategy::Backbone
    {
        walk_backbone(
            graph,
            reverse_seed,
            k,
            args.iteration,
            args.backbone_lookahead,
            branch_support,
            &right_kmers,
            args.reverse_reuse_reference_scale,
        )
    } else {
        walk_search(
            graph,
            reverse_seed,
            k,
            args.iteration,
            &right_kmers,
            args.reverse_reuse_reference_scale,
        )
    };

    let mut positions: Vec<i64> = right_positions
        .into_iter()
        .chain(left_positions)
        .filter(|position| *position > 0 && *position < 1000)
        .map(i64::from)
        .collect();
    let contig_position = if positions.len() > 1 {
        median(&positions) as i32
    } else {
        -1
    };
    positions.clear();

    let mut right = process_sides(
        right_paths,
        right_weight,
        slice_len,
        reads,
        soft_boundary,
        args.assembly_mode,
    );
    let mut left = process_sides(
        left_paths,
        left_weight,
        slice_len,
        reads,
        soft_boundary,
        args.assembly_mode,
    );
    if right.is_empty() {
        right.push(SideContig::default());
    }
    if left.is_empty() {
        left.push(SideContig::default());
    }
    let candidate_limit = if args.assembly_mode == AssemblyMode::Uce
        && args.path_strategy == PathStrategy::Backbone
    {
        1
    } else if args.assembly_mode == AssemblyMode::Uce {
        args.side_candidates
    } else if args.assembly_mode == AssemblyMode::Its2 {
        16
    } else {
        3
    };

    let seed_sequence = decode_kmer(seed, k);
    let mut candidates = Vec::new();
    for left_side in left.iter().take(candidate_limit) {
        for right_side in right.iter().take(candidate_limit) {
            let mut sequence = reverse_complement(&left_side.sequence);
            sequence.extend_from_slice(&seed_sequence);
            sequence.extend_from_slice(&right_side.sequence);
            let weight = left_side.weight + right_side.weight;
            let mut support = calculate_read_support(&sequence, slice_len, reads);
            let mut length = sequence.len();
            let mut depth = if args.assembly_mode == AssemblyMode::Uce {
                depth_stats(&sequence, k, graph)
            } else {
                (0.0, 0.0, 0.0)
            };

            if args.min_coverage > 0.0 {
                let positional_count = if args.assembly_mode == AssemblyMode::Uce {
                    support.unique_read_count
                } else {
                    support.total_read_count
                };
                let positional_span = if args.assembly_mode == AssemblyMode::Uce {
                    support.supported_bases
                } else {
                    support.supported_extent
                };
                let coverage_depth = positional_count as f64 * slice_len as f64 / 0.9;
                if positional_span == 0
                    || coverage_depth / (positional_span as f64) < args.min_coverage
                {
                    continue;
                }
                if coverage_depth / (length as f64) < args.min_coverage {
                    sequence = sequence[support.left_coord..support.right_coord].to_vec();
                    length = sequence.len();
                    support = calculate_read_support(&sequence, slice_len, reads);
                    if args.assembly_mode == AssemblyMode::Uce {
                        depth = depth_stats(&sequence, k, graph);
                    }
                }
            }

            let read_density = if length > 0 {
                support.total_read_count as f64 / length as f64
            } else {
                0.0
            };
            let support_fraction = if length > 0 {
                support.supported_extent as f64 / length as f64
            } else {
                0.0
            };
            let unique_density = if length > 0 {
                support.unique_read_count as f64 / length as f64
            } else {
                0.0
            };
            let reasons = if args.assembly_mode == AssemblyMode::Uce {
                rejection_reasons(
                    args,
                    length,
                    support.unique_read_count,
                    support.supported_bases,
                    unique_density,
                    depth.1,
                    depth.2,
                )
            } else {
                Vec::new()
            };

            candidates.push(ContigRecord {
                sequence,
                position: contig_position,
                weight,
                read_count: support.total_read_count,
                supported_span: support.supported_extent,
                flank_balance: support.flank_balance,
                read_density,
                support_fraction,
                kmer_median_depth: depth.0,
                kmer_depth_cv: depth.1,
                kmer_max_depth_ratio: depth.2,
                unique_read_count: support.unique_read_count,
                multi_mapping_read_count: support.multi_mapping_read_count,
                supported_bases: support.supported_bases,
                support_breadth: support.breadth,
                max_support_gap: support.max_gap,
                accepted: reasons.is_empty(),
                rejection_reason: reasons.join(";"),
                ..ContigRecord::default()
            });
        }
    }

    let all_kmers = right_kmers.union(&left_kmers).copied().collect();
    (candidates, all_kmers, contig_position)
}

// 候选排序按模式办事；UCE 更看证据和稳当，不是单纯挑最长。
pub fn compare_contigs(left: &ContigRecord, right: &ContigRecord, mode: AssemblyMode) -> Ordering {
    if mode == AssemblyMode::Reference {
        return left
            .read_count
            .cmp(&right.read_count)
            .then_with(|| left.weight.cmp(&right.weight))
            .then_with(|| left.sequence.cmp(&right.sequence))
            .then_with(|| left.label.cmp(&right.label));
    }

    fn score(contig: &ContigRecord) -> [f64; 9] {
        let length = contig.sequence.len();
        let unique_density = if length > 0 {
            contig.unique_read_count as f64 / length as f64
        } else {
            0.0
        };
        let density_factor = (unique_density / 0.01).min(1.0);
        let continuity_factor = 1.0 / (1.0 + contig.kmer_depth_cv);
        let repeat_factor = (10.0 / contig.kmer_max_depth_ratio.max(1.0)).min(1.0);
        let effective =
            contig.supported_bases as f64 * density_factor * continuity_factor * repeat_factor;
        let gap_fraction = if length > 0 {
            contig.max_support_gap as f64 / length as f64
        } else {
            1.0
        };
        [
            effective,
            contig.support_breadth,
            unique_density,
            -gap_fraction,
            length as f64,
            contig.unique_read_count as f64,
            contig.read_count as f64,
            contig.flank_balance,
            contig.weight as f64,
        ]
    }

    let left_score = score(left);
    let right_score = score(right);
    for (left_value, right_value) in left_score.iter().zip(right_score.iter()) {
        let ordering = left_value.total_cmp(right_value);
        if ordering != Ordering::Equal {
            return ordering;
        }
    }
    left.sequence
        .cmp(&right.sequence)
        .then_with(|| left.label.cmp(&right.label))
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::seq::encode_kmer;
    fn info(depth: i64) -> KmerInfo {
        KmerInfo {
            depth,
            position: 100,
            is_reverse: false,
            reference_weight: 0,
            fragment_support: 0,
            pe_flank_rescue: false,
        }
    }

    #[test]
    fn backbone_commits_long_arm_and_discards_sibling() {
        let k = 4;
        let entries = [
            ("AAAA", 5),
            ("AAAC", 2),
            ("AACG", 2),
            ("ACGT", 2),
            ("CGTT", 2),
            ("GTTT", 2),
            ("AAAG", 20),
            ("AAGA", 20),
            ("AGAC", 20),
        ];
        let graph: HashMap<_, _> = entries
            .into_iter()
            .map(|(sequence, depth)| (encode_kmer(sequence.as_bytes()).unwrap(), info(depth)))
            .collect();
        let seed = encode_kmer(b"AAAA").unwrap();
        let (paths, visited, _, _) =
            walk_backbone(&graph, seed, k, 100, 24, None, &HashSet::new(), 1.0);
        let extension: Vec<u8> = paths[0].bases.iter().map(|base| bits_base(*base)).collect();
        assert_eq!(extension, b"CGTTT");
        assert!(visited.contains(&encode_kmer(b"AAAC").unwrap()));
        assert!(!visited.contains(&encode_kmer(b"AAAG").unwrap()));
    }

    #[test]
    fn backbone_cycle_stops() {
        let graph: HashMap<_, _> = ["AAA", "AAC", "ACA", "CAA"]
            .into_iter()
            .map(|sequence| (encode_kmer(sequence.as_bytes()).unwrap(), info(5)))
            .collect();
        let seed = encode_kmer(b"AAA").unwrap();
        let (paths, visited, _, _) =
            walk_backbone(&graph, seed, 3, 100, 24, None, &HashSet::new(), 1.0);
        assert_eq!(paths[0].bases.len(), 3);
        assert_eq!(visited.len(), 4);
    }

    #[test]
    fn pe_fragment_support_keeps_only_independently_supported_low_depth_kmers() {
        let mut reference = HashMap::new();
        reference.insert(
            encode_kmer(b"CCC").unwrap(),
            RefKmer {
                depth: 1,
                position: 100,
                is_reverse: false,
            },
        );
        let two_fragments = vec![
            b"CCCA".to_vec(),
            b"AAAC".to_vec(),
            b"CCCA".to_vec(),
            b"AAAG".to_vec(),
        ];
        let counts = count_assemble_chunk_parallel(&two_fragments, 3, 2, true, Some(&reference));
        let kmer = encode_kmer(b"AAA").unwrap();
        let graph = build_graph_from_counts(counts, &reference, 2, 2);
        assert_eq!(graph[&kmer].depth, 2);
        assert_eq!(graph[&kmer].fragment_support, 2);

        let graph = build_graph_from_counts(
            vec![KmerCount {
                kmer,
                depth: 2,
                fragment_support: 1,
            }],
            &reference,
            2,
            2,
        );
        assert!(!graph.contains_key(&kmer));
    }

    #[test]
    fn pe_branch_support_changes_only_branch_ranking() {
        let current = encode_kmer(b"AAA").unwrap();
        let standard = encode_kmer(b"AAC").unwrap();
        let pe_supported = encode_kmer(b"AAG").unwrap();
        let mut graph = HashMap::new();
        graph.insert(current, info(5));
        graph.insert(standard, info(5));
        let mut flank = info(2);
        flank.pe_flank_rescue = true;
        flank.fragment_support = 2;
        graph.insert(pe_supported, flank);
        let blocked = HashSet::from([current]);

        assert_eq!(
            outgoing(
                &graph,
                current,
                3,
                &blocked,
                &HashSet::new(),
                None,
                false,
                None,
                1.0
            )[0]
            .kmer,
            standard
        );
        let support = HashMap::from([((current, pe_supported), 2)]);
        assert_eq!(
            outgoing(
                &graph,
                current,
                3,
                &blocked,
                &HashSet::new(),
                None,
                false,
                Some(&support),
                1.0
            )[0]
            .kmer,
            pe_supported
        );
    }

    #[test]
    fn reverse_reuse_scales_reference_bonus_but_keeps_read_depth() {
        let current = encode_kmer(b"AAC").unwrap();
        let reverse_reused = encode_kmer(b"ACA").unwrap();
        let fresh = encode_kmer(b"ACC").unwrap();
        let mut graph = HashMap::new();
        graph.insert(current, info(5));
        let mut previous_path = info(2);
        previous_path.reference_weight = 80;
        graph.insert(reverse_reused, previous_path);
        graph.insert(fresh, info(30));
        let blocked = HashSet::from([current]);
        let other_arm = HashSet::from([reverse_complement_kmer(reverse_reused, 3)]);

        assert_eq!(
            outgoing(&graph, current, 3, &blocked, &other_arm, None, false, None, 1.0)[0].kmer,
            reverse_reused
        );
        assert_eq!(
            outgoing(&graph, current, 3, &blocked, &other_arm, None, false, None, 0.25)[0].kmer,
            fresh
        );
        graph.get_mut(&reverse_reused).unwrap().depth = 50;
        assert_eq!(
            outgoing(&graph, current, 3, &blocked, &other_arm, None, false, None, 0.25)[0].kmer,
            reverse_reused
        );
    }

    #[test]
    fn narrow_count_keys_match_wide_count_keys() {
        let mut state = 0x6a09_e667_f3bc_c909_u64;
        let mut next = |bound: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % bound as u64) as usize
        };
        let genome: Vec<u8> = (0..600).map(|_| b"ACGT"[next(4)]).collect();
        let sequences: Vec<Vec<u8>> = (0..2_001)
            .map(|_| {
                let start = next(genome.len() - 40);
                let mut read = genome[start..start + 20 + next(20)].to_vec();
                for base in read.iter_mut() {
                    match next(30) {
                        0 => *base = b'N',
                        1 => *base = b"ACGT"[next(4)],
                        _ => {}
                    }
                }
                read
            })
            .collect();
        for k in [3, 11, 21, 31, 32] {
            let reference =
                crate::io_utils::build_reference_kmers(&[("r".into(), genome[..200].to_vec())], k);
            for paired in [false, true] {
                for anchors in [None, Some(&reference)] {
                    let narrow = count_chunk::<u64>(&sequences, k, 1, paired, anchors);
                    let wide = count_chunk::<u128>(&sequences, k, 1, paired, anchors);
                    assert_eq!(narrow, wide, "k={k} paired={paired}");
                    // Dynamic block claiming must not depend on the worker count.
                    for threads in [2, 3, 8] {
                        let parallel =
                            count_assemble_chunk_parallel(&sequences, k, threads, paired, anchors);
                        assert_eq!(parallel, narrow, "k={k} paired={paired} threads={threads}");
                    }
                }
            }
        }
    }

    #[test]
    fn parallel_chunk_count_matches_legacy_count() {
        let sequences = vec![
            b"ATGCCATG".to_vec(),
            b"ATGCCATG".to_vec(),
            b"TTGCCATA".to_vec(),
        ];
        let reference = HashMap::new();
        let expected = build_assemble_dictionary(&sequences, 4, &reference);
        let mut accumulator = SortedKmerCounts::default();
        for chunk in sequences.chunks(2) {
            accumulator.push(count_assemble_chunk_parallel(chunk, 4, 2, false, None));
        }
        let observed = build_graph_from_counts(accumulator.into_counts(), &reference, 0, 0);
        assert_eq!(expected.len(), observed.len());
        for (kmer, value) in expected {
            assert_eq!(value.depth, observed[&kmer].depth);
        }
    }
}
