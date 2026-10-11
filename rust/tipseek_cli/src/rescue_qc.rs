use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

const KMER: usize = 15;
const ANCHOR: usize = 25;
const MAX_KMER_OCCURRENCES: usize = 16;
const INVERTED_REPEAT_KMER: usize = 21;
const INVERTED_REPEAT_MAX_KMER_OCCURRENCES: usize = 64;
pub const MINIMUM_NEW_UNSUPPORTED_GAP: usize = 40;

#[cfg(test)]
fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

fn encode_dna_kmer(sequence: &[u8]) -> Option<u64> {
    let mut encoded = 0_u64;
    for base in sequence {
        encoded = (encoded << 2)
            | match base.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
    }
    Some(encoded)
}

fn reverse_complement_packed(mut encoded: u64, k: usize) -> u64 {
    let mut reverse = 0_u64;
    for _ in 0..k {
        reverse = (reverse << 2) | (3 - (encoded & 3));
        encoded >>= 2;
    }
    reverse
}

/// Detects an exact, long, self reverse-complement match. The 21-mer chains
/// are grouped by anti-diagonal, so a forward run at increasing positions is
/// paired with a reverse-complement run at decreasing positions. Highly
/// repetitive 21-mers are ignored to keep the guard bounded.
pub fn has_long_inverted_repeat(sequence: &str, minimum_span: usize) -> bool {
    let k = INVERTED_REPEAT_KMER;
    if minimum_span == 0 || sequence.len() < k || sequence.len() < minimum_span {
        return false;
    }
    let bytes = sequence.as_bytes();
    let encoded = (0..=bytes.len() - k)
        .map(|start| encode_dna_kmer(&bytes[start..start + k]))
        .collect::<Vec<_>>();
    let mut positions: HashMap<u64, Vec<usize>> = HashMap::new();
    for (position, kmer) in encoded.iter().enumerate() {
        let Some(kmer) = kmer else {
            continue;
        };
        let bucket = positions.entry(*kmer).or_default();
        if bucket.len() <= INVERTED_REPEAT_MAX_KMER_OCCURRENCES {
            bucket.push(position);
        }
    }

    let mut pairs = Vec::new();
    for (left, kmer) in encoded.iter().enumerate() {
        let Some(kmer) = kmer else {
            continue;
        };
        let reverse = reverse_complement_packed(*kmer, k);
        let Some(matches) = positions.get(&reverse) else {
            continue;
        };
        if matches.len() > INVERTED_REPEAT_MAX_KMER_OCCURRENCES {
            continue;
        }
        pairs.extend(
            matches
                .iter()
                .copied()
                .filter(|right| *right > left)
                .map(|right| (left + right, left, right)),
        );
    }
    pairs.sort_unstable();

    let mut chain_start = None::<(usize, usize, usize)>;
    let mut previous = None::<(usize, usize, usize)>;
    for (diagonal, left, right) in pairs {
        let consecutive = previous.is_some_and(|(old_diagonal, old_left, old_right)| {
            diagonal == old_diagonal && left == old_left + 1 && right + 1 == old_right
        });
        if !consecutive {
            chain_start = Some((diagonal, left, right));
        }
        previous = Some((diagonal, left, right));
        let Some((_, start_left, start_right)) = chain_start else {
            continue;
        };
        let span = left - start_left + k;
        if span < minimum_span {
            continue;
        }
        let first_start = start_left;
        let first_end = left + k;
        let second_start = right;
        let second_end = start_right + k;
        let overlap = first_end
            .min(second_end)
            .saturating_sub(first_start.max(second_start));
        if overlap.saturating_mul(5) <= span {
            return true;
        }
    }
    false
}

fn maximum_bracketed_zero_run(values: &[usize]) -> usize {
    let mut longest = 0;
    let mut index = 0;
    while index < values.len() {
        if values[index] != 0 {
            index += 1;
            continue;
        }
        let start = index;
        while index < values.len() && values[index] == 0 {
            index += 1;
        }
        if start > 0 && index < values.len() {
            longest = longest.max(index - start);
        }
    }
    longest
}

/// Multiplicative hash for 2-bit packed k-mers; keys are local sequence data.
#[derive(Default)]
struct PackedKmerHasher(u64);

impl Hasher for PackedKmerHasher {
    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.write_u64(u64::from(byte));
        }
    }

    fn write_u32(&mut self, value: u32) {
        self.write_u64(u64::from(value));
    }

    fn write_u64(&mut self, value: u64) {
        self.0 = (self.0.rotate_left(5) ^ value).wrapping_mul(0x9e37_79b9_7f4a_7c15);
    }

    fn finish(&self) -> u64 {
        self.0 ^ (self.0 >> 29)
    }
}

type PackedKmerMap<V> = HashMap<u32, V, BuildHasherDefault<PackedKmerHasher>>;

fn base_code(base: u8) -> Option<u32> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Visits every all-ACGT 15-mer (case-insensitive) as `(offset, packed)`.
fn for_each_gap_kmer(sequence: &[u8], mut visit: impl FnMut(usize, u32)) {
    const MASK: u32 = (1 << (2 * KMER)) - 1;
    let (mut packed, mut valid) = (0_u32, 0_usize);
    for (end, &base) in sequence.iter().enumerate() {
        match base_code(base) {
            Some(code) => {
                packed = ((packed << 2) | code) & MASK;
                valid += 1;
            }
            None => valid = 0,
        }
        if valid >= KMER {
            visit(end + 1 - KMER, packed);
        }
    }
}

fn reverse_complement_gap_kmer(mut packed: u32) -> u32 {
    let mut reverse = 0_u32;
    for _ in 0..KMER {
        reverse = (reverse << 2) | (3 - (packed & 3));
        packed >>= 2;
    }
    reverse
}

/// Read-chain support for one candidate sequence, accumulated read by read.
struct GapSupport {
    length: usize,
    positions: PackedKmerMap<Vec<usize>>,
    difference: Vec<i64>,
}

impl GapSupport {
    /// `None` when the sequence is too short to have an internal gap.
    fn new(sequence: &str) -> Option<Self> {
        let length = sequence.len();
        if length <= 2 * (KMER + ANCHOR) {
            return None;
        }
        let mut positions = PackedKmerMap::<Vec<usize>>::default();
        for_each_gap_kmer(sequence.as_bytes(), |start, kmer| {
            positions.entry(kmer).or_default().push(start);
        });
        positions.retain(|_, starts| starts.len() <= MAX_KMER_OCCURRENCES);
        Some(Self {
            length,
            positions,
            difference: vec![0; length + 1],
        })
    }

    /// Adds one oriented read given as its packed 15-mers. Each diagonal
    /// contributes once with its min/max hit, so hits can be merged in any
    /// order: consecutive hits on one diagonal are folded before sorting.
    fn add_oriented(&mut self, kmers: &[(usize, u32)], hits: &mut Vec<(i64, usize, usize)>) {
        hits.clear();
        for &(offset, kmer) in kmers {
            for &start in self.positions.get(&kmer).into_iter().flatten() {
                let diagonal = start as i64 - offset as i64;
                match hits.last_mut() {
                    Some(last) if last.0 == diagonal => {
                        last.1 = last.1.min(start);
                        last.2 = last.2.max(start);
                    }
                    _ => hits.push((diagonal, start, start)),
                }
            }
        }
        hits.sort_unstable();
        for group in hits.chunk_by(|left, right| left.0 == right.0) {
            let minimum = group.iter().map(|hit| hit.1).min().expect("nonempty group");
            let maximum = group.iter().map(|hit| hit.2).max().expect("nonempty group");
            let start = minimum + KMER + ANCHOR;
            let end = maximum.saturating_sub(ANCHOR).saturating_add(1);
            if start < end && end <= self.length {
                self.difference[start] += 1;
                self.difference[end] -= 1;
            }
        }
    }

    fn support(self) -> Vec<usize> {
        let mut support = Vec::with_capacity(self.length);
        let mut current = 0_i64;
        for delta in self.difference.into_iter().take(self.length) {
            current += delta;
            support.push(current.max(0) as usize);
        }
        support
    }

    fn maximum_gap(self) -> usize {
        let length = self.length;
        let margin = KMER + ANCHOR;
        maximum_bracketed_zero_run(&self.support()[margin..length - margin])
    }
}

/// Feeds both orientations of every read to each candidate. Read k-mers are
/// packed once and shared by all candidates. The reverse-complement window of
/// forward window `o` starts at `len - KMER - o` and is valid exactly when the
/// forward window is, so both strands come from one pass over the read.
fn accumulate_read_support(candidates: &mut [&mut GapSupport], reads: &[(String, String)]) {
    let mut forward = Vec::new();
    let mut reverse = Vec::new();
    let mut hits = Vec::new();
    for (_, read) in reads {
        let bytes = read.as_bytes();
        forward.clear();
        reverse.clear();
        for_each_gap_kmer(bytes, |offset, packed| {
            forward.push((offset, packed));
            reverse.push((
                bytes.len() - KMER - offset,
                reverse_complement_gap_kmer(packed),
            ));
        });
        for candidate in candidates.iter_mut() {
            candidate.add_oriented(&forward, &mut hits);
            candidate.add_oriented(&reverse, &mut hits);
        }
    }
}

/// Returns the longest internal interval without a coherent read chain.
/// A chain must contain matching 15-mers on one read/candidate diagonal and
/// retain 25 aligned bases on both sides of every reported boundary.
pub fn maximum_unsupported_internal_gap(sequence: &str, reads: &[(String, String)]) -> usize {
    let Some(mut support) = GapSupport::new(sequence) else {
        return 0;
    };
    accumulate_read_support(&mut [&mut support], reads);
    support.maximum_gap()
}

pub fn introduces_unsupported_internal_gap(
    before: &str,
    after: &str,
    reads: &[(String, String)],
) -> bool {
    let mut before_support = GapSupport::new(before);
    let mut after_support = GapSupport::new(after);
    let mut candidates = before_support
        .iter_mut()
        .chain(after_support.iter_mut())
        .collect::<Vec<_>>();
    accumulate_read_support(&mut candidates, reads);
    let before_gap = before_support.map_or(0, GapSupport::maximum_gap);
    let after_gap = after_support.map_or(0, GapSupport::maximum_gap);
    after_gap >= MINIMUM_NEW_UNSUPPORTED_GAP
        && after_gap.saturating_sub(before_gap) >= MINIMUM_NEW_UNSUPPORTED_GAP
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna(length: usize) -> String {
        let mut state = 0x9e37_79b9_u64;
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize] as char
            })
            .collect()
    }

    fn tiled_reads(sequence: &str) -> Vec<(String, String)> {
        (0..=sequence.len() - 120)
            .step_by(30)
            .enumerate()
            .map(|(index, start)| {
                (
                    format!("read{index}"),
                    sequence[start..start + 120].to_owned(),
                )
            })
            .collect()
    }

    /// The pre-optimization implementation, kept verbatim as an oracle.
    fn legacy_maximum_unsupported_internal_gap(
        sequence: &str,
        reads: &[(String, String)],
    ) -> usize {
        legacy_support(sequence, reads).map_or(0, |support| {
            let margin = KMER + ANCHOR;
            maximum_bracketed_zero_run(&support[margin..sequence.len() - margin])
        })
    }

    fn legacy_support(sequence: &str, reads: &[(String, String)]) -> Option<Vec<usize>> {
        if sequence.len() <= 2 * (KMER + ANCHOR) || sequence.len() < KMER {
            return None;
        }
        let sequence = sequence.to_ascii_uppercase();
        let mut positions = HashMap::<String, Vec<usize>>::new();
        for start in 0..=sequence.len() - KMER {
            let word = &sequence[start..start + KMER];
            if word
                .bytes()
                .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
            {
                positions.entry(word.to_owned()).or_default().push(start);
            }
        }
        positions.retain(|_, starts| starts.len() <= MAX_KMER_OCCURRENCES);

        let mut difference = vec![0_i64; sequence.len() + 1];
        for (_, read) in reads {
            for oriented in [read.to_ascii_uppercase(), reverse_complement(read)] {
                if oriented.len() < KMER {
                    continue;
                }
                let mut diagonals = HashMap::<i64, (usize, usize)>::new();
                for offset in 0..=oriented.len() - KMER {
                    let word = &oriented[offset..offset + KMER];
                    for start in positions.get(word).into_iter().flatten().copied() {
                        let diagonal = start as i64 - offset as i64;
                        diagonals
                            .entry(diagonal)
                            .and_modify(|range| {
                                range.0 = range.0.min(start);
                                range.1 = range.1.max(start);
                            })
                            .or_insert((start, start));
                    }
                }
                for (minimum, maximum) in diagonals.into_values() {
                    let start = minimum + KMER + ANCHOR;
                    let end = maximum.saturating_sub(ANCHOR).saturating_add(1);
                    if start < end && end <= sequence.len() {
                        difference[start] += 1;
                        difference[end] -= 1;
                    }
                }
            }
        }

        let mut support = Vec::with_capacity(sequence.len());
        let mut current = 0_i64;
        for delta in difference.into_iter().take(sequence.len()) {
            current += delta;
            support.push(current.max(0) as usize);
        }
        Some(support)
    }

    #[test]
    fn packed_gap_support_matches_legacy_string_implementation() {
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut next = |bound: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % bound as u64) as usize
        };
        for _ in 0..300 {
            // Short repeat units make k-mers recur, including past the
            // occurrence cap; lowercase and N exercise the validity rules.
            let unit: Vec<u8> = (0..1 + next(20)).map(|_| b"ACGT"[next(4)]).collect();
            let source = |length: usize, next: &mut dyn FnMut(usize) -> usize| -> String {
                (0..length)
                    .map(|i| match next(40) {
                        0 => 'N',
                        1 => (b"acgt"[next(4)]) as char,
                        2..=12 => unit[i % unit.len()] as char,
                        _ => b"ACGT"[next(4)] as char,
                    })
                    .collect()
            };
            let before = source(next(400), &mut next);
            let mut after = before.clone();
            let cut = next(after.len() + 1);
            after.insert_str(cut, &source(next(120), &mut next));
            let reads: Vec<(String, String)> = (0..next(40))
                .map(|index| {
                    let template = if next(2) == 0 { &after } else { &before };
                    let start = next(template.len() + 1);
                    let mut read: String = template[start..].chars().take(5 + next(150)).collect();
                    if next(3) == 0 {
                        read = reverse_complement(&read);
                    }
                    if next(4) == 0 {
                        read.push_str(&source(next(30), &mut next));
                    }
                    (format!("r{index}"), read)
                })
                .collect();
            for sequence in [&before, &after] {
                let packed = GapSupport::new(sequence).map(|mut support| {
                    accumulate_read_support(&mut [&mut support], &reads);
                    support.support()
                });
                assert_eq!(packed, legacy_support(sequence, &reads));
                assert_eq!(
                    maximum_unsupported_internal_gap(sequence, &reads),
                    legacy_maximum_unsupported_internal_gap(sequence, &reads)
                );
            }
            let legacy_before = legacy_maximum_unsupported_internal_gap(&before, &reads);
            let legacy_after = legacy_maximum_unsupported_internal_gap(&after, &reads);
            assert_eq!(
                introduces_unsupported_internal_gap(&before, &after, &reads),
                legacy_after >= MINIMUM_NEW_UNSUPPORTED_GAP
                    && legacy_after.saturating_sub(legacy_before) >= MINIMUM_NEW_UNSUPPORTED_GAP
            );
        }
    }

    #[test]
    fn long_inverted_repeat_detector_finds_exact_nonoverlapping_arms() {
        let arm = dna(180);
        let spacer = reverse_complement(&dna(60));
        let sequence = format!("{arm}{spacer}{}", reverse_complement(&arm));
        assert!(has_long_inverted_repeat(&sequence, 150));
        assert!(!has_long_inverted_repeat(&sequence, 200));
        assert!(!has_long_inverted_repeat(&dna(700), 150));
        assert!(!has_long_inverted_repeat("NNNNACGTNNNN", 10));
    }

    #[test]
    fn detects_a_new_internal_sequence_without_a_spanning_read_chain() {
        let before = dna(500);
        let reads = tiled_reads(&before);
        let after = format!("{}{}{}", &before[..250], dna(80), &before[250..]);
        assert_eq!(maximum_unsupported_internal_gap(&before, &reads), 0);
        assert!(maximum_unsupported_internal_gap(&after, &reads) >= 40);
        assert!(introduces_unsupported_internal_gap(&before, &after, &reads));
    }

    #[test]
    fn does_not_reject_a_contig_with_coherent_tiled_read_support() {
        let before = dna(500);
        let after = format!("{}{}", before, dna(80));
        let reads = tiled_reads(&after);
        assert!(!introduces_unsupported_internal_gap(
            &after[..500],
            &after,
            &reads
        ));
    }
}
