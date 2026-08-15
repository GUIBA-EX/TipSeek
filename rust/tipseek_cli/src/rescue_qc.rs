use std::collections::HashMap;

const KMER: usize = 15;
const ANCHOR: usize = 25;
const MAX_KMER_OCCURRENCES: usize = 16;
const INVERTED_REPEAT_KMER: usize = 21;
const INVERTED_REPEAT_MAX_KMER_OCCURRENCES: usize = 64;
pub const MINIMUM_NEW_UNSUPPORTED_GAP: usize = 40;

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

/// Returns the longest internal interval without a coherent read chain.
/// A chain must contain matching 15-mers on one read/candidate diagonal and
/// retain 25 aligned bases on both sides of every reported boundary.
pub fn maximum_unsupported_internal_gap(sequence: &str, reads: &[(String, String)]) -> usize {
    if sequence.len() <= 2 * (KMER + ANCHOR) || sequence.len() < KMER {
        return 0;
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
    let margin = KMER + ANCHOR;
    maximum_bracketed_zero_run(&support[margin..sequence.len() - margin])
}

pub fn introduces_unsupported_internal_gap(
    before: &str,
    after: &str,
    reads: &[(String, String)],
) -> bool {
    let before_gap = maximum_unsupported_internal_gap(before, reads);
    let after_gap = maximum_unsupported_internal_gap(after, reads);
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
