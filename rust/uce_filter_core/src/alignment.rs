use crate::index::{ExactSeed, UceIndex};
use crate::model::LocusId;

const MATCH_SCORE: i32 = 2;
const MISMATCH_SCORE: i32 = -3;
const GAP_OPEN: i32 = -5;
const GAP_EXTEND: i32 = -1;
const STOP: u8 = 3;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AlignmentEvidence {
    pub sequence: u32,
    pub strand: u8,
    pub score: i32,
    pub exact_seed_length: u16,
    pub query_start: u16,
    pub query_end: u16,
    pub reference_start: u32,
    pub reference_end: u32,
    pub matches: u16,
    pub mismatches: u16,
    pub gap_bases: u16,
}

impl AlignmentEvidence {
    pub fn aligned_columns(self) -> usize {
        self.matches as usize + self.mismatches as usize + self.gap_bases as usize
    }

    pub fn reference_overlap(self) -> usize {
        self.reference_end as usize - self.reference_start as usize
    }

    pub fn identity(self) -> f64 {
        let columns = self.aligned_columns();
        if columns == 0 {
            0.0
        } else {
            self.matches as f64 / columns as f64
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct TerminalEvidence {
    pub effective_window: u32,
    pub original_reference_start: u32,
    pub original_reference_end: u32,
    pub left_overhang: u16,
    pub right_overhang: u16,
    pub near_left_terminal: bool,
    pub near_right_terminal: bool,
    pub extends_left_terminal: bool,
    pub extends_right_terminal: bool,
}

pub fn terminal_evidence(
    alignment: AlignmentEvidence,
    query_length: usize,
    reference_length: usize,
    terminal_window: usize,
) -> TerminalEvidence {
    let query_left = alignment.query_start as usize;
    let query_right = query_length.saturating_sub(alignment.query_end as usize);
    let (reference_start, reference_end, left_overhang, right_overhang) = if alignment.strand == 1 {
        (
            alignment.reference_start as usize,
            alignment.reference_end as usize,
            query_left,
            query_right,
        )
    } else {
        (
            reference_length.saturating_sub(alignment.reference_end as usize),
            reference_length.saturating_sub(alignment.reference_start as usize),
            query_right,
            query_left,
        )
    };
    let effective_window = terminal_window.min((reference_length / 5).max(1));
    let near_left_terminal = reference_start <= effective_window;
    let near_right_terminal = reference_length.saturating_sub(reference_end) <= effective_window;
    TerminalEvidence {
        effective_window: effective_window.min(u32::MAX as usize) as u32,
        original_reference_start: reference_start.min(u32::MAX as usize) as u32,
        original_reference_end: reference_end.min(u32::MAX as usize) as u32,
        left_overhang: left_overhang.min(u16::MAX as usize) as u16,
        right_overhang: right_overhang.min(u16::MAX as usize) as u16,
        near_left_terminal,
        near_right_terminal,
        extends_left_terminal: near_left_terminal && left_overhang > 0,
        extends_right_terminal: near_right_terminal && right_overhang > 0,
    }
}

fn best_transition(candidates: &[(i32, u8)]) -> (i32, u8) {
    let mut best = (0, STOP);
    for &(score, state) in candidates {
        if score > best.0 || (score == best.0 && score > 0 && state < best.1) {
            best = (score, state);
        }
    }
    best
}

fn seeded_local_affine(
    query: &[u8],
    reference: &[u8],
    seed: ExactSeed,
    strand: u8,
    band: usize,
) -> Option<AlignmentEvidence> {
    if query.is_empty() || reference.is_empty() || seed.is_empty() {
        return None;
    }
    let expected_start = (seed.reference_start as usize).saturating_sub(seed.read_start as usize);
    let slice_start = expected_start.saturating_sub(band);
    let slice_end = (expected_start + query.len() + band).min(reference.len());
    if slice_start >= slice_end {
        return None;
    }
    let target = &reference[slice_start..slice_end];
    let rows = query.len() + 1;
    let cols = target.len() + 1;
    let cells = rows.checked_mul(cols)?;
    // Scores depend only on the previous row; traceback needs just two bits
    // per state. Keep the original traversal and tie-breaking order.
    let mut previous_scores = vec![[0_i32; 3]; cols];
    let mut current_scores = vec![[0_i32; 3]; cols];
    let mut trace = vec![0b11_11_11_u8; cells];
    let mut best = (0_i32, 0_usize, 0_usize, 0_u8);

    for i in 1..rows {
        for j in 1..cols {
            let idx = i * cols + j;

            let substitution = if query[i - 1] == target[j - 1] {
                MATCH_SCORE
            } else {
                MISMATCH_SCORE
            };
            let matched = best_transition(&[
                (previous_scores[j - 1][0] + substitution, 0),
                (previous_scores[j - 1][1] + substitution, 1),
                (previous_scores[j - 1][2] + substitution, 2),
            ]);

            let inserted = best_transition(&[
                (previous_scores[j][0] + GAP_OPEN, 0),
                (previous_scores[j][1] + GAP_EXTEND, 1),
                (previous_scores[j][2] + GAP_OPEN, 2),
            ]);

            let deleted = best_transition(&[
                (current_scores[j - 1][0] + GAP_OPEN, 0),
                (current_scores[j - 1][1] + GAP_OPEN, 1),
                (current_scores[j - 1][2] + GAP_EXTEND, 2),
            ]);

            current_scores[j] = [matched.0, inserted.0, deleted.0];
            trace[idx] = matched.1 | (inserted.1 << 2) | (deleted.1 << 4);
            for state in 0..3_u8 {
                let score = current_scores[j][state as usize];
                if score > best.0
                    || (score == best.0 && score > 0 && (i, j, state) < (best.1, best.2, best.3))
                {
                    best = (score, i, j, state);
                }
            }
        }
        std::mem::swap(&mut previous_scores, &mut current_scores);
    }
    if best.0 <= 0 {
        return None;
    }

    let (score, query_end, target_end, mut state) = best;
    let (mut i, mut j) = (query_end, target_end);
    let (mut matches, mut mismatches, mut gaps) = (0_usize, 0_usize, 0_usize);
    loop {
        let idx = i * cols + j;
        if state == STOP {
            break;
        }
        let previous = (trace[idx] >> (2 * state)) & 3;
        if previous == STOP {
            break;
        }
        match state {
            0 => {
                if query[i - 1] == target[j - 1] {
                    matches += 1;
                } else {
                    mismatches += 1;
                }
                i -= 1;
                j -= 1;
            }
            1 => {
                gaps += 1;
                i -= 1;
            }
            2 => {
                gaps += 1;
                j -= 1;
            }
            _ => break,
        }
        state = previous;
    }
    Some(AlignmentEvidence {
        sequence: seed.sequence,
        strand,
        score,
        exact_seed_length: seed.len().min(u16::MAX as usize) as u16,
        query_start: i.min(u16::MAX as usize) as u16,
        query_end: query_end.min(u16::MAX as usize) as u16,
        reference_start: (slice_start + j).min(u32::MAX as usize) as u32,
        reference_end: (slice_start + target_end).min(u32::MAX as usize) as u32,
        matches: matches.min(u16::MAX as usize) as u16,
        mismatches: mismatches.min(u16::MAX as usize) as u16,
        gap_bases: gaps.min(u16::MAX as usize) as u16,
    })
}

#[cfg(test)]
fn seeded_local_affine_full_matrix(
    query: &[u8],
    reference: &[u8],
    seed: ExactSeed,
    strand: u8,
    band: usize,
) -> Option<AlignmentEvidence> {
    if query.is_empty() || reference.is_empty() || seed.is_empty() {
        return None;
    }
    let expected_start = (seed.reference_start as usize).saturating_sub(seed.read_start as usize);
    let slice_start = expected_start.saturating_sub(band);
    let slice_end = (expected_start + query.len() + band).min(reference.len());
    if slice_start >= slice_end {
        return None;
    }
    let target = &reference[slice_start..slice_end];
    let rows = query.len() + 1;
    let cols = target.len() + 1;
    let cells = rows.checked_mul(cols)?;
    let mut scores = [vec![0_i32; cells], vec![0_i32; cells], vec![0_i32; cells]];
    let mut trace = [vec![STOP; cells], vec![STOP; cells], vec![STOP; cells]];
    let mut best = (0_i32, 0_usize, 0_usize, 0_u8);

    for i in 1..rows {
        for j in 1..cols {
            let idx = i * cols + j;
            let diagonal = (i - 1) * cols + j - 1;
            let above = (i - 1) * cols + j;
            let left = i * cols + j - 1;
            let substitution = if query[i - 1] == target[j - 1] {
                MATCH_SCORE
            } else {
                MISMATCH_SCORE
            };
            let matched = best_transition(&[
                (scores[0][diagonal] + substitution, 0),
                (scores[1][diagonal] + substitution, 1),
                (scores[2][diagonal] + substitution, 2),
            ]);
            scores[0][idx] = matched.0;
            trace[0][idx] = matched.1;

            let inserted = best_transition(&[
                (scores[0][above] + GAP_OPEN, 0),
                (scores[1][above] + GAP_EXTEND, 1),
                (scores[2][above] + GAP_OPEN, 2),
            ]);
            scores[1][idx] = inserted.0;
            trace[1][idx] = inserted.1;

            let deleted = best_transition(&[
                (scores[0][left] + GAP_OPEN, 0),
                (scores[1][left] + GAP_OPEN, 1),
                (scores[2][left] + GAP_EXTEND, 2),
            ]);
            scores[2][idx] = deleted.0;
            trace[2][idx] = deleted.1;

            for state in 0..3_u8 {
                let score = scores[state as usize][idx];
                if score > best.0
                    || (score == best.0 && score > 0 && (i, j, state) < (best.1, best.2, best.3))
                {
                    best = (score, i, j, state);
                }
            }
        }
    }
    if best.0 <= 0 {
        return None;
    }

    let (score, query_end, target_end, mut state) = best;
    let (mut i, mut j) = (query_end, target_end);
    let (mut matches, mut mismatches, mut gaps) = (0_usize, 0_usize, 0_usize);
    loop {
        let idx = i * cols + j;
        if state == STOP || scores[state as usize][idx] == 0 {
            break;
        }
        let previous = trace[state as usize][idx];
        match state {
            0 => {
                if query[i - 1] == target[j - 1] {
                    matches += 1;
                } else {
                    mismatches += 1;
                }
                i -= 1;
                j -= 1;
            }
            1 => {
                gaps += 1;
                i -= 1;
            }
            2 => {
                gaps += 1;
                j -= 1;
            }
            _ => break,
        }
        state = previous;
    }
    Some(AlignmentEvidence {
        sequence: seed.sequence,
        strand,
        score,
        exact_seed_length: seed.len().min(u16::MAX as usize) as u16,
        query_start: i.min(u16::MAX as usize) as u16,
        query_end: query_end.min(u16::MAX as usize) as u16,
        reference_start: (slice_start + j).min(u32::MAX as usize) as u32,
        reference_end: (slice_start + target_end).min(u32::MAX as usize) as u32,
        matches: matches.min(u16::MAX as usize) as u16,
        mismatches: mismatches.min(u16::MAX as usize) as u16,
        gap_bases: gaps.min(u16::MAX as usize) as u16,
    })
}

pub fn align_read(
    index: &UceIndex,
    query: &[u8],
    locus: LocusId,
    band: usize,
) -> Option<AlignmentEvidence> {
    let seed = index.best_exact(query, locus)?;
    let reference = index.references.get(seed.sequence as usize)?;
    seeded_local_affine(query, &reference.bases, seed, reference.strand, band)
}

/// Marker boundaries may carry substitutions that a local alignment clips.
/// Extend only into observed query bases; missing sample ends are never filled.
/// This is separate from read recruitment, whose local alignment is unchanged.
pub fn align_marker(
    index: &UceIndex,
    query: &[u8],
    locus: LocusId,
    band: usize,
) -> Option<AlignmentEvidence> {
    if query.len() > u16::MAX as usize {
        return None;
    }
    let mut alignment = align_read(index, query, locus, band)?;
    let reference = &index.references[alignment.sequence as usize].bases;
    let left = alignment.reference_start as usize;
    let right = reference.len() - alignment.reference_end as usize;
    let mut count = |a: &[u8], b: &[u8]| {
        for (&a, &b) in a.iter().zip(b) {
            if a == b {
                alignment.matches += 1;
                alignment.score += MATCH_SCORE;
            } else {
                alignment.mismatches += 1;
                alignment.score += MISMATCH_SCORE;
            }
        }
    };
    let start = alignment.query_start as usize;
    let end = alignment.query_end as usize;
    if left <= band && start >= left {
        count(&query[start - left..start], &reference[..left]);
        alignment.query_start -= left as u16;
        alignment.reference_start = 0;
    }
    if right <= band && query.len() - end >= right {
        count(
            &query[end..end + right],
            &reference[reference.len() - right..],
        );
        alignment.query_end += right as u16;
        alignment.reference_end = reference.len() as u32;
    }
    Some(alignment)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compact_trace_matches_full_matrix_including_ties_and_gaps() {
        let mut state = 146_u64;
        let mut random = |n| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(state >> 32) as usize % 4]
                })
                .collect()
        };
        for round in 0..600 {
            let reference = random(20 + round % 80);
            let start = round % (reference.len() - 12);
            let mut query = reference[start..].to_vec();
            match round % 5 {
                0 => {
                    query.extend_from_slice(&reference);
                }
                1 => {
                    query.insert(8, b'A');
                }
                2 => {
                    query.remove(8);
                }
                3 => {
                    query[8] = b'N';
                }
                _ => {
                    query = random(30);
                }
            }
            let seed = ExactSeed {
                sequence: 0,
                read_start: 0,
                read_end: 7,
                reference_start: start as u32,
                reference_end: start as u32 + 7,
            };
            let band = round % 40;
            let strand = 1 + (round % 2) as u8;
            assert_eq!(
                seeded_local_affine(&query, &reference, seed, strand, band),
                seeded_local_affine_full_matrix(&query, &reference, seed, strand, band),
                "round {round}"
            );
        }
    }

    #[test]
    fn marker_boundaries_preserve_exact_matches_and_repeats() {
        let root = std::env::temp_dir().join(format!("marker_exact_{}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        let sequence = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCG";
        let reference = root.join("ref.fasta");
        std::fs::write(
            &reference,
            format!(">ref\n{}\n", String::from_utf8_lossy(sequence)),
        )
        .unwrap();
        let index = UceIndex::build(&reference, 11).unwrap();
        for query in [
            sequence.to_vec(),
            sequence[6..24].to_vec(),
            [b"TTTT".as_slice(), sequence, b"CCCC"].concat(),
            [sequence.as_slice(), b"TTTT", sequence].concat(),
        ] {
            assert_eq!(
                align_marker(&index, &query, 0, 32),
                align_read(&index, &query, 0, 32)
            );
        }
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn local_alignment_reports_identity_and_terminal_overhang() {
        let reference = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCG";
        let query = b"TTTTACGTTGCAACGATTCGGTACC";
        let seed = ExactSeed {
            sequence: 0,
            read_start: 4,
            read_end: 25,
            reference_start: 0,
            reference_end: 21,
        };
        let alignment = seeded_local_affine(query, reference, seed, 1, 16).unwrap();
        assert_eq!(alignment.matches, 21);
        assert_eq!(alignment.mismatches, 0);
        assert_eq!(alignment.identity(), 1.0);
        let terminal = terminal_evidence(alignment, query.len(), reference.len(), 4);
        assert_eq!(terminal.left_overhang, 4);
        assert!(terminal.extends_left_terminal);
        assert!(!terminal.extends_right_terminal);
    }

    #[test]
    fn reverse_alignment_swaps_terminal_sides() {
        let alignment = AlignmentEvidence {
            sequence: 0,
            strand: 2,
            score: 20,
            exact_seed_length: 10,
            query_start: 3,
            query_end: 13,
            reference_start: 0,
            reference_end: 10,
            matches: 10,
            mismatches: 0,
            gap_bases: 0,
        };
        let terminal = terminal_evidence(alignment, 13, 100, 5);
        assert!(terminal.near_right_terminal);
        assert_eq!(terminal.right_overhang, 3);
        assert!(terminal.extends_right_terminal);
    }
}
