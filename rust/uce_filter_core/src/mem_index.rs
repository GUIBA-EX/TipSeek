//! Per-locus FM indexes for exact-match evidence.
//!
//! Recruitment remains a global rolling-k-mer lookup. Only recruited loci are
//! queried here. Each locus stores its forward and reverse-complement reference
//! segments in one FM index, separated by sentinels. Runtime state is limited to
//! symbol rank bitvectors and sparse suffix-array samples: no positional anchor
//! hash table and no full suffix array or RMQ remain resident.
//!
//! For run-k values up to 32 (the default is 19) the FM index is replaced by
//! [`RunKAnchors`]: an exact run-k table plus the forward reference text.
//! Every run-k window then costs one hash lookup, and a MEM is extended by
//! comparing eight bases per step against each occurrence of its run-k
//! prefix, instead of two rank queries per base. Both paths return identical
//! occurrence counts, strand masks and seeds.

use crate::index::{code, OrientedReference};

const ALPHABET_SIZE: usize = 6;
const SA_SAMPLE_RATE: usize = 16;
const SENTINEL: u8 = 0;
const SEPARATOR: u8 = 5;
/// Largest run-k served by [`RunKAnchors`]; longer run-k uses the FM index.
const MAX_ANCHOR_K: usize = 32;
/// Empty run-k table slot.
const EMPTY_SLOT: u32 = u32::MAX;
/// Slot flag for a run-k that occurs more than once; the value bits then
/// index `RunKAnchors::repeat_starts` instead of holding a text position.
const REPEAT_SLOT: u32 = 1 << 31;
/// Read lengths packed on the stack; longer reads use a heap buffer.
const STACK_READ_WORDS: usize = 24;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct MemSeed {
    pub sequence: u32,
    pub read_start: usize,
    pub read_end: usize,
    pub reference_start: usize,
    pub reference_end: usize,
}

impl MemSeed {
    pub fn len(self) -> usize {
        self.read_end - self.read_start
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct MemQueryProfile {
    pub index_queries: u64,
    pub run_windows: u64,
    pub matching_windows: u64,
    pub mem_starts: u64,
    pub mem_bases: u64,
}

#[derive(Clone, Copy, Debug)]
struct ReferenceSpan {
    text_start: u32,
    text_end: u32,
    sequence: u32,
    reference_start: u32,
    reference_end: u32,
    strand: u8,
}

#[derive(Debug, Default)]
pub(crate) struct LocusMemIndex {
    length: usize,
    spans: Vec<ReferenceSpan>,
    cumulative: [u32; ALPHABET_SIZE],
    symbol_bits: [Vec<u64>; ALPHABET_SIZE],
    symbol_prefix: [Vec<u32>; ALPHABET_SIZE],
    sample_bits: Vec<u64>,
    sample_prefix: Vec<u32>,
    sample_values: Vec<u32>,
    /// Strand of the reference span holding each BWT row's suffix (0 at
    /// separators). A run-k match spans only bases, so it never crosses a
    /// span boundary and its strand is that of the span it starts in.
    row_strands: Vec<u8>,
    /// Exact run-k table that replaces the FM index when run-k <= 32.
    anchors: Option<RunKAnchors>,
}

/// Exact run-k occurrences over the forward reference text of one locus.
///
/// Memory is kept near 6–12 bytes per reference base: the text is packed at
/// two bits per base, and the open-addressing table stores one `u32` per
/// slot without the run-k itself. A slot holds the text position of a run-k
/// that occurs once (the common case) or flags a group of repeated
/// occurrences. Lookups confirm a candidate slot by reading its run-k back
/// from the packed text, which is small enough to stay in cache.
#[derive(Debug, Default)]
struct RunKAnchors {
    run_k: usize,
    /// Forward reference bases of all spans, two bits each, first base in
    /// the low bits of a word; one trailing zero word keeps reads in range.
    packed: Vec<u64>,
    /// Start and end of each span in `packed`, parallel to
    /// `LocusMemIndex::spans`.
    span_starts: Vec<u32>,
    span_ends: Vec<u32>,
    /// `EMPTY_SLOT`, or `[REPEAT_SLOT flag | fingerprint | value]` where the
    /// value is a text position or, with the flag, a repeat group. The value
    /// uses the fewest bits that fit this locus; the rest of the low 31 bits
    /// hold a run-k hash fingerprint that rejects most probes without
    /// reading the packed text.
    slots: Vec<u32>,
    slot_bits: u32,
    value_bits: u32,
    /// Group `g` occurs at `repeat_positions[repeat_starts[g]..repeat_starts[g + 1]]`.
    repeat_starts: Vec<u32>,
    repeat_positions: Vec<u32>,
}

/// Text positions of one run-k.
enum Occurrences<'a> {
    One(u32),
    Many(&'a [u32]),
}

impl RunKAnchors {
    fn build(run_k: usize, span_texts: &[Vec<u8>]) -> Result<Self, String> {
        debug_assert!((1..=MAX_ANCHOR_K).contains(&run_k));
        let total = span_texts.iter().map(Vec::len).sum::<usize>();
        if total >= (REPEAT_SLOT - 1) as usize {
            return Err("locus anchor text length exceeds 2^31 bases".to_string());
        }
        let mut anchors = Self {
            run_k,
            packed: vec![0_u64; total.div_ceil(32) + 1],
            span_starts: Vec::with_capacity(span_texts.len()),
            span_ends: Vec::with_capacity(span_texts.len()),
            ..Self::default()
        };
        let mut keyed = Vec::<(u64, u32)>::with_capacity(total);
        let mut position = 0_usize;
        for span in span_texts {
            anchors.span_starts.push(position as u32);
            let mut kmer = 0_u64;
            for (offset, &value) in span.iter().enumerate() {
                anchors.packed[position / 32] |= u64::from(value) << (2 * (position % 32));
                kmer = push_base(kmer, value, run_k);
                if offset + 1 >= run_k {
                    keyed.push((kmer, (position + 1 - run_k) as u32));
                }
                position += 1;
            }
            anchors.span_ends.push(position as u32);
        }
        keyed.sort_unstable();

        let groups = keyed
            .iter()
            .enumerate()
            .filter(|&(index, entry)| index == 0 || keyed[index - 1].0 != entry.0)
            .count();
        // Load factor at most 2/3 keeps probe chains short.
        let capacity = (groups + groups / 2 + 1).next_power_of_two().max(2);
        anchors.slot_bits = capacity.trailing_zeros();
        // One spare value keeps every occupied slot distinct from EMPTY_SLOT.
        anchors.value_bits = usize::BITS - (total.max(groups) + 1).leading_zeros();
        anchors.slots = vec![EMPTY_SLOT; capacity];
        anchors.repeat_starts.push(0);
        let mut first = 0_usize;
        while first < keyed.len() {
            let kmer = keyed[first].0;
            let mut end = first + 1;
            while end < keyed.len() && keyed[end].0 == kmer {
                end += 1;
            }
            let value = if end - first == 1 {
                keyed[first].1
            } else {
                anchors
                    .repeat_positions
                    .extend(keyed[first..end].iter().map(|&(_, position)| position));
                anchors
                    .repeat_starts
                    .push(anchors.repeat_positions.len() as u32);
                REPEAT_SLOT | (anchors.repeat_starts.len() - 2) as u32
            };
            let hash = slot_hash(kmer);
            let value = value | anchors.fingerprint(hash);
            let mut slot = anchors.home_slot(hash);
            while anchors.slots[slot] != EMPTY_SLOT {
                slot = (slot + 1) & (capacity - 1);
            }
            anchors.slots[slot] = value;
            first = end;
        }
        Ok(anchors)
    }

    /// Table index from the high hash bits.
    #[inline(always)]
    fn home_slot(&self, hash: u64) -> usize {
        (hash >> (64 - self.slot_bits)) as usize
    }

    /// Low hash bits placed above the value bits of a slot.
    #[inline(always)]
    fn fingerprint(&self, hash: u64) -> u32 {
        (hash as u32) & !REPEAT_SLOT & !self.value_mask()
    }

    #[inline(always)]
    fn value_mask(&self) -> u32 {
        (1_u32 << self.value_bits) - 1
    }

    #[inline(always)]
    fn find(&self, kmer: u64) -> Option<Occurrences<'_>> {
        let mask = run_k_mask(self.run_k);
        let hash = slot_hash(kmer);
        let fingerprint = self.fingerprint(hash);
        let fingerprint_mask = !REPEAT_SLOT & !self.value_mask();
        let mut slot = self.home_slot(hash);
        loop {
            let value = self.slots[slot];
            if value == EMPTY_SLOT {
                return None;
            }
            if value & fingerprint_mask != fingerprint {
                slot = (slot + 1) & (self.slots.len() - 1);
                continue;
            }
            let occurrences = if value & REPEAT_SLOT == 0 {
                Occurrences::One(value & self.value_mask())
            } else {
                let group = (value & self.value_mask()) as usize;
                Occurrences::Many(
                    &self.repeat_positions[self.repeat_starts[group] as usize
                        ..self.repeat_starts[group + 1] as usize],
                )
            };
            let first = match occurrences {
                Occurrences::One(position) => position,
                Occurrences::Many(positions) => positions[0],
            };
            if packed_word(&self.packed, first as usize) & mask == kmer {
                return Some(occurrences);
            }
            slot = (slot + 1) & (self.slots.len() - 1);
        }
    }

    /// Index into `LocusMemIndex::spans` of the span holding `position`.
    #[inline(always)]
    fn span_index(&self, position: u32) -> usize {
        self.span_starts.partition_point(|&start| start <= position) - 1
    }
}

#[inline(always)]
fn slot_hash(kmer: u64) -> u64 {
    (kmer ^ (kmer >> 31)).wrapping_mul(0x9e37_79b9_7f4a_7c15)
}

#[inline(always)]
fn run_k_mask(run_k: usize) -> u64 {
    if run_k >= MAX_ANCHOR_K {
        u64::MAX
    } else {
        (1_u64 << (2 * run_k)) - 1
    }
}

/// Appends a base to a rolling run-k whose first base is in the low bits,
/// matching the layout of `packed_word`.
#[inline(always)]
fn push_base(kmer: u64, value: u8, run_k: usize) -> u64 {
    (kmer >> 2) | (u64::from(value) << (2 * (run_k - 1)))
}

/// Thirty-two bases starting at `position`, first base in the low bits.
/// `words` must end with a zero word past the last base.
#[inline(always)]
fn packed_word(words: &[u64], position: usize) -> u64 {
    let word = position / 32;
    let shift = 2 * (position % 32);
    if shift == 0 {
        words[word]
    } else {
        (words[word] >> shift) | (words[word + 1] << (64 - shift))
    }
}

/// Packs base codes into `words` (with one trailing zero word).
fn pack_codes(codes: &[u8], words: &mut [u64]) {
    words.fill(0);
    for (position, &value) in codes.iter().enumerate() {
        words[position / 32] |= u64::from(value) << (2 * (position % 32));
    }
}

/// Number of equal bases, at most `limit`, comparing up to 32 per step.
#[inline(always)]
fn packed_common_prefix(
    read: &[u64],
    read_start: usize,
    text: &[u64],
    text_start: usize,
    limit: usize,
) -> usize {
    let mut matched = 0_usize;
    while matched < limit {
        let step = (limit - matched).min(32);
        let mut difference =
            packed_word(read, read_start + matched) ^ packed_word(text, text_start + matched);
        if step < 32 {
            difference &= (1_u64 << (2 * step)) - 1;
        }
        if difference != 0 {
            return matched + (difference.trailing_zeros() / 2) as usize;
        }
        matched += step;
    }
    limit
}

/// Writes the 0..=3 base codes of a validated DNA read into `codes`.
pub(crate) fn read_codes_into(read: &[u8], codes: &mut Vec<u8>) {
    codes.clear();
    codes.extend(read.iter().map(|&base| code(base).unwrap_or(0)));
}

impl LocusMemIndex {
    fn build(&mut self, text: &[u8]) -> Result<(), String> {
        if text.is_empty() {
            return Ok(());
        }
        let suffixes = build_suffix_array(text)?;
        self.length = text.len();
        let words = self.length.div_ceil(64);
        self.symbol_bits = std::array::from_fn(|_| vec![0_u64; words]);
        self.sample_bits = vec![0_u64; words];
        self.sample_values.clear();

        let mut frequencies = [0_u32; ALPHABET_SIZE];
        for &symbol in text {
            frequencies[symbol as usize] += 1;
        }
        let mut preceding = 0_u32;
        for (symbol, frequency) in frequencies.into_iter().enumerate() {
            self.cumulative[symbol] = preceding;
            preceding += frequency;
        }

        for (row, &position) in suffixes.iter().enumerate() {
            let symbol = if position == 0 {
                SENTINEL
            } else {
                text[position - 1]
            };
            self.symbol_bits[symbol as usize][row / 64] |= 1_u64 << (row % 64);
            if position % SA_SAMPLE_RATE == 0 {
                self.sample_bits[row / 64] |= 1_u64 << (row % 64);
                self.sample_values.push(position as u32);
            }
        }
        self.row_strands = suffixes
            .iter()
            .map(|&position| self.span_at(position).map_or(0, |span| span.strand))
            .collect();
        self.symbol_prefix = std::array::from_fn(|symbol| {
            let mut prefix = Vec::with_capacity(words + 1);
            prefix.push(0_u32);
            for &word in &self.symbol_bits[symbol] {
                prefix.push(prefix.last().copied().unwrap() + word.count_ones());
            }
            prefix
        });
        self.sample_prefix = Vec::with_capacity(words + 1);
        self.sample_prefix.push(0_u32);
        for &word in &self.sample_bits {
            self.sample_prefix
                .push(self.sample_prefix.last().copied().unwrap() + word.count_ones());
        }
        Ok(())
    }

    #[inline(always)]
    fn rank(&self, symbol: u8, end: usize) -> usize {
        let word = end / 64;
        let offset = end % 64;
        let mut count = self.symbol_prefix[symbol as usize][word] as usize;
        if offset > 0 {
            let mask = (1_u64 << offset) - 1;
            count += (self.symbol_bits[symbol as usize][word] & mask).count_ones() as usize;
        }
        count
    }

    #[inline(always)]
    fn extend(&self, interval: (usize, usize), symbol: u8) -> (usize, usize) {
        let base = self.cumulative[symbol as usize] as usize;
        (
            base + self.rank(symbol, interval.0),
            base + self.rank(symbol, interval.1),
        )
    }

    fn symbol_at(&self, row: usize) -> u8 {
        let mask = 1_u64 << (row % 64);
        for symbol in 0..ALPHABET_SIZE {
            if self.symbol_bits[symbol][row / 64] & mask != 0 {
                return symbol as u8;
            }
        }
        unreachable!("each BWT row has exactly one symbol")
    }

    fn sampled_sa(&self, row: usize) -> Option<usize> {
        let word = row / 64;
        let offset = row % 64;
        let bit = 1_u64 << offset;
        if self.sample_bits[word] & bit == 0 {
            return None;
        }
        let preceding_mask = bit - 1;
        let sample = self.sample_prefix[word] as usize
            + (self.sample_bits[word] & preceding_mask).count_ones() as usize;
        Some(self.sample_values[sample] as usize)
    }

    fn locate(&self, mut row: usize) -> usize {
        let mut steps = 0_usize;
        loop {
            if let Some(position) = self.sampled_sa(row) {
                return (position + steps) % self.length;
            }
            let symbol = self.symbol_at(row);
            row = self.cumulative[symbol as usize] as usize + self.rank(symbol, row);
            steps += 1;
            debug_assert!(steps < SA_SAMPLE_RATE);
        }
    }

    fn span_at(&self, position: usize) -> Option<ReferenceSpan> {
        let index = self
            .spans
            .partition_point(|span| span.text_end as usize <= position);
        self.spans.get(index).copied().filter(|span| {
            span.text_start as usize <= position && position < span.text_end as usize
        })
    }

    #[cfg(test)]
    fn interval(&self, pattern: &[u8]) -> (usize, usize) {
        let mut interval = (0_usize, self.length);
        for &base in pattern {
            interval = self.extend(interval, base_symbol(base));
            if interval.0 == interval.1 {
                break;
            }
        }
        interval
    }

    fn occurrence(&self, row: usize, length: usize) -> Option<(u32, usize, usize, u8)> {
        let text_start = self.locate(row);
        let span = self.span_at(text_start)?;
        if text_start + length > span.text_end as usize {
            return None;
        }
        let reference_end = span.reference_end as usize - (text_start - span.text_start as usize);
        let reference_start = reference_end.checked_sub(length)?;
        debug_assert!(reference_start >= span.reference_start as usize);
        Some((span.sequence, reference_start, reference_end, span.strand))
    }

    fn best_occurrence(
        &self,
        interval: (usize, usize),
        length: usize,
    ) -> Option<(u32, usize, usize)> {
        (interval.0..interval.1)
            .filter_map(|row| self.occurrence(row, length))
            .map(|(sequence, start, end, _)| (sequence, start, end))
            .min_by_key(|&(sequence, start, _)| (sequence, start))
    }

    /// Run-k occurrence statistics and the longest exact match of a read.
    ///
    /// `read` holds base codes 0..=3 (see [`read_codes_into`]). For every
    /// run-k window it records the occurrence count, and the strand when the
    /// window occurs once. It returns the longest match of at least run-k
    /// bases, ties broken by (sequence, reference start, read start).
    pub fn collect(
        &self,
        read: &[u8],
        run_k: usize,
        occurrence_counts: &mut [u32],
        strand_masks: &mut [u8],
        profile: &mut MemQueryProfile,
    ) -> Option<MemSeed> {
        if self.length == 0 {
            return None;
        }
        profile.index_queries += 1;
        match self.anchors.as_ref() {
            Some(anchors) if anchors.run_k == run_k => {
                self.collect_anchored(anchors, read, occurrence_counts, strand_masks, profile)
            }
            _ => self.collect_fm(read, run_k, occurrence_counts, strand_masks, profile),
        }
    }

    fn collect_anchored(
        &self,
        anchors: &RunKAnchors,
        read: &[u8],
        occurrence_counts: &mut [u32],
        strand_masks: &mut [u8],
        profile: &mut MemQueryProfile,
    ) -> Option<MemSeed> {
        let run_k = anchors.run_k;
        let words = read.len().div_ceil(32) + 1;
        let mut stack = [0_u64; STACK_READ_WORDS];
        let mut heap = Vec::new();
        let packed_read = if words <= STACK_READ_WORDS {
            &mut stack[..words]
        } else {
            heap.resize(words, 0);
            &mut heap[..]
        };
        pack_codes(read, packed_read);
        let packed_read = &*packed_read;

        let mut kmer = read[..run_k - 1]
            .iter()
            .fold(0_u64, |kmer, &value| push_base(kmer, value, run_k));
        let mut best = None::<MemSeed>;
        for read_start in 0..occurrence_counts.len() {
            profile.run_windows += 1;
            kmer = push_base(kmer, read[read_start + run_k - 1], run_k);
            let Some(found) = anchors.find(kmer) else {
                continue;
            };
            profile.matching_windows += 1;
            let single;
            let occurrences = match found {
                Occurrences::One(position) => {
                    single = position;
                    std::slice::from_ref(&single)
                }
                Occurrences::Many(positions) => positions,
            };
            occurrence_counts[read_start] = occurrences.len() as u32;
            if occurrences.len() == 1 {
                strand_masks[read_start] = self.spans[anchors.span_index(occurrences[0])].strand;
            }
            // A match starting here can be at most this long. Once it is
            // shorter than the current best it can neither win nor tie.
            if best.is_some_and(|current| read.len() - read_start < current.len()) {
                continue;
            }
            let read_tail = read_start + run_k;
            let mut chosen = None::<(usize, u32, usize)>;
            for &position in occurrences {
                let span_index = anchors.span_index(position);
                let text_tail = position as usize + run_k;
                let limit = (read.len() - read_tail)
                    .min(anchors.span_ends[span_index] as usize - text_tail);
                let length = run_k
                    + packed_common_prefix(
                        packed_read,
                        read_tail,
                        &anchors.packed,
                        text_tail,
                        limit,
                    );
                let span = self.spans[span_index];
                let reference_start = span.reference_start as usize
                    + (position - anchors.span_starts[span_index]) as usize;
                if chosen.is_none_or(|(best_length, sequence, start)| {
                    length > best_length
                        || (length == best_length
                            && (span.sequence, reference_start) < (sequence, start))
                }) {
                    chosen = Some((length, span.sequence, reference_start));
                }
            }
            let Some((length, sequence, reference_start)) = chosen else {
                continue;
            };
            profile.mem_starts += 1;
            profile.mem_bases += length as u64;
            let candidate = MemSeed {
                sequence,
                read_start,
                read_end: read_start + length,
                reference_start,
                reference_end: reference_start + length,
            };
            if best.is_none_or(|current| candidate_wins(candidate, current)) {
                best = Some(candidate);
            }
        }
        best
    }

    fn collect_fm(
        &self,
        read: &[u8],
        run_k: usize,
        occurrence_counts: &mut [u32],
        strand_masks: &mut [u8],
        profile: &mut MemQueryProfile,
    ) -> Option<MemSeed> {
        let mut best = None::<MemSeed>;
        for read_start in 0..occurrence_counts.len() {
            profile.run_windows += 1;
            // A match starting here can be at most this long. Once it is
            // shorter than the current best it can neither win nor tie, so
            // only the run-k window statistics are still needed.
            let cannot_win = best.is_some_and(|current| read.len() - read_start < current.len());
            let mut interval = (0_usize, self.length);
            let mut length = 0_usize;
            for &value in &read[read_start..] {
                let next = self.extend(interval, value + 1);
                if next.0 == next.1 {
                    break;
                }
                interval = next;
                length += 1;
                if length == run_k {
                    let occurrences = interval.1 - interval.0;
                    profile.matching_windows += 1;
                    occurrence_counts[read_start] = occurrences.min(u32::MAX as usize) as u32;
                    if occurrences == 1 {
                        strand_masks[read_start] = self.row_strands[interval.0];
                    }
                    if cannot_win {
                        break;
                    }
                }
            }
            if length < run_k || cannot_win {
                continue;
            }
            profile.mem_starts += 1;
            profile.mem_bases += length as u64;
            let Some((sequence, reference_start, reference_end)) =
                self.best_occurrence(interval, length)
            else {
                continue;
            };
            let candidate = MemSeed {
                sequence,
                read_start,
                read_end: read_start + length,
                reference_start,
                reference_end,
            };
            if best.is_none_or(|current| candidate_wins(candidate, current)) {
                best = Some(candidate);
            }
        }
        best
    }

    pub fn symbols(&self) -> usize {
        self.length
    }
}

/// Longer seeds win; equal lengths prefer (sequence, reference start, read start).
#[inline(always)]
fn candidate_wins(candidate: MemSeed, current: MemSeed) -> bool {
    candidate.len() > current.len()
        || (candidate.len() == current.len()
            && (
                candidate.sequence,
                candidate.reference_start,
                candidate.read_start,
            ) < (
                current.sequence,
                current.reference_start,
                current.read_start,
            ))
}

pub(crate) fn build_locus_indexes(
    references: &[OrientedReference],
    locus_count: usize,
    run_k: usize,
) -> Result<Vec<LocusMemIndex>, String> {
    build_locus_indexes_with(references, locus_count, run_k, run_k <= MAX_ANCHOR_K)
}

fn build_locus_indexes_with(
    references: &[OrientedReference],
    locus_count: usize,
    run_k: usize,
    use_anchors: bool,
) -> Result<Vec<LocusMemIndex>, String> {
    let mut loci = (0..locus_count)
        .map(|_| LocusMemIndex::default())
        .collect::<Vec<_>>();
    let mut texts = (0..locus_count).map(|_| Vec::new()).collect::<Vec<_>>();
    let mut span_texts = (0..locus_count)
        .map(|_| Vec::<Vec<u8>>::new())
        .collect::<Vec<_>>();
    for (sequence, reference) in references.iter().enumerate() {
        let sequence = u32::try_from(sequence)
            .map_err(|_| "reference sequence count exceeds u32".to_string())?;
        let locus_id = reference.locus as usize;
        let locus = &mut loci[locus_id];
        let text = &mut texts[locus_id];
        let mut cursor = 0_usize;
        while cursor < reference.bases.len() {
            while cursor < reference.bases.len() && code(reference.bases[cursor]).is_none() {
                cursor += 1;
            }
            let reference_start = cursor;
            while cursor < reference.bases.len() && code(reference.bases[cursor]).is_some() {
                cursor += 1;
            }
            let reference_end = cursor;
            if reference_end - reference_start < run_k {
                continue;
            }
            if !text.is_empty() {
                text.push(SEPARATOR);
            }
            let text_start = text.len();
            text.extend(
                reference.bases[reference_start..reference_end]
                    .iter()
                    .rev()
                    .map(|&base| base_symbol(base)),
            );
            if use_anchors {
                span_texts[locus_id].push(
                    reference.bases[reference_start..reference_end]
                        .iter()
                        .map(|&base| code(base).expect("span holds only bases"))
                        .collect(),
                );
            }
            locus.spans.push(ReferenceSpan {
                text_start: u32::try_from(text_start)
                    .map_err(|_| "locus text coordinate exceeds u32".to_string())?,
                text_end: u32::try_from(text.len())
                    .map_err(|_| "locus text coordinate exceeds u32".to_string())?,
                sequence,
                reference_start: u32::try_from(reference_start)
                    .map_err(|_| "reference coordinate exceeds u32".to_string())?,
                reference_end: u32::try_from(reference_end)
                    .map_err(|_| "reference coordinate exceeds u32".to_string())?,
                strand: reference.strand,
            });
        }
    }
    for ((locus, text), spans) in loci.iter_mut().zip(&mut texts).zip(&span_texts) {
        if text.is_empty() {
            continue;
        }
        text.push(SENTINEL);
        if use_anchors {
            // The FM text length is kept as the reported index size.
            locus.length = text.len();
            locus.anchors = Some(RunKAnchors::build(run_k, spans)?);
        } else {
            locus.build(text)?;
        }
    }
    Ok(loci)
}

#[inline(always)]
fn base_symbol(base: u8) -> u8 {
    code(base).map_or(SEPARATOR, |value| value + 1)
}

/// A suffix-tree node.
///
/// The incoming edge is represented by the inclusive range `start..=end`.
/// Leaves use `usize::MAX`, meaning their edge ends at the current global
/// leaf end.
#[derive(Clone, Copy)]
struct SuffixTreeNode {
    start: usize,
    end: usize,
    suffix_link: usize,
    first_child: usize,
}

impl SuffixTreeNode {
    fn root() -> Self {
        Self {
            start: 0,
            end: 0,
            suffix_link: 0,
            first_child: usize::MAX,
        }
    }

    fn new(start: usize, end: usize) -> Self {
        Self {
            start,
            end,
            suffix_link: 0,
            first_child: usize::MAX,
        }
    }
}

/// Compact adjacency-list entry.
///
/// Children are stored in one contiguous vector rather than using a
/// HashMap or separate allocation for every node.
#[derive(Clone, Copy)]
struct ChildLink {
    child: usize,
    next: usize,
}

/// Return the symbol at an augmented-text position.
///
/// Original bytes are mapped from 0..=255 to 1..=256. Position
/// `text.len()` is a virtual sentinel with value zero, so it is:
///
/// - distinct from every byte;
/// - lexicographically smaller than every byte;
/// - not stored in a copied augmented buffer.
#[inline]
fn suffix_symbol(text: &[u8], position: usize) -> u16 {
    if position == text.len() {
        0
    } else {
        u16::from(text[position]) + 1
    }
}

#[inline]
fn suffix_edge_length(nodes: &[SuffixTreeNode], node: usize, leaf_end: usize) -> usize {
    let end = if nodes[node].end == usize::MAX {
        leaf_end
    } else {
        nodes[node].end
    };

    end - nodes[node].start + 1
}

/// Locate an outgoing edge by its first symbol.
///
/// The alphabet contains only 257 possible symbols, including the sentinel.
/// Consequently, scanning the compact child list is bounded by a constant
/// independent of the text length.
fn find_child_link(
    text: &[u8],
    nodes: &[SuffixTreeNode],
    links: &[ChildLink],
    parent: usize,
    symbol: u16,
) -> usize {
    let mut link = nodes[parent].first_child;

    while link != usize::MAX {
        let child = links[link].child;
        let child_symbol = suffix_symbol(text, nodes[child].start);

        if child_symbol == symbol {
            return link;
        }

        link = links[link].next;
    }

    usize::MAX
}

fn add_child(
    nodes: &mut [SuffixTreeNode],
    links: &mut Vec<ChildLink>,
    parent: usize,
    child: usize,
) {
    let link = links.len();

    links.push(ChildLink {
        child,
        next: nodes[parent].first_child,
    });

    nodes[parent].first_child = link;
}

fn build_suffix_array(text: &[u8]) -> Result<Vec<usize>, String> {
    if text.len() > u32::MAX as usize {
        return Err("locus FM-index text length exceeds u32".to_string());
    }

    let length = text.len();

    if length == 0 {
        return Ok(Vec::new());
    }

    // Includes one virtual terminal symbol at position `length`.
    let total_length = length
        .checked_add(1)
        .ok_or_else(|| "locus FM-index text length exceeds address space".to_string())?;

    /*
     * Build a compact suffix tree using Ukkonen's algorithm.
     *
     * Node zero is always the root. Each non-root node represents both:
     *
     * - a suffix-tree node;
     * - the incoming edge from its parent.
     *
     * Edge strings are represented only by start/end positions.
     */
    let mut nodes = Vec::with_capacity(total_length);
    let mut links = Vec::with_capacity(total_length);

    nodes.push(SuffixTreeNode::root());

    let root = 0_usize;

    let mut active_node = root;
    let mut active_edge = 0_usize;
    let mut active_length = 0_usize;
    let mut remaining_suffixes = 0_usize;
    let mut leaf_end = 0_usize;

    for position in 0..total_length {
        leaf_end = position;
        remaining_suffixes += 1;

        // The most recently created internal node whose suffix link has
        // not yet been assigned during this phase.
        let mut last_internal = usize::MAX;

        while remaining_suffixes > 0 {
            if active_length == 0 {
                active_edge = position;
            }

            let active_symbol = suffix_symbol(text, active_edge);
            let edge_link = find_child_link(text, &nodes, &links, active_node, active_symbol);

            if edge_link == usize::MAX {
                /*
                 * Rule 2: no outgoing edge begins with the current symbol.
                 * Add a new leaf.
                 */
                let leaf = nodes.len();
                nodes.push(SuffixTreeNode::new(position, usize::MAX));
                add_child(&mut nodes, &mut links, active_node, leaf);

                if last_internal != usize::MAX {
                    nodes[last_internal].suffix_link = active_node;
                    last_internal = usize::MAX;
                }
            } else {
                let next = links[edge_link].child;
                let edge_length = suffix_edge_length(&nodes, next, leaf_end);

                /*
                 * Skip/count trick: move over an entire edge rather than
                 * comparing it one symbol at a time.
                 */
                if active_length >= edge_length {
                    active_edge += edge_length;
                    active_length -= edge_length;
                    active_node = next;
                    continue;
                }

                let next_symbol = suffix_symbol(text, nodes[next].start + active_length);
                let current_symbol = suffix_symbol(text, position);

                if next_symbol == current_symbol {
                    /*
                     * Rule 3: the current symbol is already present.
                     * The current phase ends early.
                     */
                    if last_internal != usize::MAX && active_node != root {
                        nodes[last_internal].suffix_link = active_node;
                    }

                    active_length += 1;
                    break;
                }

                /*
                 * Rule 2 with edge splitting.
                 *
                 * Existing:
                 *
                 * parent ---- next
                 *
                 * Becomes:
                 *
                 * parent ---- split ---- next
                 *                 \
                 *                  leaf
                 */
                let split_start = nodes[next].start;
                let split_end = split_start + active_length - 1;

                let split = nodes.len();
                nodes.push(SuffixTreeNode::new(split_start, split_end));

                // Replace parent -> next with parent -> split.
                links[edge_link].child = split;

                let leaf = nodes.len();
                nodes.push(SuffixTreeNode::new(position, usize::MAX));
                add_child(&mut nodes, &mut links, split, leaf);

                // Shorten the beginning of the original edge and attach
                // it below the new internal node.
                nodes[next].start += active_length;
                add_child(&mut nodes, &mut links, split, next);

                if last_internal != usize::MAX {
                    nodes[last_internal].suffix_link = split;
                }

                last_internal = split;
            }

            remaining_suffixes -= 1;

            /*
             * Move to the active point for the next pending suffix.
             */
            if active_node == root && active_length > 0 {
                active_length -= 1;
                active_edge = position + 1 - remaining_suffixes;
            } else if active_node != root {
                active_node = nodes[active_node].suffix_link;
            }
        }
    }

    /*
     * A lexicographic depth-first traversal of suffix-tree leaves gives
     * the suffix array.
     *
     * Use an explicit stack to avoid overflowing the call stack on highly
     * repetitive input.
     */
    let mut suffix_array = Vec::with_capacity(length);
    let mut stack = Vec::with_capacity(nodes.len().min(1024));

    // Each stack entry is `(node, string_depth)`.
    stack.push((root, 0_usize));

    while let Some((node, depth)) = stack.pop() {
        if nodes[node].first_child == usize::MAX {
            if node != root {
                let suffix_start = total_length - depth;

                // Exclude the suffix containing only the virtual sentinel.
                if suffix_start < length {
                    suffix_array.push(suffix_start);
                }
            }

            continue;
        }

        /*
         * Child links are stored in insertion order. Sort each node's
         * outgoing edges by their leading symbol before traversal.
         *
         * A suffix-tree node has at most 257 children for byte input.
         */
        let mut children = Vec::new();
        let mut link = nodes[node].first_child;

        while link != usize::MAX {
            children.push(links[link].child);
            link = links[link].next;
        }

        children.sort_unstable_by_key(|&child| suffix_symbol(text, nodes[child].start));

        // The stack is LIFO, so push children in reverse order.
        for child in children.into_iter().rev() {
            let child_depth = depth + suffix_edge_length(&nodes, child, leaf_end);

            stack.push((child, child_depth));
        }
    }

    debug_assert_eq!(suffix_array.len(), length);

    Ok(suffix_array)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn codes(read: &[u8]) -> Vec<u8> {
        let mut values = Vec::new();
        read_codes_into(read, &mut values);
        values
    }

    fn text_index(texts: &[(&[u8], u8)]) -> LocusMemIndex {
        let mut index = LocusMemIndex::default();
        let mut text = Vec::new();
        for (sequence, &(bases, strand)) in texts.iter().enumerate() {
            if !text.is_empty() {
                text.push(SEPARATOR);
            }
            let start = text.len();
            text.extend(bases.iter().rev().map(|&base| base_symbol(base)));
            index.spans.push(ReferenceSpan {
                text_start: start as u32,
                text_end: text.len() as u32,
                sequence: sequence as u32,
                reference_start: 0,
                reference_end: bases.len() as u32,
                strand,
            });
        }
        text.push(SENTINEL);
        index.build(&text).unwrap();
        index
    }

    #[test]
    fn fm_intervals_count_repeats_and_locate_leftmost() {
        let index = text_index(&[(b"ACGTACGTACGA", 1)]);
        let interval = index.interval(b"ACG");
        assert_eq!(interval.1 - interval.0, 3);
        assert_eq!(index.best_occurrence(interval, 3), Some((0, 0, 3)));
    }

    #[test]
    fn longest_prefix_returns_leftmost_equal_occurrence() {
        let index = text_index(&[(b"TTACGTTTACGT", 1)]);
        let read = codes(b"ACGTAA");
        let mut counts = vec![0; read.len() - 3 + 1];
        let mut strands = vec![0; counts.len()];
        let seed = index
            .collect(
                &read,
                3,
                &mut counts,
                &mut strands,
                &mut MemQueryProfile::default(),
            )
            .unwrap();
        assert_eq!((seed.read_start, seed.read_end), (0, 4));
        assert_eq!((seed.reference_start, seed.reference_end), (2, 6));
    }

    /// The pre-pruning `collect`: extends every start to its full MEM.
    fn collect_unpruned(
        index: &LocusMemIndex,
        read: &[u8],
        run_k: usize,
        occurrence_counts: &mut [u32],
        strand_masks: &mut [u8],
    ) -> Option<MemSeed> {
        let mut best = None::<MemSeed>;
        for read_start in 0..occurrence_counts.len() {
            let mut interval = (0_usize, index.length);
            let mut length = 0_usize;
            for &value in &read[read_start..] {
                let next = index.extend(interval, value + 1);
                if next.0 == next.1 {
                    break;
                }
                interval = next;
                length += 1;
                if length == run_k {
                    let occurrences = interval.1 - interval.0;
                    occurrence_counts[read_start] = occurrences as u32;
                    if occurrences == 1 {
                        if let Some((_, _, _, strand)) = index.occurrence(interval.0, run_k) {
                            strand_masks[read_start] = strand;
                        }
                    }
                }
            }
            if length < run_k {
                continue;
            }
            let Some((sequence, reference_start, reference_end)) =
                index.best_occurrence(interval, length)
            else {
                continue;
            };
            let candidate = MemSeed {
                sequence,
                read_start,
                read_end: read_start + length,
                reference_start,
                reference_end,
            };
            if best.is_none_or(|current| {
                candidate.len() > current.len()
                    || (candidate.len() == current.len()
                        && (
                            candidate.sequence,
                            candidate.reference_start,
                            candidate.read_start,
                        ) < (
                            current.sequence,
                            current.reference_start,
                            current.read_start,
                        ))
            }) {
                best = Some(candidate);
            }
        }
        best
    }

    #[test]
    fn pruned_collect_matches_full_extension() {
        // Small alphabet-biased references and reads drawn from them with
        // mutations make repeats and equal-length MEM ties common.
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = |bound: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % bound as u64) as usize
        };
        for _ in 0..400 {
            let unit: Vec<u8> = (0..1 + next(12)).map(|_| b"ACGT"[next(4)]).collect();
            let references: Vec<Vec<u8>> = (0..1 + next(3))
                .map(|_| {
                    (0..20 + next(200))
                        .map(|i| {
                            if next(3) == 0 {
                                b"ACGT"[next(4)]
                            } else {
                                unit[i % unit.len()]
                            }
                        })
                        .collect()
                })
                .collect();
            let texts: Vec<(&[u8], u8)> = references
                .iter()
                .enumerate()
                .map(|(i, r)| (r.as_slice(), 1 + (i % 2) as u8))
                .collect();
            let index = text_index(&texts);
            for _ in 0..20 {
                let source = &references[next(references.len())];
                let start = next(source.len());
                let mut read: Vec<u8> = source[start..]
                    .iter()
                    .take(10 + next(60))
                    .copied()
                    .collect();
                for base in read.iter_mut() {
                    if next(15) == 0 {
                        *base = b"ACGT"[next(4)];
                    }
                }
                let run_k = 1 + next(9);
                if read.len() < run_k {
                    continue;
                }
                let read = codes(&read);
                let windows = read.len() - run_k + 1;
                let (mut counts, mut strands) = (vec![0; windows], vec![0; windows]);
                let pruned = index.collect(
                    &read,
                    run_k,
                    &mut counts,
                    &mut strands,
                    &mut MemQueryProfile::default(),
                );
                let (mut full_counts, mut full_strands) = (vec![0; windows], vec![0; windows]);
                let full =
                    collect_unpruned(&index, &read, run_k, &mut full_counts, &mut full_strands);
                assert_eq!(pruned, full);
                assert_eq!(counts, full_counts);
                assert_eq!(strands, full_strands);
            }
        }
    }

    #[test]
    fn sentinels_prevent_cross_sequence_matches() {
        let index = text_index(&[(b"AAAAC", 1), (b"GTTTT", 2)]);
        let interval = index.interval(b"ACG");
        assert_eq!(interval.0, interval.1);
    }

    fn query(
        index: &LocusMemIndex,
        read: &[u8],
        run_k: usize,
    ) -> (Option<MemSeed>, Vec<u32>, Vec<u8>, MemQueryProfile) {
        let windows = read.len() - run_k + 1;
        let (mut counts, mut strands) = (vec![0; windows], vec![0; windows]);
        let mut profile = MemQueryProfile::default();
        let seed = index.collect(read, run_k, &mut counts, &mut strands, &mut profile);
        (seed, counts, strands, profile)
    }

    #[test]
    fn anchored_collect_matches_fm_index() {
        // References drawn from short repeat units with ambiguous bases make
        // repeated run-k windows, split spans and equal-length ties common.
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut next = |bound: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % bound as u64) as usize
        };
        let mut compared = 0_usize;
        for iteration in 0..300 {
            let run_k = [1, 2, 5, 11, 19, 19, 23, 31, 32][next(9)];
            let unit: Vec<u8> = (0..1 + next(16)).map(|_| b"ACGT"[next(4)]).collect();
            let loci = 1 + next(3);
            let mut references = Vec::new();
            for locus in 0..loci {
                for _ in 0..1 + next(3) {
                    let bases: Vec<u8> = (0..run_k + next(300))
                        .map(|i| match next(40) {
                            0 => b'N',
                            1..=12 => b"ACGT"[next(4)],
                            _ => unit[i % unit.len()],
                        })
                        .collect();
                    let reverse = crate::index::reverse_complement(&bases);
                    references.push(OrientedReference {
                        locus: locus as u32,
                        strand: 1,
                        bases,
                    });
                    references.push(OrientedReference {
                        locus: locus as u32,
                        strand: 2,
                        bases: reverse,
                    });
                }
            }
            let fm = build_locus_indexes_with(&references, loci, run_k, false).unwrap();
            let anchored = build_locus_indexes_with(&references, loci, run_k, true).unwrap();
            for locus in 0..loci {
                assert_eq!(fm[locus].symbols(), anchored[locus].symbols());
            }
            for _ in 0..40 {
                let source = &references[next(references.len())];
                let start = next(source.bases.len());
                let mut read: Vec<u8> = source.bases[start..]
                    .iter()
                    .take(run_k + next(160))
                    .map(|&base| if base == b'N' { b"ACGT"[next(4)] } else { base })
                    .collect();
                for base in read.iter_mut() {
                    if next(25) == 0 {
                        *base = b"ACGTacgu"[next(8)];
                    }
                }
                if read.len() < run_k {
                    continue;
                }
                let read = codes(&read);
                let locus = next(loci);
                let (fm_seed, fm_counts, fm_strands, fm_profile) = query(&fm[locus], &read, run_k);
                let (seed, counts, strands, profile) = query(&anchored[locus], &read, run_k);
                assert_eq!(seed, fm_seed, "iteration={iteration} run_k={run_k}");
                assert_eq!(counts, fm_counts, "iteration={iteration} run_k={run_k}");
                assert_eq!(strands, fm_strands, "iteration={iteration} run_k={run_k}");
                assert_eq!(profile.mem_starts, fm_profile.mem_starts);
                assert_eq!(profile.mem_bases, fm_profile.mem_bases);
                compared += 1;
            }
        }
        assert!(compared > 5_000);
    }

    #[test]
    fn packed_prefix_matches_bytewise_comparison() {
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = |bound: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % bound as u64) as usize
        };
        for _ in 0..2_000 {
            let text: Vec<u8> = (0..1 + next(200)).map(|_| next(4) as u8).collect();
            let start = next(text.len());
            let mut read: Vec<u8> = text[start..].to_vec();
            read.truncate(next(read.len() + 1));
            if !read.is_empty() && next(2) == 0 {
                let at = next(read.len());
                read[at] = (read[at] + 1 + next(3) as u8) % 4;
            }
            read.extend((0..next(40)).map(|_| next(4) as u8));
            let mut packed_text = vec![0; text.len().div_ceil(32) + 1];
            pack_codes(&text, &mut packed_text);
            let mut packed_read = vec![0; read.len().div_ceil(32) + 1];
            pack_codes(&read, &mut packed_read);
            let read_start = next(read.len() + 1);
            let limit = (read.len() - read_start).min(text.len() - start);
            let expected = read[read_start..]
                .iter()
                .zip(&text[start..])
                .take(limit)
                .take_while(|(left, right)| left == right)
                .count();
            assert_eq!(
                packed_common_prefix(&packed_read, read_start, &packed_text, start, limit),
                expected
            );
        }
    }
}
