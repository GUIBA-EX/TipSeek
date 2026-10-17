//! Cheap protein screen exclusively for newly rescued, nucleotide-unanchored
//! graph components. A 10-aa ungapped window with >=8 matches is sufficient.
//! This is candidate screening, never model acceptance.
use super::{reverse_complement, translate_cds};
use std::collections::HashMap;

pub struct ProteinScreen<'a> {
    proteins: Vec<&'a [u8]>,
    seeds: HashMap<[u8; 3], Vec<(usize, usize)>>,
}
fn informative(window: &[u8]) -> bool {
    !window.contains(&b'X')
        && !window.contains(&b'*')
        && window
            .iter()
            .fold(0u32, |mask, &b| mask | (1 << (b - b'A')))
            .count_ones()
            >= 4
}
impl<'a> ProteinScreen<'a> {
    pub fn new(proteins: impl Iterator<Item = &'a str>) -> Self {
        let proteins: Vec<_> = proteins.map(str::as_bytes).collect();
        let mut seeds: HashMap<_, Vec<_>> = HashMap::new();
        for (id, protein) in proteins.iter().enumerate() {
            let mut informative_starts = vec![false; protein.len()];
            for (start, window) in protein.windows(10).enumerate() {
                if informative(window) {
                    informative_starts[start..start + 8].fill(true);
                }
            }
            for (pos, seed) in protein.windows(3).enumerate() {
                if informative_starts[pos] {
                    seeds
                        .entry([seed[0], seed[1], seed[2]])
                        .or_default()
                        .push((id, pos));
                }
            }
        }
        Self { proteins, seeds }
    }
    fn protein_supported(&self, target: &[u8]) -> bool {
        if target.len() < 10 {
            return false;
        }
        let mut previous = HashMap::new();
        for (pos, seed) in target.windows(3).enumerate() {
            if let Some(hits) = self.seeds.get(&[seed[0], seed[1], seed[2]]) {
                for &(id, qpos) in hits {
                    let diagonal = qpos as isize - pos as isize;
                    if let Some(last) = previous.insert((id, diagonal), pos) {
                        if pos == last || pos - last > 7 {
                            continue;
                        }
                        for start in pos.saturating_sub(7)..=last.min(target.len() - 10) {
                            let query_start = start as isize + diagonal;
                            if query_start < 0 {
                                continue;
                            }
                            let query_start = query_start as usize;
                            let Some(query) = self.proteins[id].get(query_start..query_start + 10)
                            else {
                                continue;
                            };
                            let window = &target[start..start + 10];
                            if !informative(window) || !informative(query) {
                                continue;
                            }
                            if window.iter().zip(query).filter(|(a, b)| a == b).count() >= 8 {
                                return true;
                            }
                        }
                    }
                }
            }
        }
        false
    }
    /// Screening only for observed 30-59 nt candidates independently supported
    /// by reads. Exact informative eight-residue anchors must have one placement
    /// per reference and agree on the sample strand/frame. No model is accepted here.
    pub fn supports_read_backed_short(&self, dna: &str, count: usize, covered: usize) -> bool {
        count >= 2 && covered == dna.len() && self.supports_short_exact(dna)
    }

    pub fn supports_short_exact(&self, dna: &str) -> bool {
        if !(30..60).contains(&dna.len()) || !dna.bytes().all(|b| b"ACGT".contains(&b)) {
            return false;
        }
        let reverse = reverse_complement(dna);
        let mut orientation = None;
        let mut found = false;
        for (strand, oriented) in [dna, reverse.as_str()].into_iter().enumerate() {
            for frame in 0..3 {
                let length = (oriented.len() - frame) / 3 * 3;
                let Ok(translated) = translate_cds(&oriented[frame..frame + length]) else {
                    continue;
                };
                for protein in &self.proteins {
                    let mut diagonal = None;
                    for (start, word) in translated.as_bytes().windows(8).enumerate() {
                        if !word.iter().all(|b| b"ACDEFGHIKLMNPQRSTVWY".contains(b))
                            || word
                                .iter()
                                .fold(0u32, |mask, b| mask | (1 << (b - b'A')))
                                .count_ones()
                                < 5
                        {
                            continue;
                        }
                        for (qpos, target) in protein.windows(8).enumerate() {
                            if word != target {
                                continue;
                            }
                            let d = qpos as isize - start as isize;
                            if diagonal.is_some_and(|old| old != d)
                                || orientation.is_some_and(|old| old != (strand, frame))
                            {
                                return false;
                            }
                            diagonal = Some(d);
                            orientation = Some((strand, frame));
                            found = true;
                        }
                    }
                }
            }
        }
        found
    }

    pub fn supports(&self, dna: &str) -> bool {
        let reverse = reverse_complement(dna);
        for oriented in [dna, reverse.as_str()] {
            for frame in 0..3 {
                let length = oriented.len().saturating_sub(frame) / 3 * 3;
                if length < 30 {
                    continue;
                }
                if let Ok(protein) = translate_cds(&oriented[frame..frame + length]) {
                    if self.protein_supported(protein.as_bytes()) {
                        return true;
                    }
                }
            }
        }
        false
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn short_screen_requires_exact_unique_anchor_and_independent_full_coverage() {
        let dna = "TAAATGAAAGCTTGGGATGTTTTTCCTTAA";
        let screen = ProteinScreen::new(std::iter::once("GGMKAWDVFPHH"));
        assert!(!screen.supports(dna));
        assert!(screen.supports_read_backed_short(dna, 2, dna.len()));
        assert!(screen.supports_read_backed_short(&reverse_complement(dna), 2, dna.len()));
        assert!(!screen.supports_read_backed_short(dna, 1, dna.len()));
        assert!(!screen.supports_read_backed_short(dna, 2, dna.len() - 1));
        assert!(!screen.supports_short_exact(&dna.replace("ATG", "ANG")));
        assert!(!ProteinScreen::new(std::iter::once("MKAWDVFPMKAWDVFP")).supports_short_exact(dna));
        assert!(!ProteinScreen::new(std::iter::once("MMMMMMMM"))
            .supports_short_exact(&"ATG".repeat(10)));
    }
    #[test]
    fn every_two_substitution_pattern_in_a_short_exon_passes() {
        let protein = "MALWQRSNDE";
        let index = ProteinScreen::new(std::iter::once(protein));
        for a in 0..10 {
            for b in a + 1..10 {
                let mut target = protein.as_bytes().to_vec();
                target[a] = b'G';
                target[b] = b'G';
                assert!(index.protein_supported(&target), "{a},{b}");
            }
        }
        assert!(!index.protein_supported(b"MALGGGGGGG"));
        assert!(!index.protein_supported(b"MALWQRSND"));
        assert!(!ProteinScreen::new(std::iter::once("QQQQQQQQQQ")).protein_supported(b"QQQQQQQQQQ"));
    }
    #[test]
    fn six_frame_screen_preserves_reverse_coding_fragments() {
        let dna = "ATGGCTCTGTGGCAGCGTTCTAATGATGAA";
        let index = ProteinScreen::new(std::iter::once("MALWQRSNDE"));
        assert!(index.supports(dna));
        assert!(index.supports(&format!("A{}TT", reverse_complement(dna))));
        assert!(!index.supports(&"C".repeat(120)));
    }
}
