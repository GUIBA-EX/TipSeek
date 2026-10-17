//! Exact containment is evidence redundancy, not a claim of haplotype linkage.
//! Only a unique, full-length occurrence in a longer observed contig qualifies.
use super::{reverse_complement, GeneModel, Strand};
use std::collections::{BTreeMap, BTreeSet, HashMap};

#[derive(Default)]
pub(super) struct Equivalence<'a> {
    aliases: BTreeMap<&'a str, (&'a str, usize, usize, bool)>,
}
impl<'a> Equivalence<'a> {
    pub fn new(models: &[GeneModel], dna: &'a BTreeMap<String, String>) -> Self {
        let used: BTreeSet<_> = models
            .iter()
            .filter(|m| super::exon_recovery::eligible(m))
            .map(|m| m.candidate_id.as_str())
            .collect();
        let sequences: Vec<_> = dna
            .iter()
            .filter(|(id, s)| {
                used.contains(id.as_str())
                    && s.len() >= 30
                    && s.bytes().all(|b| b"ACGT".contains(&b))
            })
            .collect();
        // A saturated seed is ignored in its entirety; never claim uniqueness
        // after searching only part of a repetitive bucket.
        let mut seeds: HashMap<&[u8], Vec<(usize, usize)>> = HashMap::new();
        for (id, (_, s)) in sequences.iter().enumerate() {
            for (pos, seed) in s.as_bytes().windows(15).enumerate() {
                let bucket = seeds.entry(seed).or_default();
                if bucket.len() <= 64 {
                    bucket.push((id, pos));
                }
            }
        }
        let mut aliases = BTreeMap::new();
        for (id, sequence) in &sequences {
            let reverse = reverse_complement(sequence);
            let mut unique = None;
            let mut ambiguous = false;
            for (reversed, query) in [(false, sequence.as_str()), (true, reverse.as_str())] {
                let Some(bucket) = seeds.get(&query.as_bytes()[..15]) else {
                    continue;
                };
                if bucket.len() > 64 {
                    ambiguous = true;
                    break;
                }
                for &(parent, start) in bucket {
                    let (name, longer) = sequences[parent];
                    if longer.len() < query.len() + 30
                        || longer.get(start..start + query.len()) != Some(query)
                    {
                        continue;
                    }
                    if unique.is_some() {
                        ambiguous = true;
                        break;
                    }
                    unique = Some((name.as_str(), start, query.len(), reversed));
                }
                if ambiguous {
                    break;
                }
            }
            if !ambiguous {
                if let Some(placement) = unique {
                    aliases.insert(id.as_str(), placement);
                }
            }
        }
        Self { aliases }
    }

    pub fn map<'b>(
        &'b self,
        mut id: &'b str,
        mut strand: Strand,
        mut positions: [usize; 3],
    ) -> (&'b str, Strand, [usize; 3]) {
        // Parent lengths strictly increase, so aliases cannot form a cycle.
        while let Some(&(parent, start, length, reversed)) = self.aliases.get(id) {
            for p in &mut positions {
                *p = start + if reversed { length - 1 - *p } else { *p };
            }
            if reversed {
                strand = if strand == Strand::Forward {
                    Strand::Reverse
                } else {
                    Strand::Forward
                };
            }
            id = parent;
        }
        (id, strand, positions)
    }
}
