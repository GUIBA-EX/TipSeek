//! Conservative paired-read audit on original contigs, not on synthetic N junctions.
//! Full-length ungapped placements only. Absence of bridging short reads is NOT
//! evidence against a long intron. This is not a general-purpose read aligner.
use super::{reverse_complement, GeneModel, Strand};
use crate::fastx::{FastxFormat, FastxReader};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufReader, Read};
use std::path::{Path, PathBuf};

const K: usize = 15;
const MIN_READ: usize = 90;
const MAX_EDIT: usize = 5;

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
struct Hit {
    contig: usize,
    start: usize,
    end: usize,
    reverse: bool,
    edits: usize,
    fragment_start: isize,
}

struct Index<'a> {
    sequences: Vec<&'a [u8]>,
    seeds: HashMap<u32, Vec<(usize, usize)>>,
}
fn key(sequence: &[u8]) -> Option<u32> {
    sequence.iter().try_fold(0, |k, b| {
        Some(
            (k << 2)
                | match b {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => return None,
                },
        )
    })
}
impl<'a> Index<'a> {
    fn new(sequences: Vec<&'a [u8]>) -> Self {
        let mut seeds: HashMap<u32, Vec<(usize, usize)>> = HashMap::new();
        for (id, sequence) in sequences.iter().enumerate() {
            for (offset, window) in sequence.windows(K).enumerate() {
                if let Some(k) = key(window) {
                    let positions = seeds.entry(k).or_default();
                    // 65 is a sentinel: repetitive seeds are never partially used.
                    if positions.len() <= 64 {
                        positions.push((id, offset));
                    }
                }
            }
        }
        Self { sequences, seeds }
    }
    #[cfg(test)]
    fn place(&self, read: &[u8]) -> Option<Hit> {
        self.place_evidence(read)
            .and_then(|(hit, partial)| (!partial).then_some(hit))
    }

    fn place_evidence(&self, read: &[u8]) -> Option<(Hit, bool)> {
        // N consumes an edit; it is neither deleted nor treated as a wildcard.
        // Six disjoint seeds still find every placement within five edits.
        if read.len() < MIN_READ || !read.iter().all(|b| b"ACGTN".contains(b)) {
            return None;
        }
        let reverse = reverse_complement(std::str::from_utf8(read).ok()?);
        let mut hits = Vec::new();
        for (is_reverse, query) in [(false, read), (true, reverse.as_bytes())] {
            let mut placements = HashSet::new();
            // Six disjoint seeds suffice: with at most five substitutions,
            // at least one must be exact. All competing placements within the
            // edit margin remain discoverable without scanning every window.
            for (block, window) in query.chunks_exact(K).take(MAX_EDIT + 1).enumerate() {
                let offset = block * K;
                if let Some(positions) = key(window).and_then(|k| self.seeds.get(&k)) {
                    // Never claim uniqueness after ignoring a repetitive seed.
                    if positions.len() > 64 {
                        return None;
                    }
                    for &(id, pos) in positions {
                        if pos >= offset {
                            placements.insert((id, pos - offset));
                        }
                    }
                }
                if placements.len() > 256 {
                    return None;
                }
            }
            for (id, start) in placements {
                let Some(target) = self.sequences[id].get(start..start + query.len()) else {
                    continue;
                };
                let edits = target
                    .iter()
                    .zip(query)
                    .filter(|(a, b)| a != b || !b"ACGT".contains(a))
                    .take(MAX_EDIT + 1)
                    .count();
                if edits <= MAX_EDIT {
                    hits.push(Hit {
                        contig: id,
                        start,
                        end: start + query.len(),
                        reverse: is_reverse,
                        edits,
                        fragment_start: start as isize,
                    });
                }
            }
        }
        hits.sort_by_key(|h| h.edits);
        let Some(best) = hits.first() else {
            return self.place_contained(read).map(|hit| (hit, true));
        };
        if best.edits >= MAX_EDIT {
            return None;
        }
        if hits.get(1).is_some_and(|h| h.edits <= best.edits + 1) {
            return None;
        }
        // At most five substitutions in >=90 bases guarantee an exact 15-mer.
        // Shorter reads cannot guarantee sensitivity/uniqueness at five edits.
        if read.len() < K * (MAX_EDIT + 1) {
            return None;
        }
        Some((best.clone(), false))
    }

    // A short observed contig can be wholly contained inside a longer read.
    // Exact >=30-nt containment, every competing placement, no pair linkage.
    fn place_contained(&self, read: &[u8]) -> Option<Hit> {
        let reverse = reverse_complement(std::str::from_utf8(read).ok()?);
        let mut placements = HashSet::new();
        for (is_reverse, query) in [(false, read), (true, reverse.as_bytes())] {
            for (block, window) in query.chunks_exact(K).enumerate() {
                let offset = block * K;
                let Some(positions) = key(window).and_then(|k| self.seeds.get(&k)) else {
                    continue;
                };
                if positions.len() > 64 {
                    return None;
                }
                for &(id, pos) in positions {
                    let target = self.sequences[id];
                    if target.len() < 30 || target.len() >= query.len() || pos > offset {
                        continue;
                    }
                    let start = offset - pos;
                    if query.get(start..start + target.len()) == Some(target) {
                        placements.insert((id, start, is_reverse));
                    }
                }
                if placements.len() > 1 {
                    return None;
                }
            }
        }
        if placements.len() != 1 {
            return None;
        }
        let (id, query_start, reverse) = placements.into_iter().next()?;
        Some(Hit {
            contig: id,
            start: 0,
            end: self.sequences[id].len(),
            reverse,
            edits: 0,
            fragment_start: -(query_start as isize),
        })
    }
}

enum Reads {
    Gm2(BufReader<File>),
    Fastq(FastxReader),
}
impl Reads {
    fn open(path: &Path) -> io::Result<Self> {
        if path.extension().is_some_and(|x| x == "gm2") {
            Ok(Self::Gm2(BufReader::new(File::open(path)?)))
        } else {
            Ok(Self::Fastq(FastxReader::open(path, FastxFormat::Fastq)?))
        }
    }
    fn next(&mut self) -> io::Result<Option<(String, Vec<u8>)>> {
        match self {
            Self::Fastq(reader) => Ok(reader.next_record()?.map(|r| {
                let name = String::from_utf8_lossy(&r.header)
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .trim_start_matches('@')
                    .to_owned();
                let name = name
                    .strip_suffix("/1")
                    .or_else(|| name.strip_suffix("/2"))
                    .unwrap_or(&name)
                    .to_owned();
                (name, r.sequence.to_ascii_uppercase())
            })),
            Self::Gm2(reader) => {
                let mut header = [0u8; 6];
                if reader.read(&mut header[..1])? == 0 {
                    return Ok(None);
                }
                reader.read_exact(&mut header[1..])?;
                let size =
                    ((header[0] as usize) << 16) | ((header[1] as usize) << 8) | header[2] as usize;
                let length = (((header[3] & 127) as usize) << 16)
                    | ((header[4] as usize) << 8)
                    | header[5] as usize;
                if size == 0 || length == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "empty GM2 pair record",
                    ));
                }
                let mut payload = vec![0; size];
                reader.read_exact(&mut payload)?;
                let mut sequence = Vec::with_capacity(length);
                let mut last = 0;
                for byte in payload {
                    let delta = byte & 127;
                    for _ in 0..if byte & 128 == 0 { 1 } else { 2 } {
                        let chunk = delta ^ last;
                        for _ in 0..((chunk >> 5) & 3) + 1 {
                            if sequence.len() < length {
                                sequence.push((chunk & 31) + 64);
                            }
                        }
                        last = chunk;
                    }
                    if sequence.len() == length {
                        return Ok(Some((String::new(), sequence)));
                    }
                }
                Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated GM2 sequence",
                ))
            }
        }
    }
}

#[derive(Default)]
pub struct PairAudit {
    pub available: bool,
    pub total_pairs: usize,
    pub unique_mapped_pairs: usize,
    pub distinct_cross_contig_pairs: usize,
    names: BTreeMap<String, usize>,
    pairs: HashMap<(usize, usize), Vec<(Hit, Hit)>>,
    single_hits: HashMap<usize, Vec<Hit>>,
    pub inputs: Vec<PathBuf>,
}
#[derive(Default, Debug)]
pub struct EdgeEvidence {
    pub concordant: usize,
    pub discordant: usize,
    pub orientation_conflicts: usize,
}
impl EdgeEvidence {
    pub fn warning(&self) -> bool {
        self.discordant >= 5 && self.discordant > self.concordant
    }
    pub fn status(&self, available: bool) -> &'static str {
        if !available {
            "reads_unavailable"
        } else if self.warning() {
            "discordant_read_warning"
        } else if self.concordant >= 2 {
            "paired_read_support"
        } else {
            "insufficient_read_evidence"
        }
    }
}
impl PairAudit {
    pub fn load(
        root: Option<&Path>,
        sample: &str,
        family: &str,
        contigs: &BTreeMap<String, (String, String)>,
    ) -> io::Result<Self> {
        let mut audit = Self::default();
        let Some(root) = root else { return Ok(audit) };
        let directory = root.join(sample).join("filtered_pe");
        let mut paths = None;
        for extension in ["gm2", "fq", "fastq", "fq.gz", "fastq.gz"] {
            let a = directory.join(format!("{family}_1.{extension}"));
            let b = directory.join(format!("{family}_2.{extension}"));
            if a.exists() != b.exists() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "missing mate in exon read audit",
                ));
            }
            if a.exists() {
                paths = Some((a, b));
                break;
            }
        }
        let Some((a, b)) = paths else {
            return Ok(audit);
        };
        let index = Index::new(contigs.values().map(|(_, s)| s.as_bytes()).collect());
        audit.names = contigs
            .values()
            .enumerate()
            .map(|(i, (name, _))| (name.clone(), i))
            .collect();
        audit.inputs = vec![a.clone(), b.clone()];
        let (mut left, mut right) = (Reads::open(&a)?, Reads::open(&b)?);
        let mut seen = HashSet::new();
        let mut single_seen = HashSet::new();
        loop {
            let pair = (left.next()?, right.next()?);
            let (a, b) = match pair {
                (None, None) => break,
                (Some(a), Some(b)) => (a, b),
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "exon audit mate counts differ",
                    ))
                }
            };
            if a.0 != b.0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "exon audit mate identifiers differ",
                ));
            }
            audit.total_pairs += 1;
            let (a, b) = (index.place_evidence(&a.1), index.place_evidence(&b.1));
            for (hit, _) in [&a, &b].into_iter().flatten() {
                if hit.edits <= 1
                    && single_seen.insert((hit.contig, hit.fragment_start, hit.end, hit.reverse))
                {
                    audit
                        .single_hits
                        .entry(hit.contig)
                        .or_default()
                        .push(hit.clone());
                }
            }
            if let (Some((a, false)), Some((b, false))) = (a, b) {
                audit.unique_mapped_pairs += 1;
                if a.contig == b.contig {
                    continue;
                }
                let (a, b) = if a.contig < b.contig { (a, b) } else { (b, a) };
                // Coordinate deduplication: PCR-identical fragments are not independent support.
                let coordinate = (
                    a.contig, a.start, a.end, a.reverse, b.contig, b.start, b.end, b.reverse,
                );
                if seen.insert(coordinate) {
                    audit.distinct_cross_contig_pairs += 1;
                    audit
                        .pairs
                        .entry((a.contig, b.contig))
                        .or_default()
                        .push((a, b));
                }
            }
        }
        audit.available = true;
        Ok(audit)
    }

    /// Individual uniquely placed mates can support an exon even when their
    /// partner cannot be placed. This does not establish a cross-contig link.
    /// Count coordinate-deduplicated reads overlapping at least 30 nt (or the
    /// whole segment when shorter), and the union of their mapped bases.
    pub fn segment_support(&self, candidate: &str, segment: super::Interval) -> (usize, usize) {
        let Some(hits) = self
            .names
            .get(candidate)
            .and_then(|id| self.single_hits.get(id))
        else {
            return (0, 0);
        };
        let mut intervals: Vec<_> = hits
            .iter()
            .filter_map(|hit| {
                let start = hit.start.max(segment.start);
                let end = hit.end.min(segment.end);
                (end.saturating_sub(start) >= 30.min(segment.len())).then_some((start, end))
            })
            .collect();
        intervals.sort_unstable();
        let mut end = segment.start;
        let mut covered = 0;
        for &(a, b) in &intervals {
            covered += b.saturating_sub(a.max(end));
            end = end.max(b);
        }
        (intervals.len(), covered)
    }
    pub fn edge(&self, left: &GeneModel, right: &GeneModel) -> EdgeEvidence {
        let mut evidence = EdgeEvidence::default();
        let (Some(&a), Some(&b)) = (
            self.names.get(&left.candidate_id),
            self.names.get(&right.candidate_id),
        ) else {
            return evidence;
        };
        let key = if a < b { (a, b) } else { (b, a) };
        let Some(pairs) = self.pairs.get(&key) else {
            return evidence;
        };
        for (first, second) in pairs {
            let (l, r) = if a < b {
                (first, second)
            } else {
                (second, first)
            };
            let within_exon = |hit: &Hit, model: &GeneModel| {
                model
                    .exons
                    .iter()
                    .any(|e| e.interval.start <= hit.start && hit.end <= e.interval.end)
            };
            if !within_exon(l, left) || !within_exon(r, right) {
                continue;
            }
            if l.reverse != (left.strand == Strand::Reverse)
                || r.reverse != (right.strand == Strand::Forward)
            {
                evidence.orientation_conflicts += 1;
                continue;
            }
            if l.edits <= 1 && r.edits <= 1 {
                evidence.concordant += 1;
            } else if (l.edits == 0 && r.edits >= 3) || (r.edits == 0 && l.edits >= 3) {
                evidence.discordant += 1;
            }
        }
        evidence
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn dna(n: usize, mut seed: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                b"ACGT"[(seed % 4) as usize]
            })
            .collect()
    }
    #[test]
    fn exact_short_contig_containment_supports_sequence_without_full_read_linkage() {
        let sequence = dna(180, 991);
        let short = &sequence[70..115];
        let index = Index::new(vec![short]);
        assert!(index.place(&sequence[..150]).is_none());
        let (hit, partial) = index.place_evidence(&sequence[..150]).unwrap();
        assert!(partial);
        assert_eq!((hit.start, hit.end, hit.fragment_start), (0, 45, -70));
        assert!(
            index
                .place_evidence(
                    &reverse_complement(std::str::from_utf8(&sequence[..150]).unwrap())
                        .into_bytes()
                )
                .unwrap()
                .0
                .reverse
        );
        assert!(Index::new(vec![short, short])
            .place_evidence(&sequence[..150])
            .is_none());
        assert!(Index::new(vec![&short[..29]])
            .place_evidence(&sequence[..150])
            .is_none());
    }
    #[test]
    fn disjoint_seeds_match_exhaustive_placement_with_substitutions() {
        let target = dna(360, 127);
        let other = dna(360, 7919);
        let index = Index::new(vec![&target, &other]);
        for length in [90, 120, 150] {
            for errors in 0..=5 {
                for shift in 0..15 {
                    let mut read = target[60..60 + length].to_vec();
                    for block in 0..errors {
                        let pos = block * K + shift;
                        read[pos] = if read[pos] == b'A' { b'C' } else { b'A' };
                    }
                    let rc = reverse_complement(std::str::from_utf8(&read).unwrap());
                    let mut brute = Vec::new();
                    for (reverse, query) in [(false, read.as_slice()), (true, rc.as_bytes())] {
                        for (contig, sequence) in [&target, &other].into_iter().enumerate() {
                            for (start, window) in sequence.windows(length).enumerate() {
                                let edits =
                                    window.iter().zip(query).filter(|(a, b)| a != b).count();
                                if edits <= MAX_EDIT {
                                    brute.push(Hit {
                                        contig,
                                        start,
                                        end: start + length,
                                        reverse,
                                        edits,
                                        fragment_start: start as isize,
                                    });
                                }
                            }
                        }
                    }
                    brute.sort_by_key(|h| h.edits);
                    let expected = brute
                        .first()
                        .filter(|best| {
                            best.edits < MAX_EDIT
                                && brute
                                    .get(1)
                                    .is_none_or(|runner| runner.edits > best.edits + 1)
                        })
                        .cloned();
                    assert_eq!(index.place(&read), expected);
                }
            }
        }
    }

    #[test]
    fn ambiguous_read_bases_consume_edits_without_inventing_unique_support() {
        let target = dna(360, 127);
        let index = Index::new(vec![&target]);
        for count in 1..=5 {
            for shift in 0..15 {
                let mut read = target[60..180].to_vec();
                for block in 0..count {
                    read[block * K + shift] = b'N';
                }
                let mapped = index.place(&read);
                if count < MAX_EDIT {
                    let hit = mapped.unwrap();
                    assert_eq!((hit.start, hit.edits), (60, count));
                    let reverse = reverse_complement(std::str::from_utf8(&read).unwrap());
                    assert!(index.place(reverse.as_bytes()).unwrap().reverse);
                } else {
                    assert!(mapped.is_none());
                }
            }
        }
        let mut read = target[60..180].to_vec();
        read[20] = b'N';
        let mut competing = target.clone();
        competing[80] = if competing[80] == b'A' { b'C' } else { b'A' };
        assert!(Index::new(vec![&target, &competing]).place(&read).is_none());
        // An N outside a known short exon may support its exact sequence,
        // but partial containment must never establish paired linkage.
        let short = &target[90..135];
        let index = Index::new(vec![short]);
        assert!(index.place_evidence(&read).unwrap().1);
        read[40] = b'N';
        assert!(index.place_evidence(&read).is_none());
    }

    #[test]
    fn mapping_is_full_length_strand_aware_and_copy_ambiguous() {
        let sequence = dna(400, 127);
        let index = Index::new(vec![&sequence]);
        let read = &sequence[80..200];
        assert_eq!(index.place(read).unwrap().start, 80);
        assert!(
            index
                .place(reverse_complement(std::str::from_utf8(read).unwrap()).as_bytes())
                .unwrap()
                .reverse
        );
        assert!(Index::new(vec![&sequence, &sequence]).place(read).is_none());
        let mut mutated = read.to_vec();
        for i in [5, 26, 60] {
            mutated[i] = if mutated[i] == b'A' { b'C' } else { b'A' };
        }
        assert_eq!(index.place(&mutated).unwrap().edits, 3);
        assert!(index.place(&read[..60]).is_none());
    }
    #[test]
    fn no_pairs_is_not_a_chimera_and_warning_requires_repeated_discordance() {
        assert_eq!(
            EdgeEvidence::default().status(true),
            "insufficient_read_evidence"
        );
        assert_eq!(EdgeEvidence::default().status(false), "reads_unavailable");
        let e = EdgeEvidence {
            concordant: 0,
            discordant: 5,
            orientation_conflicts: 0,
        };
        assert!(e.warning());
        assert!(!(EdgeEvidence { concordant: 5, ..e }).warning());
    }
    fn model(name: &str, sequence: &str) -> GeneModel {
        let length = sequence.len() / 3;
        let hit = super::super::exonerate::parse_vulgar(
            &format!(
                "vulgar: p 0 {length} . {name} 0 {} + 100 M {length} {}",
                sequence.len(),
                sequence.len()
            ),
            name,
            sequence.len(),
            &BTreeMap::from([("p".into(), length)]),
        )
        .unwrap()
        .remove(0);
        hit.normalize(name, sequence, &"A".repeat(length))
            .unwrap()
            .build(sequence, super::super::ModelConfig::default())
            .unwrap()
    }
    #[test]
    fn paired_files_deduplicate_fragments_and_require_exonic_inward_placements() {
        let root = std::env::temp_dir().join(format!("exon-pairs-{}", std::process::id()));
        let directory = root.join("s/filtered_pe");
        std::fs::create_dir_all(&directory).unwrap();
        let a = String::from_utf8(dna(600, 411)).unwrap();
        let b = String::from_utf8(dna(600, 721)).unwrap();
        let mut left = String::new();
        let mut right = String::new();
        for i in 0..11 {
            let offset = if i == 10 { 10 } else { 0 };
            left += &format!(
                "@pair{i}/1\n{}\n+\n{}\n",
                &a[150 + offset..270 + offset],
                "I".repeat(120)
            );
            right += &format!(
                "@pair{i}/2\n{}\n+\n{}\n",
                reverse_complement(&b[300 + offset..420 + offset]),
                "I".repeat(120)
            );
        }
        let r1 = directory.join("f_1.fq");
        let r2 = directory.join("f_2.fq");
        std::fs::write(&r1, &left).unwrap();
        std::fs::write(&r2, &right).unwrap();
        let contigs = BTreeMap::from([
            ("1".into(), ("a".into(), a.clone())),
            ("2".into(), ("b".into(), b.clone())),
        ]);
        let audit = PairAudit::load(Some(&root), "s", "f", &contigs).unwrap();
        assert_eq!(audit.total_pairs, 11);
        assert_eq!(audit.unique_mapped_pairs, 11);
        assert_eq!(audit.distinct_cross_contig_pairs, 2);
        assert_eq!(
            audit.segment_support("a", super::super::Interval::new(160, 280).unwrap()),
            (2, 120)
        );
        let mut ambiguous_read = a[150..270].to_owned();
        ambiguous_read.replace_range(20..21, "N");
        std::fs::write(&r1, left.replace(&a[150..270], &ambiguous_read)).unwrap();
        let with_n = PairAudit::load(Some(&root), "s", "f", &contigs).unwrap();
        assert_eq!(with_n.unique_mapped_pairs, 11);
        assert_eq!(with_n.distinct_cross_contig_pairs, 2);
        assert_eq!(
            with_n.segment_support("a", super::super::Interval::new(160, 280).unwrap()),
            (2, 120)
        );
        std::fs::write(&r1, &left).unwrap();
        let mut ma = model("a", &a);
        let mb = model("b", &b);
        assert_eq!(audit.edge(&ma, &mb).status(true), "paired_read_support");
        ma.strand = Strand::Reverse;
        assert_eq!(audit.edge(&ma, &mb).orientation_conflicts, 2);
        ma.strand = Strand::Forward;
        ma.exons[0].interval.end = 220;
        assert_eq!(audit.edge(&ma, &mb).concordant, 0);
        let unmapped = right
            .replace(&reverse_complement(&b[300..420]), &"N".repeat(120))
            .replace(&reverse_complement(&b[310..430]), &"N".repeat(120));
        std::fs::write(&r2, unmapped).unwrap();
        let single = PairAudit::load(Some(&root), "s", "f", &contigs).unwrap();
        assert_eq!(single.unique_mapped_pairs, 0);
        assert_eq!(
            single.segment_support("a", super::super::Interval::new(160, 280).unwrap()),
            (2, 120)
        );
        std::fs::write(&r2, right.replace("@pair0/2", "@wrong/2")).unwrap();
        assert!(PairAudit::load(Some(&root), "s", "f", &contigs).is_err());
        std::fs::write(&r2, "").unwrap();
        assert!(PairAudit::load(Some(&root), "s", "f", &contigs).is_err());
        std::fs::remove_dir_all(root).unwrap();
    }
    #[test]
    fn repeated_discordance_is_counted_only_for_exonic_inward_pairs() {
        let a = String::from_utf8(dna(600, 51)).unwrap();
        let b = String::from_utf8(dna(600, 91)).unwrap();
        let mut audit = PairAudit {
            available: true,
            names: BTreeMap::from([("a".into(), 0), ("b".into(), 1)]),
            ..PairAudit::default()
        };
        audit.pairs.insert(
            (0, 1),
            (0..5)
                .map(|i| {
                    (
                        Hit {
                            contig: 0,
                            start: 100 + i,
                            end: 220 + i,
                            reverse: false,
                            edits: 0,
                            fragment_start: 0,
                        },
                        Hit {
                            contig: 1,
                            start: 200 + i,
                            end: 320 + i,
                            reverse: true,
                            edits: 3,
                            fragment_start: 0,
                        },
                    )
                })
                .collect(),
        );
        assert!(audit.edge(&model("a", &a), &model("b", &b)).warning());
    }
    #[test]
    fn gm2_sequence_decoder_handles_packed_deltas_and_truncation() {
        // A single chunk A followed by C, stored using the high-bit XOR shorthand.
        let path = std::env::temp_dir().join(format!("exon-reads-{}.gm2", std::process::id()));
        std::fs::write(&path, [0, 0, 2, 0, 0, 3, 1, 130]).unwrap();
        let mut reader = Reads::open(&path).unwrap();
        assert_eq!(reader.next().unwrap().unwrap().1, b"ACA");
        assert!(reader.next().unwrap().is_none());
        std::fs::write(&path, [0, 0, 1, 0, 0, 3, 1]).unwrap();
        assert!(Reads::open(&path).unwrap().next().is_err());
        std::fs::remove_file(path).unwrap();
    }
}
