use crate::assembly::{
    add_pe_branch_support, add_read_slices, assemble_seed, branch_edges, build_graph_from_counts,
    compare_contigs, count_assemble_chunk_parallel, filter_and_weight_graph, Its2PathEvidence,
    SortedKmerCounts,
};
use crate::hash::{HashMap, HashSet};
use crate::io_utils::{
    find_filtered, for_each_sequence_chunk, has_interleaved_pairs, load_or_build_reference_kmers,
    minimum_sequence_length, read_fasta, read_linked_reads_with_quality, LinkedReads,
};
use crate::model::{Args, AssemblyMode, ContigRecord, GraphFormat, LocusResult, LocusTask};
use crate::seq::{calculate_auto_k, encode_kmer, reverse_complement, reverse_complement_kmer};
use crate::unitig::write_graphs;
use std::fs::{self, File, OpenOptions};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

pub const SUMMARY_HEADER: &str = "locus,status,accepted,rejection_reason,selected_contig_length,read_supported_span,slice_supported_bases,slice_support_breadth,max_slice_support_gap,read_count,unique_read_count,multi_mapping_read_count,read_density,unique_read_density,support_fraction,flank_balance,kmer_median_depth,kmer_depth_cv,kmer_max_depth_ratio,candidate_count,low_quality";

#[derive(Default)]
pub struct ProfileStats {
    processed_loci: AtomicU64,
    completed_loci: AtomicU64,
    failed_loci: AtomicU64,
    skipped_loci: AtomicU64,
    reads_decode_ns: AtomicU64,
    kmer_count_ns: AtomicU64,
    graph_filter_ns: AtomicU64,
    seed_path_ns: AtomicU64,
    qc_write_ns: AtomicU64,
}

impl ProfileStats {
    fn add(counter: &AtomicU64, started: Instant) {
        counter.fetch_add(
            started.elapsed().as_nanos().min(u64::MAX as u128) as u64,
            Ordering::Relaxed,
        );
    }
    fn add_duration(counter: &AtomicU64, elapsed: Duration) {
        counter.fetch_add(
            elapsed.as_nanos().min(u64::MAX as u128) as u64,
            Ordering::Relaxed,
        );
    }
    fn record_reads(&self, started: Instant) {
        Self::add(&self.reads_decode_ns, started);
    }
    fn record_reads_duration(&self, elapsed: Duration) {
        Self::add_duration(&self.reads_decode_ns, elapsed);
    }
    fn record_kmers(&self, started: Instant) {
        Self::add(&self.kmer_count_ns, started);
    }
    fn record_kmers_duration(&self, elapsed: Duration) {
        Self::add_duration(&self.kmer_count_ns, elapsed);
    }
    fn record_graph(&self, started: Instant) {
        Self::add(&self.graph_filter_ns, started);
    }
    fn record_path(&self, started: Instant) {
        Self::add(&self.seed_path_ns, started);
    }
    fn record_qc(&self, started: Instant) {
        Self::add(&self.qc_write_ns, started);
    }
    fn record_locus(&self, result: &LocusResult) {
        self.processed_loci.fetch_add(1, Ordering::Relaxed);
        if result.skipped {
            self.skipped_loci.fetch_add(1, Ordering::Relaxed);
        } else if matches!(result.status.as_str(), "success" | "low quality") {
            self.completed_loci.fetch_add(1, Ordering::Relaxed);
        } else {
            self.failed_loci.fetch_add(1, Ordering::Relaxed);
        }
    }
}

pub fn write_assembly_profile(output: &Path, profile: &ProfileStats) -> io::Result<()> {
    let mut out = BufWriter::new(File::create(output.join("assembly_profile.tsv"))?);
    writeln!(
        out,
        "stage	seconds	processed_loci	completed_loci	failed_loci	skipped_loci"
    )?;
    let processed = profile.processed_loci.load(Ordering::Relaxed);
    let completed = profile.completed_loci.load(Ordering::Relaxed);
    let failed = profile.failed_loci.load(Ordering::Relaxed);
    let skipped = profile.skipped_loci.load(Ordering::Relaxed);
    for (name, value) in [
        ("reads_decode", &profile.reads_decode_ns),
        ("kmer_count", &profile.kmer_count_ns),
        ("graph_filter_weight", &profile.graph_filter_ns),
        ("seed_path_search", &profile.seed_path_ns),
        ("contig_qc_write", &profile.qc_write_ns),
    ] {
        writeln!(
            out,
            "{name}	{:.6}	{processed}	{completed}	{failed}	{skipped}",
            value.load(Ordering::Relaxed) as f64 / 1e9
        )?;
    }
    out.flush()
}

// 多线程写同一个日志得排队，不然字儿搅一块就没法看了。
pub fn log_line(output: &Path, lock: &Mutex<()>, message: &str) {
    let _guard = lock.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    if let Ok(mut log) = OpenOptions::new()
        .create(true)
        .append(true)
        .open(output.join("log.txt"))
    {
        let _ = writeln!(log, "{message}");
    }
    println!("{message}");
}

fn remove_if_exists(path: &Path) {
    if path.is_file() {
        let _ = fs::remove_file(path);
    }
}

// 某个位点重跑或失败时，只扫它自己的旧产物，别误伤邻居。
fn clean_locus_outputs(args: &Args, key: &str) {
    remove_if_exists(&args.output.join("results").join(format!("{key}.fasta")));
    remove_if_exists(
        &args
            .output
            .join("results")
            .join(format!("{key}.its2_support.tsv")),
    );
    remove_if_exists(&args.output.join("contigs_all").join(format!("{key}.fasta")));
    remove_if_exists(
        &args
            .output
            .join("contigs_all_low")
            .join(format!("{key}.fasta")),
    );
}

fn format_header(contig: &ContigRecord, mode: AssemblyMode, prefix: &str) -> String {
    let mut header = format!(
        ">{prefix}_{}_{}_{}_{}_{}_span_{}",
        contig.sequence.len(),
        contig.seed_count,
        contig.position,
        contig.weight,
        contig.read_count,
        contig.supported_span
    );
    if mode == AssemblyMode::Uce {
        header.push_str(&format!("_supported_{}", contig.supported_bases));
    }
    if mode == AssemblyMode::Its2 {
        header.push_str(&format!(
            "_fragments_{}_paired_{}_diagnostic_{}_em_{:.3}_abundance_{:.6}",
            contig.fragment_support,
            contig.paired_fragment_support,
            contig.diagnostic_fragment_support,
            contig.em_fragment_support,
            contig.em_abundance
        ));
    }
    if mode == AssemblyMode::Its2 && !contig.label.is_empty() {
        header = format!(">{} {}", contig.label, &header[1..]);
    }
    header.push_str(&format!("_balance_{:.3}", contig.flank_balance));
    header
}

// 把候选按等级落盘，标题里带够指标，后续人工复核不抓瞎。
fn write_contigs(
    path: &Path,
    contigs: &[ContigRecord],
    mode: AssemblyMode,
    prefix: &str,
) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for contig in contigs {
        writeln!(writer, "{}", format_header(contig, mode, prefix))?;
        writeln!(writer, "{}", String::from_utf8_lossy(&contig.sequence))?;
    }
    writer.flush()
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Its2Placement {
    start: usize,
    end: usize,
    reverse: bool,
    read_start: isize,
    read_end: isize,
    mismatch: Option<usize>,
}

// Check the entire read/candidate overlap. Retain alternative placements so
// repeats cannot manufacture positional coverage or phase links.
fn its2_placements(mate: &[u8], candidate: &[u8], slice_len: usize) -> Vec<Its2Placement> {
    let mut placements = Vec::new();
    if mate.is_empty() || candidate.is_empty() || slice_len == 0 {
        return placements;
    }
    let minimum_overlap = slice_len.min(mate.len()).min(candidate.len());
    let reverse = reverse_complement(mate);
    for (strand, oriented) in [mate, reverse.as_slice()].into_iter().enumerate() {
        for offset in -(mate.len() as isize - minimum_overlap as isize)
            ..=(candidate.len() - minimum_overlap) as isize
        {
            let read_start = (-offset).max(0) as usize;
            let target_start = offset.max(0) as usize;
            let overlap = (mate.len() - read_start).min(candidate.len() - target_start);
            if oriented[read_start..read_start + overlap]
                .iter()
                .zip(&candidate[target_start..target_start + overlap])
                .all(|(&base, &target)| matches!(base, b'A' | b'C' | b'G' | b'T') && base == target)
            {
                placements.push(Its2Placement {
                    start: target_start,
                    end: target_start + overlap,
                    reverse: strand == 1,
                    read_start: offset,
                    read_end: offset + mate.len() as isize,
                    mismatch: None,
                });
            }
        }
    }

    // A short terminal match must not compete with a longer placement.
    if let Some(longest) = placements.iter().map(|p| p.end - p.start).max() {
        placements.retain(|p| p.end - p.start == longest);
    }
    placements
}

// Every overlap of at least 2*k-1 bases contains one complete read-grid seed.
// Seeds nominate offsets only; the entire overlap still has to match exactly.
struct Its2PlacementIndex<'a> {
    candidates: &'a [ContigRecord],
    seeds: HashMap<&'a [u8], Vec<(usize, usize)>>,
    minimum: usize,
    k: usize,
}

impl<'a> Its2PlacementIndex<'a> {
    fn new(candidates: &'a [ContigRecord], minimum: usize) -> Self {
        let k = minimum.div_ceil(2).clamp(1, 16);
        let mut seeds: HashMap<&[u8], Vec<(usize, usize)>> = HashMap::new();
        for (i, candidate) in candidates.iter().enumerate() {
            for (position, seed) in candidate.sequence.windows(k).enumerate() {
                if seed.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                    seeds.entry(seed).or_default().push((i, position));
                }
            }
        }
        Self {
            candidates,
            seeds,
            minimum,
            k,
        }
    }

    #[cfg(test)]
    fn placements(&self, mate: &[u8]) -> Vec<Vec<Its2Placement>> {
        self.place(mate, None)
    }

    fn place(&self, mate: &[u8], quality: Option<&[u8]>) -> Vec<Vec<Its2Placement>> {
        let quality = quality.filter(|q| q.len() == mate.len() && q.iter().any(|&v| v <= 12));
        let minimum = self.minimum;
        let mut result = vec![Vec::new(); self.candidates.len()];
        if self.minimum == 0 || mate.is_empty() {
            return result;
        }
        if mate.len() < minimum {
            return self
                .candidates
                .iter()
                .map(|c| its2_placements(mate, &c.sequence, self.minimum))
                .collect();
        }
        let reverse = reverse_complement(mate);
        for (strand, oriented) in [mate, reverse.as_slice()].into_iter().enumerate() {
            let mut offsets = vec![Vec::new(); self.candidates.len()];
            for (i, seed) in oriented.windows(self.k).enumerate().step_by(self.k) {
                if let Some(hits) = self.seeds.get(seed) {
                    for &(candidate, position) in hits {
                        offsets[candidate].push(position as isize - i as isize);
                    }
                }
            }
            for (i, offsets) in offsets.iter_mut().enumerate() {
                let candidate = &self.candidates[i].sequence;
                if candidate.len() < minimum {
                    continue;
                }
                offsets.sort_unstable();
                offsets.dedup();
                for &offset in offsets.iter() {
                    let a = (-offset).max(0) as usize;
                    let b = offset.max(0) as usize;
                    let overlap = (mate.len() - a).min(candidate.len() - b);
                    let mut mismatch = None;
                    let valid = overlap >= minimum
                        && oriented[a..a + overlap]
                            .iter()
                            .zip(&candidate[b..b + overlap])
                            .enumerate()
                            .all(|(j, (&base, &target))| {
                                if !matches!(base, b'A' | b'C' | b'G' | b'T') {
                                    return false;
                                }
                                if base == target {
                                    return true;
                                }
                                let read_position = if strand == 0 {
                                    a + j
                                } else {
                                    mate.len() - 1 - a - j
                                };
                                if mismatch.is_some()
                                    || quality.is_none_or(|q| q[read_position] > 12)
                                {
                                    return false;
                                }
                                mismatch = Some(b + j);
                                true
                            });
                    if valid && (mismatch.is_none() || overlap >= 75) {
                        result[i].push(Its2Placement {
                            start: b,
                            end: b + overlap,
                            reverse: strand == 1,
                            read_start: offset,
                            read_end: offset + mate.len() as isize,
                            mismatch,
                        });
                    }
                }
            }
        }
        for (i, placements) in result.iter_mut().enumerate() {
            if self.candidates[i].sequence.len() < minimum {
                *placements = its2_placements(mate, &self.candidates[i].sequence, self.minimum);
            } else if let Some(longest) = placements.iter().map(|p| p.end - p.start).max() {
                placements.retain(|p| p.end - p.start == longest);
            }
        }
        result
    }
}

#[cfg(test)]
fn mate_matches_candidate(mate: &[u8], candidate: &[u8], slice_len: usize) -> bool {
    !its2_placements(mate, candidate, slice_len).is_empty()
}

fn inward_pair(a: Its2Placement, b: Its2Placement) -> bool {
    if a.reverse == b.reverse {
        return false;
    }
    let (forward, reverse) = if a.reverse { (b, a) } else { (a, b) };
    forward.start <= reverse.start && forward.end <= reverse.end
}

// For equal-length alternatives the distinguishing sites are explicit. For
// length-changing alternatives conservatively include the differing block and
// its flanks, until a gapped alignment can narrow the phase requirement.
fn its2_variant_sites(index: usize, candidates: &[ContigRecord]) -> Vec<usize> {
    let sequence = &candidates[index].sequence;
    let mut sites = std::collections::BTreeSet::new();
    for (other_index, other) in candidates.iter().enumerate() {
        if other_index == index || other.sequence.is_empty() || sequence.is_empty() {
            continue;
        }
        if sequence.len() == other.sequence.len() {
            sites.extend(
                sequence
                    .iter()
                    .zip(&other.sequence)
                    .enumerate()
                    .filter_map(|(pos, (a, b))| (a != b).then_some(pos)),
            );
        } else {
            let prefix = sequence
                .iter()
                .zip(&other.sequence)
                .take_while(|(a, b)| a == b)
                .count();
            let suffix = sequence[prefix..]
                .iter()
                .rev()
                .zip(other.sequence[prefix..].iter().rev())
                .take_while(|(a, b)| a == b)
                .count();
            sites.extend(
                prefix.saturating_sub(1)..(sequence.len() - suffix + 1).min(sequence.len()),
            );
        }
    }
    sites.into_iter().collect()
}

// Trim to the aligned marker span, including a partial span. The retained sequence is
// always a slice of the assembly, including sample insertions; no bases are filled.
fn trim_its2_reference_span(sequence: &mut Vec<u8>, index: &uce_filter_core::index::UceIndex) {
    if sequence.len() > u16::MAX as usize {
        return;
    }
    let Some(alignment) = uce_filter_core::alignment::align_marker(index, sequence, 0, 32) else {
        return;
    };
    if alignment.identity() < 0.8 {
        return;
    }
    let core = sequence[alignment.query_start as usize..alignment.query_end as usize].to_vec();
    *sequence = if alignment.strand == 2 {
        reverse_complement(&core)
    } else {
        core
    };
}

fn append_equivalence_member(members: &mut String, label: &str) {
    if label.is_empty() || members.split(";").any(|member| member == label) {
        return;
    }
    if !members.is_empty() {
        members.push(';');
    }
    members.push_str(label);
}

fn prefer_its2_candidate(
    existing: &ContigRecord,
    candidate: &ContigRecord,
    mode: AssemblyMode,
) -> bool {
    if existing.label.is_empty() != candidate.label.is_empty() {
        return !candidate.label.is_empty();
    }
    compare_contigs(existing, candidate, mode).is_lt()
}

// ITS2 候选分组、配对支持和 EM 丰度都在这里补齐。
#[cfg(test)]
fn annotate_its2_candidates(
    candidates: &mut [ContigRecord],
    fragments: &[Vec<Vec<u8>>],
    slice_len: usize,
) {
    annotate_its2_candidates_with_quality(candidates, fragments, None, slice_len);
}

fn annotate_its2_candidates_with_quality(
    candidates: &mut [ContigRecord],
    fragments: &[Vec<Vec<u8>>],
    qualities: Option<&LinkedReads>,
    slice_len: usize,
) {
    annotate_its2_evidence_with_quality(candidates, fragments, qualities, slice_len, false);

    for candidate in candidates
        .iter_mut()
        .filter(|c| !c.its2_full_reference_span)
    {
        candidate.accepted = false;
        if candidate.supported_bases == candidate.sequence.len() {
            candidate.rejection_reason = "partial_reference_span".into();
        }
    }
    // Unsupported extensions remain in the diagnostic output, but must not
    // prevent a fully covered marker from obtaining diagnostic support.
    let covered: Vec<usize> = candidates
        .iter()
        .enumerate()
        .filter(|(_, c)| {
            c.its2_full_reference_span
                && c.fragment_support >= 2
                && c.supported_bases == c.sequence.len()
        })
        .map(|(i, _)| i)
        .collect();
    let mut complete = Vec::new();
    for &index in &covered {
        let sequence = &candidates[index].sequence;
        let contained = covered.iter().any(|&other| {
            candidates[other].sequence.len() > sequence.len()
                && !sequence.is_empty()
                && candidates[other]
                    .sequence
                    .windows(sequence.len())
                    .any(|s| s == sequence)
        });
        if contained {
            candidates[index].accepted = false;
            candidates[index].rejection_reason = "contained_model".into();
        } else {
            complete.push(index);
        }
    }
    if !complete.is_empty() {
        let mut viable: Vec<_> = complete.iter().map(|&i| candidates[i].clone()).collect();
        annotate_its2_evidence_with_quality(&mut viable, fragments, qualities, slice_len, true);
        for (index, candidate) in complete.into_iter().zip(viable) {
            candidates[index] = candidate;
        }
    }
}

type Its2FragmentPosition = Vec<(isize, isize, bool)>;

// Compare only identical alignment geometry: clipping and indels need a gapped
// likelihood model and must not acquire an accidental length advantage here.
fn its2_likelihood_weights(
    indices: &[usize],
    candidates: &[ContigRecord],
    placements: &[Vec<Vec<Its2Placement>>],
    fragment: &[Vec<u8>],
    quality: Option<&Vec<Vec<u8>>>,
) -> Vec<f64> {
    let neutral = || vec![1.0; indices.len()];
    if indices.len() <= 1 {
        return neutral();
    }
    let Some(quality) = quality else {
        return neutral();
    };
    let Some(&first) = indices.first() else {
        return Vec::new();
    };
    let geometry = |p: &Its2Placement| (p.start, p.end, p.reverse, p.read_start, p.read_end);
    for &i in indices {
        if candidates[i].sequence.len() != candidates[first].sequence.len()
            || placements[i].iter().zip(&placements[first]).any(|(a, b)| {
                a.len() > 1 || b.len() > 1 || a.first().map(geometry) != b.first().map(geometry)
            })
        {
            return neutral();
        }
    }
    // A read pair is one fragment. Overlapping mates do not provide two
    // independent observations; conflicting bases contribute no likelihood.
    let mut observations = HashMap::<usize, (u8, u8, bool)>::new();
    for (mate, ps) in placements[first].iter().enumerate() {
        let Some(p) = ps.first() else { continue };
        let Some(q) = quality
            .get(mate)
            .filter(|q| q.len() == fragment[mate].len())
        else {
            return neutral();
        };
        for site in p.start..p.end {
            if indices
                .iter()
                .all(|&i| candidates[i].sequence[site] == candidates[first].sequence[site])
            {
                continue;
            }
            let offset = (site as isize - p.read_start) as usize;
            let offset = if p.reverse {
                fragment[mate].len() - 1 - offset
            } else {
                offset
            };
            let base = fragment[mate][offset];
            let base = if p.reverse {
                match base {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'N',
                }
            } else {
                base
            };
            let entry = observations.entry(site).or_insert((base, q[offset], false));
            entry.2 |= entry.0 != base;
            entry.1 = entry.1.max(q[offset]);
        }
    }
    let mut scores = vec![0.0_f64; indices.len()];
    // Ordered accumulation keeps output reproducible across randomized hash seeds.
    let mut observations: Vec<_> = observations.into_iter().collect();
    observations.sort_unstable_by_key(|x| x.0);
    for (site, (base, q, conflict)) in observations {
        if conflict || base == b'N' {
            continue;
        }
        let error = 10_f64.powf(-(q.min(60) as f64) / 10.0).min(0.75);
        let penalty = ((error / 3.0) / (1.0 - error)).ln();
        for (score, &i) in scores.iter_mut().zip(indices) {
            if candidates[i].sequence[site] != base {
                *score += penalty;
            }
        }
    }
    let best = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    scores.into_iter().map(|s| (s - best).exp()).collect()
}

// Each edge requires two different mapped fragment positions. Connectivity,
// rather than an obligatory adjacent-site chain, establishes the phase block.
fn its2_phase_component_count(
    sites: usize,
    links: &HashMap<(usize, usize), HashSet<Its2FragmentPosition>>,
) -> usize {
    let mut roots: Vec<_> = (0..sites).collect();
    for (&(a, b), evidence) in links {
        if evidence.len() < 2 {
            continue;
        }
        let old = roots[b];
        let new = roots[a];
        for root in &mut roots {
            if *root == old {
                *root = new;
            }
        }
    }
    roots.sort_unstable();
    roots.dedup();
    roots.len()
}

fn its2_phase_base(
    position: usize,
    placements: &[Vec<Its2Placement>],
    fragment: &[Vec<u8>],
    quality: Option<&Vec<Vec<u8>>>,
) -> bool {
    placements.iter().enumerate().any(|(mate, p)| {
        if p.len() != 1 {
            return false;
        }
        let p = p[0];
        if position < p.start || position >= p.end || p.mismatch == Some(position) {
            return false;
        }
        let read_position = (position as isize - p.read_start) as usize;
        let read_position = if p.reverse {
            fragment[mate].len() - 1 - read_position
        } else {
            read_position
        };
        quality
            .and_then(|q| q.get(mate))
            .filter(|q| !q.is_empty())
            .is_none_or(|q| q.get(read_position).is_some_and(|&q| q >= 20))
    })
}

#[cfg(test)]
fn annotate_its2_evidence(
    candidates: &mut [ContigRecord],
    fragments: &[Vec<Vec<u8>>],
    slice_len: usize,
    check_phase: bool,
) {
    annotate_its2_evidence_with_quality(candidates, fragments, None, slice_len, check_phase);
}

fn annotate_its2_evidence_with_quality(
    candidates: &mut [ContigRecord],
    fragments: &[Vec<Vec<u8>>],
    qualities: Option<&LinkedReads>,
    slice_len: usize,
    check_phase: bool,
) {
    let mut compatibilities = Vec::new();
    let mut diagnostic_positions: Vec<HashSet<Its2FragmentPosition>> =
        (0..candidates.len()).map(|_| HashSet::new()).collect();
    let mut coverage: Vec<Vec<u64>> = candidates
        .iter()
        .map(|c| vec![0; c.sequence.len()])
        .collect();
    let sites: Vec<_> = (0..candidates.len())
        .map(|i| {
            if check_phase {
                its2_variant_sites(i, candidates)
            } else {
                Vec::new()
            }
        })
        .collect();
    let mut links: Vec<HashMap<(usize, usize), HashSet<Its2FragmentPosition>>> =
        (0..sites.len()).map(|_| HashMap::new()).collect();
    let mut site_positions: Vec<Vec<HashSet<Its2FragmentPosition>>> = sites
        .iter()
        .map(|s| vec![HashSet::new(); s.len()])
        .collect();
    for candidate in candidates.iter_mut() {
        candidate.phase_components = None;
        candidate.min_variant_site_support = None;
        candidate.read_count = 0;
        candidate.unique_read_count = 0;
        candidate.multi_mapping_read_count = 0;
        candidate.fragment_support = 0;
        candidate.paired_fragment_support = 0;
        candidate.diagnostic_fragment_support = 0;
        candidate.quality_rescued_fragment_support = 0;
        candidate.em_fragment_support = f64::NAN;
        candidate.em_abundance = f64::NAN;
        candidate.accepted = false;
        candidate.rejection_reason = "no_fragment_support".into();
    }
    // Repeated fragments contribute their original counts, but need alignment
    // only once. Borrow the existing sequences rather than copying read data.
    let mut unique_fragments = HashMap::new();
    for (i, fragment) in fragments.iter().enumerate() {
        let quality = qualities.map(|q| &q[i]);
        *unique_fragments.entry((fragment, quality)).or_insert(0_u64) += 1;
    }
    let mut unique_fragments: Vec<_> = unique_fragments.into_iter().collect();
    unique_fragments.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let indexed_candidates = candidates.to_vec();
    let placement_index = Its2PlacementIndex::new(&indexed_candidates, slice_len);
    for ((fragment, quality), multiplicity) in unique_fragments {
        let mut indices = Vec::new();
        let by_mate: Vec<_> = fragment
            .iter()
            .enumerate()
            .map(|(j, mate)| {
                // Low-quality alternatives remain shared compatibility, not corrections
                // or unique evidence. Verify exact and tolerated offsets in one pass.
                placement_index.place(mate, quality.and_then(|q| q.get(j)).map(Vec::as_slice))
            })
            .collect();
        let all_placements: Vec<Vec<Vec<Its2Placement>>> = (0..candidates.len())
            .map(|i| {
                by_mate
                    .iter()
                    .map(|placements| placements[i].clone())
                    .collect()
            })
            .collect();
        let paired_alternative = all_placements.iter().enumerate().any(|(i, p)| {
            candidates[i].its2_full_reference_span
                && p.len() == 2
                && p[0]
                    .iter()
                    .any(|&a| p[1].iter().any(|&b| inward_pair(a, b)))
        });
        for (index, candidate) in candidates.iter_mut().enumerate() {
            let placements = &all_placements[index];
            if paired_alternative && placements.iter().any(Vec::is_empty) {
                continue;
            }
            if placements.iter().all(Vec::is_empty) {
                continue;
            }
            let both = placements.len() == 2 && placements.iter().all(|p| !p.is_empty());
            if both
                && !placements[0]
                    .iter()
                    .any(|&a| placements[1].iter().any(|&b| inward_pair(a, b)))
            {
                continue;
            }
            indices.push(index);
            candidate.fragment_support += multiplicity;
            if placements.iter().flatten().any(|p| p.mismatch.is_some()) {
                candidate.quality_rescued_fragment_support += multiplicity;
            }
            candidate.read_count +=
                placements.iter().filter(|p| !p.is_empty()).count() as u64 * multiplicity;
            candidate.unique_read_count +=
                placements.iter().filter(|p| p.len() == 1).count() as u64 * multiplicity;
            candidate.multi_mapping_read_count +=
                placements.iter().filter(|p| p.len() > 1).count() as u64 * multiplicity;
            let unique: Vec<_> = placements
                .iter()
                .filter(|p| p.len() == 1)
                .map(|p| p[0])
                .collect();
            if both && unique.len() == 2 {
                candidate.paired_fragment_support += multiplicity;
            }
            let covers = |position| {
                unique.iter().any(|p| {
                    p.start <= position && position < p.end && p.mismatch != Some(position)
                })
            };
            for (position, depth) in coverage[index].iter_mut().enumerate() {
                if covers(position) {
                    *depth += multiplicity;
                }
            }
            let observed: Vec<_> = sites[index]
                .iter()
                .enumerate()
                .filter(|(_, site)| its2_phase_base(**site, placements, fragment, quality))
                .map(|(i, _)| i)
                .collect();
            if !observed.is_empty() {
                let mut signature: Vec<_> = unique
                    .iter()
                    .map(|p| (p.read_start, p.read_end, p.reverse))
                    .collect();
                signature.sort_unstable();
                for (offset, &a) in observed.iter().enumerate() {
                    site_positions[index][a].insert(signature.clone());
                    for &b in &observed[offset + 1..] {
                        let evidence = links[index].entry((a, b)).or_default();
                        if evidence.len() < 2 {
                            evidence.insert(signature.clone());
                        }
                    }
                }
            }
        }
        if indices.len() == 1 {
            let index = indices[0];
            let placements = &all_placements[index];
            let high_quality = sites[index].iter().all(|&site| {
                let observed = placements
                    .iter()
                    .flatten()
                    .any(|p| p.start <= site && site < p.end);
                !observed || its2_phase_base(site, placements, fragment, quality)
            });
            if high_quality && placements.iter().all(|p| p.len() <= 1) {
                let mut positions: Vec<_> = placements
                    .iter()
                    .flatten()
                    .map(|p| (p.read_start, p.read_end, p.reverse))
                    .collect();
                positions.sort_unstable();
                diagnostic_positions[index].insert(positions);
            }
        }
        if check_phase && !indices.is_empty() {
            let weights =
                its2_likelihood_weights(&indices, candidates, &all_placements, fragment, quality);
            compatibilities.push((indices, weights, multiplicity));
        }
    }
    for (index, candidate) in candidates.iter_mut().enumerate() {
        if check_phase {
            candidate.phase_components = Some(its2_phase_component_count(
                sites[index].len(),
                &links[index],
            ));
            candidate.min_variant_site_support = site_positions[index]
                .iter()
                .map(|positions| positions.len())
                .min();
        }
        candidate.diagnostic_fragment_support = diagnostic_positions[index].len() as u64;
        candidate.supported_bases = coverage[index].iter().filter(|&&depth| depth >= 2).count();
        candidate.support_breadth = if candidate.sequence.is_empty() {
            0.0
        } else {
            candidate.supported_bases as f64 / candidate.sequence.len() as f64
        };
        let mut gap = 0;
        candidate.max_support_gap = 0;
        for depth in &coverage[index] {
            gap = if *depth < 2 { gap + 1 } else { 0 };
            candidate.max_support_gap = candidate.max_support_gap.max(gap);
        }
        let start = coverage[index].iter().position(|&depth| depth >= 2);
        let end = coverage[index].iter().rposition(|&depth| depth >= 2);
        candidate.supported_span = start.zip(end).map_or(0, |(a, b)| b - a + 1);
        let length = candidate.sequence.len().max(1) as f64;
        candidate.support_fraction = candidate.supported_span as f64 / length;
        candidate.read_density = candidate.read_count as f64 / length;
        let middle = coverage[index].len() / 2;
        let left: u64 = coverage[index][..middle].iter().sum();
        let right: u64 = coverage[index][middle..].iter().sum();
        candidate.flank_balance = left.min(right) as f64 / left.max(right).max(1) as f64;
        candidate.rejection_reason = if candidate.fragment_support < 2 {
            "low_fragment_support"
        } else if candidate.supported_bases != candidate.sequence.len() {
            "incomplete_read_coverage"
        } else if candidate.phase_components.is_some_and(|count| count > 1) {
            "unresolved_phase_links"
        } else if candidate.diagnostic_fragment_support < 2 {
            "unresolved_competing_paths"
        } else {
            ""
        }
        .into();
        candidate.accepted = candidate.rejection_reason.is_empty();
    }
    if !check_phase || compatibilities.is_empty() {
        return;
    }
    let mut abundance: Vec<f64> = candidates
        .iter()
        .map(|candidate| candidate.fragment_support.max(1) as f64)
        .collect();
    let normalizer: f64 = abundance.iter().sum();
    for value in &mut abundance {
        *value /= normalizer;
    }
    let mut expected = vec![0.0; candidates.len()];
    for _ in 0..64 {
        expected.fill(0.0);
        for (compatible, weights, multiplicity) in &compatibilities {
            let denominator: f64 = compatible
                .iter()
                .zip(weights)
                .map(|(index, weight)| abundance[*index] * weight)
                .sum();
            if denominator > 0.0 {
                for (index, weight) in compatible.iter().zip(weights) {
                    expected[*index] +=
                        *multiplicity as f64 * abundance[*index] * weight / denominator;
                }
            }
        }
        let total: f64 = expected.iter().sum();
        if total > 0.0 {
            for (value, count) in abundance.iter_mut().zip(expected.iter()) {
                *value = *count / total;
            }
        }
    }
    for (candidate, (count, value)) in candidates
        .iter_mut()
        .zip(expected.into_iter().zip(abundance))
    {
        candidate.em_fragment_support = if candidate.accepted { count } else { f64::NAN };
        candidate.em_abundance = if candidate.accepted { value } else { f64::NAN };
    }
}

fn its2_status(contig: &ContigRecord) -> &'static str {
    if contig.accepted {
        "READ_SUPPORTED"
    } else if matches!(
        contig.rejection_reason.as_str(),
        "unresolved_phase_links" | "unresolved_competing_paths" | "contained_model"
    ) {
        "UNPHASED"
    } else if matches!(
        contig.rejection_reason.as_str(),
        "incomplete_read_coverage" | "partial_reference_span"
    ) {
        "PARTIAL"
    } else {
        "LOW_SUPPORT"
    }
}

// 把 ITS2 的支持证据另存表，主 FASTA 不塞得乱七八糟。
fn write_its2_support(path: &Path, contigs: &[ContigRecord]) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "variant\tequivalence_members\tlength\tfragment_support\tpaired_fragment_support\tdiagnostic_fragment_support\tem_fragment_support\tem_abundance\tstatus\tsupported_bases\tsupport_breadth\tmax_support_gap\trejection_reason\tquality_rescued_fragment_support\tphase_components\tmin_variant_site_support")?;
    for (index, contig) in contigs.iter().enumerate() {
        let label = if contig.label.is_empty() {
            format!("assembly_{}", index + 1)
        } else {
            contig.label.clone()
        };
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}",
            label,
            if contig.equivalence_members.is_empty() {
                label.clone()
            } else {
                contig.equivalence_members.clone()
            },
            contig.sequence.len(),
            contig.fragment_support,
            contig.paired_fragment_support,
            contig.diagnostic_fragment_support,
            if contig.em_fragment_support.is_finite() {
                format!("{:.6}", contig.em_fragment_support)
            } else {
                "NA".into()
            },
            if contig.em_abundance.is_finite() {
                format!("{:.8}", contig.em_abundance)
            } else {
                "NA".into()
            },
            its2_status(contig),
            contig.supported_bases,
            contig.support_breadth,
            contig.max_support_gap,
            contig.rejection_reason,
            contig.quality_rescued_fragment_support,
            contig
                .phase_components
                .map_or_else(|| "NA".into(), |v| v.to_string()),
            contig
                .min_variant_site_support
                .map_or_else(|| "NA".into(), |v| v.to_string())
        )?;
    }
    writer.flush()
}

pub fn run_manifest(args: &Args, tasks: &[LocusTask]) -> io::Result<String> {
    let mut hasher = DefaultHasher::new();
    if args.assembly_mode == AssemblyMode::Its2 {
        "its2-evidence-v16".hash(&mut hasher);
    }
    args.kmer_size.hash(&mut hasher);
    args.kmer_min.hash(&mut hasher);
    args.kmer_max.hash(&mut hasher);
    args.error_limit.hash(&mut hasher);
    args.iteration.hash(&mut hasher);
    args.min_coverage.to_bits().hash(&mut hasher);
    args.soft_boundary.hash(&mut hasher);
    args.assembly_mode.hash(&mut hasher);
    args.side_candidates.hash(&mut hasher);
    args.path_strategy.hash(&mut hasher);
    args.backbone_lookahead.hash(&mut hasher);
    args.reverse_reuse_reference_scale
        .to_bits()
        .hash(&mut hasher);
    args.max_contig_length.hash(&mut hasher);
    args.min_read_density.to_bits().hash(&mut hasher);
    args.density_check_min_length.hash(&mut hasher);
    args.max_depth_cv.to_bits().hash(&mut hasher);
    args.max_depth_ratio.to_bits().hash(&mut hasher);
    for task in tasks {
        task.key.hash(&mut hasher);
        let reference = fs::metadata(&task.reference_path)?;
        reference.len().hash(&mut hasher);
        reference
            .modified()
            .ok()
            .and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
            .map(|time| time.as_nanos())
            .unwrap_or(0)
            .hash(&mut hasher);
        match find_filtered(&args.output, &task.key) {
            Some((path, fasta)) => {
                path.hash(&mut hasher);
                fasta.hash(&mut hasher);
                let filtered = fs::metadata(path)?;
                filtered.len().hash(&mut hasher);
                filtered
                    .modified()
                    .ok()
                    .and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
                    .map(|time| time.as_nanos())
                    .unwrap_or(0)
                    .hash(&mut hasher);
            }
            None => false.hash(&mut hasher),
        }
    }
    Ok(format!(
        "GM2_ASSEMBLER_RUN_MANIFEST_V1\n{:016x}\n",
        hasher.finish()
    ))
}

pub fn process_locus(
    args: &Args,
    task: &LocusTask,
    completed: &HashSet<String>,
    log_lock: &Arc<Mutex<()>>,
    profile: Option<&ProfileStats>,
) -> LocusResult {
    // 单个位点出岔子就收拾它自己的输出，别把整批组装拖下水。
    let result = match process_locus_inner(args, task, completed, log_lock, profile) {
        Ok(result) => result,
        Err(error) => {
            clean_locus_outputs(args, &task.key);
            log_line(
                &args.output,
                log_lock,
                &format!("error assembling {}: {error}", task.key),
            );
            LocusResult::failure(&task.key, "error")
        }
    };
    if let Some(profile) = profile {
        profile.record_locus(&result);
    }
    result
}

fn process_locus_inner(
    args: &Args,
    task: &LocusTask,
    completed: &HashSet<String>,
    log_lock: &Arc<Mutex<()>>,
    profile: Option<&ProfileStats>,
) -> io::Result<LocusResult> {
    // 一个 locus 从过滤 reads 到结果文件，全流程在这旮沓串起来。
    let best_path = args
        .output
        .join("results")
        .join(format!("{}.fasta", task.key));
    let all_path = args
        .output
        .join("contigs_all")
        .join(format!("{}.fasta", task.key));
    let low_path = args
        .output
        .join("contigs_all_low")
        .join(format!("{}.fasta", task.key));

    if completed.contains(&task.key)
        && best_path.is_file()
        && fs::metadata(&best_path).is_ok_and(|metadata| metadata.len() > 0)
    {
        return Ok(LocusResult {
            key: task.key.clone(),
            status: "skipped".to_string(),
            skipped: true,
            ..LocusResult::default()
        });
    }

    let Some((filtered_path, fasta)) = find_filtered(&args.output, &task.key) else {
        clean_locus_outputs(args, &task.key);
        return Ok(LocusResult::failure(&task.key, "no filtered file"));
    };
    let paired_fragments =
        args.assembly_mode == AssemblyMode::Uce && has_interleaved_pairs(&filtered_path, fasta)?;
    let count_chunk_size = if paired_fragments {
        (args.read_chunk_size.max(2) / 2) * 2
    } else {
        args.read_chunk_size
    };
    let reads_started = Instant::now();
    let (linked_fragments, linked_qualities) = if args.assembly_mode == AssemblyMode::Its2 {
        read_linked_reads_with_quality(&filtered_path, fasta)?
    } else {
        (Vec::new(), Vec::new())
    };
    let minimum = minimum_sequence_length(&filtered_path, fasta, args.read_chunk_size)?;
    let slice_len = minimum.map_or(0, |length| ((length as f64 * 0.9) as usize).max(1));
    if slice_len == 0 {
        clean_locus_outputs(args, &task.key);
        log_line(
            &args.output,
            log_lock,
            &format!("No reads were obtained for gene {}", task.key),
        );
        return Ok(LocusResult::failure(&task.key, "no reads"));
    }

    let mut reads = HashMap::new();
    let fixed_k = args.kmer_size != 0;
    if !fixed_k {
        for_each_sequence_chunk(&filtered_path, fasta, args.read_chunk_size, |chunk| {
            add_read_slices(&mut reads, chunk, slice_len);
            Ok(())
        })?;
    }
    if let Some(profile) = profile {
        profile.record_reads(reads_started);
    }
    if !fixed_k && reads.is_empty() {
        clean_locus_outputs(args, &task.key);
        log_line(
            &args.output,
            log_lock,
            &format!("No reads were obtained for gene {}", task.key),
        );
        return Ok(LocusResult::failure(&task.key, "no reads"));
    }

    let kmer_started = Instant::now();
    let reference_records = read_fasta(&task.reference_path)?;
    let reference_sequences: Vec<Vec<u8>> = reference_records
        .iter()
        .map(|(_, sequence)| sequence.clone())
        .collect();
    let current_k = if fixed_k {
        args.kmer_size
    } else {
        // 用户没指定 k 就按 reads 和参考自动挑，别死抱默认值。
        calculate_auto_k(
            &reference_sequences,
            &reads,
            slice_len,
            args.kmer_min,
            args.kmer_max,
            args.error_limit,
        )
    };
    if current_k == 0 || current_k > 63 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("assembly k-mer must be between 1 and 63, got {current_k}"),
        ));
    }

    let soft_boundary = if args.soft_boundary < 0 {
        slice_len / 2
    } else {
        args.soft_boundary as usize
    };
    log_line(
        &args.output,
        log_lock,
        &format!("Use k={current_k} for assembling gene {}.", task.key),
    );
    log_line(
        &args.output,
        log_lock,
        &format!("Assembling {} {} / {}", task.key, task.ordinal, task.total),
    );

    let reference = load_or_build_reference_kmers(
        &task.reference_path,
        &reference_records,
        current_k,
        args.reference_cache_dir.as_deref(),
    )?;
    let mut count_runs = SortedKmerCounts::default();
    if fixed_k {
        if let Some(profile) = profile {
            profile.record_kmers(kmer_started);
        }
        for_each_sequence_chunk(&filtered_path, fasta, count_chunk_size, |chunk| {
            let reads_chunk_started = Instant::now();
            add_read_slices(&mut reads, chunk, slice_len);
            if let Some(profile) = profile {
                profile.record_reads_duration(reads_chunk_started.elapsed());
            }
            let kmers_chunk_started = Instant::now();
            count_runs.push(count_assemble_chunk_parallel(
                chunk,
                current_k,
                args.kmer_count_threads,
                paired_fragments,
                paired_fragments.then_some(&reference),
            ));
            if let Some(profile) = profile {
                profile.record_kmers_duration(kmers_chunk_started.elapsed());
            }
            Ok(())
        })?;
        if reads.is_empty() {
            clean_locus_outputs(args, &task.key);
            log_line(
                &args.output,
                log_lock,
                &format!("No reads were obtained for gene {}", task.key),
            );
            return Ok(LocusResult::failure(&task.key, "no reads"));
        }
    } else {
        for_each_sequence_chunk(&filtered_path, fasta, count_chunk_size, |chunk| {
            count_runs.push(count_assemble_chunk_parallel(
                chunk,
                current_k,
                args.kmer_count_threads,
                paired_fragments,
                paired_fragments.then_some(&reference),
            ));
            Ok(())
        })?;
        if let Some(profile) = profile {
            profile.record_kmers(kmer_started);
        }
    }
    let graph_started = Instant::now();
    let min_fragment_support = u32::from(paired_fragments) * 2;
    let counts = count_runs.into_counts();
    // Both candidate graphs share one k-mer count. The core graph supplies a
    // safe no-PE path; the assisted graph only needs to beat it on the normal
    // UCE evidence ranking to become the reported contig.
    let mut core_graph = build_graph_from_counts(counts.clone(), &reference, args.error_limit, 0);
    filter_and_weight_graph(&mut core_graph, args.error_limit, task.reference_count, 0);
    let mut graph =
        build_graph_from_counts(counts, &reference, args.error_limit, min_fragment_support);
    filter_and_weight_graph(
        &mut graph,
        args.error_limit,
        task.reference_count,
        min_fragment_support,
    );
    let mut branch_support = branch_edges(&graph, current_k);
    if paired_fragments && !branch_support.is_empty() {
        for_each_sequence_chunk(&filtered_path, fasta, count_chunk_size, |chunk| {
            add_pe_branch_support(chunk, current_k, &reference, &mut branch_support);
            Ok(())
        })?;
    }
    let (gfa, dot) = match args.graph_format {
        GraphFormat::None => (false, false),
        GraphFormat::Gfa => (true, false),
        GraphFormat::Dot => (false, true),
        GraphFormat::Both => (true, true),
    };
    if gfa || dot {
        write_graphs(
            &args.output.join("assembly_graphs"),
            &task.key,
            &graph,
            current_k,
            gfa,
            dot,
        )?;
    }
    if let Some(profile) = profile {
        profile.record_graph(graph_started);
    }
    if graph.len() < 3 {
        clean_locus_outputs(args, &task.key);
        log_line(
            &args.output,
            log_lock,
            "Could not get enough reads from filter.",
        );
        return Ok(LocusResult::failure(
            &task.key,
            "insufficient genomic kmers",
        ));
    }

    let path_started = Instant::now();
    let mut seeds: Vec<(u128, i64, i32, i64)> = graph
        .iter()
        .filter(|(_, value)| value.position > 1 && value.position < 1000 && !value.is_reverse)
        .map(|(kmer, value)| (*kmer, value.depth, value.position, value.reference_weight))
        .collect();
    seeds.sort_by(|left, right| {
        right
            .3
            .cmp(&left.3)
            .then_with(|| right.1.cmp(&left.1))
            .then_with(|| left.0.cmp(&right.0))
    });
    if seeds.is_empty() {
        clean_locus_outputs(args, &task.key);
        log_line(&args.output, log_lock, "Could not get enough seeds.");
        return Ok(LocusResult::failure(&task.key, "no seed"));
    }

    // Prevent ITS2 paths from expanding through unrelated rDNA flanks before
    // combining sides. Ends are observed graph k-mers, not reference-filled bases.
    let terminal_kmers: Option<HashSet<u128>> =
        (args.assembly_mode == AssemblyMode::Its2).then(|| {
            reference_sequences
                .iter()
                .filter(|s| s.len() >= current_k)
                .flat_map(|sequence| {
                    [
                        encode_kmer(&sequence[sequence.len() - current_k..]),
                        encode_kmer(&sequence[..current_k])
                            .map(|kmer| reverse_complement_kmer(kmer, current_k)),
                    ]
                    .into_iter()
                    .flatten()
                })
                .collect()
        });
    let path_evidence = terminal_kmers.map(|mut ends| {
        // A shorter panel entry may end inside another complete ITS2. Such a
        // shared internal k-mer must not truncate every longer sample path.
        for sequence in &reference_sequences {
            if sequence.len() > current_k + 1 {
                for window in sequence[1..sequence.len() - 1].windows(current_k) {
                    if let Some(kmer) = encode_kmer(window) {
                        ends.remove(&kmer);
                        ends.remove(&reverse_complement_kmer(kmer, current_k));
                    }
                }
            }
        }
        let longest = reference_sequences.iter().map(Vec::len).max().unwrap_or(0);
        Its2PathEvidence::new(&reads, current_k, ends)
            .with_max_extension(longest.saturating_add(2 * current_k))
    });
    let original_seed_count = seeds.len();
    let seed_set: HashSet<u128> = seeds.iter().map(|seed| seed.0).collect();
    let mut accepted = Vec::new();
    let mut rejected = Vec::new();
    let mut core_accepted = Vec::new();

    // The legacy single-contig stopping rule can discard a smaller disconnected
    // ITS2 component once the dominant component accounts for half the seeds.
    while !seeds.is_empty()
        && (args.assembly_mode == AssemblyMode::Its2 || seeds.len() > original_seed_count / 2)
    {
        let seed = seeds[0].0;
        let (mut candidates, used_kmers, _) = assemble_seed(
            args,
            &reads,
            slice_len,
            &graph,
            seed,
            current_k,
            soft_boundary,
            paired_fragments.then_some(&branch_support),
            path_evidence.as_ref(),
        );
        let mut core_start = None;
        let core_used_kmers = if paired_fragments {
            let (mut core_candidates, core_used_kmers, _) = assemble_seed(
                args,
                &reads,
                slice_len,
                &core_graph,
                seed,
                current_k,
                soft_boundary,
                None,
                None,
            );
            core_start = Some(candidates.len());
            candidates.append(&mut core_candidates);
            Some(core_used_kmers)
        } else {
            None
        };
        let seed_coverage = core_used_kmers.as_ref().unwrap_or(&used_kmers);
        let used_seed_count = seed_set.intersection(seed_coverage).count();
        for candidate in &mut candidates {
            candidate.seed_count = used_seed_count;
        }

        seeds.retain(|candidate| {
            !seed_coverage.contains(&candidate.0)
                && !seed_coverage.contains(&reverse_complement_kmer(candidate.0, current_k))
        });
        for (candidate_index, candidate) in candidates.into_iter().enumerate() {
            let is_core_candidate = core_start.is_some_and(|start| candidate_index >= start);
            if args.assembly_mode == AssemblyMode::Uce {
                if candidate.accepted {
                    if is_core_candidate {
                        core_accepted.push(candidate.clone());
                    }
                    accepted.push(candidate);
                } else {
                    rejected.push(candidate);
                }
            } else if args.assembly_mode == AssemblyMode::Its2
                || candidate.read_count as usize * slice_len > candidate.sequence.len()
            {
                accepted.push(candidate);
            } else {
                rejected.push(candidate);
            }
        }
    }

    if let Some(profile) = profile {
        profile.record_path(path_started);
    }
    let qc_started = Instant::now();
    if args.assembly_mode == AssemblyMode::Its2 {
        // References guide recruitment and annotation; they are not sample assemblies.
        let marker_index = uce_filter_core::index::UceIndex::build(&task.reference_path, current_k)
            .map_err(io::Error::other)?;
        let normalize_started = Instant::now();
        let mut unique: HashMap<Vec<u8>, ContigRecord> = HashMap::new();
        let mut normalized: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
        for mut candidate in accepted.drain(..).chain(rejected.drain(..)) {
            let original = std::mem::take(&mut candidate.sequence);
            candidate.sequence = normalized
                .entry(original)
                .or_insert_with_key(|sequence| {
                    let mut trimmed = sequence.clone();
                    trim_its2_reference_span(&mut trimmed, &marker_index);
                    trimmed
                })
                .clone();
            if candidate.equivalence_members.is_empty() {
                candidate.equivalence_members = candidate.label.clone();
            }
            match unique.get_mut(&candidate.sequence) {
                Some(existing) => {
                    let mut members = existing.equivalence_members.clone();
                    append_equivalence_member(&mut members, &candidate.label);
                    if prefer_its2_candidate(existing, &candidate, args.assembly_mode) {
                        candidate.equivalence_members = members;
                        *existing = candidate;
                    } else {
                        existing.equivalence_members = members;
                    }
                }
                None => {
                    unique.insert(candidate.sequence.clone(), candidate);
                }
            }
        }
        // Equal read compatibility does not imply equal sequences. Preserve alternatives
        // until read linkage can resolve them; only identical sequences are collapsed.
        drop(normalized);
        let mut candidates: Vec<ContigRecord> = unique.into_values().collect();
        if args.profile {
            eprintln!(
                "ITS2 normalization {:.3}s",
                normalize_started.elapsed().as_secs_f64()
            );
        }
        let span_started = Instant::now();
        for candidate in &mut candidates {
            candidate.its2_full_reference_span =
                uce_filter_core::alignment::align_marker(&marker_index, &candidate.sequence, 0, 32)
                    .is_some_and(|alignment| {
                        alignment.reference_start == 0
                            && alignment.reference_end as usize
                                == marker_index.references[alignment.sequence as usize]
                                    .bases
                                    .len()
                            && alignment.identity() >= 0.8
                    });
        }
        if args.profile {
            eprintln!(
                "ITS2 span annotation {:.3}s",
                span_started.elapsed().as_secs_f64()
            );
        }
        let evidence_started = Instant::now();
        candidates.sort_by(|a, b| a.sequence.cmp(&b.sequence));
        for (index, candidate) in candidates.iter_mut().enumerate() {
            candidate.label = format!("its2_{:06}", index + 1);
            if candidate.equivalence_members.is_empty() {
                candidate.equivalence_members = candidate.label.clone();
            }
        }
        // Terminally clipped reads can cover a valid marker boundary without
        // spanning the central read slice used by the legacy assembly score.
        annotate_its2_candidates_with_quality(
            &mut candidates,
            &linked_fragments,
            Some(&linked_qualities),
            current_k,
        );
        if args.profile {
            eprintln!(
                "ITS2 read evidence {:.3}s",
                evidence_started.elapsed().as_secs_f64()
            );
        }
        accepted = candidates
            .iter()
            .filter(|candidate| candidate.accepted)
            .cloned()
            .collect();
        rejected = candidates
            .into_iter()
            .filter(|candidate| !candidate.accepted)
            .collect();
    }

    let mut low_quality = accepted.is_empty();
    if low_quality && args.assembly_mode == AssemblyMode::Reference {
        accepted = rejected.clone();
        low_quality = accepted.is_empty();
    }
    let pool = if accepted.is_empty() {
        &rejected
    } else {
        &accepted
    };
    let Some(mut best) = pool
        .iter()
        .max_by(|left, right| compare_contigs(left, right, args.assembly_mode))
        .cloned()
    else {
        clean_locus_outputs(args, &task.key);
        log_line(
            &args.output,
            log_lock,
            "Insufficient reads coverage, unable to build contigs.",
        );
        return Ok(LocusResult::failure(&task.key, "no contigs"));
    };
    if paired_fragments {
        if let Some(core_best) = core_accepted
            .iter()
            .max_by(|left, right| compare_contigs(left, right, args.assembly_mode))
        {
            // PE may improve evidence-based ranking, but it must never shorten
            // the final locus relative to its independently assembled core path.
            if best.sequence.len() < core_best.sequence.len() {
                best = core_best.clone();
            }
        }
    }

    accepted.sort_by(|left, right| compare_contigs(right, left, args.assembly_mode));
    rejected.sort_by(|left, right| compare_contigs(right, left, args.assembly_mode));

    if args.assembly_mode == AssemblyMode::Its2 {
        write_contigs(&best_path, &accepted, args.assembly_mode, "its2_variant")?;
        write_contigs(&all_path, &accepted, args.assembly_mode, "its2_variant")?;
        let mut support = accepted.clone();
        support.extend(rejected.clone());
        write_its2_support(
            &args
                .output
                .join("results")
                .join(format!("{}.its2_support.tsv", task.key)),
            &support,
        )?;
    } else if args.assembly_mode == AssemblyMode::Reference || !low_quality {
        write_contigs(
            &best_path,
            std::slice::from_ref(&best),
            args.assembly_mode,
            "contig",
        )?;
        write_contigs(&all_path, &accepted, args.assembly_mode, "contig")?;
    } else {
        remove_if_exists(&best_path);
        remove_if_exists(&all_path);
    }
    if matches!(args.assembly_mode, AssemblyMode::Uce | AssemblyMode::Its2) && !rejected.is_empty()
    {
        write_contigs(
            &low_path,
            &rejected,
            args.assembly_mode,
            "low_support_contig",
        )?;
    } else {
        remove_if_exists(&low_path);
    }

    let accepted_locus = args.assembly_mode == AssemblyMode::Reference || !low_quality;
    let length = best.sequence.len();
    if let Some(profile) = profile {
        profile.record_qc(qc_started);
    }
    Ok(LocusResult {
        key: task.key.clone(),
        status: if low_quality {
            "low quality".to_string()
        } else {
            "success".to_string()
        },
        value: best.read_count,
        accepted: accepted_locus,
        rejection_reason: if low_quality {
            best.rejection_reason.clone()
        } else {
            String::new()
        },
        selected_contig_length: length,
        read_supported_span: best.supported_span,
        slice_supported_bases: best.supported_bases,
        slice_support_breadth: best.support_breadth,
        max_slice_support_gap: best.max_support_gap,
        read_count: best.read_count,
        unique_read_count: best.unique_read_count,
        multi_mapping_read_count: best.multi_mapping_read_count,
        read_density: best.read_density,
        unique_read_density: if length > 0 {
            best.unique_read_count as f64 / length as f64
        } else {
            0.0
        },
        support_fraction: best.support_fraction,
        flank_balance: best.flank_balance,
        kmer_median_depth: best.kmer_median_depth,
        kmer_depth_cv: best.kmer_depth_cv,
        kmer_max_depth_ratio: best.kmer_max_depth_ratio,
        candidate_count: accepted.len() + rejected.len(),
        low_quality,
        skipped: false,
    })
}

pub fn read_result_dict(path: &Path) -> io::Result<HashMap<String, (String, u64)>> {
    let mut results = HashMap::new();
    if !path.is_file() {
        return Ok(results);
    }
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() >= 3 && !fields[0].is_empty() {
            results.insert(
                fields[0].to_string(),
                (fields[1].to_string(), fields[2].parse().unwrap_or(0)),
            );
        }
    }
    Ok(results)
}

pub fn write_result_dict(path: &Path, results: &HashMap<String, (String, u64)>) -> io::Result<()> {
    let mut keys: Vec<&String> = results.keys().collect();
    keys.sort();
    let mut writer = BufWriter::new(File::create(path)?);
    for key in keys {
        let (status, value) = &results[key];
        writeln!(writer, "{key},{status},{value},")?;
    }
    writer.flush()
}

pub fn read_summary_lines(path: &Path) -> io::Result<HashMap<String, String>> {
    let mut rows = HashMap::new();
    if !path.is_file() {
        return Ok(rows);
    }
    for (index, line) in BufReader::new(File::open(path)?).lines().enumerate() {
        let line = line?;
        if index == 0 {
            continue;
        }
        let key = line.split(',').next().unwrap_or("");
        if !key.is_empty() {
            rows.insert(key.to_string(), line);
        }
    }
    Ok(rows)
}

fn rounded(value: f64, decimals: usize) -> String {
    let mut text = format!("{value:.decimals$}");
    while text.contains('.') && text.ends_with('0') {
        text.pop();
    }
    if text.ends_with('.') {
        text.pop();
    }
    if text == "-0" {
        "0".to_string()
    } else {
        text
    }
}

fn rounded_float(value: f64, decimals: usize) -> String {
    let mut text = rounded(value, decimals);
    if !text.contains('.') {
        text.push_str(".0");
    }
    text
}

// 汇总行把状态和关键指标压成一条，统计脚本拿来就能用。
pub fn summary_line(result: &LocusResult) -> String {
    format!(
        "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
        result.key,
        result.status,
        u8::from(result.accepted),
        result.rejection_reason,
        result.selected_contig_length,
        result.read_supported_span,
        result.slice_supported_bases,
        rounded(result.slice_support_breadth, 6),
        result.max_slice_support_gap,
        result.read_count,
        result.unique_read_count,
        result.multi_mapping_read_count,
        rounded(result.read_density, 6),
        rounded(result.unique_read_density, 6),
        rounded(result.support_fraction, 3),
        rounded(result.flank_balance, 3),
        rounded_float(result.kmer_median_depth, 3),
        rounded_float(result.kmer_depth_cv, 3),
        rounded_float(result.kmer_max_depth_ratio, 3),
        result.candidate_count,
        u8::from(result.low_quality),
    )
}

// 汇总按 locus 名排序再写，重复运行也不至于顺序乱飘。
pub fn write_summary(path: &Path, rows: &HashMap<String, String>) -> io::Result<()> {
    let mut keys: Vec<&String> = rows.keys().collect();
    keys.sort();
    let mut writer = BufWriter::new(File::create(path)?);
    write!(writer, "{SUMMARY_HEADER}\r\n")?;
    for key in keys {
        write!(writer, "{}\r\n", rows[key])?;
    }
    writer.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn candidate(label: &str, sequence: &[u8]) -> ContigRecord {
        ContigRecord {
            label: label.to_string(),
            equivalence_members: label.to_string(),
            sequence: sequence.to_vec(),
            its2_full_reference_span: true,
            ..ContigRecord::default()
        }
    }

    #[test]
    fn phase_graph_accepts_nonadjacent_bridges_but_not_duplicates() {
        let one = vec![(0, 100, false)];
        let two = vec![(1, 101, false)];
        let mut links = HashMap::new();
        links.insert((0, 2), HashSet::from([one.clone(), two.clone()]));
        links.insert((1, 2), HashSet::from([one.clone()]));
        assert_eq!(its2_phase_component_count(3, &links), 2);
        links.get_mut(&(1, 2)).unwrap().insert(one);
        assert_eq!(its2_phase_component_count(3, &links), 2);
        links.get_mut(&(1, 2)).unwrap().insert(two);
        assert_eq!(its2_phase_component_count(3, &links), 1);
        assert_eq!(its2_phase_component_count(4, &links), 2);
        assert_eq!(its2_phase_component_count(0, &HashMap::new()), 0);
    }

    #[test]
    fn likelihood_uses_quality_without_double_counting_mate_overlap() {
        let candidates = vec![candidate("a", b"ACGT"), candidate("b", b"ATGT")];
        let p = Its2Placement {
            start: 0,
            end: 4,
            reverse: false,
            read_start: 0,
            read_end: 4,
            mismatch: None,
        };
        let placements = vec![vec![vec![p], vec![p]], vec![vec![p], vec![p]]];
        let fragment = vec![b"ACGT".to_vec(), b"ACGT".to_vec()];
        let q = vec![vec![12; 4], vec![12; 4]];
        let w = its2_likelihood_weights(&[0, 1], &candidates, &placements, &fragment, Some(&q));
        let e = 10_f64.powf(-1.2);
        assert!((w[1] - (e / 3.0) / (1.0 - e)).abs() < 1e-12);
        assert_eq!(w[0], 1.0);
        let high = vec![vec![40; 4], vec![40; 4]];
        assert!(
            its2_likelihood_weights(&[0, 1], &candidates, &placements, &fragment, Some(&high))[1]
                < w[1]
        );
        let conflicting = vec![b"ACGT".to_vec(), b"ATGT".to_vec()];
        assert_eq!(
            its2_likelihood_weights(&[0, 1], &candidates, &placements, &conflicting, Some(&q)),
            vec![1.0, 1.0]
        );
        assert_eq!(
            its2_likelihood_weights(&[0, 1], &candidates, &placements, &fragment, None),
            vec![1.0, 1.0]
        );
        let zero = vec![vec![0; 4], vec![0; 4]];
        assert_eq!(
            its2_likelihood_weights(&[0, 1], &candidates, &placements, &fragment, Some(&zero)),
            vec![1.0, 1.0]
        );
        let reverse_reads: Vec<_> = fragment.iter().map(|r| reverse_complement(r)).collect();
        let mut reverse_placements = placements.clone();
        for candidate in &mut reverse_placements {
            for mate in candidate {
                mate[0].reverse = true;
            }
        }
        let reverse_quality = vec![vec![40, 40, 12, 40], vec![40, 40, 12, 40]];
        assert_eq!(
            its2_likelihood_weights(
                &[0, 1],
                &candidates,
                &reverse_placements,
                &reverse_reads,
                Some(&reverse_quality)
            ),
            w
        );
        let mut clipped = placements;
        clipped[1][0][0].end = 3;
        assert_eq!(
            its2_likelihood_weights(&[0, 1], &candidates, &clipped, &fragment, Some(&q)),
            vec![1.0, 1.0]
        );
    }

    #[test]
    fn quality_rescue_preserves_error_positions_and_orientation() {
        let sequence = phase_fixture().remove(0).sequence;
        let candidates = vec![candidate("a", &sequence)];
        let index = Its2PlacementIndex::new(&candidates, 21);
        let mut read = sequence.clone();
        read[30] = if read[30] == b'A' { b'C' } else { b'A' };
        let mut q = vec![40; read.len()];
        q[30] = 12;
        let p = index.place(&read, Some(&q));
        assert_eq!(p[0].len(), 1);
        assert_eq!(p[0][0].mismatch, Some(30));
        let reverse = reverse_complement(&read);
        q.reverse();
        let p = index.place(&reverse, Some(&q));
        assert_eq!(p[0].len(), 1);
        assert!(p[0][0].reverse);
        assert_eq!(p[0][0].mismatch, Some(30));
        assert!(index.place(&read, Some(&vec![40; read.len()]))[0].is_empty());
        assert!(index.place(&read[..74], Some(&[12; 74]))[0].is_empty());
        read[35] = if read[35] == b'A' { b'C' } else { b'A' };
        assert!(index.place(&read, Some(&vec![12; read.len()]))[0].is_empty());
        read[35] = b'N';
        assert!(index.place(&read, Some(&vec![12; read.len()]))[0].is_empty());
    }

    #[test]
    fn quality_seed_index_matches_exhaustive_error_search() {
        let mut state = 719_u64;
        let mut random = |n| {
            (0..n)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(state >> 32) as usize % 4]
                })
                .collect::<Vec<_>>()
        };
        for round in 0..200 {
            let sequence = random(180);
            let candidates = vec![candidate("a", &sequence)];
            let index = Its2PlacementIndex::new(&candidates, 21);
            let start = round % 70;
            let mut read = sequence[start..start + 90].to_vec();
            let at = round % 90;
            read[at] = if read[at] == b'A' { b'C' } else { b'A' };
            let mut quality = vec![40; 90];
            quality[at] = if round % 3 == 0 { 40 } else { 12 };
            if round % 2 == 0 {
                read = reverse_complement(&read);
                quality.reverse();
            }
            let mut expected = Vec::new();
            let reverse = reverse_complement(&read);
            for (strand, q) in [read.as_slice(), reverse.as_slice()]
                .into_iter()
                .enumerate()
            {
                for offset in -69_isize..=159 {
                    let a = (-offset).max(0) as usize;
                    let b = offset.max(0) as usize;
                    let length = (90 - a).min(180 - b);
                    let errors: Vec<_> = (0..length)
                        .filter(|&j| q[a + j] != sequence[b + j])
                        .collect();
                    let allowed = errors.is_empty()
                        || (length >= 75
                            && errors.len() == 1
                            && quality[if strand == 0 {
                                a + errors[0]
                            } else {
                                89 - a - errors[0]
                            }] <= 12);
                    if allowed {
                        expected.push(Its2Placement {
                            start: b,
                            end: b + length,
                            reverse: strand == 1,
                            read_start: offset,
                            read_end: offset + 90,
                            mismatch: errors.first().map(|j| b + j),
                        });
                    }
                }
            }
            if let Some(longest) = expected.iter().map(|p| p.end - p.start).max() {
                expected.retain(|p| p.end - p.start == longest);
            }
            assert_eq!(
                index.place(&read, Some(&quality))[0],
                expected,
                "round {round}"
            );
        }
    }

    #[test]
    fn quality_rescue_never_fills_the_mismatching_base() {
        let sequence = phase_fixture().remove(0).sequence;
        let mut read = sequence.clone();
        read[30] = if read[30] == b'A' { b'C' } else { b'A' };
        let mut q = vec![40; read.len()];
        q[30] = 12;
        let mut candidates = vec![candidate("a", &sequence)];
        annotate_its2_candidates_with_quality(
            &mut candidates,
            &vec![vec![read]; 2],
            Some(&vec![vec![q]; 2]),
            21,
        );
        assert_eq!(candidates[0].quality_rescued_fragment_support, 2);
        assert_eq!(candidates[0].supported_bases, sequence.len() - 1);
        assert!(!candidates[0].accepted);
    }

    #[test]
    fn low_quality_alleles_cannot_manufacture_a_minor_haplotype() {
        let sequence = phase_fixture().remove(0).sequence;
        let mut minor = sequence.clone();
        minor[30] = if minor[30] == b'A' { b'C' } else { b'A' };
        let mut candidates = vec![candidate("major", &sequence), candidate("minor", &minor)];
        let mut fragments = Vec::new();
        let mut qualities = Vec::new();
        for (seq, allele_q) in [(&sequence, 40), (&minor, 12)] {
            for flank in 0..2 {
                let mut read = vec![b'T'; flank];
                read.extend(seq);
                let mut q = vec![40; read.len()];
                q[30 + flank] = allele_q;
                fragments.push(vec![read]);
                qualities.push(vec![q]);
            }
        }
        let mut without_quality = candidates.clone();
        annotate_its2_candidates(&mut without_quality, &fragments, 21);
        assert!(
            without_quality[1].accepted,
            "exact-only evidence misses low-quality allele uncertainty"
        );
        annotate_its2_candidates_with_quality(&mut candidates, &fragments, Some(&qualities), 21);
        assert!(candidates[0].accepted);
        assert!(candidates[0].quality_rescued_fragment_support > 0);
        assert!(!candidates[1].accepted);
        assert_eq!(candidates[1].diagnostic_fragment_support, 0);
        // The same minor allele with high-quality independent observations survives.
        for q in &mut qualities {
            q[0].fill(40);
        }
        annotate_its2_candidates_with_quality(&mut candidates, &fragments, Some(&qualities), 21);
        assert!(candidates.iter().all(|c| c.accepted));
    }

    #[test]
    fn quality_rescue_cannot_establish_a_bridge_through_its_error() {
        let read = vec![b'A'; 100];
        let quality = vec![vec![40; 100]];
        let p = Its2Placement {
            start: 0,
            end: 100,
            read_start: 0,
            read_end: 100,
            reverse: false,
            mismatch: Some(30),
        };
        assert!(!its2_phase_base(
            30,
            &[vec![p]],
            &[read.clone()],
            Some(&quality)
        ));
        assert!(its2_phase_base(31, &[vec![p]], &[read], Some(&quality)));
    }

    #[test]
    fn its2_index_matches_exhaustive_placements() {
        let mut state = 918273_u64;
        let mut random = |length: usize| -> Vec<u8> {
            (0..length)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(state >> 32) as usize % 4]
                })
                .collect()
        };
        for round in 0..400 {
            let sequence = random(30 + round % 90);
            let mut repeated = sequence.clone();
            repeated.extend_from_slice(&sequence);
            let candidates = vec![
                candidate("a", &sequence),
                candidate("repeat", &repeated),
                candidate("short", &sequence[..3]),
                candidate("empty", b""),
                candidate("N", b"ACGTNNNNACGT"),
            ];
            let minimum = round % 35;
            let index = Its2PlacementIndex::new(&candidates, minimum);
            let start = round % sequence.len();
            let mut read = sequence[start..].to_vec();
            if round % 3 == 0 {
                read.splice(0..0, random(13));
            }
            if round % 4 == 0 {
                read.extend(random(19));
            }
            if round % 5 == 0 {
                let middle = read.len() / 2;
                read[middle] = b'N';
            }
            if round % 2 == 0 {
                read = reverse_complement(&read);
            }
            let observed = index.placements(&read);
            for (i, c) in candidates.iter().enumerate() {
                assert_eq!(
                    observed[i],
                    its2_placements(&read, &c.sequence, minimum),
                    "round {round}, candidate {i}"
                );
            }
        }
    }

    #[test]
    fn its2_matching_checks_bases_outside_the_seed() {
        assert!(!mate_matches_candidate(b"TTTACGACCC", b"GGGACGAGGG", 4));
        assert!(mate_matches_candidate(b"TTTACGACCC", b"ATTTACGACCCA", 4));
        assert!(mate_matches_candidate(b"GGGTCGTAAA", b"ATTTACGACCCA", 4));
        // A read extending past an assembled end still supplies overlap evidence.
        assert!(mate_matches_candidate(b"TTTACGACCC", b"ACGACCCA", 4));
        assert!(!mate_matches_candidate(b"NNNN", b"NNNN", 4));
        assert!(!mate_matches_candidate(b"ACGT", b"ACGT", 0));
    }

    fn phase_fixture() -> Vec<ContigRecord> {
        let base =
            b"ACGTTGCAAGTCGATCGTACCTGAGCTAGTCGACCTAGGCTACGATGCTAGCTTGACGATCGTACGAGTCCAGTGCATGA";
        let mut candidates = Vec::new();
        for (index, (left, right)) in [(b'A', b'A'), (b'A', b'C'), (b'C', b'A'), (b'C', b'C')]
            .into_iter()
            .enumerate()
        {
            let mut sequence = base.to_vec();
            sequence[12] = left;
            sequence[65] = right;
            candidates.push(candidate(&format!("v{index}"), &sequence));
        }
        candidates
    }

    #[test]
    fn its2_local_coverage_does_not_establish_phase() {
        let mut candidates = phase_fixture();
        let mut fragments = Vec::new();
        for c in &candidates {
            for _ in 0..2 {
                fragments.push(vec![c.sequence[..40].to_vec()]);
                fragments.push(vec![c.sequence[40..].to_vec()]);
            }
        }
        annotate_its2_candidates(&mut candidates, &fragments, 16);
        for c in &candidates {
            assert_eq!(c.support_breadth, 1.0);
            assert!(!c.accepted);
            assert_eq!(its2_status(c), "UNPHASED");
            assert!(c.em_abundance.is_nan());
        }
    }

    #[test]
    fn its2_shared_middle_allele_does_not_resolve_competing_paths() {
        let mut candidates = phase_fixture();
        let mut insertion = candidates[0].sequence.clone();
        insertion.insert(40, b'G');
        candidates.push(candidate("insertion", &insertion));
        let mut fragments = Vec::new();
        for c in &candidates[..4] {
            for _ in 0..2 {
                fragments.push(vec![c.sequence[..40].to_vec()]);
                fragments.push(vec![c.sequence[40..].to_vec()]);
                fragments.push(vec![c.sequence[20..60].to_vec()]);
            }
        }
        annotate_its2_candidates(&mut candidates, &fragments, 16);
        for c in &candidates[..4] {
            assert_eq!(c.support_breadth, 1.0);
            assert_eq!(c.diagnostic_fragment_support, 0);
            assert!(!c.accepted, "shared local links must not phase {}", c.label);
            assert_eq!(its2_status(c), "UNPHASED");
        }
    }

    #[test]
    fn its2_valid_pairs_link_distant_variants() {
        let mut candidates = phase_fixture();
        let mut fragments = Vec::new();
        for c in &candidates {
            for shift in 0..2 {
                fragments.push(vec![
                    c.sequence[..40 + shift].to_vec(),
                    reverse_complement(&c.sequence[40 - shift..]),
                ]);
            }
        }
        annotate_its2_candidates(&mut candidates, &fragments, 16);
        for c in &candidates {
            assert!(c.accepted, "{}: {}", c.label, c.rejection_reason);
            assert_eq!(c.paired_fragment_support, 2);
            assert_eq!(c.diagnostic_fragment_support, 2);
            assert_eq!(its2_status(c), "READ_SUPPORTED");
            assert!((c.em_abundance - 0.25).abs() < 1e-6);
        }
    }

    #[test]
    fn its2_local_hits_cannot_accept_uncovered_sequence() {
        let mut candidates = vec![phase_fixture().remove(0)];
        let fragments = vec![vec![candidates[0].sequence[..40].to_vec()]; 3];
        annotate_its2_candidates(&mut candidates, &fragments, 16);
        assert!(!candidates[0].accepted);
        assert_eq!(its2_status(&candidates[0]), "PARTIAL");
        assert!(candidates[0].max_support_gap > 0);
        assert_eq!(candidates[0].phase_components, None);
        assert_eq!(candidates[0].min_variant_site_support, None);
    }

    #[test]
    fn its2_uncovered_extension_cannot_block_supported_marker() {
        let core = phase_fixture().remove(0).sequence;
        let mut extended = core.clone();
        extended.extend_from_slice(b"GGGTTTAAACCCGGGTTTAAACCC");
        let mut candidates = vec![candidate("core", &core), candidate("extension", &extended)];
        let fragments = vec![
            vec![core.clone()],
            vec![core[..60].to_vec()],
            vec![core[20..].to_vec()],
        ];
        annotate_its2_candidates(&mut candidates, &fragments, 16);
        assert!(candidates[0].accepted);
        assert_eq!(candidates[0].diagnostic_fragment_support, 3);
        assert!(!candidates[1].accepted);
        assert_eq!(its2_status(&candidates[1]), "PARTIAL");
    }

    #[test]
    fn phase_audit_distinguishes_duplicates_low_quality_and_independent_bridges() {
        let all = phase_fixture();
        let mut candidates = vec![all[0].clone(), all[2].clone(), all[3].clone()];
        let sequence = candidates[0].sequence.clone();
        let fragments = vec![
            vec![sequence.clone()],
            vec![sequence.clone()],
            vec![sequence[1..].to_vec()],
        ];
        let mut qualities: LinkedReads = fragments
            .iter()
            .map(|f| f.iter().map(|r| vec![40; r.len()]).collect())
            .collect();
        qualities[2][0][64] = 12;
        annotate_its2_evidence_with_quality(
            &mut candidates,
            &fragments,
            Some(&qualities),
            16,
            true,
        );
        assert_eq!(candidates[0].fragment_support, 3);
        assert_eq!(candidates[0].min_variant_site_support, Some(1));
        assert_eq!(candidates[0].phase_components, Some(2));
        assert!(!candidates[0].accepted);
        qualities[2][0][64] = 40;
        annotate_its2_evidence_with_quality(
            &mut candidates,
            &fragments,
            Some(&qualities),
            16,
            true,
        );
        assert_eq!(candidates[0].min_variant_site_support, Some(2));
        assert_eq!(candidates[0].phase_components, Some(1));
        assert!(candidates[0].accepted);
        annotate_its2_evidence_with_quality(
            &mut candidates,
            &fragments,
            Some(&qualities),
            16,
            false,
        );
        assert_eq!(candidates[0].phase_components, None);
        assert_eq!(candidates[0].min_variant_site_support, None);
    }

    #[test]
    fn its2_duplicate_bridge_cannot_manufacture_phase_links() {
        let all = phase_fixture();
        // Only the first allele distinguishes the target from both competitors.
        // Local diagnostic reads are not independent bridges to the second site.
        let mut candidates = vec![all[0].clone(), all[2].clone(), all[3].clone()];
        let sequence = candidates[0].sequence.clone();
        let fragments = vec![
            vec![sequence.clone()],
            vec![sequence.clone()],
            vec![sequence[..35].to_vec()],
            vec![sequence[1..36].to_vec()],
        ];
        annotate_its2_evidence(&mut candidates, &fragments, 16, true);
        assert!(candidates[0].diagnostic_fragment_support >= 2);
        assert_eq!(candidates[0].supported_bases, sequence.len());
        assert!(!candidates[0].accepted);
        assert_eq!(candidates[0].rejection_reason, "unresolved_phase_links");
    }

    #[test]
    fn its2_duplicate_fragments_cannot_manufacture_diagnostic_support() {
        let core = phase_fixture().remove(0).sequence;
        let mut candidates = vec![candidate("core", &core)];
        annotate_its2_candidates(&mut candidates, &vec![vec![core]; 10], 16);
        assert_eq!(candidates[0].fragment_support, 10);
        assert_eq!(candidates[0].diagnostic_fragment_support, 1);
        assert!(!candidates[0].accepted);
    }

    #[test]
    fn its2_contained_fragment_does_not_compete_with_full_marker() {
        let core = phase_fixture().remove(0).sequence;
        let fragments = vec![
            vec![core[..60].to_vec()],
            vec![core[20..].to_vec()],
            vec![core.clone()],
        ];
        let mut candidates = vec![
            candidate("core", &core),
            candidate("fragment", &core[20..60]),
        ];
        annotate_its2_candidates(&mut candidates, &fragments, 16);
        assert!(candidates[0].accepted);
        assert!(!candidates[1].accepted);
        assert_eq!(candidates[1].rejection_reason, "contained_model");
    }

    #[test]
    fn partial_reference_model_cannot_shadow_a_full_marker() {
        let core = phase_fixture().remove(0).sequence;
        let fragments = vec![
            vec![core[..60].to_vec()],
            vec![core[20..].to_vec()],
            vec![core.clone()],
        ];
        let mut partial = candidate("partial", &core[20..60]);
        partial.its2_full_reference_span = false;
        let mut candidates = vec![candidate("full", &core), partial];
        annotate_its2_candidates(&mut candidates, &fragments, 16);
        assert!(candidates[0].accepted);
        assert!(!candidates[1].accepted);
        assert_eq!(candidates[1].rejection_reason, "partial_reference_span");
    }

    #[test]
    fn its2_terminal_reads_need_not_span_a_central_read_slice() {
        let core = phase_fixture().remove(0).sequence;
        let left = [vec![b'A'; 100], core[..51].to_vec()].concat();
        let right = [core[29..].to_vec(), vec![b'C'; 100]].concat();
        let fragments = vec![
            vec![left.clone()],
            vec![left],
            vec![right.clone()],
            vec![right],
        ];
        let mut candidates = vec![candidate("terminal", &core)];
        annotate_its2_candidates(&mut candidates, &fragments, 21);
        assert!(candidates[0].accepted);
        assert_eq!(candidates[0].supported_bases, core.len());
    }

    #[test]
    fn its2_same_strand_mates_do_not_establish_links() {
        let mut candidates = vec![phase_fixture().remove(0)];
        let fragments = vec![
            vec![
                candidates[0].sequence[..40].to_vec(),
                candidates[0].sequence[40..].to_vec()
            ];
            2
        ];
        annotate_its2_candidates(&mut candidates, &fragments, 16);
        assert_eq!(candidates[0].fragment_support, 0);
        assert!(!candidates[0].accepted);
    }

    #[test]
    fn its2_labelled_reference_beats_anonymous_graph_candidate() {
        let anonymous = candidate("", b"AAAATATA");
        let labelled = candidate("C40", b"AAAATATA");
        assert!(prefer_its2_candidate(
            &anonymous,
            &labelled,
            AssemblyMode::Its2
        ));
        assert!(!prefer_its2_candidate(
            &labelled,
            &anonymous,
            AssemblyMode::Its2
        ));
        let mut members = "C40".to_string();
        append_equivalence_member(&mut members, "C40e");
        append_equivalence_member(&mut members, "C40");
        assert_eq!(members, "C40;C40e");
    }
}
