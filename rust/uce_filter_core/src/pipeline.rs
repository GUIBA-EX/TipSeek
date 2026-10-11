use crate::alignment::{align_read, terminal_evidence};
use crate::evidence::{collect_runs_stats, infer_orientation};
use crate::index::{
    default_verify_k, ExactSeed, IndexProfile, ReadEvidenceScratch, RecruitScratch, UceIndex,
};
use crate::model::{default_spill_path, Candidate, Fragment, FragmentBank, LocusId};
use crate::selection::{choose_auto, choose_legacy, selected};
use gm2_tools::fastx::{gzip_backend_name, open_input, FastxRecord};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::thread::{self, JoinHandle};
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct Config {
    pub references: PathBuf,
    pub recruit_references: Option<PathBuf>,
    pub read1: PathBuf,
    pub read2: PathBuf,
    pub output: PathBuf,
    pub kmer_size: usize,
    pub verify_kmer_size: Option<usize>,
    pub step: usize,
    pub max_locus_count: usize,
    pub retain_loci_file: Option<PathBuf>,
    pub minimum_alignment_overlap: usize,
    pub minimum_alignment_identity: f64,
    pub min_depth: i64,
    pub max_depth: i64,
    pub max_size_mb: i64,
    pub max_fragments: u64,
    pub memory_limit_mib: u64,
    pub selection_auto: bool,
    pub reference_is_contig: bool,
    pub alignment_shadow: bool,
    pub shadow_per_locus: usize,
    pub shadow_band: usize,
    pub terminal_window: usize,
    pub profile: bool,
}

#[derive(Clone, Debug, Default)]
pub struct RunSummary {
    pub fragments_read: u64,
    pub fragments_retained_once: usize,
    pub fragment_bases_retained_once: u64,
    pub assignments: usize,
    pub loci_written: usize,
    pub fragment_memory_bytes: u64,
    pub fragment_spill_bytes: u64,
    pub evidence_scratch_bytes: u64,
    pub candidate_memory_bytes: u64,
    pub shadow_sampled_assignments: usize,
    pub shadow_aligned_mates: usize,
    pub shadow_seconds: f64,
    pub index_seconds: f64,
    pub scan_seconds: f64,
    pub selection_seconds: f64,
    pub output_seconds: f64,
    pub decode_seconds: f64,
    pub recruit_seconds: f64,
    pub evidence_seconds: f64,
    pub store_seconds: f64,
    pub index_profile: IndexProfile,
    pub elapsed_seconds: f64,
}

const LOCUS_BUFFER_BYTES: usize = 64 * 1024;
const DECODE_CHUNK_BYTES: usize = 1024 * 1024;
const DECODE_BUFFERS_PER_MATE: usize = 2;

fn passes_alignment_filter(
    index: &UceIndex,
    locus: LocusId,
    read1: &[u8],
    read2: &[u8],
    band: usize,
    minimum_overlap: usize,
    minimum_identity: f64,
) -> bool {
    if minimum_overlap == 0 && minimum_identity == 0.0 {
        return true;
    }
    [read1, read2].into_iter().any(|read| {
        align_read(index, read, locus, band).is_some_and(|alignment| {
            alignment.reference_overlap() >= minimum_overlap
                && alignment.identity() >= minimum_identity
        })
    })
}

fn retain_locus_mask(path: &Path, index: &UceIndex) -> Result<Vec<bool>, String> {
    let requested = fs::read_to_string(path)
        .map_err(|error| {
            format!(
                "Unable to read retain-loci file '{}': {error}",
                path.display()
            )
        })?
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect::<BTreeSet<_>>();
    let known = index
        .loci
        .iter()
        .enumerate()
        .map(|(id, locus)| (locus.name.as_str(), id))
        .collect::<BTreeMap<_, _>>();
    let unknown = requested
        .iter()
        .filter(|name| !known.contains_key(name.as_str()))
        .cloned()
        .collect::<Vec<_>>();
    if !unknown.is_empty() {
        return Err(format!(
            "retain-loci file contains unknown loci: {}",
            unknown.join(", ")
        ));
    }
    let mut mask = vec![false; index.loci.len()];
    for name in requested {
        mask[known[name.as_str()]] = true;
    }
    Ok(mask)
}

struct LocusOutput {
    path: PathBuf,
    buffer: Vec<u8>,
}

struct OutputRouter {
    outputs: Vec<Option<LocusOutput>>,
}

struct FragmentRoutes {
    offsets: Vec<u32>,
    loci: Vec<LocusId>,
}

const NO_CANDIDATE: u32 = u32::MAX;
const CANDIDATE_MIN_GROWTH: usize = 4096;

/// Append-only candidate arena. Per-locus linked indices preserve input order
/// without thousands of independently growing `Vec<Candidate>` allocations.
struct LocusCandidateStore {
    heads: Vec<u32>,
    tails: Vec<u32>,
    next: Vec<u32>,
    candidates: Vec<Candidate>,
}

impl LocusCandidateStore {
    fn new(locus_count: usize) -> Self {
        Self {
            heads: vec![NO_CANDIDATE; locus_count],
            tails: vec![NO_CANDIDATE; locus_count],
            next: Vec::new(),
            candidates: Vec::new(),
        }
    }

    fn push(&mut self, locus: LocusId, candidate: Candidate) -> Result<(), String> {
        if self.candidates.len() == self.candidates.capacity() {
            let growth = (self.candidates.capacity() / 4).max(CANDIDATE_MIN_GROWTH);
            self.candidates.reserve_exact(growth);
        }
        if self.next.len() == self.next.capacity() {
            let growth = (self.next.capacity() / 4).max(CANDIDATE_MIN_GROWTH);
            self.next.reserve_exact(growth);
        }
        let node = u32::try_from(self.candidates.len())
            .map_err(|_| "too many UCE locus candidates".to_string())?;
        let locus = locus as usize;
        let tail = self.tails[locus];
        if tail == NO_CANDIDATE {
            self.heads[locus] = node;
        } else {
            self.next[tail as usize] = node;
        }
        self.tails[locus] = node;
        self.next.push(NO_CANDIDATE);
        self.candidates.push(candidate);
        Ok(())
    }

    fn copy_locus(&self, locus: usize, output: &mut Vec<Candidate>) {
        output.clear();
        let mut node = self.heads[locus];
        while node != NO_CANDIDATE {
            let index = node as usize;
            output.push(self.candidates[index]);
            node = self.next[index];
        }
    }

    fn allocated_bytes(&self) -> usize {
        (self.heads.capacity() + self.tails.capacity() + self.next.capacity())
            * std::mem::size_of::<u32>()
            + self.candidates.capacity() * std::mem::size_of::<Candidate>()
    }
}

impl FragmentRoutes {
    fn from_pairs(fragment_count: usize, pairs: &[(u32, LocusId)]) -> Result<Self, String> {
        if pairs.len() > u32::MAX as usize {
            return Err("too many UCE locus assignments".to_string());
        }
        let mut offsets = vec![0_u32; fragment_count + 1];
        for &(fragment, _) in pairs {
            offsets[fragment as usize + 1] += 1;
        }
        for fragment in 0..fragment_count {
            offsets[fragment + 1] += offsets[fragment];
        }
        let mut cursors = offsets[..fragment_count].to_vec();
        let mut loci = vec![0_u32; pairs.len()];
        for &(fragment, locus) in pairs {
            let cursor = &mut cursors[fragment as usize];
            loci[*cursor as usize] = locus;
            *cursor += 1;
        }
        Ok(Self { offsets, loci })
    }

    fn get(&self, fragment: u32) -> &[LocusId] {
        let fragment = fragment as usize;
        &self.loci[self.offsets[fragment] as usize..self.offsets[fragment + 1] as usize]
    }
}

struct ShadowWriter {
    writer: BufWriter<File>,
    band: usize,
    terminal_window: usize,
    sampled_assignments: usize,
    aligned_mates: usize,
    elapsed_seconds: f64,
    loci: Vec<ShadowLocusStats>,
}

#[derive(Clone, Debug, Default)]
struct ShadowLocusStats {
    sampled_assignments: u64,
    aligned_mates: u64,
    both_mates_aligned: u64,
    linked_mate_only: u64,
    multi_locus_assignments: u64,
    identity_sum: f64,
    identity_ge_90: u64,
    overlap_sum: u64,
    near_terminal_assignments: u64,
    extends_left_assignments: u64,
    extends_right_assignments: u64,
    covered_bins: u64,
}

impl ShadowWriter {
    fn new(
        path: &Path,
        band: usize,
        terminal_window: usize,
        locus_count: usize,
    ) -> Result<Self, String> {
        let mut writer = BufWriter::new(File::create(path).map_err(|e| e.to_string())?);
        writeln!(
            writer,
            "locus\tfragment_ordinal\tmate\tselected_locus_count\tstatus\treference_id\tstrand\treference_length\tterminal_window\tquery_length\tquery_start\tquery_end\treference_start\treference_end\texact_seed_length\tscore\tmatches\tmismatches\tgap_bases\toverlap\tidentity\tleft_overhang\tright_overhang\tnear_left_terminal\tnear_right_terminal\textends_left_terminal\textends_right_terminal"
        )
        .map_err(|e| e.to_string())?;
        Ok(Self {
            writer,
            band,
            terminal_window,
            sampled_assignments: 0,
            aligned_mates: 0,
            elapsed_seconds: 0.0,
            loci: vec![ShadowLocusStats::default(); locus_count],
        })
    }

    fn write_fragment(
        &mut self,
        index: &UceIndex,
        locus: LocusId,
        locus_name: &str,
        fragment: &Fragment,
        selected_locus_count: usize,
    ) -> Result<(), String> {
        self.sampled_assignments += 1;
        let mut mapped_mates = 0_u64;
        let mut near_terminal = false;
        let mut extends_left = false;
        let mut extends_right = false;
        for (mate, record) in [(1, &fragment.r1), (2, &fragment.r2)] {
            let started = Instant::now();
            let alignment = align_read(index, &record.sequence, locus, self.band);
            self.elapsed_seconds += started.elapsed().as_secs_f64();
            let Some(alignment) = alignment else {
                writeln!(
                    self.writer,
                    "{locus_name}\t{}\t{mate}\t{selected_locus_count}\tunmapped\tNA\tNA\tNA\tNA\t{}\tNA\tNA\tNA\tNA\t0\tNA\t0\t0\t0\t0\t0.000000\t0\t0\t0\t0\t0\t0",
                    fragment.ordinal,
                    record.sequence.len(),
                )
                .map_err(|e| e.to_string())?;
                continue;
            };
            self.aligned_mates += 1;
            mapped_mates += 1;
            let reference = &index.references[alignment.sequence as usize];
            let terminal = terminal_evidence(
                alignment,
                record.sequence.len(),
                reference.bases.len(),
                self.terminal_window,
            );
            near_terminal |= terminal.near_left_terminal || terminal.near_right_terminal;
            extends_left |= terminal.extends_left_terminal;
            extends_right |= terminal.extends_right_terminal;
            let stats = &mut self.loci[locus as usize];
            stats.aligned_mates += 1;
            stats.identity_sum += alignment.identity();
            stats.identity_ge_90 += u64::from(alignment.identity() >= 0.90);
            stats.overlap_sum += alignment.reference_overlap() as u64;
            let reference_length = reference.bases.len().max(1);
            let first_bin =
                (terminal.original_reference_start as usize * 64 / reference_length).min(63);
            let last_position = (terminal.original_reference_end as usize).saturating_sub(1);
            let last_bin = (last_position * 64 / reference_length).min(63);
            for bin in first_bin..=last_bin {
                stats.covered_bins |= 1_u64 << bin;
            }
            writeln!(
                self.writer,
                "{locus_name}\t{}\t{mate}\t{selected_locus_count}\taligned\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}",
                fragment.ordinal,
                alignment.sequence,
                alignment.strand,
                reference.bases.len(),
                terminal.effective_window,
                record.sequence.len(),
                alignment.query_start,
                alignment.query_end,
                terminal.original_reference_start,
                terminal.original_reference_end,
                alignment.exact_seed_length,
                alignment.score,
                alignment.matches,
                alignment.mismatches,
                alignment.gap_bases,
                alignment.reference_overlap(),
                alignment.identity(),
                terminal.left_overhang,
                terminal.right_overhang,
                terminal.near_left_terminal as u8,
                terminal.near_right_terminal as u8,
                terminal.extends_left_terminal as u8,
                terminal.extends_right_terminal as u8,
            )
            .map_err(|e| e.to_string())?;
        }
        let stats = &mut self.loci[locus as usize];
        stats.sampled_assignments += 1;
        stats.both_mates_aligned += u64::from(mapped_mates == 2);
        stats.linked_mate_only += u64::from(mapped_mates == 1);
        stats.multi_locus_assignments += u64::from(selected_locus_count > 1);
        stats.near_terminal_assignments += u64::from(near_terminal);
        stats.extends_left_assignments += u64::from(extends_left);
        stats.extends_right_assignments += u64::from(extends_right);
        Ok(())
    }

    fn flush(&mut self) -> Result<(), String> {
        self.writer.flush().map_err(|e| e.to_string())
    }

    fn write_summary(&self, path: &Path, index: &UceIndex) -> Result<(), String> {
        let mut writer = BufWriter::new(File::create(path).map_err(|e| e.to_string())?);
        writeln!(
            writer,
            "locus\tsampled_assignments\taligned_mates\tboth_mates_aligned\tlinked_mate_only\tmulti_locus_assignments\tmean_identity\tidentity_ge_0.90\tmean_overlap\tnear_terminal_assignments\textends_left_assignments\textends_right_assignments\tcovered_bins_64"
        )
        .map_err(|e| e.to_string())?;
        for (locus, stats) in index.loci.iter().zip(&self.loci) {
            if stats.sampled_assignments == 0 {
                continue;
            }
            let mean_identity = if stats.aligned_mates == 0 {
                0.0
            } else {
                stats.identity_sum / stats.aligned_mates as f64
            };
            let mean_overlap = if stats.aligned_mates == 0 {
                0.0
            } else {
                stats.overlap_sum as f64 / stats.aligned_mates as f64
            };
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{:.3}\t{}\t{}\t{}\t{}",
                locus.name,
                stats.sampled_assignments,
                stats.aligned_mates,
                stats.both_mates_aligned,
                stats.linked_mate_only,
                stats.multi_locus_assignments,
                mean_identity,
                stats.identity_ge_90,
                mean_overlap,
                stats.near_terminal_assignments,
                stats.extends_left_assignments,
                stats.extends_right_assignments,
                stats.covered_bins.count_ones(),
            )
            .map_err(|e| e.to_string())?;
        }
        writer.flush().map_err(|e| e.to_string())
    }
}

fn evenly_sample(ids: &[u32], limit: usize) -> Vec<u32> {
    if limit == 0 || ids.len() <= limit {
        return ids.to_vec();
    }
    if limit == 1 {
        return vec![ids[ids.len() / 2]];
    }
    (0..limit)
        .map(|i| ids[i * (ids.len() - 1) / (limit - 1)])
        .collect()
}

impl OutputRouter {
    fn new(locus_count: usize) -> Self {
        Self {
            outputs: (0..locus_count).map(|_| None).collect(),
        }
    }

    fn register(&mut self, locus: usize, path: PathBuf) -> Result<(), String> {
        File::create(&path).map_err(|e| e.to_string())?;
        self.outputs[locus] = Some(LocusOutput {
            path,
            buffer: Vec::with_capacity(LOCUS_BUFFER_BYTES + 1024),
        });
        Ok(())
    }

    fn write_fragment(&mut self, locus: usize, fragment: &Fragment) -> Result<(), String> {
        let should_flush = {
            let output = self.outputs[locus]
                .as_mut()
                .ok_or_else(|| format!("locus {locus} was not registered"))?;
            write_record(&mut output.buffer, &fragment.r1)?;
            write_record(&mut output.buffer, &fragment.r2)?;
            output.buffer.len() >= LOCUS_BUFFER_BYTES
        };
        if should_flush {
            self.flush_one(locus)?;
        }
        Ok(())
    }

    fn flush_one(&mut self, locus: usize) -> Result<(), String> {
        let Some(output) = self.outputs[locus].as_mut() else {
            return Ok(());
        };
        if output.buffer.is_empty() {
            return Ok(());
        }
        let mut file = OpenOptions::new()
            .append(true)
            .open(&output.path)
            .map_err(|e| e.to_string())?;
        file.write_all(&output.buffer).map_err(|e| e.to_string())?;
        output.buffer.clear();
        Ok(())
    }

    fn flush(&mut self) -> Result<(), String> {
        for locus in 0..self.outputs.len() {
            self.flush_one(locus)?;
        }
        Ok(())
    }
}

fn write_record(out: &mut impl Write, record: &FastxRecord) -> Result<(), String> {
    out.write_all(&record.header).map_err(|e| e.to_string())?;
    out.write_all(b"\n").map_err(|e| e.to_string())?;
    out.write_all(&record.sequence).map_err(|e| e.to_string())?;
    out.write_all(b"\n").map_err(|e| e.to_string())?;
    out.write_all(if record.plus.is_empty() {
        b"+\n"
    } else {
        &record.plus
    })
    .map_err(|e| e.to_string())?;
    if !record.plus.is_empty() {
        out.write_all(b"\n").map_err(|e| e.to_string())?;
    }
    out.write_all(&record.quality).map_err(|e| e.to_string())?;
    out.write_all(b"\n").map_err(|e| e.to_string())?;
    Ok(())
}

type DecodeMessage = Result<Option<Vec<u8>>, String>;

/// Pull-based view over one bounded background decompressor. Two reusable
/// chunks allow inflate to run ahead while keeping memory bounded.
struct BackgroundReader {
    filled: Option<Receiver<DecodeMessage>>,
    empty: Option<SyncSender<Vec<u8>>>,
    current: Vec<u8>,
    offset: usize,
    finished: bool,
    worker: Option<JoinHandle<()>>,
}

impl BackgroundReader {
    fn open(path: &Path, mate: &str) -> Result<Self, String> {
        let (filled_tx, filled_rx) = mpsc::sync_channel(DECODE_BUFFERS_PER_MATE);
        let (empty_tx, empty_rx) = mpsc::sync_channel(DECODE_BUFFERS_PER_MATE);
        for _ in 0..DECODE_BUFFERS_PER_MATE {
            empty_tx
                .send(vec![0_u8; DECODE_CHUNK_BYTES])
                .map_err(|error| error.to_string())?;
        }
        let path = path.to_path_buf();
        let worker = thread::Builder::new()
            .name(format!("uce-decode-{mate}"))
            .spawn(move || {
                gm2_tools::platform::prefer_performance_cores();
                let mut input = match open_input(&path) {
                    Ok(input) => input,
                    Err(error) => {
                        let _ = filled_tx.send(Err(error.to_string()));
                        return;
                    }
                };
                while let Ok(mut buffer) = empty_rx.recv() {
                    buffer.resize(DECODE_CHUNK_BYTES, 0);
                    match input.read(&mut buffer) {
                        Ok(0) => {
                            let _ = filled_tx.send(Ok(None));
                            return;
                        }
                        Ok(length) => {
                            buffer.truncate(length);
                            if filled_tx.send(Ok(Some(buffer))).is_err() {
                                return;
                            }
                        }
                        Err(error) => {
                            let _ = filled_tx.send(Err(error.to_string()));
                            return;
                        }
                    }
                }
            })
            .map_err(|error| error.to_string())?;
        Ok(Self {
            filled: Some(filled_rx),
            empty: Some(empty_tx),
            current: Vec::new(),
            offset: 0,
            finished: false,
            worker: Some(worker),
        })
    }

    fn next_chunk(&mut self) -> std::io::Result<bool> {
        if !self.current.is_empty() {
            let mut buffer = std::mem::take(&mut self.current);
            buffer.clear();
            if let Some(empty) = &self.empty {
                let _ = empty.send(buffer);
            }
        }
        let message = self
            .filled
            .as_ref()
            .ok_or_else(|| std::io::Error::other("background decoder is closed"))?
            .recv()
            .map_err(|_| std::io::Error::other("background decoder stopped unexpectedly"))?;
        match message {
            Ok(Some(buffer)) => {
                self.current = buffer;
                self.offset = 0;
                Ok(true)
            }
            Ok(None) => {
                self.finished = true;
                Ok(false)
            }
            Err(error) => Err(std::io::Error::new(std::io::ErrorKind::InvalidData, error)),
        }
    }
}

impl Read for BackgroundReader {
    fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
        if output.is_empty() || self.finished {
            return Ok(0);
        }
        if self.offset == self.current.len() && !self.next_chunk()? {
            return Ok(0);
        }
        let length = output.len().min(self.current.len() - self.offset);
        output[..length].copy_from_slice(&self.current[self.offset..self.offset + length]);
        self.offset += length;
        Ok(length)
    }
}

impl Drop for BackgroundReader {
    fn drop(&mut self) {
        self.filled.take();
        self.empty.take();
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}

/// UCE-only FASTQ reader. Unlike the shared general-purpose reader, it fills
/// caller-owned records so non-candidate fragments retain their allocations.
struct FastqScratchReader {
    input: BufReader<BackgroundReader>,
}

impl FastqScratchReader {
    fn open(path: &Path, mate: &str) -> Result<Self, String> {
        Ok(Self {
            input: BufReader::with_capacity(
                DECODE_CHUNK_BYTES,
                BackgroundReader::open(path, mate)?,
            ),
        })
    }

    fn read_line_into(&mut self, target: &mut Vec<u8>) -> Result<bool, String> {
        target.clear();
        if self
            .input
            .read_until(b'\n', target)
            .map_err(|e| e.to_string())?
            == 0
        {
            return Ok(false);
        }
        while matches!(target.last(), Some(b'\n' | b'\r')) {
            target.pop();
        }
        Ok(true)
    }

    fn next_record_into(&mut self, record: &mut FastxRecord) -> Result<bool, String> {
        record.header.clear();
        record.sequence.clear();
        record.plus.clear();
        record.quality.clear();
        loop {
            if !self.read_line_into(&mut record.header)? {
                return Ok(false);
            }
            if !record.header.is_empty() {
                break;
            }
        }
        if record.header.first() != Some(&b'@') {
            return Err("malformed FASTQ record".to_string());
        }
        if !self.read_line_into(&mut record.sequence)? {
            return Err("truncated FASTQ sequence".to_string());
        }
        if !self.read_line_into(&mut record.plus)? || record.plus.first() != Some(&b'+') {
            return Err("malformed or truncated FASTQ plus line".to_string());
        }
        if !self.read_line_into(&mut record.quality)? {
            return Err("truncated FASTQ quality".to_string());
        }
        if record.quality.len() != record.sequence.len() {
            return Err("FASTQ sequence and quality lengths differ".to_string());
        }
        Ok(true)
    }
}

fn next_pair_into(
    r1: &mut FastqScratchReader,
    r2: &mut FastqScratchReader,
    first: &mut FastxRecord,
    second: &mut FastxRecord,
) -> Result<bool, String> {
    let first_present = r1.next_record_into(first)?;
    let second_present = r2.next_record_into(second)?;
    match (first_present, second_present) {
        (false, false) => Ok(false),
        (true, true) if paired_read_id(&first.header) == paired_read_id(&second.header) => Ok(true),
        (true, true) => Err("paired input files contain mismatched read identifiers".to_string()),
        _ => Err("paired input files contain different numbers of records".to_string()),
    }
}

/// Normalizes the identifier before whitespace and an optional /1 or /2
/// suffix, covering both common FASTQ paired-end header conventions.
fn paired_read_id(header: &[u8]) -> &[u8] {
    let header = header.strip_prefix(b"@").unwrap_or(header);
    let token_end = header
        .iter()
        .position(|base| base.is_ascii_whitespace())
        .unwrap_or(header.len());
    let token = &header[..token_end];
    token
        .strip_suffix(b"/1")
        .or_else(|| token.strip_suffix(b"/2"))
        .unwrap_or(token)
}

#[inline]
fn keep_linked_pair(orient1: u8, orient2: u8) -> bool {
    !((1..=2).contains(&orient1) && orient1 == orient2 || orient1 == 0 && orient2 == 0)
}

#[derive(Clone, Copy, Debug, Default)]
struct FastEvidence {
    covered_bins: u64,
    max_exact: u16,
    left_extension: u16,
    right_extension: u16,
    terminal_mask: u8,
    aligned_mates: u8,
}

#[derive(Clone, Copy, Debug)]
struct MateEvidence {
    best: Option<ExactSeed>,
    orientation: u8,
}

fn add_exact_evidence(
    index: &UceIndex,
    seed: ExactSeed,
    query_length: usize,
    terminal_window: usize,
    evidence: &mut FastEvidence,
) {
    evidence.max_exact = evidence
        .max_exact
        .max(seed.len().min(u16::MAX as usize) as u16);
    evidence.aligned_mates = evidence.aligned_mates.saturating_add(1);
    let reference = &index.references[seed.sequence as usize];
    let reference_length = reference.bases.len().max(1);
    let (original_start, original_end, left_overhang, right_overhang) = if reference.strand == 1 {
        (
            seed.reference_start as usize,
            seed.reference_end as usize,
            seed.read_start as usize,
            query_length.saturating_sub(seed.read_end as usize),
        )
    } else {
        (
            reference_length.saturating_sub(seed.reference_end as usize),
            reference_length.saturating_sub(seed.reference_start as usize),
            query_length.saturating_sub(seed.read_end as usize),
            seed.read_start as usize,
        )
    };
    let first_bin = (original_start * 64 / reference_length).min(63);
    let last_position = original_end.saturating_sub(1);
    let last_bin = (last_position * 64 / reference_length).min(63);
    for bin in first_bin..=last_bin {
        evidence.covered_bins |= 1_u64 << bin;
    }
    let effective_window = terminal_window.min(reference_length / 5);
    if original_start <= effective_window && left_overhang > 0 {
        evidence.terminal_mask |= 1;
        evidence.left_extension = evidence
            .left_extension
            .max(left_overhang.min(u16::MAX as usize) as u16);
    }
    if reference_length.saturating_sub(original_end) <= effective_window && right_overhang > 0 {
        evidence.terminal_mask |= 2;
        evidence.right_extension = evidence
            .right_extension
            .max(right_overhang.min(u16::MAX as usize) as u16);
    }
}

fn constrain_evidence(
    evidence: &mut Vec<(LocusId, FastEvidence)>,
    max_locus_count: usize,
    retain_loci: Option<&[bool]>,
) {
    if max_locus_count > 0 && evidence.len() > max_locus_count {
        evidence.clear();
    } else if let Some(retain_loci) = retain_loci {
        evidence.retain(|(locus, _)| retain_loci[*locus as usize]);
    }
}

pub fn run(config: &Config) -> Result<RunSummary, String> {
    gm2_tools::platform::prefer_performance_cores();
    let started = Instant::now();
    let index_started = Instant::now();
    let verify_k = config
        .verify_kmer_size
        .unwrap_or_else(|| default_verify_k(config.kmer_size));
    let index = UceIndex::build_split_with_verify_k(
        config
            .recruit_references
            .as_deref()
            .unwrap_or(&config.references),
        &config.references,
        config.kmer_size,
        verify_k,
    )?;
    let index_seconds = index_started.elapsed().as_secs_f64();
    let retain_loci = config
        .retain_loci_file
        .as_deref()
        .map(|path| retain_locus_mask(path, &index))
        .transpose()?;
    eprintln!(
        "UCEFilter index: {} loci, {} FM-index symbols, k={}, run-k={} ({:.3}s)",
        index.loci.len(),
        index.exact_index_symbols(),
        index.k,
        index.run_k,
        index_seconds
    );
    if index.requires_panel_recruitment_expansion() {
        eprintln!("UCEFilter recruitment: subset gate with full-panel candidate expansion");
    }
    eprintln!("UCEFilter gzip backend: {}", gzip_backend_name());
    eprintln!(
        "UCEFilter CPU: {}",
        gm2_tools::platform::cpu_info().summary()
    );
    if config.profile {
        eprintln!(
            "UCEFilter paired decode: 2 background workers, {} x {} KiB buffers per mate",
            DECODE_BUFFERS_PER_MATE,
            DECODE_CHUNK_BYTES / 1024,
        );
    }
    let mut reader1 = FastqScratchReader::open(&config.read1, "r1")?;
    let mut reader2 = FastqScratchReader::open(&config.read2, "r2")?;
    // Parse one complete paired record before replacing a prior result. This
    // catches malformed leading FASTQ records and paired-file mix-ups without
    // losing that fragment from the actual filtering pass.
    let initial_decode_started = config.profile.then(Instant::now);
    let mut r1 = FastxRecord::default();
    let mut r2 = FastxRecord::default();
    let mut pending_pair = next_pair_into(&mut reader1, &mut reader2, &mut r1, &mut r2)?;
    let initial_decode_seconds =
        initial_decode_started.map_or(0.0, |started| started.elapsed().as_secs_f64());
    // Do not replace a previous result until the references and both inputs
    // have passed their initial open/index checks.
    fs::create_dir_all(&config.output).map_err(|e| e.to_string())?;
    let filtered = config.output.join("filtered");
    if filtered.exists() {
        fs::remove_dir_all(&filtered).map_err(|e| e.to_string())?;
    }
    fs::create_dir_all(&filtered).map_err(|e| e.to_string())?;
    let memory_limit = config.memory_limit_mib.saturating_mul(1024 * 1024);
    let mut bank = FragmentBank::new(memory_limit, default_spill_path(&config.output));
    let mut locus_candidates = LocusCandidateStore::new(index.loci.len());
    let mut coarse_counts = vec![0_u64; index.loci.len()];
    let mut ordinal = 0_u64;
    let mut recruited_loci = RecruitScratch::default();
    let mut read_evidence = ReadEvidenceScratch::default();
    let mut mate1_evidence = Vec::<MateEvidence>::new();
    let mut evidence = Vec::<(LocusId, FastEvidence)>::new();
    let mut decode_seconds = initial_decode_seconds;
    let mut recruit_seconds = 0.0_f64;
    let mut evidence_seconds = 0.0_f64;
    let mut store_seconds = 0.0_f64;
    let mut index_profile = IndexProfile::default();
    let scan_started = Instant::now();
    loop {
        if config.max_fragments > 0 && ordinal >= config.max_fragments {
            break;
        }
        let stage_started = config.profile.then(Instant::now);
        let pair_available = if pending_pair {
            pending_pair = false;
            true
        } else {
            next_pair_into(&mut reader1, &mut reader2, &mut r1, &mut r2)?
        };
        if let Some(started) = stage_started {
            decode_seconds += started.elapsed().as_secs_f64();
        }
        if !pair_available {
            break;
        }
        let stage_started = config.profile.then(Instant::now);
        recruited_loci.begin(index.loci.len());
        index.recruit(
            &r1.sequence,
            config.step,
            &mut recruited_loci,
            config.profile.then_some(&mut index_profile),
        );
        index.recruit(
            &r2.sequence,
            config.step,
            &mut recruited_loci,
            config.profile.then_some(&mut index_profile),
        );
        recruited_loci.sort();
        if recruited_loci.loci().is_empty() {
            if let Some(started) = stage_started {
                recruit_seconds += started.elapsed().as_secs_f64();
            }
            ordinal += 1;
            continue;
        }
        if index.requires_panel_recruitment_expansion() {
            recruited_loci.begin(index.loci.len());
            index.expand_recruitment_to_panel(
                &r1.sequence,
                config.step,
                &mut recruited_loci,
                config.profile.then_some(&mut index_profile),
            );
            index.expand_recruitment_to_panel(
                &r2.sequence,
                config.step,
                &mut recruited_loci,
                config.profile.then_some(&mut index_profile),
            );
            recruited_loci.sort();
        }
        let loci = recruited_loci.loci();
        for &locus in loci {
            coarse_counts[locus as usize] += 2;
        }
        if let Some(started) = stage_started {
            recruit_seconds += started.elapsed().as_secs_f64();
        }
        let stage_started = config.profile.then(Instant::now);
        index.read_evidence(
            &r1.sequence,
            loci,
            &mut read_evidence,
            config.profile.then_some(&mut index_profile),
        );
        mate1_evidence.clear();
        for i in 0..loci.len() {
            mate1_evidence.push(MateEvidence {
                best: read_evidence.best(i),
                orientation: infer_orientation(collect_runs_stats(read_evidence.orientation(i))),
            });
        }
        index.read_evidence(
            &r2.sequence,
            loci,
            &mut read_evidence,
            config.profile.then_some(&mut index_profile),
        );
        evidence.clear();
        for (i, &locus) in loci.iter().enumerate() {
            let mate1 = mate1_evidence[i];
            let orient2 = infer_orientation(collect_runs_stats(read_evidence.orientation(i)));
            if !keep_linked_pair(mate1.orientation, orient2) {
                continue;
            }
            if !passes_alignment_filter(
                &index,
                locus,
                &r1.sequence,
                &r2.sequence,
                config.shadow_band,
                config.minimum_alignment_overlap,
                config.minimum_alignment_identity,
            ) {
                continue;
            }
            let mut fast = FastEvidence::default();
            if let Some(seed) = mate1.best {
                add_exact_evidence(
                    &index,
                    seed,
                    r1.sequence.len(),
                    config.terminal_window,
                    &mut fast,
                );
            }
            if let Some(seed) = read_evidence.best(i) {
                add_exact_evidence(
                    &index,
                    seed,
                    r2.sequence.len(),
                    config.terminal_window,
                    &mut fast,
                );
            }
            evidence.push((locus, fast));
        }
        constrain_evidence(
            &mut evidence,
            config.max_locus_count,
            retain_loci.as_deref(),
        );
        if let Some(started) = stage_started {
            evidence_seconds += started.elapsed().as_secs_f64();
        }
        let stage_started = config.profile.then(Instant::now);
        if !evidence.is_empty() {
            let fragment_bases = (r1.sequence.len() + r2.sequence.len()) as u32;
            let fragment_id = bank.insert(Fragment {
                ordinal,
                r1: std::mem::take(&mut r1),
                r2: std::mem::take(&mut r2),
            })?;
            let locus_count = evidence.len().min(u16::MAX as usize) as u16;
            for &(locus, fast) in &evidence {
                locus_candidates.push(
                    locus,
                    Candidate {
                        fragment_id,
                        fragment_bases,
                        covered_bins: fast.covered_bins,
                        max_exact: fast.max_exact,
                        left_extension: fast.left_extension,
                        right_extension: fast.right_extension,
                        locus_count,
                        terminal_mask: fast.terminal_mask,
                        aligned_mates: fast.aligned_mates,
                    },
                )?;
            }
        }
        if let Some(started) = stage_started {
            store_seconds += started.elapsed().as_secs_f64();
        }
        ordinal += 1;
        if ordinal.is_multiple_of(1_048_576) {
            eprintln!("UCEFilter handled {} Mi fragments", ordinal / 1_048_576);
        }
    }
    let scan_seconds = scan_started.elapsed().as_secs_f64();
    let evidence_scratch_bytes = read_evidence.allocated_bytes()
        + mate1_evidence.capacity() * std::mem::size_of::<MateEvidence>()
        + evidence.capacity() * std::mem::size_of::<(LocusId, FastEvidence)>();
    let candidate_memory_bytes = locus_candidates.allocated_bytes();
    let selection_started = Instant::now();
    let mut loci_written = 0_usize;
    let mut assignments = 0_usize;
    let mut route_pairs = Vec::<(u32, LocusId)>::new();
    let mut shadow_routes = config.alignment_shadow.then(|| {
        (0..bank.len())
            .map(|_| Vec::new())
            .collect::<Vec<Vec<LocusId>>>()
    });
    let mut router = OutputRouter::new(index.loci.len());
    let mut summary = BufWriter::new(
        File::create(config.output.join("uce_filter_summary.tsv")).map_err(|e| e.to_string())?,
    );
    if config.selection_auto {
        writeln!(
            summary,
            "locus\tcoarse_reads\trun_fragments\tselected_fragments\tminimum_exact\tthinning_interval\tselection_mode\teligible_fragments\tcovered_bins_64\ttarget_core_fragments\tleft_terminal_candidates\tright_terminal_candidates\tselected_left_terminal\tselected_right_terminal"
        )
        .map_err(|e| e.to_string())?;
    } else {
        writeln!(
            summary,
            "locus\tcoarse_reads\trun_fragments\tselected_fragments\tminimum_exact\tthinning_interval"
        )
        .map_err(|e| e.to_string())?;
    }
    let mut candidates = Vec::<Candidate>::new();
    for (locus_id, locus) in index.loci.iter().enumerate() {
        locus_candidates.copy_locus(locus_id, &mut candidates);
        let decision = choose_legacy(
            &candidates,
            locus.effective_length,
            config.kmer_size,
            config.min_depth,
            config.max_depth,
            config.max_size_mb,
        );
        let (selected_ids, auto_details) = if config.selection_auto {
            let automatic = choose_auto(
                &candidates,
                &decision,
                locus.effective_length,
                config.reference_is_contig,
            );
            let details = (
                automatic.mode,
                automatic.eligible_fragments,
                automatic.covered_bins,
                automatic.target_core_fragments,
                automatic.left_terminal_candidates,
                automatic.right_terminal_candidates,
                automatic.selected_left_terminal,
                automatic.selected_right_terminal,
            );
            (automatic.selected_ids, Some(details))
        } else {
            let mut passing = 0_usize;
            let ids = candidates
                .iter()
                .filter(|candidate| selected(candidate, &decision, &mut passing))
                .map(|candidate| candidate.fragment_id)
                .collect();
            (ids, None)
        };
        let minimum_exact = decision
            .minimum_exact
            .map_or_else(|| "all".to_string(), |v| v.to_string());
        let thinning_interval = decision
            .thinning_interval
            .map_or_else(|| "none".to_string(), |v| v.to_string());
        if let Some((
            mode,
            eligible,
            covered_bins,
            target_core,
            left_terminal_candidates,
            right_terminal_candidates,
            selected_left_terminal,
            selected_right_terminal,
        )) = auto_details
        {
            writeln!(
                summary,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                locus.name,
                coarse_counts[locus_id],
                candidates.len(),
                selected_ids.len(),
                minimum_exact,
                thinning_interval,
                mode.as_str(),
                eligible,
                covered_bins,
                target_core,
                left_terminal_candidates,
                right_terminal_candidates,
                selected_left_terminal,
                selected_right_terminal,
            )
            .map_err(|e| e.to_string())?;
        } else {
            writeln!(
                summary,
                "{}\t{}\t{}\t{}\t{}\t{}",
                locus.name,
                coarse_counts[locus_id],
                candidates.len(),
                selected_ids.len(),
                minimum_exact,
                thinning_interval,
            )
            .map_err(|e| e.to_string())?;
        }
        if selected_ids.is_empty() {
            continue;
        }
        loci_written += 1;
        assignments += selected_ids.len();
        router.register(locus_id, filtered.join(format!("{}.fq", locus.name)))?;
        if let Some(shadow_routes) = shadow_routes.as_mut() {
            for id in evenly_sample(&selected_ids, config.shadow_per_locus) {
                shadow_routes[id as usize].push(locus_id as LocusId);
            }
        }
        for id in selected_ids {
            route_pairs.push((id, locus_id as LocusId));
        }
    }
    drop(candidates);
    drop(locus_candidates);
    summary.flush().map_err(|e| e.to_string())?;
    let routes = FragmentRoutes::from_pairs(bank.len(), &route_pairs)?;
    drop(route_pairs);
    let selection_seconds = selection_started.elapsed().as_secs_f64();
    let fragment_memory_bytes = bank.memory_bytes();
    let fragment_spill_bytes = bank.spill_bytes();
    let mut shadow_writer = if config.alignment_shadow {
        Some(ShadowWriter::new(
            &config.output.join("alignment_shadow.tsv"),
            config.shadow_band,
            config.terminal_window,
            index.loci.len(),
        )?)
    } else {
        None
    };
    let output_started = Instant::now();
    bank.stream_in_order(|id, fragment| {
        for &locus in routes.get(id) {
            router.write_fragment(locus as usize, fragment)?;
        }
        if let (Some(shadow_routes), Some(shadow_writer)) =
            (shadow_routes.as_ref(), shadow_writer.as_mut())
        {
            for &locus in &shadow_routes[id as usize] {
                shadow_writer.write_fragment(
                    &index,
                    locus,
                    &index.loci[locus as usize].name,
                    fragment,
                    routes.get(id).len(),
                )?;
            }
        }
        Ok(())
    })?;
    router.flush()?;
    if let Some(shadow_writer) = shadow_writer.as_mut() {
        shadow_writer.flush()?;
        shadow_writer.write_summary(&config.output.join("alignment_shadow_summary.tsv"), &index)?;
    }
    let output_seconds = output_started.elapsed().as_secs_f64();
    let mut counts = BufWriter::new(
        File::create(config.output.join("ref_reads_count_dict.txt")).map_err(|e| e.to_string())?,
    );
    for (locus, count) in index.loci.iter().zip(coarse_counts) {
        if count > 0 {
            writeln!(counts, "{},{}", locus.name, count).map_err(|e| e.to_string())?;
        }
    }
    counts.flush().map_err(|e| e.to_string())?;
    let (shadow_sampled_assignments, shadow_aligned_mates, shadow_seconds) =
        shadow_writer.map_or((0, 0, 0.0), |writer| {
            (
                writer.sampled_assignments,
                writer.aligned_mates,
                writer.elapsed_seconds,
            )
        });
    Ok(RunSummary {
        fragments_read: ordinal,
        fragments_retained_once: bank.len(),
        fragment_bases_retained_once: bank.bases(),
        assignments,
        loci_written,
        fragment_memory_bytes,
        fragment_spill_bytes,
        evidence_scratch_bytes: evidence_scratch_bytes as u64,
        candidate_memory_bytes: candidate_memory_bytes as u64,
        shadow_sampled_assignments,
        shadow_aligned_mates,
        shadow_seconds,
        index_seconds,
        scan_seconds,
        selection_seconds,
        output_seconds,
        decode_seconds,
        recruit_seconds,
        evidence_seconds,
        store_seconds,
        index_profile,
        elapsed_seconds: started.elapsed().as_secs_f64(),
    })
}

pub fn output_exists(path: &Path) -> bool {
    path.join("filtered").is_dir() && path.join("ref_reads_count_dict.txt").is_file()
}

#[cfg(test)]
mod tests {
    use super::{
        add_exact_evidence, constrain_evidence, evenly_sample, keep_linked_pair, next_pair_into,
        paired_read_id, run, BackgroundReader, Config, FastEvidence, FastqScratchReader,
        FragmentRoutes, LocusCandidateStore, DECODE_CHUNK_BYTES,
    };
    use crate::index::{reverse_complement, UceIndex};
    use crate::model::Candidate;
    use gm2_tools::fastx::FastxRecord;
    use std::fs;
    use std::io::{Read, Write};

    #[test]
    fn paired_end_decision_preserves_legacy_linked_mate_semantics() {
        assert!(!keep_linked_pair(0, 0));
        assert!(keep_linked_pair(1, 0));
        assert!(keep_linked_pair(0, 2));
        assert!(!keep_linked_pair(1, 1));
        assert!(!keep_linked_pair(2, 2));
        assert!(keep_linked_pair(1, 2));
        assert!(keep_linked_pair(3, 0));
        assert!(keep_linked_pair(3, 3));
    }

    #[test]
    fn panel_ambiguity_is_rejected_before_the_retain_locus_gate() {
        let mut evidence = vec![(0, FastEvidence::default()), (1, FastEvidence::default())];
        constrain_evidence(&mut evidence, 1, Some(&[true, false]));
        assert!(evidence.is_empty());

        let mut evidence = vec![(0, FastEvidence::default())];
        constrain_evidence(&mut evidence, 1, Some(&[true, false]));
        assert_eq!(evidence.len(), 1);

        let mut evidence = vec![(1, FastEvidence::default())];
        constrain_evidence(&mut evidence, 1, Some(&[true, false]));
        assert!(evidence.is_empty());
    }

    #[test]
    fn subset_recruitment_expands_to_the_panel_before_ambiguity_rejection() {
        let root = std::env::temp_dir().join(format!(
            "uce-filter-subset-panel-ambiguity-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&root);
        let references = root.join("references");
        let recruit_references = root.join("recruit-references");
        fs::create_dir_all(&references).unwrap();
        fs::create_dir_all(&recruit_references).unwrap();
        let reference =
            b"ACGTTGCAACGATTCGGTACCATGCAAGTTCGATCGGATCCGTAACCGGTTAGCTACGATGCTAGGCTTACCGAT";
        let fasta = format!(">ref\n{}\n", String::from_utf8_lossy(reference));
        fs::write(references.join("target.fa"), &fasta).unwrap();
        fs::write(references.join("other.fa"), &fasta).unwrap();
        fs::write(recruit_references.join("target.fa"), &fasta).unwrap();
        let r1 = &reference[..25];
        let r2 = reverse_complement(&reference[35..60]);
        let fastq = |mate: usize, sequence: &[u8]| {
            format!(
                "@pair/{}\n{}\n+\n{}\n",
                mate,
                String::from_utf8_lossy(sequence),
                "I".repeat(sequence.len())
            )
        };
        let read1 = root.join("r1.fq");
        let read2 = root.join("r2.fq");
        fs::write(&read1, fastq(1, r1)).unwrap();
        fs::write(&read2, fastq(2, &r2)).unwrap();
        let retain = root.join("retain.txt");
        fs::write(&retain, "target\n").unwrap();
        let output = root.join("output");
        let summary = run(&Config {
            references,
            recruit_references: Some(recruit_references),
            read1,
            read2,
            output: output.clone(),
            kmer_size: 21,
            verify_kmer_size: Some(19),
            step: 1,
            max_locus_count: 1,
            retain_loci_file: Some(retain),
            minimum_alignment_overlap: 0,
            minimum_alignment_identity: 0.0,
            min_depth: 50,
            max_depth: 768,
            max_size_mb: 6,
            max_fragments: 0,
            memory_limit_mib: 16,
            selection_auto: true,
            reference_is_contig: false,
            alignment_shadow: false,
            shadow_per_locus: 64,
            shadow_band: 32,
            terminal_window: 150,
            profile: false,
        })
        .unwrap();
        assert_eq!(summary.loci_written, 0);
        assert!(!output.join("filtered/target.fq").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn sensitive_pass_recovers_a_pair_that_fast_k31_cannot_recruit() {
        let root =
            std::env::temp_dir().join(format!("uce-filter-sensitive-pass-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        let references = root.join("references");
        fs::create_dir_all(&references).unwrap();
        let reference =
            b"ACGTTGCAACGATTCGGTACCATGCAAGTTCGATCGGATCCGTAACCGGTTAGCTACGATGCTAGGCTTACCGAT";
        fs::write(
            references.join("locus.fa"),
            format!(">ref\n{}\n", String::from_utf8_lossy(reference)),
        )
        .unwrap();
        let r1 = &reference[..25];
        let r2 = reverse_complement(&reference[35..60]);
        let fastq = |mate: usize, sequence: &[u8]| {
            format!(
                "@pair/{}\n{}\n+\n{}\n",
                mate,
                String::from_utf8_lossy(sequence),
                "I".repeat(sequence.len())
            )
        };
        let read1 = root.join("r1.fq");
        let read2 = root.join("r2.fq");
        fs::write(&read1, fastq(1, r1)).unwrap();
        fs::write(&read2, fastq(2, &r2)).unwrap();
        let mut config = Config {
            references: references.clone(),
            recruit_references: None,
            read1,
            read2,
            output: root.join("fast"),
            kmer_size: 31,
            verify_kmer_size: None,
            step: 4,
            max_locus_count: 0,
            retain_loci_file: None,
            minimum_alignment_overlap: 0,
            minimum_alignment_identity: 0.0,
            min_depth: 50,
            max_depth: 768,
            max_size_mb: 6,
            max_fragments: 0,
            memory_limit_mib: 16,
            selection_auto: true,
            reference_is_contig: false,
            alignment_shadow: false,
            shadow_per_locus: 64,
            shadow_band: 32,
            terminal_window: 150,
            profile: false,
        };
        let fast = run(&config).unwrap();
        assert_eq!(fast.loci_written, 0);
        assert!(!config.output.join("filtered/locus.fq").exists());

        let retain = root.join("retain.txt");
        fs::write(&retain, "locus\n").unwrap();
        let recruit_references = root.join("recruit-references");
        fs::create_dir_all(&recruit_references).unwrap();
        fs::copy(
            references.join("locus.fa"),
            recruit_references.join("locus.fa"),
        )
        .unwrap();
        config.output = root.join("sensitive");
        config.recruit_references = Some(recruit_references);
        config.kmer_size = 21;
        config.verify_kmer_size = Some(19);
        config.step = 1;
        config.max_locus_count = 1;
        config.retain_loci_file = Some(retain);
        let sensitive = run(&config).unwrap();
        assert_eq!(sensitive.loci_written, 1);
        assert!(config.output.join("filtered/locus.fq").is_file());

        config.output = root.join("alignment-filtered");
        config.minimum_alignment_overlap = 26;
        let alignment_filtered = run(&config).unwrap();
        assert_eq!(alignment_filtered.loci_written, 0);
        assert!(!config.output.join("filtered/locus.fq").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn paired_read_ids_normalize_common_fastq_conventions() {
        assert_eq!(paired_read_id(b"@read-42/1"), b"read-42");
        assert_eq!(paired_read_id(b"@read-42/2"), b"read-42");
        assert_eq!(paired_read_id(b"@read-42 1:N:0:ACGT"), b"read-42");
    }

    #[test]
    fn background_reader_preserves_recycled_chunks_and_eof() {
        let root = std::env::temp_dir().join(format!(
            "uce-filter-background-reader-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let path = root.join("input.fastq");
        let expected = (0..DECODE_CHUNK_BYTES * 3 + 17)
            .map(|index| (index % 251) as u8)
            .collect::<Vec<_>>();
        fs::write(&path, &expected).unwrap();
        let mut reader = BackgroundReader::open(&path, "test-bytes").unwrap();
        let mut observed = Vec::new();
        reader.read_to_end(&mut observed).unwrap();
        assert_eq!(observed, expected);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn background_reader_propagates_open_errors() {
        let path = std::env::temp_dir().join(format!(
            "uce-filter-missing-background-input-{}",
            std::process::id()
        ));
        let mut reader = BackgroundReader::open(&path, "test-error").unwrap();
        let mut observed = Vec::new();
        assert!(reader.read_to_end(&mut observed).is_err());
    }

    #[test]
    fn mismatched_paired_fastq_identifiers_fail() {
        let root = std::env::temp_dir().join(format!("uce-filter-pair-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let r1_path = root.join("r1.fq");
        let r2_path = root.join("r2.fq");
        fs::write(&r1_path, b"@read-a/1\nACGT\n+\n!!!!\n").unwrap();
        fs::write(&r2_path, b"@read-b/2\nACGT\n+\n!!!!\n").unwrap();
        let mut reader1 = FastqScratchReader::open(&r1_path, "test-r1").unwrap();
        let mut reader2 = FastqScratchReader::open(&r2_path, "test-r2").unwrap();
        let mut r1 = FastxRecord::default();
        let mut r2 = FastxRecord::default();
        assert!(next_pair_into(&mut reader1, &mut reader2, &mut r1, &mut r2).is_err());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn bounded_shadow_sampling_is_deterministic_and_spans_input_order() {
        let ids: Vec<u32> = (0..100).collect();
        assert_eq!(evenly_sample(&ids, 4), vec![0, 33, 66, 99]);
        assert_eq!(evenly_sample(&ids[..3], 4), vec![0, 1, 2]);
        assert_eq!(evenly_sample(&ids, 1), vec![50]);
    }

    #[test]
    fn compact_routes_preserve_fragment_and_locus_order() {
        let routes = FragmentRoutes::from_pairs(4, &[(2, 7), (0, 3), (2, 9), (3, 1)]).unwrap();
        assert_eq!(routes.get(0), &[3]);
        assert!(routes.get(1).is_empty());
        assert_eq!(routes.get(2), &[7, 9]);
        assert_eq!(routes.get(3), &[1]);
    }

    #[test]
    fn compact_candidate_store_preserves_per_locus_input_order() {
        let candidate = |fragment_id| Candidate {
            fragment_id,
            fragment_bases: 300,
            covered_bins: 1,
            max_exact: 31,
            left_extension: 0,
            right_extension: 0,
            locus_count: 1,
            terminal_mask: 0,
            aligned_mates: 1,
        };
        let mut store = LocusCandidateStore::new(3);
        store.push(2, candidate(4)).unwrap();
        store.push(0, candidate(5)).unwrap();
        store.push(2, candidate(8)).unwrap();
        store.push(1, candidate(9)).unwrap();
        let mut scratch = Vec::new();
        store.copy_locus(2, &mut scratch);
        assert_eq!(
            scratch
                .iter()
                .map(|value| value.fragment_id)
                .collect::<Vec<_>>(),
            vec![4, 8]
        );
        store.copy_locus(0, &mut scratch);
        assert_eq!(
            scratch
                .iter()
                .map(|value| value.fragment_id)
                .collect::<Vec<_>>(),
            vec![5]
        );
    }

    #[test]
    fn evidence_layout_remains_compact() {
        assert!(std::mem::size_of::<FastEvidence>() <= 16);
    }

    #[test]
    fn candidate_arena_limits_unused_capacity() {
        let candidate = Candidate {
            fragment_id: 0,
            fragment_bases: 300,
            covered_bins: 1,
            max_exact: 31,
            left_extension: 0,
            right_extension: 0,
            locus_count: 1,
            terminal_mask: 0,
            aligned_mates: 1,
        };
        let mut store = LocusCandidateStore::new(1);
        for _ in 0..70_000 {
            store.push(0, candidate).unwrap();
        }
        assert!(store.allocated_bytes() < 3 * 1024 * 1024);
    }

    #[test]
    fn exact_evidence_normalizes_reverse_coordinates_and_terminal_sides() {
        let root = std::env::temp_dir().join(format!("uce-filter-pipeline-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let sequence = b"ACGTTGCAACGATTCGGTACCATGCAAGTTCGATCGGATCCGTAACCGGTT";
        let mut out = fs::File::create(root.join("locus.fa")).unwrap();
        writeln!(out, ">ref\n{}", String::from_utf8_lossy(sequence)).unwrap();
        drop(out);
        let index = UceIndex::build(&root, 16).unwrap();
        let read = crate::index::reverse_complement(&sequence[..32]);
        let seed = index.best_exact(&read, 0).unwrap();
        let mut evidence = FastEvidence::default();
        add_exact_evidence(&index, seed, read.len() + 5, 150, &mut evidence);
        assert_ne!(evidence.covered_bins, 0);
        assert_eq!(evidence.aligned_mates, 1);
        assert_ne!(evidence.terminal_mask, 0);
        fs::remove_dir_all(root).unwrap();
    }
}
