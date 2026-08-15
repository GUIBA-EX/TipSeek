use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use uce_filter_core::alignment::align_read;
use uce_filter_core::index::{RecruitScratch, UceIndex};

pub(crate) const FALLBACK_CONTIG_MIN_PROBE_COVERAGE: f64 = 0.80;
pub(crate) const FALLBACK_CONTIG_MIN_PROBE_IDENTITY: f64 = 0.80;
pub(crate) const FALLBACK_CONTIG_MIN_LENGTH: usize = 200;
const FALLBACK_CONTIG_OTHER_SCORE_RATIO: f64 = 0.95;
const FALLBACK_CONTIG_ALIGNMENT_BAND: usize = 64;

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ContigProbeEvidence {
    pub(crate) locus: String,
    pub(crate) sequence_length: usize,
    pub(crate) accepted: bool,
    pub(crate) reason: &'static str,
    pub(crate) core_anchor_status: &'static str,
    pub(crate) structural_review: &'static str,
    pub(crate) long_inverted_repeat: bool,
    pub(crate) maximum_unsupported_internal_gap: usize,
    pub(crate) target_score: i32,
    pub(crate) target_probe_coverage: f64,
    pub(crate) target_identity: f64,
    pub(crate) best_other_locus: String,
    pub(crate) best_other_score: i32,
    pub(crate) near_tie_other_loci: usize,
}

impl ContigProbeEvidence {
    pub(crate) fn unavailable(locus: &str, reason: &'static str) -> Self {
        Self {
            locus: locus.to_owned(),
            sequence_length: 0,
            accepted: false,
            reason,
            core_anchor_status: reason,
            structural_review: "none",
            long_inverted_repeat: false,
            maximum_unsupported_internal_gap: 0,
            target_score: 0,
            target_probe_coverage: 0.0,
            target_identity: 0.0,
            best_other_locus: String::new(),
            best_other_score: 0,
            near_tie_other_loci: 0,
        }
    }

    pub(crate) fn apply_structure_checks(
        &mut self,
        long_inverted_repeat: bool,
        maximum_unsupported_internal_gap: usize,
    ) {
        if !self.accepted {
            return;
        }
        self.long_inverted_repeat = long_inverted_repeat;
        self.maximum_unsupported_internal_gap = maximum_unsupported_internal_gap;
        if long_inverted_repeat {
            self.accepted = false;
            self.reason = "long_inverted_repeat";
            self.core_anchor_status = "structure_rejected";
        } else if maximum_unsupported_internal_gap >= crate::rescue_qc::MINIMUM_NEW_UNSUPPORTED_GAP
        {
            self.core_anchor_status = "anchored_with_review";
            self.structural_review = "internal_gap_ge40";
        } else {
            self.core_anchor_status = "anchored";
        }
    }
}

fn alignment_metrics(index: &UceIndex, sequence: &[u8], locus: u32) -> Option<(i32, f64, f64)> {
    let alignment = align_read(index, sequence, locus, FALLBACK_CONTIG_ALIGNMENT_BAND)?;
    let probe_length = index
        .references
        .get(alignment.sequence as usize)?
        .bases
        .len();
    if probe_length == 0 {
        return None;
    }
    Some((
        alignment.score,
        alignment.reference_overlap() as f64 / probe_length as f64,
        alignment.identity(),
    ))
}

pub(crate) fn evaluate_contig_probe_support_parallel(
    references: &Path,
    contigs: &BTreeMap<String, String>,
    workers: usize,
) -> Result<Vec<ContigProbeEvidence>, String> {
    if contigs.is_empty() {
        return Ok(Vec::new());
    }
    let index = UceIndex::build_split_with_verify_k(references, references, 11, 11)?;
    let locus_ids = index
        .loci
        .iter()
        .enumerate()
        .map(|(id, locus)| (locus.name.as_str(), id as u32))
        .collect::<BTreeMap<_, _>>();
    let contigs = contigs.iter().collect::<Vec<_>>();
    for &(locus, _) in &contigs {
        if !locus_ids.contains_key(locus.as_str()) {
            return Err(format!(
                "UCE fallback probe gate cannot find reference locus '{locus}'"
            ));
        }
    }
    let worker_count = workers.max(1).min(contigs.len());
    if worker_count == 1 {
        return Ok(evaluate_contig_probe_chunk(&index, &locus_ids, &contigs));
    }
    let next_contig = AtomicUsize::new(0);
    thread::scope(|scope| {
        let mut handles = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            handles.push(scope.spawn(|| {
                let mut recruited = RecruitScratch::default();
                let mut rows = Vec::new();
                loop {
                    let row_index = next_contig.fetch_add(1, Ordering::Relaxed);
                    let Some(&(locus, sequence)) = contigs.get(row_index) else {
                        break;
                    };
                    rows.push((
                        row_index,
                        evaluate_contig_probe_row(
                            &index,
                            &locus_ids,
                            locus,
                            sequence,
                            &mut recruited,
                        ),
                    ));
                }
                rows
            }));
        }
        let mut indexed_evidence = Vec::with_capacity(contigs.len());
        for handle in handles {
            let rows = handle
                .join()
                .map_err(|_| "UCE fallback probe gate worker panicked".to_owned())?;
            indexed_evidence.extend(rows);
        }
        indexed_evidence.sort_unstable_by_key(|(row_index, _)| *row_index);
        Ok(indexed_evidence.into_iter().map(|(_, row)| row).collect())
    })
}

fn evaluate_contig_probe_chunk(
    index: &UceIndex,
    locus_ids: &BTreeMap<&str, u32>,
    contigs: &[(&String, &String)],
) -> Vec<ContigProbeEvidence> {
    let mut recruited = RecruitScratch::default();
    contigs
        .iter()
        .map(|&(locus, sequence)| {
            evaluate_contig_probe_row(index, locus_ids, locus, sequence, &mut recruited)
        })
        .collect()
}

fn evaluate_contig_probe_row(
    index: &UceIndex,
    locus_ids: &BTreeMap<&str, u32>,
    locus: &str,
    sequence: &str,
    recruited: &mut RecruitScratch,
) -> ContigProbeEvidence {
    let target_locus = locus_ids[locus];
    if sequence.len() < FALLBACK_CONTIG_MIN_LENGTH {
        return ContigProbeEvidence {
            locus: locus.to_owned(),
            sequence_length: sequence.len(),
            accepted: false,
            reason: "contig_length_below_200",
            core_anchor_status: "probe_rejected",
            structural_review: "none",
            long_inverted_repeat: false,
            maximum_unsupported_internal_gap: 0,
            target_score: 0,
            target_probe_coverage: 0.0,
            target_identity: 0.0,
            best_other_locus: String::new(),
            best_other_score: 0,
            near_tie_other_loci: 0,
        };
    }
    let Some((target_score, target_probe_coverage, target_identity)) =
        alignment_metrics(index, sequence.as_bytes(), target_locus)
    else {
        return ContigProbeEvidence {
            locus: locus.to_owned(),
            sequence_length: sequence.len(),
            accepted: false,
            reason: "no_target_probe_alignment",
            core_anchor_status: "probe_rejected",
            structural_review: "none",
            long_inverted_repeat: false,
            maximum_unsupported_internal_gap: 0,
            target_score: 0,
            target_probe_coverage: 0.0,
            target_identity: 0.0,
            best_other_locus: String::new(),
            best_other_score: 0,
            near_tie_other_loci: 0,
        };
    };
    let target_failure = if target_probe_coverage < FALLBACK_CONTIG_MIN_PROBE_COVERAGE {
        Some("target_probe_coverage_below_0.80")
    } else if target_identity < FALLBACK_CONTIG_MIN_PROBE_IDENTITY {
        Some("target_probe_identity_below_0.80")
    } else {
        None
    };
    if let Some(reason) = target_failure {
        return ContigProbeEvidence {
            locus: locus.to_owned(),
            sequence_length: sequence.len(),
            accepted: false,
            reason,
            core_anchor_status: "probe_rejected",
            structural_review: "none",
            long_inverted_repeat: false,
            maximum_unsupported_internal_gap: 0,
            target_score,
            target_probe_coverage,
            target_identity,
            best_other_locus: String::new(),
            best_other_score: 0,
            near_tie_other_loci: 0,
        };
    }
    recruited.begin(index.loci.len());
    index.recruit(sequence.as_bytes(), 1, recruited, None);
    let mut best_other = None::<(u32, i32)>;
    let mut near_tie_other_loci = 0_usize;
    for &candidate in recruited.loci() {
        if candidate == target_locus {
            continue;
        }
        let Some((score, coverage, identity)) =
            alignment_metrics(index, sequence.as_bytes(), candidate)
        else {
            continue;
        };
        if best_other.is_none_or(|(_, best_score)| score > best_score) {
            best_other = Some((candidate, score));
        }
        if coverage >= FALLBACK_CONTIG_MIN_PROBE_COVERAGE
            && identity >= FALLBACK_CONTIG_MIN_PROBE_IDENTITY
            && score as f64 >= FALLBACK_CONTIG_OTHER_SCORE_RATIO * target_score as f64
        {
            near_tie_other_loci += 1;
        }
    }
    let (accepted, reason) = if near_tie_other_loci > 0 {
        (false, "near_tie_other_locus")
    } else {
        (true, "pass")
    };
    let (best_other_locus, best_other_score) = best_other.map_or_else(
        || (String::new(), 0),
        |(id, score)| (index.loci[id as usize].name.clone(), score),
    );
    ContigProbeEvidence {
        locus: locus.to_owned(),
        sequence_length: sequence.len(),
        accepted,
        reason,
        core_anchor_status: if accepted {
            "probe_pass"
        } else {
            "probe_rejected"
        },
        structural_review: "none",
        long_inverted_repeat: false,
        maximum_unsupported_internal_gap: 0,
        target_score,
        target_probe_coverage,
        target_identity,
        best_other_locus,
        best_other_score,
        near_tie_other_loci,
    }
}

pub(crate) fn fallback_recruited_loci(path: &Path) -> Result<BTreeSet<String>, String> {
    let text = fs::read_to_string(path).map_err(|error| {
        format!(
            "Unable to read UCE recruitment audit '{}': {error}",
            path.display()
        )
    })?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("UCE recruitment audit '{}' is empty", path.display()))?
        .split('\t')
        .collect::<Vec<_>>();
    let locus_column = header
        .iter()
        .position(|name| *name == "locus")
        .ok_or_else(|| {
            format!(
                "UCE recruitment audit '{}' has no locus column",
                path.display()
            )
        })?;
    let pass_column = header
        .iter()
        .position(|name| *name == "final_pass")
        .ok_or_else(|| {
            format!(
                "UCE recruitment audit '{}' has no final_pass column",
                path.display()
            )
        })?;
    let required_columns = locus_column.max(pass_column) + 1;
    let mut loci = BTreeSet::new();
    for (offset, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < required_columns {
            return Err(format!(
                "UCE recruitment audit '{}' row {} has too few columns",
                path.display(),
                offset + 2
            ));
        }
        if fields[pass_column] == "fallback" {
            loci.insert(fields[locus_column].to_owned());
        }
    }
    Ok(loci)
}

pub(crate) fn write_contig_probe_audit(
    path: &Path,
    evidence: &[ContigProbeEvidence],
) -> Result<(), String> {
    let mut writer = BufWriter::new(fs::File::create(path).map_err(|error| {
        format!(
            "Unable to write UCE fallback contig probe audit '{}': {error}",
            path.display()
        )
    })?);
    writeln!(
        writer,
        "locus	sequence_length	accepted	reason	core_anchor_status	structural_review	long_inverted_repeat	maximum_unsupported_internal_gap	target_score	target_probe_coverage	target_identity	best_other_locus	best_other_score	near_tie_other_loci"
    )
    .map_err(|error| error.to_string())?;
    for row in evidence {
        writeln!(
            writer,
            "{}	{}	{}	{}	{}	{}	{}	{}	{}	{:.6}	{:.6}	{}	{}	{}",
            row.locus,
            row.sequence_length,
            u8::from(row.accepted),
            row.reason,
            row.core_anchor_status,
            row.structural_review,
            u8::from(row.long_inverted_repeat),
            row.maximum_unsupported_internal_gap,
            row.target_score,
            row.target_probe_coverage,
            row.target_identity,
            row.best_other_locus,
            row.best_other_score,
            row.near_tie_other_loci,
        )
        .map_err(|error| error.to_string())?;
    }
    writer.flush().map_err(|error| error.to_string())
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct RecruitPass {
    pub(crate) name: &'static str,
    pub(crate) kmer_size: String,
    pub(crate) step: String,
    pub(crate) verify_kmer_size: Option<String>,
    pub(crate) max_locus_count: usize,
    pub(crate) minimum_alignment_overlap: Option<String>,
    pub(crate) minimum_alignment_identity: Option<String>,
}

impl RecruitPass {
    pub(crate) fn fast(kmer_size: &str, step: &str) -> Self {
        Self {
            name: "fast",
            kmer_size: kmer_size.to_owned(),
            step: step.to_owned(),
            verify_kmer_size: None,
            max_locus_count: 0,
            minimum_alignment_overlap: None,
            minimum_alignment_identity: None,
        }
    }

    pub(crate) fn fallback(
        kmer_size: &str,
        step: &str,
        verify_kmer_size: &str,
        minimum_alignment_overlap: &str,
        minimum_alignment_identity: &str,
    ) -> Self {
        Self {
            name: "fallback",
            kmer_size: kmer_size.to_owned(),
            step: step.to_owned(),
            verify_kmer_size: Some(verify_kmer_size.to_owned()),
            max_locus_count: 1,
            minimum_alignment_overlap: Some(minimum_alignment_overlap.to_owned()),
            minimum_alignment_identity: Some(minimum_alignment_identity.to_owned()),
        }
    }
}

pub(crate) fn read_selected_fragments(path: &Path) -> Result<BTreeMap<String, u64>, String> {
    let text = fs::read_to_string(path).map_err(|error| {
        format!(
            "Unable to read UCEFilter summary '{}': {error}",
            path.display()
        )
    })?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("UCEFilter summary '{}' is empty", path.display()))?
        .split('\t')
        .collect::<Vec<_>>();
    let locus_column = header
        .iter()
        .position(|name| *name == "locus")
        .ok_or_else(|| format!("UCEFilter summary '{}' has no locus column", path.display()))?;
    let selected_column = header
        .iter()
        .position(|name| *name == "selected_fragments")
        .ok_or_else(|| {
            format!(
                "UCEFilter summary '{}' has no selected_fragments column",
                path.display()
            )
        })?;
    let required_columns = locus_column.max(selected_column) + 1;
    let mut selected = BTreeMap::new();
    for (offset, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < required_columns {
            return Err(format!(
                "UCEFilter summary '{}' row {} has too few columns",
                path.display(),
                offset + 2
            ));
        }
        let locus = fields[locus_column];
        let count = fields[selected_column].parse::<u64>().map_err(|_| {
            format!(
                "UCEFilter summary '{}' row {} has invalid selected_fragments",
                path.display(),
                offset + 2
            )
        })?;
        if selected.insert(locus.to_owned(), count).is_some() {
            return Err(format!(
                "UCEFilter summary '{}' contains duplicate locus '{}'",
                path.display(),
                locus
            ));
        }
    }
    Ok(selected)
}

pub(crate) fn unresolved_loci(selected: &BTreeMap<String, u64>) -> BTreeSet<String> {
    selected
        .iter()
        .filter(|(_, count)| **count == 0)
        .map(|(locus, _)| locus.clone())
        .collect()
}

pub(crate) fn write_locus_allowlist(path: &Path, loci: &BTreeSet<String>) -> Result<(), String> {
    let mut writer = BufWriter::new(fs::File::create(path).map_err(|error| {
        format!(
            "Unable to create UCE fallback locus list '{}': {error}",
            path.display()
        )
    })?);
    for locus in loci {
        writeln!(writer, "{locus}").map_err(|error| error.to_string())?;
    }
    writer.flush().map_err(|error| error.to_string())
}

fn reference_locus_name(path: &Path) -> Option<String> {
    let logical_path = if path
        .extension()
        .and_then(|value| value.to_str())
        .is_some_and(|value| value.eq_ignore_ascii_case("gz"))
    {
        PathBuf::from(path.file_stem()?)
    } else {
        path.to_path_buf()
    };
    let extension = logical_path.extension()?.to_str()?.to_ascii_lowercase();
    if !matches!(extension.as_str(), "fa" | "fas" | "fasta") {
        return None;
    }
    logical_path
        .file_stem()
        .and_then(|value| value.to_str())
        .map(str::to_owned)
}

pub(crate) fn materialize_recruit_reference_subset(
    source: &Path,
    destination: &Path,
    loci: &BTreeSet<String>,
) -> Result<(), String> {
    if !source.is_dir() {
        return Err(format!(
            "UCE auto recruitment requires a reference directory, found '{}'",
            source.display()
        ));
    }
    fs::create_dir_all(destination).map_err(|error| {
        format!(
            "Unable to create UCE fallback reference subset '{}': {error}",
            destination.display()
        )
    })?;
    let mut found = BTreeSet::new();
    for entry in fs::read_dir(source).map_err(|error| {
        format!(
            "Unable to read UCE reference directory '{}': {error}",
            source.display()
        )
    })? {
        let entry = entry.map_err(|error| error.to_string())?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let Some(locus) = reference_locus_name(&path) else {
            continue;
        };
        if !loci.contains(&locus) {
            continue;
        }
        if !found.insert(locus.clone()) {
            return Err(format!(
                "UCE reference directory '{}' has multiple files for locus '{locus}'",
                source.display()
            ));
        }
        let file_name = path
            .file_name()
            .ok_or_else(|| format!("UCE reference path '{}' has no filename", path.display()))?;
        let target = destination.join(file_name);
        fs::copy(&path, &target).map_err(|error| {
            format!(
                "Unable to stage UCE fallback reference '{}' as '{}': {error}",
                path.display(),
                target.display()
            )
        })?;
    }
    let missing = loci.difference(&found).cloned().collect::<Vec<_>>();
    if !missing.is_empty() {
        return Err(format!(
            "UCE reference directory '{}' has no file-stem match for {} unresolved locus/loci: {}",
            source.display(),
            missing.len(),
            missing.join(",")
        ));
    }
    Ok(())
}

fn read_recruit_counts(path: &Path) -> Result<BTreeMap<String, u64>, String> {
    if !path.is_file() {
        return Ok(BTreeMap::new());
    }
    let text = fs::read_to_string(path).map_err(|error| {
        format!(
            "Unable to read recruit counts '{}': {error}",
            path.display()
        )
    })?;
    let mut counts = BTreeMap::new();
    for (offset, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let (locus, raw_count) = line.split_once(',').ok_or_else(|| {
            format!(
                "Recruit counts '{}' row {} is malformed",
                path.display(),
                offset + 1
            )
        })?;
        let count = raw_count.parse::<u64>().map_err(|_| {
            format!(
                "Recruit counts '{}' row {} has an invalid count",
                path.display(),
                offset + 1
            )
        })?;
        counts.insert(locus.to_owned(), count);
    }
    Ok(counts)
}

fn write_recruit_counts(path: &Path, counts: &BTreeMap<String, u64>) -> Result<(), String> {
    let mut writer = BufWriter::new(fs::File::create(path).map_err(|error| {
        format!(
            "Unable to write recruit counts '{}': {error}",
            path.display()
        )
    })?);
    for (locus, count) in counts {
        writeln!(writer, "{locus},{count}").map_err(|error| error.to_string())?;
    }
    writer.flush().map_err(|error| error.to_string())
}

pub(crate) fn merge_fallback_outputs(
    sample_dir: &Path,
    fallback_dir: &Path,
    unresolved: &BTreeSet<String>,
    fallback_selected: &BTreeMap<String, u64>,
) -> Result<BTreeSet<String>, String> {
    let recovered = unresolved
        .iter()
        .filter(|locus| fallback_selected.get(*locus).copied().unwrap_or(0) > 0)
        .cloned()
        .collect::<BTreeSet<_>>();
    let destination_filtered = sample_dir.join("filtered");
    fs::create_dir_all(&destination_filtered).map_err(|error| error.to_string())?;
    for locus in &recovered {
        let source = fallback_dir.join("filtered").join(format!("{locus}.fq"));
        if !source.is_file() {
            return Err(format!(
                "UCE fallback selected locus '{locus}' but did not write '{}'",
                source.display()
            ));
        }
        let destination = destination_filtered.join(format!("{locus}.fq"));
        if destination.exists() {
            return Err(format!(
                "UCE fallback refused to overwrite existing fast-pass output '{}'",
                destination.display()
            ));
        }
        fs::copy(&source, &destination).map_err(|error| {
            format!(
                "Unable to merge UCE fallback output '{}' into '{}': {error}",
                source.display(),
                destination.display()
            )
        })?;
    }

    let main_counts_path = sample_dir.join("ref_reads_count_dict.txt");
    let fallback_counts_path = fallback_dir.join("ref_reads_count_dict.txt");
    let mut main_counts = read_recruit_counts(&main_counts_path)?;
    let fallback_counts = read_recruit_counts(&fallback_counts_path)?;
    for locus in &recovered {
        let count = fallback_counts.get(locus).copied().ok_or_else(|| {
            format!("UCE fallback recovered locus '{locus}' without a recruit count")
        })?;
        main_counts.insert(locus.clone(), count);
    }
    write_recruit_counts(&main_counts_path, &main_counts)?;
    Ok(recovered)
}

pub(crate) fn write_recruit_audit(
    path: &Path,
    fast: &RecruitPass,
    fallback: &RecruitPass,
    fast_selected: &BTreeMap<String, u64>,
    fallback_selected: &BTreeMap<String, u64>,
    unresolved: &BTreeSet<String>,
    recovered: &BTreeSet<String>,
) -> Result<(), String> {
    let fallback_verify = fallback.verify_kmer_size.as_deref().unwrap_or("default");
    let fallback_overlap = fallback
        .minimum_alignment_overlap
        .as_deref()
        .unwrap_or("off");
    let fallback_identity = fallback
        .minimum_alignment_identity
        .as_deref()
        .unwrap_or("off");
    let mut writer = BufWriter::new(fs::File::create(path).map_err(|error| {
        format!(
            "Unable to write UCE recruit audit '{}': {error}",
            path.display()
        )
    })?);
    writeln!(
        writer,
        "locus	fast_k	fast_step	fast_selected_fragments	fallback_attempted	fallback_k	fallback_step	fallback_verify_k	fallback_min_alignment_overlap	fallback_min_alignment_identity	panel_unique_locus_only	fallback_selected_fragments	final_pass"
    )
    .map_err(|error| error.to_string())?;
    for (locus, fast_count) in fast_selected {
        let attempted = unresolved.contains(locus);
        let fallback_count = fallback_selected.get(locus).copied().unwrap_or(0);
        let final_pass = if *fast_count > 0 {
            fast.name
        } else if recovered.contains(locus) {
            fallback.name
        } else {
            "unresolved"
        };
        writeln!(
            writer,
            "{locus}	{}	{}	{fast_count}	{}	{}	{}	{fallback_verify}	{fallback_overlap}	{fallback_identity}	{}	{fallback_count}	{final_pass}",
            fast.kmer_size,
            fast.step,
            u8::from(attempted),
            fallback.kmer_size,
            fallback.step,
            u8::from(fallback.max_locus_count == 1),
        )
        .map_err(|error| error.to_string())?;
    }
    writer.flush().map_err(|error| error.to_string())
}

pub(crate) fn preserve_summary(source: &Path, destination: &Path) -> Result<PathBuf, String> {
    fs::copy(source, destination).map_err(|error| {
        format!(
            "Unable to preserve UCEFilter summary '{}' as '{}': {error}",
            source.display(),
            destination.display()
        )
    })?;
    Ok(destination.to_path_buf())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_TEST_DIR: AtomicU64 = AtomicU64::new(0);

    fn test_root(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "gm2_uce_recruit_{name}_{}_{}",
            std::process::id(),
            NEXT_TEST_DIR.fetch_add(1, Ordering::Relaxed)
        ))
    }

    #[test]
    fn summary_columns_are_resolved_by_name() {
        let root = test_root("summary");
        fs::create_dir_all(&root).unwrap();
        let summary = root.join("summary.tsv");
        fs::write(
            &summary,
            "selected_fragments\tlocus\tcoarse_reads\n0\tuce-1\t0\n3\tuce-2\t8\n",
        )
        .unwrap();
        let selected = read_selected_fragments(&summary).unwrap();
        assert_eq!(selected["uce-1"], 0);
        assert_eq!(selected["uce-2"], 3);
        assert_eq!(unresolved_loci(&selected), BTreeSet::from(["uce-1".into()]));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn fallback_merge_is_locus_scoped_and_updates_counts() {
        let root = test_root("merge");
        let sample = root.join("sample");
        let fallback = root.join("fallback");
        fs::create_dir_all(sample.join("filtered")).unwrap();
        fs::create_dir_all(fallback.join("filtered")).unwrap();
        fs::write(
            sample.join("ref_reads_count_dict.txt"),
            "uce-1,2\nuce-2,8\n",
        )
        .unwrap();
        fs::write(
            fallback.join("ref_reads_count_dict.txt"),
            "uce-1,14\nuce-3,7\n",
        )
        .unwrap();
        fs::write(fallback.join("filtered/uce-1.fq"), b"@r\nAC\n+\nII\n").unwrap();
        fs::write(fallback.join("filtered/uce-3.fq"), b"@r\nGT\n+\nII\n").unwrap();
        let unresolved = BTreeSet::from(["uce-1".into(), "uce-3".into()]);
        let selected = BTreeMap::from([("uce-1".into(), 4), ("uce-3".into(), 0)]);
        let recovered = merge_fallback_outputs(&sample, &fallback, &unresolved, &selected).unwrap();
        assert_eq!(recovered, BTreeSet::from(["uce-1".into()]));
        assert!(sample.join("filtered/uce-1.fq").is_file());
        assert!(!sample.join("filtered/uce-3.fq").exists());
        assert_eq!(
            fs::read_to_string(sample.join("ref_reads_count_dict.txt")).unwrap(),
            "uce-1,14\nuce-2,8\n"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn fallback_recruit_reference_contains_only_unresolved_loci() {
        let root = test_root("reference_subset");
        let source = root.join("source");
        let subset = root.join("subset");
        fs::create_dir_all(&source).unwrap();
        fs::write(source.join("uce-1.fasta"), b">uce-1\nACGT\n").unwrap();
        fs::write(source.join("uce-2.fa"), b">uce-2\nTGCA\n").unwrap();
        fs::write(source.join("uce-3.fasta.gz"), b"compressed fixture").unwrap();
        fs::write(source.join("notes.txt"), b"not a reference locus\n").unwrap();

        materialize_recruit_reference_subset(
            &source,
            &subset,
            &BTreeSet::from(["uce-2".into(), "uce-3".into()]),
        )
        .unwrap();
        assert!(!subset.join("uce-1.fasta").exists());
        assert_eq!(
            fs::read_to_string(subset.join("uce-2.fa")).unwrap(),
            ">uce-2\nTGCA\n"
        );
        assert!(!subset.join("notes.txt").exists());
        assert_eq!(
            fs::read(subset.join("uce-3.fasta.gz")).unwrap(),
            b"compressed fixture"
        );
        assert!(materialize_recruit_reference_subset(
            &source,
            &root.join("missing"),
            &BTreeSet::from(["uce-4".into()])
        )
        .is_err());
        fs::remove_dir_all(root).unwrap();
    }

    fn deterministic_sequence(mut state: u32, length: usize) -> String {
        (0..length)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                b"ACGT"[((state >> 29) & 3) as usize] as char
            })
            .collect()
    }

    #[test]
    fn final_contig_probe_gate_requires_target_coverage_and_rejects_panel_ties() {
        let root = test_root("contig_probe_gate");
        let references = root.join("references");
        fs::create_dir_all(&references).unwrap();
        let target = deterministic_sequence(1, 120);
        let partial = deterministic_sequence(2, 120);
        let ambiguous = deterministic_sequence(3, 120);
        let short = deterministic_sequence(4, 120);
        for (locus, sequence) in [
            ("target", target.as_str()),
            ("partial", partial.as_str()),
            ("ambiguous", ambiguous.as_str()),
            ("ambiguous-copy", ambiguous.as_str()),
            ("short", short.as_str()),
        ] {
            fs::write(
                references.join(format!("{locus}.fasta")),
                format!(">{locus}_probe\n{sequence}\n"),
            )
            .unwrap();
        }
        let contigs = BTreeMap::from([
            (
                "target".into(),
                format!(
                    "{}{target}{}",
                    deterministic_sequence(10, 40),
                    deterministic_sequence(11, 40)
                ),
            ),
            (
                "partial".into(),
                format!(
                    "{}{}{}",
                    deterministic_sequence(12, 80),
                    &partial[..57],
                    deterministic_sequence(13, 80)
                ),
            ),
            (
                "ambiguous".into(),
                format!(
                    "{}{ambiguous}{}",
                    deterministic_sequence(14, 40),
                    deterministic_sequence(15, 40)
                ),
            ),
            (
                "short".into(),
                format!(
                    "{}{short}{}",
                    deterministic_sequence(16, 20),
                    deterministic_sequence(17, 20)
                ),
            ),
        ]);
        let evidence = evaluate_contig_probe_support_parallel(&references, &contigs, 1).unwrap();
        let parallel = evaluate_contig_probe_support_parallel(&references, &contigs, 3).unwrap();
        let overprovisioned =
            evaluate_contig_probe_support_parallel(&references, &contigs, 64).unwrap();
        assert_eq!(parallel, evidence);
        assert_eq!(overprovisioned, evidence);
        let by_locus = evidence
            .iter()
            .map(|row| (row.locus.as_str(), row))
            .collect::<BTreeMap<_, _>>();
        assert!(by_locus["target"].accepted);
        assert_eq!(by_locus["target"].reason, "pass");
        assert!(!by_locus["partial"].accepted);
        assert_eq!(
            by_locus["partial"].reason,
            "target_probe_coverage_below_0.80"
        );
        assert!(by_locus["partial"].best_other_locus.is_empty());
        assert!(!by_locus["ambiguous"].accepted);
        assert_eq!(by_locus["ambiguous"].reason, "near_tie_other_locus");
        assert!(!by_locus["short"].accepted);
        assert_eq!(by_locus["short"].reason, "contig_length_below_200");
        assert_eq!(by_locus["short"].target_score, 0);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn fallback_loci_are_read_from_audit_columns_by_name() {
        let root = test_root("audit_loci");
        fs::create_dir_all(&root).unwrap();
        let audit = root.join("audit.tsv");
        fs::write(
            &audit,
            "final_pass\tignored\tlocus\nfast\tx\tuce-1\nfallback\tx\tuce-2\nunresolved\tx\tuce-3\n",
        )
        .unwrap();
        assert_eq!(
            fallback_recruited_loci(&audit).unwrap(),
            BTreeSet::from(["uce-2".into()])
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn provisional_core_structure_checks_reject_inversion_but_review_gap() {
        let mut clean = ContigProbeEvidence::unavailable("uce-clean", "fixture");
        clean.accepted = true;
        clean.reason = "pass";
        clean.core_anchor_status = "probe_pass";
        clean.apply_structure_checks(false, 0);
        assert!(clean.accepted);
        assert_eq!(clean.core_anchor_status, "anchored");
        assert_eq!(clean.structural_review, "none");

        let mut gap = clean.clone();
        gap.locus = "uce-gap".into();
        gap.apply_structure_checks(false, crate::rescue_qc::MINIMUM_NEW_UNSUPPORTED_GAP);
        assert!(gap.accepted);
        assert_eq!(gap.core_anchor_status, "anchored_with_review");
        assert_eq!(gap.structural_review, "internal_gap_ge40");

        let mut inverted = clean.clone();
        inverted.locus = "uce-inverted".into();
        inverted.apply_structure_checks(true, 0);
        assert!(!inverted.accepted);
        assert_eq!(inverted.reason, "long_inverted_repeat");
        assert_eq!(inverted.core_anchor_status, "structure_rejected");

        let root = test_root("core_audit");
        fs::create_dir_all(&root).unwrap();
        let audit = root.join("audit.tsv");
        write_contig_probe_audit(&audit, &[clean, gap, inverted]).unwrap();
        let text = fs::read_to_string(&audit).unwrap();
        assert!(text
            .lines()
            .next()
            .unwrap()
            .contains("maximum_unsupported_internal_gap"));
        assert!(text.contains("uce-gap	0	1	pass	anchored_with_review	internal_gap_ge40"));
        fs::remove_dir_all(root).unwrap();
    }
}
