# 3. Profiling

[中文版本](profiling_ZH.md) · [project overview](../README_EN.md) · [command-line guide](../manual/EN_US/command_line.md) · [output guide](../manual/EN_US/output.md)

By default, profiling is an assembly-free workflow for recovering reference-level evidence for any amplicon marker from WGS or metagenomic reads. It reports compatibility with reference sequences, not organismal abundance. Optional `--profile-copy-mode host|symbiont` adds ITS2 multi-copy assembly and read-support checks.

## Workflow

```text
TipSeek k-mer recruitment → Themisto pseudoalignment → reference-level support
```

Without copy mode, it runs neither `refilter`, `assemble`, `combine`, nor `tree`. Recruitment and Themisto both use k=21.

## Inputs

Always required:

- One `.fa` or `.fasta` marker library passed to `-r`.
- A sample table with WGS or metagenomic reads.
- `themisto` on `PATH`.

Optional:

- `--profile-group-map`, a two-column TSV: `reference_id<TAB>group`. It is an optional annotation column in the reference-level result.
- `--profile-decoy`, a FASTA of plausible non-target sequences.

The reference ID is the first whitespace-delimited FASTA-header field. A supplied map must cover every reference exactly once; duplicate rows are allowed only when they assign the same group.

## Run and cache

```bash
cli/tipseek profiling \
  -f samples.tsv -r marker_reference.fasta \
  -o output -p 8
```

Add `--profile-group-map marker_groups.tsv` only when the `group` annotation column is useful. Themisto indices use a content-addressed cache. Use the common `--reference-cache-dir` option to share it. To rebuild, use a fresh cache directory. When calling `marker_profile` directly, use a separate cache or `--force-rebuild` after changing the reference, group map, or decoy.

The only profiling-specific options are `--profile-copy-mode`, `--profile-group-map`, and `--profile-decoy`. Recruitment/assembly/index k=21, pseudoalignment threshold=0.8, relevant k-mer fraction=0.5, and index memory=2 GB retain the validated defaults. Index memory is per concurrent build, not a process memory cap. Use `-p` to control concurrent samples. Former profiling tuning, tool-path and cache flags are rejected; remove them and use `PATH` and `--reference-cache-dir` where needed. `GM2_THEMISTO` is no longer read. For experimental non-default pseudoalignment, invoke `marker_profile` directly; its outputs are not equivalent to the validated integrated defaults.

Parameter simplification is checked against the prior defaults on the same three WGS samples in both modes; see the [comparison protocol and results](development/profiling-parameter-validation.md). Shared cache keys are independent of the output directory and change when reference inputs change.

## Outputs and interpretation

Every sample writes `marker_profile/`.

### Primary result: `marker_reference_support.tsv`

One row is written for every hit reference sequence:

- `hit_queries`: number of queries compatible with that reference.
- `fractional_queries`: shared-query support; a query with N candidates contributes `1/N` to each, so it is never counted N times.
- `singleton_queries`: queries compatible with this reference alone.
- `ambiguity_status`: `has_singleton_support` or `shared_only`.

This is evidence of compatibility with a reference sequence, not a claim that the reference is uniquely present, and not an organismal abundance.

`marker_qc.tsv` records pseudoalignment and run parameters. `marker_reference_metadata.tsv` in the index cache records the Themisto color-to-reference mapping and optional group annotation.

## QC and calibration

Check pseudoaligned queries, each reference’s fractional and singleton support, and the decoy-reference evidence before interpreting a result. Use negative controls, mixtures, and depth-matched downsampling to evaluate reference libraries and the fixed pseudoalignment defaults. Keep the reference library, optional annotation map, k-mer size, thresholds, and decoy strategy fixed for cross-sample comparisons.

See [Filter](filter_EN.md), the [output guide](../manual/EN_US/output.md), and the [command-line guide](../manual/EN_US/command_line.md).

## ITS2 multi-copy recovery

`--profile-copy-mode host|symbiont` adds ITS2 graph assembly to the same profiling workflow. Both modes use the same evidence checks; choose the target reference with `-r`, and provide competing/non-target sequences with `--profile-decoy`. The mode records the intended biological interpretation and does not infer a sequence's origin by itself. Omitting the option retains assembly-free profiling.

Copy recovery recruits against the combined target/decoy panel, retains linked mates and shared reference matches, and disables abundance thinning. Explicit input read limits still apply. References provide recruitment and graph seeds, never substitute for sample sequence. Distinct assembled sequences are retained separately. Read support checks full overlaps with conservative base-quality tolerance, pair orientation and links between distinguishing sites. Sequencing-error tolerance and biological ITS boundaries still require further validation; `READ_SUPPORTED` means the candidate sequence has supporting evidence, not that every genomic copy or the complete ITS region was recovered.

Per sample, `marker_profile/its2_candidates.fasta` contains accepted and unresolved candidates. `marker_profile/its2_candidates.tsv` joins stable candidate IDs, coverage, phase/acceptance status, rejection reason, copy mode and independent source annotations. `PARTIAL` indicates insufficient read coverage or a partial reference span; `UNPHASED` indicates missing links between distinguishing sites or unresolved competing paths. EM estimates are `NA` for unresolved candidates and must not be interpreted as cell abundance.

Source annotations are relative to the supplied panel: `target_only`, `decoy_only`, `ambiguous_source`, `unclassified`, or `target_match_no_decoy`. In particular, reads support does not resolve an ambiguous source, and no-decoy matches do not establish exclusion of other origins. Graphs and rejected candidates are retained for inspection. An initial benchmark on three real WGS samples is complete, as described below; recovery of all low-frequency copies is not established.

At >=80% alignment identity to the target reference, ITS2 recovery trims to the aligned sample interval, including partial spans, before recalculating reads support. Terminal substitutions are evaluated using observed sample bases. Sample insertions are retained; partial spans are not padded. The report includes reference-relative boundaries and alignment identity; `its2_reference_spans.fasta` contains oriented aligned sample intervals, including partial/unresolved candidates. These annotations describe the supplied marker reference and do not prove complete genomic rDNA copies.

Copy recovery currently requires distinct paired R1 and R2 files; supplying one file twice is rejected because it cannot establish paired evidence.

Graph extension uses observed read contexts, falling back to a shorter observed context when longer context coverage is absent. A reference endpoint that also occurs internally in another reference does not stop all paths. Extension on each side is bounded by the longest reference length plus `2 × k`; unusually large insertions may remain partial. ITS2 visits all remaining graph seeds, avoiding loss of smaller disconnected components after half the seeds have been visited. Only fully read-covered, full-reference-span candidates enter final phasing and competition; partial-reference models cannot divert paired evidence from full models. Models strictly contained in longer covered full candidates remain `UNPHASED`; all rejected candidates remain available for inspection. Terminal reads use the existing k-mer length as the minimum exact overlap, with the entire overlap and pair orientation still checked.

Diagnostic and variant-site bridging support are deduplicated by mapped read endpoints and orientation, so repeated positions cannot manufacture independent phasing evidence. This is an approximation, not UMI molecule counting; raw compatible fragment counts are retained. Combinations without joint evidence remain `UNPHASED`; amplicon controls never select assembly models.

Reference and candidate coordinates are zero-based, half-open intervals; `reference_record` is the 1-based record ordinal in the reference index. Outputs distinguish sequence variants, not the number of genomic copies with identical sequences.

### Initial WGS validation

Three WGS samples from [Denis et al.](https://doi.org/10.1016/j.cub.2026.04.057), NC7_3, NC2_12 and NC1_14, were compared with independent replicate SymPortal amplicons. Recruitment uses independent Symb-SHIN/GenBank panels; amplicons evaluate outputs only and never select assemblies. Denominators below are distinct sequences in the replicate union.

| Sample | Complete amplicon sequences contained in candidates: before → after | Read-supported, source-resolved, full-reference-span matches: before → after | Symbiont elapsed: before → after |
| --- | --- | --- | --- |
| NC7_3 | 7/19 → 8/19 | 3/19 → 1/19 | 20.55 → 28.56 s |
| NC2_12 | 9/21 → 10/21 | 2/21 → 0/21 | 29.30 → 36.19 s |
| NC1_14 | 10/16 → 10/16 | 1/16 → 1/16 | 35.47 → 37.64 s |

More complete candidates do not establish more reliably phased haplotypes. Broader competition exposes previously missing compatible combinations, and mapped-position deduplication prevents duplicate reads from manufacturing independent diagnostic evidence. Strict acceptance consequently decreases; recovery of all low-frequency copies is not established. NC7_3 C3ba and C3bd have local compatible fragments at zero-based sites 96 and 218, but none covering both sites. Several complete NC2_12 models remain indistinguishable. Unresolved sequences are retained without using amplicon truth to force phasing.

Evaluation-only recruitment against each sample's amplicons found 154, 221 and 207 fragments, all already present in the independently recruited pools. This assay found no additional panel-related recruitment loss, but cannot exclude shared recruiter blind spots involving errors or N bases. These evaluation pools never enter assembly.

The NC1_14 2 nt candidate is a GC/CG repeat expansion with equivalent insertion boundaries. Four event-spanning reads represent two distinct starts and both orientations; nearby instrument coordinates within each pair suggest optical duplicates rather than four independent templates. High-quality local insertion evidence does not establish a complete biological haplotype. The candidate remains `UNPHASED`, without an amplicon match, and needs independent-library or long-read confirmation. Amplicon absence alone does not prove a false positive.

Host mode yields 27, 17, 18 read-supported, source-resolved, full-reference-span candidates, taking 38.63, 17.41, 18.16 s. There is no matched host ITS truth, so these are not host recall estimates. All workflows start from the same raw WGS reads with `-p 1`; download time is excluded. Amplicon agreement is not recall against a complete biological copy inventory or a cell-abundance estimate.

Broader search increases total runtime; this revision does not claim an overall speedup. Three alternating runs on the same NC7_3 symbiont recruited pool measured 11.69 s for accelerated overlap search versus 11.45 s for the original, and 11.28 versus 11.35 s for an exact-match DP bypass. Both produced byte-identical outputs, but neither showed sufficient reliable benefit to retain. The earlier coverage-first phasing optimization remains: 32.76 → 9.59 s median on NC7_3 host recruited reads, with byte-identical outputs across three alternating pairs.

### Sequencing-error sensitivity and phasing audit

An evaluation-only audit compared exact overlaps of at least 21 nt with one Q≤12 substitution in overlaps of at least 75 nt. Mismatching bases never contribute coverage; variable positions among equal-length control amplicons were excluded from error rescue. Controls are used for evaluation only, never production assembly or phasing.

Control sequences with at least two compatible fragments per base changed from 5→6, 12→12 and 11→11; only NC7_3 C3st gained complete compatible coverage. Allowing up to two Q≤12 substitutions in overlaps of at least 100 nt added no further complete sequences. This shared-read coverage assay does not enforce unique placement, pair orientation or independent phasing, so these counts are not haplotype recall. Unrestricted one-mismatch matching admits many Q40 differences: potential new C3ba bridges depend on a Q40 mismatch at site 119. Such evidence cannot simply be treated as sequencing error. Production exact matching remains unchanged, with no new parameter.

Adjacent-site bridges now count distinct mapped fragment positions rather than duplicate multiplicity. Two local diagnostic positions do not substitute for two independent bridge positions. A regression test covers adequate local diagnosis with a single duplicated bridge. Models without sufficient joint evidence remain unresolved.

A shared seed index nominates candidate offsets before full-overlap verification, preserving repeat ambiguity and terminal matching. Four hundred cases containing short sequences, repeats, reverse complements, N bases and clipped ends match the exhaustive implementation exactly. Three alternating pairs on the NC7_3 symbiont recruited pool reduced median elapsed time from 11.39 to 10.36 s (about 9%), with byte-identical FASTA, support tables, summaries and GFA. Peak memory stayed near 147 MiB. This performance experiment isolates indexing before the bridge-deduplication change.

Full symbiont workflows were rerun from raw WGS with `-p 1`, excluding download:

| Sample | Complete amplicon sequences in candidates: before→after | Strict accepted matches: before→after | Full workflow elapsed: before→after |
| --- | --- | --- | --- |
| NC7_3 | 8/19 → 8/19 | 1 → 1 | 28.56 → 27.57 s |
| NC2_12 | 10/21 → 10/21 | 0 → 0 | 36.19 → 34.94 s |
| NC1_14 | 10/16 → 10/16 | 1 → 1 | 37.64 → 35.48 s |

Candidate sequences are unchanged across all six host/symbiont workflows. Bridge deduplication changes one source-unclassified NC7_3 host-mode candidate from `READ_SUPPORTED` to `UNPHASED`; strict symbiont match counts are unchanged.

Reliable low-frequency haplotype recall did not increase in this revision; production error tolerance still needs independent mixture calibration. All 144 tests and static checks pass. Full-workflow timings are single measurements; repeated evidence for speedup comes from the separate alternating benchmark.

### Quality-aware compatible-evidence rescue

Original Phred+33 qualities are retained with paired FASTQ sequences; quality-free FASTA remains exact-only. Rescue permits one Q≤12 substitution in an overlap of at least 75 nt, without indel or N tolerance and without rewriting reads or candidate sequences. Mismatching bases never count as coverage, diagnosis or bridges. Matching bases used to distinguish models or bridge variant sites require Q≥20. Exact and tolerated models remain jointly compatible; tolerance does not erase competing alleles. These internal conservative rules add no user parameter and are not a trained DADA2/WGS error model.

Exact and tolerated offsets share the existing seed-index pass. Identical reads with different qualities are evaluated separately; diagnostic and bridge evidence remains deduplicated by mapped positions. `quality_rescued_fragment_support` reports raw compatible fragment multiplicity involving tolerated placements, not independent molecules or unique haplotype evidence.

NC7_3 C3st improves from 257/263 to 263/263 covered bases, closing the 6 bp gap and changing `PARTIAL` to `UNPHASED`. It has six compatible rescue fragment counts but only one distinct diagnostic position. An independent audit reproduces production counts: among five full competing models, zero-based edges 96–146 and 146–218 each have only one independent high-quality bridge, below the requirement of two. The two key low-quality reads also exactly match other models, so they cannot force acceptance of C3st.

| Sample | Complete amplicon sequences in candidates | Strict accepted matches: before→after | Symbiont workflow elapsed: before→after |
| --- | --- | --- | --- |
| NC7_3 | 8/19 | 1 → 1 | 27.57 → 27.55 s |
| NC2_12 | 10/21 | 0 → 0 | 34.94 → 34.79 s |
| NC1_14 | 10/16 | 1 → 1 | 35.48 → 35.43 s |

No sample gains an accepted source-resolved full-reference candidate absent from its amplicon control. This is not a measured biological false-positive rate: amplicon absence alone does not establish false positivity. A targeted synthetic negative with two distinct Q12 error-allele positions yields one false acceptance under exact-only evidence and zero with quality-aware evidence. A true minor allele with Q40 observations remains accepted. Broader mixture calibration is still needed.

All 150 tests and static checks pass, including 200 comparisons of tolerated indexed placements against exhaustive search. Three alternating pairs on the same NC7_3 recruited pool give median times of 10.37→10.34 s, within timing variation; no additional speedup is claimed. Repeat outputs are stable within each version. Full-workflow times above are single measurements excluding download.

Method references: [Aspat GraftM ITS2](https://github.com/hde08/Aspat_WesternPacific_Popgen/blob/master/Scripts/Symbionts/GraftM_ITS2.sh) and [Symb-SHIN](https://github.com/hisatakeishida/Symb-SHIN/blob/main/C_Read-based_alignment.md) inform WGS marker recruitment/classification; [DADA2](https://pmc.ncbi.nlm.nih.gov/articles/PMC4927377/) motivates quality-dependent error assessment; [weebill architecture](https://github.com/wwood/weebill/blob/208eec3cba3aae60198b200d84dd2224b1d84ae9/ARCHITECTURE.md) informs staged screening and shared evidence. Genome-level sparse sketches, ANI and abundance models are not directly used for short multicopy ITS2 phasing. None is added as a dependency.

### Quality likelihood and phase connectivity (v15)

Inspired by [Lorikeet's candidate likelihood calculation](https://github.com/rhysnewell/Lorikeet/blob/master/src/pair_hmm/pair_hmm_likelihood_calculation_engine.rs) and [Floria's quality-aware haplotyping](https://github.com/bluenote-1577/floria), this stage reuses existing candidate placements and EM without new CLI parameters or dependencies. It does not implement their complete PairHMM, clustering or flow algorithms.

For equal-length candidates with identical unique placement geometry per mate, relative Phred-based base likelihoods now weight EM responsibilities. Common terms cancel. Overlapping mates count each position once using the higher quality; discordant bases contribute no likelihood. Different lengths, clipping, ambiguous placements or missing qualities retain equal weights, avoiding unmodelled indel/alignment-length bias. EM does not determine acceptance and is not cell abundance.

Phasing now requires a connected graph of distinguishing sites, rather than a supported edge at every adjacent pair. Nonadjacent links can contribute, but each edge still requires two distinct mapped fragment positions; acceptance still requires two independent diagnostic positions and complete coverage. Each edge stores at most two positions to bound memory. The coverage-only first pass no longer stores an unused EM compatibility matrix.

All six host/symbiont workflows completed from raw WGS reads. Candidate sequences and acceptance statuses match v14. Five models (four host, one symbiont) changed rejection reason from missing phase links to unresolved competing paths; they remain unresolved.

| Sample | Complete amplicon sequence candidates: before→after | Strictly accepted: before→after | Symbiont workflow: before→after |
| --- | --- | --- | --- |
| NC7_3 | 8/19 → 8/19 | 1 → 1 | 27.55 → 27.44 s |
| NC2_12 | 10/21 → 10/21 | 0 → 0 | 34.79 → 34.87 s |
| NC1_14 | 10/16 → 10/16 | 1 → 1 | 35.43 → 35.92 s |

C3st retains 263/263 bp coverage, 73 compatible fragments and one independent diagnostic position. An independent audit finds three distinct bridges across nonadjacent sites 96–218, but only one each across 96–146 and 146–218. Site 146 remains disconnected at the two-position threshold, so the model remains `UNPHASED`. A uniform-prior likelihood audit gives unequal weights to only one compatible fragment; the C3st fractional score changes from 22.00 under equal compatibility to 21.68. This is neither independent support nor production EM abundance. Quality modelling cannot create missing phase information.

No accepted, source-resolved, complete symbiont candidate lacks an amplicon match. This is not a measured biological false-positive rate. The targeted low-quality false-SNP control has zero false acceptances in both v14 and v15, while the high-quality true minor allele survives. Independent host amplicon truth is unavailable. All 152 tests passed, followed by an additional passing reverse-strand quality-position assertion and clean static checks. Workflow timings above are single observations.

Three alternating timing pairs on the same NC7_3 recruited pool gave medians of 10.35→10.35 s, with stable repeated outputs within each version; no speedup was demonstrated. Cross-version output differs in quality-weighted EM estimates, while candidate sequences and acceptance statuses remain unchanged.

### Raw-read audit of C3st site 146 (v16)

Site 146 is zero-based (147 in one-based coordinates). An independent scan of all 21,143,328 NC7_3 WGS pairs used every 15-mer and reverse complement from the existing WGS candidate, retaining both mates when either matched. Amplicon sequences did not drive production recruitment or assembly. Of 2,842 seed-positive pairs, 2,689 were absent from the recruitment pool. Gapped local alignment requiring at least 75 nt, 95% identity and coverage of the site retained 72 reads, all already recruited: 70 carried A and two carried the C3st G. No missing independent bridge was found under this audit criterion; a short seed hit alone is not evidence of missed recruitment.

Both G observations have Q40 and map to `[76,227)`. Their other mates lie outside the marker and have identical sequences. Their sequencing cluster coordinates on the same lane/tile are `(37723,8872)` and `(37731,8886)`, consistent with nearby duplicates, although no UMI proves template identity. Both are compatible only with C3st among the five final complete models: competing models did not remove them. A Q12 substitution at site 168 in one read is already rescued by the existing error-tolerant placement rule. They still represent one independent mapped position.

Two result columns are appended without new CLI parameters: `phase_components` counts components connected by edges with two distinct mapped positions; `min_variant_site_support` reports the minimum high-quality distinct-position support across distinguishing sites. Candidates not assessed in complete-model phasing receive `NA`; no distinguishing sites yields `0` and `NA`, respectively. Connectivity does not replace independent diagnostic support. The coverage-only pass now avoids unused fragment-position signatures. Acceptance thresholds and likelihood allocation are unchanged.

Methodological references remain [Lorikeet's quality likelihoods and preservation of original reads](https://github.com/rhysnewell/Lorikeet/blob/master/src/pair_hmm/pair_hmm_likelihood_calculation_engine.rs) and [Floria's fragment/site quality evidence](https://github.com/bluenote-1577/floria). Neither tool is a new dependency; missing phase evidence is not replaced by forced base correction.

All 153 tests and static checks passed, including regressions preventing duplicate positions or low-quality bridges from supplying independent evidence. The independent raw-read scanner also matched naive seed searching on 500 synthetic pairs. All six raw-WGS workflows completed. Candidate sequences and every pre-existing TSV field match the baseline: acceptance, rejection reasons and EM estimates are unchanged.

| Sample | Complete amplicon sequence candidates | Strictly accepted: before→after | Symbiont workflow: before→after |
| --- | --- | --- | --- |
| NC7_3 | 8/19 | 1 → 1 | 27.44 → 27.58 s |
| NC2_12 | 10/21 | 0 → 0 | 34.87 → 34.88 s |
| NC1_14 | 10/16 | 1 → 1 | 35.92 → 35.96 s |

C3st reports `phase_components=2` and `min_variant_site_support=1`, remaining `UNPHASED` with 263/263 bp coverage. No accepted, source-resolved, full-reference symbiont candidate lacks an amplicon match; this does not establish a zero biological false-positive rate. The targeted false-SNP control has zero false acceptances in v15 and v16; the high-quality true minor allele remains accepted. Independent host amplicon truth is unavailable. Workflow timings are single observations and exclude this raw-read audit.

Three alternating pairs on the same NC7_3 recruited pool gave medians of 10.36→10.39 s (about +0.3%), with stable repeated outputs within each version. This stage adds diagnostics and demonstrates no speedup. Cross-version file differences are the added columns; existing fields and candidate sequences remain unchanged.

### ITS2 alignment performance

Stage timing on NC7_3 symbiont reads identified reference-boundary normalization (~8.6 s) and full-span checks (~1.2 s) as the bottlenecks; read evidence took ~0.03 s. Affine-gap local alignment now keeps two score rows and packs the three traceback states into one byte per cell. Traversal, scoring and tie-breaking remain unchanged, as do candidate, acceptance and phase rules. Existing `--profile` output includes these three substages; no new parameters are introduced.

All 154 tests and static checks passed, including 600 field-by-field comparisons against the previous full-matrix algorithm covering indels, repeats, ties and different reference windows. The old matrix remains only as a test oracle.

The final version was rerun from raw reads for all six host/symbiont workflows across three samples. Recruited FASTQ, candidate FASTA/TSV and reference-span FASTA are byte-identical; strict acceptance and amplicon matches are unchanged. Full-workflow timings below are single observations.

| Sample | Host: before→after | Symbiont: before→after |
| --- | --- | --- |
| NC7_3 | 30.23 → 29.17 s | 27.58 → 27.34 s |
| NC2_12 | 16.75 → 16.56 s | 34.88 → 34.35 s |
| NC1_14 | 17.53 → 17.50 s | 35.96 → 35.55 s |

Three alternating single-threaded timing pairs on the same NC7_3 symbiont recruited pool gave medians of **10.39→10.12 s**, a **2.6% elapsed-time reduction**. FASTA, support tables, GFA and summary files are byte-identical across versions and repeats. Recovery and false-acceptance thresholds are unchanged.
