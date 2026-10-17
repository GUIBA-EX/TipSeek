# Exonerate integration: stages 1–3

The optional `--gene-exonerate-shadow` flag reaches the exon annotation component through the normal TipSeek CLI. The annotation component also accepts `--exonerate-shadow` directly. Both flags take no value and invoke `exonerate-rs` through the inherited PATH. Miniprot remains authoritative. Default runs create no shadow files and invoke no Exonerate process.

This initial integration uses a separate executable to validate the interface without adding a library dependency or changing model selection. Both per-candidate annotation and attempted padded joins are covered. Execution failures and invalid output are recorded independently from no-hits. Each alignment has a 60-second timeout. Arguments are passed without a shell; input snapshots and version output remain available after temporary workflow cleanup.

The parser validates query/target identity, oriented endpoints, sequence bounds, supported operations, and summed advances. It reconstructs CDS intervals from vulgar paths, excluding splice signals/introns and retaining split-codon fragments. All reported intervals are ascending zero-based half-open coordinates. The comparison deliberately requires exact coordinates; it does not equate different scoring systems or synthesize miniprot PAF evidence.

## Validation

- 297 workspace tests pass, including reverse-strand split codons, malformed paths, timeouts, failed executables, CLI forwarding/exon-only validation, and byte-identical authoritative outputs with a differing or missing shadow executable.
- Strict workspace Clippy and formatting checks pass.
- The existing 56-gene Patiria pectinifera dataset was run with and without shadow mode using miniprot 0.18: 56 intact candidate comparisons; 112 fragmented candidate comparisons plus 30 attempted padded joins.
- All 392 compared intact output files and 495 fragmented output files were byte-identical between modes. This includes model/segment tables, CDS/proteins, GFF3, unresolved outputs, reference translation audit and padding decisions.
- There were no Exonerate execution/parse failures. Of the intact report rows, 86 had strict coordinate differences. Fragmented rows included 139 coordinate differences, 30 reported-hit matches, 5 Exonerate-only hits and 5 both-no-hit records. Rows are not biological locus counts or unique comparison calls.

Evidence, exact commands, binary/input SHA-256 hashes, output locations and timings are in [exonerate-shadow-validation.json](exonerate-shadow-validation.json). The initial report retains its historical explicit-path commands. The current valueless flags and inherited-PATH lookup were revalidated on the same dataset in [exonerate-shadow-path-validation.json](exonerate-shadow-path-validation.json), again with all 887 authoritative output files unchanged and no Exonerate execution/parse failures. This rerun prepended the existing release binary directory to the test process PATH because the desktop shell did not already expose `exonerate-rs`. Raw run data remains at the local output location recorded there. This verifies non-interference and instrumentation, not annotation correctness or speed superiority. It does not rerun read recruitment or assembly.

## Stage 2: shared evidence and structural QC

`AlignmentEvidence` records the engine, protein length and interval, score, frameshifts, reported stops and large query gaps. The miniprot adapter retains PAF/CIGAR semantics, and the common model builder now requires alignment evidence rather than a PAF object. Missing miniprot evidence still produces the legacy `missing_paf` flag and cannot become resolve-eligible.

The Exonerate adapter validates protein-to-genome operation advances, translates M operations and joined S split codons, and derives CDS phases in transcript order on either strand. F operations remain structural failures. G query gaps of at least 10 amino acids use the existing large-gap flag. Target insertions remain in the CDS. Terminal stops are retained only when present within the reported coding path; an adjacent unaligned stop is not silently appended. Internal stops, phase errors, ambiguous codons, coverage and splice sites use the same QC implementation as miniprot.

Identity is the exact amino-acid match fraction among aligned M/S pairs, excluding gaps and treating translated X as non-identical. This definition includes split codons and is explicitly distinct from either engine's native percent identity. Substitution-matrix positives are unavailable and recorded as `.`; no miniprot PAF, rank or positive fraction is fabricated. Scores retain their engine provenance and are not cross-engine rankings. The original Exonerate `upstream/src/c4/alignment.c` and `upstream/src/model/frameshift.c` informed the split-codon and frameshift handling.

Each comparison directory now includes `normalized_models.tsv`, with structural status, eligibility, identity definition, phases, frameshifts, large query gaps, CDS, translation and QC flags. Eligibility in this audit is diagnostic only: normalized Exonerate models never enter selection, accepted output or padding decisions. An adapter error is recorded in this table independently of coordinate diagnostics; an execution failure is recorded in the existing status files and comparison report.

Stage-2 real-data results and the comparison against saved pre-refactor outputs are in [exonerate-evidence-validation.json](exonerate-evidence-validation.json). Unit tests cover both split-codon phases on both strands, repeated loci, partial query offsets, missing evidence, malformed operations, gaps, terminal/internal stops and frameshifts.

## Earlier selection validation (before strict boundary QC)

The replay tool `rust/gm2_tools/examples/validate_exonerate_selection.rs` builds and selects each engine's models separately from the retained inputs and alignments. It never compares scores between engines. Run it with a shadow output directory and a destination TSV; the adjacent `.models.tsv` includes per-model state, selection, eligibility, coordinates and translation. The current replay uses the default 0.20/0.80 coverage thresholds, matching the recorded dataset runs.

Two regression tests reproduced selector defects before the fixes: a third nearly tied distinct CDS was hidden by two identical leading CDS models, and a translated internal stop without a reported stop count could outrank a clean model. Selection now checks all direct competitors for ambiguity and counts translated internal-stop state as a structural defect. The ambiguity scan is linear within each group. Existing locus grouping and tie thresholds remain unchanged.

| Saved input set | Comparisons | Miniprot accepted | Exonerate accepted | Exact CDS matches | CDS matches after terminal-stop removal |
| --- | ---: | ---: | ---: | ---: | ---: |
| Intact | 56 | 50 | 44 | 2 | 31 |
| Fragmented, including attempted padded joins | 142 | 112 | 101 | 36 | 75 |

Counts are models/comparison calls, not unique genes. A match requires the same protein, strand and overlapping genomic locus. Stop-normalized matching removes at most one translated terminal stop codon and compares the remaining nucleotides; it is not proof of exon-boundary or full candidate-set equivalence. The replay includes only the miniprot-driven padded joins previously attempted, so these results do not validate an independent Exonerate padding workflow.

Concrete intact discrepancies include `ppect_mg004` and `ppect_mg007`: Exonerate begins three nucleotides later and ends three nucleotides earlier than miniprot, so removing a final stop alone cannot equate them. At `ppect_mg010`, the Exonerate model fails the shared noncanonical-splice gate while the miniprot model passes. These observations are backend disagreements, not an accuracy ranking without reference adjudication.

[Selection summary and input hashes](exonerate-selection-validation.json), [intact comparisons](exonerate-selection-intact.tsv), [fragmented comparisons](exonerate-selection-fragmented.tsv), and the adjacent model tables preserve the evidence. The fixed binary was also rerun with and without shadow mode: all 887 authoritative files matched both the no-shadow output and saved pre-fix output; see [pipeline regression](exonerate-selection-pipeline-validation.json). All 293 workspace tests, strict Clippy, and formatting pass. One initial test run encountered an intermittent OS `ExecutableFileBusy` error launching a temporary fake aligner; the subsequent full run passed.

These results do not justify switching the default backend or claiming speed superiority. Start/stop policy, splice discrepancies, missing candidate models and the independent padded-join workflow still require validation.

## Source-reference adjudication and boundary corrections

Miniprot is not the truth oracle. The source GFF, extracted genomic targets and nucleotide references were checked independently of either aligner's model predictions. For all 56 intact candidates, the saved target equals the dataset genomic target and the source GFF CDS extracted from that target equals the reference CDS. These references come from the same annotation used to derive query proteins, so this is annotation concordance rather than independent biological validation.

Both engines exactly reconstruct the source CDS (ignoring one terminal stop) and source intron boundaries for 28 of 56 families. Each has 24 same-strand models missing at least one source microexon of at most 9 nt. Source annotation can itself be imperfect; matching another aligner or passing structural QC is not sufficient evidence of correctness. See [reference audit](exonerate-reference-audit.tsv) and [audit metadata](exonerate-reference-audit.json). The reproducible tool is `scripts/audit_exonerate_reference.py`.

For `ppect_mg004`, the source first exon is 6 nt at 0–6; both engines miss it. Miniprot instead starts at 5010 with CCA (P), despite aligning a reference beginning with M, while Exonerate starts at query residue 1 and target 5013. For `ppect_mg007`, both miss the 4 nt first exon, with miniprot starting at 233 with ACA (T) and Exonerate starting at 236. These are not merely removable terminal stops. For `ppect_mg010`, the source includes a 4 nt internal exon at 10880–10884. Both engines skip it; miniprot reports a canonical alternative with an insertion, while Exonerate reports a noncanonical alternative. Canonical motifs alone do not validate the miniprot model.

The original Exonerate executable was rerun on all 198 saved comparison inputs using the recorded arguments. 197 stdout files match Rust exactly. The remaining case, `ppect_mg046`, differs in a rejected low-coverage opposite-strand hit; its high-scoring source-strand path agrees. Original exhaustive search changes that short-hit result again. At `ppect_mg010`, original exhaustive search retains the same incorrect source-strand microexon-skipping path, so the disagreement cannot be corrected simply by reproducing the original algorithm or removing heuristic search. [Parity hashes](exonerate-original-parity.tsv) and the [boundary validation report](exonerate-boundary-validation.json) retain the evidence. Original `src/model/modeltype.c` constructs protein2genome as a LOCAL model; `protein2genome.c` inserts the phased-intron graph, and `intron.c` applies scored splice transitions. Neither local alignment nor a high score guarantees a complete annotated gene.

The actionable QC defect was an over-permissive completeness label. A model is now `terminal_partial` if its aligned query omits any reference terminal residues, even within the former tolerance. `query_start_truncated` and `query_end_truncated` expose which end is missing. If the aligned query starts at zero, its reference starts with M, and the reconstructed protein does not, `reference_start_mismatch` is recorded and an otherwise complete model becomes partial. Existing structural failures remain failures. No bases or exon coordinates are invented; partial models can still be resolve-eligible under the existing policy, but cannot validate a supposedly complete padded join. `complete` remains a structural/reference-coverage category, not a biological correctness claim; it now requires both reference ends regardless of the configured coverage floor.

On the real dataset this reclassifies 13 intact miniprot models from complete to partial. Fragmented validation accepts 23 padded joins instead of 29; six incomplete joins are rejected and their source fragments retained. Miniprot raw alignments remain byte-identical. The freshly rebuilt workflow, with and without shadow mode, produces 880 byte-identical authoritative files (392 intact and 488 fragmented). Output counts differ from earlier runs because the stricter padded-join gate intentionally changes accepted derived models. All 293 workspace tests, strict Clippy and formatting pass. New selection replay tables are [intact](exonerate-boundary-intact.tsv) and [fragmented](exonerate-boundary-fragmented.tsv).

Microexon reconstruction and dedicated validation are discontinued at the user’s request. The source candidates are predictions from ANNEVO, including low-scoring short CDS features; they are not established biological truth. The historical concordance tables are retained as audit evidence, not as targets for future correction. No exon is automatically restored to match that annotation.

## Integration hardening and performance

Malformed PAF frameshift/stop counts and alignment scores previously fell back to zero or missing values. Invalid numeric tags now produce parse warnings instead of supplying apparently clean alignment evidence, and query/target intervals extending beyond their reported sequence lengths are rejected. A regression test confirms that a valid GFF paired with malformed PAF cannot become resolve-eligible. Missing optional numeric tags retain their existing defaults.

The common model builder also rejects overlapping CDS features on either strand. Previously such features could duplicate coding bases during extraction and still pass structural QC. These guards apply equally to the common evidence path; they do not alter valid alignment scores or exon coordinates. Review of original Exonerate `c4/alignment.c` and `model/frameshift.c` informed preservation of explicit frameshift evidence and monotonic alignment-path semantics; miniprot PAF parsing is a TipSeek-specific input check, not an Exonerate scoring change.

External Exonerate completion polling now starts at 1 ms and caps at 2 ms, reducing the delay after short alignments finish. The same 60-second alignment timeout and 10-second version timeout remain. Longer-running children consequently require more polling wakeups; no CPU or universal throughput improvement is claimed from wall-clock measurements alone. An initial backoff capped at 10 ms produced no meaningful end-to-end gain and was replaced before final benchmarking.

All 296 workspace tests, strict Clippy and formatting pass. The regression tests first reproduced overlapping-CDS acceptance and silent malformed-tag fallback, then passed after the fixes. Timeout, failed-process, reverse-strand and original-output checks remain in the suite. The reproducible benchmark is `scripts/benchmark_exon_shadow.py`; it alternates before/after binaries for three repetitions, verifies executable hashes, and compares all authoritative files plus saved Exonerate stdout and normalized models. See [benchmark commands, hashes and results](exonerate-hardening-benchmark.json).

Final real-data medians (three alternating repetitions):

| Dataset/mode | Before (s) | After (s) | Time reduction |
| --- | ---: | ---: | ---: |
| intact / baseline | 1.566 | 1.566 | 0.0% |
| intact / shadow | 2.568 | 2.367 | 7.8% |
| fragmented / baseline | 3.870 | 3.870 | -0.0% |
| fragmented / shadow | 6.023 | 5.473 | 9.1% |

All 880 authoritative outputs (392 intact, 488 fragmented) match across before/after and shadow-on/off runs; saved Exonerate stdout and normalized model tables also match. This measures the optional comparison workflow, not a change to the aligner algorithm. Baseline fluctuations of this size do not establish a default-workflow speedup.


## Stage 3: independent exon backend

`--gene-exon-backend exonerate-rs` selects the standalone backend; the component accepts `--backend exonerate-rs`. Miniprot remains the default. Exonerate is called by the bare executable name through PATH, with no dependency on miniprot and no Exonerate executable-path option. Selecting it together with shadow mode is rejected before annotation output is removed. Candidate annotation and padded reannotation both use the same selected backend.

Normalized Exonerate evidence is carried explicitly on raw models, with engine provenance propagated into built models. It does not fabricate PAF records. The existing model selector and complete-padded-join gates run on these models. Unknown substitution-matrix positives are exported as `.`. Engine/version and identity definitions are recorded in `manifest/backend.tsv`; raw vulgar records, input snapshots, arguments, stderr and exit status remain under `raw_exonerate/`. Missing executables fail during the version preflight, preserving existing output. A successful empty alignment and a failed process produce distinct candidate states; malformed alignment evidence fails the task rather than being accepted or silently called a no-hit.

The integration test supplies only a fake `exonerate-rs` in PATH, with no miniprot available, and verifies a complete independently validated padded join, CDS output, unavailable-positive export, no-hit/failure separation, malformed-output rejection and preservation of existing output when the executable is missing. CLI tests cover forwarding, invalid backend names and incompatible shadow flags. All 297 workspace tests, strict Clippy and formatting pass.

Real-data validation uses 56 intact families and their fragmented candidates, with independently chosen padded-join attempts for each backend. Three alternating repetitions are deterministic per backend. All 880 miniprot authoritative output files remain identical to the prior implementation. Exonerate uses no miniprot execution or output. Timing excludes upstream read recruitment/assembly and measures annotation of existing candidates on this machine:

| Dataset | Miniprot median | Exonerate median | Speed ratio | Accepted models (miniprot / Exonerate) | Validated padded joins |
| --- | ---: | ---: | ---: | ---: | ---: |
| Intact | 1.566 s | 0.765 s | 2.05× | 50 / 44 | — |
| Fragmented | 3.920 s | 1.515 s | 2.59× | 61 / 61 | 23 / 17 |

The fragmented accepted models cover 54 families for each backend; this does not establish identical sequences or loci. Intact accepted families are 50 versus 44. Exonerate attempted 24 padded joins, while miniprot attempted 30, so the end-to-end speed comparison includes different computational work. Best-path (`subopt=no`) search, different scoring and terminal-stop conventions still limit candidate-set equivalence. Faster execution and equal counts are not accuracy claims. Parallel Exonerate calls and independent biological validation are not part of this stage.

[Commands, binary hashes, per-run times and functional counts](exonerate-independent-validation.json) are retained; reproduce with `scripts/validate_exon_backends.py`. The test process prepended the existing release binary directory to PATH because the desktop shell did not already expose exonerate-rs. No path is passed to the workflow for that executable.

## Further validation

Broaden testing of repeat loci, candidate-set completeness, backend-specific padded joins and resource use before considering any default-backend change. Microexon reconstruction remains out of scope at the user's request.

A representative intact case (`ppect_mg001`) has matching internal exon boundaries but miniprot retains a terminal stop codon; Exonerate also reports a short reverse-strand hit. These are concrete reasons to investigate differences before treating this backend as interchangeable.

## Stage 3 follow-up: admission, padded CDS and process completion

The two backends share QC but produce different paths. The current real-data audit explains eight families with different padded-join outcomes:

- `mg010`, `mg012`, `mg019`, `mg020`, `mg041`: Exonerate source paths contain noncanonical introns and cannot enter padding.
- `mg045`: Exonerate's source path contains an internal stop.
- `mg018`: Exonerate's padded model misses the last reference residue; miniprot's reaches it.
- `mg043`: miniprot's padded model misses the first eight residues; Exonerate's reaches the reference start.

There are 16 joins accepted by both engines. Their CDS sequences are all identical after removing the terminal stop codon retained by miniprot. No other CDS differences remain in these shared joins. This comparison does not establish biological correctness. No microexon work was performed, and no QC threshold was relaxed to equalize backend counts.

Original C Exonerate replay of all 192 current saved inputs matches 191 Rust vulgar outputs exactly. The sole mismatch is the previously recorded low-coverage reverse-strand hit in intact `mg046`; it is not an accepted model. Original `src/comparison/match.c` (`Match_1_3_split_score_func`) translates target triplets and scores them against the query, while vulgar operations describe advances rather than stop-defect tags. Thus an `M` or `S` operation cannot itself establish a stop-free CDS.

A new regression exposed an adapter bug: a stop codon aligned to a normal reference residue at the alignment endpoint was stripped from generic internal-stop checking, allowing a `complete` model. The Exonerate adapter now counts such stops in `M` and split-codon `S` pairs and counts stops in target insertions (`G`). A stop explicitly paired with the reference terminal `*` remains supported; no stop is appended outside the path. Tests cover both strands, split codons, inserted stops, and rejection of a padded model without superseding the clean source fragments. This bug was reproduced with controlled fixtures; it did not change accepted real-data sequences in this dataset.

`manifest/padded_models.tsv` now retains every built padded model's reference interval, state, selection/eligibility, QC flags and all failed padding gates. Previously these rejected models disappeared behind the generic `padded_model_not_validated` reason. The table is diagnostic and does not alter selection. No-hit/build/process failures still use the existing warnings and fragment outcome records.

On Linux, the bounded child runner now waits for a pidfd exit notification rather than timer polling. The existing polling timeout fallback remains for systems without pidfd support. Timeout kill/reap and failure handling tests pass. Five alternating before/after repetitions of independent Exonerate annotation produced:

| Dataset | Before median | After median | Wall-time reduction |
| --- | ---: | ---: | ---: |
| Intact | 0.751 s | 0.707 s | 5.9% |
| Fragmented | 1.491 s | 1.372 s | 8.0% |

Accepted CDS, protein and supercontig outputs are byte-identical before/after. CPU time did not decrease (intact 0.693→0.707 s; fragmented 1.347→1.369 s); this optimization reduces completion latency, not alignment computation. Three independent functional repetitions remain deterministic, including the new diagnostic table. All 880 pre-existing miniprot output files match the previous version. Accepted models remain 50/44 intact and 61/61 fragmented (miniprot/Exonerate); padded joins remain 23/17. All 298 workspace tests and strict Clippy pass.

Reports and reproduction tools:

- [Functional checks and backend timings](exonerate-padding-validation.json): `scripts/validate_exon_backends.py`.
- [Per-family admissions, padded rejection gates, shared CDS comparison and original replay](exonerate-padding-differences.json): `scripts/audit_exon_backend_differences.py VALIDATION.json REPORT.json ORIGINAL_EXONERATE`.
- [Paired performance measurements, CPU times and binary hashes](exonerate-padding-performance.json): `scripts/benchmark_exon_backend_update.py BEFORE AFTER DATASET REPORT.json`.

These measurements cover annotation of the saved real candidates on this machine, not read recruitment/assembly or a universal runtime guarantee. The workflow continues to invoke the bare `exonerate-rs` name from PATH.

## Stage 3 follow-up: incomplete splice signals and output errors

Further review found two reproducible edge-case defects, rather than evidence that the remaining real-data backend differences should be normalized away:

1. An intron only 1–3 nt long was classified `unavailable` and escaped the `noncanonical_splice` rejection. Such a gap cannot contain two non-overlapping dinucleotide splice signals. Models now carry `invalid_splice_length`, fail resolve/padding eligibility, and incur a structural defect during competition. Regression tests cover all three lengths on both strands. Original Exonerate's `src/model/intron.c` models the splice signals separately and defaults `minintron` to 30; this fix does **not** impose that engine-specific default on miniprot, only rejects an incomplete signal pair.
2. FASTA appends used `BufWriter` without checking the final flush. A small write to a failing destination could return success even though the sequence was lost. `/dev/full` reproduces that failure on Linux. FASTA and normalized GFF writers now propagate flush errors, and annotation flushes its model/manifest writers before successful cleanup.

Exon and intron records are now appended in one batch per model and output type, preserving order and the absence of empty intron files. This reduces the number of append opens implied by the real selected models: Exonerate intact 586→110, fragmented 578→154; miniprot intact 636→111, fragmented 601→142. These counts are derived from the selected segment records, not system-call tracing.

All 300 workspace tests and strict Clippy pass. Three repetitions per backend/dataset remain deterministic. All **1,862 pre-existing authoritative files** across both datasets and both backends are byte-identical to the preceding version. The edge-case splice defect is reproduced by fixtures, not observed in this real dataset. Backend admissions remain 50/44 intact, 61/61 fragmented; padded joins remain 23/17. All 16 shared joins remain CDS-identical after terminal-stop normalization.

The per-family audit now includes defective source intron coordinates and donor/acceptor motifs. Remaining paths contain combinations such as GA–AG and GT–TT; these are actual aligned target bases, not a reversed-coordinate extraction error. Original replay still matches 191/192 outputs, with the same unaccepted intact `mg046` short reverse hit as the sole discrepancy. These observations do not establish the biological status of noncanonical splice sites and do not justify relaxing QC.

Five alternating before/after repetitions on the standalone Exonerate workflow:

| Dataset | Before median | After median | Interpretation |
| --- | ---: | ---: | --- |
| Intact | 0.712 s | 0.696 s | 2.2% lower wall time in this run |
| Fragmented | 1.371 s | 1.368 s | Essentially unchanged; 0.3% is too small to claim a reliable gain |

CPU medians were 0.706→0.694 s intact and 1.365→1.365 s fragmented. This is an output-writing optimization; alignment computation and candidate paths are unchanged. Full-pipeline speedups or universal runtime guarantees are not established.

- [Functional regression and per-backend timings](exonerate-splice-io-validation.json).
- [Remaining per-family differences, splice coordinates, padded gates and original replay](exonerate-splice-io-differences.json).
- [Paired performance and derived append-open counts](exonerate-splice-io-performance.json).

Use the same validation, audit and benchmark scripts documented above. Exonerate continues to be resolved by its bare executable name in the current environment; microexon reconstruction remains excluded.
