# Changelog

## Unreleased

## v1.6.2 — TipSeek identity and evidence-bounded UCE recovery

- Renamed the public project to TipSeek and changed the main command from `geneminer2` to `tipseek`; documentation and release links now use the TipSeek identity.
- Decoupled UCE recruitment k from exact-match verification k inside UCEFilter, while leaving the existing component defaults unchanged.
- Changed UCE-mode defaults to k=23/step=4 automatic recruitment followed by one evidence-constrained rescue round. Original assembly keeps k=31/step=4 without rescue; explicit `-kf`, `--uce-recruit-mode`, `--uce-rescue-rounds`, and the new `--no-uce-rescue-reads` override the UCE defaults.
- Strengthened `--uce-recruit-mode auto` by expanding gated fragments against the complete panel before ambiguity and unresolved-locus output gates, parallelizing final probe checks, and preserving all pass and rejection evidence in auditable TSV sidecars and archived outputs.
- Made sensitive-pass contigs explicit provisional cores before rescue. The core gate now audits every fallback locus, rejects a core that already contains an exact long inverted repeat, records internal read-chain gaps as review evidence, and lets only anchored cores without review flags seed the existing whole-contig and terminal rescue path. A gap alone remains review-only because same-individual calibration found a correct low-coverage collinear locus with that signal; the candidate is retained without rescue extension.
- Added the TipSeek UCE manuscript and supplementary methods, including the coral benchmark across TipSeek R0/R1/R2, GeneMiner2, and SPAdes + PHYLUCE configurations.

## v1.6.1 — Automatic UCE recruitment and manuscript draft

- Added the opt-in `--uce-recruit-mode auto`, which retries only unresolved UCE loci with a sensitive k=21/step=1 recruitment gate and retains auditable pass and rejection evidence.
- Renamed the public project to TStools while retaining the `geneminer2` entry point and historical compatibility identifiers.
- Added the initial Chinese UCE manuscript draft for collaborative revision and provenance tracking.

## v1.6.0 — Rescue QC and adaptive mitochondrial completion

- Added per-locus UCE rescue guards that revert newly introduced long inverted repeats and unsupported internal gaps.
- Added explicit one- or two-round rescue control and corresponding per-locus status reporting.
- Fixed mitochondrial adaptive checkpoints so deeper read budgets continue until the same validated circular genome is recovered at two consecutive depths.
- Fixed the nightly FASTX fuzz workflow and expanded regression coverage.

## v1.5.9 — UCE scheduling and rescue correctness

- Made the default shared CPU budget `-p auto`: it detects physical cores permitted by affinity or cpuset constraints, then caps that value using cgroup or scheduler limits. Explicit integer values still override auto-detection, while UCEFilter's single supported compute worker remains separate from recruitment's 1--2-unit reservation.
- Hardened UCE raw-read rescue: rescue assembly now uses fixed k=21, preserves only active loci in the rescue reference, and retains evidence from prior rounds. Every terminal extension is independently checked for length, breadth, gap, fragment, and bridge support; failed extensions revert without discarding accepted evidence.
- Added regression coverage for automatic CPU budgeting, fixed-k rescue assembly, active-locus rescue references, terminal-evidence acceptance, cross-round evidence preservation, and per-locus rescue-status reporting.

## v1.5.7 — UCE rescue I/O and PanRefV2.2 graph evidence

- Avoided byte-for-byte sample-directory copies during UCE rescue. Rescue now stages work outside the sample directory and atomically moves the settled round back, preserving rollback and output layout while reducing slow-disk I/O.
- Extended PanRefV2 with quality-masked bounded read ledgers and sequential candidate spools, sparse per-sample unitig evidence, conservative complete sample-backbone paths, stable backbone IDs/coordinates, and explicit bubble QC.
- Resolved acyclic PanRefV2 graphs globally using only observed accepted-read transitions, preserved GFA path orientation, added an adaptive k=25 retry for supported unresolved graphs, and made unitig-edge construction near-linear.

## v1.5.4 — UCEFilter FM-index performance pipeline

- Replaced positional-anchor verification with per-locus FM indexes that report maximal exact matches, retaining repeat-aware orientation handling while substantially reducing UCEFilter index memory.
- Added bounded parallel R1/R2 gzip decode, with zlib-ng used automatically when available; decode overlaps recruitment without changing paired-read ordering or filtering output.
- Compacted evidence and per-locus candidate storage, added memory reporting, and updated the UCE scheduler to budget one compute worker plus two decode workers per sample.
- Validated byte-identical full Octocoral output against the prior storage implementation, alongside Rust and workflow regression suites.

## v1.5.3 — Conservative mitochondrial rescue efficiency

- Kept the existing global adaptive-stop semantics while avoiding a read de Bruijn graph when no terminal mate link reaches the established support threshold; this cannot remove an otherwise admissible component bridge.
- Made mitochondrial rescue seeds reverse-complement-aware and auditable. Contigs containing ambiguous bases contribute their unambiguous ACGT segments; duplicate and non-informative low-complexity seed segments are excluded without using reference similarity as a filter.
- Added `mito_rescue_seeds.tsv` provenance and regression coverage for seed decisions and the read-graph short circuit.

## v1.4 — MainFilter I/O and gzip backend acceleration

- Kept `MainFilterNew` output file handles open for the life of the run instead of reopening on every buffer flush, and raised the process file-descriptor limit at startup to accommodate large per-locus output sets.
- Replaced `String`-based FASTQ/FASTA line reading with byte-level parsing, removing per-line UTF-8 validation and a redundant copy; enlarged the gzip and file input buffers to cut read syscall frequency.
- Added runtime detection of zlib-ng via `dlopen`/`dlsym`: uses its SIMD-accelerated gzip decompression when present in the environment, and falls back transparently to system zlib otherwise with no build-time dependency change.
- Verified byte-identical filtering output across all output modes and both gzip backends on real UCE target-capture data; documented measured gains in `docs/development/mainfilter-performance.md`.
- Added a joint mito rescue round: after the first UCE-style assembly pass, all retained contigs become sample-specific seeds that are combined with the GenBank baits into one rescue reference, then recruited and reassembled together with the original paired reads.
- Replaced the near-identical-consensus adaptive-stop heuristic with an exact, cut- and strand-independent circular sequence comparison; adaptive stages now retain partial assemblies across read-depth increases instead of requiring circularity at every stage, and reuse the immutable GenBank-derived reference cache across stages.
- Added unit tests for the joint rescue reference builder, the exact circular comparator, and the adaptive stage state machine.

## v1.3.2 — Rust consensus generation

- Replaced the production `build_consensus` executable with a Rust implementation while preserving its command name and the existing consensus stage in the main CLI.
- Added SAM and gzip-SAM parsing, CIGAR-aware base/deletion/insertion handling, primary-alignment filtering, multi-reference IUPAC consensus FASTA output, and compatible `.sam.gz` default prefixes.
- Replaced the optional Python matplotlib/SciPy mutation-density plot with a portable Rust bitmap PNG implementation; the consensus production path no longer requires Python plotting dependencies.
- Added Rust unit, gzip/multi-reference/indel/primary-alignment, PNG, and legacy-output compatibility regression tests.

## v1.3 — MainFilter canonical index and safe cache

- Made Rust `MainFilterNew` use canonical k-mers for strand-invariant exact recruitment, retaining the `-gr` flag only as a compatibility alias. Recruitment and per-locus output remain equivalent for identical inputs and parameters.
- Replaced inline/shared multi-locus hit storage with compact packed locus postings, reducing the memory cost of multi-reference dictionaries.
- Added reusable dictionary format v3 with canonical-policy metadata and a SHA-256 fingerprint of reference content. Stale, legacy, or mismatched dictionaries are rebuilt rather than silently reused.
- Added cache-invalidation, randomized lookup, end-to-end, and old-versus-new byte-equivalence coverage for MainFilter outputs.

## v1.2 — Gene QC hardening and documentation refresh

- Added two-pass QC to `gene-resolve`: pre-alignment filtering by translated candidate length and distinct-sample occupancy, followed by post-alignment occupancy and effective-codon-site checks. The new `occupancy_qc.tsv` records both decisions and rejection reasons per family.
- Added `--gene-min-aa-length` and `--gene-min-effective-codon-sites`; occupancy is now explicitly counted by distinct samples, so multiple candidates from one sample cannot inflate a family’s taxon coverage.
- Modernized the Chinese and English README and workflow guides with task-oriented navigation, concise input/output/boundary summaries, dedicated bilingual mitochondrial guides, and a testing-stage support notice.
- Removed two unused legacy Python utilities.

## v1.1 — Gene-family recovery and conservative repeatome analysis

- Added the `gene` workflow for multi-species bait-defined nuclear gene families: candidate recovery, protein-guided miniprot annotation, conservative strict one-to-one clade resolution, and strict or multicopy ASTER2 species-tree inputs with provenance and QC.
- Added the standalone `te` workflow for reference-free short-read repeatome analysis: taxon-balanced discovery, exact-equivalence curation, read-supported conservative annotation, and candidate-read RPM quantification without repeated MainFilter runs.
- Added explicit, reproducible manifests and bounded evidence rules for both workflows; TE annotation never merges repeat equivalence groups and reports weak evidence as unresolved rather than forcing a TE family.
- Added synchronized Chinese and English README, manual, and dedicated Gene/TE workflow documentation; the release build now includes the new Rust executables.

## v1.0 — Reference-level profiling, pan-reference population analysis, and mitochondrial workflow

- Reworked marker profiling into a general reference-level evidence workflow: GeneMiner2 recruitment plus Themisto pseudoalignment now reports per-reference hit counts, fractional shared-query support, and singleton support. Removed the mSWEEP dependency, group-abundance output, and associated CLI options; an optional group map is annotation only.
- Added the experimental `panref` population engine, which builds a graph-backed cohort reference from UCE baits and sample reads while retaining the existing pseudo-reference engine and staged population outputs.
- Added the `mito` workflow for annotated GenBank references: one mitochondrial read pool, UCE-style assembly, overlap merging, mate-link/read-graph validation, and conservative circularity reporting.
- Reorganized Chinese and English READMEs, manuals, profiling/population chapters, and mitochondrial documentation; added ignore rules for generated Cython artifacts.

## v0.9.3 — Deterministic MainFilter I/O optimization

- Avoided retaining FASTQ headers, `+` lines, and duplicate normalized text buffers in default GM2 output and scan-only modes; text-output modes remain byte-compatible.
- On the DK40 target-capture benchmark (one million read pairs), reduced default GM2 filtering time by about 7–9% while preserving all 4,466 GM2 files and the read-count report byte-for-byte.
- Polished the MainFilter performance note with explicit compatibility boundaries, benchmark scope, and release-validation requirements.


## v0.9.2 — Four-chapter documentation

- Reorganized user documentation into Filter, Assembler, Profiling, and Population chapters in both Chinese and English.
- Consolidated UCE workflow and assembler rationale into the Assembler chapter; added the first dedicated generic marker-profiling chapter.
- Reduced duplicated workflow and QC prose in the command-line manuals while retaining command and option references.
- Moved the MainFilter performance note to `docs/development/` and removed superseded workflow documents.

## v0.9.1 — Marker profiling hardening and mode clarification

- Made marker profiling fully group-map driven: dynamic reporting groups, exact reference-to-group coverage checks, content-addressed Themisto cache keys, safer output handling, and expanded QC.
- Made `--profile-kmer-size` apply consistently to both GeneMiner2 recruitment and Themisto pseudoalignment; compute immutable profiling cache inputs once per run rather than once per sample.
- Renamed the public conventional assembly mode from `reference` to `original`. `original` is for exon, SCO, and nuclear or mitochondrial marker recovery; `uce` is for UCE recovery from genome-skimming or target-capture data.
- Updated Chinese and English READMEs, manuals, output descriptions, and assembler documentation for profiling and the `original` / `uce` split.


## v0.8 — Original-Rust default and assembler validation

- Made `original-rust` the `reference + auto` default and renamed the user-facing UCE Rust implementation selector from `rust` to `uce-rust`.

- Documented a fixed-parameter, 40-locus single-thread comparison between the upstream Python assembler and `main_assembler-original-rust`: identical locus status, 38/39 identical best-contig sequences, and a documented remaining difference at `v1__uce-1200`; reference mode now defaults to `original-rust`, while the upstream Python implementation remains available as `original` for strict comparison.

- Restored the byte-identical upstream GeneMiner2 Python assembler for reference-mode fallback and reproducibility.
- Removed the UCE-aware Python fallback source, executable, CLI option, build target, and mode routing; UCE and ITS2 now fail directly when Rust assembly is unavailable or fails.
- Made the upstream original assembler the direct default for reference mode; Rust reference assembly now requires explicit selection.
- Added a Chinese algorithm note comparing the upstream and Rust GeneMiner2 assemblers and separating the contributions of MaSuRCA, SPAdes, and Sparrowhawk from features not adopted.
- Added a versioned binary reference k-mer cache for `main_assembler-original-rust`, with reference identity validation, corrupt-cache rebuilds, and atomic replacement.

## v0.7.2 — Documentation structure and readable Rust internals

- Reorganized the Chinese and English READMEs around mode selection, installation, quick start, and primary outputs.
- Added synchronized UCE and Population workflow guides covering assembly guardrails, rescue fallback, pseudo-reference validation, staged execution, and required QC.
- Added concise Chinese Northeast-dialect comments throughout the Rust MainFilter, Refilter, and Assembler without changing behavior.

## v0.7.1 — Repository cleanup

- Removed the unreferenced population pseudo-reference comparison helper.
- Made `clean` remove Python bytecode caches and made `distclean` remove all generated PyInstaller spec files.

## v0.7 — ITS2 assembly and Rust utility migration

- Added ITS2 multi-candidate assembly with paired-fragment compatibility, equivalence groups, diagnostic support, and EM abundance estimates; ITS2 now remains strictly Rust-only on failure.
- Reimplemented alignment cleanup, sequence merging, reference trimming, and UCE statistics as readable Rust utilities while preserving their command-line contracts.
- Removed the unused MUSCLE integration and obsolete validation-only helper scripts.
- Fixed single-end statistics, rescue scheduling after sample failures, deterministic gene-tree ordering, and top-level CLI error handling.
- Synchronized Chinese and English README, command-line, and output documentation with the current CLI; removed obsolete GUI-era console output and local debug artifacts.

## v0.65 — MainFilter deterministic lookup optimization

- Optimized the Rust primary filter's short-k-mer scan with a DNA lookup table, modulo-free probe scheduling, and `AHashMap` k-mer lookup.
- Kept filtering semantics, command-line options, and cache/output formats unchanged; documented byte-level output verification and the decisions not to add threads, LRU output handles, or low-gain hash alternatives.

## v0.6 — Scalable UCE graph assembly

- Stream filtered reads in bounded batches and count k-mers through parallel, sorted per-batch aggregation.
- Compress non-branching UCE backbone paths into unitigs, retaining bounded decisions only at graph junctions.
- Add optional compact GFA and DOT assembly-graph output via `--assembler-graph-format`.
- Add `--assembler-read-chunk-size` and `--assembler-kmer-count-threads`, while preserving the unmodified Python fallback.

## v0.5 — Rust UCE assembly and reusable population analysis

- Added the high-performance Rust UCE assembler with compact rolling k-mers, a bounded non-backtracking backbone path strategy, reference caching, and parallel per-locus assembly.
- Made Rust assembly the default via `--assembler-implementation auto`; failed or unavailable Rust runs now clean incomplete outputs and retry the unmodified Git-baseline Python assembler.
- Retained strict `uce-rust` and direct `original` assembler modes for reproducibility and diagnosis.
- Added fixed external cohort-reference support, checked resume stages (`mapping`, `calling`, and `selection`), and per-stage variant-count QC to the population workflow.
- Added reusable tools for summarizing UCE validation runs and comparing a population pseudo-reference with an external reference.
- Updated Chinese and English command-line/output documentation and regression coverage.

## v0.4 — Population analysis

- Added the Rust `population` workflow: cohort-reference construction, uniform minibwa mapping, joint bcftools variant calling, and one representative SNP per UCE.
- Added SqCL-inspired longest-eligible-contig reference selection, with a read-support-first alternative and per-sample reference-contribution diagnostics.
- Added all-SNP, one-SNP-per-UCE, and LD-pruned VCF/PLINK panels with PCA for each panel.
- Added automated ADMIXTURE K-range analysis, cross-validation summaries, status reporting, and retained logs.
- Added mapping-rate, coverage-breadth, depth, sample-name, and reference-provenance quality-control reports.
- Added real-tool integration tests for minibwa, samtools, bcftools, PLINK, and ADMIXTURE.
- Updated Chinese and English usage and output documentation.

## v0.3 — Rust primary filter

- Reimplemented the primary read filter in Rust while retaining the original command-line and cache compatibility behavior.

## v0.2 — UCE assembly validation and rescue

- Added UCE assembly guardrails, read-support validation, and controlled raw-read rescue.

## v0.1 — UCE workflow foundations

- Added sequence-integrity fixes and the initial UCE-focused command-line workflow.
