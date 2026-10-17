# Exon/intron annotation

[中文](exon_ZH.md) · [gene workflow](gene_EN.md) · [command-line guide](../manual/EN_US/command_line.md)

`--assembly-mode exon` extends the default gene-family recovery workflow with protein-guided exon/intron annotation. By default it derives protein references directly from coding nucleotide baits, so separate `.faa` files are not required; external proteins remain an optional override. This mode is not for UCEs.

## Inputs and command

Each `.fa`, `.fas`, or `.fasta` file under `-r` defines one family and may contain coding nucleotide baits from several species, such as CDS or concatenated exons that preserve coding frame.

```text
family_reference/
├── geneA.fasta
└── geneB.fa
```

These nucleotide references are translated automatically with the standard nuclear genetic code. If a bait contains introns, is non-coding, uses a non-standard genetic code, or cannot be interpreted as coding sequence, supply an external protein directory instead. Each override must use the matching family stem and `.faa` extension. The directory may be incomplete: families with `.faa` files use them, while the remaining families still use automatic translation.

Run exon mode as a complete workflow; do not add stage subcommands such as `filter` or `assemble`.

```bash
cli/tipseek --assembly-mode exon \
  -f samples.tsv \
  -r family_reference \
  -o exon_output \
  -p 8
```

To override selected families with external proteins:

```bash
cli/tipseek --assembly-mode exon \
  -f samples.tsv -r family_reference \
  --gene-protein-reference family_proteins \
  -o exon_output -p 8
```

Annotation requires miniprot 0.18 or later. TipSeek first writes recovered candidates to `exon_output/gene/`, then writes structural results to `exon_output/exon/`.

## Method

1. The default gene workflow recruits, refilters, and assembles reads, then summarizes candidates across samples.
2. Without a matching `.faa`, TipSeek examines all six reading frames equally. A terminal stop is allowed and removed from the protein reference; any frame with an internal stop is excluded rather than repaired or masked.
3. A sequence is accepted directly only when it has one stop-free frame or one frame has uniquely strongest CDS endpoint support from an initial methionine and/or terminal stop. The longest such translation becomes the family anchor. Remaining ambiguous sequences are accepted only when one frame has a unique positive amino-acid 3-mer Dice match to that anchor. With no trustworthy anchor, a tied match, or no positive match, TipSeek emits no protein for that sequence and records `ambiguous_reading_frame`; a matching `.faa` is then required if the family has no other usable reference. Every candidate and decision is recorded in `manifest/reference_translation.tsv`.
4. Miniprot runs independently on each candidate contig. This prevents one partial candidate from suppressing a complementary candidate through target competition.
5. Embedded PAF supplies protein coordinates, coverage, score, identity, protein CIGAR, frameshifts, and in-frame stops. GFF3 supplies exon/CDS coordinates, strand, and phase.
6. Overlapping models compete deterministically by structural defects, protein coverage, alignment score, identity, splice evidence, and miniprot rank. Non-overlapping models remain independent.
7. Selected models are classified, written to structured audit tables, and routed to resolve-eligible or unresolved output.

`N` and other IUPAC bases are preserved in genomic outputs, so reported coordinates always refer to the candidate sequence used for annotation.

## Model acceptance

Only selected `complete`, `terminal_partial`, and `internal_partial` models can enter `cds/` and later `gene-resolve`. An eligible model must have:

- no frameshift or internal stop;
- explicit, continuous CDS phase;
- CDS length divisible by three;
- translation without `X`;
- no unrecognized noncanonical splice;
- no unresolved model competition.

A single terminal stop codon is retained in the annotated CDS. It and the corresponding terminal `*` are removed together only when protein and codon alignments are prepared. Other states—including `low_coverage`, `frameshifted`, `internal_stop`, `phase_inconsistent`, `ambiguous_cds`, and `ambiguous_model`—remain under `unresolved/` and in the manifests.

## Short exons and multi-contig recovery

Exon mode searches remaining reference seeds and recovers unbranched paths from supported, nucleotide-unanchored read-graph components, retaining candidates of at least 30 nt; ordinary gene mode retains its existing behavior. Miniprot uses `-L 10 -n 2 --outc 0` to retain short evidence. Low-coverage fragments are not thereby accepted for resolve.

Reference-ordered chains can include two or more contigs and internal fragments covering at least 10 aa with ≥80% identity and no structural defects or query gaps. The same contig cannot re-enter a chain. Adjacent identities must differ by no more than 15 percentage points. Overlapping CDS must agree base-for-base at the reference coordinates; there is no fixed overlap-length cap. Indel-bearing overlaps without a reliable coding-coordinate correspondence are not trimmed. Near-equal alternative chains are ambiguous. Search is bounded to 128 nodes per protein and 64 paths per endpoint; truncated searches remain unresolved. Equivalent contig orders inferred through different reference orthologs are deduplicated.

`fragment_cds/` contains explicitly unvalidated reference-guided CDS hypotheses, excluded from automatic resolve. Missing internal coding residues are represented by three Ns per amino acid; these are not intron padding. Only coordinate-compatible CDS fragments are exported.

For a unique chain meeting coverage and terminal requirements, the selected backend reannotates oriented source contigs separated by `--gene-fragment-padding` Ns (default 100). Acceptance requires one complete forward model, every padding interval wholly inside an intron and outside exons, and coding support from every source fragment. An alignment that skips the middle source is rejected. Accepted CDS has no synthetic coding bases; genomic padding remains in supercontigs/intron outputs. Missing splice flanks or unresolved copies can still prevent a join.

`manifest/fragment_chains.tsv` records all source members, reference coverage, coding gaps and trimming. `fragment_groups.tsv` retains first/last endpoints. `physical_linkage=unverified` applies even after structural reannotation: inferred order does not prove physical linkage, copy identity or the real intron length. Padding set to zero disables stitching output while retaining diagnostics.

Coverage now subtracts query-only gaps (`I` in miniprot; query-only `G` in Exonerate). Models spanning both reference termini but missing internal reference residues are `internal_partial`, not complete, and cannot pass complete padded-join validation. This is a reference-relative distinction, not proof of a biological exon deletion.

## Paired-read audit and unanchored rescue

Only new unanchored components (`assembler_seed_position=-2`) undergo a six-frame protein screen: one ungapped 10-aa window with at least eight matches and at least four distinct residues in each window. Rejected candidates remain in the input and receive `protein_screen_rejected`; anchored candidates bypass this screen. This conservatively excludes isolated low-complexity or shorter continuous coding evidence, not proof that it is noncoding.

Annotation automatically reads the workflow's `sample/filtered_pe/family_1,2` files. `--gene-reads-root` selects another workflow directory (underlying annotate option: `--reads-root`). Paired GM2 and FASTQ, including gzip, are supported. FASTQ names and counts must agree; GM2 lacks original names and relies on the producing workflow's preserved mate order. Missing files are explicitly unavailable; missing one mate or inconsistent pairing is an error.

A family-wide index audits original contigs once. Only unique, full-length ungapped placements of reads ≥90 nt with ≤4 substitutions and a ≥2-substitution advantage over the runner-up count. Repetitive or overflowing searches remain unassigned. Coordinate-identical pairs are deduplicated. Both ends must be entirely exonic and inward after orienting the source fragments. Two distinct low-mismatch pairs give `paired_read_support`, which is not proof of physical linkage. No bridge is merely insufficient evidence, particularly across long introns.

Following HybPiper's asymmetric mismatch principle, at least five distinct pairs with one exact mate and another having 3–4 substitutions, exceeding low-mismatch support, trigger `discordant_read_warning`. Such chains are withheld from CDS hypotheses and automatic padding. This limited audit does not model base quality, indels or insert length; it is not a general-purpose read mapper or a false-merge-rate estimator. Orientation conflicts are reported separately.

Both backends use a GeneMiner2-style SCO task queue. `-p` controls concurrent SCO workers (at most 64); a worker handles one entire SCO and immediately claims another on completion. Individual alignment processes remain single-threaded. As the SCO queue drains, remaining families can borrow idle slots for candidate jobs within the shared CPU budget; the two concurrency limits are not multiplied. Workers own isolated directories; outputs merge in stable sample, family and candidate order. `manifest/execution.tsv` records the scheduling mode, and the log reports completed SCO tasks.

See `read_audit.tsv`, `fragment_read_support.tsv`, and `completion_audit.tsv` under `manifest/`. Internal join flanks are limited to 120 observed bases around each model, retaining outer termini; `join_sources.tsv` records the retained original coordinates. Full structural reannotation remains mandatory.

## Main options

| Option | Default | Meaning |
| --- | ---: | --- |
| `--gene-protein-reference` | auto-derived | Optional `.faa` directory; matching families override automatic translation |
| `--gene-miniprot` | `miniprot` | miniprot executable |
| `--gene-max-intron` | `50000` | Maximum intron length accepted by miniprot |
| `--gene-min-model-coverage` | `0.20` | Minimum protein coverage for a normal partial model; padding trials alone allow a 1-aa discretization margin |
| `--gene-complete-coverage` | `0.80` | Minimum coverage for a complete model; both protein ends are also required |
| `--gene-fragment-padding` | `100` | Ns inserted between contigs for unique-chain validation; `0` disables joining |
| `--gene-reads-root` | workflow output | Workflow root containing per-sample recruited paired reads |
| `--gene-flank` | `0` | Observed bases added on each side in `genes_flanked/`; never N-padded |

## Outputs

```text
exon_output/
├── gene/                         # recovered candidates and family summaries
└── exon/
    ├── cds/                      # resolve-eligible CDS; no padded N in exons
    ├── proteins/                 # one-to-one with cds/
    ├── exons/                    # one record per supported CDS exon
    ├── introns/                  # one record per predicted intron
    ├── genes/                    # observed single-contig gene spans
    ├── genes_flanked/
    ├── supercontigs/             # selected spans; validated joins may contain intronic N
    ├── gff3/                     # normalized gene/mRNA/exon/CDS/intron hierarchy
    ├── raw_miniprot/             # per-candidate output and padded validation runs
    ├── models/
    │   ├── gene_models.tsv
    │   └── gene_segments.tsv
    ├── manifest/
    │   ├── reference_translation.tsv # automatic frame/strand selection audit
    │   ├── derived_proteins/     # automatically derived proteins actually used
    │   └── ...                   # candidate, ID, warning, fragment, and command audits
    └── unresolved/
```

Use `manifest/reference_translation.tsv` for reference source and frame selection, `models/gene_models.tsv` for model-level status and eligibility, `models/gene_segments.tsv` for coordinates and splice classes, and `manifest/fragment_groups.tsv` for every attempted padded join. Continue with:

```bash
cli/tipseek gene-resolve --gene-input exon_output/exon -o gene_resolved -p 8
```

See the [gene workflow guide](gene_EN.md) for resolve QC and strict or multicopy species-tree inference.

## Optional Exonerate comparison (stages 1–2)

Pass `--gene-exonerate-shadow` in exon mode to compare each candidate and each attempted N-padded sequence using `exonerate-rs` from the current environment’s PATH; no executable path is required. The option is off by default; Exonerate is not bundled. Miniprot and the existing QC remain the sole source of accepted models, CDS, and padding decisions.

The comparison runs `protein2genome --exhaustive no --subopt no --score 1 --revcomp yes`, using `--gene-max-intron` and a 60-second timeout per alignment. Failed, timed-out or invalid runs are audited as failures, not no-hits. Workflow I/O errors still propagate.

`manifest/exonerate_shadow.tsv` records raw miniprot model counts, Exonerate hit counts, matching counts, scores, protein/genomic intervals, CDS intervals and frameshift operation counts. Intervals are zero-based, half-open and genomically ascending, including reverse-strand hits. `stage` distinguishes `candidate` and `padded_join`. `raw_exonerate_shadow/` retains version output, arguments, input snapshots, stdout/stderr, status and elapsed time after temporary workflow files are removed.

`reported_hits_match` means every reported Exonerate hit has a miniprot raw model with the same protein interval, strand, genomic span and CDS intervals. It does not establish equivalent candidate sets, CDS phases, identity, structural QC or accepted models: this Exonerate mode reports only best paths, while miniprot can report multiple models. Frameshift operation counts are not assumed to have miniprot's semantics. Coordinate differences need investigation and do not identify a winner. The stage-2 adapter also writes `normalized_models.tsv` in each comparison directory, applying the shared structural QC to reconstructed CDS phases and translated proteins. Identity explicitly counts matching aligned M/S amino-acid pairs, excluding gaps; substitution-matrix positives are unavailable (`.`). These diagnostic models do not enter selection. Independent backend selection remains a later stage.

Complete status now requires both reference-protein ends to be covered. Even a small terminal omission produces `terminal_partial` and `query_start_truncated` or `query_end_truncated`. A reconstructed non-M start aligned to the origin of an M-starting reference is marked `reference_start_mismatch`; such partial models cannot validate a complete padded join. Complete is a structural/reference-coverage category, not proof of biological correctness, and miniprot is not treated as ground truth.

## Independent Exonerate backend

Add `--gene-exon-backend exonerate-rs` to an exon command to use Exonerate for candidate annotation, model selection and padded-join reannotation. The executable is resolved as `exonerate-rs` through the current PATH, without an executable-path argument or a miniprot dependency. Miniprot remains the default. This option cannot be combined with `--gene-exonerate-shadow`.

```bash
cli/tipseek --assembly-mode exon -f reads.tsv -r references -o exon_output -p 8 \
  --gene-exon-backend exonerate-rs
```

For existing candidates, the component accepts `gene_workflow annotate --backend exonerate-rs --input gene_output/gene --nucleotide-reference references --out exon_output/exon`.

Both backends use shared structural QC and padding gates. `raw_exonerate/` preserves vulgar output, input snapshots, arguments, exit status and stderr; `manifest/backend.tsv` records engine/version and identity semantics. Standard model, CDS, protein, GFF3 and resolve input layouts are preserved. Exonerate identity counts aligned M/S amino-acid pairs; unavailable substitution-matrix positives are exported as `.` rather than fabricated PAF evidence.

The current mode uses `protein2genome --exhaustive no --subopt no --revcomp yes`, a 60-second per-alignment timeout and bounded SCO-level parallelism (`-p`, at most 64 concurrent SCO tasks); each Exonerate process remains independent. An unavailable executable fails before old output is removed. Failed candidates are marked `exonerate_failed`, distinct from no-hits; malformed output fails the task. Best-path search does not guarantee recovery of repeated loci, and candidate sets, scores and terminal handling are not interchangeable with miniprot. Timing and accepted-model counts are not accuracy measures.

Per-model padded reannotation diagnostics are retained in `manifest/padded_models.tsv`, including reference intervals, model state, QC flags and failed join checks. These distinguish missing terminal residues, splice defects and failures to bridge the N padding. Backend CDS outputs may differ by a retained terminal stop codon (3 nt); a file mismatch alone does not establish a splice difference.


## One recovered sequence per SCO

`exon/sco/<SCO>.fasta` contains one sequence per sample in a fixed reference coding coordinate system. Missing bases and unresolved disagreements are N; reference bases are never copied into gaps. Structurally clean short exon evidence survives whole-gene competition. This is an unphased evidence composite, not a validated single-copy CDS or physical cross-contig linkage.

Use `--gene-sco-reference reference_id` (component: `--sco-reference reference_id`) to choose a reference record within each SCO. Matching must be unique. Without it, the backbone is selected across the entire cohort by identity-weighted unique residue support, with a penalty for atypical reference lengths. Duplicate or contained hits cannot inflate this score. Evidence against other references is projected only through unambiguous, locally anchored protein coordinates; insertions and repeat-shift alternatives remain unmapped. An original sample codon assigned to multiple backbone positions masks those positions as N. Both backends retain residue-to-CDS maps, including query gaps, target insertions and split codons. Inserted sample codons are omitted from the fixed-coordinate composite but retained in original models. A terminal-stop coordinate present in the translated nucleotide reference is filled only from an observed sample stop, otherwise NNN. Automatically reverse-translated or shifted references use the selected coding translation's coordinates; consult `reference_translation.tsv`.

`manifest/sco_recovery.tsv` reports observed, missing and conflicting bases and unmappable models. `exon_candidates.tsv` and `exon_candidate_sources.tsv` preserve coordinates, source models, boundary evidence and read coverage without separate candidate-exon FASTA exports. Individual uniquely mapped mates support exon sequence even when their partner is unmapped; they do not establish linkage. Legacy model/GFF/exon outputs remain diagnostic; use `sco/` for the requested per-SCO result.


Reuse saved alignments without invoking either backend:

```bash
gene_workflow recover-sco --annotation output/exon --input output/gene \
  --out recovered_sco --threads 16
```

Add `--sco-backbone-manifest output/exon/manifest/sco_recovery.tsv` to preserve each original backbone for a coordinate-matched control. Recovery reads the original candidate DNA, saved alignments and cached reference proteins, and fails on missing required caches. It does not recruit reads, assemble contigs or rerun alignments. `manifest/reference_projection.tsv` audits coordinate mappings; `sco_selection_control.tsv` reports the reference-selection-only control. Use the candidate DNA associated with the original annotation.

Redundant evidence across contigs is compared only when an entire observed contig (at least 30 nt) occurs exactly once in another observed contig at least 30 nt longer, including reverse-complement containment. Inferior models are removed only when the existing identity, informative extra context and complete physical-codon containment requirements also hold. Base variants, repeat placements and near ties remain unresolved. This does not establish same-copy linkage.

Read support allows N in reads but counts every N as a mismatch, never a deletion or wildcard. Full-length placement, uniqueness and the runner-up edit margin remain required. Exact containment of a short contig inside a read supports only its local sequence, not a paired link.

Add `--reads-root output` to cache replay to recompute paired-read chain evidence in `chain_audit/<sample>/<SCO>.tsv`, including ambiguity, search limits and supported edges. This diagnostic does not insert chain hypotheses or unknown bases into the SCO composite. Normal annotation plans chains before global model competition; repeated planning on previously selected models recomputes competition within each reference.

When protein alignment omits the stop, SCO recovery can inspect the immediately adjacent in-frame triplet on the same original contig. This requires an exact, unique ten-residue terminal anchor containing at least five distinct amino acids, consecutive projected coordinates and no terminal coding gap. Only observed stop codons are used, only if the nucleotide reference contains that coordinate; disagreements remain N. There is no downstream search, cross-contig stop inference or reference-base filling.

An observed 30–59 nt candidate failing the ordinary protein screen can receive a supplemental exact eight-residue screen when at least two coordinate-deduplicated, uniquely placed reads jointly cover its entire contig. The anchor must contain at least five distinct amino acids; repeat placements and conflicting strand/frame evidence fail. Passing only permits backend alignment: existing model and projection checks still apply. This rescue is disabled without reads, does not infer unobserved microexons, and never treats short-contig coverage as paired linkage.

Chain filtering stops comparing later paths once overflow of the retained prefix is certain. Candidates are sorted and survivors are append-only, so the first 64 paths, search-limit flag and ambiguity decisions are unchanged. This optimization does not relax conflict or splice checks.
