# 2. Assembler

[中文版本](assembler_ZH.md)

Assembly turns recruited locus reads into contigs. Choose the workflow by biological target, not by sequencer or library label.

| Workflow | Appropriate target | Default backend | Result |
|---|---|---|---|
| `--assembly-mode original` | Exons, SCOs, nuclear markers, and mitochondrial markers | `original-rust` | Reference-guided contigs, trimmed by default |
| `--assembly-mode uce` | UCEs from genome-skimming or target-capture data | `uce-rust` | UCE cores plus read-supported flanks |

## Backends

`original-rust` is the deterministic Rust compatibility backend for `original`; `original` remains an alias for compatibility. `uce-rust` is the only backend for `uce`.

All backends use reference-positioned seeds, read-k-mer support, bidirectional extension, and read-slice validation.

## Original workflow

The original algorithm extends the highest-weight edge first and keeps alternatives on a stack for backtracking. It retains up to three candidates per side before combining and scoring them. This is appropriate for shorter reference-guided targets such as exons and organellar markers.

```bash
cli/tipseek -f samples.tsv -r references -o output -p 8 \
  --assembly-mode original
```

`original-rust` can reuse a versioned binary reference-k-mer cache with `--reuse-reference-cache`.

## UCE read routes: `ucefilter` and main + re

The default UCE route is:

```text
raw paired FASTQ → ucefilter → filtered/ → uce-rust → contig
```

`ucefilter` is the UCE-specific fused read stage. In one scan of the original FASTQ it performs broad recruitment, complete paired-fragment retention, orientation and maximum-exact-match evidence, and per-locus depth/position selection. It therefore replaces the generic `MainFilter → filtered_pe/ → refilter → filtered/` (**main + re**) route for default UCE. Use `filter assemble`; a separate `refilter` stage is neither needed nor valid afterwards.

Only `--legacy-uce-filter` restores the main + re compatibility route, chiefly for historical comparison or diagnosis. There `filtered_pe/` is broad recruitment and `filtered/` is the assembly input. Default `ucefilter` writes `filtered/` directly and does not create `filtered_pe/`.

UCE mode enables one rescue round by default after primary assembly has accepted contigs; `--no-uce-rescue-reads` disables it and `--uce-rescue-rounds 2` requests the second round. Rescue borrows the core baiting-and-iterative-mapping idea from MITObim and uses a fixed k=21. Round 1 recruits with original references plus accepted contigs; round 2 uses terminal windows only for loci still growing. The previous contig must remain intact, and each added side independently requires at least 30 bp, 85% breadth, a maximum 30-bp gap, two independent fragments, and one fragment bridging the frozen-core boundary; unsupported sides or loci revert individually. Before a round is finally accepted, its pre- and post-rescue sequences are also compared. A newly introduced exact inverted repeat of at least 150 bp restores that locus from the round backup with status `reverted_inverted_repeat`; `--uce-rescue-inverted-repeat-min-bp 0` disables this guard. The experimental `--uce-rescue-reverse-reuse-reference-scale` changes only the reference bonus for a node whose reverse complement is already present in either assembly arm; its default `1.0` disables scaling and never changes read depth. Rescue is constrained extension of a primary result, never reference gap filling.

## UCE selection, assembly, and QC

UCE mode uses the fused Rust `ucefilter` by default. Rolling-k-mer recruitment, run-k orientation verification, maximum-exact-match evidence, and adaptive per-locus selection are completed during one scan of the original reads, followed by direct output to `filtered/`. Weak evidence is first removed with the existing dynamic exact-match threshold. Reservoir reduction is bypassed when a locus has fewer than 512 candidates, estimated depth at most 160x, or exact seeds spanning fewer than 48 of 64 reference bins. A saturated locus retains at least 60% of its eligible core with reference-position diversity. Reads crossing bait or contig edges are stratified by overhang length so that a small set at each length preserves the overlap ladder from core to flank. It creates neither GM2 nor per-locus duplicated candidate FASTQ; each complete PE fragment is stored once and selected atomically. TipSeek follows the original GeneMiner2 resource-budget scheduler: `-p` is shared between ready filter, refilter, and assembly jobs; completed samples immediately enter the next stage queue. With `-p 1`, stages run serially. With larger values, assembly/refilter jobs receive a bounded 2--6-thread share and each recruitment job reserves 1--2 budget threads. The fragment-bank limit is internal and automatic: it uses the lower of Linux `MemAvailable` and the process cgroup's unused memory, reserves half for the system and other workflow buffers, divides the remainder across possible concurrent UCEFilter jobs, and caps each sample at 4 GiB. It falls back to 512 MiB when memory accounting is unavailable. The `backbone` strategy makes one bounded look-ahead decision at each bubble, commits the winning branch, and does not backtrack. See "UCE read routes" above for `--uce-rescue-reads`.

CPU parallelism defaults to `-p auto`: TipSeek counts physical cores allowed by affinity/cpuset, then applies cgroup or scheduler limits; an explicit integer overrides detection. The scheduler keeps the original 1--2-unit recruitment reservation separate from UCEFilter's supported single compute worker.

`--uce-alignment-shadow` is an opt-in evidence mode. After adaptive selection it deterministically samples at most 64 fragments per locus, locates a reference window with maximum-exact-match seeds, and runs local affine-gap alignment. Raw per-mate evidence and per-locus summaries are written to `alignment_shadow.tsv` and `alignment_shadow_summary.tsv`. Initial bait, whole-contig rescue, and terminal rescue evidence are preserved by round. Short-bait terminals are design boundaries; only assembly-contig terminals are suitable for extension assessment.

Inspect `uce_assembly_summary.csv`, `uce_rescue_summary.csv`, `uce_rescue_rounds.csv`, `contigs_all_low/`, and downstream alignments; per-side length, breadth, gap, fragment, and bridge evidence is retained in the round audit. Important guardrails include unique-read density, supported breadth, largest unsupported gap, k-mer depth CV, and max/median depth ratio. Use `--uce-path-strategy search` only for sensitivity comparison with legacy branch enumeration.

File definitions are in the [output guide](../manual/EN_US/output.md); options are in the [command-line guide](../manual/EN_US/command_line.md).
