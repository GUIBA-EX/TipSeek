# 4. Population

[中文版本](population_ZH.md)

`population` creates one cohort coordinate system for multiple diploid UCE samples, then produces a joint SNP set, PCA, and optional ADMIXTURE results. It reports unphased genotypes; it is not a phasing or per-locus gene-tree workflow.

## Input and quick start

Each sample needs its completed UCE assembly (`uce_assembly_summary.csv` and accepted contigs in `results/`) plus the original R1/R2 files listed in the sample table. For `population`, optional `population` and `batch` fields may follow R2; leave the R2 field empty for single-end samples with metadata. The workflow requires minibwa, samtools, bcftools, and PLINK 1.9; ADMIXTURE is optional.

```bash
cli/tipseek population \
  -f samples.tsv -r baits_by_locus -o output -p 32 \
  --assembly-mode uce --engine panrefv2
```

`baits_by_locus` contains one FASTA per locus. Use `pseudoref` for accepted-contig representative references, or `--population-reference-fasta FILE` with `pseudoref` for a fixed external reference.

## Reference engines

| Engine | Purpose |
| --- | --- |
| `pseudoref` | Default accepted-contig representative reference; supports `sqcl-longest` and `supported`. |
| `panref` | Legacy experimental implementation. |
| `panrefv2` | Recommended streaming local-graph UCE reference builder. |

PanRefV2 builds a frozen paired-read core per locus, then admits only bounded single-mate rescue reads connected to that core. Each input FASTQ is decoded once; strict-core and accepted-rescue evidence remain bounded in memory, while one sequential temporary candidate spool per sample bounds RAM independently of the number of candidate reads. Low-quality bases (Phred < 20) are excluded from graph k-mers and duplicate strict-core or accepted-rescue pairs are collapsed. The temporary spool is removed after rescue selection. PanRefV2.2 also writes sparse per-sample unitig color evidence and conservatively reports a sample backbone path only when that sample supports every canonical node and transition; it writes no per-locus recruited FASTQ and never scaffolds across UCE loci.

Backbone order: bait agreement > minimum supporting samples across the path > observed graph-edge read/PE support > minimum per-sample depth stability across the path > length. Only transitions with accepted-ledger evidence participate in global path extension. Graph edges are resolved globally when the local graph is acyclic; cyclic graphs use the safe local fallback. `population_graph.gfa` records every emitted backbone as a direction-correct `P` path. An unresolved supported graph is retried at k=25; its graph, path, sparse evidence, and `assembly_k` record remain reproducible alongside the standard k=31 graph.

## Stages and outputs

```text
reference → mapping → calling → selection
```

The stages build the reference, map/QC reads, jointly call/filter variants, then create all-SNP, one-SNP-per-UCE, and LD-pruned panels with PCA and optional ADMIXTURE. Use `--population-start-at` and `--population-stop-after` only with validated outputs.

PanRefV2 writes `population/reference/panrefv2/`:

- `index_metadata.tsv` — bait minimizer summary;
- `recruitment_summary.tsv` — strong, candidate rescue, candidate-spool bytes, accepted rescue, and ambiguous pairs;
- `population_graph.gfa` — locus-partitioned local graphs;
- `unitig_color_evidence.tsv` — sparse accepted read depth for each locus, unitig, and sample;
- `bubble_qc.tsv` — conservative graph-ambiguity QC: only isolated, non-branching alternatives that rejoin one later backbone node are `simple_bubble`; all other branches remain non-allelic `complex_branch` or `terminal_branch`;
- `backbone_manifest.tsv` — SHA-256 stable IDs for every resolved backbone, mapping the original locus, optional FASTA record, GFA path, orientation, k, length, sequence fingerprint, and node count;
- `backbone_coordinates.tsv` — each backbone node as a 0-based half-open interval in final sequence order, with its emitted GFA orientation;
- `sample_backbone_paths.tsv` — a path string only for complete, evidence-supported canonical sample paths; `partial` and `no_coverage` rows never assert a haplotype path;
- `locus_summary.tsv` — locus status and backbone evidence.

Only `pass` loci enter `population_reference.fasta` by default. `short`, `low_sample_support`, `low_coverage`, `complex`, and `no_core` remain in the QC report. Use `--population-panrefv2-include-low-confidence` only for diagnostics.

## QC

Default filters are MAPQ 20, base quality 20, genotype DP 5, GQ 20, site QUAL 20, call rate 0.8, and MAC 2. Check PanRefV2 locus status, mapping QC, variant QC, panel counts, PCA agreement across the three panels, and ADMIXTURE CV error before interpretation. Individual QC is written to `structure/qc/individuals.imiss`, `individuals.het`, and `individuals.genome` (pairwise PI_HAT); use `--skip-relatedness-qc` to omit the quadratic relatedness calculation for large cohorts.
