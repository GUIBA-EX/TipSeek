<p align="center">
  <img src="docs/images/tipseek_logo.png" alt="TipSeek logo" width="760">
</p>

# TipSeek

[![CI](https://github.com/GUIBA-EX/TipSeek/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/TipSeek/actions/workflows/ci.yml)
[![CodeQL](https://github.com/GUIBA-EX/TipSeek/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/TipSeek/actions/workflows/codeql.yml)
[![Latest release](https://img.shields.io/github/v/release/GUIBA-EX/TipSeek?display_name=tag)](https://github.com/GUIBA-EX/TipSeek/releases/latest)
[![MSRV: 1.87](https://img.shields.io/badge/MSRV-1.87-orange)](rust-toolchain.toml)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

**[中文](README.md)** · [Changelog](CHANGELOG.md) · [Report an issue](https://github.com/GUIBA-EX/TipSeek/issues)

TipSeek is a Rust-native toolkit for short-read recovery and analysis. Through one entry point, it performs reference-guided read recruitment, target assembly, evidence reporting, and cohort analysis for genome skimming, target capture, nuclear gene families, UCEs, animal mitochondria, RAD augmentation, and reference-free repeatomes. Release builds require no Python runtime.

## Workflows

In the table below, `tipseek` refers to the built `cli/tipseek` executable.

| Goal | Command | Main output |
| --- | --- | --- |
| Recover nuclear gene-family candidates | `tipseek` (default `gene` mode) | Per-sample candidates and cross-sample family summaries |
| Annotate exon/intron structure | `tipseek --assembly-mode exon` | Gene candidates, CDS, exons, introns, and supercontigs |
| Recover UCE cores and read-supported flanks | `tipseek --assembly-mode uce` | UCE contigs, recovery summaries, and per-locus evidence |
| Recover animal mitochondria | `tipseek mito` | Closed, linear, or ambiguous structure call |
| Measure marker support | `tipseek profiling` | Read support for each reference sequence |
| Analyse UCE population data | `tipseek population` | Cohort reference, VCF, PCA, and related outputs |
| Add WGS samples to a RAD matrix | `tipseek rad-probe` → `tipseek rad` → `tipseek rad-validate` | Independent-arm recovery and strict matrix |
| Analyse a reference-free repeatome | `tipseek te` | Repeat library, annotation, and RPM |

All workflows share input validation, CPU scheduling, and run-state records while applying task-specific biological criteria. The UCE path uses fragment-aware hierarchical recruitment, independent core and terminal evidence budgets, PE-supported dual-graph assembly, and reversible per-locus rescue. Other workflows run only the steps needed for their inference target.

## Installation

See the [command-line guide](manual/EN_US/command_line.md) for complete dependencies and platform notes. Build from source:

```bash
git clone https://github.com/GUIBA-EX/TipSeek.git
cd TipSeek
cargo run -p xtask -- build
cli/tipseek -h
```

Build artifacts are written to `cli/` together with `SHA256SUMS` and `SBOM.spdx.json`.

## Quick start

The sample manifest is tab-delimited, with one `sample_id  R1  [R2]` record per line. Absolute paths are recommended.

```text
sample_1<TAB>/data/sample_1_R1.fastq.gz<TAB>/data/sample_1_R2.fastq.gz
sample_2<TAB>/data/sample_2_R1.fastq.gz<TAB>/data/sample_2_R2.fastq.gz
```

Each FASTA file in the reference directory defines one gene family or UCE locus.

```bash
# Default: recover nuclear gene-family candidates
cli/tipseek -f samples.tsv -r family_references -o gene_out -p auto

# Derive protein references from nucleotide baits, then annotate exon/intron structure
cli/tipseek --assembly-mode exon \
  -f samples.tsv -r family_references \
  -o exon_out -p auto
```

To recover UCE contigs only:

```bash
cli/tipseek filter assemble \
  -f samples.tsv \
  -r uce_references \
  -o uce_out \
  -p auto \
  --assembly-mode uce
```

UCE mode defaults to k=23, step=4, `auto` recruitment, and one evidence-constrained rescue round. Use `--no-uce-rescue-reads` to disable rescue or `--uce-rescue-rounds 2` to request a second round. See the [command-line guide](manual/EN_US/command_line.md#73-assembly-and-uce-options) for all options.

Inspect these UCE outputs first:

- `<sample>/uce_assembly_summary.csv`: per-locus recovery status and support metrics;
- `<sample>/results/`: final accepted contigs;
- `<sample>/uce_recruit_passes.tsv` and `<sample>/uce_recruit_contig_probe_gate.tsv`: recruitment source, probe gates, and candidate state;
- `<sample>/uce_rescue_rounds.csv` and `<sample>/uce_rescue_summary.csv`: acceptance, trimming, or rollback in each rescue round.

## Evidence and interpretation boundaries

- TipSeek assembly and rescue are governed by read evidence. UCE rescue never fills gaps from the reference, and review-only cores cannot seed rescue.
- Candidate counts in `gene` and `exon` are assembly evidence, not allele counts or biological copy-number estimates.
- `mito` targets ordinary single-circular animal mitochondria. Short reads cannot reliably determine the copy number of exact repeats longer than the insert size, so unresolved cases remain linear or ambiguous.
- `profiling` reports compatibility between reads and references; it is not species identification or abundance estimation.
- RAD R1/R2 are independent restriction-site arms. WGS recovery alone does not demonstrate allele dropout; use the two-arm checks from `rad-validate`.

## Reproducible runs

- `workflow_manifest.tsv` records the command, version, key options, reference and manifest SHA-256 values, and input-read metadata.
- `workflow_status.tsv` atomically records `succeeded` or `failed`. `--resume` reuses output only when the inputs, options, and successful state match exactly.
- `--workflow-profile` records timing and I/O without changing the analysis. `--cleanup-dry-run` previews removable intermediate files before cleanup.

## Documentation

| Topic | English | 中文 |
| --- | --- | --- |
| Installation, inputs, and options | [Command-line guide](manual/EN_US/command_line.md) | [命令行指南](manual/ZH_CN/command_line.md) |
| Output directories and tables | [Output reference](manual/EN_US/output.md) | [输出说明](manual/ZH_CN/output.md) |
| Filtering and caches | [Filter](docs/filter_EN.md) | [Filter](docs/filter_ZH.md) |
| Gene, exon, and UCE assembly | [Assembler](docs/assembler_EN.md) | [Assembler](docs/assembler_ZH.md) |
| Mitochondria | [Mito](docs/mitochondria_EN.md) | [Mito](docs/mitochondria_CN.md) |
| Gene, exon, RAD, and TE | [Gene](docs/gene_EN.md) · [Exon](docs/exon_EN.md) · [RAD](docs/rad_EN.md) · [TE](docs/te_EN.md) | [Gene](docs/gene_ZH.md) · [Exon](docs/exon_ZH.md) · [RAD](docs/rad_CN.md) · [TE](docs/te_ZH.md) |
| Population and profiling | [Population](docs/population_EN.md) · [Profiling](docs/profiling_EN.md) | [Population](docs/population_ZH.md) · [Profiling](docs/profiling_ZH.md) |

## Citation and licence

Cite the current software release as:

```bibtex
@software{TipSeek,
  author    = {XIA, Fei and TANG, Zizhen and XU, Yan},
  title     = {TipSeek: Reference-Guided Short-Read Recovery and Analysis},
  year      = {2026},
  version   = {1.6.3},
  url       = {https://github.com/GUIBA-EX/TipSeek},
  publisher = {GitHub}
}
```

When using the GeneMiner2-derived gene assembler, also cite: Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26:e70111. [doi:10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

TipSeek is released under [GPL-3.0-or-later](LICENSE). See [NOTICE](NOTICE) for the provenance of third-party and ported code.
