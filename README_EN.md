# TipSeek

[![CI](https://github.com/GUIBA-EX/TipSeek/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/TipSeek/actions/workflows/ci.yml)
[![CodeQL](https://github.com/GUIBA-EX/TipSeek/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/TipSeek/actions/workflows/codeql.yml)
[![Latest release](https://img.shields.io/github/v/release/GUIBA-EX/TipSeek?display_name=tag)](https://github.com/GUIBA-EX/TipSeek/releases/latest)
[![MSRV: 1.87](https://img.shields.io/badge/MSRV-1.87-orange)](rust-toolchain.toml)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

**[中文](README.md)** · [Changelog](CHANGELOG.md) · [Report an issue](https://github.com/GUIBA-EX/TipSeek/issues)

TipSeek is a Rust-native toolkit for recovering and analysing short reads. It recruits reads against references, then performs task-specific assembly, evidence quantification, or cohort analysis for genome skimming, target capture, UCEs, mitochondria, nuclear gene families, RAD augmentation, and reference-free repeatomes. The single entry point is `tipseek`; no Python runtime is required.

![TipSeek workflow](docs/images/summary_EN.png)

## Workflows

| Goal | Command | Main output |
| --- | --- | --- |
| Recover exons, SCOs, or other markers | `filter assemble` | Reference-guided contigs |
| Recover UCE cores and read-supported flanks | `filter assemble --assembly-mode uce` | UCE contigs, recovery summary, and per-locus evidence |
| Recover animal mitochondria | `mito` | Closed, linear, or ambiguous structure call |
| Measure marker support | `profiling` | Read support for each reference |
| Analyse UCE population data | `population` | Cohort reference, VCF, PCA, and related outputs |
| Recover nuclear gene families | `gene` | Family candidates, copy states, and resolution inputs |
| Add WGS samples to a RAD matrix | `rad-probe` → `rad` → `rad-validate` | Independent-arm recovery and strict matrix |
| Analyse a reference-free repeatome | `te` | Repeat library, annotation, and RPM |

The workflows share input handling, parallel scheduling, and run-state records while retaining task-specific evidence models. The UCE path uses fragment-aware hierarchical recruitment, separate core and terminal evidence budgets, PE-supported dual-graph assembly, and reversible per-locus rescue. Other workflows run only the steps required for their inference target.

## Installation

See the [command-line guide](manual/EN_US/command_line.md) for complete dependencies. Build from source:

```bash
git clone https://github.com/GUIBA-EX/TipSeek.git
cd TipSeek
cargo run -p xtask -- build
cli/tipseek -h
```

Build artifacts are written to `cli/` together with `SHA256SUMS` and `SBOM.spdx.json`.

## Minimal UCE example

The sample manifest is tab-delimited, with one `sample_id  R1  [R2]` record per line. Each FASTA file in the reference directory represents one locus or bait.

```text
sample_1<TAB>/data/sample_1_R1.fastq.gz<TAB>/data/sample_1_R2.fastq.gz
sample_2<TAB>/data/sample_2_R1.fastq.gz<TAB>/data/sample_2_R2.fastq.gz
```

```bash
cli/tipseek filter assemble \
  -f samples.tsv \
  -r uce_references \
  -o uce_out \
  -p auto \
  --assembly-mode uce
```

UCE mode defaults to k=23, step=4, `auto` recruitment, and one evidence-constrained rescue round. Use `--no-uce-rescue-reads` to disable rescue or `--uce-rescue-rounds 2` to request a second round. See the [command-line guide](manual/EN_US/command_line.md#73-assembly-and-uce-options) for all options and legacy-path reproduction.

Inspect these outputs first:

- `uce_assembly_summary.csv`: recovery status for each sample and locus;
- `uce_contigs/`: final accepted candidate sequences;
- `uce_recruit_passes.tsv` and `uce_recruit_contig_probe_gate.tsv`: recruitment source, probe gates, and candidate state;
- `uce_rescue_rounds.csv` and `uce_rescue_summary.csv`: acceptance, trimming, or rollback in each rescue round.

## Evidence and interpretation boundaries

- TipSeek assembly and rescue are governed by read evidence; UCE rescue never fills gaps from the reference. Candidates may be accepted, flagged for review, or rejected, and review-only cores cannot seed rescue.
- The `original` assembly mode handles conventional markers and retains the GeneMiner2 baseline path. Interpret TipSeek UCE, population, and other workflows through their task-specific QC tables.
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
| Conventional and UCE assembly | [Assembler](docs/assembler_EN.md) | [Assembler](docs/assembler_ZH.md) |
| Mitochondria | [Mito](docs/mitochondria_EN.md) | [Mito](docs/mitochondria_CN.md) |
| Gene, RAD, and TE | [Gene](docs/gene_EN.md) · [RAD](docs/rad_EN.md) · [TE](docs/te_EN.md) | [Gene](docs/gene_ZH.md) · [RAD](docs/rad_CN.md) · [TE](docs/te_ZH.md) |
| Population and profiling | [Population](docs/population_EN.md) · [Profiling](docs/profiling_EN.md) | [Population](docs/population_ZH.md) · [Profiling](docs/profiling_ZH.md) |

## Citation and licence

Cite the current software release as:

```bibtex
@software{TipSeek,
  author    = {XIA, Fei and TANG, Zizhen and XU, Yan},
  title     = {TipSeek: Reference-Guided Short-Read Recovery and Analysis},
  year      = {2026},
  version   = {1.6.2},
  url       = {https://github.com/GUIBA-EX/TipSeek},
  publisher = {GitHub}
}
```

For method provenance or analyses using the `original` baseline, also cite: Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26:e70111. [doi:10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

TipSeek is released under [GPL-3.0-or-later](LICENSE). See [NOTICE](NOTICE) for the provenance of third-party and ported code.
