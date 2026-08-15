# TipSeek

[![CI](https://github.com/GUIBA-EX/TipSeek/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/TipSeek/actions/workflows/ci.yml)
[![CodeQL](https://github.com/GUIBA-EX/TipSeek/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/TipSeek/actions/workflows/codeql.yml)
[![Latest release](https://img.shields.io/github/v/release/GUIBA-EX/TipSeek?display_name=tag)](https://github.com/GUIBA-EX/TipSeek/releases/latest)
[![MSRV: 1.87](https://img.shields.io/badge/MSRV-1.87-orange)](rust-toolchain.toml)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

**[English](README_EN.md)** · [更新记录](CHANGELOG.md) · [报告问题](https://github.com/GUIBA-EX/TipSeek/issues)

TipSeek 是 Rust 原生的短 reads 恢复与分析工具。它以参考序列招募 reads，再按任务完成组装、证据量化或群体分析，覆盖 genome skimming、target capture、UCE、线粒体、核基因家族、RAD 补充和无参考 repeatome。统一入口为 `tipseek`，运行时不依赖 Python。

![TipSeek 工作流](docs/images/summary_ZH.png)

## 工作流

| 目标 | 命令 | 主要输出 |
| --- | --- | --- |
| Exon、SCO 或其他 marker 恢复 | `filter assemble` | 参考引导 contigs |
| UCE core 与 reads 支持的侧翼 | `filter assemble --assembly-mode uce` | UCE contigs、恢复汇总和逐 locus 证据 |
| 动物线粒体恢复 | `mito` | closed、linear 或 ambiguous 结构判定 |
| Marker 支持度评估 | `profiling` | 每条参考的 reads 支持 |
| UCE 群体分析 | `population` | cohort reference、VCF、PCA 等 |
| 核基因家族恢复 | `gene` | family 候选、copy 状态和解析输入 |
| WGS 补充 RAD 矩阵 | `rad-probe` → `rad` → `rad-validate` | 双 arm 恢复和严格矩阵 |
| 无参考 repeatome | `te` | repeat library、注释和 RPM |

各工作流共享输入、并行调度和运行状态记录，但保留各自的证据模型。UCE 路径采用 fragment-aware 分层招募、核心与末端的分离预算、PE-supported 双图组装和逐 locus 可逆 rescue；其他任务只运行与其推断目标相关的步骤。

## 安装

完整依赖见[命令行指南](manual/ZH_CN/command_line.md)。从源码构建：

```bash
git clone https://github.com/GUIBA-EX/TipSeek.git
cd TipSeek
cargo run -p xtask -- build
cli/tipseek -h
```

构建结果位于 `cli/`，并包含 `SHA256SUMS` 和 `SBOM.spdx.json`。

## UCE 最小示例

样本表为 tab 分隔文本，每行为 `sample_id  R1  [R2]`；参考目录中的每个 FASTA 文件代表一个 locus 或 bait。

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

UCE 模式默认使用 k=23、step=4、`auto` 招募和一轮 evidence-constrained rescue。可用 `--no-uce-rescue-reads` 关闭 rescue，或用 `--uce-rescue-rounds 2` 请求第二轮；完整参数和旧路径复现方式见[命令行指南](manual/ZH_CN/command_line.md#73-组装与-uce)。

首先检查：

- `uce_assembly_summary.csv`：每个样本和 locus 的恢复状态；
- `uce_contigs/`：最终接受的候选序列；
- `uce_recruit_passes.tsv` 和 `uce_recruit_contig_probe_gate.tsv`：招募来源、probe 门控和候选状态；
- `uce_rescue_rounds.csv` 和 `uce_rescue_summary.csv`：每轮接纳、裁切或回滚结果。

## 证据与结果边界

- TipSeek 的组装与 rescue 以 reads 证据为准；UCE rescue 不使用参考序列填补缺口。候选可标记为 accepted、review 或 rejected，review-only core 不参与 rescue。
- `original` 组装模式用于常规 marker，也保留 GeneMiner2 基线复现路径；TipSeek 的 UCE、群体和其他工作流应按各自的 QC 表解释。
- `mito` 面向常规单环动物线粒体。超过 insert size 的完全重复不能由短 reads 可靠确定拷贝数，结果会保留为 linear 或 ambiguous。
- `profiling` 报告 reads 与参考序列的相容性，不等同于物种鉴定或丰度估计。
- RAD 的 R1/R2 是独立限制性位点 arms；WGS 恢复本身不能证明 allele dropout，应以 `rad-validate` 的双 arm 检查为准。

## 可复现运行

- `workflow_manifest.tsv` 记录命令、版本、关键参数、参考与样本表 SHA-256，以及输入 reads 元数据。
- `workflow_status.tsv` 原子记录 `succeeded` 或 `failed`；`--resume` 只在输入、参数和成功状态完全一致时返回已有结果。
- `--workflow-profile` 仅记录时间与 I/O，不改变分析；`--cleanup-dry-run` 可在删除可再生中间文件前生成审核清单。

## 文档

| 内容 | 中文 | English |
| --- | --- | --- |
| 安装、输入与参数 | [命令行指南](manual/ZH_CN/command_line.md) | [Command-line guide](manual/EN_US/command_line.md) |
| 输出目录与结果表 | [输出说明](manual/ZH_CN/output.md) | [Output reference](manual/EN_US/output.md) |
| Filter 与缓存 | [Filter](docs/filter_ZH.md) | [Filter](docs/filter_EN.md) |
| 常规与 UCE 组装 | [Assembler](docs/assembler_ZH.md) | [Assembler](docs/assembler_EN.md) |
| 线粒体 | [Mito](docs/mitochondria_CN.md) | [Mito](docs/mitochondria_EN.md) |
| Gene、RAD、TE | [Gene](docs/gene_ZH.md) · [RAD](docs/rad_CN.md) · [TE](docs/te_ZH.md) | [Gene](docs/gene_EN.md) · [RAD](docs/rad_EN.md) · [TE](docs/te_EN.md) |
| Population 与 profiling | [Population](docs/population_ZH.md) · [Profiling](docs/profiling_ZH.md) | [Population](docs/population_EN.md) · [Profiling](docs/profiling_EN.md) |

## 引用与许可

引用当前软件版本：

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

方法来源或 `original` 基线相关分析另请引用：Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26:e70111. [doi:10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

TipSeek 以 [GPL-3.0-or-later](LICENSE) 发布；第三方与移植代码的来源见 [NOTICE](NOTICE)。
