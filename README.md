# TStools（原 GeneMiner2-UCE）

[![CI](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/ci.yml)
[![CodeQL](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/codeql.yml)
[![Fuzz smoke](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/fuzz-smoke.yml/badge.svg?branch=master)](https://github.com/GUIBA-EX/GeneMiner2-UCE/actions/workflows/fuzz-smoke.yml)
[![MSRV: 1.87](https://img.shields.io/badge/MSRV-1.87-orange)](rust-toolchain.toml)
[![Rust edition: 2021](https://img.shields.io/badge/Rust%20edition-2021-orange)](Cargo.toml)
[![Dependency policy: cargo-deny](https://img.shields.io/badge/dependency%20policy-cargo--deny-blue)](deny.toml)
[![SBOM: SPDX 2.3](https://img.shields.io/badge/SBOM-SPDX%202.3-blueviolet)](rust/xtask/src/main.rs)
[![Release integrity: SHA-256](https://img.shields.io/badge/release%20integrity-SHA--256-blueviolet)](rust/xtask/src/main.rs)
[![Latest release](https://img.shields.io/github/v/release/GUIBA-EX/GeneMiner2-UCE?display_name=tag)](https://github.com/GUIBA-EX/GeneMiner2-UCE/releases/latest)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

**[English](README_EN.md)** · [更新记录](CHANGELOG.md) · [报告问题](https://github.com/GUIBA-EX/GeneMiner2-UCE/issues)

TStools（原 GeneMiner2-UCE）是面向短 reads 的参考引导恢复工具：先以参考序列招募 reads，再按任务完成组装、证据量化或群体分析。它覆盖 genome skimming、target capture、UCE、线粒体、核基因家族、RAD 补充和无参考 repeatome；生产工作流为 Rust 原生实现，不需要 Python 运行时。

> **与原版 GeneMiner2 的关系：** GeneMiner2 是本项目的算法来源与结果兼容基线；TStools 不是对上游的维护型 fork，而是独立演进的工作流工具箱。`cli/geneminer2`、部分历史输出名、仓库地址和论文中的旧名称仅为兼容标识，并不表示两者功能或实现完全相同。

![TStools 流程](docs/images/summary_ZH.png)

## 与原版 GeneMiner2 的主要区别

| 维度 | 原版 GeneMiner2 | 当前 TStools（v1.5.8） |
| --- | --- | --- |
| 项目定位 | 基因组浅层测序数据的基因恢复算法与原始工作流 | 在兼容基线之上独立演进的短 reads 恢复与分析工具箱 |
| 核心实现 | 上游实现 | Rust 原生生产路径；不依赖 Python 运行时 |
| reads 招募 | 原始招募语义 | canonical 双链 2-bit k-mer、内容校验的参考缓存及有界流式 I/O；保持招募语义与历史输出兼容，同时降低 CPU、内存和 I/O 开销 |
| 常规组装 | 上游算法基线 | 默认确定性的 `original-rust`；保留上游 `original` 路线用于严格对照与复现 |
| UCE 组装 | 非专用的一般组装路径 | `ucefilter → uce-rust` 完成招募、成对 fragments 保留、方向/精确匹配证据和逐 locus 选择；默认一轮 rescue 仅接受 reads 支持的延伸，绝不以参考补洞 |
| 工作流范围 | 常规基因恢复 | 另含线粒体、marker profiling、UCE 群体、核基因家族、RAD 矩阵补充和无参考 repeatome |
| 结果解释 | 以恢复序列为主 | 闭环、RAD 严格矩阵和群体图路径都要求显式证据；输出 QC、provenance 与审计记录 |

因此，若目标是复现上游基线，可显式选择 `original` 并保留完整输入与参数；若目标是 UCE、群体或扩展工作流，应按 TStools 的专用模式和 QC 规则解释结果，不能将其直接等同于原版 GeneMiner2 的输出。

算法、性能与适用边界见 [Filter](docs/filter_ZH.md)、[Assembler](docs/assembler_ZH.md) 和 [Population](docs/population_ZH.md)。

## 安装与最小示例

安装 Rust/Cargo 与所需生物信息学软件后，在仓库根目录运行：

```bash
cargo run -p xtask -- build
```

统一入口是 `cli/geneminer2`。样本表使用 tab 分隔的 `sample_id  R1  [R2]`；参考目录中每个 FASTA 文件代表一个 locus 或 bait。

```bash
# UCE：原始 paired reads → 选择性招募 → UCE 组装
cli/geneminer2 filter assemble \
  -f samples.tsv -r uce_references -o uce_out -p 8 \
  --assembly-mode uce
```

首先查看 `uce_out/uce_assembly_summary.csv` 和 `uce_out/uce_contigs/`。UCE 默认使用 `kf=23`、`step=4`、`auto` 敏感招募和一轮受证据约束的 rescue；rescue 不会用参考序列虚构缺口。可用 `--no-uce-rescue-reads` 关闭 rescue，或用 `--uce-rescue-rounds 2` 显式运行第二轮。

如需重现旧的单遍 k=31、无 rescue 路径，请显式覆盖默认值：

```bash
cli/geneminer2 filter assemble \
  -f samples.tsv -r uce_references -o uce_out -p auto \
  --assembly-mode uce -kf 31 --uce-recruit-mode fast --no-uce-rescue-reads
```

`auto` 先以默认 k=23、step=4 执行快速招募，只对快速阶段没有选中 fragments 的 locus 再扫描一次 FASTQ（fallback 默认 k=21、step=1、独立验证 k=19）。敏感阶段先用未恢复 locus 子集作粗招募门控；只对命中该子集的 fragments 用完整 probe 面板扩展候选并检查多 locus 歧义，再要求 read pair 至少一端与目标 probe/reference 的局部比对达到 45 bp、80% identity，最后只合并唯一支持某个未恢复 locus 的 reads；随后只组装一次。仅由敏感阶段恢复的 provisional core 还必须达到 200 bp，并与目标 probe 达到至少 80% coverage、80% identity，不能存在得分接近的其他 probe locus，也不能包含至少 150 bp 的长倒置重复，才会锚定为该 locus。锚定 core 中至少 40 bp 的内部 reads 链缺口会被明确标记为复核项，但不会单独硬拒绝；同个体基因组校准显示，低覆盖的真实共线位点也可能出现该信号。默认的一轮 rescue 只允许无复核项的锚定 core 作为样本自身 bait 进入 whole-contig rescue 与 guard；review-only core 保留为候选，但不会被 rescue 扩展。该策略扩大的是候选 reads 招募范围，不代表候选位点已经被证明正确。

自动敏感招募会保留 `uce_filter_summary.fast.tsv`；实际发生第二阶段时还会生成 `uce_filter_summary.fallback.tsv`。`uce_recruit_passes.tsv` 逐 locus 记录快速阶段、敏感阶段和最终来源，`uce_recruit_contig_probe_gate.tsv` 覆盖所有 fallback locus，记录 assembler、probe、锚定、倒置重复和内部 gap 证据；被拒绝的结果移入 `fallback_probe_rejected/`，不会静默删除。启用 rescue 时，后续状态继续写入 `uce_rescue_rounds.csv` 和 `uce_rescue_summary.csv`。

## 选择命令

| 目标 | 命令 | 主要结果 |
| --- | --- | --- |
| 常规 exon、SCO 或核 marker | `filter assemble` | 参考引导 contig |
| UCE core 与 reads 支持的侧翼 | `filter assemble --assembly-mode uce` | UCE contig 与 QC 汇总 |
| 线粒体 | `mito` | closed、linear 或明确的结构歧义判定 |
| marker 的 reads 支持 | `profiling` | 每条参考的支持证据 |
| UCE 群体遗传 | `population` | 公共伪参考、VCF、PCA 等 |
| 核基因家族 | `gene` | family 候选、copy 状态与后续解析输入 |
| WGS 补充既有 RAD 矩阵 | `rad-probe` → `rad` → `rad-validate` | 双 arm 恢复与严格矩阵 |
| 无参考 repeatome | `te` | 保守 repeat library、注释和 RPM |

默认的常规组装使用确定性的 `original-rust`；UCE 使用专用的 `ucefilter` 与 `uce-rust`。两条路线都是 Rust 实现。

## 常用命令

```bash
# 常规 marker；original 为默认模式
cli/geneminer2 filter assemble -f samples.tsv -r references -o marker_out -p 8

# 单环动物线粒体；需要带注释的 GenBank 参考
cli/geneminer2 mito -f samples.tsv -o mito_out -p 8 \
  --mito-genbank mitochondrial_reference.gb

# 核基因家族
cli/geneminer2 gene -f samples.tsv -r family_baits -o gene_out -p 8

# 不组装，只评估 reads 对 marker 的支持
cli/geneminer2 profiling -f samples.tsv -r marker_reference.fasta -o profile_out -p 8

# UCE 组装完成后的群体流程
cli/geneminer2 population -f samples.tsv -r uce_references -o population_out -p 8 \
  --assembly-mode uce --engine panrefv2

# 无参考 repeatome；使用独立样本表
cli/geneminer2 te -f te_samples.tsv -o te_out -p 32
```

## 结果的边界

- UCE 默认将广泛招募与每 locus 的 reads 选择合并为一次 FASTQ 扫描；`--legacy-uce-filter` 仅用于对照与诊断。
- `mito` 只适用于常规单环动物线粒体。串联重复或超过 insert size 的完全重复不能由短 reads 可靠定拷贝数，结果会保留为 linear 或 ambiguous。
- `profiling` 是参考相容性证据，不是物种鉴定或丰度估计。
- RAD 中 R1/R2 是独立限制性位点 arm；WGS 恢复不直接证明 allele dropout。请以 `rad-validate` 的双 arm 检查为准。
- `--cleanup-intermediates` 只在同次完整流程成功后删除可再生的过滤 reads；先加 `--cleanup-dry-run` 可生成 `cleanup_preview.tsv` 审核候选路径和字节数，且不会删除任何文件。最终 contig、汇总、原始 reads 和参考始终保留。

## 文档

| 主题 | 文档 |
| --- | --- |
| 安装、参数、外部依赖 | [命令行指南](manual/ZH_CN/command_line.md) |
| 输出目录与结果表 | [输出说明](manual/ZH_CN/output.md) |
| 过滤与缓存 | [Filter](docs/filter_ZH.md) |
| 常规/UCE 组装 | [Assembler](docs/assembler_ZH.md) |
| 线粒体、Gene、RAD、TE | [Mito](docs/mitochondria_CN.md) · [Gene](docs/gene_ZH.md) · [RAD](docs/rad_CN.md) · [TE](docs/te_ZH.md) |
| 群体与 profiling | [Population](docs/population_ZH.md) · [Profiling](docs/profiling_ZH.md) |

`--workflow-profile` 会写入仅记录时间与 I/O 的 `workflow_profile.tsv`，不会改变分析结果。

每次标准工作流都会在输出根目录原子写入 `workflow_manifest.tsv`，记录 CLI 版本、命令、参考 SHA-256、样本表 SHA-256、关键参数以及输入 reads 的路径/大小/修改时间，用于复现和审计。

`--resume` 是保守的整次工作流恢复：仅当现有 manifest 与本次输入/参数完全一致且 `workflow_status.tsv` 为成功时才作为无操作成功返回；任何失败或不匹配都会拒绝，且不会覆盖旧状态或跳过部分 stage。

输出目录一旦创建，CLI 会在结束时原子写入 `workflow_status.tsv`；其中 `state` 为 `succeeded` 或 `failed`，失败时附带错误摘要，便于批处理系统识别部分输出不可消费。

`cargo run -p xtask -- build` 会在 `cli/` 中生成 `SHA256SUMS` 与 `SBOM.spdx.json`，用于发布校验与二进制物料清单。FASTX fuzz 目标位于 `fuzz/`，仅由手动或每周的受限 CI smoke 任务运行（1,000 次输入、最多 60 秒）。

## 引用

请引用：Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111. [doi:10.1111/1755-0998.70111](https://doi.org/10.1111/1755-0998.70111)

```bibtex
@software{TStools,
  author  = {XIA, Fei and TANG, Zizhen and XU, Yan},
  title   = {TStools (formerly GeneMiner2-UCE): Reference-Guided Short-Read Recovery for UCE, Mitochondrial, Gene-Family, and RAD Workflows},
  year    = {2026},
  version = {1.5.8},
  url     = {https://github.com/GUIBA-EX/GeneMiner2-UCE},
  publisher = {GitHub},
  note    = {GPL-3.0-or-later licensed software}
}
```

项目以 [GPL-3.0-or-later](LICENSE) 发布；第三方与移植代码的来源边界见 [NOTICE](NOTICE)。
