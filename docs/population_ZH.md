# 4. Population

[English](population_EN.md) · [项目总览](../README.md) · 前置章节：[Filter](filter_ZH.md)、[Assembler](assembler_ZH.md) · [命令行指南](../manual/ZH_CN/command_line.md) · [输出说明](../manual/ZH_CN/output.md)

`population` 将多个已完成的二倍体 UCE 样本统一到一套参考坐标，输出联合 SNP、PCA 和可选 ADMIXTURE。它产生未定相基因型，不替代单倍型或单 locus 基因树分析。

## 输入与运行

每个样本需要保留：`uce_assembly_summary.csv`、`results/`中的接受 contig，以及样本表中的原始 R1/R2。样本表可在 R2 后附加可选的 `population` 和 `batch` 列。运行依赖 minibwa、samtools、bcftools 和 PLINK 1.9；ADMIXTURE 缺失时，其余结果仍会完成。

已完成 UCE 组装时，默认 `pseudoref` 直接复用已接受的 contig 与支持统计，不重新招募或建图：

```bash
cli/tipseek population \
  -f samples.tsv -r baits_by_locus -o output -p 32 \
  --assembly-mode uce --engine pseudoref
```

`baits_by_locus`中每个 FASTA 文件代表一个 locus。默认参考策略仍是 `sqcl-longest`；可用 `--population-reference-strategy supported` 单独比较支持优先策略。需要从 reads 构建局部图参考时再选择 `--engine panrefv2`。PanRefV2 不接受 `--population-reference-fasta`；若已有固定参考，使用 `--engine pseudoref --population-reference-fasta FILE`。

## 参考引擎

| 引擎 | 用途 |
| --- | --- |
| `pseudoref` | 默认；从已接受 contig 选择每个 locus 的代表序列。可用 `--population-reference-strategy sqcl-longest|supported`。 |
| `panref` | 旧实验实现，仅为兼容保留。 |
| `panrefv2` | 可选的流式局部图参考构建。 |

PanRefV2使用唯一 bait minimizer 将read pair归入locus：先以双端共同命中建立冻结核心图，再只接收连接核心图的单端救援read。每个输入FASTQ只解码一次；严格核心与实际接收的救援read在内存中仍受每样本/位点上限约束，每个样本使用一个顺序临时 candidate spool，使内存不随候选read数增长。Phred < 20 的碱基不进入图的 k-mer，重复的严格核心或实际接收的救援 pair 会折叠。救援筛选完成后自动删除。PanRefV2.2 还会从已接收的 ledger 写出稀疏的每样本 unitig 颜色证据；仅当一个样本支持 canonical backbone 的全部节点和转换边时，才保守地写出该样本路径；它不写按locus拆分的FASTQ，不跨UCE locus连接序列。

每条backbone按以下顺序选择：

```text
bait一致性 > 路径最小支持样本数 > 有观测证据的图边read/PE支持 > 路径最小样本深度稳定性 > 长度
```

无环局部图只使用已接收 ledger 实际支持的转换边进行全局路径解析；有环图安全回退到局部策略。`population_graph.gfa`为每个输出 backbone 写入方向正确的 `P` 路径。若支持图仍无法解析，会以 k=25 重建，并同时保留其图、路径、稀疏证据和 `locus_summary.tsv` 的 `assembly_k` 标识。

## 阶段与结果

```text
reference → mapping → calling → selection
```

- `reference`：构建或复制公共参考；
- `mapping`：minibwa 与 samtools 统一比对、去重复和QC；
- `calling`：bcftools 联合检测与过滤；
- `selection`：生成 all-SNP、每UCE一个SNP和LD-pruned面板，运行PCA与可选ADMIXTURE。

用 `--population-start-at` 和 `--population-stop-after` 重启或截停阶段；后续阶段只复用已检查的已有输出。

PanRefV2额外输出到 `population/reference/panrefv2/`：

- `index_metadata.tsv`：bait minimizer索引摘要；
- `recruitment_summary.tsv`：强证据、候选救援、candidate spool 字节数、实际救援和歧义pair；
- `population_graph.gfa`：按locus分区的局部DBG；
- `unitig_color_evidence.tsv`：每个 locus、unitig、样本的稀疏已接收 read 深度；
- `bubble_qc.tsv`：保守的图歧义 QC；只有隔离、无再次分叉且汇合到同一后续 backbone 节点的替代路径才标为 `simple_bubble`，其余分支保留为不作等位推断的 `complex_branch` 或 `terminal_branch`；
- `backbone_manifest.tsv`：每个已解析 backbone 的 SHA-256 稳定 ID；关联原始 locus、可选 FASTA record、GFA 路径、方向、k、长度、序列指纹和节点数；
- `backbone_coordinates.tsv`：每个 backbone node 在最终序列顺序中的 0-based 半开区间，以及其输出 GFA 方向；
- `sample_backbone_paths.tsv`：只有完整且有证据支持的 canonical 样本路径才写出 path；`partial` 与 `no_coverage` 行不会断言单倍型路径；
- `locus_summary.tsv`：状态、长度、unitig数、样本支持和backbone路径证据。

默认只有`pass` locus写入`population_reference.fasta`。`short`、`low_sample_support`、`low_coverage`、`complex`和`no_core`保留在统计中但不参与后续比对。仅作诊断时可加入：

```bash
--population-panrefv2-include-low-confidence
```

## 默认过滤与QC

默认阈值：MAPQ 20、base quality 20、genotype DP 5、GQ 20、site QUAL 20、call rate 0.8、MAC 2。低DP/GQ基因型先设为缺失，再执行位点过滤。

解释群体结果前至少检查：

- `reference/panrefv2/locus_summary.tsv`：恢复数、短locus和复杂图；
- `mapping/mapping_qc.tsv`：mapping rate、breadth和depth；
- `variants/variant_qc.tsv`：各过滤阶段SNP数；
- `structure/panel_summary.tsv`：三套面板的位点数；
- `structure/qc/individuals.imiss`、`individuals.het`、`individuals.genome`：个体缺失率、杂合度与成对 PI_HAT；样本很多时可用 `--skip-relatedness-qc` 跳过最后一项。
- 三套面板PCA是否给出一致的主要结构；以及ADMIXTURE不同K的CV误差。

主解释通常使用每UCE一个SNP的面板，并以all-SNP和LD-pruned结果作敏感性检查。

## 缓存校验与资源预算

`population/cache/` 保存成功完成的 reference、mapping、calling 和 filtering 阶段记录。恢复前检查样本身份与顺序、原始 reads 文件身份、UCE 汇总和结果内容、参考 SHA-256、相关参数以及外部工具文件指纹。原始 reads 和大型结果采用规范路径、文件大小、修改时间以及 Unix 文件系统身份/ctime 校验，避免为恢复再完整扫描 reads；这不是 reads 内容的密码学校验。参考和参考输入使用内容 SHA-256。

旧输出没有完成记录，或输入、参数、输出不匹配时，拒绝复用并提示重跑对应阶段。开始重写阶段前撤销旧记录，成功后才原子写入新记录。保持参考输入不变时：

- 修改 DP、GQ、QUAL、call rate 或 MAC：用 `--population-start-at calling`，复用已校验的原始联合检测结果，再做过滤。
- 修改 MAPQ 或 base quality：同样从 `calling` 开始，但重新联合检测。
- 仅修改 LD、PCA/ADMIXTURE 相关设置：用 `--population-start-at selection`，保留已验证的比对和过滤结果。
- 修改参考、样本集合或 reads：从 `reference` 重跑；不把旧 BAM 套用到新参考。

改变线程数不使生物学结果缓存失效。比对沿用样本队列，将同时运行的 mapper 与 samtools sort 合并计入线程预算；单线程时先写临时 SAM，再排序。sort 显式设置每线程内存上限，样本并发结合可用内存及 Linux cgroup-v2 限制收紧；这不是对外部程序全部内存的硬限制。实际安排写入 `mapping/resource_budget.tsv`。基因型检测仍使用原始 reads，不使用构图限额后的 reads 替代。

临时 BAM 不压缩，以减少反复压缩的 CPU 开销；至少有两个流水线线程预算时，fixmate 经坐标排序后直接通过管道送入去重，省去完整的中间坐标排序 BAM；内存不足时，排序仍可能使用各样本独立的临时文件。最终 BAM 仍压缩；单线程仍采用顺序临时文件。最终 BAM 写出时同步生成 `.bam.bai`，省去单独建索引时的一次完整 BAM 读取。samtools 的 `sort`/`markdup` 需要支持 `--write-index`（已用 1.16.1 验证）。实测范围和取舍见[真实 reads 验证](development/population-integration-validation.md)。
