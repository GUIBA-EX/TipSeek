# 输出文件

本文档说明当前命令行版本会生成的输出文件。旧 GUI、细胞器组装和 Windows 辅助脚本相关输出已经不属于本 CLI fork 的支持范围。

## 单样本目录

输入样本表中的每个样本都会在输出目录下生成一个同名文件夹。

**uce_filter_summary.tsv**：逐 locus 自动选择审计表。`selection_mode` 为 `pass-through`、`core`、`rescue` 或 `legacy-fallback`；同时报告合格/已选 fragment、64-bin breadth、核心目标以及左右端候选与保留数。

**filtered**：默认融合 UCEFilter 直接写出的每-locus interleaved paired-end FASTQ。UCE 模式不再生成 `ucefilter_candidates`；**filtered_pe** 仅保留给 legacy MainFilter/GM2 回退路径。

**filtered**：进一步过滤后保留的 reads。UCE 模式下，只要 paired-end 的任一端通过 locus 过滤，整对 reads 都会被保留。

**alignment_shadow.tsv**：仅在 `--uce-alignment-shadow` 下生成的逐 mate 内部比对证据，包括 identity、overlap、linked-mate、terminal 和参考坐标。**alignment_shadow_summary.tsv** 是对应的 per-locus 计数与 64-bin breadth 汇总。初始结果另存为 `alignment_shadow_initial*.tsv`；rescue 证据保存在对应 `uce_rescue_round_N/`。

**large_files**：进一步过滤时超过深度限制或文件大小限制的 reads。只有产生这类文件时才会出现。

**results**：主组装结果。每个 locus 的最佳 contig 写为 `<locus>.fasta`。

**contigs_all**：组装器评估过的候选 contigs。

**contigs_all_low**：UCE 模式下保留的低支持延伸候选 contigs，仅用于检查，不会提升到主 `results` 目录。

**consensus**：可选 `consensus` 命令生成的一致性序列。

**blast**：可选 `trim` 命令生成的参考切齐结果。

**log.txt**：单样本日志文件。

**result_dict.txt**：每个 locus 的组装状态和 reads 支持摘要。

**ref_reads_count_dict.txt**：第一轮过滤阶段每个 locus 分配到的 reads 数量。

**uce_assembly_summary.csv**：UCE 模式下的单样本逐 locus 组装摘要，包含状态、最佳 contig 长度、reads 支持跨度、read count、read density、支持比例、侧翼平衡度、k-mer 深度指标、候选数和低质量标记。

**uce_rescue_summary.csv**：UCE 模式默认 rescue 阶段生成的单样本首轮至最终结果摘要；使用 `--no-uce-rescue-reads` 时不生成。

**uce_rescue_rounds.csv**：逐轮、逐 locus 记录 active/revert/terminal-side 决策、长度和 unique-read 增量；第二轮还记录左右新增长度、breadth、最大 gap、fragment 数、跨旧 core 边界 fragment 数及是否接受。`reverted_inverted_repeat` 表示该轮新引入达到当前阈值的精确倒置重复；`reverted_unsupported_internal_gap` 表示 rescue 新引入至少 40 bp、没有同一条 read 上连贯 k-mer 链支持的内部区间。两种状态都会逐 locus 恢复该轮备份。

**assembly_graphs**：仅在使用 `--assembler-graph-format gfa`、`dot` 或 `both` 时生成的逐 locus 压缩组装图；默认不生成。

## Marker profiling 输出

**marker_profile/marker_reference_support.tsv**：每条命中的参考序列一行。`hit_queries` 是相容 query 总数；`fractional_queries` 将共享 query 平分给候选；`singleton_queries` 是只与该参考相容的 query 数。

**marker_profile/marker_qc.tsv**：该样本的伪比对与运行参数汇总；计数单位是单条 FASTA/FASTQ query record。

**marker_profile/marker_reference_metadata.tsv**：本次 profiling 所用的 reference ID、Themisto color 与可选 group 注释映射。

## 合并输出

**combined_results**：按 locus 合并不同样本恢复序列后的文件。

**combined_results/aligned**：`combine` 阶段生成的多序列比对结果。

**combined_trimed**：比对列过滤后的逐 locus alignment。默认由 trimAl 生成；使用 `--alignment-filter alifilter` 时由 AliFilter 生成；使用 `--alignment-filter none` 或 `--no-trimal` 时不会生成。

**combined_results.fasta**：比对列过滤前的串联 alignment。

**combined_trimed.fasta**：比对列过滤后的串联 alignment。

**failed_gene_trees.tsv**：溯祖树流程中单基因树构建失败的 locus 列表。

**failed_samples.tsv**：当任一样本在 filter、refilter、assemble 或 UCE rescue 阶段失败时生成。命令行流程会在写出该文件后停止，避免下游步骤在样本结果不完整时继续运行。

**Coalescent.tree**：溯祖流程生成的物种树。

**Concatenation.tree**：串联流程生成的系统树。

## UCE 专用输出

**uce_contigs**：UCE 组装模式生成的 phyluce 兼容 contig 输出。每个样本一个 `*.contigs.fasta` 文件。`sample_name_map.tsv` 记录 TipSeek 样本名（包括保留的兼容 GeneMiner2 旧标识）与 phyluce 安全样本名之间的映射关系。

**uce_rescue_summary.csv**：跨样本合并后的首轮至最终 rescue 摘要。UCE 模式默认生成，使用 `--no-uce-rescue-reads` 时不生成。

**uce_rescue_rounds.csv**：跨样本合并后的逐轮 rescue 审计表。

## Population 输出

`population` 子命令在已有输出目录下创建 `population/`：

- `sample_manifest.tsv`：原始样本名、TipSeek 内部目录名（包括保留的兼容 GeneMiner2 旧名称）、VCF 样本名、reads 路径和 SE/PE 布局的对应关系。
- `reference/population_reference.fasta`：所有样本统一 mapping 使用的公共 UCE 参考。使用 `--population-reference-fasta` 时为复制后的外部参考。
- `reference/population_reference_provenance.tsv`：仅 `pseudoref` 内部构建时，每个 locus 的来源样本、选择策略、候选数、长度及 reads 支持指标。
- `reference/reference_contribution.tsv`：仅 `pseudoref` 内部构建时，每个样本贡献的参考 loci 数及比例，用于检查公共参考是否由少数样本主导。
- `reference/locus_name_map.tsv`：仅 `pseudoref` 内部构建时原始 locus 名和 VCF 安全名称的对应关系。
- `reference/reference_source.tsv`：使用固定外部参考时，记录源文件与复制后的参考路径。
- `mapping/<sample>.bam` 和索引：minibwa 统一 mapping 后经 samtools 处理的 BAM。
- `mapping/mapping_qc.tsv`：每个样本的 mapped/properly-paired reads、mapping rate、覆盖广度和平均深度。
- `variants/cohort.raw.bcf`：联合检测的原始 cohort BCF。
- `variants/cohort.biallelic.snps.vcf.gz`：拆分并保留的双等位 SNP。
- `variants/cohort.genotype_filtered.vcf.gz`：低 DP/GQ 基因型设为缺失后的 VCF。
- `variants/cohort.tagged.vcf.gz`：补充群体统计标签后的 VCF。
- `variants/cohort.filtered.vcf.gz`：经过 QUAL、call rate 和 MAC 过滤的分析起点。
- `variants/variant_qc.tsv`：raw calling、双等位、基因型过滤、标签补充和位点过滤各阶段的位点数与文件路径。
- `structure/all_snps.vcf.gz` 和同名前缀 PLINK/PCA 文件：保留同一 UCE 内多个 SNP 的高敏感度面板。
- `structure/one_snp_per_uce.vcf.gz`、`population.{bed,bim,fam}` 和 `population_pca.*`：每个 UCE 按 call rate、QUAL 和 MAC 选择一个代表 SNP 的主面板。
- `structure/ld_pruned.vcf.gz` 和同名前缀 PLINK/PCA 文件：对全部 SNP 做 LD pruning 的敏感性面板。
- `structure/selected_snps.tsv`：每个 UCE 被选中代表 SNP 的位置和筛选统计。
- `structure/panel_summary.tsv`：三种面板的 SNP 数量、路径和建议用途。
- `structure/admixture/K<K>.log`、`population.<K>.Q` 和 `population.<K>.P`：各 K 的 ADMIXTURE 日志、个体遗传成分和祖先群体等位基因频率。
- `structure/admixture/cv_errors.tsv`：各 K 的交叉验证误差以及最低误差 K 标记。
- `structure/admixture/status.tsv`：ADMIXTURE 的 `complete`、`skipped`、`unavailable` 或 `failed` 状态。

PCA 和 ADMIXTURE 的主解释应优先使用每个 UCE 一个 SNP 的面板，并与 all-SNP 和 LD-pruned PCA 比较。若样本 mapping rate、coverage breadth 明显偏低，应先排查缺失数据；仅在内部构建参考时，才根据 `reference_contribution.tsv` 的参考来源不均衡排查参考偏倚。

PanRefV2还会写入 `reference/panrefv2/`：`index_metadata.tsv`、`recruitment_summary.tsv`、`population_graph.gfa` 和 `locus_summary.tsv`。这些文件记录bait minimizer、read招募、局部图及每个locus的QC；默认只有状态为`pass`的locus进入公共参考。

## 统计输出

`stats` 子命令会从已有输出目录生成类似 HybPiper 的 UCE 恢复统计。

**uce_stats.tsv**：样本级统计表，汇总过滤 reads、成功 loci、低质量 loci、主要失败类型、按参考长度比例统计的恢复数量、rescue 回退数量、总恢复碱基数、contig 长度、reads 支持跨度和 read density。

**uce_locus_stats.tsv**：locus 级统计表，汇总参考长度、样本占有率、恢复长度、reads 支持、侧翼平衡度和候选 contig 数量。

**uce_seq_lengths.tsv**：sample-by-locus contig 长度矩阵，包含根据参考目录计算的 `MeanLength` 行。

**uce_read_counts.tsv**：sample-by-locus 矩阵，记录最佳 contig 的 read count。

**uce_filtered_read_counts.tsv**：sample-by-locus 矩阵，根据 `ref_reads_count_dict.txt` 记录第一轮过滤阶段分配到各 locus 的 reads 数量。

**uce_rescue_stats.tsv**：由已有 `uce_rescue_summary.csv` 转换得到的 rescue 明细表。

**uce_recovery_heatmap.png** 和 **uce_read_counts_heatmap.png**：`stats` 子命令生成的热图；如果环境中没有 `pandas`、`seaborn` 或 `matplotlib`，或使用了 `--stats-no-heatmap`，则不会生成。
