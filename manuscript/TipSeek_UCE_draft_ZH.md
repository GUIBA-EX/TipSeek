# TipSeek：片段感知、证据约束的超保守元件恢复

# TipSeek: fragment-aware, evidence-bounded recovery of ultraconserved elements from short reads

**Original Paper | Genome analysis**

**Authors:** X<sup>1,*</sup>

<sup>1</sup> X

<sup>*</sup>Corresponding author. X. E-mail: X.

## Abstract

**Motivation:** 超保守元件（ultraconserved elements，UCEs）已成为动物系统发育、物种界定和分类修订的重要基因组标记。其保守核心用于跨类群识别同源 locus，而变异更丰富的侧翼为浅层系统发育推断提供主要信息。对于低深度短 reads，通用参考比对可以招募参考已表示区域的 reads，却不能独立重建参考边界外的样本特异侧翼，并且保守核心可能使同一 fragment 同时支持多个 loci。现有全局组装和目标招募流程因此分别面临吞吐量与侧翼恢复的权衡。TipSeek 将 UCE recovery 建模为连续的 fragment-level evidence propagation：仅对未恢复 loci 提高招募敏感度，在完整参考面板上核验 locus 支持，为饱和核心和稀疏末端分别分配证据预算，并通过 core/PE-assisted 双图与逐 locus 回滚控制延伸。

**Results:** 在 11 个珊瑚 WGS 样本和 3,023-locus OCTO-V2 UCE/exon 面板上，TipSeek 默认配置获得 25,933 个 sample–locus recoveries，候选序列中位长度为 771 bp；GeneMiner2 *k*=23 获得 26,528 个 recoveries，中位长度为 281 bp；SPAdes + PHYLUCE 获得 24,231 个 recoveries，中位长度为 2,072 bp。TipSeek 第一轮 rescue 同时提高跨样本恢复完整度和候选长度，第二轮则主要延长已有候选。因此，*k*=23 与一轮 rescue 被设为当前默认配置。

**Availability and implementation:** TipSeek 以 Rust 实现，采用 GPL-3.0-or-later 许可证，源代码见 https://github.com/GUIBA-EX/TipSeek。

**Keywords:** ultraconserved elements; target-restricted assembly; genome skimming; paired-end reads; phylogenomics; coral

## 1 Introduction

超保守元件（ultraconserved elements，UCEs）已广泛用于动物系统发育、物种界定和分类修订。其高度保守的核心可在较深进化尺度上锚定对应的同源 locus，而共同恢复的侧翼通常积累更多变异，是解析近缘种和种内谱系关系的主要信息来源（Faircloth et al. 2012, Smith et al. 2014, Erickson et al. 2021）。因此，UCE recovery 的目标不只是检测到 probe 或保守核心，而是从每个样本中恢复 locus 归属明确、由 reads 支持并包含足够侧翼的连续序列。

当前非模式动物的 UCE 系统发育研究大量依赖 Illumina 低深度 WGS 或 genome-skimming 数据。通用短读比对器可以快速招募与参考足够相似、且位于参考已表示区域内的 reads，但单独比对并不能重建参考边界以外的样本特异序列。随着样本与参考之间的分化增加，变异较快的侧翼 reads 也比保守核心 reads 更容易发生比对脱落和参考偏倚。与此同时，不同 UCE loci 的保守区可能共享短序列，使同一 read 或 paired fragment 同时命中多个 loci；仅依赖单条 read 的最高比对分数不足以稳定恢复直系同源位点。低深度短 reads 的 UCE recovery 因而通常需要全局组装，或在参考引导招募后进行逐 locus 组装，并利用 paired-end linkage 与局部重叠关系从核心继续恢复侧翼（Hahn et al. 2013, Allen et al. 2015）。

现有工作流主要采用 assembly-first 或 recruitment-first 两类路径。PHYLUCE 通常先使用 SPAdes 对全样本 reads 构建 contigs，再按 probe 识别和提取目标序列（Bankevich et al. 2012, Faircloth 2016）；该路径能够利用全局重叠关系恢复较长 contigs，并已用于从珊瑚 genome-skimming 数据中提取 OCTO-V2 UCE 和 exon loci（Quattrini et al. 2024），但同时需要处理大量非目标 reads，增加计算、存储和中间文件开销，并限制大批量浅层测序样本的分析吞吐量。GeneMiner 和 GeneMiner2 先招募目标相关 reads，再进行逐 locus 组装，从而缩小搜索空间并提高目标恢复吞吐量（Xie et al. 2024, Yu et al. 2026）。对于 UCE 数据，这一路径还需要处理一个特有的不平衡：保守核心可迅速积累高深度、跨 locus 共享的证据，而连接核心与侧翼的 fragments 数量少、覆盖低且位置分散。不同流程如何扩大招募、处理多 locus 命中、压缩核心冗余和保留末端重叠链，会共同改变恢复量、候选长度和计算成本（Bossert et al. 2024）。

TipSeek 针对这一核心—侧翼不平衡，将 UCE recovery 表述为连续的 fragment-level evidence propagation。分层招募只对未恢复 loci 提高敏感度，并在完整参考面板上重新核验每个 fragment 的 locus 支持；结构感知的证据预算优先压缩已饱和的核心覆盖，同时保护具有左右末端延伸潜力的 fragments；参考锚定的 core graph 和 PE-assisted graph 再结合 paired-end 分支支持、局部深度和有限前瞻生成候选，rescue 则按 locus 独立接纳、裁切或回滚新增序列。本次发布聚焦 UCE recovery 模块，描述其 fragment 招募、panel-wide locus verification、核心—末端证据选择和可逆延伸机制，并在 11 个珊瑚 WGS 样本上比较 TipSeek、GeneMiner2 和 PHYLUCE 的恢复完整度、候选长度与计算资源，同时通过珊瑚真值数据检验 locus assignment 和侧翼恢复的正确性。

## 2 Materials and methods

### 2.1 Fragment-aware hierarchical recruitment and panel-wide locus verification

TipSeek 将同一 read pair 视为一个 fragment。任一 mate 通过粗招募后，R1 和 R2 作为不可分割的单元在共享 fragment bank 中仅存储一次；两个 mates 共同提供 locus、链方向、参考位置和末端延伸证据。UCEFilter 使用 canonical rolling k-mers、blocked Bloom filter 和精确哈希表完成粗招募，并为每个 locus 的正向与反向参考建立 FM-index，以最长精确匹配和 run-based orientation 核验候选（Bloom 1970, Ferragina and Manzini 2000, Putze et al. 2009）。因此，k-mer 命中只负责产生候选集合，最终保留由 fragment 层面的双端证据决定。

默认 `auto` 模式首先执行快速招募，仅对没有 selected fragments 的 loci 再扫描 FASTQ。敏感阶段使用的 locus 子集只承担粗招募门控；任一 fragment 通过该门控后，其候选集合会重新扩展到完整 probe/reference 面板，再进行方向、精确匹配和可选局部比对核验。该步骤避免因预先缩小参考集合而将共享核心 reads 人为判为 locus 特异证据。fallback-only contig 还需通过 probe coverage、identity、近似并列 locus、长倒置重复和内部无支持 gap 检查，才能成为 rescue-eligible core；未通过接纳条件的候选记录为 ambiguous、review 或 rejected，不进入后续延伸。完整阈值见 Supplementary Methods S1。

TipSeek 不以单条 read 的最高比对分数强制完成唯一 locus 分配，而是为每个 fragment–locus 组合建立确定性的证据向量。对每个候选 locus，程序分别提取两个 mates 的最佳精确 seed，记录参考覆盖区间、最长精确匹配长度（`max_exact`）、左右末端延伸量和直接命中的 mate 数；run-based orientation 用于排除两个 mates 均无有效方向证据或具有相同明确链方向的配对。若启用局部比对门控，还要求至少一个 mate 达到设定的参考 overlap 和 identity。完成全 panel 核验后，`locus_count` 记录该 fragment 同时通过核验的 locus 数。候选 locus 数超过配置上限时，该 fragment 被整体排除；未超过上限的多 locus fragment 可继续为相应 loci 提供候选证据，而不是被任意指定给单一 best hit。

在后续的饱和 locus 选择中，每个 fragment–locus 候选按字典序排序：首先优先保留 `locus_count` 较小、即 locus 特异性更高的 fragments，其次依次比较更长的 `max_exact`、更多直接命中的 mates 和稳定的 fragment ID。该顺序同时用于核心和末端候选的质量比较，末端候选另优先考虑延伸量。由此，多 locus 命中、双端支持和确定性 tie breaking 在招募与证据选择之间保持一致。

### 2.2 Structural fragment labeling and topology-aware evidence budgeting

TipSeek 中的 core-supporting 和 terminal-supporting evidence 不是两个互斥的 read 集合。通过 locus 核验的 fragment–locus 候选首先作为参考锚定的核心证据；terminal 是附加在同一候选上的左右端结构标签。对每个 mate，程序根据其在该 locus 上的最佳精确 seed 计算相对于匹配参考序列的原始方向位置。设参考长度为 *L*<sub>ref</sub>、用户设置的末端窗口为 *W*，有效末端窗口为 *W*<sub>eff</sub> = min(*W*, *L*<sub>ref</sub>/5)。当 seed 的参考起点位于左侧 *W*<sub>eff</sub> 内，且 read 在 seed 左侧仍有碱基时，该 fragment 被标记为 left-terminal；当 seed 的参考终点位于右侧 *W*<sub>eff</sub> 内，且 read 在 seed 右侧仍有碱基时，被标记为 right-terminal。左右 `extension` 分别取两个 mates 中对应未包含于精确 seed 的最大碱基数。没有 terminal 标签的候选仅提供核心覆盖；具有 terminal 标签的候选同时提供核心锚定和末端延伸证据，完整 paired fragment 中未直接命中的 mate 仍被保留以维持 core-to-flank linkage。

证据压缩只在参考跨度充分且覆盖已饱和的 loci 上启动。设通过既有 exact-match 门控的 fragment 数为 *N*、有效参考长度为 *L*<sub>eff</sub>、平均 fragment 碱基数为 *b̄*。当 *N* 不少于 512、估计覆盖深度高于 160×，且 exact seeds 覆盖至少 48/64 个参考区间时，核心目标预算为

`N_core = min[N, max(512, ceil(80L_eff / b̄), ceil(0.60N))]`。

不满足饱和条件的 loci 保留全部合格 fragments；若既有深度或文件大小限制已触发，则沿用动态 exact-match selector。对于进入自动压缩的 loci，核心候选按 `locus_count`、`max_exact`、直接命中 mate 数和稳定 fragment ID 排序，并先以 64-bin 配额在参考范围内分散保留，再按同一质量顺序补足 `N_core`，避免高深度局部堆积取代全参考覆盖。

左右 terminal candidates 分别建立保护队列。每侧首先按 extension 由长到短排序，再依次比较更小的 `locus_count`、更长的 `max_exact`、更多直接命中的 mates 和稳定 fragment ID；每个精确 extension 长度最多保留 4 个 fragments，每侧最多保留 768 个。核心选择与左右末端保护结果最终按 fragment ID 取并集，因此同一 fragment 即使同时携带核心和末端标签也只写入一次。该策略在压缩饱和核心冗余的同时，保留从参考核心通向左右侧翼的低频重叠阶梯。

### 2.3 PE-supported dual-graph assembly and reversible rescue

经核心预算与末端保护合并后的去重 fragments 进入 `uce-rust`。程序从同一组 k-mer counts 构建两张参考锚定的加权 de Bruijn 图（Idury and Waterman 1995）。core graph 使用常规深度和参考证据；PE-assisted graph 还可接纳由至少两个独立 fragments 支持的低深度非参考 k-mers。paired-fragment support 只累计到真实分支边，并与参考连续性、局部深度和有限前瞻共同决定路径。两张图独立产生候选；如果 PE-assisted 候选短于已通过 QC 的 core 候选，最终结果保留 core path。

默认 rescue 使用已接受 core 与原参考重新招募 reads，并按 locus 独立执行 accept、trim 或 revert。该步骤采用 target-restricted iterative assembly 的 baiting 思路（Hahn et al. 2013, Allen et al. 2015）。第一轮依据 unique-read density、内部无支持 gap 和长倒置重复控制 whole-contig extension；可选第二轮只处理仍在增长的 loci，并要求新增末端同时满足支持 breadth、gap、独立 fragment 和 core-to-extension bridge 条件。任一 locus 的失败只回滚该 locus。TipSeek 将 fast/fallback 来源、full-panel locus support 汇总、fragment 预算、terminal candidate 计数、probe gate、rescue 轮次和回滚结果写入结构化汇总与逐轮审计表，使最终候选可追溯至招募来源和接纳规则。完整门控参数见 Supplementary Methods S1，整体证据流见 Figure 1。

![Figure 1. TipSeek algorithm](figures/figure1_algorithm.svg)

**Figure 1. TipSeek 的 fragment-aware、evidence-bounded UCE recovery 流程。** Paired reads 以 fragment 为单位经过快速招募和 unresolved-only fallback；通过子集门控的 fragments 返回完整 probe/reference 面板进行 locus 核验，并记录多 locus 支持、mate linkage、方向和参考位置。通过核验的候选进一步标记为 core-only、left-terminal 和/或 right-terminal evidence；核心预算与左右末端保护结果去重合并后进入 core graph 和 PE-assisted graph。rescue 对每个 locus 独立执行 accept、trim 或 revert，虚线表示逐 locus 回滚。

### 2.4 Software implementation

TipSeek 以 Rust 实现，统一入口为 `tipseek`。本文分析使用的 UCE 路径为 `UCEFilter fast pass → unresolved-only fallback with full-panel verification → topology-aware selection → uce-rust → one-round rescue`。参考引导过滤、FM-index 和加权 de Bruijn 图是算法基础；TipSeek 的方法贡献是 fragment-level 状态传递、panel-wide multi-locus evidence、核心与末端的分离预算、PE/core 双候选路径及逐 locus 可逆延伸。完整参数和复现命令见 Supplementary Methods S2。

### 2.5 Benchmark data and comparison design

benchmark 使用 You et al.（2026）公开的 12 个八放珊瑚 WGS 样本中的 11 个（CRR2698935、CRR2698936、CRR2698938–CRR2698946）；CRR2698937 因 X 未纳入预设分析集。参考面板为 OCTO-V2，由 29,181 条 probes 组成，共靶向 3,023 个 loci，包括 1,337 个 UCE loci 和 1,686 个 exon loci（Erickson et al. 2021）。主比较包括 GeneMiner2 *k*=23、TipSeek 默认配置（*k*=23、`auto`、一轮 rescue）以及 SPAdes X + PHYLUCE X genome-harvesting。TipSeek 内部比较包括 *k*=23 下的零、一和两轮 rescue，以及 *k*=31、`auto`、一轮 rescue。完整版本、命令和运行环境见 Supplementary Methods S2。

主要指标为 sample–locus recoveries、每样本候选数、面板独立位点数、shared loci、候选序列中位长度、wall time 和峰值 RSS。一个样本在一个 locus 获得被相应工作流接受的候选，计为一次 sample–locus recovery；面板独立位点指至少在一个样本中恢复的 loci，shared loci 分别按全部 11 个样本和至少 9 个样本恢复统计；长度中位数由全部 accepted sample–locus candidates 计算。GeneMiner2 与 TipSeek 使用同一批样本直接计时。PHYLUCE 的主分析报告恢复量和候选长度，其资源记录不参与比较。

### 2.6 Coral ground-truth validation

为区分候选长度与序列正确性，使用具有公开参考基因组的珊瑚 X 建立真值数据集。首先将 3,023-locus 面板与参考基因组比对，仅保留唯一定位且包含完整目标区间的 loci；随后以固定随机种子 X 模拟与实测数据一致的 paired-end reads，并设置低、中和高三个测序深度（X、X 和 X）。TipSeek、GeneMiner2 和 SPAdes + PHYLUCE 使用与 WGS benchmark 相同的参数运行。每条 accepted candidate 依据其最佳参考位置和覆盖范围判定 locus assignment，并统计正确 locus assignment 比例、序列一致性、嵌合或错误延伸比例及正确侧翼恢复长度。真值定义、模拟命令和逐 locus 结果见 Supplementary Methods S3 和 Supplementary Table S3。

## 3 Results

### 3.1 Recovery completeness, sequence length and resource use

GeneMiner2 获得 26,528 个 sample–locus recoveries，其中 1,653 个 loci 在全部 11 个样本中恢复；TipSeek 默认配置分别为 25,933 和 1,471。GeneMiner2 比 TipSeek 多获得 595 个 sample–locus recoveries（2.24%），TipSeek 的候选序列中位长度为 771 bp，是 GeneMiner2 的 2.74 倍。TipSeek 和 GeneMiner2 的 wall time 分别为 28.02 和 10.57 min，峰值 RSS 分别为 6.25 和 0.51 GiB。

SPAdes + PHYLUCE 获得 24,231 个 sample–locus recoveries，其中 1,084 个 loci 在全部样本中恢复，候选序列中位长度为 2,072 bp。相较 TipSeek，PHYLUCE 少获得 1,702 个 sample–locus recoveries，候选序列中位长度增加 1,301 bp。PHYLUCE 的比较限于恢复量和候选长度。

**Table 1. 三个工作流的五个配置在 11 个珊瑚 WGS 样本上的结果。**

| 配置 | Sample–locus recoveries | 每样本平均/中位数（范围） | 面板独立位点 | Shared loci（11/11） | Shared loci（≥9/11） | 长度中位数（bp） | Wall time（min） | 峰值 RSS（GiB） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| GeneMiner2 *k*=23 | 26,528 | 2,411.64/2,384（2,364–2,546） | 2,909 | 1,653 | 2,250 | 281 | 10.57 | 0.51 |
| TipSeek *k*=23，R0 | 25,101 | 2,281.91/2,262（2,173–2,464） | 2,881 | 1,125 | 2,050 | 439 | 6.88 | 0.72 |
| TipSeek *k*=23，R1 | 25,933 | 2,357.55/2,341（2,259–2,500） | 2,885 | 1,471 | 2,134 | 771 | 28.02 | 6.25 |
| TipSeek *k*=23，R2 | 25,933 | 2,357.55/2,341（2,259–2,500） | 2,885 | 1,471 | 2,134 | 990 | 61.90 | 7.14 |
| SPAdes + PHYLUCE genome-harvesting | 24,231 | 2,202.82/2,224（1,985–2,330） | 2,718 | 1,084 | 1,984 | 2,072 | — | — |

*Note:* R0、R1 和 R2 分别表示零、一和两轮 rescue。PHYLUCE 缺少同次、同核的完整计时与内存记录，因此不比较资源消耗。

TipSeek 内部比较确定了默认运行点。在一轮 rescue 下，*k*=23 比 *k*=31 增加 3,848 个 sample–locus recoveries（17.42%），候选序列中位长度分别为 771 和 752 bp。*k*=23 从 R0 到 R1 增加 832 个 sample–locus recoveries，并将中位长度从 439 bp 增至 771 bp；R2 不再增加 loci，将中位长度提高到 990 bp，同时把 wall time 从 28.02 min 增至 61.90 min。因此，*k*=23 和一轮 rescue 被设为当前默认值：第一轮同时增加恢复量和候选长度，第二轮主要延伸已有候选。完整六组结果见 Supplementary Table S1。

该 benchmark 直接量化恢复量、候选长度和计算资源。候选长度表示各工作流恢复的连续序列范围，其 locus assignment、碱基一致性和正确侧翼长度由珊瑚真值数据进一步检验。TipSeek 同时为每个候选保留从 fragment 招募、full-panel locus support 和证据预算到 rescue 接纳或回滚的结构化记录。Figure 2 汇总了五个配置的恢复量、shared loci、候选长度及 TipSeek 与 GeneMiner2 的资源运行点。

![Figure 2. Coral benchmark](figures/figure2_benchmark.svg)

**Figure 2. 三个工作流在 11 个珊瑚 WGS 样本上的恢复结果与计算权衡。** (A) sample–locus recoveries；(B) 在全部 11 个样本中恢复的 shared loci；(C) accepted candidates 的序列长度中位数；(D) GeneMiner2 与 TipSeek 配置在 wall time–峰值 RSS 平面上的运行点。白点和橙色外圈标记当前默认 TipSeek R1。SPAdes + PHYLUCE 因缺少同次、同核的完整记录而不进入资源面板。

### 3.2 Sequence correctness on the coral ground-truth dataset

在珊瑚真值数据集的低、中和高深度条件下，三个工作流分别获得 X、X 和 X 个可评分候选。TipSeek 默认配置的正确 locus assignment 比例为 X，候选序列一致性为 X，嵌合或错误延伸比例为 X，正确侧翼恢复长度中位数为 X bp；GeneMiner2 和 SPAdes + PHYLUCE 的对应结果见 Supplementary Table S3。该结果用于检验更长候选是否来自目标 locus 的受支持侧翼，而不将候选长度本身解释为准确性。

## 4 Discussion

本研究将 UCE recovery 与普通参考比对区分为两个不同任务。短读比对可以高效识别参考已表示且与样本足够相似的区域，但 UCE 系统发育分析还要求恢复变异更快的侧翼，并在保守核心可能跨 locus 共享的条件下维持正确的位点归属。TipSeek 因此将 mapper-like seed recruitment 作为入口，而不是终点：完整 panel 核验负责描述每个 fragment 的候选 locus 范围，paired-end linkage 和局部重叠随后用于逐 locus 组装，末端证据则推动路径越过 probe/reference 边界进入样本自身的侧翼。

三个工作流体现了不同的算法取舍。GeneMiner2 获得最多的 sample–locus recoveries 和 shared loci，但候选序列中位长度最短；SPAdes + PHYLUCE 通过全局组装获得最长候选，但恢复量和跨样本完整度较低；TipSeek R1 位于二者之间。具体而言，相较 GeneMiner2，TipSeek R1 少恢复 595 个 sample–locus recoveries、24 个面板独立位点和 182 个 11/11 shared loci，但候选序列中位长度增加 490 bp。相较 PHYLUCE，TipSeek R1 多恢复 1,702 个 sample–locus recoveries、167 个面板独立位点和 387 个 11/11 shared loci，而候选序列中位长度短 1,301 bp。该结果表明，recruitment-first、assembly-first 和 fragment-aware target-restricted 路径分别占据恢复广度、序列延伸和计算成本之间的不同运行区间（Bossert et al. 2024）。

R0 到 R1 的变化说明第一轮 rescue 主要补全已有面板位点在不同样本中的恢复，并同步延长其侧翼。R1 比 R0 增加 832 个 sample–locus recoveries，而面板独立位点仅增加 4 个；与此同时，11/11 shared loci 增加 346 个，≥9/11 shared loci 增加 84 个，候选序列中位长度由 439 bp 增至 771 bp。第一轮 rescue 因而不是简单扩大可检出的参考 locus 范围，而是利用已接受 core 重新招募与之连接的 fragments，补回部分样本中的缺失结果并继续延伸 reads 支持的序列。

第二轮 rescue 的作用更集中于延伸已有候选。R2 与 R1 的 sample–locus recoveries、面板独立位点和 shared loci 完全相同，仅将候选序列中位长度由 771 bp 提高到 990 bp，同时使 wall time 由 28.02 min 增至 61.90 min、峰值 RSS 由 6.25 GiB 增至 7.14 GiB。因此，R0 适合作为资源优先的快速运行点，R1 同时改善跨样本完整度和序列长度，R2 则适用于优先扩展已恢复 loci 的分析；*k*=23 与一轮 rescue 据此作为当前默认配置。

TipSeek 的主要算法特征是 fragment evidence 在招募、选择、图路径和 rescue 之间连续传递。对于同时命中多个 loci 的 fragment，程序不依赖单一 best-hit 分数强制唯一归属，而是在完整 panel 上记录候选 locus 数，并在证据排序中优先选择 locus 特异性更高、精确匹配更长且双端支持更充分的 fragments。core 与 terminal 也不是互斥类别：terminal 标签在核心锚定证据之上标记左右端延伸潜力，并通过独立保护队列进入最终 fragment 并集。由此，高深度核心可以被压缩，而连接 core 与 flank 的低频 fragments 不会随深度截断一起丢失。第一轮 rescue 对 shared loci 和候选长度的同步提升，以及第二轮只继续增加长度的结果，与这一分阶段 evidence-bounded recovery 机制一致。

珊瑚真值数据进一步检验完整 panel 核验和末端延伸是否能够恢复正确的 locus 与侧翼。在低、中和高深度条件下，TipSeek 的 locus assignment、序列一致性、错误延伸和正确侧翼恢复结果分别为 X；GeneMiner2 和 SPAdes + PHYLUCE 的对应结果为 X。该验证与实测 WGS benchmark 共同用于确定默认配置的运行点，并区分“恢复得更长”与“恢复了正确侧翼”这两个指标。

## Acknowledgements

X

## Author contributions

X（按 CRediT taxonomy 填写）。

## Supplementary material

Supplementary Material 包含算法阈值、软件版本、完整命令、运行环境、样本与探针元数据、六组实测 benchmark 结果及珊瑚真值验证结果。

## Conflict of interests

X

## Funding

X

## Data availability

TipSeek 源代码见 https://github.com/GUIBA-EX/TipSeek，并采用 GPL-3.0-or-later 许可证。珊瑚 benchmark reads 为 CRR2698935、CRR2698936、CRR2698938–CRR2698946，来自 BioProject PRJCA057506（You et al. 2026）。OCTO-V2 bait set 归档于 https://doi.org/10.6084/m9.figshare.12061038。

## Ethics statement

本研究分析公开测序数据，未产生新的生物样本或测序数据，伦理审批不适用。

## AI disclosure

OpenAI Codex 用于协助手稿结构重组和语言编辑。作者核查了全部技术描述、数据、引用与结论，并对稿件内容负责。

## References

Allen JM, Huang DI, Cronk QC, Johnson KP. aTRAM—automated target restricted assembly method: a fast method for assembling loci across divergent taxa from next-generation sequencing data. *BMC Bioinformatics* 2015;16:98. https://doi.org/10.1186/s12859-015-0515-2

Bankevich A, Nurk S, Antipov D, et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. *Journal of Computational Biology* 2012;19:455–477. https://doi.org/10.1089/cmb.2012.0021

Bloom BH. Space/time trade-offs in hash coding with allowable errors. *Communications of the ACM* 1970;13:422–426. https://doi.org/10.1145/362686.362692

Bossert S, Pauly A, Danforth BN, Orr MC, Murray EA. Lessons from assembling UCEs: a comparison of common methods and the case of *Clavinomia* (Halictidae), an Old World member of the tribe Dieunomiini. *Molecular Ecology Resources* 2024;24:e13925. https://doi.org/10.1111/1755-0998.13925

Erickson KL, Pentico A, Quattrini AM, McFadden CS. New approaches to species delimitation and population structure of anthozoans: two case studies of octocorals using ultraconserved elements and exons. *Molecular Ecology Resources* 2021;21:78–92. https://doi.org/10.1111/1755-0998.13241

Faircloth BC. PHYLUCE is a software package for the analysis of conserved genomic loci. *Bioinformatics* 2016;32:786–788. https://doi.org/10.1093/bioinformatics/btv646

Faircloth BC, McCormack JE, Crawford NG, Harvey MG, Brumfield RT, Glenn TC. Ultraconserved elements anchor thousands of genetic markers spanning multiple evolutionary timescales. *Systematic Biology* 2012;61:717–726. https://doi.org/10.1093/sysbio/sys004

Ferragina P, Manzini G. Opportunistic data structures with applications. In: *Proceedings of the 41st Annual Symposium on Foundations of Computer Science*. Redondo Beach, CA: IEEE Computer Society, 2000, 390–398. https://doi.org/10.1109/SFCS.2000.892127

Hahn C, Bachmann L, Chevreux B. Reconstructing mitochondrial genomes directly from genomic next-generation sequencing reads—a baiting and iterative mapping approach. *Nucleic Acids Research* 2013;41:e129. https://doi.org/10.1093/nar/gkt371

Idury RM, Waterman MS. A new algorithm for DNA sequence assembly. *Journal of Computational Biology* 1995;2:291–306. https://doi.org/10.1089/cmb.1995.2.291

Putze F, Sanders P, Singler J. Cache-, hash-, and space-efficient Bloom filters. *ACM Journal of Experimental Algorithmics* 2009;14:4.4–4.18. https://doi.org/10.1145/1498698.1594230

Quattrini AM, McCartin LJ, Easton EE, et al. Skimming genomes for systematics and DNA barcodes of corals. *Ecology and Evolution* 2024;14:e11254. https://doi.org/10.1002/ece3.11254

Smith BT, Harvey MG, Faircloth BC, Glenn TC, Brumfield RT. Target capture and massively parallel sequencing of ultraconserved elements for comparative studies at shallow evolutionary time scales. *Systematic Biology* 2014;63:83–95. https://doi.org/10.1093/sysbio/syt061

Xie P, Guo Y, Teng Y, Zhou W, Yu Y. GeneMiner: a tool for extracting phylogenetic markers from next-generation sequencing data. *Molecular Ecology Resources* 2024;24:e13924. https://doi.org/10.1111/1755-0998.13924

You L, Xia F, Liu X. A new gorgonian *Paraplexaura binyuani* sp. nov. (Cnidaria, Octocorallia, Acanthogorgiidae) from the Huaguang Atoll, Xisha Islands, South China Sea. *Diversity* 2026;18:166. https://doi.org/10.3390/d18030166

Yu X, Tang Z, Zhang Z, Song Y, He H, Shi Y, Hou J, Yu Y. GeneMiner2: accurate and automated recovery of genes from genome skimming data. *Molecular Ecology Resources* 2026;26:e70111. https://doi.org/10.1111/1755-0998.70111
