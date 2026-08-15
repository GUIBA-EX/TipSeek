# Supplementary Material for “TipSeek: fragment-aware, evidence-bounded recovery of ultraconserved elements from short reads”

## Supplementary Methods S1. Algorithm parameters

### S1.1 Fragment-aware hierarchical recruitment

快速阶段以 recruitment *k*=23、step=4 和 verification *k*=11 扫描 paired FASTQ。命中任一 mate 后，完整 read pair 作为一个 fragment 存入共享 fragment bank。每个 locus 的正反向参考使用独立 FM-index；maximal exact match 和 run-based orientation 提供参考位置、方向、末端 overhang 及多 locus 证据。

默认 `auto` 模式仅对快速阶段没有 selected fragments 的 loci 再扫描 FASTQ，参数为 recruitment *k*=21、step=1 和 verification *k*=19。通过 unresolved-locus 粗门控的 fragments 返回完整 probe 面板重新扩展候选。新增 read pair 需以至少 45 bp、80% identity 的局部比对唯一支持一个未恢复 locus。

fallback-only contig 首先作为 provisional core 接受完整 probe 面板复核。候选长度需至少为 200 bp，目标 probe coverage 和 identity 均需至少为 80%。如果另一 locus 在相同 coverage/identity 条件下达到目标得分的 95%，该候选标记为 ambiguous。至少 150 bp 的精确长倒置重复触发结构拒绝；最大内部无支持 gap 达到 40 bp 的候选标记为 `anchored_with_review`，保留输出但不进入 rescue。

### S1.2 Topology-aware evidence budgeting

自适应证据预算在以下条件同时满足时启动：候选数不少于 512、估计深度高于 160×、精确 seeds 覆盖至少 48/64 个参考区间。其他 loci 保留全部合格 fragments；预先触发文件大小保护时使用 legacy selector。

设合格 fragment 数为 *N*、有效参考长度为 *L*<sub>eff</sub>、平均 fragment 碱基数为 *b̄*，核心预算为

`N_core = min[N, max(512, ceil(80L_eff / b̄), ceil(0.60N))]`。

核心候选依次按较少的 locus assignments、较长 maximal exact match、更多 aligned mates 和稳定 fragment ID 排序，再以 64-bin 配额维持参考跨度。terminal fragments 按左右侧和 overhang 长度分层，每个长度档最多保留 4 个 fragments，每侧总上限为 768。核心和末端集合去重合并后写入组装输入。

### S1.3 PE-supported dual-graph assembly

`uce-rust` 使用同一组 k-mer counts 构建 core graph 和 PE-assisted graph。core graph 使用常规深度及参考证据；PE-assisted graph 可接纳深度不高于默认错误阈值的非参考 k-mers，但每个此类 k-mer 需由至少两个独立 fragments 支持。paired-fragment support 仅累计到真实分支边。默认 backbone 策略对每个分支执行最多 24 步前瞻。

两张图从相同 seeds 独立生成候选。候选依次按有效支持碱基数、support breadth、unique-read density、无支持 gap 比例、长度、unique/total read 数、flank balance 和图权重排序。缺少 unique-read 或位置支持的候选被拒绝。PE-assisted 最优候选短于已通过 QC 的 core-graph 候选时，最终结果保留 core path。

### S1.4 Reversible rescue

UCE 模式默认以固定 *k*=21 执行一轮 rescue。第一轮将长度至少为 60 bp 的已接受 contigs 与原参考共同作为 bait，重新招募 reads 并组装。对于 rescue 前已接受的 loci，新结果需保持至少 50% 的原 unique-read density，且不得新增至少 40 bp 的内部无支持 gap。新增至少 150 bp 的精确长倒置重复触发该 locus 回滚。

可选第二轮仅处理第一轮后仍在增长的 loci，并使用 350 bp 末端窗口构建 baits。每侧延伸需至少为 30 bp、支持 breadth 至少为 85%、最大无支持 gap 不超过 30 bp，并由至少两个独立 fragments 支持，其中至少一个 fragment 跨越第一轮 core 与新增序列的边界。未达到条件的一侧修剪至第一轮 core。

rescue 对每个 locus 独立执行 accept、trim 或 revert。fast/fallback 来源、fragment 预算、provisional-core 判定、rescue 轮次和回滚结果分别写入 TSV/CSV 文件；probe 拒绝和结构拒绝结果保留在独立目录中。

## Supplementary Methods S2. Reproducibility

### S2.1 Software and archived release

- TipSeek version：X
- Git commit：X
- Archived release DOI：X
- Source repository：X
- License：GPL-3.0-or-later
- SPAdes version：X
- PHYLUCE version：X
- GeneMiner2 version and commit：X

### S2.2 Input data

- Samples：CRR2698935、CRR2698936、CRR2698938–CRR2698946
- Excluded from the prespecified analysis set：CRR2698937（reason：X）
- Probe panel：3,023 loci
- BioProject：X
- Probe-set permanent URL/DOI：X
- Input FASTQ checksums：X

### S2.3 Compute environment

- Operating system：X
- CPU model：X
- Physical cores allocated：72
- RAM：X
- Compression backend：zlib-ng
- File system/storage：X
- Thread and memory limits for each workflow：X

### S2.4 Commands and outputs

- TipSeek *k*=23, R0：X
- TipSeek *k*=23, R1：X
- TipSeek *k*=23, R2：X
- TipSeek *k*=31, R1：X
- GeneMiner2 *k*=23：X
- SPAdes + PHYLUCE genome-harvesting：X
- Benchmark output archive：X
- Figure/table generation scripts：X

## Supplementary Table S1. Complete coral benchmark

| Configuration | Sample–locus recoveries | Mean per sample | Loci recovered in all 11 samples | Median length (bp) | Wall time (min) | Peak RSS (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GeneMiner2 *k*=23 | 26,528 | 2,411.64 | 1,653 | 281 | 10.57 | 0.51 |
| TipSeek *k*=23, R0 | 25,101 | 2,281.91 | 1,125 | 439 | 6.88 | 0.72 |
| TipSeek *k*=23, R1 | 25,933 | 2,357.55 | 1,471 | 771 | 28.02 | 6.25 |
| TipSeek *k*=23, R2 | 25,933 | 2,357.55 | 1,471 | 990 | 61.90 | 7.14 |
| TipSeek *k*=31, R1 | 22,085 | 2,007.73 | 833 | 752 | 19.84 | 2.83 |
| SPAdes + PHYLUCE genome-harvesting | 24,231 | 2,202.82 | 1,084 | 2,072 | — | — |

R0, R1 and R2 denote zero, one and two rescue rounds, respectively. PHYLUCE resource values are not reported because the complete workflow lacks a same-run, same-core timing and memory record.
