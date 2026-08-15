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

- TipSeek release：v1.6.2（https://github.com/GUIBA-EX/TipSeek/releases/tag/v1.6.2）
- Benchmark Git commit：d46ab9d
- Archived release DOI：X
- Source repository：https://github.com/GUIBA-EX/TipSeek
- License：GPL-3.0-or-later
- SPAdes version：X
- PHYLUCE version：X
- GeneMiner2 version and commit：X

### S2.2 Input data

- Samples：You et al.（2026）公开的 12 个八放珊瑚 WGS 样本中的 11 个（CRR2698935、CRR2698936、CRR2698938–CRR2698946）
- Excluded from the prespecified analysis set：CRR2698937（reason：X）
- Probe panel：OCTO-V2，29,181 probes targeting 3,023 loci（Erickson et al. 2021）
- Probe composition：1,337 UCE loci and 1,686 exon loci
- Probe/reference sequence version：OCTO-V2（实际下载文件名与校验和：X）
- BioProject：PRJCA057506（You et al. 2026）
- Probe-set permanent URL/DOI：https://doi.org/10.6084/m9.figshare.12061038
- Input FASTQ checksums：X
- Per-sample taxon, library and sequencing metadata：Supplementary Table S2

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

## Supplementary Methods S3. Coral ground-truth validation

### S3.1 Reference genome and truth-locus definition

- Coral species：X
- Reference assembly accession and version：X
- Reference assembly DOI/URL：X
- Probe-to-reference aligner, version and command：X
- Unique-locus criteria：X
- Number of uniquely mapped truth loci：X

仅将满足预设唯一定位、覆盖和 identity 条件的 probe targets 纳入评分。每个 truth locus 保留目标区间及其两侧 X bp 参考序列，用于分别评估 locus assignment、目标区间一致性和侧翼延伸。

### S3.2 Read simulation and workflow execution

- Simulator and version：X
- Fixed random seed：X
- Read length：X bp
- Insert-size distribution：X
- Sequencing depths：X、X、X
- Read-error model：X
- Simulated FASTQ checksums：X
- TipSeek, GeneMiner2 and SPAdes + PHYLUCE commands：X

三个工作流沿用实测 WGS benchmark 的核心参数。除方法本身产生的随机步骤外，输入 reads、truth loci 和评分脚本完全相同。

### S3.3 Accuracy metrics

accepted candidate 首先与完整 truth reference 比对。最佳命中属于预期 locus 且满足 X coverage 和 X identity 时，记为正确 locus assignment。目标区间的碱基一致性按对齐的匹配碱基比例计算；candidate 中跨越两个非相邻 truth intervals 或包含超过 X bp 的高一致性非目标片段时，记为嵌合或错误延伸。正确侧翼恢复长度为目标区间之外、仍连续比对到该 truth locus 两侧的碱基数。所有阈值在查看工作流间结果前固定。

## Supplementary Table S1. Complete coral benchmark

| Configuration | Sample–locus recoveries | Mean/median per sample (range) | Distinct panel loci | Shared loci (11/11) | Shared loci (≥9/11) | Median length (bp) | Wall time (min) | Peak RSS (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| GeneMiner2 *k*=23 | 26,528 | 2,411.64/2,384 (2,364–2,546) | 2,909 | 1,653 | 2,250 | 281 | 10.57 | 0.51 |
| TipSeek *k*=23, R0 | 25,101 | 2,281.91/2,262 (2,173–2,464) | 2,881 | 1,125 | 2,050 | 439 | 6.88 | 0.72 |
| TipSeek *k*=23, R1 | 25,933 | 2,357.55/2,341 (2,259–2,500) | 2,885 | 1,471 | 2,134 | 771 | 28.02 | 6.25 |
| TipSeek *k*=23, R2 | 25,933 | 2,357.55/2,341 (2,259–2,500) | 2,885 | 1,471 | 2,134 | 990 | 61.90 | 7.14 |
| TipSeek *k*=31, R1 | 22,085 | 2,007.73/2,008 (1,871–2,112) | 2,787 | 833 | 1,699 | 752 | 19.84 | 2.83 |
| SPAdes + PHYLUCE genome-harvesting | 24,231 | 2,202.82/2,224 (1,985–2,330) | 2,718 | 1,084 | 1,984 | 2,072 | — | — |

R0, R1 and R2 denote zero, one and two rescue rounds, respectively. PHYLUCE resource values are not reported because the complete workflow lacks a same-run, same-core timing and memory record.

## Supplementary Table S2. Coral WGS sample metadata

| Run accession | Taxon | BioSample | BioProject | Library layout | Read length (bp) | Read pairs | Total bases | Analysis status |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | --- |
| CRR2698935 | X | X | PRJCA057506 | X | X | X | X | Included |
| CRR2698936 | X | X | PRJCA057506 | X | X | X | X | Included |
| CRR2698937 | X | X | PRJCA057506 | X | X | X | X | Excluded: X |
| CRR2698938 | X | X | PRJCA057506 | X | X | X | X | Included |
| CRR2698939 | X | X | PRJCA057506 | X | X | X | X | Included |
| CRR2698940 | X | X | PRJCA057506 | X | X | X | X | Included |
| CRR2698941 | X | X | PRJCA057506 | X | X | X | X | Included |
| CRR2698942 | X | X | PRJCA057506 | X | X | X | X | Included |
| CRR2698943 | X | X | PRJCA057506 | X | X | X | X | Included |
| CRR2698944 | X | X | PRJCA057506 | X | X | X | X | Included |
| CRR2698945 | X | X | PRJCA057506 | X | X | X | X | Included |
| CRR2698946 | X | X | PRJCA057506 | X | X | X | X | Included |

## Supplementary Table S3. Coral ground-truth accuracy benchmark

| Workflow | Depth | Scored candidates | Correct locus assignment (%) | Target identity (%) | Chimera/incorrect extension (%) | Correct flank length, median (bp) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GeneMiner2 *k*=23 | Low (X) | X | X | X | X | X |
| TipSeek *k*=23, R1 | Low (X) | X | X | X | X | X |
| SPAdes + PHYLUCE | Low (X) | X | X | X | X | X |
| GeneMiner2 *k*=23 | Medium (X) | X | X | X | X | X |
| TipSeek *k*=23, R1 | Medium (X) | X | X | X | X | X |
| SPAdes + PHYLUCE | Medium (X) | X | X | X | X | X |
| GeneMiner2 *k*=23 | High (X) | X | X | X | X | X |
| TipSeek *k*=23, R1 | High (X) | X | X | X | X | X |
| SPAdes + PHYLUCE | High (X) | X | X | X | X | X |
