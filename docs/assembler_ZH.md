# 2. Assembler

[English version](assembler_EN.md)

Assembler 将招募后的逐 locus reads 组装为 contig。应按生物学目标选择流程，而不是按测序平台或文库名称选择。

| 工作流 | 适用目标 | 默认后端 | 结果 |
|---|---|---|---|
| `--assembly-mode original` | exon、SCO、核 marker、线粒体 marker | `original-rust` | 参考引导 contig，默认按参考裁切 |
| `--assembly-mode uce` | genome skimming 或 target capture 中的 UCE | `uce-rust` | UCE core 及有 read 支持的 flank |

## 后端

`original-rust` 是 `original` 的确定性 Rust 兼容后端；为兼容旧命令，`original` 仍是其别名。`uce-rust` 是 `uce` 的唯一后端。

所有后端均使用参考位置明确的 seed、read k-mer 支持、双向延伸和 read-slice 验证。

## Original 工作流

原版算法优先延伸权重最高的边，将其他分支压栈以便回溯。每侧最多保留三条候选，再组合和评分。该行为适用于 exon 与细胞器 marker 等较短的参考引导目标。

```bash
cli/geneminer2 -f samples.tsv -r references -o output -p 8 \
  --assembly-mode original
```

`original-rust` 可用 `--reuse-reference-cache` 复用带版本的二进制参考 k-mer cache。

## UCE reads 路径：`ucefilter` 与 main + re

默认 UCE 路径为：

```text
原始 paired FASTQ → ucefilter → filtered/ → uce-rust → contig
```

`ucefilter` 是 UCE 专用的融合 reads 阶段：它在一次原始 FASTQ 扫描中完成宽松招募、完整 paired fragment 保留、方向和 maximum-exact-match 证据、以及按 locus 的深度/位置选择。因此它在默认 UCE 中替代通用的 `MainFilter → filtered_pe/ → refilter → filtered/`（**main + re**）路径；默认命令应使用 `filter assemble`，不需要也不能单独把 `refilter` 当作后续步骤。

只有 `--legacy-uce-filter` 会恢复 main + re 兼容路径，主要用于历史结果对照或诊断。其 `filtered_pe/` 是宽松招募候选，`filtered/` 才是组装输入；默认 `ucefilter` 直接写 `filtered/`，不产生前者。

UCE 模式在主组装接受 contig 后默认执行一轮 rescue；`--no-uce-rescue-reads` 可关闭，`--uce-rescue-rounds 2` 可显式请求第二轮。rescue 借鉴 MITObim 的 baiting-and-iterative-mapping 核心思路，组装固定使用 k=21：第一轮用原参考与已接受 contig 共同招募，第二轮只对仍增长的 locus 使用两端窗口。旧 contig 必须完整保留；每侧新增区域独立要求至少 30 bp、85% breadth、最大 gap 30 bp、2 个独立 fragment 和 1 个跨旧 core 边界 fragment；不通过即按侧或按 locus 逐一回退。每轮最终接受前还会比较 rescue 前后序列；若本轮新引入至少 150 bp 的精确倒置重复，则以 `reverted_inverted_repeat` 状态恢复该 locus 的本轮备份。`--uce-rescue-inverted-repeat-min-bp 0` 可关闭该检查。实验性 `--uce-rescue-reverse-reuse-reference-scale` 只缩放与任一组装臂已走路径互为反向互补节点的 reference bonus，默认 `1.0`（不缩放），不改变 reads depth。rescue 是对主结果的受约束延伸，不是参考补洞。

## UCE 选择、组装与 QC

UCE 模式默认通过融合的 Rust `ucefilter` 招募和选择 reads：rolling k-mer 粗招募、run-k 方向验证、最大精确匹配和逐 locus 自动选择在同一次原始 reads 扫描中完成，最终直接写入 `filtered/`。自动选择先沿用动态 exact-match 阈值去除弱证据；候选少于 512、估计深度不超过 160×或只覆盖少于 48/64 个参考区间的 locus 不做 reservoir 压缩。饱和 locus 至少保留 60% 合格核心，并按参考位置分散选择；跨 bait/contig 边缘的 reads 按 overhang 长度分层，每个长度保留少量高质量 PE fragment，避免切断 core 到 flank 的重叠阶梯。它不产生 GM2 或按 locus 重复的候选 FASTQ；完整 PE fragment 只保存一次并作为一个单位保留。TStools 沿用 GeneMiner2 原版的资源预算调度：`-p` 由已就绪的 filter、refilter 和 assemble 任务共享；样品完成一个阶段后立即进入下一阶段队列。`-p 1` 时各阶段串行；更大数值时，assemble/refilter 任务获得受限的 2--6 线程份额，每个招募任务预留 1--2 个预算线程。fragment bank 上限为内部自动策略：取 Linux `MemAvailable` 与进程 cgroup 剩余内存的较小值，预留一半给系统及工作流其他缓冲区，再按可能并发的 UCEFilter 样品数均分，并将每样品限制在 4 GiB；无法获得内存信息时回退到 512 MiB。`backbone` 策略在每个气泡作一次有限前瞻决策，提交胜出分支而不回溯。`--uce-rescue-reads` 的细节见前文"UCE reads 路径"一节。

CPU 并行度默认使用 `-p auto`：TStools 统计 affinity/cpuset 允许的物理核心，再应用 cgroup 或作业调度器上限；显式整数可覆盖自动检测。调度器保留原版每个招募任务 1--2 个预算单位的设计，并将该预算与 UCEFilter 当前支持的单个主体计算 worker 分开。

`--uce-alignment-shadow` 是默认关闭的证据模式：在自动选择完成后确定性、均匀地抽取每 locus 最多 64 个 fragment，用最大精确匹配定位参考窗口，再执行局部 affine-gap 比对。原始逐 mate 证据和 per-locus 汇总分别写入 `alignment_shadow.tsv` 与 `alignment_shadow_summary.tsv`。初始 bait、whole-contig rescue 和 terminal rescue 分轮保存；短 bait terminal 只代表设计边界，只有 assembly contig 上的 terminal 才可用于评价延伸。

应检查 `uce_assembly_summary.csv`、`uce_rescue_summary.csv`、`uce_rescue_rounds.csv`、`contigs_all_low/` 和下游比对；逐侧长度、breadth、gap、fragment 与 bridge 证据都写入逐轮审计表。主要护栏包括 unique-read density、支持 breadth、最大无支持缺口、k-mer depth CV 与 max/median depth ratio。`--uce-path-strategy search` 只应用于与旧分支枚举的敏感性比较。

文件字段见[输出说明](../manual/ZH_CN/output.md)，参数见[命令行指南](../manual/ZH_CN/command_line.md)。
