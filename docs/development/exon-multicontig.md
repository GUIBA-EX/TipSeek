# 长内含子与多 contig 恢复验证（2026-10-11）

本阶段改善了长内含子背景下短编码片段的保留和参考顺序恢复，但没有证明已解决完整基因重建。HX50 的 12 个困难 SCO 家族仍没有通过完整 CDS 验收的模型。参考引导的拼接单独输出，不能等同于同一拷贝的真实连接。

## 实现与参考逻辑

本地参考了 Captus 的 `captus/bioformats.py` 中 `parse_model`、`calc_prot_len_matched`，以及 DNA 提取分支的 `greedy_assembly_partial_hits`；参考了 HybPiper 的 `hybpiper/exonerate_hits.py` 中 `_remove_subsumed_hits`、`_trim_overlapping_hits`、`_create_stitched_contig`。借鉴其参考坐标排序、重叠去重、缺失编码区补 N 和区分 CDS 拼接与基因组结构的原则，独立实现 Rust 逻辑。没有把 Captus 的 DNA 分支误作蛋白模型，也没有实现 HybPiper 的 reads 回贴嵌合检查。

- exon 模式不再在消耗约一半参考种子后停止组装；保留高、低支持 contig 供后续注释筛选。每轮移除起始种子保证前进。普通 gene 模式维持原行为。
- exon 分类保留至少 30 nt 的候选，不再用参考长度的三分之一过滤短片段。miniprot 以 `-L 10 -n 2 --outc 0` 保留短证据；最低 10 aa、80% identity 和结构质控仍约束串联来源。没有恢复微外显子功能。
- 在每条参考蛋白内先比较模型，再搜索两个及以上 contig 的单调链，避免各片段选择不同物种参考后无法串联。同一 contig 不可重复进入链；重叠必须有一致 CDS，最多 12 aa 且不超过较短片段的三分之一；相邻 identity 差不超过 15 个百分点。
- 竞争路径和搜索截断不自动接纳。每参考最多 128 节点、每终点最多 64 路径；先剪掉可证明冗余的子路径，再施加上限。未来仍可出现的 contig 会参与剪枝判断，避免错误丢弃可扩展路径。十片段和重复 contig 回归测试覆盖此逻辑。
- `fragment_cds/` 保存 `reference_ordered_unvalidated` 假说，内部编码缺口按每 aa 三个 N 填补，不进入自动 resolve。跨参考得到的相同 contig 顺序去重。
- 唯一链满足覆盖及两端条件后，才尝试对定向 contig 间补 N 并重新注释。所有补 N 边界都必须落在内含子内，且每个来源必须保留编码支持；跳过中间片段的比对不能接纳。来源模型只有在验证成功后才被替代。
- `fragment_chains.tsv` 记录全体来源、编码缺口和重叠裁剪。即使结构复核通过，`physical_linkage` 仍为 `unverified`，不能据此声称拷贝同源或实际内含子长度已确认。
- 修复了到达参考两端却跳过内部残基的假完整模型。miniprot 的 query-only `I` 和 Exonerate 的 query-only `G` 从覆盖率扣除，模型标为 `internal_partial`，不能通过完整补 N 验证。这是参考相对差异，不是该物种真实缺失的证据。
- `exonerate-rs` 仍作为独立后端，通过当前环境的可执行程序名调用，不写入绝对路径。

## HX50：从真实 reads 开始

使用用户已有的海星 SCO 数据：`/home/yasuda2/staruce/test/01_fastp/HX50_{1,2}.clean.fq.gz`，样本 HX50（Choriaster granulatus）。这是实际测序清洗 reads，区别于旧测试中的 ART 模拟 reads。

先按 P. pectinifera 注释选出 12 个 SCO 家族：每家族至少两个 ≥20 kb 内含子，最大参考内含子约 25.9–71.1 kb。选择发生在观察结果前。每家族保留原来的多物种参考。运行实际 reads 招募、重过滤、组装、分类和注释；最终二进制重新运行两后端注释。旧版本复用完全相同的已招募 reads，以控制招募差异，运行旧组装、分类与 miniprot 注释。

**限制：P. pectinifera 不是 HX50 同物种、同个体真值。** 下表表示参考相对找回率。不能据此断言 HX50 本身有同样长度的内含子。覆盖率使用通过 identity/结构质控的证据并集，允许短模型和未在跨参考竞争中胜出的证据；它不证明这些片段来自一个拷贝。query 缺口不计覆盖；跨剪接密码子的 miniprot U/V 保守地不计入覆盖。

| 指标 | 旧版 miniprot | 新版 miniprot | 新版 exonerate-rs |
| --- | ---: | ---: | ---: |
| 保留候选 contig | 22 | 215 | 215 |
| ≥80% 覆盖的参考外显子 | 1/185 | 56/185 | 57/185 |
| 参考外显子找回率 | 0.54% | 30.27% | 30.81% |
| 每家族 CDS 证据并集平均覆盖 | 5.17% | 42.81% | 44.50% |
| 完整、可接纳 CDS 模型 | 0 | 0 | 0 |
| 单独输出的未验证 CDS 假说 | 0 | 14 | 12 |
| 注释耗时（秒，单次） | 0.714 | 6.15 | 81.02 |

新版 miniprot 有 117 条歧义链、20 条搜索超限记录；exonerate-rs 分别为 119、24。这些记录按参考蛋白产生，不能解释为独立基因或真实错误连接的数量。

完整 CDS 数是已选中、可进入 resolve 且状态为 complete 的模型数。参考顺序 CDS 假说可能只覆盖局部，不能计为完整 CDS。没有同个体组装或验证连接的真值，因此真实错误合并率记录为 `null`；歧义拒绝数量不是错误合并率。

## 已知来源片段对照

对同一批 12 个 P. pectinifera 基因，从真实参考基因组按 ≥20 kb 内含子截成 3–4 个片段，并保留每侧 120 nt 边缘。这是控制性基因组片段对照，**不是从 reads 组装的真实测序测试**，参考注释本身也不是实验验证真值。

最终 miniprot 输出 9 个参考顺序 CDS 假说，exonerate-rs 输出 5 个及 1 个歧义结果；分别尝试 5、3 个基因组补 N 拼接，均被完整结构要求拒绝。保留这些失败，避免把内部参考缺口、移码或非典型剪接掩饰为成功。三 contig、45/10/45 aa 的独立集成测试验证可接受的完整链；跳过中间片段、冲突重叠、竞争拷贝与重复 contig 有负对照。

## 性能与验证边界

全量 reads 招募日志记录 95.729 秒。旧组装 0.865 秒，本次扩大种子搜索的组装日志记录 0.936 秒；这是单次观察，不能当作稳定性能比。候选从 22 增加到 215，扩大了注释工作量。旧候选集 miniprot 注释约 0.714 秒；新候选集耗时见上表。两后端线程参数不同，且计时期间存在并行验证任务，因此这些时间仅用于记录，不能当作后端公平速度排名。本阶段没有满足“所有情形只能比旧版快”的目标；不能把冗余路径剪枝等同于端到端已加速。

312 项 workspace 测试通过（串行测试线程）、Clippy 全目标零警告、release 构建成功、验证脚本编译检查通过。此前一次并发测试遇到临时脚本执行 `Text file busy`，串行完整重跑通过；未据此修改生产程序或宣称已定位该瞬时竞争。

## 复现与记录

- [真实 reads 指标、逐家族结果、输入信息和二进制/源码摘要](exon-multicontig-hx50-validation.json)
- [控制性片段结果与执行参数](exon-multicontig-fragment-controls.json)
- `scripts/prepare_starfish_sco_reads_validation.py`：建立 panel，保留原多物种参考。
- `scripts/run_starfish_sco_baseline.py`：旧版本共享 reads 基线。
- `scripts/report_starfish_sco_recovery.py`：同口径参考外显子/蛋白残基覆盖统计。
- `scripts/validate_exon_fragment_controls.py`：已知来源多 contig 对照。

原始和中间结果位于 `/tmp/tipseek-hx50-multicontig`、`/tmp/tipseek-sco-fragment-controls-verified`；这些临时目录不属于持久交付，汇总报告已复制到本目录。脚本目前针对这台机器的海星数据路径，复现到其他环境时需替换输入路径。

从 reads 开始的执行参数：

```sh
GM2_COMPONENT_DIR=/tmp/tipseek-hx50-multicontig/components \
  tipseek --assembly-mode exon \
  -f /tmp/tipseek-hx50-multicontig/samples.tsv \
  -r /tmp/tipseek-hx50-multicontig/references \
  -o /tmp/tipseek-hx50-multicontig/after \
  -p 8 -kf 21 -s 1 -sb unlimited --min-depth 2 \
  --gene-miniprot miniprot
```

本次组件目录只用于验证新构建；生产调用无须固定组件路径。最终独立后端注释执行 `gene_workflow annotate --backend exonerate-rs --input .../after/gene --nucleotide-reference .../references --out .../verified-exonerate --threads 1`，默认从 PATH 找到 `exonerate-rs`。完整参数见报告中的计时记录和旧版本记录。

后续主要缺口是 reads 支持的拷贝区分/嵌合检查，以及仍未组装或未通过结构验收的编码区。当前阶段改善了碎片找回和透明报告，没有把这些未解决问题转成自动接受。
