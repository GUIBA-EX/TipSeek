# 线粒体工作流：常规环形动物线粒体

[English](mitochondria_EN.md) · [项目总览](../README.md) · [命令行指南](../manual/ZH_CN/command_line.md) · [输出说明](../manual/ZH_CN/output.md)

`mito` 使用现有 MainFilter、refilter 与 Rust UCE assembler，从带注释的 GenBank 参考招募 reads，并只在样本 reads 支持时报告闭环。

| 从什么开始 | 得到什么 | 主要边界 |
| --- | --- | --- |
| 带注释 GenBank + 样本 reads | 有 read 支持的闭环或 partial 线粒体序列 | 仅适用于常规单条环形动物线粒体 |

## 适用范围

本流程只面向**常规的单条环形动物线粒体基因组**：基因顺序与参考大体可比，不存在明显多分子结构、大片段重排或需要分相报告的强异质性。它不适合植物或真菌复杂线粒体、分裂或多染色体线粒体、严重重排，或需要专门处理异质性与 NUMT 的研究。

## 最小命令

```bash
cli/tipseek mito \
  -f samples.tsv \
  -o mito_output \
  -r mitochondrial_reference.gb
```

`--mito-max-reads 320` 限制自适应流程的最大输入；`--no-mito-adaptive-stop` 改为一次性过滤。
`--cleanup-intermediates` 仅在最终闭环确认后删除 `filtered/mitochondrion.fq` 与本次样本的自适应 stage 副本，并在输出根目录留下 `cleanup_manifest.tsv`；默认不删除。配合 `--cleanup-dry-run` 时仅生成 `cleanup_preview.tsv`，不删除文件。成功后也清理 `filtered_pe/` 候选 reads；失败或部分组装保留中间文件。

## 流程

```text
GenBank gene/rRNA/tRNA + genome/tile baits
→ MainFilter 招募 paired reads
→ 单一 mitochondrial read pool
→ refilter → Rust UCE 联合图组装
→ finalize 并判定闭环
→ 未成环时：contigs_all 作为样本特异性 seed，再招募一次并重组装（已成环则跳过）
→ 严格 overlap、GFA 唯一路径与 mate-link 连接
→ 原 k 失败时，局部尝试 k−10 / k / k+10 的唯一图路径
→ junction-spanning reads 验证闭环
```

所有 bait 写入一个线粒体 locus。参考仅用于招募与 seed；最终序列不按参考坐标拼接或补洞。mate link 只决定候选邻接和方向，断点碱基必须由同一 filtered read pool 的唯一路径恢复；不能恢复时保持断裂，不插入 `N`。多 k 仅在原 k 不存在可接受路径时运行，且不能绕过唯一性或 junction-read 条件。

`mito` 默认启用无限延伸与 GFA 图输出；若参考较远或覆盖较低，建议显式使用较敏感的 `-kf 17–25 -s 1`。每个自适应深度仍受 `--mito-max-reads` 限制；需要扫描完整文库时应将其设为不小于输入量。

## 成功标准

闭环结果必须同时满足：单一 component、无 `N`、每个连接都来自严格 overlap 或唯一 GFA/read-graph bridge、首尾闭合，并有至少 `--mito-min-junction-support` 条 reads 跨越最终 junction。junction 支持取沿接缝平铺的一小段 k-mer 各位置跨越深度的**最小值**——单个侥幸命中的 k-mer 不能替代一致的跨越覆盖；每条 read 允许在任一链上 ≤1 个错配，一个测序错误不再丢弃真实的跨越 read；若某接缝 k-mer 在 contig 内部也出现，则不计为闭环证据。条件不足时保留 partial 输出，但命令返回失败。

## 参数尽量使用默认值

常规运行只需 `-f samples.tsv -r reference.gb -o output`。线程数自动检测，reads 深度自动递增，配对保留、补救回退和闭环验证自动执行，无需额外开关。内部拼接阈值沿用已验证的默认值。

仅在需要时使用通用选项 `-p`、`--cleanup-intermediates` 和 `--cleanup-dry-run`。旧的 `--mito-genbank` 及 mito 专家参数仍兼容既有脚本，但不需要在新命令中显式填写；`-r` 与旧参考参数不能指定不同文件。

## 输出

- `<sample>/mito/mitochondrial_assembly.fasta`：闭环序列或 partial components（经审计的原始组装，从不旋转）。
- `<sample>/mito/mitochondrial_standardized.fasta`：仅对已验证闭环写出，旋转到可复现的基因起点（优先 tRNA-Phe，否则用保守基因兜底）并置于该基因的编码链，便于跨样本比较。只对已有碱基重排或逆相补——不添加任何参考碱基——FASTA 头记录 anchor、strand、偏移与错配数。无法可靠定位锚定基因时不写出。
- `<sample>/mito/mitochondrial_assembly_summary.tsv`：兼容的粗粒度状态、`resolution_reason`、连接、junction 支持,以及 `ambiguous_bases` / `ambiguous_per_kb`（primary contig 每 1000bp 的 N 数）。
- `<sample>/mito/mitochondrial_evidence.json`：结构、图、mate-link 与 junction 的机器可读证据。
- `<sample>/mito/mitochondrial_feature_evidence.tsv`：参考 feature 的 canonical 21-mer 相似度；仅说明与参考的精确锚点共享，不能据此判断远缘样本的基因存在或缺失；`translation_status=not_checked`，不是 CDS 注释或翻译判定。
- `<sample>/mito/mitochondrial_mate_links.tsv`：已接受的 read-supported links。
- `.gm2_mito_reference/metadata/mitochondrial_genes.tsv`：bait 元数据；坐标为 0-based 半开区间，`segments_0_half_open` 保留跨 origin 或 `join(...)` feature 的全部片段。

## 重跑与补救结果

二次过滤在任一端通过时保留双端 reads，保证后续图连接与闭环验证接收正确的配对。补救组装只有通过 reads 支持的闭环验证后才替换首轮部分结果；为空、仍为部分组装或运行失败时恢复首轮结果。

每次最终判定都会替换本模块管理的输出：无 contig 时摘要和 JSON 均记录 `no_contigs`，非预期错误记录 `failed`，旧标准化序列和特征证据会清除。要求闭环但未成环时保留本次部分结果，同时返回错误。单次模式将 `--max-reads` 传给招募阶段，单位为每个输入文件 1,048,576 条记录，`0` 表示全部 reads。
