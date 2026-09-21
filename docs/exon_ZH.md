# Exon/intron 注释

[English](exon_EN.md) · [Gene 工作流](gene_ZH.md) · [命令行指南](../manual/ZH_CN/command_line.md)

`--assembly-mode exon` 在默认 gene family 恢复流程之后增加蛋白引导的 exon/intron 注释。默认直接从 coding 核酸 bait 推导蛋白参考，不要求另备 `.faa`；外部蛋白可作为可选覆盖。该模式不用于 UCE。

## 输入与命令

`-r` 目录中每个 `.fa`、`.fas` 或 `.fasta` 文件定义一个 family，可包含多个物种的 CDS、拼接 exon 或其他保持 coding frame 的核酸 bait。

```text
family_reference/
├── geneA.fasta
└── geneB.fa
```

默认流程使用标准核遗传密码表自动翻译这些核酸参考。若 bait 包含 intron、属于非编码序列、使用非标准遗传密码或不能解释为 coding sequence，应提供外部蛋白目录；文件必须以同名 family stem 加 `.faa` 命名。目录可以不完整：有 `.faa` 的 family 使用外部蛋白，其余 family 仍自动翻译。

exon 模式是一条完整工作流；不要再添加 `filter`、`assemble` 等阶段子命令。

```bash
cli/tipseek --assembly-mode exon \
  -f samples.tsv \
  -r family_reference \
  -o exon_output \
  -p 8
```

可选的外部蛋白覆盖为：

```bash
cli/tipseek --assembly-mode exon \
  -f samples.tsv -r family_reference \
  --gene-protein-reference family_proteins \
  -o exon_output -p 8
```

注释要求 miniprot 0.18 或更高版本。TipSeek 先将恢复候选写入 `exon_output/gene/`，再将结构结果写入 `exon_output/exon/`。

## 实现方法

1. 默认 gene 工作流完成 reads 招募、refilter、组装及跨样本候选汇总。
2. 未提供同名 `.faa` 时，TipSeek 平等检查正反两条链的六个 reading frame。允许并移除正常末端终止；任何含内部终止的 frame 都直接排除，不修补也不屏蔽后继续使用。
3. 只有一种无内部终止 frame，或某个 frame 从起始甲硫氨酸和/或末端终止获得唯一最强的 CDS 端点支持时，才直接接受该序列；其中最长的可信翻译作为 family anchor。其余歧义序列只有在某个 frame 对该 anchor 获得唯一且大于零的氨基酸 3-mer Dice 相似度时才接受。没有可信 anchor、最佳结果并列或完全不相似时，不为该序列输出蛋白，并记录 `ambiguous_reading_frame`；若该 family 没有其他可用参考，则需要提供同名 `.faa`。所有候选和最终选择写入 `manifest/reference_translation.tsv`。
4. 对每条候选 contig 独立运行 miniprot，避免一条 partial candidate 通过 target competition 压制另一条互补 candidate。
5. 嵌入 PAF 提供蛋白坐标、覆盖度、score、identity、protein CIGAR、frameshift 与 in-frame stop；GFF3 提供 exon/CDS 坐标、方向和 phase。
6. 重叠模型按结构缺陷、蛋白覆盖度、alignment score、identity、剪接证据与 miniprot rank 确定性竞争；不重叠模型保持独立。
7. 对选中模型分类，写入结构化审计表，并分别进入 resolve-eligible 或 unresolved 输出。

基因组输出保留 `N` 与其他 IUPAC 碱基，因此所有坐标始终对应实际用于注释的候选序列。

## 模型接纳条件

只有已选中的 `complete` 与 `terminal_partial` 模型可以进入 `cds/` 和后续 `gene-resolve`。合格模型必须满足：

- 无 frameshift 或内部终止；
- CDS phase 明确且连续；
- CDS 长度为 3 的倍数；
- 翻译结果不含 `X`；
- 不含未识别的非经典剪接；
- 不存在未解决的模型竞争。

注释 CDS 保留一个正常末端终止密码子；仅在准备蛋白与密码子比对时，才将它与对应末端 `*` 一并移除。其他状态，包括 `low_coverage`、`frameshifted`、`internal_stop`、`phase_inconsistent`、`ambiguous_cds` 和 `ambiguous_model`，保留在 `unresolved/` 与 manifest 中。

## 经验证的补 N

只有在恰好存在一对已选中的 `terminal_partial` 模型、二者互补覆盖同一蛋白两端，且该蛋白没有已选完整模型时，才尝试补 N。候选对还必须满足：重叠不超过较短模型的 10%，合并后新增覆盖不少于蛋白长度的 20%，总覆盖达到 `--gene-complete-coverage`，并在蛋白两端各 2%（至少 3 aa）的容差内到达首尾。

TipSeek 按模型方向排列两条来源 contig，在中间插入 `--gene-fragment-padding` 个 N（默认 100），再对派生序列重新运行 miniprot。仅当重注释得到唯一、完整、可进入 resolve 的正向模型，而且整段补 N 区完全位于预测 intron 内且不与任何 exon 重叠时，才接受拼接。

这一步依赖来源 contig 在断点一侧保留足够的真实 intron 与剪接边界证据。若两个片段在 exon 边界精确截断并完全缺失 intronic splice flank，重注释通常不能把 N 区可靠识别为 intron，因此不会接受拼接。补 N 也不会重建未恢复的末端 exon 或其他 coding sequence。

因此：

- 接受的 CDS 与蛋白不含人工碱基；
- 补 N 仅保留在 `supercontigs/`、intron 输出、GFF3 与审计表中；
- 来源 partial 模型转为仅审计状态 `superseded_by_padded_join`；
- 有多种组合或重注释失败时不拼接；
- 补 N 只验证结构桥接，不估计真实 intron 长度。

在真实海星 *Patiria pectinifera* 基因测试中，3 个含真实 2000-nt intron 的完整候选均在不提供 `.faa` 时恢复了正确 CDS、方向、exon/intron 边界与 phase。保留断点两侧各 120 nt 真实 intronic flank 的互补片段中，满足上述门槛的一对通过 100-N 重注释并得到与真实 CDS 完全一致的结果；另外两组分别因蛋白区间重叠过多和末端片段低于覆盖门槛而保持未拼接。该测试验证的是保守接纳逻辑，而不是保证 reads 组装一定能恢复一对互补末端片段。

使用 `--gene-fragment-padding 0` 可关闭此步骤。

## 主要参数

| 参数 | 默认值 | 含义 |
| --- | ---: | --- |
| `--gene-protein-reference` | 自动推导 | 可选 `.faa` 目录；同名 family 覆盖自动翻译 |
| `--gene-miniprot` | `miniprot` | miniprot 可执行文件 |
| `--gene-max-intron` | `50000` | miniprot 接受的最大 intron 长度 |
| `--gene-min-model-coverage` | `0.20` | partial model 最低蛋白覆盖度 |
| `--gene-complete-coverage` | `0.80` | complete model 最低覆盖度；仍要求覆盖蛋白两端 |
| `--gene-fragment-padding` | `100` | 唯一双片段验证插入的 N 数；`0` 关闭拼接 |
| `--gene-flank` | `0` | `genes_flanked/` 两侧加入的实测碱基数；不会补 N |

## 输出

```text
exon_output/
├── gene/                         # 恢复候选与 family 汇总
└── exon/
    ├── cds/                      # 可进入 resolve 的 CDS；exon 中无补入的 N
    ├── proteins/                 # 与 cds/ 一一对应
    ├── exons/                    # 每个受支持 CDS exon 一条记录
    ├── introns/                  # 每个预测 intron 一条记录
    ├── genes/                    # 实测单 contig gene span
    ├── genes_flanked/
    ├── supercontigs/             # 已选 span；验证拼接可含 intron N
    ├── gff3/                     # 标准化 gene/mRNA/exon/CDS/intron 层级
    ├── raw_miniprot/             # 逐候选结果与补 N 验证运行
    ├── models/
    │   ├── gene_models.tsv
    │   └── gene_segments.tsv
    ├── manifest/
    │   ├── reference_translation.tsv # 自动 frame/strand 选择审计
    │   ├── derived_proteins/     # 实际使用的自动推导蛋白
    │   └── ...                   # candidate、ID、warning、fragment 与命令审计
    └── unresolved/
```

参考来源和 frame 选择查看 `manifest/reference_translation.tsv`；模型状态与 eligibility 查看 `models/gene_models.tsv`，坐标与 splice class 查看 `models/gene_segments.tsv`，所有补 N 尝试查看 `manifest/fragment_groups.tsv`。后续运行：

```bash
cli/tipseek gene-resolve --gene-input exon_output/exon -o gene_resolved -p 8
```

resolve QC 以及 strict/multicopy 物种树推断见 [gene 工作流说明](gene_ZH.md)。
