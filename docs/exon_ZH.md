# Exon/intron 注释

[English](exon_EN.md) · [Gene 工作流](gene_ZH.md) · [命令行指南](../manual/ZH_CN/command_line.md)

`--assembly-mode exon` 在默认 gene family 恢复流程之后增加蛋白引导的 exon/intron 注释。它适用于具有同名蛋白参考的核基因家族，不用于 UCE。结果包括保留坐标的基因模型、可进入 resolve 的 CDS 与蛋白，以及被拒绝或存在歧义模型的完整审计记录。

## 输入与命令

核酸参考目录与蛋白参考目录必须使用相同的 family stem。核酸 bait 可包含多个物种；对应蛋白文件必须使用 `.faa` 扩展名。

```text
family_reference/
├── geneA.fasta
└── geneB.fa

family_proteins/
├── geneA.faa
└── geneB.faa
```

exon 模式是一条完整工作流；不要再添加 `filter`、`assemble` 等阶段子命令。

```bash
cli/tipseek --assembly-mode exon \
  -f samples.tsv \
  -r family_reference \
  --gene-protein-reference family_proteins \
  -o exon_output \
  -p 8
```

必需输入为样本表、核酸 family bait 与同名蛋白参考。注释要求 miniprot 0.18 或更高版本。TipSeek 先将恢复候选写入 `exon_output/gene/`，再将结构结果写入 `exon_output/exon/`。

## 实现方法

1. 默认 gene 工作流完成 reads 招募、refilter、组装及跨样本候选汇总。
2. 对每条候选 contig 独立运行 miniprot，避免一条 partial candidate 通过 target competition 压制另一条互补 candidate。
3. 嵌入 PAF 提供蛋白坐标、覆盖度、score、identity、protein CIGAR、frameshift 与 in-frame stop；GFF3 提供 exon/CDS 坐标、方向和 phase。
4. 重叠模型按结构缺陷、蛋白覆盖度、alignment score、identity、剪接证据与 miniprot rank 确定性竞争；不重叠模型保持独立。
5. 对选中模型分类，写入结构化审计表，并分别进入 resolve-eligible 或 unresolved 输出。

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

只有在恰好存在一对已选中的 `terminal_partial` 模型、二者互补覆盖同一蛋白两端，且该蛋白没有已选完整模型时，才尝试补 N。

TipSeek 按模型方向排列两条来源 contig，在中间插入 `--gene-fragment-padding` 个 N（默认 100），再对派生序列重新运行 miniprot。仅当重注释得到唯一、完整、可进入 resolve 的正向模型，而且整段补 N 区完全位于预测 intron 内且不与任何 exon 重叠时，才接受拼接。

因此：

- 接受的 CDS 与蛋白不含人工碱基；
- 补 N 仅保留在 `supercontigs/`、intron 输出、GFF3 与审计表中；
- 来源 partial 模型转为仅审计状态 `superseded_by_padded_join`；
- 有多种组合或重注释失败时不拼接；
- 补 N 只验证结构桥接，不估计真实 intron 长度。

使用 `--gene-fragment-padding 0` 可关闭此步骤。

## 主要参数

| 参数 | 默认值 | 含义 |
| --- | ---: | --- |
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
    ├── manifest/                 # candidate、ID、warning、fragment 与命令审计
    └── unresolved/
```

模型状态与 eligibility 查看 `models/gene_models.tsv`，坐标与 splice class 查看 `models/gene_segments.tsv`，所有补 N 尝试查看 `manifest/fragment_groups.tsv`。后续运行：

```bash
cli/tipseek gene-resolve --gene-input exon_output/exon -o gene_resolved -p 8
```

resolve QC 以及 strict/multicopy 物种树推断见 [gene 工作流说明](gene_ZH.md)。
