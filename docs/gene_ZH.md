# Gene 与 exon 工作流

默认 `gene` 模式面向多物种 bait 定义的核基因家族，恢复并保留样本内候选 contig；`exon`、`gene-resolve` 与 `gene-tree` 再增加结构注释和系统发育解析。候选数是组装观察，**不是**等位基因或真实拷贝数结论。

## 快速开始

每个 `family_reference/*.fasta` 是一个 family，可含多个物种的 bait；`family_proteins/` 为同名 family 的蛋白 FASTA。

```bash
# 仅恢复候选（默认 gene 模式）
cli/tipseek -f samples.tsv -r family_reference -o gene_output -p 8

# 恢复候选后继续进行蛋白引导的 exon/intron 注释
cli/tipseek --assembly-mode exon \
  -f samples.tsv -r family_reference \
  --gene-protein-reference family_proteins -o exon_output -p 8

# 对齐、建树、选择 strict 一对一子树
cli/tipseek gene-resolve --gene-input exon_output/exon -o gene_resolved -p 8

# strict pseudo-SCO 物种树
cli/tipseek gene-tree --gene-input gene_resolved -o species_strict -p 8 \
  --gene-species-mode strict --gene-aster astral

# 多拷贝家族物种树
cli/tipseek gene-tree --gene-input gene_resolved -o species_multi -p 8 \
  --gene-species-mode multicopy --gene-aster astral
```

默认 `gene` 模式需要 `-f/-r/-o`。`--assembly-mode exon` 还需要 `--gene-protein-reference` 和 miniprot。resolve 需要 MAFFT 与 IQ-TREE；物种树需要 ASTER2 `astral`。exon 模式是默认工作流的扩展，不属于 UCE 模式。

## 模式与后处理命令

| 调用 | 输入 | 作用 | 主要输出 |
| --- | --- | --- | --- |
| `tipseek` | reads + family bait | 招募、refilter、组装及候选汇总 | `<output>/gene/` |
| `tipseek --assembly-mode exon` | reads + family bait + 蛋白参考 | 默认 gene recovery 后继续进行 miniprot ≥0.18 结构注释与模型 QC | `<output>/{gene,exon}/` |
| `gene-resolve` | `<exon-output>/exon/` | 蛋白 MSA、密码子回译、gene tree、无根树上一对一子树选择 | `gene_resolved/` |
| `gene-tree` | `gene_resolved/` | ASTER2 strict 或 multicopy 物种树 | 物种树与 provenance |

## 注释与结构 QC

exon 模式对每条候选 contig 独立运行 miniprot，并保存 `--gff` 产生的 GFF3 与嵌入 PAF。这样可避免一条 partial contig 在 miniprot 的靶序列竞争中压制另一条互补 partial contig。PAF 提供蛋白覆盖度、alignment score、protein CIGAR、frameshift 与 in-frame stop；GFF3 提供 CDS 坐标和 phase。候选中的 `N` 和其他 IUPAC 碱基会保留，因此坐标仍对应原始 contig。

同一 contig 上的模型仅在直接达到重叠阈值时竞争，并按结构缺陷、蛋白覆盖度、alignment score、identity、剪接位点和 miniprot rank 依次择优；被淘汰的桥接模型不会把两个原本独立的位点合并为一个竞争组。主要状态包括 `complete`、`terminal_partial`、`low_coverage`、`frameshifted`、`internal_stop`、`phase_inconsistent`、`ambiguous_cds` 与 `ambiguous_model`。

默认只有 `complete` 和 `terminal_partial` 且无 frameshift、无内部终止、具有明确且连续的 phase、CDS 长度为 3 的倍数、翻译不含 `X`、没有未识别非经典剪接、没有竞争歧义的模型进入 `cds/` 和后续 resolve。注释 CDS 保留一个正常的末端终止密码子；仅在准备蛋白与密码子比对时，才将该密码子及其翻译出的 `*` 成对移除。其他模型仍保存在 `unresolved/` 和结构化 manifest。

如果恰好有一对已选中的 `terminal_partial` 模型互补覆盖同一蛋白的两端，exon 模式会按模型方向排列两条完整候选 contig，默认在中间加入 100 个 `N`，再仅对该派生序列运行一次 miniprot。只有重注释得到唯一、完整、可进入 resolve 的正向模型，而且整段补 N 区完全位于预测内含子内、与外显子无重叠时才接受。因而 CDS 与蛋白不含人工碱基；补 N 的区段只保留在 `supercontigs/`、内含子输出、GFF3 和审计表中。两条来源 partial 仅作审计保留，标记为 `superseded_by_padded_join`。有多种组合或重注释失败时不拼接。该结构桥接不用于推断真实内含子长度。

主要 exon 参数：

- `--gene-max-intron`：最大内含子长度，默认 50,000 bp。
- `--gene-min-model-coverage`：partial model 最低蛋白覆盖度，默认 0.20。
- `--gene-complete-coverage`：完整模型最低覆盖度；仍需覆盖蛋白两端，默认 0.80。
- `--gene-fragment-padding`：上述唯一两片段验证所插入的 `N` 数，默认 100；设为 0 可关闭拼接。
- `--gene-flank`：`genes_flanked/` 两侧实际序列长度，默认 0；不会补 N。

## Resolve 与 QC

`gene-resolve` 默认仅做快速 ML tree。它在树推断前后执行两轮保守 QC：pre-alignment 仅保留可翻译且长度不少于 `--gene-min-aa-length`（默认 30 aa）的候选，并按**不同样本数**检查 `--gene-min-taxa`；post-alignment 在 MAFFT/TAPER 与密码子回译后再次检查样本占有率，并要求至少 `--gene-min-effective-codon-sites`（默认 30）个有效 codon 位点。`--gene-ufboot` 只能为 `0`（默认）或 `≥1000`；后者才在 `tree_selection_qc.tsv` 中提供可用的 branch support。可选参数：

```bash
cli/tipseek gene-resolve --gene-input exon_output/exon -o gene_resolved -p 8 \
  --gene-outgroup outgroups.tsv \
  --gene-taper /path/to/correction_multi.jl --gene-julia julia \
  --gene-ufboot 1000
```

- `--gene-outgroup`：TSV/CSV 第一列为 outgroup sample ID；要求其在 gene tree 上单系。
- `--gene-taper`：在蛋白 MSA 后运行 TAPER；异常、重复或缺失 header 的输出会被拒绝并写入 unresolved。
- `occupancy_qc.tsv`：每个 family 的 pre/post 保留候选数、不同样本数、中位长度、阈值和拒绝原因；同一样本的多个候选只计一次占有率。
- `family_qc.tsv`：仅为通过 post-alignment 检查的对齐 QC（`alignment_pass`），不是整条 resolve 成功标志。
- `tree_selection_qc.tsv`：每个 strict 子树的候选占有率、多候选样本数与 branch support。
- `resolve_manifest.tsv`：每个 family 的最终 resolved/unresolved 原因。

## 输出与解释

```text
gene_output/gene/               # 默认 gene 模式输出
├── family_summary.tsv
├── family_count_matrix.tsv
├── pseudo_sco/
└── multiple_candidate_families/

exon_output/exon/               # --assembly-mode exon 的结构注释输出
├── cds/                         # 仅 resolve-eligible CDS；验证拼接的外显子中无 N
├── proteins/                    # 与 cds/ 一一对应
├── exons/                       # 每个 protein-supported CDS exon 一条记录
├── introns/                     # 每个 intron 一条记录
├── genes/                       # 仅单 contig、含实测 intron 的基因区间
├── genes_flanked/
├── supercontigs/                # 已选模型；验证通过的双 contig 可含内含子 N
├── gff3/                        # 标准化 gene/mRNA/exon/CDS/intron
├── raw_miniprot/                # 每条候选的 GFF3+PAF/stderr，含补 N 验证运行
├── models/{gene_models.tsv,gene_segments.tsv}
├── manifest/
└── unresolved/

gene_resolved/
├── resolved_1to1/                 # 每个 strict 子树的 CDS 与审计 tree
├── unresolved_multicandidate/     # 多拷贝、冲突或失败 family
├── astral_input/resolved_1to1.trees
├── astralpro_input/{multicopy.trees,leaf_to_species.tsv}
├── occupancy_qc.tsv
├── family_qc.tsv
├── tree_selection_qc.tsv
└── resolve_manifest.tsv
```

strict 路线将每个已选子树规范为**每样本一条叶**并交给 ASTER2。multicopy 路线保留完整 gene tree，同时用 `leaf_to_species.tsv` 将候选叶映射到样本。`gene-tree` 在输出目录写入 `gene_tree_provenance.tsv`，记录命令、输入和 SHA-256。
