# Gene 子命令

`gene` 面向多物种 bait 定义的核基因家族。它恢复并保留样本内候选 contig，再以蛋白注释和 gene tree 将可靠的一对一子树与多拷贝/歧义家族分开。候选数是组装观察，**不是**等位基因或真实拷贝数结论。

| 从什么开始 | 得到什么 | 主要边界 |
| --- | --- | --- |
| reads + 多物种 family bait | 候选 contig、严格一对一子树或多拷贝树输入 | 不直接声明真实拷贝数 |

## 快速开始

每个 `family_reference/*.fasta` 是一个 family，可含多个物种的 bait；`family_proteins/` 为同名 family 的蛋白 FASTA。

```bash
# 1. 从 reads 恢复候选；固定使用 original-rust
cli/tipseek gene -f samples.tsv -r family_reference -o gene_output -p 8

# 2. 蛋白引导注释
cli/tipseek gene-annotate --gene-input gene_output/gene \
  --gene-protein-reference family_proteins -o gene_annotation -p 8

# 3. 对齐、建树、选择 strict 一对一子树
cli/tipseek gene-resolve --gene-input gene_annotation -o gene_resolved -p 8

# 4a. strict pseudo-SCO 物种树
cli/tipseek gene-tree --gene-input gene_resolved -o species_strict -p 8 \
  --gene-species-mode strict --gene-aster astral

# 4b. 多拷贝家族物种树
cli/tipseek gene-tree --gene-input gene_resolved -o species_multi -p 8 \
  --gene-species-mode multicopy --gene-aster astral
```

`gene` 需要 `-f/-r/-o`；其余三个子命令只需要 `--gene-input` 和 `-o`。注释需要 miniprot；resolve 需要 MAFFT 与 IQ-TREE；物种树需要 ASTER2 `astral`。

## 子命令

| 子命令 | 输入 | 作用 | 主要输出 |
| --- | --- | --- | --- |
| `gene` | reads + family bait | 招募、refilter、original-rust 组装及候选汇总 | `gene/` |
| `gene-annotate` | `gene/` + 蛋白参考 | miniprot 提取 CDS、exon、intron、supercontig | `gene_annotation/` |
| `gene-resolve` | `gene_annotation/` | 蛋白 MSA、密码子回译、gene tree、无根树上一对一子树选择 | `gene_resolved/` |
| `gene-tree` | `gene_resolved/` | ASTER2 strict 或 multicopy 物种树 | 物种树与 provenance |

## Resolve 与 QC

`gene-resolve` 默认仅做快速 ML tree。它在树推断前后执行两轮保守 QC：pre-alignment 仅保留可翻译且长度不少于 `--gene-min-aa-length`（默认 30 aa）的候选，并按**不同样本数**检查 `--gene-min-taxa`；post-alignment 在 MAFFT/TAPER 与密码子回译后再次检查样本占有率，并要求至少 `--gene-min-effective-codon-sites`（默认 30）个有效 codon 位点。`--gene-ufboot` 只能为 `0`（默认）或 `≥1000`；后者才在 `tree_selection_qc.tsv` 中提供可用的 branch support。可选参数：

```bash
cli/tipseek gene-resolve --gene-input gene_annotation -o gene_resolved -p 8 \
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
gene_output/gene/
├── family_summary.tsv
├── family_count_matrix.tsv
├── pseudo_sco/
└── multiple_candidate_families/

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
