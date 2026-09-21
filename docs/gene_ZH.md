# Gene family 恢复与系统发育解析

[English](gene_EN.md) · [Exon 注释](exon_ZH.md) · [命令行指南](../manual/ZH_CN/command_line.md)

默认 `gene` 模式根据多物种核酸 bait 恢复核基因家族候选 contig。候选数是组装证据，**不是**等位基因数或真实生物学拷贝数估计。

## 候选恢复

`family_reference/` 中每个 `.fa` 或 `.fasta` 文件定义一个 family，可包含多个物种。

```bash
cli/tipseek \
  -f samples.tsv \
  -r family_reference \
  -o gene_output \
  -p 8
```

TipSeek 对每个样本和 family 完成 reads 招募、refilter 与候选组装，并将跨样本汇总写入 `gene_output/gene/`。该模式不分配 exon 坐标，也不直接判断单拷贝状态。

如需蛋白引导的 CDS、exon、intron 与 supercontig，应使用 `--assembly-mode exon`；详见 [exon 注释说明](exon_ZH.md)。`gene-resolve` 消费结构化的 `exon/` 输出，而不是未注释的 `gene/` 目录。

## 工作流边界

| 调用 | 输入 | 作用 | 主要输出 |
| --- | --- | --- | --- |
| `tipseek` | reads + 核酸 family bait | 候选恢复与跨样本汇总 | `<output>/gene/` |
| `tipseek --assembly-mode exon` | reads + coding 核酸 family bait；蛋白参考可选 | 候选恢复与结构注释 | `<output>/{gene,exon}/` |
| `tipseek gene-resolve` | `<exon-output>/exon/` | MSA、密码子回译、gene tree 与一对一子树选择 | `gene_resolved/` |
| `tipseek gene-tree` | `gene_resolved/` | strict 或 multicopy ASTER2 物种树 | 物种树与 provenance |

## Resolve 与 QC

`gene-resolve` 只接收 exon 注释标记为 resolve-eligible 的模型，随后再执行两轮 QC：

- pre-alignment：仅保留可翻译且长度不少于 `--gene-min-aa-length`（默认 30 aa）的候选；`--gene-min-taxa` 按**不同样本数**计算；
- post-alignment：在 MAFFT、可选 TAPER 与密码子回译后再次检查占有率，并要求至少 `--gene-min-effective-codon-sites` 个有效 codon 位点（默认 30）。

默认构建不含 branch support 的快速 ML gene tree。`--gene-ufboot` 只能为 `0` 或至少 `1000`；只有后者能为 `tree_selection_qc.tsv` 提供 branch support。

```bash
cli/tipseek gene-resolve \
  --gene-input exon_output/exon \
  -o gene_resolved \
  -p 8 \
  --gene-outgroup outgroups.tsv \
  --gene-taper /path/to/correction_multi.jl \
  --gene-julia julia \
  --gene-ufboot 1000
```

- `--gene-outgroup`：TSV/CSV 第一列列出 outgroup sample ID；要求其在 gene tree 上单系。
- `--gene-taper`：蛋白 MSA 后的可选 masking；异常、重复或缺失 header 会将该 family 拒绝到 unresolved。
- `occupancy_qc.tsv`：每个 family 的保留候选数、不同样本占有率、中位长度、阈值与拒绝原因。
- `family_qc.tsv`：通过 post-alignment 检查的 alignment QC，不是整体成功标志。
- `tree_selection_qc.tsv`：每个 strict clade 的占有率、多候选样本数与 branch support。
- `resolve_manifest.tsv`：每个 family 最终的 resolved 或 unresolved 状态。

## 物种树

```bash
# strict pseudo-SCO：每个样本一条已选叶
cli/tipseek gene-tree --gene-input gene_resolved -o species_strict -p 8 \
  --gene-species-mode strict --gene-aster astral

# multicopy：完整 gene tree 加叶到样本映射
cli/tipseek gene-tree --gene-input gene_resolved -o species_multi -p 8 \
  --gene-species-mode multicopy --gene-aster astral
```

两条路线均需要 ASTER2 `astral`。strict 路线将每个已选子树规范为每样本一条叶；multicopy 路线保留完整 gene tree，并提供 `leaf_to_species.tsv`。`gene_tree_provenance.tsv` 记录命令、输入路径与 SHA-256。

## 主要输出

```text
gene_output/gene/
├── family_summary.tsv
├── family_count_matrix.tsv
├── pseudo_sco/
└── multiple_candidate_families/

gene_resolved/
├── resolved_1to1/                 # 每个 strict clade 的 CDS 与审计 tree
├── unresolved_multicandidate/     # multicopy、冲突或失败 family
├── astral_input/resolved_1to1.trees
├── astralpro_input/{multicopy.trees,leaf_to_species.tsv}
├── occupancy_qc.tsv
├── family_qc.tsv
├── tree_selection_qc.tsv
└── resolve_manifest.tsv
```
