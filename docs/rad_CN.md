# RAD：用 WGS 补充已有 ipyrad loci

`rad` 用 WGS reads 补充 ipyrad `.loci` 中缺失的新样本；它不从缺失的 RAD reads 伪造等位基因。

| 起点 | 产物 | 边界 |
| --- | --- | --- |
| 已完成的 `.loci` 或已拆样本的 paired RAD reads + 新 WGS reads | 独立 R1/R2 arm、恢复状态、验证后的严格矩阵 | 不推断 R1/R2 中间插入区；不直接证明 allele dropout |

## 推荐路线

```bash
# 从已完成 ipyrad .loci 建立可复用 probe
cli/tipseek rad-probe --ipyrad-loci assembly.loci -o rad_probe

# 为新 WGS 样本恢复独立的两条 arm
cli/tipseek rad --rad-probe rad_probe/rad_reference \
  -f wgs_samples.tsv -o rad_out -p 8

# 写入独立验证后的 strict 矩阵
cli/tipseek rad-validate --rad-probe rad_probe/rad_reference \
  --rad-recovery rad_out/rad_recovery -o rad_validate_out
```

## 输入与 probe 构建

优先使用 `--ipyrad-loci FILE` 读取已完成的 `.loci`。也可用 `--ipyrad-params FILE` 调用 ipyrad（默认步骤 `1--7`），或用 `--rad-denovo` 从已拆样本的 paired RAD reads 建立候选 probe。

`rad-probe` 不负责拆样本、识别酶切位点或去接头。原生 de novo 模式只建立保守 probe，不替代 ipyrad 的完整聚类和矩阵构建。

```text
paired_rad_samples.tsv
sample_id<TAB>R1.fastq.gz<TAB>R2.fastq.gz
```

de novo 模式用 paired-arm 多 seed 找候选，再以 R1、R2 全长距离确认 locus；代表 bait 始终选择真实观测 read。默认要求 `k=31`、深度 `3`、至少 `2` 个样本、arm 长度 `60 bp`、每条 arm 最大编辑距离 `3`。仅在已知 R2 从第二酶切端开始时使用 `--rad-overhang-r2`。

## 恢复与验证

`rad` 只接受 probe 中不存在的新 WGS 样本。样本重复、R1/R2 不配对或 arm FASTA 损坏时会直接停止。

默认流程只使用 k31 招募，速度最快。先按默认参数运行；恢复不足时再依次尝试：

1. `--rad-linked-recruitment`：把一个 arm 命中的 paired fragment 限量提供给同 locus 的另一 arm；默认上限为 256，可用 `--rad-link-max-fragments` 调整。
2. `--rad-fallback-kmers 25`：仅对 k31 未命中的 fragment 使用短 k-mer。它更慢，也更容易招募非特异 reads。

`rad-validate` 默认要求目标覆盖度至少 `0.80`、identity 至少 `0.90`，且本 locus 得分比最佳外源 locus 高至少 `5%`。只有 R1、R2 都通过的样本才进入 strict matrix。完整 contig 保留在 `rad_recovery/`，矩阵中只写目标区间。

## 输出与解释

- `rad_reference/arms/`：每 locus 的 R1/R2 多等位 bait FASTA；
- `rad_recovery/`：每个 WGS sample 的招募、refilter 与组装结果；
- `rad_matrix/rad_sample_locus.tsv`：sample × locus 的 arm 与联合恢复状态；
- `rad_matrix/recovered_arms/`：保留有证据单 arm 的未比对 FASTA；
- `rad_matrix/paired_arms/`：两个 arm 都恢复的 WGS sample，尚未 validate；
- `rad_validated/rad_validation.tsv`：逐 sample × locus × arm 的验证指标和状态；
- `rad_validated/strict_arms/`：原始 RAD bait 加上两个 arm 都验证通过的 WGS sample。

R1/R2 永远是独立观察，流程不会桥接中间未测序区。所有 arm FASTA 都未比对；应在确定缺失数据策略后再进行 MSA。`rad-validate` 不重新招募或组装 reads，也不会修改原始 `rad_recovery/` 或 `rad_matrix/`。

`rad_missing_wgs_recovered` 只表示输入 RAD 矩阵缺失的 locus 已由 WGS 恢复；它本身不是 restriction-site allele dropout 的直接证据。若需作该解释，还应有酶切位点和跨 locus WGS 证据。
