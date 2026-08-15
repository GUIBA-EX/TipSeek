# 3. Profiling

[English version](profiling_EN.md)

Profiling 是从 WGS 或 metagenome reads 中免组装恢复任意扩增子 marker 的参考序列级证据流程。它报告与参考序列的相容证据，不生成 contig，也不是生物量比例。

## 流程

```text
TipSeek k-mer 招募 → Themisto 伪比对 → 参考序列级证据
```

不会运行 `refilter`、`assemble`、`combine` 或 `tree`。`--profile-kmer-size` 同时设定招募与 Themisto 的 k-mer，必须是 15–31 的奇数。

## 输入

始终需要：

- 通过 `-r` 提供一个 `.fa` 或 `.fasta` marker 参考库。
- 含 WGS 或 metagenome reads 的样本表。
- `themisto`，放在 `PATH` 或以 `--profile-themisto` 指定。

可选：

- `--profile-group-map` 两列 TSV：`reference_id<TAB>group`。它是参考序列级结果中的可选注释列。
- `--profile-decoy`：可能的非目标序列 FASTA。

reference ID 是 FASTA 标题第一个空白前字段。若提供 map，它必须恰好覆盖每条参考；重复行仅在 group 相同时允许。

## 运行与 cache

```bash
cli/tipseek profiling \
  -f samples.tsv -r marker_reference.fasta \
  -o output -p 8
```

只有需要给参考序列附加类别时，才加 `--profile-group-map marker_groups.tsv`。Themisto index 使用内容寻址 cache。用 `--profile-index-dir` 在多次运行间共享；只有需要重建时才加 `--profile-force-rebuild`。直接调用 `marker_profile` 时，若改了参考、group map 或 decoy，应换用新 cache 或加 `--force-rebuild`。

## 输出与解释

每个样本写入 `marker_profile/`。

### 主结果：`marker_reference_support.tsv`

每条命中的参考序列一行：

- `hit_queries`：与该参考相容的 query 数。
- `fractional_queries`：共享 query 的分数化支持；一条 query 有 N 个候选时，每个候选得 `1/N`，不会被重复计 N 次。
- `singleton_queries`：只与该参考相容的 query 数。
- `ambiguity_status`：`has_singleton_support` 或 `shared_only`。

它表示“与某参考序列相容”的证据，不代表该参考必然唯一存在，也不等于生物量丰度。

`marker_qc.tsv` 记录伪比对及运行参数。`marker_reference_metadata.tsv` 记录 Themisto color 到参考序列及可选 group 注释的映射。

## QC 与校准

解释结果前，应检查伪比对 query 数、各参考的分数化和单例支持，以及 decoy 参考的证据。应以匹配 marker 与测序深度的阴性对照、混合样和 downsampling 选择参考库与伪比对阈值。跨样本比较时，固定参考库、可选注释 map、k-mer、阈值与 decoy 策略。

参见[Filter](filter_ZH.md)、[输出说明](../manual/ZH_CN/output.md)和[命令行指南](../manual/ZH_CN/command_line.md)。
