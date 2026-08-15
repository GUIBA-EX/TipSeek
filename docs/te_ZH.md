# TE / repeatome

`tipseek te` 面向 genome-skimming 与 WGS 短读长，输出可复现的重复证据与相对丰度。它不宣称完成 TE 注释、插入位点检测，也不能从 capture 数据推断全基因组拷贝数。

## 从这里开始

```bash
cli/tipseek te -f te_samples.tsv -o te_out -p 16
```

```text
taxon_id  sample_id  read1  read2
Taxon_A   A01        A01_R1.fq.gz  A01_R2.fq.gz
Taxon_B   B01        B01.fq.gz
```

`read2` 可省略。只有拥有已分类、最好是近缘物种的库时才使用 `--te-library library.fa`；FASTA 标题采用 `name#Class/Subclass`。

## 一个流程，两条证据通道

```text
reads → discover → curate → annotate → quantify
                 └→ interspersed（可选恢复通道）
```

主通道识别精确重复单元（EQ），保留 paired-read 证据，组装短而受支持的 fragment，并用全部 eligible reads 定量；它刻意保持保守。

`interspersed` 通道不要求 reads 唯一归属 EQ，而使用共享候选 reads 建立稀疏 minimizer-overlap component，再对每个 component 联合组装。需要恢复非串联重复 consensus 时运行：

```bash
cli/tipseek te -f te_samples.tsv -o te_out -p 16 --te-stage interspersed
```

## 将类别读作证据，而不是结论

| 输出类别 | 含义 |
| --- | --- |
| `simple_repeat` | 短周期 motif |
| `tandem_repeat_candidate` / `satellite_candidate` | 重复阵列；短读长不能证明染色体位置 |
| `foldback_like_DNA` | 长且有 reads 支持的倒置重复候选；不等于已证明可转座 |
| `interspersed_repeat_candidate` | 非周期 component，仍需结构或同源证据 |
| `unknown_repeat` | 证据不足；应保留而非过度分类 |

没有 Dfam 或 library 命中，**不等于**不是 repeat。对非模式刺胞动物，这通常只是尚无足以支持家族名称的证据。

## 最值得看的产物

```text
03_annotate/annotation_evidence.tsv    类别、支持度、period、倒置重复分数
03_annotate/repeat_families.tsv        保守相似性分组；EQ 本身不合并
03_interspersed/clusters.tsv           overlap component 与联合组装结构
03_interspersed/consensus.fasta        可供外部注释的 consensus 候选
04_quantify/repeat_signal.tsv          每样本 EQ 丰度与覆盖
04_quantify/repeat_landscape.tsv       reads 到 consensus 的 divergence proxy
05_compare/repeat_superfamilies.tsv    shared、taxon-shared、sample-specific family
```

`signal_rpm` 是相对 read signal。`estimated_genome_fraction` 只是 read-fraction proxy；只有可比较的随机 WGS 文库才能近似解释为 genome fraction，UCE off-target reads 绝不可如此解释。

## 实际解释原则

Tandem 与 foldback 调用首先是结构观察。用 `consensus.fasta` 再接 Dfam、蛋白域 HMM 或近缘珊瑚自建库。只有具有末端结构或蛋白域证据时，才命名为 LINE、LTR 或 DNA transposon。稳定保留 `unknown` 是正确结果，不是失败。
