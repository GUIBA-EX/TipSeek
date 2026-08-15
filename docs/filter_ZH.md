# MainFilter：reads 招募

[English](filter_EN.md) · [输出说明](../manual/ZH_CN/output.md) · [命令行指南](../manual/ZH_CN/command_line.md)

MainFilter 是 TipSeek 的通用 reads 招募层。它用参考序列中的 k-mer 扫描原始 reads，把命中 read（或完整 paired-end fragment）分配到相应 locus，供后续 `refilter` 和组装使用。它不是比对器，也不直接判定 contig 正确性、marker 丰度或群体结构。默认 UCE 不使用本页的两步路径，而使用融合的 `ucefilter`；见 [Assembler](assembler_ZH.md)。

## 在工作流中的位置

```text
FASTQ/FASTA + 每 locus 一个参考 FASTA
              ↓
          MainFilter
              ↓
      filtered_pe/（第一次、宽松招募）
              ↓
          refilter（去歧义、深度/大小限制）
              ↓
          filtered/ → assembler
```

在 `original`、gene 与 `--legacy-uce-filter` 路径中，只要一对 reads 中任一 mate 命中，该 paired fragment 的两端都会保留。核心区命中可因此携带另一端的侧翼变异信息。默认 UCE 的 paired-fragment 规则由 `ucefilter` 在同一次扫描中执行。profiling 只使用首次招募，不运行 `refilter`。

## 输入与参考

- 样本表是 tab 分隔的 `sample<TAB>R1<TAB>R2`；单端样本可省略 R2。每个非注释行必须恰为两列或三列，规范化后的样本名必须唯一，且列出的 reads 文件必须存在。
- 输入可为 FASTA、FASTQ 或 gzip 压缩 FASTQ；同一次运行的输入格式必须一致。
- `original` 和 `uce` 使用参考目录：每个 locus 一个 `.fa`/`.fasta`，文件主名就是 locus 名，不能重复。
- profiling 使用单个 marker FASTA 库，而不是 locus 目录。

`-kf` 是参考与 reads 共享的 k-mer 长度；`-s` 是 read 上的采样步长。每个 read 的末端窗口始终额外检查，避免尾端因步长被遗漏。含 N 或其他非 A/C/G/T/U 字符的窗口不会参与匹配。

## 匹配与分配语义

MainFilter 将一个 k-mer 与其反向互补统一为 canonical key，因此一次查询就能覆盖两条链；`-gr` 仅保留为兼容参数。一个 read 可以命中多个 locus，首次招募阶段会写入每个命中 locus；`refilter` 才负责进一步处理这种歧义。

`ref_reads_count_dict.txt` 表示“首次被招募到该 locus 的 reads 数”，不是最终组装覆盖度，也不是基因/等位基因拷贝数。

## 缓存

同一套未变参考可通过 `--reuse-reference-cache` 复用 k-mer 字典。字典保存 canonical 策略、k-mer 长度、locus 名和参考内容 SHA-256；参考内容、k 或格式版本不匹配时会拒绝旧缓存并自动重建，不会静默复用。

扫描步长不改变参考 k-mer 字典本身，只改变 reads 的探测位置。共享缓存目录应只由使用相同参考与 `-kf` 的任务共用。

## 输出与检查

- `ref_reads_count_dict.txt`：首次招募的逐 locus read 数。
- `filtered_pe/`：首次招募 reads；refilter 成功后通常可清理。
- `filtered/`：refilter 后交给组装器的 reads。
- `large_files/`：仅在 refilter 的深度或文件大小限制触发时出现。

某个 locus 没有 reads 并不等于生物学上缺失。应同时检查参考距离、`-kf`、`-s`、测序深度和 `ref_reads_count_dict.txt`。若大量 locus 同时为零，优先检查参考目录、样本表路径和输入 reads 格式。

## Rust MainFilter 相比历史实现

历史 MainFilter 基线由上游维护，本 Rust-only 仓库不再包含其源码。Rust 版保持命令行语义、canonical 双链招募规则和输出格式不变；优化目标是减少每个 read 的 CPU、分配和 I/O 开销，而不是以更大 `-s` 换取速度。

| 部分 | 历史路径的主要成本 | Rust MainFilter 的做法 |
| --- | --- | --- |
| k-mer 扫描 | 解释器循环、对象/字符串处理 | 2-bit rolling k-mer；`k≤32` 用 `u64`，`33–64` 用 `u128`，更长 k-mer 才用字节键。 |
| 双链匹配 | 分别处理正反链 | canonical key，一次索引和一次查询覆盖正反链。 |
| 多 locus 命中 | 多个容器对象与间接访问 | 多命中 locus ID 压入连续 `packed_hits`，查询时直接取切片。 |
| 参考重复使用 | 每个样本重复构建/解析参考 | 带内容 SHA-256 校验的二进制字典缓存；失效时安全重建。 |
| FASTA/FASTQ 读取 | `String`、UTF-8 校验和大小写规范化拷贝 | 字节级 `read_until`，复用行/记录缓冲；DNA 查表同时接受大小写。 |
| gzip 输入 | 默认小缓冲导致较多 FFI/系统调用 | 外层 reader 与 zlib `gzbuffer` 均使用 1 MiB 缓冲；运行时有 zlib-ng 时自动使用，否则回退系统 zlib。 |
| 分 locus 输出 | 频繁打开文件、全局 flush 抖动 | 常驻文件句柄、每 locus 缓冲、高/低水位批量 flush、有限 buffer pool 与大缓冲回收。 |
| 编码与日志 | 临时对象和频繁小写入 | 文本/GM2 编码复用 scratch buffer；日志使用缓冲写入。 |

实际收益取决于 reads、参考大小、k 和文件系统。在 DK40 UCE 的真实测试中，`k=33`、3,579 loci、1 Mi read pairs 的分 locus 输出从 18.35 s 降至 3.72 s，峰值内存从约 375 MiB 降至约 186 MiB；所有输出文件及 read-count 均逐字节一致。更详细的开发验证记录见[性能说明](development/mainfilter-performance.md)。

## 使用原则

1. 先用适合类群的参考与保守的 `-kf`/`-s` 完成准确性验证，再讨论速度。
2. 多样本运行时复用参考缓存；不要手动复制或混用不同参考生成的缓存文件。
3. 将 MainFilter 看作宽松招募，不将其计数直接解释为丰度或缺失。
4. 正式 gene/UCE 工作流使用正常筛选输出；内部仅扫描模式不产生可供组装的 reads，不是推荐分析路径。
