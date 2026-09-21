# MainFilter 性能优化与兼容性说明

[English](mainfilter-performance_EN.md)

本文记录当前 Rust `MainFilterNew` 与上游历史 Haxe/C++ 基线的真实 UCE 对比、性能设计与兼容性边界。它只讨论 reads 招募的工程性能；不将运行时间等同于 UCE recovery、组装准确性或群体遗传推断质量。

## 结论：真实 Haxe/C++ 对比

历史 MainFilter 基线为 Haxe 编译的 C++/hxcpp 实现，源码仅在上游仓库保留。本次以该上游源码重新编译历史二进制，并与当前 Rust 二进制在完全相同的真实 UCE 数据、参考、参数和输出模式下比较。

| 情形 | 历史 Haxe/C++ | 当前 Rust | 改善 |
| --- | ---: | ---: | ---: |
| 首次建 cache + 过滤 | 20.32 s | 3.52 s | **82.7% 更快**（5.8×） |
| 复用 cache + 过滤 | 18.35 s | 3.72 s | **79.7% 更快**（4.9×） |
| 首次建 cache 的峰值 RSS | 381 MiB | 187 MiB | **51.0% 更低** |
| 复用 cache 的峰值 RSS | 375 MiB | 186 MiB | **50.3% 更低** |

过滤阶段本身由 Haxe/C++ 的 18.59 s（首次）/18.06 s（复用）降至 Rust 的 3.42 s/3.57 s。这里的 wall time 包含 gzip 解压、FASTQ 解析、k-mer 招募、按 locus 输出与计数写入。

## 基准协议

### 实现

- 历史基线：上游 `sculab/GeneMiner2` 的 MainFilter 源码，使用 Haxe 4.3.3、hxcpp 4.3.2、`g++ -O2` 和项目 README 指定的编译选项重新构建。
- 候选实现：本仓库 `rust/main_filter_new` 的 release `MainFilterNew`。
- 两版各自创建并加载自己的 dictionary cache；两种 cache 格式不同，不能交叉复用。

### 固定输入

- reads：`validation/reads/fastq/DK40__SRR29729138_{1,2}.fastq.gz`
- reference：`validation/references/bivalve_uce_2k_v1_loci/`，共 3,579 个 locus。
- 参数：`-kf 33 -s 3 -m_reads 1 -m 0 -gr`。
- `-m_reads 1` 表示每个输入文件最多处理 `2^20` 个记录，即 1 Mi read pairs。
- 输出：`filtered_pe/`，每 locus 一个文本 read 文件。

每个二进制分别运行两次：第一次无 cache，第二次加载本实现第一次生成的 cache。计时使用 `/usr/bin/time`，报告 wall time、user/system time 与最大 RSS。

## 结果一致性

两版的招募与输出通过如下检查：

| 检查 | 结果 |
| --- | --- |
| 非空 locus 输出文件数 | 两版均为 2,228 |
| 输出总字节数 | 两版均为 460,812,556 bytes |
| `filtered_pe/` 逐文件比较 | `diff -qr` 完全一致 |
| `ref_reads_count_dict.txt` 数值 | 完全一致 |

最后一项的原始文件不能直接 `cmp`：历史 Haxe 版按无序哈希表顺序写入，且每行格式是 `locus,count,`（末尾多一个逗号）；Rust 按参考顺序写入 `locus,count`。去掉 Haxe 的尾随逗号并排序后，2,228 个 locus 的 count 完全一致。这是格式与行序差异，不是 read recruitment 差异。

## Rust 版的优化

| 部分 | 历史 Haxe/C++ 路径的主要成本 | 当前 Rust 的做法 | 语义边界 |
| --- | --- | --- | --- |
| k-mer 表示 | 长模式与通用容器开销 | 2-bit rolling key：`u64`（`k≤32`）、`u128`（`33–64`）；仅 `k>64` 用字节键。 | canonical key、歧义碱基和 read 末端检查规则保持。 |
| 双链匹配 | 正反链 k-mer 处理 | canonical k-mer，一次查询覆盖正反链。 | `-gr` 仅保留为兼容参数。 |
| 多 locus 命中 | 分散容器与间接访问 | `ReferenceHits` + 连续 `packed_hits`；结构实测为 12 bytes。 | locus ID 集合不变。 |
| 索引查找 | 通用散列表与对象访问 | `AHashMap` 的原生整数/字节 key。 | 哈希表迭代顺序不作为输出语义。 |
| dictionary cache | 旧格式、较多逐项 I/O | v4 cache、参考内容 SHA-256、4 MiB `BufWriter`、加载时直写 `packed_hits`。 | cache 不匹配、损坏或过期时安全重建。 |
| FASTA/FASTQ | Haxe 字符串与对象路径 | 字节级 `read_until`、复用行与记录缓冲，不进行全 read 大写复制。 | 正常 FASTA/FASTQ 输出字节不变。 |
| gzip | 历史 `GzipReader` 路径 | 外层 reader 与 zlib `gzbuffer` 均为 1 MiB；构建时通过 `pkg-config zlib-ng` 自动直连原生 zlib-ng，未检测到时保留运行时 zlib-ng 探测并最终回退系统 zlib。 | 解压得到的字节不变；zlib-ng 仍是可选构建依赖。 |
| 输出 | buffer flush 与反复文件操作 | 常驻句柄、每 locus 缓冲、64 MiB 高水位/32 MiB 低水位批量 flush。 | 每个文件的记录内容与顺序不变。 |
| 内存回收 | 热 locus 峰值容量易常驻 | 有界 buffer pool；大于 1 MiB 的 buffer 在 flush 后释放。 | 不影响已写或后续记录。 |
| 编码/日志 | 临时对象与小写入 | 文本/GM2 scratch buffer 复用、缓冲日志。 | 编码格式与日志语义不变。 |

## cache 行为

dictionary cache 只依赖参考内容、locus 名、k-mer 长度与索引格式；扫描步长 `-s` 不改变参考 dictionary，只改变 reads 上的探测位置。共享 cache 只能用于相同参考与相同 `-kf` 的工作。

Rust cache 当前为 v4。旧 Rust cache、历史 Haxe cache、损坏 cache 或参考内容发生变化时都会被拒绝并重建；不允许静默读错索引。

## 不作为正式工作流性能结论的内容

- **scan-only mode**：不产生可组装的 reads，仅用于内部诊断，不作为 gene/UCE 正式性能结论。
- **增大 `step`**：会改变采样密度与潜在检出率，属于分析参数改变，而不是无损优化。
- **正式分 locus 输出的并发写入**：多线程向同一组 locus 文件写入会引入合并、顺序和 I/O 争用问题；当前正式路径优先保留单写入者与逐字节验证。
- **SSHash/GGCAT 等索引替换**：会改变静态索引与 cache 设计，只有证实现有索引是主要瓶颈后才应独立实现并重新验证。

## 修改后的验证要求

修改 MainFilter 后至少执行：

```bash
cargo test --manifest-path rust/main_filter_new/Cargo.toml
cargo clippy --manifest-path rust/main_filter_new/Cargo.toml -- -D warnings
cargo build --release --manifest-path rust/main_filter_new/Cargo.toml
```

凡涉及招募逻辑、k-mer 编码、cache、解析或输出，还必须用固定真实双端样本进行旧/新二进制比较：

```bash
diff -qr old_run/filtered_pe new_run/filtered_pe
# 对 count 文件：先统一尾随逗号格式并按 locus 排序，再比较数值。
```

只有在输出一致、测试通过且真实数据计时有稳定收益时，才应发布候选二进制。
