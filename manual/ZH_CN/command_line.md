# TStools 命令行指南

TStools（原 GeneMiner2-UCE）是面向 UCE 和其他系统发育标记的命令行流程。本仓库不包含 GUI、内置演示数据或旧版图形界面文档。

## 1. 从源码构建

### 1.1 构建依赖

完整构建需要：

- C/C++ 编译器和 [zlib](https://zlib.net/)；
- [Rust 和 Cargo](https://www.rust-lang.org/tools/install/)；

Rust/Cargo 是当前完整构建的必需依赖，用于编译主 reads 过滤器、二次过滤器、assembler、population 以及 alignment cleanup、sequence merge、reference trimming 和 UCE statistics 工具。

Ubuntu 可先安装系统依赖：

```bash
sudo apt install build-essential zlib1g zlib1g-dev
```

### 1.2 运行时依赖

系统发育流程按所选步骤调用部分外部程序：

```bash
conda install -c bioconda \
  aster blast clustalo fasttree iqtree mafft magicblast miniprot minimap2 raxml-ng trimal veryfasttree
```

可选工具：

- `--alignment-filter alifilter` 需要单独安装 AliFilter，并确保 `AliFilter` 位于 `PATH`。
- `population` 需要 `minibwa`、samtools、bcftools 和 PLINK 1.9。
- ADMIXTURE 用于遗传成分分析。若未安装，population 仍会完成公共参考、VCF、PLINK 面板和 PCA，并将状态记录为 `unavailable`。

population 外部程序默认从 `PATH` 查找，也可通过 `--population-minibwa`、`--population-samtools`、`--population-bcftools`、`--population-plink` 和 `--population-admixture` 分别指定。

### 1.3 下载与构建

```bash
git clone --depth 1 https://github.com/GUIBA-EX/GeneMiner2-UCE.git
cd GeneMiner2-UCE
cargo run -p xtask -- build
```

构建后的统一入口为：

```bash
cli/geneminer2 -h
```

源码更新后应重新运行相同的 Cargo 命令。

## 2. 输入文件

多数 reads 恢复命令需要以下三个参数；`gene-annotate`、`gene-resolve` 和 `gene-tree` 改用 `--gene-input` 与 `-o`，不需要 `-f/-r`：

- `-f FILE`：tab 分隔的样本表；
- `-r DIR`：参考序列目录；
- `-o DIR`：输出目录。

### 2.1 样本表

每行表示一个样本。单端数据使用两列，双端数据使用三列：

```text
Sample_A	/data/reads/Sample_A_R1.fq.gz	/data/reads/Sample_A_R2.fq.gz
Sample_B	/data/reads/Sample_B_R1.fq.gz	/data/reads/Sample_B_R2.fq.gz
Sample_C	/data/reads/Sample_C_R1.fq.gz	/data/reads/Sample_C_R2.fq.gz
```

建议使用绝对路径。样本名应保持唯一；`population` 可在 R2 后追加可选的 `population` 与 `batch` 列，单端样本若填写元数据则应保留空 R2 列；population 会另行生成内部样本名和 VCF 样本名的对应表。

### 2.2 参考序列目录

每个 locus 使用一个 FASTA 文件，文件名即 locus 名；同一文件可包含一条或多条参考序列：

```text
references/
  uce-0001.fasta
  uce-0002.fasta
```

## 3. 子命令与默认流程

可在同一条命令中按执行顺序列出一个或多个子命令：

| 子命令 | 功能 |
| --- | --- |
| `filter` | 根据参考 k-mer 从原始数据招募 reads |
| `profiling` | 一次招募 marker reads，免组装估计 group-level marker 信号 |
| `refilter` | 进一步分配和过滤每个 locus 的 reads |
| `assemble` | 使用 wDBG 组装目标序列 |
| `gene` | 恢复核基因家族候选 contig |
| `gene-annotate` | 使用蛋白参考进行 miniprot 注释 |
| `gene-resolve` | 比对、建基因树并解析严格一对一子树 |
| `gene-tree` | 从 strict 或 multicopy gene trees 推断物种树 |
| `te` | 从短读长数据发现、整理、注释并定量保守 repeatome 单元 |
| `population` | 构建公共 UCE 参考并生成群体 SNP、PCA 和 ADMIXTURE 结果 |
| `consensus` | 在杂合位点生成一致性序列 |
| `trim` | 按参考序列裁切侧翼 |
| `combine` | 合并样本、比对、清理和过滤 alignment |
| `tree` | 构建溯祖树或串联树 |
| `stats` | 汇总 UCE 恢复统计并可选生成热图 |

不显式指定子命令时：

- `--assembly-mode original`（默认）用于 exon、SCO 及核/线粒体 marker 的参考引导恢复，运行 `filter refilter assemble trim combine tree`；
- `--assembly-mode uce` 用于从 genome skimming 或 target capture 恢复 UCE，运行 `filter assemble combine tree`；融合 UCEFilter 已包含 refilter 语义，并跳过 `trim`，避免新恢复的 UCE 侧翼再次被裁回参考范围；
- `profiling` 先做一次招募，再由 Themisto 伪比对并输出参考序列级支持；不组装，也不运行下游系统发育步骤。

默认 original 模式示例：

```bash
cli/geneminer2 \
  -f /home/user/project/samples.tsv \
  -r /home/user/project/references \
  -o /home/user/project/output \
  -p 8
```

## 4. UCE 组装与 marker profiling

### 4.1 UCE

UCE 模式用于从 genome skimming 或 target capture reads 恢复 UCE core 及有 read 支持的 flank。默认融合 `ucefilter` 在一次扫描中完成招募、run-k 验证和逐 locus 自动选择，低深度 locus 原样通过，饱和 locus 才压缩冗余核心并保留 bait/contig 两端的 overhang 阶梯，直接写最终 per-locus FASTQ；完整 read pair 作为一个单位保留，不写 GM2 或候选 FASTQ。默认流程跳过 `trim`，避免新恢复的 flank 被裁回参考范围。

可选 `--uce-alignment-shadow` 只收集有界的内部比对证据，不改变自动 reads 选择；默认关闭。

```bash
cli/geneminer2 \
  -f samples.tsv -r references -o output -p 8 \
  --assembly-mode uce
```

组装策略、后端选择、rescue、cache 和 QC 见[Assembler 章节](../../docs/assembler_ZH.md)。本手册只保留参数定义；输出字段见[输出文件说明](output.md)。

### 4.2 Marker profiling

Profiling 执行一次招募、Themisto 伪比对并输出参考序列级支持，不组装。它需要一个 `.fa` 或 `.fasta` marker 参考库；`--profile-group-map` 可选，仅增加 group 注释列。

```bash
cli/geneminer2 profiling \
  -f samples.tsv -r marker_reference.fasta \ -o output -p 8
```

输入、decoy、cache、QC 与定量解释见[Profiling 章节](../../docs/profiling_ZH.md)。

### 4.3 Gene 家族恢复与解析

`gene` 是独立完整流程：每个 bait FASTA 为一个 family，可含多个物种。它固定使用 `original-rust`，输出候选而不直接宣称单拷贝。

```bash
cli/geneminer2 gene -f samples.tsv -r family_reference -o gene_output -p 8
cli/geneminer2 gene-annotate --gene-input gene_output/gene \
  --gene-protein-reference family_proteins -o gene_annotation -p 8
cli/geneminer2 gene-resolve --gene-input gene_annotation -o gene_resolved -p 8
```

`gene-resolve` 需要 MAFFT 与 IQ-TREE；可用 `--gene-taper correction_multi.jl` 做 masking。它先按不同样本数和 `--gene-min-aa-length`（默认 30 aa）做 pre-alignment QC，再以 `--gene-min-effective-codon-sites`（默认 30）和占有率做 post-alignment QC；详情见 `occupancy_qc.tsv`。`--gene-ufboot` 只能为 `0`（默认）或 `≥1000`。`family_qc.tsv` 是通过 post-alignment QC 的对齐统计，`tree_selection_qc.tsv` 记录 strict 子树和占有率。

```bash
cli/geneminer2 gene-tree --gene-input gene_resolved -o species_strict -p 8 \
  --gene-species-mode strict --gene-aster astral
cli/geneminer2 gene-tree --gene-input gene_resolved -o species_multi -p 8 \
  --gene-species-mode multicopy --gene-aster astral
```

strict 使用每树一条/样本的 ASTRAL 输入；multicopy 同时传入候选到样本的映射。两条路线均需要 ASTER2 `astral`。

### 4.4 TE / repeatome

`te` 是独立完整流程，不需要 `-r`。它使用四列样本表：`taxon_id sample_id read1 read2`（单端时省略第四列），并按 `discover → curate → annotate → quantify` 运行。

```bash
cli/geneminer2 te -f te_samples.tsv -o te_output -p 32
```

`--te-library` 可选地提供 `name#Class/Subclass` 格式的已分类 TE FASTA；注释不会合并 EQ 或替代完整 TE 注释。输出解释、阈值与重跑规则见 [TE / repeatome 章节](../../docs/te_ZH.md)。

## 5. Population 群体遗传分析

### 5.1 适用范围与示例

`population` 拿多个已完成 UCE 组装的样本及其原始 reads，构建公共伪参考、联合 VCF、PCA 和 ADMIXTURE 输入。每个样本必须保留 `uce_assembly_summary.csv`、`results/` 中已接受的 contig，以及样本表所列 reads。

```bash
cli/geneminer2 population \
  -f /home/user/project/samples.tsv \
  -r /home/user/project/references -o output -p 8 \
  --assembly-mode uce \
  --population-admixture-k-min 2 --population-admixture-k-max 6
```

公共伪参考策略、分阶段重启、SNP 面板和必查 QC 见[Population 流程说明](../../docs/population_ZH.md)。个体缺失率、杂合度和成对 PI_HAT 写入 `population/structure/qc/`；大队列可用 `--skip-relatedness-qc` 跳过亲缘关系计算。仅在对应阶段产物已经通过校验时使用 `--population-start-at`。

## 6. 分步运行示例

只基于已有结果重建串联树：

```bash
cli/geneminer2 tree \
  -f samples.tsv -r references -o output \
  -m concatenation --phylo-program iqtree
```

生成 consensus，按参考裁切并合并：

```bash
cli/geneminer2 consensus trim combine \
  -f samples.tsv -r references -o output \
  -c 0.75 -ts consensus -tm all -tr 0.5 \
  -cs trimmed
```

从已有 UCE 结果生成统计表但不绘图：

```bash
cli/geneminer2 stats \
  -f samples.tsv -r references -o output \
  --stats-no-heatmap
```

## 7. 参数参考

以下列出主要公开参数及当前默认值。运行 `cli/geneminer2 -h` 可查看同版本的完整帮助。

### 7.1 通用输入与并行

| 参数 | 说明 |
| --- | --- |
| `-f FILE` | 样本表，必需 |
| `-r DIR` | 参考序列目录，必需 |
| `-o DIR` | 输出目录，必需 |
| `-p INT\|auto` | 共享 CPU 预算，默认 `auto`。自动模式统计 affinity/cpuset 允许的物理核心，并受 cgroup 或作业调度器上限约束；整数值可覆盖自动检测。 |

### 7.2 Reads 过滤与 re-filtering

| 参数 | 说明 |
| --- | --- |
| `-kf INT` | 过滤 k-mer 大小；UCE 模式默认 `23`，其他模式默认 `31` |
| `-s, --step-size INT` | reads 扫描步长，默认 `4` |
| `--max-reads INT` | 每个文件最多处理的 reads 数，单位为百万；`0` 表示不限 |
| `--reuse-reference-cache` | 复用带输入指纹的 reference k-mer index；显式使用 `original-rust` 时也启用带版本和 k 校验的 assembler 二进制 cache |
| `--reference-cache-dir DIR` | reference cache 目录；默认 `output/.gm2_reference_cache`，且必须与上一参数同时使用 |
| `--depth-low-water-mark INT` | 低于该深度时尝试放宽条件招募更多 reads，默认 `50` |
| `--depth-limit INT` | re-filtering 处理的最高深度，默认 `768` |
| `--file-size-limit INT` | re-filtering 文件大小限制，默认 `6` |

### 7.3 组装与 UCE

| 参数 | 说明 |
| --- | --- |
| `-ka INT` | assembly k-mer 大小；默认 `0`，表示自动估算 |
| `--min-ka INT` / `--max-ka INT` | 自动估算范围，默认 `21` / `51` |
| `-e, --error-threshold INT` | k-mer 错误阈值，默认 `2` |
| `-sb, --soft-boundary VALUE` | 软边界：整数、`auto` 或 `unlimited`；默认 `auto` |
| `-i, --search-depth INT` | 搜索深度，默认 `4096` |
| `--min-coverage INT` | contig 最低 read depth，默认 `0` |
| `--assembler-implementation MODE` | `auto`（默认）在 original 模式使用 `original-rust`，在 uce 模式使用 `uce-rust`；`original` 与 `original-rust` 均使用单线程、确定性的 Rust 原版兼容实现；uce 不会回退至其他实现 |
| `--assembler-read-chunk-size INT` | Rust assembler 每批读取的 reads 数，默认 `8192` |
| `--assembler-kmer-count-threads INT` | 每个 locus 的 k-mer 排序和计数线程；默认 `0`，表示自动分配 |
| `--assembler-graph-format MODE` | 可选组装图输出：`none`（默认）、`gfa`、`dot` 或 `both` |
| `--assembly-mode MODE` | `original` 或 `uce`；默认 `original` |
| `--assembly-mode uce` | 默认使用 k=23/step=4、自动敏感招募、一轮受控 rescue，以及固定的 backbone/QC 安全设置 |
| `--uce-recruit-mode fast\|auto` | 招募策略；UCE 模式默认 `auto`，其他模式默认 `fast` |
| `--uce-rescue-reads` | 显式启用固定 k=21 的受控 rescue；为兼容旧命令而保留，UCE 模式已默认启用 |
| `--no-uce-rescue-reads` | 关闭 UCE 模式默认的 rescue 阶段 |
| `--uce-rescue-rounds 1\|2` | rescue 轮数；默认 `1` |
| `--uce-rescue-inverted-repeat-min-bp INT` | 若本轮 rescue 新引入至少该长度的精确倒置重复，则逐 locus 回退；默认 `150`，`0` 关闭 |
| `--uce-rescue-reverse-reuse-reference-scale FLOAT` | 实验性：缩放与任一组装臂已走路径互为反向互补节点的 reference bonus；范围 `0`--`1`，默认 `1.0`（关闭缩放），不改变 reads depth |

### 7.4 Population

| 参数 | 说明 |
| --- | --- |
| `--engine MODE` | `pseudoref`（默认）、旧 `panref` 或 `panrefv2`；后两者使用 `-r` 指定的按locus bait目录 |
| `--population-panrefv2-include-low-confidence` | 将 `short` 或 `low_sample_support` PanRefV2 locus 也写入参考；默认仅保留 `pass` |
| `--population-reference-strategy MODE` | 仅 `pseudoref` 使用：`sqcl-longest`（默认）或 `supported` |
| `--population-reference-fasta FILE` | 使用固定外部 FASTA 作为公共参考；文件会复制到 `population/reference/`，不生成按样本的参考贡献统计 |
| `--population-min-mapq INT` / `--population-min-baseq INT` | 联合检测最低 MAPQ / base quality，默认 `20` / `20` |
| `--population-min-dp INT` / `--population-min-gq INT` | 低于阈值的基因型设为缺失，默认 `5` / `20` |
| `--population-min-qual FLOAT` | 最低位点 QUAL，默认 `20` |
| `--population-min-call-rate FLOAT` | 最低非缺失基因型比例，默认 `0.8` |
| `--population-min-mac INT` | 最低 minor allele count，默认 `2` |
| `--population-ld-window INT` / `--population-ld-step INT` | LD pruning 窗口和步长，默认 `50` / `5` SNP |
| `--population-ld-r2 FLOAT` | LD pruning 的 r² 阈值，默认 `0.2` |
| `--population-admixture-k-min INT` / `--population-admixture-k-max INT` | ADMIXTURE K 范围，默认 `2` / `6`；最大 K 不超过样本数 |
| `--population-admixture-cv INT` | ADMIXTURE CV 折数，默认 `10`，不超过样本数 |
| `--population-start-at STAGE` | 从 `reference`（默认）、`mapping`、`calling` 或 `selection` 开始；后 3 者分别复用已检查的参考、BAM 或过滤 VCF |
| `--population-stop-after STAGE` | 在 `reference`、`mapping`、`calling` 或 `selection` 后停止；默认 `selection` |
| `--population-skip-mark-duplicates` | 跳过 samtools duplicate marking |
| `--population-skip-plink` | 不生成 PLINK、PCA、LD-pruned 或 ADMIXTURE 结果 |
| `--population-skip-admixture` | 生成 PLINK 和 PCA，但不运行 ADMIXTURE |
| `--population-minibwa PATH` | minibwa 可执行文件，默认 `minibwa` |
| `--population-samtools PATH` | samtools 可执行文件，默认 `samtools` |
| `--population-bcftools PATH` | bcftools 可执行文件，默认 `bcftools` |
| `--population-plink PATH` | PLINK 1.9 可执行文件，默认 `plink` |
| `--population-admixture PATH` | ADMIXTURE 可执行文件，默认 `admixture` |

### 7.5 Consensus、裁切与合并

| 参数 | 说明 |
| --- | --- |
| `-c, --consensus-threshold FLOAT` | consensus 阈值，默认 `0.75` |
| `-ts, --trim-source SOURCE` | `assembly` 或 `consensus`；默认使用上一步结果 |
| `-tm, --trim-mode MODE` | `all`、`longest`、`terminal` 或 `isoform`；默认 `terminal` |
| `-tr, --trim-retention FLOAT` | 参考裁切的保留长度比例，默认 `0` |
| `-cs, --combine-source SOURCE` | `assembly`、`consensus` 或 `trimmed`；默认使用上一步结果 |
| `-cd, --clean-difference FLOAT` | alignment 最大可接受成对差异，默认 `1.0` |
| `-cn, --clean-sequences INT` | alignment 最少序列数，默认 `0` |
| `--msa-program PROGRAM` | `clustalo` 或 `mafft`；默认 `mafft` |
| `--msa-threads INT` | 每个 MSA 任务的线程数，默认 `1`，且不能大于 `-p` |
| `--alignment-filter PROGRAM` | `trimal`、`alifilter` 或 `none`；默认 `trimal` |
| `--filter-processes INT` | 并行 alignment-filter 任务上限，默认等于 `-p` |
| `--alifilter-model MODEL` | AliFilter 模型名或 `model.json` 路径 |
| `--strict-combine-errors` | 任一 locus 的 alignment、清理或过滤失败时立即停止 |
| `--no-alignment` | 跳过多序列比对 |
| `--no-trimal` | 已弃用；等同于 `--alignment-filter none` |

### 7.6 建树与统计

| 参数 | 说明 |
| --- | --- |
| `-m, --tree-method METHOD` | `coalescent` 或 `concatenation`；默认 `coalescent` |
| `-b, --bootstrap INT` | bootstrap 次数，默认 `1000` |
| `--phylo-program PROGRAM` | `raxmlng`、`iqtree`、`fasttree` 或 `veryfasttree`；默认 `fasttree` |
| `--stats-no-heatmap` | 不生成 UCE 统计热图 |
| `--stats-count-input-reads` | 统计输入 FASTQ reads 以填写 `InputReads` 和 `PctFiltered`；大数据集上较慢 |

比例类参数使用 `0–1` 小数。对正式数据修改过滤阈值前，建议先保留默认值完成小规模测试，并同时检查中间 VCF、mapping QC 和 locus 占有率。
