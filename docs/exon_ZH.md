# Exon/intron 注释

[English](exon_EN.md) · [Gene 工作流](gene_ZH.md) · [命令行指南](../manual/ZH_CN/command_line.md)

`--assembly-mode exon` 在默认 gene family 恢复流程之后增加蛋白引导的 exon/intron 注释。默认直接从 coding 核酸 bait 推导蛋白参考，不要求另备 `.faa`；外部蛋白可作为可选覆盖。该模式不用于 UCE。

## 输入与命令

`-r` 目录中每个 `.fa`、`.fas` 或 `.fasta` 文件定义一个 family，可包含多个物种的 CDS、拼接 exon 或其他保持 coding frame 的核酸 bait。

```text
family_reference/
├── geneA.fasta
└── geneB.fa
```

默认流程使用标准核遗传密码表自动翻译这些核酸参考。若 bait 包含 intron、属于非编码序列、使用非标准遗传密码或不能解释为 coding sequence，应提供外部蛋白目录；文件必须以同名 family stem 加 `.faa` 命名。目录可以不完整：有 `.faa` 的 family 使用外部蛋白，其余 family 仍自动翻译。

exon 模式是一条完整工作流；不要再添加 `filter`、`assemble` 等阶段子命令。

```bash
cli/tipseek --assembly-mode exon \
  -f samples.tsv \
  -r family_reference \
  -o exon_output \
  -p 8
```

可选的外部蛋白覆盖为：

```bash
cli/tipseek --assembly-mode exon \
  -f samples.tsv -r family_reference \
  --gene-protein-reference family_proteins \
  -o exon_output -p 8
```

注释要求 miniprot 0.18 或更高版本。TipSeek 先将恢复候选写入 `exon_output/gene/`，再将结构结果写入 `exon_output/exon/`。

## 实现方法

1. 默认 gene 工作流完成 reads 招募、refilter、组装及跨样本候选汇总。
2. 未提供同名 `.faa` 时，TipSeek 平等检查正反两条链的六个 reading frame。允许并移除正常末端终止；任何含内部终止的 frame 都直接排除，不修补也不屏蔽后继续使用。
3. 只有一种无内部终止 frame，或某个 frame 从起始甲硫氨酸和/或末端终止获得唯一最强的 CDS 端点支持时，才直接接受该序列；其中最长的可信翻译作为 family anchor。其余歧义序列只有在某个 frame 对该 anchor 获得唯一且大于零的氨基酸 3-mer Dice 相似度时才接受。没有可信 anchor、最佳结果并列或完全不相似时，不为该序列输出蛋白，并记录 `ambiguous_reading_frame`；若该 family 没有其他可用参考，则需要提供同名 `.faa`。所有候选和最终选择写入 `manifest/reference_translation.tsv`。
4. 对每条候选 contig 独立运行 miniprot，避免一条 partial candidate 通过 target competition 压制另一条互补 candidate。
5. 嵌入 PAF 提供蛋白坐标、覆盖度、score、identity、protein CIGAR、frameshift 与 in-frame stop；GFF3 提供 exon/CDS 坐标、方向和 phase。
6. 重叠模型按结构缺陷、蛋白覆盖度、alignment score、identity、剪接证据与 miniprot rank 确定性竞争；不重叠模型保持独立。
7. 对选中模型分类，写入结构化审计表，并分别进入 resolve-eligible 或 unresolved 输出。

基因组输出保留 `N` 与其他 IUPAC 碱基，因此所有坐标始终对应实际用于注释的候选序列。

## 模型接纳条件

只有已选中的 `complete`、`terminal_partial` 与 `internal_partial` 模型可以进入 `cds/` 和后续 `gene-resolve`。合格模型必须满足：

- 无 frameshift 或内部终止；
- CDS phase 明确且连续；
- CDS 长度为 3 的倍数；
- 翻译结果不含 `X`；
- 不含未识别的非经典剪接；
- 不存在未解决的模型竞争。

注释 CDS 保留一个正常末端终止密码子；仅在准备蛋白与密码子比对时，才将它与对应末端 `*` 一并移除。其他状态，包括 `low_coverage`、`frameshifted`、`internal_stop`、`phase_inconsistent`、`ambiguous_cds` 和 `ambiguous_model`，保留在 `unresolved/` 与 manifest 中。

## 短外显子与多 contig 恢复

exon 模式会继续搜索剩余参考种子，并从已有 reads 的高支持图节点恢复无参考种子的非分叉片段，保留至少 30 nt 的候选；普通 gene 模式保留原有搜索和长度门槛。miniprot 采用 `-L 10 -n 2 --outc 0` 保留短命中，低覆盖片段不会因此自动进入 resolve。

片段先在同一蛋白参考下比较，再按参考坐标构建有向无环图。参与串联的片段至少覆盖 10 aa、identity ≥80%，且无移码、内部终止、未知编码碱基、内部 query 缺口或剪接缺陷。可组合两个及以上 contig，包括位于蛋白中间的短片段；同一 contig 不能重复进入一条链。相邻片段 identity 差不超过 15 个百分点；重叠 CDS 必须按参考坐标逐碱基一致才可去重，不再使用固定的 12 aa 重叠上限；含插入/缺失、无法对应编码坐标的重叠不裁剪。

存在覆盖接近的替代链时标记 `ambiguous`。每个参考最多 128 个节点、每个末端最多 64 条路径；搜索被截断时标记 `search_limit`，不据此自动接纳。跨不同参考物种、但 contig 顺序和方向相同的链会合并记录。这些限制控制组合规模，也意味着复杂重复区域可能仍未解决。

两种输出分别表达不同证据：

- `fragment_cds/`：按参考排序的 CDS 串联假设，标记 `reference_ordered_unvalidated`，不会自动进入 resolve。内部缺失编码区按缺失 aa × 3 补 N；这些 N 不代表内含子。无法可靠按 query 坐标切割的片段不导出该假设。
- `supercontigs/` 与正常 `cds/`：对唯一、覆盖足够且到达参考两端的链，在每两个已定向来源 contig 之间插入 `--gene-fragment-padding` 个 N（默认 100），调用所选后端重新注释。只有唯一完整、结构可靠的正向模型可以接纳；每段 N 都必须位于内含子内、不能与 exon 重叠，且每个来源片段必须保留编码支持，不能跳过中间片段。

`manifest/fragment_chains.tsv` 保留来源成员、片段数、query 覆盖、编码缺口、重叠去重长度和状态；`physical_linkage=unverified` 明确表示未证明这些 contig 在真实基因组中的物理连接。结构重注释通过也不等同于同一拷贝或读段连接证明。`fragment_groups.tsv` 保留拼接首尾端点，完整成员以 `fragment_chains.tsv` 为准。

缺少剪接边界、存在拷贝冲突或尚未组装出的外显子，都可能导致拼接失败。补 N 不估计真实内含子长度，也不会恢复不存在的编码序列。`--gene-fragment-padding 0` 关闭拼接输出，诊断表仍保留候选链。

## CDS 完整度

覆盖率扣除 miniprot 的 query 插入（`I`）或 Exonerate 的 query-only `G`。到达参考两端但内部有缺口的模型标为 `internal_partial`，不再误标 `complete`，也不能通过完整补 N 验证。该标记描述相对参考的差异，不代表该物种一定缺失了外显子。


## reads 支持与新增片段预筛

无参考核苷酸种子的图片段不会跨分叉猜测连接，来源在 `candidate_assembly.tsv` 中标为 `assembler_seed_position=-2`。仅此类新增候选先做六读框蛋白预筛：至少一个 10 aa 窗口有 8 个一致残基，双方窗口均至少含 4 种氨基酸。未通过的候选记录为 `protein_screen_rejected`，原序列保留在 gene 候选文件；原有参考锚定候选不受该预筛影响。这会保守地放弃只有低复杂度或不足 10 aa 连续编码证据的独立救援片段，不等同于证明该片段非编码。

默认从工作流输出下的 `样本/filtered_pe/家族_1、_2` 读取配对 reads；分开运行注释时可用 `--gene-reads-root` 指向原工作流目录（底层 `gene_workflow annotate --reads-root`）。支持本流程成对输出的 GM2，以及匹配名称/数量的双端 FASTQ，包括 gzip；GM2 不含原始名称，只能检查数量并依赖生产步骤保留的配对顺序。没有文件时记为 `reads_unavailable`，不会假装验证通过；缺失单侧或 FASTQ 配对不一致时报错。

检查在一个家族的所有原始候选上只建立一次索引，不将 reads 比对到人工 N 边界。当前采用保守的全长、无 indel 比对：reads 至少 90 nt，最多 4 个替换，且与次优位置至少差 2 个替换。重复种子或搜索超限不判断唯一性。相同位置和方向的 read pair 去重；用于连接判断的两端都必须完全处于各来源的外显子内，并在定向片段中相向。至少两个不同位置的低错配 pair 记为 `paired_read_support`，仍不是物理连接证明。没有桥接短 reads 时记为证据不足，不因超长内含子无法跨越就判错。

借鉴 HybPiper 的不对称错配检查：一端全匹配、另一端有 3–4 个替换，累计至少 5 个去重 pair 且多于低错配支持时，记为 `discordant_read_warning`，阻止该链导出和自动补 N 接纳。该告警不是错误合并率；当前不建模碱基质量、indel 或插入片段长度，不是完整 reads 比对器。方向冲突另行报告。

两种后端均采用 GeneMiner2 式的 SCO 任务池：`-p` 控制同时处理的 SCO 数量（最多 64），每个任务处理一个完整 SCO，完成后立即领取下一个；每次比对使用一个线程；队列收尾时，剩余 SCO 可在共享 CPU 总预算内借用空闲核处理候选，不把 SCO 并发数与候选并发数相乘。各任务使用独立目录，最终按样本、SCO 和候选顺序合并。实际调度记录在 `manifest/execution.tsv`，完成进度定期写入日志。

新增 `manifest/read_audit.tsv`、`fragment_read_support.tsv` 和 `completion_audit.tsv`。内部拼接保留模型两侧最多 120 nt 真实侧翼，外部末端保持原有序列；实际裁剪坐标记录在 `join_sources.tsv`。补 N 后仍需通过完整模型和全部来源的结构复核。

## 主要参数

| 参数 | 默认值 | 含义 |
| --- | ---: | --- |
| `--gene-protein-reference` | 自动推导 | 可选 `.faa` 目录；同名 family 覆盖自动翻译 |
| `--gene-miniprot` | `miniprot` | miniprot 可执行文件 |
| `--gene-max-intron` | `50000` | miniprot 接受的最大 intron 长度 |
| `--gene-min-model-coverage` | `0.20` | 常规 partial model 的最低蛋白覆盖度；仅补 N 尝试允许 1-aa 离散容差 |
| `--gene-complete-coverage` | `0.80` | complete model 最低覆盖度；仍要求覆盖蛋白两端 |
| `--gene-fragment-padding` | `100` | 每个相邻 contig 之间插入的 N 数；`0` 关闭拼接 |
| `--gene-reads-root` | 当前工作流输出 | 含各样本配对招募 reads 的工作流目录 |
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
    ├── manifest/
    │   ├── reference_translation.tsv # 自动 frame/strand 选择审计
    │   ├── derived_proteins/     # 实际使用的自动推导蛋白
    │   └── ...                   # candidate、ID、warning、fragment 与命令审计
    └── unresolved/
```

参考来源和 frame 选择查看 `manifest/reference_translation.tsv`；模型状态与 eligibility 查看 `models/gene_models.tsv`，坐标与 splice class 查看 `models/gene_segments.tsv`，所有补 N 尝试查看 `manifest/fragment_groups.tsv`。后续运行：

```bash
cli/tipseek gene-resolve --gene-input exon_output/exon -o gene_resolved -p 8
```

resolve QC 以及 strict/multicopy 物种树推断见 [gene 工作流说明](gene_ZH.md)。

## 可选 Exonerate 对照（第一、二阶段）

增加 `--gene-exonerate-shadow` 可对每条候选及已尝试的补 N 派生序列运行独立的 protein2genome 比对。此选项仅在 exon 模式可用；未指定时不调用 Exonerate，也不产生对照文件。默认通过当前环境的 PATH 查找 `exonerate-rs`，无需指定路径；当前不随 TipSeek 打包。

```bash
cli/tipseek --assembly-mode exon -f samples.tsv -r family_reference \
  -o exon_output -p 8 \
  --gene-exonerate-shadow
```

这一阶段仍完全由 miniprot 及现有质控决定 CDS、模型接纳和补 N 结果。对照仅使用 Exonerate 的 `protein2genome --exhaustive no --subopt no --score 1 --revcomp yes`，最大 intron 沿用 `--gene-max-intron`。每次比对最长 60 秒，进程失败、超时或无效输出写入审计，不作为“无命中”，也不替代主模型。文件系统写入等工作流错误仍会正常报错。

新增输出：

- `manifest/exonerate_shadow.tsv`：逐蛋白汇总 miniprot 原始模型数量、Exonerate 命中数量、匹配数量及各 Exonerate 命中的分数、坐标、CDS 区间和移码操作数。`stage` 区分 `candidate` 与 `padded_join`；所有区间均为 0-based、左闭右开，负链区间也按基因组坐标升序表示。
- `raw_exonerate_shadow/`：可执行文件版本响应、每次运行的参数、输入快照、stdout、stderr、执行状态与耗时，便于在 `.work` 清理后复查。

`reported_hits_match` 只表示每条**已报告的 Exonerate 命中**都能找到蛋白区间、方向、基因组跨度和 CDS 区间相同的 miniprot 原始模型，不表示双方候选集合、CDS phase、identity、结构质量或最终接纳结果等价。Exonerate 当前只报告最佳路径，miniprot 可报告多个模型；移码操作数也不保证与 miniprot 计数定义一致。`coordinate_difference` 是待核查差异，不自动判定哪一方正确。执行失败、仅一侧命中和双方无命中分别报告。

现已统一带引擎来源的比对证据与结构检查；后续阶段将验证模型竞争和指标定义，再开放独立后端选择。

第二阶段另在每个对照目录写入 `normalized_models.tsv`，通过统一结构检查生成 CDS 阅读框、翻译蛋白及质量状态。identity 明确按 M/S 对齐氨基酸对计算，排除缺口；替换矩阵 positive 比例缺失时记为 `.`。这些模型只用于审计，不参与模型选择和补 N 决策。

完整状态现在要求比对覆盖参考蛋白两端；少量末端缺失也会标为 `terminal_partial`，并记录 `query_start_truncated` 或 `query_end_truncated`。参考以 M 开头而重建蛋白起始不一致时，记录 `reference_start_mismatch`；这类部分模型不能验证完整补 N 拼接。`complete` 是结构与参考覆盖分类，不是生物学准确性的保证，miniprot 也不作为真值。

## 独立 Exonerate 后端

在 exon 命令中添加 `--gene-exon-backend exonerate-rs`，即可让 Exonerate 独立负责候选注释、模型选择及补 N 重注释。程序直接从当前环境 PATH 调用 `exonerate-rs`，不接受或要求 Exonerate 路径；此模式不需要 miniprot。未指定后端时仍使用 miniprot。此选项不能与 `--gene-exonerate-shadow` 同时使用。

```bash
cli/tipseek --assembly-mode exon -f reads.tsv -r references -o exon_output -p 8 \
  --gene-exon-backend exonerate-rs
```

组件也可直接处理已有候选：

```bash
gene_workflow annotate --backend exonerate-rs --input gene_output/gene \
  --nucleotide-reference references --out exon_output/exon
```

两种后端使用共同的结构质控和拼接验证。`raw_exonerate/` 保存原始 vulgar、输入快照、参数、退出状态和 stderr；`manifest/backend.tsv` 记录引擎、版本及 identity 口径，`command_provenance.tsv` 记录实际调用。Exonerate 的 identity 按 M/S 对齐氨基酸对计算，未知 positive 值记为 `.`，不会伪造 PAF。常规模型、CDS、蛋白、GFF3 和 resolve 输入布局保持一致。

当前采用 `protein2genome --exhaustive no --subopt no --revcomp yes`，每次比对有 60 秒超时；`-p` 用于 SCO 任务池并行，最多 64 个 SCO 同时注释；SCO 内的 Exonerate 比对逐个独立执行，单次比对保持单线程。程序不可用时会在清理旧输出前报错；单个候选运行失败标记为 `exonerate_failed`，无命中单独记录，无法解析的输出会终止任务。当前最佳路径搜索不保证找全重复位点，其候选集、打分和末端处理与 miniprot 不等价。后端速度和接纳数量不能用来代替准确性验证。

补 N 复注释的逐模型诊断保存在 `manifest/padded_models.tsv`，包含参考蛋白覆盖区间、模型状态、质量标记和未通过的拼接判定条件；可据此区分缺失首尾残基、剪接缺陷和跨 N 区域失败。两种后端的 CDS 可能因是否保留末端终止密码子而相差 3 nt，不能仅凭文件不一致认定剪接结果不同。


## 按 SCO 汇总外显子

主要恢复结果为 `exon/sco/<SCO>.fasta`：每个样本、每个 SCO 一条序列，按固定参考编码坐标排列；缺失或不同候选间不能消除的碱基冲突用 `N`。不会用参考碱基补缺，也不会因整条 CDS 未通过模型选择而丢弃结构合格的短片段。该序列是参考坐标上的证据汇总，不代表已验证的单一拷贝、分相 CDS 或真实跨 contig 连接。

`--gene-sco-reference reference_id` 固定该 SCO 文件内的参考记录；底层对应 `gene_workflow annotate --sco-reference reference_id`。未指定时按有效命中的唯一残基覆盖与一致性加权选择参考，并限制异常参考长度的影响；同一个 SCO 在整个样本组中使用同一参考。匹配必须唯一，不匹配时报错。其他参考的命中通过保守蛋白坐标映射投影到该参考：插入、重复区多解和弱锚点不填入样本碱基，同一原始密码子对应多个位置时填 N。选择参考时重复或被包含的命中不会重复增加支持度，冲突证据仍保留为 N。

两后端保留逐残基到样本 CDS 的映射，插入不会推移下游坐标，参考缺失位置保留 N，跨剪接位点的密码子仍占一个参考密码子。参考编码区中的样本插入不进入固定长度汇总，但原始模型保留。参考末端含终止密码子时保留对应三个位点，只填入实际观察到的样本终止密码子，否则用 NNN；不会复制参考终止密码子。若参考经过自动反向或移框翻译，坐标以所选编码翻译为准，详见 `reference_translation.tsv`。

审计文件：`manifest/sco_recovery.tsv` 记录参考、实测/缺失/冲突碱基和无法映射的模型数；`exon_candidates.tsv` 与 `exon_candidate_sources.tsv` 保留坐标、模型来源、边界与 reads 支持，不另拆分交付外显子 FASTA。单端唯一比对也计入外显子支持，但不等于双端跨 contig 证据。内部剪接边界和未验证的比对终点分开记录。旧的模型、GFF、`exons/` 等目录保留用于诊断，主要汇总结果请使用 `sco/`。


复用已有比对结果修复 SCO 汇总，无需调用任一比对后端：

```bash
gene_workflow recover-sco --annotation output/exon --input output/gene \
  --out recovered_sco --threads 16
```

使用 `--sco-backbone-manifest output/exon/manifest/sco_recovery.tsv` 固定每个 SCO 原有参考，便于比较投影前后的恢复率。该操作读取保存的原始比对、候选 DNA 和参考蛋白，不重新招募 reads、组装或运行比对。缓存不完整时会报错；应提供原注释对应的候选 DNA。`manifest/reference_projection.tsv` 记录映射状态，`sco_selection_control.tsv` 提供仅修复参考选择的对照。

对跨 contig 的冗余证据，只在短 contig 全长（至少 30 nt）唯一、精确地包含在长至少 30 nt 的另一条实测 contig 中时比较竞争模型，支持反向互补。只有原有的高一致性、额外信息量和完整密码子包含条件同时满足时才剔除弱模型；不同碱基、重复位置和近似并列仍保留歧义。该处理不证明同一拷贝的物理连接。

reads 支持检查允许 read 含 N，但每个 N 消耗一个不匹配，不能删除或充当通配符。仍需完整比对、唯一定位以及与次优位置的差距；短 contig 在 read 中的精确包含仅提供局部序列支持，不作为双端连接证据。

缓存重放可加 `--reads-root output`，从原有成对 reads 重算拼接链证据，在 `chain_audit/<sample>/<SCO>.tsv` 输出有无歧义、搜索限制和双端支持。该诊断不会把候选链直接写进 SCO 序列，也不会填补未知碱基。常规流程先规划链再做全局模型竞争；对已经筛选过的模型重复规划时，会在单一参考内重新计算竞争状态。

蛋白比对未包含终止密码子时，SCO 汇总可从同一原始 contig 紧邻的同框三碱基恢复样本终止密码子：末端必须有唯一、完全一致且含至少五种氨基酸的 10 aa 锚点，投影连续，没有末端编码缺口。不会向下游搜索、跨 contig 寻找终止密码子或复制参考序列；不同实测终止密码子的分歧仍为 N。只有参考 CDS 本来包含终止密码子时才输出这些坐标。

对原先未通过普通蛋白筛选的 30–59 nt 实测候选，若至少两条坐标去重、唯一定位的 reads 合计覆盖整条 contig，可用含至少五种氨基酸的精确 8 aa 锚点补充筛选。重复放置或方向/阅读框不一致会拒绝；仍需运行所选后端并通过原有模型和参考投影检查。没有 reads 时不启用这条补救路径，短 contig 覆盖不等于跨 contig 配对支持。该机制不推断不存在的微外显子。

拼接图的路径筛选在确认超过保留上限时停止比较剩余路径；因候选已排序且保留集合只追加，该优化保持前 64 条路径、搜索受限标记和歧义判定不变，不通过放宽冲突或剪接条件来提速。
