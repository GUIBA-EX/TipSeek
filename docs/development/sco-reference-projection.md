# HX50：参考选择和跨参考证据投影修复

已复用两后端全库保存的比对结果修复 SCO 汇总。数据来自 HX50 真实双端 reads 的既有组装候选，参考为 `/home/yasuda2/starsco/geneminer_sco_refs/fasta`，共 7,006 个 SCO。未重新招募 reads、组装或运行比对。所有有候选的 6,513 个 SCO 均已重建；493 个 SCO 在原始候选输入中不存在，仍无输出。

## 参考原始代码后的改动

- Captus `captus/bioformats.py` 的蛋白模型处理使用覆盖长度、identity 与结构质量加权选择参考；HybPiper `hybpiper/distribute_targets.py` 根据命中支持选择每个基因的目标参考。TipSeek 改为 identity 平方加权的唯一残基支持，并限制异常参考长度的影响，替代按蛋白 ID 字典序选择首条参考。
- HybPiper `hybpiper/exonerate_hits.py::_remove_subsumed_hits` 去除包含命中，后续处理修剪重叠并对 query 缺口补 N。这里在选择参考时将重叠/重复残基支持折叠为一次，并按原始 contig、链方向、三个原始碱基坐标去重实际密码子。按用户要求保留矛盾证据为 N；不会用最高分命中覆盖冲突，也不会复制参考碱基。
- 整个样本组内一个 SCO 使用同一参考坐标。多样本注释合并后重新统一 SCO 汇总，串行与并行保持一致。保存的恢复清单可冻结每个 SCO 的参考，供后续重建复用。
- 其他参考通过蛋白全局坐标映射投影，只有所有最低代价路径均同意、局部锚点足够的位置可用。参考插入、重复区位移多解、弱同源和资源上限失败均保留未映射。原始样本密码子若对应多个参考位置，这些位置全部用 N。
- `gene_workflow recover-sco` 从原始比对重建逐残基到 CDS 的映射。归一化模型 TSV 没有完整映射，不能用区间跨度乘三替代。读取原注释对应的候选 DNA、ID 清单和冻结的参考蛋白；缺失必要缓存时明确报错。新注释也保存外部蛋白参考快照。

这不是直接照搬两个软件的流程：保留了用户要求的固定参考顺序、每样本每 SCO 一条序列、缺失和冲突填 N。未进行推测性的微外显子补齐。

## 全库结果

所有阈值按非 N 样本碱基占所选参考编码坐标的比例计算，分母均为 7,006 个参考家族。计数表示有证据的坐标，不表示已验证的完整 CDS、单一拷贝或跨 contig 连接。

| 指标 | miniprot 原结果 | miniprot 修复后 | exonerate-rs 原结果 | exonerate-rs 修复后 |
| --- | ---: | ---: | ---: | ---: |
| 有样本碱基的 SCO | 4,905 | 6,330 | 5,193 | 6,430 |
| 至少 50% 坐标非 N | 956 | 2,211 | 962 | 2,234 |
| 至少 80% 坐标非 N | 194 | 499 | 183 | 465 |
| 至少 95% 坐标非 N | 60 | 128 | 61 | 108 |
| 整条汇总序列无 N | 2 | 3 | 0 | 0 |
| 输出存在但全 N | 1,608 | 183 | 1,320 | 83 |
| 无候选、无输出 | 493 | 493 | 493 | 493 |
| 实测非 N 碱基数 | 2,972,111 | 5,380,936 | 3,058,769 | 5,420,951 |

参考选择是主要瓶颈。仅修复选择、不做跨参考投影时，≥80% 的 SCO 为 510 / 477；投影后为 499 / 465，因为新增证据也揭示了竞争坐标与碱基冲突。不能将更多命中一律算作更多确定序列。

为排除参考切换造成的分母变化，另将每个 SCO 固定为原参考：miniprot ≥80% 从 194 增至 307，exonerate-rs 从 183 增至 291。对应实测碱基分别从 2,972,111 增至 4,527,989、从 3,058,769 增至 4,584,852。该对照逐条验证参考 ID 和序列长度与原结果相同。

## 长内含子面板与局限

12 个既定 SCO、185 个 ppect 外显子坐标区间，固定 ppect 参考；两后端在 PATH 中没有可执行比对程序时仍能重建成功。

| 指标 | miniprot 原 / 修复 | exonerate-rs 原 / 修复 |
| --- | ---: | ---: |
| 有样本碱基的外显子区间 | 140 / 160 | 148 / 162 |
| 整个区间无 N | 18 / 18 | 15 / 12 |
| 实测碱基 / 32,298 参考碱基 | 21,603 / 24,835 | 22,753 / 24,509 |

exonerate-rs 的完整区间数下降是新冲突被掩码后的结果。不能借删除竞争证据来报告虚假的完整度。参考来自其他物种，不是真实 HX50 外显子边界标注；也没有同一个体、同一拷贝的连接真值，不能报告生物学错误合并率。本阶段没有重新改善 reads 招募或建立新的跨 contig reads 连接。全部外显子完整恢复仍未解决。

## 性能与验证

在相同 HX50 miniprot 全库候选和保存比对上，16 个工作线程、顺序进行的暖缓存测试：初版投影实现 21.04 秒、峰值 RSS 1,085,528 KiB；优化后 18.28 秒、876,356 KiB，分别减少 13.1% 和 19.3%。通过窄带起步、按代价上界扩展、避免克隆全组模型及密码子键中的候选字符串实现；扩大计算范围的条件保证坐标歧义判断不变。更合理的资源分配额外恢复了 1,036 个可映射碱基。

实际全库重建用时：miniprot 约 18 秒，exonerate-rs 约 19 秒。上述性能是保存比对后的 SCO 重建，不是 reads 到组装的整条流程提速。原全库注释耗时分别为 1,192.70 与 4,294.77 秒，复用比对避免再次支付该成本。

341 项工作区测试通过；全目标 Clippy（警告作为错误）、发布构建、格式和差异检查通过。覆盖支持驱动选择、重复支持去重、异参考密码子投影、参考插入、同一密码子竞争位置掩码、负链、剪接处分裂密码子、多样本统一坐标、串并行一致、缺失缓存报错，以及窄带与完整矩阵的最优路径对照。逐条校验全部 SCO 的样本记录数、参考 ID、长度、字符集和实测/缺失/冲突计数。

## 文件与复现

- 完整记录：[sco-reference-projection.json](sco-reference-projection.json)。两个后端的独立审计为 `sco-reference-projection-miniprot.json`、`sco-reference-projection-exonerate.json`。
- 最终序列：`/home/yasuda2/tipseek_HX50_sco_projection_20261011/miniprot-optimized/sco/` 和 `exonerate-rs-optimized/sco/`。
- 对照：同目录的 `miniprot-fixed-optimized/`、`exonerate-rs-fixed-optimized/`，保留原参考 ID。
- `manifest/sco_recovery.tsv` 给出参考、覆盖计数、跨参考模型、竞争位置和选择原因；`reference_projection.tsv` 给出每个源参考的映射状态；`sco_selection_control.tsv` 为选择修复单独对照；`sco_replay.tsv` 记录缓存解析与原失败/预筛选排除数量。

```bash
gene_workflow recover-sco \
  --annotation /home/yasuda2/tipseek_HX50_all7006_familyqueue_20261011/miniprot \
  --input /home/yasuda2/tipseek_HX50_all7006_20261011/workflow/gene \
  --out NEW_RECOVERY_DIRECTORY --threads 16
```

固定原参考时追加 `--sco-backbone-manifest ORIGINAL_ANNOTATION/manifest/sco_recovery.tsv`。完整命令、冻结的二进制校验值和时间记录见数据目录的 `optimized-status.json`；审计可使用 `scripts/report_sco_projection.py` 重算。未推送或发布本阶段代码。
