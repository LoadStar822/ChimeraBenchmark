# Benchmark 结果说明

本目录存放 benchmark 的输出与汇总：
- `results/builds`：数据库构建（build）
- `results/classify`：分类/分箱（per-read）
- `results/profile`：丰度画像（profile）
- `results/reports`：额外报告

默认评估 rank：`species`、`genus`。

## 默认参数公平性收口状态

- 公开 benchmark 口径已经收紧为：除统一 `threads=32` 外，推断参数保持软件默认。
- 公开评估统一使用官方 NCBI taxonomy 快照 `results/taxonomy/ncbi_20260408/`。
- `ganon` 结果已经按新的默认参数合同刷新。
- `bracken` 结果已经按新的端到端 `kraken2 -> bracken` 合同刷新。
- 当前仓库里已落盘的 `taxor` 结果，仍是收口前生成的历史结果；刷新后才能作为论文主表结果引用。
- `results/classify/README.md` 与 `results/profile/README.md` 当前仍会展示 `taxor` 的历史结果，因为它们是从现有 run 目录自动汇总出来的。

## 软件版本

- kraken2: 2.1.3
- bracken: 2.9
- centrifuger: 1.1.0-r291
- ganon2: 2.1.0
- sylph: 0.8.1
- taxor: 0.1.3（SeqAn 3.4.0-rc.1）

## 结果分区说明

- `classify`：只展示具备 **per-read** 指标的工具（例如 ganon2/Chimera/Taxor）。像 sylph 这类仅输出丰度的工具不会出现在 classify 表中。
- `profile`：展示 **profile** 指标；profile-only 工具只出现在该表中。

注：本项目配置文件里实验名仍可能写作 `ganon`，但实际运行的软件为 **ganon2**，结果表中统一记作 **ganon2**。

## Bracken 说明（默认参数 + 端到端 pipeline）

Bracken 仅做 abundance 重估（不输出 per-read 分类），因此只出现在 `results/profile/README.md`。

为保证“默认参数”的可比性：
- Bracken 运行使用其默认参数（read_len=100、level=S、threshold=10；仅线程数用于加速）。
- 公开 benchmark 口径按普通学术用户的使用方式执行完整 `kraken2 -> bracken` pipeline，不再复用外部已有 `kraken2.report`。
- 因此 Bracken 的 runtime / memory 必须解释为端到端 pipeline 成本，而不是 Bracken 单步成本。

Bracken 需要在 Kraken2 DB 内生成 `database100mers.kmer_distrib`（通过 `kmer2read_distr` + `generate_kmer_distribution.py`，参数与 bracken 默认一致：k=35、read_len=100）。
由于我们通常会在 Kraken2 build 完成后清理 DB 中的 fasta 库文件以节省空间，Bracken build 时会临时从 CAMIRefseq 的 `target.tsv` 重新生成 `library/added/library.fna`，生成完分布文件后再删除该 fasta。

## Centrifuger 说明（默认参数 + CAMI profile 输出）

- classify：使用 `centrifuger` 默认参数（`-k 1` 等），仅指定 `-t 32` 与 `-x`（DB 前缀）以及输入 reads。
- profile：用 `centrifuger-quant` 基于 classify 的输出生成 profile；为便于评估，使用 `--output-format 2` 让其输出 CAMI 格式（不改变算法与数值，仅改变输出格式）。

## Taxor 说明（默认参数 + 无 benchmark 特殊修补）

- 公开 benchmark 口径不再对 Taxor 做工具级静默跳过。
- 公开 benchmark 口径不再对 Taxor 的中间输出执行 `fix_search` 修补。
- 若 Taxor 在默认流程下因资源或输出格式问题失败，应如实记录为失败或排除，而不是在 benchmark 层静默改写。

## 公开评估口径

本仓库对外只使用一套主口径：
- 所有公开结果都按目标层级判断
- 默认层级为 `species`、`genus`

- `per-read`：看每条 read 或 contig 在指定层级上有没有分对
- `profile`：看整个样本在指定层级上的组成，以及整棵 taxonomy tree 上的距离像不像真值

对 `species/genus` 的判断规则是：
- 先把真值 taxid 沿 taxonomy 向上提升到目标层级
- 再把预测 taxid 也提升到同一层级
- 两者相同就记为正确，不同就记为错误
- 这里看的是“能否提升到目标层级”，不是“原始 rank 是否恰好等于目标层级”

例子：
- 真值是某个 `species`
- 预测到了这个 `species` 下的某个 strain
- 在 `species` 列里，这条预测记为正确
- 真值是某个 `strain/subspecies`
- 只要它能上提到某个 `species`，就进入 `species` 分母
- 真值若只有 `family/order/class`，则不进入 `species/genus` 分母

## Per-read 指标

每个层级都只统计能提升到该层级的 truth 条目。

主表解读规则：
- `F1` 不能脱离 `Truth Mapped Rate / Pred Mapped Rate` 单独引用
- 若某数据集的 mapped rate 明显低于 `1`，说明 taxonomy 映射仍有覆盖缺口，headline 只能解释为“可映射子集上的结果”

计数规则：
- `FN`：没有预测，或预测为空
- `TP`：真值和预测提升到该层级后相同
- `FP + FN`：有预测，但提升到该层级后不相同

公式：
- `Precision = TP / (TP + FP)`
- `Recall = TP / (TP + FN)`
- `F1 = 2 * Precision * Recall / (Precision + Recall)`

补充字段：
- `Classified Rate`：有预测的 read 占比
- `Unclassified Rate`：没有预测的 read 占比
- `Truth Mapped Rate`：真值可提升到该层级的占比；例如 `strain/subspecies` 若能上提到该 rank，则记入 mapped
- `Pred Mapped Rate`：预测可提升到该层级的占比

## Profile 指标

Profile 主表使用 OPAL core：
- `Completeness`：在目标层级上，真值里出现的 taxa 有多少被预测到了
- `Purity`：在目标层级上，预测到的 taxa 里有多少真的在真值里
- `L1 Norm`：把 truth 和 prediction 都折叠到目标层级、再归一化到总和为 `1` 后的 `L1` 距离，范围 `0..2`
- `Weighted UniFrac`：沿 taxonomy tree 比较两边的累积质量差异；本仓库使用 OPAL 默认的树边长度 `1 / depth`，报告未归一化版本

公开结果只统计工具原生输出的 profile 文件：
- 有原生 profile 输出：直接按该输出计算
- 没有原生 profile 输出：不出现在 `profile/README.md`，也不会用 classify 结果回推 abundance

按 `species/genus` 的处理方式：
- 先把 truth 和 prediction 都提升并汇总到目标层级
- `Completeness` 和 `Purity` 只按该层级是否出现来计数
- `L1 Norm` 按该层级的相对丰度分布计数

计数规则：
- `TP`：该层级的 taxon 同时出现在 truth 和 prediction 中
- `FP`：该层级的 taxon 只出现在 prediction 中
- `FN`：该层级的 taxon 只出现在 truth 中

公式：
- `Completeness = TP / (TP + FN)`
- `Purity = TP / (TP + FP)`
- `L1 Norm = sum_i |p_i - t_i|`

解读方式：
- `Completeness`、`Purity` 越大越好
- `L1 Norm`、`Weighted UniFrac` 越小越好

## 结果表字段说明

通用资源字段：
- `Elapsed (s)`：一次 run 的 wall time（秒），来自 `meta.json`
- `Max RSS (GB)`：最大常驻内存（GB），来自 `/usr/bin/time -v`

结果表：
- `classify/README.md`：展示 per-read 的公开主口径结果
- `profile/README.md`：展示 profile 的公开主口径结果

原始 JSON：
- `metrics.json` 里可能保留一些历史辅助审计字段
- benchmark README、对外结果表和论文主表默认不使用这些辅助字段

备注：
- truth_profile 若为百分比（总和 > 1.5），会自动除以 100 再参与计算。
