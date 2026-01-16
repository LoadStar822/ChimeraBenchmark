# Benchmark 结果说明

本目录存放 benchmark 的输出与汇总：
- `results/builds`：数据库构建（build）
- `results/classify`：分类/分箱（per-read）
- `results/profile`：丰度/检出（abundance/presence）
- `results/reports`：额外报告

默认评估 rank：`species`、`genus`。

## 软件版本

- ganon2: 2.1.0
- sylph: 0.8.1
- taxor: 0.1.3（SeqAn 3.4.0-rc.1）

## 结果分区说明

- `classify`：只展示具备 **per-read** 指标的工具（例如 ganon2/Chimera/Taxor）。像 sylph 这类仅输出丰度的工具不会出现在 classify 表中。
- `profile`：展示 **abundance/presence** 指标；profile-only 工具只出现在该表中。

注：本项目配置文件里实验名仍可能写作 `ganon`，但实际运行的软件为 **ganon2**，结果表中统一记作 **ganon2**。

## Taxor 运行限制（资源）

Taxor 在 read 型数据集上生成极大的 `taxor_search.tsv` 中间文件，并触发高内存/磁盘压力（甚至被系统 kill），影响服务器稳定性。

因此目前 Taxor 结果仅在 CAMI contigs 数据集上报告：
- `cami2-marine-long-sample0`
- `cami2-marine-long`（包含 sample_0..9 的 contigs，合并为一次 run）
- `cami2-marine-short`（包含 sample_0..9 的 contigs，合并为一次 run）

由于资源限制，以下数据集暂不运行 Taxor：
- `atcc-illumina`
- `atcc-hifi`
- `zymo-gridion-even`
- `zymo-gridion-log`
- `zymo-promethion-even`
- `zymo-promethion-log`

为保护服务器，同时不将 `taxor_search.tsv` 作为长期结果保留（需要时请重跑 search 步骤）。

## 评估口径（desc / exact）

本项目的指标口径对齐 `/mnt/sda/tianqinzhong/code/Chimera/bench`：
- **descendant-aware（desc）**：允许预测落在真值分类单元的后代（更细粒度）时仍计为正确；用于结果表（默认展示）。
- **exact**：不做 descendant 折叠/放宽（或更严格的 rank 对齐），作为辅助诊断，保留在 `metrics.json` 的 `exact_` 前缀字段中。

## Per-read 指标（rank: species/genus）

记号：
- `N`：truth reads 总数（分母始终使用 `N`，便于不同 rank 横向比较）
- `g`：truth taxid
- `p`：pred taxid（`unclassified`/缺失视为 no-call，记作 `p = None`）
- `lift(x, r)`：将 taxid `x` 沿 taxonomy 向上提升到 rank `r`（若不存在则返回 `None`）

对每个 rank `r`：
- 仅对 **truth 可提升到该 rank** 的 reads 计入该 rank 的评估（即 support 过滤）：
  - `g_r = lift(g, r)`，若 `g_r` 为 `None`，则该 read 在 rank `r` 上跳过

计数（对齐 bench 的 FP/FN 处理方式）：
- `FN`：`p = None`
- `TP`（desc）：`p` 在 taxonomy 上是 `g_r` 的后代（等价于 `g_r ∈ lineage(p)`）
- `TP`（exact）：`lift(p, r) == g_r`
- `FP + FN`：`p` 有预测但不满足 TP 条件（错误预测记作同时贡献一次 FP 与一次 FN）

公式：
- `classified_rate = (# p != None) / N`
- `unclassified_rate = (# p == None) / N`
- `TruthMappedRate(r) = (# support reads at rank r) / N`
- `PredMappedRate(r) = (# support reads with p != None) / N`
- `Precision = TP / (TP + FP)`
- `Recall = TP / (TP + FN)`
- `F1 = 2 * Precision * Recall / (Precision + Recall)`

说明：
- 结果表中的 `Precision/Recall/F1` 对应 **desc** 口径；`metrics.json` 中的 `exact_` 前缀字段对应 exact 口径。

## Abundance 指标（rank: species/genus）

输入：
- truth 向量 `t_i` 与预测向量 `p_i`（可为 counts 或 fractions）
- 我们会统一重归一化到百分比，使得：`sum_i t_i = 100`、`sum_i p_i = 100`（对齐 bench）
- Presence 阈值 `tau`（默认 `tau = 0`，表示 “>0 即认为 present”）

Presence 计数（对每个 taxon i）：
- `TP`：`t_i > tau` 且 `p_i > tau`
- `FP`：`t_i <= tau` 且 `p_i > tau`
- `FN`：`t_i > tau` 且 `p_i <= tau`

Presence 指标：
- `Precision = TP / (TP + FP)`
- `Recall = TP / (TP + FN)`
- `F1 = 2 * Precision * Recall / (Precision + Recall)`

距离指标（百分比点）：
- `L1 = sum_i |p_i - t_i|`（范围 `0..200`）
- `TV = L1 / 200`
- `Bray-Curtis = L1 / 200`（在该归一化下与 TV 相同）

desc/exact 说明：
- **desc**：会先将预测折叠到 truth taxon（最近真值祖先）再评估（对齐 bench 的 `collapse_to_truth_ancestors` 思路）。
- **exact**：不做上述折叠，直接对齐到 rank 后评估（保存在 `exact_` 前缀字段）。

## 结果表字段说明

通用资源字段：
- `Elapsed (s)`：一次 run 的 wall time（秒），来自 `meta.json`
- `Max RSS (GB)`：最大常驻内存（GB），来自 `/usr/bin/time -v`

表格字段：
- `classify/README.md`：展示 per-read 的 `Precision/Recall/F1 (species/genus)`（desc 口径）
- `profile/README.md`：展示 abundance 的 `Presence*` 与 `L1/TV/Bray-Curtis (species/genus)`（desc 口径）

metrics.json 常见字段：
- `per_read_*`：per-read 指标（desc）；`exact_per_read_*`：exact 版本
- `abundance_*`、`presence_*`：abundance/presence 指标（desc）；`exact_abundance_*`、`exact_presence_*`：exact 版本
- `truth_profile_*`：truth_profile 的条目数/质量覆盖统计（用于诊断 truth 与 taxonomy 映射是否完备）

备注：
- truth_profile 若为百分比（总和 > 1.5），会自动除以 100 再参与计算。
