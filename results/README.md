# Benchmark Results

This directory stores benchmark outputs:
- results/builds
- results/classify
- results/profile
- results/reports

Metrics are computed at ranks: species and genus.

## 软件版本

- ganon: 2.1.0
- sylph: 0.8.1

## 结果分区说明

- **classify**：只展示具备 per-read 指标的工具（例如 ganon/Chimera）。像 sylph 这类仅输出丰度的工具不会出现在 classify 表中。当前表格默认展示 **UNK 指标**。
- **profile**：展示丰度/检出类指标。当前表格默认展示 **UNK 指标**；closed-set 与 exact 口径保留在 `metrics.json` 中用于诊断。

## 工具类型与指标适配

工具大致分三类（示例）：
- **只输出丰度（profile‑only）**：只能做 abundance/presence 指标，不产生 per-read 指标。示例：**sylph**。
- **只输出 per-read（per-read‑only）**：有 per-read 指标，但没有直接的 abundance 输出。
- **两者都有（dual）**：同时产出 per-read 与 abundance。示例：**ganon**、**Chimera**。

对 per-read‑only 工具，理论上可以把每条 read 的预测按 taxon 汇总成“伪丰度”再计算 abundance 指标，
但由于读长/覆盖/多重比对等因素，这种转换往往偏差较大，结果**大概率不够好**。
因此仅建议作为**辅助诊断**使用，正式对比优先使用工具原生输出的 abundance/profile。

## 评估口径与 OOD（UNK / closed-set / exact）

默认口径为 **descendant-aware**：
- per-read：若预测 taxid 是真值 taxid（在该 rank 的 taxon）**的后代**（含自身），则计为 TP。
- abundance：先把预测分配到**最近的真值祖先**（truth 中出现的 taxon），再计算 presence 与距离指标。

同时会在 `metrics.json` 中保留 **exact**（严格等价）指标，字段前缀为 `exact_`。

**UNK 指标（后缀 `_unk`）**：
把映射失败的质量统一并入一个 `__UNMAPPED__` bucket（taxid = -1），用于显式反映 OOD/覆盖缺口。
结果表默认展示 UNK 指标，便于快速诊断 coverage 问题。

Coverage/映射率（按 rank）：
- `truth_mapped_mass_{rank}` / `truth_unmapped_mass_{rank}`
- `truth_mapped_mass_rate_{rank} = mapped / (mapped + unmapped)`
- `pred_mapped_mass_{rank}` / `pred_unmapped_mass_{rank}`
- `pred_mapped_mass_rate_{rank} = mapped / (mapped + unmapped)`

Per-read 覆盖率（按 rank）：
- `per_read_truth_mapped_rate_{rank} = (# true_rank 可映射) / N`
- `per_read_pred_mapped_rate_{rank} = (# pred_taxid 非空) / N`

备注：
- 若在实验配置中提供 `coverage_target_tsv` + `coverage_nodes_dmp`，
  覆盖判定会基于 **target.tsv 中的 taxid 及其在 nodes.dmp 上溯后的 rank 集合**；
  用于保证不同工具（同一 DB）覆盖率口径一致。

## Per-read metrics (species/genus)

Definitions:
- N = total reads
- g = gold taxon at target rank
- p = predicted taxon id (any rank)
- p is unclassified if missing or labeled "unclassified"

Counts (per rank r, descendant-aware):
- TP: p is descendant of g (including g)
- FP: p is classified and not descendant of g
- FN: p is unclassified or not descendant of g

Exact counts (strict, `exact_`):
- TP: p mapped to rank r equals g
- FP: p mapped to rank r != g and p is classified
- FN: p unclassified or p mapped to rank r != g

Formulas:
- classified_rate = (# p != unclassified) / N
- unclassified_rate = (# p == unclassified) / N
- Precision (P) = TP / (TP + FP)
- Recall (R) = TP / (TP + FN)
- F1 = 2 * P * R / (P + R)

Notes:
- closed-set per-read 指标仅对 **truth 可映射** 的 reads 计算；
- `_unk` 后缀的 per-read 指标会把无法映射的 truth/pred 统一为 UNK 参与计算；
- `exact_` 指标为严格等价口径，用于诊断误差来源。

## Abundance metrics (species/genus)

Definitions:
- Taxa union i = 1..k at rank r
- Truth vector t_i, prediction vector p_i
- If inputs are counts, normalize so sum_i t_i = 1 and sum_i p_i = 1
- Presence threshold tau, default tau = 0

Presence counts:
- TP: t_i > tau and p_i > tau
- FP: t_i <= tau and p_i > tau
- FN: t_i > tau and p_i <= tau

Presence formulas:
- Precision (P) = TP / (TP + FP)
- Recall (R) = TP / (TP + FN)
- F1 = 2 * P * R / (P + R)

Distance metrics:
- L1 = sum_i |p_i - t_i|
- TV = 0.5 * L1
- Bray-Curtis (BC) = sum_i |p_i - t_i| / (sum_i p_i + sum_i t_i)
  - If normalized (sum = 1), BC = 0.5 * L1

Notes:
- descendant-aware 指标会先将预测折叠到 truth taxon（最近真值祖先）再评估；
- closed-set 指标只在可映射 taxa 上评估；
- `_unk` 后缀的 abundance/presence 指标会把未映射质量并入 UNK bucket；
- `exact_` 指标为严格等价口径，用于诊断误差来源。

## 指标字段说明（按结果表列名）

### 通用资源指标
- `Elapsed (s)`：该次运行的 wall time（秒），取自 meta.json。
- `Max RSS (GB)`：该次运行的最大常驻内存（GB），取自 `/usr/bin/time -v`。

### Per-read 表（classify）
- `Total Reads`：总 reads 数量。
- `Classified Reads`：被分配到 taxon 的 reads 数量。
- `Unclassified Reads`：未分配（unclassified）的 reads 数量。
- `Classified Rate`：`Classified Reads / Total Reads`。
- `Unclassified Rate`：`Unclassified Reads / Total Reads`。
- `Truth Mapped Rate (rank)`：真值在该 rank 可映射的 reads 比例。
- `Pred Mapped Rate (rank)`：预测在该 rank 可映射的 reads 比例。
- `Precision/Recall/F1 (rank)`：descendant-aware per-read 指标（只在 truth 可映射 reads 上计算）。
- `Precision/Recall/F1 (rank, UNK)`：UNK 版本 per-read 指标（将无法映射的 truth/pred 统一视为 UNK 参与计算）。

### Abundance 表（profile）
- `Truth Mapped Rate (rank)`：真值在该 rank 可映射的丰度质量比例。
- `Pred Mapped Rate (rank)`：预测在该 rank 可映射的丰度质量比例。
- `Presence Precision/Recall/F1 (rank)`：descendant-aware presence 指标。
- `L1/TV/Bray-Curtis (rank)`：descendant-aware 丰度距离指标。
- `Presence Precision/Recall/F1 (rank, UNK)`：UNK 版本 presence 指标。
- `L1/TV/Bray-Curtis (rank, UNK)`：UNK 版本丰度距离指标。

### 额外统计（metrics.json 中可见）
- `truth_profile_species_total / mapped / unmapped / unmapped_rate`：truth_profile 的物种条目统计。
- `truth_profile_mass_total / mapped / unmapped / mapped_rate`：truth_profile 的丰度质量覆盖。
- `truth_mapped_mass_{rank} / truth_unmapped_mass_{rank}`：按 rank 的真值质量覆盖。
- `pred_mapped_mass_{rank} / pred_unmapped_mass_{rank}`：按 rank 的预测质量覆盖。

### 备注
- truth_profile 若为百分比（总和 > 1.5），会自动除以 100 再参与计算。
