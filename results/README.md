# Benchmark Results

This directory stores benchmark outputs:
- results/builds
- results/classify
- results/profile
- results/reports

Metrics are computed at ranks: species and genus.

## 软件版本

- ganon2: 2.1.0
- sylph: 0.8.1

## 结果分区说明

- **classify**：只展示具备 per-read 指标的工具（例如 ganon2/Chimera）。像 sylph 这类仅输出丰度的工具不会出现在 classify 表中。
- **profile**：展示丰度/检出类指标（abundance/presence）。像 sylph 这类 profile-only 工具只出现在 profile 表中。

## 工具类型与指标适配

工具大致分三类（示例）：
- **只输出丰度（profile‑only）**：只能做 abundance/presence 指标，不产生 per-read 指标。示例：**sylph**。
- **只输出 per-read（per-read‑only）**：有 per-read 指标，但没有直接的 abundance 输出。
- **两者都有（dual）**：同时产出 per-read 与 abundance。示例：**ganon2**、**Chimera**。

注：本项目配置文件里实验名仍可能写作 `ganon`，但实际运行的软件为 **ganon2**，结果表中统一记作 **ganon2**。

对 per-read‑only 工具，理论上可以把每条 read 的预测按 taxon 汇总成“伪丰度”再计算 abundance 指标，
但由于读长/覆盖/多重比对等因素，这种转换往往偏差较大，结果**大概率不够好**。
因此仅建议作为**辅助诊断**使用，正式对比优先使用工具原生输出的 abundance/profile。

## 评估口径（descendant-aware / exact）

本项目的主要评估口径对齐 `/mnt/sda/tianqinzhong/code/Chimera/bench`：
- **descendant-aware（desc）**：允许预测落在真值分类单元的后代（更细粒度）时仍计为正确；用于主表展示。
- **exact**：不做 descendant 折叠（或更严格的 rank 对齐），作为辅助诊断指标，保留在 `metrics.json` 的 `exact_` 前缀字段中。

exact 与 desc 的差异在 abundance/profile 中最明显；per-read（rank 对齐）中 desc 与 exact 往往趋同。

## Per-read metrics（rank: species/genus）

Definitions:
- N = total reads
- g_r = gold taxon lifted to target rank r（若无法 lift 到该 rank，则该 read 在 rank r 上不计入评估）
- p_r = predicted taxon lifted to target rank r（若无法 lift 到该 rank，视为错误预测）
- “unclassified”/缺失 视为未分类（no call）

Counts（per rank r，对齐 bench 的 rank-aware 口径）：
- TP: p_r == g_r
- FN: p 缺失/未分类
- FP+FN: p 有预测但 p_r 为空（无法 lift 到 rank r）或 p_r != g_r（错误预测）

Formulas:
- classified_rate = (# p != unclassified) / N
- unclassified_rate = (# p == unclassified) / N
- Precision (P) = TP / (TP + FP)
- Recall (R) = TP / (TP + FN)
- F1 = 2 * P * R / (P + R)

Notes:
- per-read 指标仅对 **truth 可 lift 到目标 rank** 的 reads 计入该 rank 的评估（即 support 过滤）。
- `exact_` 前缀字段为辅助诊断口径（保持与 bench 一致的 FP/FN 处理）。

## Abundance metrics（rank: species/genus）

Definitions:
- Taxa union i = 1..k at rank r
- Truth vector t_i, prediction vector p_i
- Inputs may be counts or fractions; we renormalize to percentages so **sum_i t_i = 100** and **sum_i p_i = 100**
  - This matches the convention used in `/mnt/sda/tianqinzhong/code/Chimera/bench`
- Presence threshold tau, default tau = 0 (tau=0 means "present if > 0")

Presence counts:
- TP: t_i > tau and p_i > tau
- FP: t_i <= tau and p_i > tau
- FN: t_i > tau and p_i <= tau

Presence formulas:
- Precision (P) = TP / (TP + FP)
- Recall (R) = TP / (TP + FN)
- F1 = 2 * P * R / (P + R)

Distance metrics:
- L1 distance (pct points) = sum_i |p_i - t_i|  (range: 0..200)
- Bray-Curtis (BC) = L1 / 200
- Total variation distance (TV) = L1 / 200 (same as BC under the % normalization)

Notes:
- descendant-aware（desc）会先将预测折叠到 truth taxon（最近真值祖先）再评估（bench 的 `collapse_to_truth_ancestors` 思路）。
- `exact_` 前缀字段为不做 descendant 折叠的辅助诊断口径。

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
- `Precision/Recall/F1 (rank)`：rank-aware per-read 指标（对齐 bench 的 lift-to-rank 口径）。

### Abundance 表（profile）
- `Presence Precision/Recall/F1 (rank)`：descendant-aware presence 指标。
- `L1/TV/Bray-Curtis (rank)`：descendant-aware 丰度距离指标（其中 `L1` 为 percent-points 距离 0..200，`TV/BC` 为 0..1）。

### 额外统计（metrics.json 中可见）
- `truth_profile_species_total / mapped / unmapped / unmapped_rate`：truth_profile 的物种条目统计。
- `truth_profile_mass_total / mapped / unmapped / mapped_rate`：truth_profile 的丰度质量覆盖。
- `truth_mapped_mass_{rank} / truth_unmapped_mass_{rank}`：按 rank 的真值质量覆盖。
- `pred_mapped_mass_{rank} / pred_unmapped_mass_{rank}`：按 rank 的预测质量覆盖。

### 备注
- truth_profile 若为百分比（总和 > 1.5），会自动除以 100 再参与计算。
