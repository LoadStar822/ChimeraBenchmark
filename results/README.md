# Benchmark Results

This directory stores benchmark outputs:
- results/builds
- results/classify
- results/profile
- results/reports

Metrics are computed at ranks: species and genus.

## Per-read metrics (species/genus)

Definitions:
- N = total reads
- g = gold taxon at target rank
- p = predicted taxon at target rank
- p is unclassified if missing or labeled "unclassified"

Counts (per rank r):
- TP: p == g
- FP: p != g and p is not unclassified
- FN: p == unclassified or p != g

Formulas:
- classified_rate = (# p != unclassified) / N
- unclassified_rate = (# p == unclassified) / N
- Precision (P) = TP / (TP + FP)
- Recall (R) = TP / (TP + FN)
- F1 = 2 * P * R / (P + R)

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
