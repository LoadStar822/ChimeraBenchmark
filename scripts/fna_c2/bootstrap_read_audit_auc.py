#!/usr/bin/env python3
from __future__ import annotations

import argparse
import bisect
import csv
import random
import statistics
from collections import Counter
from pathlib import Path


ITERATIONS = 10_000
SEED = 20_260_710
THRESHOLDS = (0, 500, 1000, 2000, 5000)
COHORTS = ("2015 Yu et al", "2018 Wirbel et al", "YachidaS_2019")
TOOLS = ("Chimera", "Centrifuger", "Kraken2_LF01")
METRICS = {
    "raw_candidate_reads_per_million": "candidate_reads_per_million",
    "audited_c2_supported_reads_per_million": "c2_supported_reads_per_million",
    "audited_strict_c2_best_reads_per_million": "strict_c2_best_reads_per_million",
}
FIELDS = [
    "tool",
    "metric",
    "threshold_expected_c2_reads",
    "n_positive",
    "n_negative",
    "positive_mean",
    "negative_mean",
    "positive_median",
    "negative_median",
    "auc_c2_positive_vs_zero",
    "auc_ci_low",
    "auc_ci_high",
    "chimera_minus_tool_auc",
    "difference_ci_low",
    "difference_ci_high",
    "difference_p_value",
    "bootstrap_iterations",
    "bootstrap_seed",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def auc(positive: list[float], negative: list[float]) -> float:
    if not positive or not negative:
        raise ValueError("AUC requires non-empty positive and negative groups")
    ordered_negative = sorted(negative)
    wins = 0.0
    for value in positive:
        lower = bisect.bisect_left(ordered_negative, value)
        upper = bisect.bisect_right(ordered_negative, value)
        wins += lower + 0.5 * (upper - lower)
    return wins / (len(positive) * len(negative))


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def build_rows(
    manifest_rows: list[dict[str, str]],
    audit_rows: list[dict[str, str]],
    *,
    iterations: int = ITERATIONS,
    seed: int = SEED,
) -> list[dict[str, object]]:
    if iterations <= 0:
        raise ValueError("bootstrap iterations must be positive")
    sample_order = [row["sample"] for row in manifest_rows]
    if len(sample_order) != len(set(sample_order)):
        raise ValueError("sample manifest contains duplicate samples")
    manifest = {row["sample"]: row for row in manifest_rows}
    if {row["cohort"] for row in manifest_rows} != set(COHORTS):
        raise ValueError("sample manifest cohort set differs from the frozen three-cohort design")

    keys = [(row["sample"], row["tool"]) for row in audit_rows]
    duplicates = [key for key, count in Counter(keys).items() if count != 1]
    if duplicates:
        raise ValueError(f"duplicate audit tool/sample rows: {duplicates[:5]}")
    expected_keys = {(sample, tool) for sample in sample_order for tool in TOOLS}
    if set(keys) != expected_keys:
        raise ValueError(
            f"audit matrix mismatch: missing={len(expected_keys - set(keys))} "
            f"extra={len(set(keys) - expected_keys)}"
        )

    sample_info: dict[str, dict[str, str | float]] = {}
    values = {tool: {metric: {} for metric in METRICS} for tool in TOOLS}
    for row in audit_rows:
        sample = row["sample"]
        source = manifest[sample]
        info = {
            "role": source["role"],
            "cohort": source["cohort"],
            "expected": float(source["expected_c2_reads_at_depth"]),
        }
        if row["role"] != info["role"] or row["cohort"] != info["cohort"]:
            raise ValueError(f"audit/manifest metadata mismatch for {sample}/{row['tool']}")
        sample_info[sample] = info
        for metric, field in METRICS.items():
            values[row["tool"]][metric][sample] = float(row[field])

    pools: dict[int, tuple[dict[str, list[str]], dict[str, list[str]]]] = {}
    for threshold in THRESHOLDS:
        positive = {cohort: [] for cohort in COHORTS}
        negative = {cohort: [] for cohort in COHORTS}
        for sample in sample_order:
            info = sample_info[sample]
            cohort = str(info["cohort"])
            if info["role"] == "paper_c2_positive" and float(info["expected"]) >= threshold:
                positive[cohort].append(sample)
            elif info["role"] == "paper_zero_fna":
                negative[cohort].append(sample)
        pools[threshold] = positive, negative

    observed: dict[tuple[str, str, int], float] = {}
    distributions: dict[tuple[str, str, int], list[float]] = {}
    for threshold, (positive, negative) in pools.items():
        positive_ids = [sample for cohort in COHORTS for sample in positive[cohort]]
        negative_ids = [sample for cohort in COHORTS for sample in negative[cohort]]
        for metric in METRICS:
            for tool in TOOLS:
                key = tool, metric, threshold
                observed[key] = auc(
                    [values[tool][metric][sample] for sample in positive_ids],
                    [values[tool][metric][sample] for sample in negative_ids],
                )
                distributions[key] = []

    generator = random.Random(seed)
    for _ in range(iterations):
        sampled_negative = {
            cohort: generator.choices(pools[0][1][cohort], k=len(pools[0][1][cohort]))
            for cohort in COHORTS
        }
        for threshold, (positive, _) in pools.items():
            positive_ids = [
                sample
                for cohort in COHORTS
                for sample in generator.choices(positive[cohort], k=len(positive[cohort]))
            ]
            negative_ids = [sample for cohort in COHORTS for sample in sampled_negative[cohort]]
            for metric in METRICS:
                for tool in TOOLS:
                    distributions[(tool, metric, threshold)].append(
                        auc(
                            [values[tool][metric][sample] for sample in positive_ids],
                            [values[tool][metric][sample] for sample in negative_ids],
                        )
                    )

    output: list[dict[str, object]] = []
    for tool in TOOLS:
        for metric in METRICS:
            for threshold in THRESHOLDS:
                positive_ids = [
                    sample for cohort in COHORTS for sample in pools[threshold][0][cohort]
                ]
                negative_ids = [sample for cohort in COHORTS for sample in pools[0][1][cohort]]
                positive_values = [values[tool][metric][sample] for sample in positive_ids]
                negative_values = [values[tool][metric][sample] for sample in negative_ids]
                key = tool, metric, threshold
                boot = distributions[key]
                row: dict[str, object] = {
                    "tool": tool,
                    "metric": metric,
                    "threshold_expected_c2_reads": threshold,
                    "n_positive": len(positive_values),
                    "n_negative": len(negative_values),
                    "positive_mean": statistics.fmean(positive_values),
                    "negative_mean": statistics.fmean(negative_values),
                    "positive_median": statistics.median(positive_values),
                    "negative_median": statistics.median(negative_values),
                    "auc_c2_positive_vs_zero": observed[key],
                    "auc_ci_low": quantile(boot, 0.025),
                    "auc_ci_high": quantile(boot, 0.975),
                    "chimera_minus_tool_auc": None,
                    "difference_ci_low": None,
                    "difference_ci_high": None,
                    "difference_p_value": None,
                    "bootstrap_iterations": iterations,
                    "bootstrap_seed": seed,
                }
                if tool != "Chimera":
                    chimera_key = "Chimera", metric, threshold
                    differences = [
                        chimera - comparator
                        for chimera, comparator in zip(distributions[chimera_key], boot)
                    ]
                    tail_low = (sum(value <= 0.0 for value in differences) + 1) / (iterations + 1)
                    tail_high = (sum(value >= 0.0 for value in differences) + 1) / (iterations + 1)
                    row.update(
                        {
                            "chimera_minus_tool_auc": observed[chimera_key] - observed[key],
                            "difference_ci_low": quantile(differences, 0.025),
                            "difference_ci_high": quantile(differences, 0.975),
                            "difference_p_value": min(1.0, 2.0 * min(tail_low, tail_high)),
                        }
                    )
                output.append(row)
    return output


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    field: (
                        ""
                        if row.get(field) is None
                        else f"{row[field]:.6g}"
                        if field == "difference_p_value"
                        else f"{row[field]:.6f}"
                        if isinstance(row.get(field), float)
                        else row.get(field, "")
                    )
                    for field in FIELDS
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recompute stratified bootstrap AUCs for the frozen Fna C2 read audit."
    )
    parser.add_argument("--sample-manifest", required=True, type=Path)
    parser.add_argument("--audit-table", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--bootstrap", type=int, default=ITERATIONS)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()

    rows = build_rows(
        read_tsv(args.sample_manifest),
        read_tsv(args.audit_table),
        iterations=args.bootstrap,
        seed=args.seed,
    )
    write_tsv(args.out, rows)
    print(f"wrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
