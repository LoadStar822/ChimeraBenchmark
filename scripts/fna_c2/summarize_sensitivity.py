#!/usr/bin/env python3
from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import random
import re
import statistics
from collections import defaultdict
from pathlib import Path


DEFAULT_THRESHOLDS = (0, 500, 1000, 2000, 5000)
COUNT_METRICS = {
    "strict_c2_best_reads_per_million": "strict_c2_best_reads",
    "c2_supported_reads_per_million": "c2_supported_reads",
}


def auc_with_ties(positive: list[float], negative: list[float]) -> float:
    if not positive or not negative:
        raise ValueError("AUC requires both positive and negative values")
    ordered = sorted(negative)
    score = 0.0
    for value in positive:
        lower = bisect.bisect_left(ordered, value)
        upper = bisect.bisect_right(ordered, value)
        score += lower + 0.5 * (upper - lower)
    return score / (len(positive) * len(negative))


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("percentile requires values")
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def stratified_bootstrap_auc(
    positive: list[float],
    negative: list[float],
    *,
    replicates: int,
    seed: int,
) -> tuple[float, float]:
    if replicates <= 0:
        raise ValueError("bootstrap replicates must be positive")
    rng = random.Random(seed)
    values: list[float] = []
    for _ in range(replicates):
        sampled_positive = [positive[rng.randrange(len(positive))] for _ in positive]
        sampled_negative = [negative[rng.randrange(len(negative))] for _ in negative]
        values.append(auc_with_ties(sampled_positive, sampled_negative))
    return percentile(values, 0.025), percentile(values, 0.975)


def input_reads(row: dict[str, str]) -> int:
    raw = row.get("input_reads") or row.get("input_reads_used")
    if raw:
        return int(float(raw))
    match = re.search(r"head(\d+)k", row.get("sample", ""))
    if match:
        return int(match.group(1)) * 1000
    raise ValueError(f"missing input read depth for {row.get('sample', '<unknown>')}")


def metric_value(row: dict[str, str], count_field: str) -> float:
    depth = input_reads(row)
    return float(row[count_field]) / depth * 1_000_000.0 if depth else 0.0


def expected_c2_reads(row: dict[str, str]) -> float:
    return float(row["paper_fna_c2_pct"]) / 100.0 * input_reads(row)


def stable_seed(*values: object) -> int:
    digest = hashlib.sha256("\x1f".join(map(str, values)).encode()).digest()
    return int.from_bytes(digest[:8], "big")


def build_auc_rows(
    rows: list[dict[str, str]],
    *,
    tool: str,
    audit_scheme: str,
    thresholds: list[int] | tuple[int, ...],
    bootstrap_replicates: int,
) -> list[dict[str, object]]:
    negatives = [row for row in rows if row["role"] == "paper_zero_fna"]
    output: list[dict[str, object]] = []
    for metric, count_field in COUNT_METRICS.items():
        negative_values = [metric_value(row, count_field) for row in negatives]
        for threshold in thresholds:
            positives = [
                row
                for row in rows
                if row["role"] == "paper_c2_positive"
                and expected_c2_reads(row) >= threshold
            ]
            positive_values = [metric_value(row, count_field) for row in positives]
            auc = auc_with_ties(positive_values, negative_values) if positive_values and negative_values else None
            if auc is not None and bootstrap_replicates:
                ci_low, ci_high = stratified_bootstrap_auc(
                    positive_values,
                    negative_values,
                    replicates=bootstrap_replicates,
                    seed=stable_seed(tool, audit_scheme, metric, threshold),
                )
            else:
                ci_low = ci_high = None
            output.append(
                {
                    "audit_scheme": audit_scheme,
                    "tool": tool,
                    "metric": metric,
                    "threshold_expected_c2_reads": threshold,
                    "n_positive": len(positive_values),
                    "n_negative": len(negative_values),
                    "auc": auc,
                    "auc_ci_low": ci_low,
                    "auc_ci_high": ci_high,
                }
            )
    return output


def read_rows(paths: list[Path]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    seen: set[str] = set()
    for path in paths:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sample = row["sample"]
                if sample in seen:
                    raise ValueError(f"duplicate sample across input tables: {sample}")
                seen.add(sample)
                output.append(row)
    return output


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    field: (
                        ""
                        if row.get(field) is None
                        else f"{row[field]:.6f}"
                        if isinstance(row.get(field), float)
                        else row.get(field, "")
                    )
                    for field in fields
                }
            )


def parse_old_tables(values: list[str]) -> dict[str, list[Path]]:
    output: dict[str, list[Path]] = defaultdict(list)
    for value in values:
        if "=" not in value:
            raise ValueError(f"old table must be TOOL=PATH: {value}")
        tool, raw_path = value.split("=", 1)
        path = Path(raw_path)
        if not path.exists():
            raise FileNotFoundError(path)
        output[tool].append(path)
    return dict(output)


def count_change_rows(
    old_rows: list[dict[str, str]],
    new_rows: list[dict[str, str]],
    *,
    tool: str,
) -> list[dict[str, object]]:
    old_by_sample = {row["sample"]: row for row in old_rows}
    new_by_sample = {row["sample"]: row for row in new_rows}
    if old_by_sample.keys() != new_by_sample.keys():
        raise ValueError(f"old/new sample mismatch for {tool}")
    output: list[dict[str, object]] = []
    for metric, count_field in COUNT_METRICS.items():
        changes = [
            metric_value(new_by_sample[sample], count_field)
            - metric_value(old_by_sample[sample], count_field)
            for sample in old_by_sample
        ]
        absolute = [abs(value) for value in changes]
        output.append(
            {
                "tool": tool,
                "metric": metric,
                "samples": len(changes),
                "median_new_minus_old": statistics.median(changes),
                "median_absolute_change": statistics.median(absolute),
                "max_absolute_change": max(absolute),
            }
        )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize the Fna C2 reference/direct-alignment sensitivity run.")
    parser.add_argument("--new-dir", action="append", required=True, type=Path)
    parser.add_argument("--old-table", action="append", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--bootstrap-replicates", type=int, default=10_000)
    args = parser.parse_args()

    old_paths = parse_old_tables(args.old_table)
    auc_rows: list[dict[str, object]] = []
    change_rows: list[dict[str, object]] = []
    tools = ("chimera", "centrifuger", "kraken2")
    for tool in tools:
        new_rows = read_rows([path / f"{tool}_candidate_identity_summary.tsv" for path in args.new_dir])
        old_rows = read_rows(old_paths[tool])
        selected_samples = {row["sample"] for row in new_rows}
        old_selected = [row for row in old_rows if row["sample"] in selected_samples]
        if len(old_selected) != len(new_rows):
            raise ValueError(f"old table missing selected samples for {tool}")
        auc_rows.extend(
            build_auc_rows(
                old_selected,
                tool=tool,
                audit_scheme="original_combined_panel_candidate",
                thresholds=DEFAULT_THRESHOLDS,
                bootstrap_replicates=args.bootstrap_replicates,
            )
        )
        auc_rows.extend(
            build_auc_rows(
                new_rows,
                tool=tool,
                audit_scheme="exact_dedup_per_clade_candidate",
                thresholds=DEFAULT_THRESHOLDS,
                bootstrap_replicates=args.bootstrap_replicates,
            )
        )
        change_rows.extend(count_change_rows(old_selected, new_rows, tool=tool))

    direct_rows = read_rows([path / "direct_alignment_sample_summary.tsv" for path in args.new_dir])
    auc_rows.extend(
        build_auc_rows(
            direct_rows,
            tool="direct_minimap2",
            audit_scheme="exact_dedup_per_clade_all_reads",
            thresholds=DEFAULT_THRESHOLDS,
            bootstrap_replicates=args.bootstrap_replicates,
        )
    )

    auc_fields = [
        "audit_scheme",
        "tool",
        "metric",
        "threshold_expected_c2_reads",
        "n_positive",
        "n_negative",
        "auc",
        "auc_ci_low",
        "auc_ci_high",
    ]
    write_tsv(args.out_dir / "sensitivity_auc.tsv", auc_rows, auc_fields)
    write_tsv(
        args.out_dir / "reference_count_change.tsv",
        change_rows,
        [
            "tool",
            "metric",
            "samples",
            "median_new_minus_old",
            "median_absolute_change",
            "max_absolute_change",
        ],
    )

    runtime_by_sample: dict[str, float] = defaultdict(float)
    for path in args.new_dir:
        with (path / "runtime.tsv").open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                runtime_by_sample[row["sample"]] += float(row["seconds"])
    runtime_rows = [
        {"sample": sample, "direct_alignment_seconds": seconds}
        for sample, seconds in sorted(runtime_by_sample.items())
    ]
    write_tsv(
        args.out_dir / "runtime_per_sample.tsv",
        runtime_rows,
        ["sample", "direct_alignment_seconds"],
    )
    aggregate = [
        {
            "samples": len(runtime_rows),
            "mean_seconds": statistics.fmean(runtime_by_sample.values()),
            "median_seconds": statistics.median(runtime_by_sample.values()),
            "min_seconds": min(runtime_by_sample.values()),
            "max_seconds": max(runtime_by_sample.values()),
        }
    ]
    write_tsv(
        args.out_dir / "runtime_summary.tsv",
        aggregate,
        ["samples", "mean_seconds", "median_seconds", "min_seconds", "max_seconds"],
    )


if __name__ == "__main__":
    main()
