#!/usr/bin/env python3
from __future__ import annotations

import argparse
import bisect
import csv
import math
import statistics
from collections import Counter
from pathlib import Path


DEFAULT_THRESHOLDS = (0, 500, 1000, 2000, 5000)
TOOLS = {
    "centrifuger_default": "Centrifuger",
    "chimera_default": "Chimera",
    "kraken2_safe_lf01": "Kraken2_LF01",
    "sylph_default": "sylph",
}
OVERALL_FIELDS = [
    "tool",
    "unit",
    "n",
    "n_c2_positive",
    "n_zero_fna",
    "c2_positive_mean",
    "zero_fna_mean",
    "c2_positive_median",
    "zero_fna_median",
    "positive_minus_zero_mean",
    "auc_c2_positive_vs_zero",
    "pearson_paper_c2",
    "spearman_paper_c2",
]
FLOOR_FIELDS = [
    "tool",
    "unit",
    "positive_rule",
    "expected_c2_reads_min",
    "n_pos",
    "n_zero",
    "pos_mean",
    "zero_mean",
    "pos_median",
    "zero_median",
    "delta_mean",
    "auc_c2_positive_vs_zero",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def auc_with_ties(positive: list[float], negative: list[float]) -> float:
    if not positive or not negative:
        raise ValueError("AUC requires non-empty positive and negative groups")
    ordered_negative = sorted(negative)
    wins = 0.0
    for value in positive:
        lower = bisect.bisect_left(ordered_negative, value)
        upper = bisect.bisect_right(ordered_negative, value)
        wins += lower + 0.5 * (upper - lower)
    return wins / (len(positive) * len(negative))


def _average_ranks(values: list[float]) -> list[float]:
    ordered = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end][1] == ordered[start][1]:
            end += 1
        rank = (start + 1 + end) / 2.0
        for index, _ in ordered[start:end]:
            ranks[index] = rank
        start = end
    return ranks


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        raise ValueError("correlation requires equal vectors with at least two values")
    left_mean = statistics.fmean(left)
    right_mean = statistics.fmean(right)
    numerator = sum(
        (left_value - left_mean) * (right_value - right_mean)
        for left_value, right_value in zip(left, right)
    )
    left_scale = math.sqrt(sum((value - left_mean) ** 2 for value in left))
    right_scale = math.sqrt(sum((value - right_mean) ** 2 for value in right))
    if left_scale == 0.0 or right_scale == 0.0:
        raise ValueError("correlation is undefined for a constant vector")
    return numerator / (left_scale * right_scale)


def spearman(left: list[float], right: list[float]) -> float:
    return pearson(_average_ranks(left), _average_ranks(right))


def validate_inputs(
    manifest_rows: list[dict[str, str]],
    signal_rows: list[dict[str, str]],
) -> dict[str, dict[str, str]]:
    manifest = {row["sample"]: row for row in manifest_rows}
    if len(manifest) != len(manifest_rows):
        raise ValueError("sample manifest contains duplicate samples")
    if len(manifest) != 760:
        raise ValueError(f"expected 760 frozen samples, found {len(manifest)}")

    keys = [(row["sample"], row["tool"]) for row in signal_rows]
    duplicates = [key for key, count in Counter(keys).items() if count != 1]
    if duplicates:
        raise ValueError(f"duplicate signal tool/sample rows: {duplicates[:5]}")
    expected = {(sample, tool) for sample in manifest for tool in TOOLS.values()}
    observed = set(keys)
    if observed != expected:
        raise ValueError(
            f"signal matrix mismatch: missing={len(expected - observed)} "
            f"extra={len(observed - expected)}"
        )
    for row in signal_rows:
        source = manifest[row["sample"]]
        for field in ("cohort", "condition", "role", "input_reads_used"):
            if row[field] != source[field]:
                raise ValueError(f"signal/manifest {field} mismatch for {row['sample']}/{row['tool']}")
    return manifest


def _summary(values: list[float]) -> tuple[float, float]:
    if not values:
        raise ValueError("summary requires at least one value")
    return statistics.fmean(values), statistics.median(values)


def build_tables(
    manifest_rows: list[dict[str, str]],
    signal_rows: list[dict[str, str]],
    thresholds: tuple[int, ...] = DEFAULT_THRESHOLDS,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    manifest = validate_inputs(manifest_rows, signal_rows)
    by_tool = {
        output_tool: [row for row in signal_rows if row["tool"] == source_tool]
        for output_tool, source_tool in TOOLS.items()
    }
    overall_rows: list[dict[str, object]] = []
    floor_rows: list[dict[str, object]] = []
    for tool in TOOLS:
        rows = by_tool[tool]
        units = {row["signal_unit"] for row in rows}
        if len(units) != 1:
            raise ValueError(f"{tool} contains mixed signal units: {sorted(units)}")
        unit = next(iter(units))
        values = [float(row["fna_c2_signal"]) for row in rows]
        paper_values = [float(manifest[row["sample"]]["paper_fna_c2_percent"]) for row in rows]
        positive = [
            float(row["fna_c2_signal"])
            for row in rows
            if row["role"] == "paper_c2_positive"
        ]
        negative = [
            float(row["fna_c2_signal"])
            for row in rows
            if row["role"] == "paper_zero_fna"
        ]
        positive_mean, positive_median = _summary(positive)
        negative_mean, negative_median = _summary(negative)
        overall_rows.append(
            {
                "tool": tool,
                "unit": unit,
                "n": len(rows),
                "n_c2_positive": len(positive),
                "n_zero_fna": len(negative),
                "c2_positive_mean": positive_mean,
                "zero_fna_mean": negative_mean,
                "c2_positive_median": positive_median,
                "zero_fna_median": negative_median,
                "positive_minus_zero_mean": positive_mean - negative_mean,
                "auc_c2_positive_vs_zero": auc_with_ties(positive, negative),
                "pearson_paper_c2": pearson(values, paper_values),
                "spearman_paper_c2": spearman(values, paper_values),
            }
        )
        for threshold in thresholds:
            selected_positive = [
                float(row["fna_c2_signal"])
                for row in rows
                if row["role"] == "paper_c2_positive"
                and float(row["expected_c2_reads_at_depth"]) >= threshold
            ]
            selected_mean, selected_median = _summary(selected_positive)
            floor_rows.append(
                {
                    "tool": tool,
                    "unit": unit,
                    "positive_rule": "all_pos" if threshold == 0 else f">={threshold}_expected_reads",
                    "expected_c2_reads_min": threshold,
                    "n_pos": len(selected_positive),
                    "n_zero": len(negative),
                    "pos_mean": selected_mean,
                    "zero_mean": negative_mean,
                    "pos_median": selected_median,
                    "zero_median": negative_median,
                    "delta_mean": selected_mean - negative_mean,
                    "auc_c2_positive_vs_zero": auc_with_ties(selected_positive, negative),
                }
            )
    return overall_rows, floor_rows


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    field: f"{row[field]:.6f}" if isinstance(row.get(field), float) else row[field]
                    for field in fields
                }
            )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Recompute frozen Fna C2 sample-level summaries and detection-floor AUCs."
    )
    parser.add_argument("--sample-manifest", required=True, type=Path)
    parser.add_argument("--signal-table", required=True, type=Path)
    parser.add_argument("--overall-out", required=True, type=Path)
    parser.add_argument("--floor-out", required=True, type=Path)
    args = parser.parse_args(argv)

    overall_rows, floor_rows = build_tables(
        read_tsv(args.sample_manifest),
        read_tsv(args.signal_table),
    )
    write_tsv(args.overall_out, overall_rows, OVERALL_FIELDS)
    write_tsv(args.floor_out, floor_rows, FLOOR_FIELDS)
    print(f"wrote {len(overall_rows)} overall rows and {len(floor_rows)} floor rows")


if __name__ == "__main__":
    main()
