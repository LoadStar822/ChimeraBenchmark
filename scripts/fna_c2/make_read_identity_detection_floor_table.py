#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


DEFAULT_THRESHOLDS = (0, 500, 1000, 2000, 5000)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def auc_with_ties(positive: list[float], negative: list[float]) -> float | None:
    if not positive or not negative:
        return None

    score = 0.0
    for pos in positive:
        for neg in negative:
            if pos > neg:
                score += 1.0
            elif pos == neg:
                score += 0.5
    return score / (len(positive) * len(negative))


def fmt(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6f}"


def mean(values: list[float]) -> float | None:
    return statistics.fmean(values) if values else None


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def weighted_rate(rows: list[dict[str, str]], numerator_field: str) -> float | None:
    denom = sum(int(row["selected_reads"]) for row in rows)
    if not denom:
        return None
    return sum(int(row[numerator_field]) for row in rows) / denom


def build_rows(sample_rows: list[dict[str, str]], thresholds: list[int]) -> list[dict[str, str | int]]:
    cohorts = ["all", *sorted({row["cohort"] for row in sample_rows})]
    out: list[dict[str, str | int]] = []

    for cohort in cohorts:
        cohort_rows = sample_rows if cohort == "all" else [row for row in sample_rows if row["cohort"] == cohort]
        zero_rows = [row for row in cohort_rows if row["role"] == "paper_zero_fna"]
        zero_supported = [float(row["c2_supported_rate"]) for row in zero_rows]
        zero_strict = [float(row["strict_c2_best_rate"]) for row in zero_rows]

        for threshold in thresholds:
            positive_rows = [
                row
                for row in cohort_rows
                if row["role"] == "paper_c2_positive"
                and float(row["expected_c2_reads_at_depth"]) >= threshold
            ]
            pos_supported = [float(row["c2_supported_rate"]) for row in positive_rows]
            pos_strict = [float(row["strict_c2_best_rate"]) for row in positive_rows]

            out.append(
                {
                    "cohort": cohort,
                    "threshold_expected_c2_reads": threshold,
                    "n_positive": len(positive_rows),
                    "n_zero": len(zero_rows),
                    "pos_mean_c2_supported_rate": fmt(mean(pos_supported)),
                    "pos_median_c2_supported_rate": fmt(median(pos_supported)),
                    "pos_weighted_c2_supported_rate": fmt(weighted_rate(positive_rows, "c2_supported_reads")),
                    "zero_mean_c2_supported_rate": fmt(mean(zero_supported)),
                    "zero_median_c2_supported_rate": fmt(median(zero_supported)),
                    "zero_weighted_c2_supported_rate": fmt(weighted_rate(zero_rows, "c2_supported_reads")),
                    "auc_supported_rate_pos_vs_zero": fmt(auc_with_ties(pos_supported, zero_supported)),
                    "pos_mean_strict_c2_best_rate": fmt(mean(pos_strict)),
                    "zero_mean_strict_c2_best_rate": fmt(mean(zero_strict)),
                    "pos_weighted_unmapped_rate": fmt(weighted_rate(positive_rows, "unmapped_reads")),
                    "zero_weighted_unmapped_rate": fmt(weighted_rate(zero_rows, "unmapped_reads")),
                }
            )
    return out


def write_tsv(path: Path, rows: list[dict[str, str | int]]) -> None:
    if not rows:
        raise ValueError(f"no rows to write: {path}")

    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build Fna C2 read-identity detection-floor tables from Chimera C2 read audit summaries."
    )
    parser.add_argument("--sample-summary", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--threshold",
        action="append",
        type=int,
        dest="thresholds",
        default=[],
        help="Minimum paper-expected C2 reads for positives. Can be repeated.",
    )
    args = parser.parse_args()

    thresholds = args.thresholds or list(DEFAULT_THRESHOLDS)
    write_tsv(args.out, build_rows(read_tsv(args.sample_summary), thresholds))


if __name__ == "__main__":
    main()
