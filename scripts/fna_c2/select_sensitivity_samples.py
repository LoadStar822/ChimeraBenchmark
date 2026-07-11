#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable


def evenly_spaced(
    rows: list[dict[str, str]],
    count: int,
    *,
    key: Callable[[dict[str, str]], object],
) -> list[dict[str, str]]:
    if count < 0:
        raise ValueError("selection count must be non-negative")
    if count > len(rows):
        raise ValueError(f"cannot select {count} rows from {len(rows)} candidates")
    if count == 0:
        return []
    ordered = sorted(rows, key=lambda row: (key(row), row["sample"]))
    if count == 1:
        return [ordered[len(ordered) // 2]]
    indices = [round(i * (len(ordered) - 1) / (count - 1)) for i in range(count)]
    if len(set(indices)) != len(indices):
        raise AssertionError("evenly spaced selection produced duplicate indices")
    return [ordered[index] for index in indices]


def select_cohort(rows: list[dict[str, str]], *, n_positive: int) -> list[dict[str, str]]:
    positives = [row for row in rows if row["role"] == "paper_c2_positive"]
    zero_fna = [row for row in rows if row["role"] == "paper_zero_fna"]
    selected_positive = evenly_spaced(
        positives,
        n_positive,
        key=lambda row: float(row["paper_fna_c2_pct"]),
    )

    condition_counts = Counter(row["condition"] for row in selected_positive)
    selected_zero: list[dict[str, str]] = []
    for condition, count in sorted(condition_counts.items()):
        candidates = [row for row in zero_fna if row["condition"] == condition]
        selected_zero.extend(evenly_spaced(candidates, count, key=lambda row: row["sample"]))

    selected: list[dict[str, str]] = []
    for order, row in enumerate(selected_positive, start=1):
        selected.append(
            {
                **row,
                "sensitivity_group": "paper_c2_positive_abundance_quantile",
                "selection_order": str(order),
            }
        )
    for order, row in enumerate(sorted(selected_zero, key=lambda item: item["sample"]), start=1):
        selected.append(
            {
                **row,
                "sensitivity_group": "condition_matched_zero_fna",
                "selection_order": str(order),
            }
        )
    return selected


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Select reference-sensitivity samples without reading tool outputs: "
            "published C2 abundance quantiles plus condition-matched zero-Fna controls."
        )
    )
    parser.add_argument("--manifest", action="append", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--positive-per-cohort", type=int, default=10)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    output: list[dict[str, str]] = []
    by_manifest: dict[Path, list[str]] = defaultdict(list)
    for manifest in args.manifest:
        by_cohort: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in read_tsv(manifest):
            by_cohort[row["cohort"]].append(row)
        for cohort, cohort_rows in sorted(by_cohort.items()):
            for row in select_cohort(cohort_rows, n_positive=args.positive_per_cohort):
                output.append({"source_manifest": str(manifest), **row})
                by_manifest[manifest].append(row["sample"])

    fields = [
        "source_manifest",
        "sample",
        "role",
        "cohort",
        "condition",
        "paper_fna_c2_pct",
        "sensitivity_group",
        "selection_order",
    ]
    with (args.out_dir / "sensitivity_sample_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in output)

    for manifest, samples in by_manifest.items():
        path = args.out_dir / f"sample_list.{safe_token(manifest.stem)}.txt"
        path.write_text("".join(f"{sample}\n" for sample in samples))


if __name__ == "__main__":
    main()
