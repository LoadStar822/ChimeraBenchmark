#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))

from analyze_crc_association import (  # noqa: E402
    COHORT_FIELDS,
    META_FIELDS,
    _analyze_prepared_metrics,
    analyze_published_fna_c2,
    load_frozen_metadata,
    load_paper_metadata,
    match_metadata_rows,
    write_tsv,
)


RAW_TOOL_LABELS = {
    "chimera_default": "Chimera",
    "centrifuger_default": "Centrifuger",
    "kraken2_safe_lf01": "Kraken2_LF01",
    "sylph_default": "sylph",
}
TOOL_LABELS = tuple(RAW_TOOL_LABELS.values())
UNIT_SPECS = {
    "reads_per_million": ("raw_c2_reads_per_million", 1_000_000.0),
    "sequence_abundance_pct": ("raw_c2_sequence_abundance_percent", 100.0),
}


def infer_cohort(sample: str) -> str:
    if sample.startswith("2015_Yu_"):
        return "2015 Yu et al"
    if sample.startswith("2018_Wirbel_"):
        return "2018 Wirbel et al"
    if sample.startswith("DRR"):
        return "YachidaS_2019"
    raise ValueError(f"cannot infer cohort from sample: {sample}")


def read_signal_tables(paths: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for path in paths:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = set(reader.fieldnames or [])
            required = {"tool", "sample", "condition"}
            missing = required - fields
            if missing:
                raise ValueError(f"{path} is missing columns: {sorted(missing)}")
            normalized = {"signal_unit", "fna_c2_signal"} <= fields
            legacy = {"unit", "Fna_C2_signal"} <= fields
            if not normalized and not legacy:
                raise ValueError(
                    f"{path} requires signal_unit/fna_c2_signal or unit/Fna_C2_signal"
                )
            for row in reader:
                raw_tool = row["tool"]
                tool = RAW_TOOL_LABELS.get(raw_tool, raw_tool)
                if tool not in TOOL_LABELS:
                    raise ValueError(f"unsupported tool in signal table: {raw_tool}")
                canonical = {
                    **row,
                    "tool": tool,
                    "unit": row.get("unit") or row["signal_unit"],
                    "Fna_C2_signal": row.get("Fna_C2_signal") or row["fna_c2_signal"],
                    "cohort": row.get("cohort") or infer_cohort(row["sample"]),
                }
                key = tool, row["sample"]
                if key in seen:
                    raise ValueError(f"duplicate tool/sample signal row: {key}")
                seen.add(key)
                rows.append(canonical)
    return rows


def prepare_tool_signal(
    tool: str,
    metadata_rows: list[dict[str, object]],
    sample_rows: list[dict[str, str]],
) -> tuple[list[dict[str, object]], str]:
    units = {row["unit"] for row in sample_rows}
    if len(units) != 1:
        raise ValueError(f"{tool}: mixed signal units: {sorted(units)}")
    unit = next(iter(units))
    if unit not in UNIT_SPECS:
        raise ValueError(f"{tool}: unsupported signal unit: {unit}")
    metric, denominator = UNIT_SPECS[unit]

    matched = sorted(
        match_metadata_rows(metadata_rows, sample_rows),
        key=lambda row: (str(row["cohort"]), str(row["paper_sample"])),
    )
    if len(matched) != len(sample_rows):
        raise ValueError(
            f"{tool}: matched {len(matched)} metadata rows for {len(sample_rows)} signals"
        )
    matched_samples = [str(row["sample"]) for row in matched]
    if len(set(matched_samples)) != len(matched_samples):
        raise ValueError(f"{tool}: signal samples are not unique after metadata matching")

    prepared: list[dict[str, object]] = []
    for row in matched:
        condition = str(row["condition"]).lower()
        audit_condition = str(row["audit_condition"]).lower()
        sex = str(row["sex"]).lower()
        if condition not in {"disease", "control"} or condition != audit_condition:
            raise ValueError(f"{tool}: condition mismatch for {row['sample']}")
        if sex not in {"female", "male"}:
            raise ValueError(f"{tool}: unsupported sex {sex!r} for {row['sample']}")
        signal = float(str(row["Fna_C2_signal"]))
        fraction = signal / denominator
        if not 0.0 <= fraction <= 1.0:
            raise ValueError(f"{tool}: signal outside fraction range for {row['sample']}")
        prepared.append(
            {
                **row,
                "case": int(condition == "disease"),
                "male": int(sex == "male"),
                metric: signal,
                "sample_level_c2_signal_fraction": fraction,
            }
        )
    return prepared, metric


QC_FIELDS = [
    "tool",
    "signal_unit",
    "paper_metadata_samples",
    "signal_samples",
    "matched_samples",
    "unique_matched_samples",
    "status",
]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Compare sample-level Fna C2 signals against CRC phenotype across three cohorts."
    )
    metadata_group = parser.add_mutually_exclusive_group(required=True)
    metadata_group.add_argument("--supplement-xlsx", type=Path)
    metadata_group.add_argument("--metadata-table", type=Path)
    parser.add_argument("--signal-table", required=True, action="append", type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20_260_710)
    args = parser.parse_args(argv)

    metadata_rows = (
        load_paper_metadata(args.supplement_xlsx)
        if args.supplement_xlsx
        else load_frozen_metadata(args.metadata_table)
    )
    signal_rows = read_signal_tables(args.signal_table)
    by_tool: dict[str, list[dict[str, str]]] = {}
    for row in signal_rows:
        by_tool.setdefault(row["tool"], []).append(row)
    unknown = set(by_tool) - set(TOOL_LABELS)
    if unknown:
        raise ValueError(f"unsupported tools in signal tables: {sorted(unknown)}")
    missing_tools = set(TOOL_LABELS) - set(by_tool)
    if missing_tools:
        raise ValueError(f"missing tools in signal tables: {sorted(missing_tools)}")

    cohort_rows: list[dict[str, object]] = []
    meta_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    for label in TOOL_LABELS:
        rows = by_tool[label]
        prepared, metric = prepare_tool_signal(label, metadata_rows, rows)
        tool_cohort, tool_meta = _analyze_prepared_metrics(
            label,
            prepared,
            [(metric, "sample_level_c2_signal_fraction")],
            bootstrap_iterations=args.bootstrap,
            seed=args.seed,
        )
        cohort_rows.extend(tool_cohort)
        meta_rows.extend(tool_meta)
        qc_rows.append(
            {
                "tool": label,
                "signal_unit": rows[0]["unit"],
                "paper_metadata_samples": len(metadata_rows),
                "signal_samples": len(rows),
                "matched_samples": len(prepared),
                "unique_matched_samples": len({str(row["sample"]) for row in prepared}),
                "status": "pass",
            }
        )

    paper_cohort, paper_meta = analyze_published_fna_c2(
        metadata_rows,
        bootstrap_iterations=args.bootstrap,
        seed=args.seed,
    )
    cohort_rows.extend(paper_cohort)
    meta_rows.extend(paper_meta)

    write_tsv(args.out_dir / "crc_sample_signal_association_cohort.tsv", cohort_rows, COHORT_FIELDS)
    write_tsv(args.out_dir / "crc_sample_signal_association_meta.tsv", meta_rows, META_FIELDS)
    write_tsv(args.out_dir / "crc_sample_signal_association_qc.tsv", qc_rows, QC_FIELDS)
    print(
        f"wrote sample-level CRC association for {len(TOOL_LABELS)} tools, "
        f"the published reference, and {len(metadata_rows)} samples"
    )


if __name__ == "__main__":
    main()
