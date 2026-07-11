from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Iterable

from .config import load_yaml_dir
from .core.results_readme import (
    ABUNDANCE_MAIN_COLUMNS,
    PER_READ_MAIN_COLUMNS,
    _aggregate_collection_records,
    _collect_runs_with_supplements,
    _display_dataset,
    _has_current_profile_metrics,
    _has_per_read_metrics,
    _parse_build_readme_rows,
)


FORMAL_CLASSIFY_DATASETS = (
    "cami-strain-madness-long",
    "cami-strain-madness-short",
    "cami2-marine-long",
    "cami2-marine-short",
    "prjna637878-supported19-single-read",
)
FORMAL_CLASSIFY_TOOLS = ("centrifuger", "chimera", "kraken2")
FORMAL_PROFILE_DATASETS = (
    "atcc-hifi",
    "atcc-illumina",
    "cami-strain-madness-long",
    "cami-strain-madness-short",
    "cami2-marine-long",
    "cami2-marine-short",
    "prjna637878-supported19-single-read",
    "zymo-gridion-even",
    "zymo-gridion-log",
    "zymo-promethion-even",
    "zymo-promethion-log",
)
SUPPLEMENTARY_PROFILE_TOOLS = ("bracken", "centrifuger", "chimera", "sylph")
FORMAL_BUILD_TOOLS = ("centrifuger", "chimera", "kraken2", "sylph")
TOOL_VERSIONS = {
    "bracken": "2.9",
    "centrifuger": "1.1.0-r291",
    "chimera": "1.7",
    "kraken2": "2.1.3",
    "sylph": "0.8.1",
}
FNA_TOOL_VERSIONS = {
    "Chimera": "1.6.3",
    "Centrifuger": "v1.0.12-10-gdf5bd32",
    "Kraken2_LF01": "2.1.3",
    "sylph": "0.8.1",
}
FNA_SIGNAL_TOOLS = ("Centrifuger", "Chimera", "Kraken2_LF01", "sylph")
FNA_AUDIT_TOOLS = ("Centrifuger", "Chimera", "Kraken2_LF01")


def _validate_complete_matrix(
    rows: list[dict],
    *,
    datasets: Iterable[str],
    tools: Iterable[str],
    label: str,
) -> None:
    expected = {(dataset, tool) for dataset in datasets for tool in tools}
    keys = [(str(row["dataset"]), str(row["tool"])) for row in rows]
    actual = set(keys)
    duplicates = sorted(key for key, count in Counter(keys).items() if count != 1)
    missing = sorted(expected - actual)
    extra = sorted(actual - expected)
    if duplicates or missing or extra:
        raise ValueError(
            f"{label} matrix mismatch: missing={missing} extra={extra} duplicates={duplicates}"
        )


def _validate_tool_sample_coverage(
    rows: list[dict[str, str]],
    *,
    expected_samples: set[str],
    expected_tools: Iterable[str],
    label: str,
) -> None:
    expected_tool_set = set(expected_tools)
    actual_tools = {row["tool"] for row in rows}
    keys = [(row["sample"], row["tool"]) for row in rows]
    duplicate_keys = sorted(key for key, count in Counter(keys).items() if count != 1)
    coverage = {
        tool: {row["sample"] for row in rows if row["tool"] == tool}
        for tool in actual_tools
    }
    incomplete = {
        tool: {
            "missing": sorted(expected_samples - coverage.get(tool, set())),
            "extra": sorted(coverage.get(tool, set()) - expected_samples),
        }
        for tool in expected_tool_set | actual_tools
        if coverage.get(tool, set()) != expected_samples
    }
    if actual_tools != expected_tool_set or duplicate_keys or incomplete:
        raise ValueError(
            f"{label} tool/sample coverage mismatch: "
            f"tools={sorted(actual_tools)} expected_tools={sorted(expected_tool_set)} "
            f"duplicates={duplicate_keys} incomplete={incomplete}"
        )


def _metric_row(record: dict, columns: list[tuple[str, str]]) -> dict[str, object]:
    metrics = record.get("metrics") or {}
    row: dict[str, object] = {
        "dataset": _display_dataset(record) or "",
        "tool": record.get("tool") or "",
        "db": record.get("db_name") or "",
    }
    for _, key in columns:
        row[key] = metrics.get(key)
    return row


def build_classify_rows(records: list[dict], dataset_samples: dict[str, int]) -> list[dict]:
    aggregated = _aggregate_collection_records(records, PER_READ_MAIN_COLUMNS)
    rows = []
    for record in aggregated:
        dataset = _display_dataset(record)
        tool = record.get("tool")
        if dataset not in FORMAL_CLASSIFY_DATASETS or tool not in FORMAL_CLASSIFY_TOOLS:
            continue
        if record.get("dataset_collection"):
            observed = int(record.get("sample_count") or 0)
            expected = dataset_samples[dataset]
            if observed != expected:
                raise ValueError(
                    f"classify {dataset}/{tool}: expected {expected} completed samples, found {observed}"
                )
        row = _metric_row(record, PER_READ_MAIN_COLUMNS)
        row.update(
            {
                "samples": dataset_samples[dataset],
                "status": "complete",
                "paper_scope": "main",
                "tool_version": TOOL_VERSIONS[str(tool)],
            }
        )
        rows.append(row)
    return sorted(
        rows,
        key=lambda row: (
            FORMAL_CLASSIFY_DATASETS.index(str(row["dataset"])),
            str(row["tool"]),
        ),
    )


def build_profile_rows(records: list[dict], dataset_samples: dict[str, int]) -> list[dict]:
    aggregated = _aggregate_collection_records(records, ABUNDANCE_MAIN_COLUMNS)
    rows = []
    for record in aggregated:
        dataset = _display_dataset(record)
        tool = record.get("tool")
        if (
            not dataset
            or dataset not in FORMAL_PROFILE_DATASETS
            or tool not in SUPPLEMENTARY_PROFILE_TOOLS
        ):
            continue
        if record.get("dataset_collection"):
            observed = int(record.get("sample_count") or 0)
            expected = dataset_samples[dataset]
            if observed != expected:
                raise ValueError(
                    f"profile {dataset}/{tool}: expected {expected} completed samples, found {observed}"
                )
        row = _metric_row(record, ABUNDANCE_MAIN_COLUMNS)
        row.update(
            {
                "samples": dataset_samples[dataset],
                "status": "complete",
                "paper_scope": "supplementary",
                "tool_version": TOOL_VERSIONS[str(tool)],
            }
        )
        rows.append(row)
    return sorted(rows, key=lambda row: (str(row["dataset"]), str(row["tool"])))


def build_classify_sample_rows(records: list[dict]) -> list[dict]:
    rows = []
    for record in records:
        dataset = _display_dataset(record)
        tool = record.get("tool")
        if dataset != "prjna637878-supported19-single-read" or tool not in FORMAL_CLASSIFY_TOOLS:
            continue
        row = _metric_row(record, PER_READ_MAIN_COLUMNS)
        row.update(
            {
                "sample_id": record.get("sample_id") or record.get("dataset") or "",
                "status": "complete",
                "paper_scope": "main",
                "tool_version": TOOL_VERSIONS[str(tool)],
            }
        )
        rows.append(row)
    return sorted(rows, key=lambda row: (str(row["sample_id"]), str(row["tool"])))


def _write_tsv(path: Path, rows: Iterable[dict], fields: list[str]) -> None:
    materialized = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(materialized)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _real_manifest_rows(results_root: Path) -> list[dict]:
    real_root = results_root / "real/fna_c2_crc3_head3m"
    signal_path = real_root / "sample_level_signals.tsv"
    audit_path = real_root / "read_audit_sample_metrics.tsv"
    if not signal_path.exists() or not audit_path.exists():
        return []

    signal_rows = _read_tsv(signal_path)
    audit_rows = _read_tsv(audit_path)
    signal_samples = {row["sample"] for row in signal_rows}
    audit_samples = {row["sample"] for row in audit_rows}
    if signal_samples != audit_samples:
        raise ValueError("Fna sample-signal and read-audit sample sets differ")
    _validate_tool_sample_coverage(
        signal_rows,
        expected_samples=signal_samples,
        expected_tools=FNA_SIGNAL_TOOLS,
        label="Fna sample-level signal",
    )
    _validate_tool_sample_coverage(
        audit_rows,
        expected_samples=signal_samples,
        expected_tools=FNA_AUDIT_TOOLS,
        label="Fna read audit",
    )
    sample_count = len(signal_samples)

    rows = []
    for task, source_rows, scope, result_file in (
        (
            "real_sample_signal",
            signal_rows,
            "supplementary",
            "results/real/fna_c2_crc3_head3m/sample_level_signals.tsv",
        ),
        (
            "real_read_audit",
            audit_rows,
            "main",
            "results/real/fna_c2_crc3_head3m/read_audit_sample_metrics.tsv",
        ),
    ):
        for tool in sorted({row["tool"] for row in source_rows}):
            policy = "default_except_threads_and_output_paths"
            if tool == "Kraken2_LF01":
                policy += "; db_load_factor_0.1"
            rows.append(
                {
                    "task": task,
                    "dataset": "fna-c2-crc3-head3m",
                    "tool": tool,
                    "db": "frozen_fna_near_neighbor",
                    "samples": sample_count,
                    "status": "complete",
                    "paper_scope": scope,
                    "tool_version": FNA_TOOL_VERSIONS[tool],
                    "parameter_policy": policy,
                    "result_file": result_file,
                }
            )
    return rows


def _dataset_sample_counts(config_root: Path) -> dict[str, int]:
    datasets = load_yaml_dir(config_root / "datasets")
    counts = {}
    for name, config in datasets.items():
        catalog_samples = (config.get("catalog") or {}).get("samples")
        if catalog_samples is not None:
            counts[name] = int(catalog_samples)
        elif isinstance(config.get("samples"), list):
            counts[name] = len(config["samples"])
        elif isinstance(config.get("sample_ids"), list):
            counts[name] = len(config["sample_ids"])
        elif config.get("group") == "cami2-marine" and isinstance(config.get("reads"), list):
            counts[name] = len(config["reads"])
        else:
            counts[name] = 1
    return counts


def _build_rows(results_root: Path) -> list[dict]:
    parsed = _parse_build_readme_rows(results_root / "builds/README.md")
    rows = []
    for (tool, db), values in parsed.items():
        if tool not in FORMAL_BUILD_TOOLS or db != "cami_refseq":
            continue
        max_rss_kb = float(values[3]) if values[3] else None
        rows.append(
            {
                "tool": tool,
                "db": db,
                "elapsed_seconds": float(values[2]) if values[2] else None,
                "max_rss_gb": max_rss_kb / (1024 * 1024) if max_rss_kb is not None else None,
                "db_size": values[4],
                "started_at": values[5],
                "finished_at": values[6],
                "status": "complete",
                "paper_scope": "supplementary",
                "tool_version": TOOL_VERSIONS[tool],
            }
        )
    return sorted(rows, key=lambda row: str(row["tool"]))


def write_paper_tables(*, config_root: Path, results_root: Path) -> dict[str, int]:
    sample_counts = _dataset_sample_counts(config_root)
    raw_records = _collect_runs_with_supplements(results_root / "classify")
    classify_records = [record for record in raw_records if _has_per_read_metrics(record.get("metrics") or {})]
    profile_records = [
        record for record in raw_records if _has_current_profile_metrics(record.get("metrics") or {})
    ]

    classify_rows = build_classify_rows(classify_records, sample_counts)
    classify_sample_rows = build_classify_sample_rows(classify_records)
    profile_rows = build_profile_rows(profile_records, sample_counts)
    build_rows = _build_rows(results_root)

    _validate_complete_matrix(
        classify_rows,
        datasets=FORMAL_CLASSIFY_DATASETS,
        tools=FORMAL_CLASSIFY_TOOLS,
        label="classify",
    )
    _validate_complete_matrix(
        profile_rows,
        datasets=FORMAL_PROFILE_DATASETS,
        tools=SUPPLEMENTARY_PROFILE_TOOLS,
        label="profile",
    )
    _validate_complete_matrix(
        [{"dataset": "cami_refseq", **row} for row in build_rows],
        datasets=("cami_refseq",),
        tools=FORMAL_BUILD_TOOLS,
        label="build",
    )
    prjna_samples = {
        str(row["sample_id"])
        for row in classify_sample_rows
        if row["dataset"] == "prjna637878-supported19-single-read"
    }
    if len(prjna_samples) != sample_counts["prjna637878-supported19-single-read"]:
        raise ValueError(
            "PRJNA single-read sample coverage mismatch: "
            f"expected={sample_counts['prjna637878-supported19-single-read']} "
            f"found={len(prjna_samples)}"
        )
    _validate_tool_sample_coverage(
        [
            {"sample": str(row["sample_id"]), "tool": str(row["tool"])}
            for row in classify_sample_rows
        ],
        expected_samples=prjna_samples,
        expected_tools=FORMAL_CLASSIFY_TOOLS,
        label="PRJNA single-read classify",
    )

    common_fields = ["dataset", "tool", "db", "samples"]
    classify_fields = common_fields + [key for _, key in PER_READ_MAIN_COLUMNS] + [
        "status",
        "paper_scope",
        "tool_version",
    ]
    profile_fields = common_fields + [key for _, key in ABUNDANCE_MAIN_COLUMNS] + [
        "status",
        "paper_scope",
        "tool_version",
    ]
    sample_fields = ["dataset", "sample_id", "tool", "db"] + [
        key for _, key in PER_READ_MAIN_COLUMNS
    ] + ["status", "paper_scope", "tool_version"]
    build_fields = [
        "tool",
        "db",
        "elapsed_seconds",
        "max_rss_gb",
        "db_size",
        "started_at",
        "finished_at",
        "status",
        "paper_scope",
        "tool_version",
    ]

    _write_tsv(results_root / "classify/summary.tsv", classify_rows, classify_fields)
    _write_tsv(results_root / "classify/sample_metrics.tsv", classify_sample_rows, sample_fields)
    _write_tsv(results_root / "profile/summary.tsv", profile_rows, profile_fields)
    _write_tsv(results_root / "builds/summary.tsv", build_rows, build_fields)

    manifest_rows = []
    for task, rows, result_file in (
        ("classify", classify_rows, "results/classify/summary.tsv"),
        ("profile", profile_rows, "results/profile/summary.tsv"),
        ("build", build_rows, "results/builds/summary.tsv"),
    ):
        for row in rows:
            manifest_rows.append(
                {
                    "task": task,
                    "dataset": row.get("dataset", "cami_refseq"),
                    "tool": row["tool"],
                    "db": row["db"],
                    "samples": row.get("samples", ""),
                    "status": row["status"],
                    "paper_scope": row["paper_scope"],
                    "tool_version": row["tool_version"],
                    "parameter_policy": "default_except_threads_and_output_paths",
                    "result_file": result_file,
                }
            )
    manifest_rows.extend(_real_manifest_rows(results_root))
    _write_tsv(
        results_root / "paper_run_manifest.tsv",
        manifest_rows,
        [
            "task",
            "dataset",
            "tool",
            "db",
            "samples",
            "status",
            "paper_scope",
            "tool_version",
            "parameter_policy",
            "result_file",
        ],
    )
    return {
        "build": len(build_rows),
        "classify": len(classify_rows),
        "classify_sample": len(classify_sample_rows),
        "profile": len(profile_rows),
        "manifest": len(manifest_rows),
    }
