from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List


PER_READ_COLUMNS = [
    ("Elapsed (s)", "run_elapsed_seconds"),
    ("Max RSS (GB)", "resource_max_rss_gb"),
    ("Total Reads", "total_reads"),
    ("Classified Reads", "classified_reads"),
    ("Unclassified Reads", "unclassified_reads"),
    ("Classified Rate", "per_read_classified_rate"),
    ("Unclassified Rate", "per_read_unclassified_rate"),
    ("Precision (species)", "per_read_precision_species"),
    ("Recall (species)", "per_read_recall_species"),
    ("F1 (species)", "per_read_f1_species"),
    ("Precision (genus)", "per_read_precision_genus"),
    ("Recall (genus)", "per_read_recall_genus"),
    ("F1 (genus)", "per_read_f1_genus"),
]

ABUNDANCE_COLUMNS = [
    ("Elapsed (s)", "run_elapsed_seconds"),
    ("Max RSS (GB)", "resource_max_rss_gb"),
    ("Presence Precision (species)", "presence_precision_species"),
    ("Presence Recall (species)", "presence_recall_species"),
    ("Presence F1 (species)", "presence_f1_species"),
    ("L1 (species)", "abundance_l1_species"),
    ("TV (species)", "abundance_tv_species"),
    ("Bray-Curtis (species)", "abundance_bc_species"),
    ("Presence Precision (genus)", "presence_precision_genus"),
    ("Presence Recall (genus)", "presence_recall_genus"),
    ("Presence F1 (genus)", "presence_f1_genus"),
    ("L1 (genus)", "abundance_l1_genus"),
    ("TV (genus)", "abundance_tv_genus"),
    ("Bray-Curtis (genus)", "abundance_bc_genus"),
]


def _format_value(value):
    if value is None:
        return ""
    if isinstance(value, float):
        text = f"{value:.6f}"
        return text.rstrip("0").rstrip(".")
    return str(value)


def _collect_runs(root: Path) -> List[Dict]:
    records = []
    if not root.exists():
        return records
    for meta_path in root.rglob("meta.json"):
        try:
            meta = json.loads(meta_path.read_text())
        except json.JSONDecodeError:
            continue
        metrics_path = meta_path.parent / "metrics.json"
        metrics = {}
        if metrics_path.exists():
            try:
                metrics = json.loads(metrics_path.read_text())
            except json.JSONDecodeError:
                metrics = {}
        if meta.get("elapsed_seconds") is not None:
            metrics.setdefault("run_elapsed_seconds", meta.get("elapsed_seconds"))
        resource = meta.get("resource", {})
        max_rss_kb = resource.get("max_rss_kb")
        if max_rss_kb:
            metrics.setdefault("resource_max_rss_gb", max_rss_kb / (1024 * 1024))
        records.append(
            {
                "exp": meta.get("exp"),
                "tool": meta.get("tool"),
                "dataset": meta.get("dataset"),
                "db_name": meta.get("db_name"),
                "metrics": metrics,
            }
        )
    return records


def _append_table(lines: list[str], title: str, records: list[Dict], columns: list[tuple[str, str]]):
    lines.append(f"### {title}")
    lines.append("")
    header = ["Tool", "DB"] + [label for label, _ in columns]
    lines.append("| " + " | ".join(header) + " |")
    lines.append("| " + " | ".join(["---"] * len(header)) + " |")
    for rec in sorted(records, key=lambda r: r.get("tool") or ""):
        metrics = rec.get("metrics", {})
        row = [rec.get("tool") or "", rec.get("db_name") or ""]
        row += [_format_value(metrics.get(key)) for _, key in columns]
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")


def write_classify_readme(root: Path) -> None:
    records = _collect_runs(root)
    lines = ["# Classify Results", "", "Auto-generated. Do not edit.", ""]
    if not records:
        (root / "README.md").write_text("\n".join(lines) + "\n")
        return

    grouped: Dict[str, List[Dict]] = {}
    for rec in records:
        grouped.setdefault(rec.get("dataset") or "unknown", []).append(rec)

    for dataset in sorted(grouped.keys()):
        lines.append(f"## Dataset: {dataset}")
        lines.append("")
        _append_table(lines, "Per-read Metrics", grouped[dataset], PER_READ_COLUMNS)
    (root / "README.md").write_text("\n".join(lines) + "\n")


def write_profile_readme(profile_root: Path, runs_root: Path) -> None:
    records = _collect_runs(runs_root)
    lines = ["# Profile Results", "", "Auto-generated. Do not edit.", ""]
    if not records:
        (profile_root / "README.md").write_text("\n".join(lines) + "\n")
        return

    grouped: Dict[str, List[Dict]] = {}
    for rec in records:
        grouped.setdefault(rec.get("dataset") or "unknown", []).append(rec)

    for dataset in sorted(grouped.keys()):
        lines.append(f"## Dataset: {dataset}")
        lines.append("")
        _append_table(lines, "Abundance Metrics", grouped[dataset], ABUNDANCE_COLUMNS)
    (profile_root / "README.md").write_text("\n".join(lines) + "\n")


def write_builds_readme(root: Path) -> None:
    lines = ["# Build Results", "", "Auto-generated. Do not edit.", ""]
    if not root.exists():
        (root / "README.md").write_text("\n".join(lines) + "\n")
        return

    rows = []
    for meta_path in root.rglob("meta.json"):
        try:
            meta = json.loads(meta_path.read_text())
        except json.JSONDecodeError:
            continue
        resource = meta.get("resource", {})
        rows.append(
            {
                "tool": meta.get("tool"),
                "db_name": meta.get("db_name") or meta_path.parent.name,
                "elapsed_seconds": meta.get("elapsed_seconds"),
                "max_rss_kb": resource.get("max_rss_kb"),
                "started_at": meta.get("started_at"),
                "finished_at": meta.get("finished_at"),
            }
        )

    build_columns = [
        ("Tool", "tool"),
        ("DB Name", "db_name"),
        ("Elapsed Seconds", "elapsed_seconds"),
        ("Max RSS (KB)", "max_rss_kb"),
        ("Started At", "started_at"),
        ("Finished At", "finished_at"),
    ]
    header = ["Tool"] + [label for label, _ in build_columns[1:]]
    lines.append("| " + " | ".join(header) + " |")
    lines.append("| " + " | ".join(["---"] * len(header)) + " |")
    for row in sorted(rows, key=lambda r: (r.get("tool") or "", r.get("db_name") or "")):
        values = [_format_value(row.get(key)) for _, key in build_columns]
        lines.append("| " + " | ".join(values) + " |")

    (root / "README.md").write_text("\n".join(lines) + "\n")
