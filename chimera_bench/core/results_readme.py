from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List


CLASSIFY_COLUMNS = [
    "total_reads",
    "classified_reads",
    "unclassified_reads",
    "per_read_classified_rate",
    "per_read_unclassified_rate",
    "per_read_precision_species",
    "per_read_recall_species",
    "per_read_f1_species",
    "per_read_precision_genus",
    "per_read_recall_genus",
    "per_read_f1_genus",
]

PROFILE_COLUMNS = [
    "presence_precision_species",
    "presence_recall_species",
    "presence_f1_species",
    "abundance_l1_species",
    "abundance_tv_species",
    "abundance_bc_species",
    "presence_precision_genus",
    "presence_recall_genus",
    "presence_f1_genus",
    "abundance_l1_genus",
    "abundance_tv_genus",
    "abundance_bc_genus",
]


def _format_value(value):
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
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
        records.append(
            {
                "exp": meta.get("exp"),
                "tool": meta.get("tool"),
                "dataset": meta.get("dataset"),
                "metrics": metrics,
            }
        )
    return records


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
        header = ["tool"] + CLASSIFY_COLUMNS
        lines.append("| " + " | ".join(header) + " |")
        lines.append("| " + " | ".join(["---"] * len(header)) + " |")
        for rec in sorted(grouped[dataset], key=lambda r: r.get("tool") or ""):
            metrics = rec.get("metrics", {})
            row = [rec.get("tool") or ""]
            row += [_format_value(metrics.get(col)) for col in CLASSIFY_COLUMNS]
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")

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
        header = ["tool"] + PROFILE_COLUMNS
        lines.append("| " + " | ".join(header) + " |")
        lines.append("| " + " | ".join(["---"] * len(header)) + " |")
        for rec in sorted(grouped[dataset], key=lambda r: r.get("tool") or ""):
            metrics = rec.get("metrics", {})
            row = [rec.get("tool") or ""]
            row += [_format_value(metrics.get(col)) for col in PROFILE_COLUMNS]
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")

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

    header = ["tool", "db_name", "elapsed_seconds", "max_rss_kb", "started_at", "finished_at"]
    lines.append("| " + " | ".join(header) + " |")
    lines.append("| " + " | ".join(["---"] * len(header)) + " |")
    for row in sorted(rows, key=lambda r: (r.get("tool") or "", r.get("db_name") or "")):
        lines.append(
            "| "
            + " | ".join(_format_value(row.get(col)) for col in header)
            + " |"
        )

    (root / "README.md").write_text("\n".join(lines) + "\n")
