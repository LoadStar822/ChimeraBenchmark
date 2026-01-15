from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List


TOOL_DISPLAY_NAMES = {
    # configs/experiments still use "ganon", but the actual software is ganon2.
    "ganon": "ganon2",
}

PER_READ_COLUMNS = [
    ("Elapsed (s)", "run_elapsed_seconds"),
    ("Max RSS (GB)", "resource_max_rss_gb"),
    ("Total Reads", "total_reads"),
    ("Classified Reads", "classified_reads"),
    ("Unclassified Reads", "unclassified_reads"),
    ("Classified Rate", "per_read_classified_rate"),
    ("Unclassified Rate", "per_read_unclassified_rate"),
    ("Truth Mapped Rate (species)", "per_read_truth_mapped_rate_species"),
    ("Pred Mapped Rate (species)", "per_read_pred_mapped_rate_species"),
    ("Precision (species)", "per_read_precision_species"),
    ("Recall (species)", "per_read_recall_species"),
    ("F1 (species)", "per_read_f1_species"),
    ("Truth Mapped Rate (genus)", "per_read_truth_mapped_rate_genus"),
    ("Pred Mapped Rate (genus)", "per_read_pred_mapped_rate_genus"),
    ("Precision (genus)", "per_read_precision_genus"),
    ("Recall (genus)", "per_read_recall_genus"),
    ("F1 (genus)", "per_read_f1_genus"),
]

PER_READ_UNK_COLUMNS = [
    ("Elapsed (s)", "run_elapsed_seconds"),
    ("Max RSS (GB)", "resource_max_rss_gb"),
    ("Precision (species, UNK)", "per_read_precision_species_unk"),
    ("Recall (species, UNK)", "per_read_recall_species_unk"),
    ("F1 (species, UNK)", "per_read_f1_species_unk"),
    ("Precision (genus, UNK)", "per_read_precision_genus_unk"),
    ("Recall (genus, UNK)", "per_read_recall_genus_unk"),
    ("F1 (genus, UNK)", "per_read_f1_genus_unk"),
]

ABUNDANCE_COLUMNS = [
    ("Elapsed (s)", "run_elapsed_seconds"),
    ("Max RSS (GB)", "resource_max_rss_gb"),
    ("Truth Mapped Rate (species)", "truth_mapped_mass_rate_species"),
    ("Pred Mapped Rate (species)", "pred_mapped_mass_rate_species"),
    ("Presence Precision (species)", "presence_precision_species"),
    ("Presence Recall (species)", "presence_recall_species"),
    ("Presence F1 (species)", "presence_f1_species"),
    ("L1 (species)", "abundance_l1_species"),
    ("TV (species)", "abundance_tv_species"),
    ("Bray-Curtis (species)", "abundance_bc_species"),
    ("Truth Mapped Rate (genus)", "truth_mapped_mass_rate_genus"),
    ("Pred Mapped Rate (genus)", "pred_mapped_mass_rate_genus"),
    ("Presence Precision (genus)", "presence_precision_genus"),
    ("Presence Recall (genus)", "presence_recall_genus"),
    ("Presence F1 (genus)", "presence_f1_genus"),
    ("L1 (genus)", "abundance_l1_genus"),
    ("TV (genus)", "abundance_tv_genus"),
    ("Bray-Curtis (genus)", "abundance_bc_genus"),
]

ABUNDANCE_UNK_COLUMNS = [
    ("Elapsed (s)", "run_elapsed_seconds"),
    ("Max RSS (GB)", "resource_max_rss_gb"),
    ("Presence Precision (species, UNK)", "presence_precision_species_unk"),
    ("Presence Recall (species, UNK)", "presence_recall_species_unk"),
    ("Presence F1 (species, UNK)", "presence_f1_species_unk"),
    ("L1 (species, UNK)", "abundance_l1_species_unk"),
    ("TV (species, UNK)", "abundance_tv_species_unk"),
    ("Bray-Curtis (species, UNK)", "abundance_bc_species_unk"),
    ("Presence Precision (genus, UNK)", "presence_precision_genus_unk"),
    ("Presence Recall (genus, UNK)", "presence_recall_genus_unk"),
    ("Presence F1 (genus, UNK)", "presence_f1_genus_unk"),
    ("L1 (genus, UNK)", "abundance_l1_genus_unk"),
    ("TV (genus, UNK)", "abundance_tv_genus_unk"),
    ("Bray-Curtis (genus, UNK)", "abundance_bc_genus_unk"),
]

PER_READ_MAIN_COLUMNS = [
    ("Elapsed (s)", "run_elapsed_seconds"),
    ("Max RSS (GB)", "resource_max_rss_gb"),
    ("Precision (species)", "per_read_precision_species"),
    ("Recall (species)", "per_read_recall_species"),
    ("F1 (species)", "per_read_f1_species"),
    ("Precision (genus)", "per_read_precision_genus"),
    ("Recall (genus)", "per_read_recall_genus"),
    ("F1 (genus)", "per_read_f1_genus"),
]

ABUNDANCE_MAIN_COLUMNS = [
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
        db_name = meta.get("db_name")
        if not db_name:
            db_path = meta.get("db")
            if db_path:
                db_name = Path(db_path).name
        tool = meta.get("tool")
        if isinstance(tool, str):
            tool = TOOL_DISPLAY_NAMES.get(tool, tool)
        records.append(
            {
                "exp": meta.get("exp"),
                "tool": tool,
                "dataset": meta.get("dataset"),
                "db_name": db_name,
                "metrics": metrics,
            }
        )
    return records


def _split_md_row(line: str) -> list[str]:
    return [c.strip() for c in line.strip().strip("|").split("|")]


def _parse_readme_rows(readme_path: Path) -> tuple[list[str], dict[str, dict[tuple[str, str], list[str]]]]:
    """Parse an existing results README into per-dataset rows.

    Returns:
      - dataset_order: datasets in appearance order
      - rows_by_dataset: dataset -> {(tool, db): cells}
    """
    if not readme_path.exists():
        return [], {}
    text = readme_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    dataset_order: list[str] = []
    rows_by_dataset: dict[str, dict[tuple[str, str], list[str]]] = {}
    current_dataset: str | None = None
    header: list[str] | None = None
    for line in text:
        if line.startswith("## Dataset:"):
            current_dataset = line.split(":", 1)[1].strip()
            if current_dataset and current_dataset not in dataset_order:
                dataset_order.append(current_dataset)
            header = None
            continue
        if current_dataset is None:
            continue
        if line.startswith("| Tool | DB |"):
            header = _split_md_row(line)
            continue
        if header is None:
            continue
        if not line.startswith("|"):
            continue
        row = _split_md_row(line)
        if not row or set(row) == {"---"}:
            continue
        # tolerate mismatch by padding/truncating to header length
        if len(row) < len(header):
            row = row + [""] * (len(header) - len(row))
        elif len(row) > len(header):
            row = row[: len(header)]
        tool = row[0].strip() if len(row) > 0 else ""
        db = row[1].strip() if len(row) > 1 else ""
        if not tool:
            continue
        rows_by_dataset.setdefault(current_dataset, {})[(tool, db)] = row
    return dataset_order, rows_by_dataset


def _parse_build_readme_rows(readme_path: Path) -> dict[tuple[str, str], list[str]]:
    if not readme_path.exists():
        return {}
    text = readme_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    header: list[str] | None = None
    rows: dict[tuple[str, str], list[str]] = {}
    for line in text:
        if line.startswith("| Tool |"):
            header = _split_md_row(line)
            continue
        if header is None:
            continue
        if not line.startswith("|"):
            continue
        row = _split_md_row(line)
        if not row or set(row) == {"---"}:
            continue
        if len(row) < len(header):
            row = row + [""] * (len(header) - len(row))
        elif len(row) > len(header):
            row = row[: len(header)]
        tool = row[0].strip() if len(row) > 0 else ""
        db = row[1].strip() if len(row) > 1 else ""
        if not tool:
            continue
        tool = TOOL_DISPLAY_NAMES.get(tool, tool)
        row[0] = tool
        rows[(tool, db)] = row
    return rows


def _format_float(value: float) -> str:
    return _format_value(float(value))


def _try_float(text: str) -> float | None:
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _normalize_abundance_l1_from_tv(
    header: list[str], rows: dict[tuple[str, str], list[str]]
) -> None:
    """Ensure L1 columns are percent-points (0..200), matching TV/BC columns."""
    try:
        l1_species = header.index("L1 (species, UNK)")
        tv_species = header.index("TV (species, UNK)")
    except ValueError:
        l1_species = tv_species = None
    try:
        l1_genus = header.index("L1 (genus, UNK)")
        tv_genus = header.index("TV (genus, UNK)")
    except ValueError:
        l1_genus = tv_genus = None

    for row in rows.values():
        if l1_species is not None and tv_species is not None:
            tv = _try_float(row[tv_species])
            if tv is not None:
                row[l1_species] = _format_float(200.0 * tv)
        if l1_genus is not None and tv_genus is not None:
            tv = _try_float(row[tv_genus])
            if tv is not None:
                row[l1_genus] = _format_float(200.0 * tv)


def _merge_rows(
    *,
    existing_order: list[str],
    existing_rows: dict[str, dict[tuple[str, str], list[str]]],
    records: list[Dict],
    columns: list[tuple[str, str]],
) -> tuple[list[str], dict[str, dict[tuple[str, str], list[str]]], list[str]]:
    header = ["Tool", "DB"] + [label for label, _ in columns]
    dataset_order = list(existing_order)
    rows_by_dataset = {k: dict(v) for k, v in existing_rows.items()}

    for rec in records:
        dataset = rec.get("dataset")
        tool = rec.get("tool")
        if not dataset or not tool:
            continue
        dataset = str(dataset)
        tool = str(tool)
        if dataset not in dataset_order:
            dataset_order.append(dataset)
        rows = rows_by_dataset.setdefault(dataset, {})

        db = rec.get("db_name") or ""
        if db:
            row_key = (tool, str(db))
        else:
            candidates = [k for k in rows.keys() if k[0] == tool]
            if len(candidates) == 1:
                row_key = candidates[0]
                db = row_key[1]
            else:
                row_key = (tool, "")

        metrics = rec.get("metrics", {})
        row = [tool, str(db)] + [_format_value(metrics.get(metric_key)) for _, metric_key in columns]
        if len(row) < len(header):
            row = row + [""] * (len(header) - len(row))
        rows[row_key] = row

    return dataset_order, rows_by_dataset, header


def _has_per_read_metrics(metrics: dict) -> bool:
    for key in metrics.keys():
        if key.startswith("per_read_"):
            return True
    return False


def _has_abundance_metrics(metrics: dict) -> bool:
    for key in metrics.keys():
        if key.startswith("presence_") or key.startswith("abundance_"):
            return True
    return False


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
    records = [r for r in _collect_runs(root) if _has_per_read_metrics(r.get("metrics", {}))]
    readme_path = root / "README.md"
    lines = ["# Classify Results", "", "Auto-generated. Do not edit.", ""]
    if not records:
        readme_path.write_text("\n".join(lines) + "\n")
        return

    header = ["Tool", "DB"] + [label for label, _ in PER_READ_MAIN_COLUMNS]
    datasets = sorted({r.get("dataset") for r in records if r.get("dataset")})
    for dataset in datasets:
        rows = [r for r in records if r.get("dataset") == dataset]
        lines.append(f"## Dataset: {dataset}")
        lines.append("")
        lines.append("### Per-read Metrics")
        lines.append("")
        lines.append("| " + " | ".join(header) + " |")
        lines.append("| " + " | ".join(["---"] * len(header)) + " |")
        for rec in sorted(rows, key=lambda r: (r.get("tool") or "", r.get("db_name") or "")):
            metrics = rec.get("metrics", {})
            row = [rec.get("tool") or "", rec.get("db_name") or ""]
            row += [_format_value(metrics.get(key)) for _, key in PER_READ_MAIN_COLUMNS]
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")
    readme_path.write_text("\n".join(lines) + "\n")


def write_profile_readme(profile_root: Path, runs_root: Path) -> None:
    records = [r for r in _collect_runs(runs_root) if _has_abundance_metrics(r.get("metrics", {}))]
    readme_path = profile_root / "README.md"

    lines = ["# Profile Results", "", "Auto-generated. Do not edit.", ""]
    if not records:
        readme_path.write_text("\n".join(lines) + "\n")
        return

    header = ["Tool", "DB"] + [label for label, _ in ABUNDANCE_MAIN_COLUMNS]
    datasets = sorted({r.get("dataset") for r in records if r.get("dataset")})
    for dataset in datasets:
        rows = [r for r in records if r.get("dataset") == dataset]
        lines.append(f"## Dataset: {dataset}")
        lines.append("")
        lines.append("### Abundance Metrics")
        lines.append("")
        lines.append("| " + " | ".join(header) + " |")
        lines.append("| " + " | ".join(["---"] * len(header)) + " |")
        for rec in sorted(rows, key=lambda r: (r.get("tool") or "", r.get("db_name") or "")):
            metrics = rec.get("metrics", {})
            row = [rec.get("tool") or "", rec.get("db_name") or ""]
            row += [_format_value(metrics.get(key)) for _, key in ABUNDANCE_MAIN_COLUMNS]
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")
    readme_path.write_text("\n".join(lines) + "\n")


def write_builds_readme(root: Path) -> None:
    readme_path = root / "README.md"
    lines = ["# Build Results", "", "Auto-generated. Do not edit.", ""]
    if not root.exists():
        readme_path.write_text("\n".join(lines) + "\n")
        return

    header = ["Tool", "DB Name", "Elapsed Seconds", "Max RSS (KB)", "Started At", "Finished At"]
    rows_by_key = _parse_build_readme_rows(readme_path)
    for meta_path in root.rglob("meta.json"):
        try:
            meta = json.loads(meta_path.read_text())
        except json.JSONDecodeError:
            continue
        if meta.get("return_code") != 0:
            continue
        resource = meta.get("resource", {})
        tool = meta.get("tool")
        if isinstance(tool, str):
            tool = TOOL_DISPLAY_NAMES.get(tool, tool)
        if not tool:
            continue
        db_name = meta.get("db_name") or meta_path.parent.name
        rows_by_key[(str(tool), str(db_name))] = [
            str(tool),
            str(db_name),
            _format_value(meta.get("elapsed_seconds")),
            _format_value(resource.get("max_rss_kb")),
            _format_value(meta.get("started_at")),
            _format_value(meta.get("finished_at")),
        ]

    lines.append("| " + " | ".join(header) + " |")
    lines.append("| " + " | ".join(["---"] * len(header)) + " |")
    for _, row in sorted(rows_by_key.items(), key=lambda r: (r[0][0], r[0][1])):
        if len(row) < len(header):
            row = row + [""] * (len(header) - len(row))
        elif len(row) > len(header):
            row = row[: len(header)]
        lines.append("| " + " | ".join(row) + " |")

    readme_path.write_text("\n".join(lines) + "\n")
