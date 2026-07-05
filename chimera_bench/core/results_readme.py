from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, List

from .metrics import METRIC_VERSION


TOOL_DISPLAY_NAMES = {
    # configs/experiments still use "ganon", but the actual software is ganon2.
    "ganon": "ganon2",
}
TOOL_DIR_NAMES = {
    # reverse mapping for looking up run dirs on disk
    "ganon2": "ganon",
}
SUPPLEMENTAL_MAIN_RESULT_RUNS = [
    Path("results/reruns/prjna_single_read_alltools_20260513_163846"),
    Path("results/reruns/chimera_1_7_20260701_003612"),
]
BUILD_README_EXCLUDED_TOOLS = {"bracken"}

BUILD_COLUMNS = [
    "Tool",
    "DB Name",
    "Elapsed Seconds",
    "Max RSS (KB)",
    "DB Size",
    "Started At",
    "Finished At",
]

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
    ("Completeness (species)", "completeness_species"),
    ("Purity (species)", "purity_species"),
    ("L1 Norm (species)", "l1_norm_species"),
    ("Completeness (genus)", "completeness_genus"),
    ("Purity (genus)", "purity_genus"),
    ("L1 Norm (genus)", "l1_norm_genus"),
    ("Weighted UniFrac", "weighted_unifrac"),
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

ABUNDANCE_MAIN_COLUMNS = [
    ("Elapsed (s)", "run_elapsed_seconds"),
    ("Max RSS (GB)", "resource_max_rss_gb"),
    ("Completeness (species)", "completeness_species"),
    ("Purity (species)", "purity_species"),
    ("L1 Norm (species)", "l1_norm_species"),
    ("Completeness (genus)", "completeness_genus"),
    ("Purity (genus)", "purity_genus"),
    ("L1 Norm (genus)", "l1_norm_genus"),
    ("Weighted UniFrac", "weighted_unifrac"),
]

PUBLIC_METRIC_ALIASES = {
    "per_read_classified_rate": "exact_per_read_classified_rate",
    "per_read_unclassified_rate": "exact_per_read_unclassified_rate",
    "per_read_truth_mapped_rate_species": "exact_per_read_truth_mapped_rate_species",
    "per_read_truth_mapped_rate_genus": "exact_per_read_truth_mapped_rate_genus",
    "per_read_pred_mapped_rate_species": "exact_per_read_pred_mapped_rate_species",
    "per_read_pred_mapped_rate_genus": "exact_per_read_pred_mapped_rate_genus",
    "per_read_precision_species": "exact_per_read_precision_species",
    "per_read_recall_species": "exact_per_read_recall_species",
    "per_read_f1_species": "exact_per_read_f1_species",
    "per_read_precision_genus": "exact_per_read_precision_genus",
    "per_read_recall_genus": "exact_per_read_recall_genus",
    "per_read_f1_genus": "exact_per_read_f1_genus",
}

CLASSIFY_PUBLIC_NOTE = [
    "本表按 `species/genus` 层级判断分类是否正确：先将真值和预测提升到对应层级，再比较是否一致。",
    "例如：真值为某个 species，预测到该 species 下的 strain 时，`species` 列记为正确。",
    "原始真值若为 `strain/subspecies/isolate`，只要能提升到对应 rank，就进入该 rank 分母；原始真值若只有 `family/order/class`，则不进入更细层级分母。",
    "Truth Mapped Rate / Pred Mapped Rate 会同时展示；F1 必须结合映射率一起解读，不能脱离分母单独引用。",
    "只有满足默认参数设置的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。",
]

PROFILE_PUBLIC_NOTE = [
    "本表按 OPAL core 计算 profile 结果：`species/genus` 两层分别报告 Completeness、Purity、L1 Norm。",
    "只统计工具原生输出的 profile 文件；没有原生 profile 的工具不会出现在本表里。",
    "只包含当前 profile 评估版本的结果；旧版结果不进入公开表。",
    "只有满足默认参数设置的 run 才能作为论文主表结果引用；部分工具的历史结果仍需刷新，请以 `results/README.md` 的状态节为准。",
    "Weighted UniFrac 是基于整棵 taxonomy tree 的全局距离；Completeness、Purity 越大越好，L1 Norm、Weighted UniFrac 越小越好。",
]

COLLECTION_CLASSIFY_NOTE = (
    "集合聚合（collection aggregate）行的 `Samples` 为已完成 sample 数；"
    "`Elapsed` 为 sample 运行时间总和，`Max RSS` 为最大值，read 数为总和；"
    "率值和准确率为 sample 算术平均。"
)

COLLECTION_PROFILE_NOTE = (
    "集合聚合（collection aggregate）行的 `Samples` 为已完成 sample 数；"
    "`Elapsed` 为 sample 运行时间总和，`Max RSS` 为最大值；"
    "profile 指标为 sample 算术平均。"
)

SUM_METRIC_KEYS = {
    "run_elapsed_seconds",
    "total_reads",
    "classified_reads",
    "unclassified_reads",
}

MAX_METRIC_KEYS = {
    "resource_max_rss_gb",
}


def _format_value(value):
    if value is None:
        return ""
    if isinstance(value, float):
        text = f"{value:.6f}"
        return text.rstrip("0").rstrip(".")
    return str(value)


def _apply_public_metric_aliases(metrics: dict) -> dict:
    aliased = dict(metrics)
    for public_key, preferred_key in PUBLIC_METRIC_ALIASES.items():
        preferred_value = metrics.get(preferred_key)
        if preferred_value is not None:
            aliased[public_key] = preferred_value
    return aliased


def _collect_runs(root: Path) -> List[Dict]:
    records = []
    if not root.exists():
        return records
    for meta_path in root.rglob("meta.json"):
        try:
            meta = json.loads(meta_path.read_text())
        except json.JSONDecodeError:
            continue
        if meta.get("return_code") not in {None, 0}:
            continue
        metrics_path = meta_path.parent / "metrics.json"
        metrics = {}
        if metrics_path.exists():
            try:
                metrics = json.loads(metrics_path.read_text())
            except json.JSONDecodeError:
                metrics = {}
        metrics = _apply_public_metric_aliases(metrics)
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
                "dataset_collection": meta.get("dataset_collection"),
                "display_dataset": meta.get("display_dataset"),
                "sample_id": meta.get("sample_id"),
                "db_name": db_name,
                "metrics": metrics,
            }
        )
    return records


def _supplemental_classify_roots(root: Path) -> list[Path]:
    try:
        is_main_classify = root.resolve() == Path("results/classify").resolve()
    except OSError:
        is_main_classify = False
    if not is_main_classify:
        return []
    return [run_root / "classify" for run_root in SUPPLEMENTAL_MAIN_RESULT_RUNS if (run_root / "classify").exists()]


def _supplemental_build_roots(root: Path) -> list[Path]:
    try:
        is_main_builds = root.resolve() == Path("results/builds").resolve()
    except OSError:
        is_main_builds = False
    if not is_main_builds:
        return []
    return [run_root / "builds" for run_root in SUPPLEMENTAL_MAIN_RESULT_RUNS if (run_root / "builds").exists()]


def _record_result_key(rec: Dict) -> tuple[str, str, str, str] | None:
    tool = rec.get("tool")
    dataset = rec.get("dataset_collection") or rec.get("display_dataset") or rec.get("dataset")
    if not tool or not dataset:
        return None
    sample_id = ""
    if rec.get("dataset_collection"):
        sample_id = str(rec.get("sample_id") or rec.get("dataset") or "")
    return (
        str(dataset),
        sample_id,
        str(tool),
        str(rec.get("db_name") or ""),
    )


def _prefer_later_records(records: list[Dict]) -> list[Dict]:
    ordered_keys: list[tuple[str, str, str, str]] = []
    rows_by_key: dict[tuple[str, str, str, str], Dict] = {}
    passthrough: list[Dict] = []
    for rec in records:
        key = _record_result_key(rec)
        if key is None:
            passthrough.append(rec)
            continue
        if key not in rows_by_key:
            ordered_keys.append(key)
        rows_by_key[key] = rec
    return passthrough + [rows_by_key[key] for key in ordered_keys]


def _collect_runs_with_supplements(root: Path) -> list[Dict]:
    records = _collect_runs(root)
    for supplemental_root in _supplemental_classify_roots(root):
        records.extend(_collect_runs(supplemental_root))
    return _prefer_later_records(records)


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
        header_map = {header[i].strip(): row[i].strip() for i in range(min(len(header), len(row)))}
        tool = header_map.get("Tool", row[0].strip() if len(row) > 0 else "")
        db = header_map.get("DB Name", header_map.get("DB", row[1].strip() if len(row) > 1 else ""))
        if not tool:
            continue
        tool = TOOL_DISPLAY_NAMES.get(tool, tool)
        if tool in BUILD_README_EXCLUDED_TOOLS:
            continue
        normalized = [
            tool,
            db,
            header_map.get("Elapsed Seconds", ""),
            header_map.get("Max RSS (KB)", ""),
            header_map.get("DB Size", header_map.get("DB Size (GiB)", "")),
            header_map.get("Started At", ""),
            header_map.get("Finished At", ""),
        ]
        rows[(tool, db)] = normalized
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


def _dir_size_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    if path.is_file():
        try:
            return path.stat().st_size
        except OSError:
            return 0
    total = 0
    for root, _dirs, files in os.walk(path, followlinks=False):
        root_path = Path(root)
        for name in files:
            file_path = root_path / name
            try:
                total += file_path.stat().st_size
            except OSError:
                continue
    return total


def _format_bytes(size_bytes: int | None) -> str:
    if size_bytes is None or size_bytes < 0:
        return ""
    units = ["B", "KiB", "MiB", "GiB", "TiB", "PiB"]
    value = float(size_bytes)
    unit = units[0]
    for candidate in units[1:]:
        if value < 1024.0:
            break
        value /= 1024.0
        unit = candidate
    if unit == "B":
        return f"{int(value)} {unit}"
    if value >= 100:
        return f"{value:.0f} {unit}"
    if value >= 10:
        return f"{value:.1f} {unit}"
    return f"{value:.2f} {unit}"


def _resolve_build_db_size(meta_path: Path, meta: dict) -> str:
    run_dir = meta_path.parent
    size_candidates: list[int] = []
    db_dir = run_dir / "DB"
    if db_dir.exists():
        size_candidates.append(_dir_size_bytes(db_dir))

    outputs = meta.get("outputs", {}) or {}
    candidates: list[Path] = []
    for key in ("db_file", "db_prefix"):
        value = outputs.get(key)
        if isinstance(value, str) and value.strip():
            candidates.append(Path(value))
    db_prefix = meta.get("db_prefix")
    if isinstance(db_prefix, str) and db_prefix.strip():
        candidates.append(Path(db_prefix))

    for candidate in candidates:
        path = candidate if candidate.is_absolute() else (run_dir / candidate)
        if path.exists():
            size_candidates.append(_dir_size_bytes(path))
            continue
        if path.suffix != ".imcf":
            imcf_path = path.with_suffix(".imcf")
            if imcf_path.exists():
                size_candidates.append(_dir_size_bytes(imcf_path))

    if size_candidates:
        return _format_bytes(max(size_candidates))
    return ""


def _resolve_build_db_size_from_run_dir(root: Path, tool: str, db_name: str) -> str:
    tool_dir = TOOL_DIR_NAMES.get(tool, tool)
    db_dir = root / tool_dir / db_name / "DB"
    if db_dir.exists():
        return _format_bytes(_dir_size_bytes(db_dir))
    return ""


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
        if key.startswith("completeness_") or key.startswith("purity_") or key.startswith("l1_norm_"):
            return True
        if key == "weighted_unifrac":
            return True
    return False


def _has_current_profile_metrics(metrics: dict) -> bool:
    if not _has_abundance_metrics(metrics):
        return False
    return metrics.get("profile_metric_version") == METRIC_VERSION


def _numeric_value(value) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _aggregate_metric(records: list[Dict], metric_key: str) -> float | int | None:
    values = []
    for rec in records:
        value = _numeric_value((rec.get("metrics") or {}).get(metric_key))
        if value is not None:
            values.append(value)
    if not values:
        return None
    if metric_key in SUM_METRIC_KEYS:
        total = sum(values)
        return int(total) if total.is_integer() else total
    if metric_key in MAX_METRIC_KEYS:
        return max(values)
    return sum(values) / len(values)


def _aggregate_collection_records(records: list[Dict], columns: list[tuple[str, str]]) -> list[Dict]:
    """Collapse sample-level runs into one display row per collection/tool/DB."""
    output: list[Dict] = []
    groups: dict[tuple[str, str, str], list[Dict]] = {}
    metric_keys = [key for _, key in columns]

    for rec in records:
        collection = rec.get("dataset_collection")
        if not collection:
            item = dict(rec)
            item["display_dataset"] = rec.get("dataset")
            output.append(item)
            continue
        tool = rec.get("tool")
        if not tool:
            continue
        key = (str(collection), str(tool), str(rec.get("db_name") or ""))
        groups.setdefault(key, []).append(rec)

    for (collection, tool, db_name), grouped in groups.items():
        display_dataset = next(
            (str(rec["display_dataset"]) for rec in grouped if rec.get("display_dataset")),
            collection,
        )
        sample_ids = sorted(
            {
                str(rec.get("sample_id") or rec.get("dataset"))
                for rec in grouped
                if rec.get("sample_id") or rec.get("dataset")
            }
        )
        metrics = {metric_key: _aggregate_metric(grouped, metric_key) for metric_key in metric_keys}
        output.append(
            {
                "exp": grouped[0].get("exp"),
                "tool": tool,
                "dataset": collection,
                "display_dataset": display_dataset,
                "dataset_collection": collection,
                "sample_count": len(sample_ids) or len(grouped),
                "sample_ids": sample_ids,
                "db_name": db_name,
                "metrics": metrics,
            }
        )

    return output


def _display_dataset(rec: Dict) -> str | None:
    dataset = rec.get("display_dataset") or rec.get("dataset_collection") or rec.get("dataset")
    return str(dataset) if dataset else None


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
    records = [r for r in _collect_runs_with_supplements(root) if _has_per_read_metrics(r.get("metrics", {}))]
    records = _aggregate_collection_records(records, PER_READ_MAIN_COLUMNS)
    readme_path = root / "README.md"
    lines = ["# Classify Results", "", "Auto-generated. Do not edit.", ""]
    if not records:
        readme_path.write_text("\n".join(lines) + "\n")
        return

    datasets = sorted({dataset for r in records if (dataset := _display_dataset(r))})
    for dataset in datasets:
        rows = [r for r in records if _display_dataset(r) == dataset]
        has_collection_rows = any(r.get("dataset_collection") for r in rows)
        header = ["Tool", "DB"]
        if has_collection_rows:
            header.append("Samples")
        header.extend([label for label, _ in PER_READ_MAIN_COLUMNS])
        lines.append(f"## Dataset: {dataset}")
        lines.append("")
        lines.append("### Per-read Metrics")
        lines.append("")
        lines.extend(CLASSIFY_PUBLIC_NOTE)
        if has_collection_rows:
            lines.append(COLLECTION_CLASSIFY_NOTE)
        lines.append("")
        lines.append("| " + " | ".join(header) + " |")
        lines.append("| " + " | ".join(["---"] * len(header)) + " |")
        for rec in sorted(rows, key=lambda r: (r.get("tool") or "", r.get("db_name") or "")):
            metrics = rec.get("metrics", {})
            row = [rec.get("tool") or "", rec.get("db_name") or ""]
            if has_collection_rows:
                row.append(_format_value(rec.get("sample_count")))
            row += [_format_value(metrics.get(key)) for _, key in PER_READ_MAIN_COLUMNS]
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")
    readme_path.write_text("\n".join(lines) + "\n")


def write_profile_readme(profile_root: Path, runs_root: Path) -> None:
    records = [
        r
        for r in _collect_runs_with_supplements(runs_root)
        if _has_current_profile_metrics(r.get("metrics", {}))
    ]
    records = _aggregate_collection_records(records, ABUNDANCE_MAIN_COLUMNS)
    readme_path = profile_root / "README.md"

    lines = ["# Profile Results", "", "Auto-generated. Do not edit.", ""]
    if not records:
        readme_path.write_text("\n".join(lines) + "\n")
        return

    datasets = sorted({dataset for r in records if (dataset := _display_dataset(r))})
    for dataset in datasets:
        rows = [r for r in records if _display_dataset(r) == dataset]
        has_collection_rows = any(r.get("dataset_collection") for r in rows)
        header = ["Tool", "DB"]
        if has_collection_rows:
            header.append("Samples")
        header.extend([label for label, _ in ABUNDANCE_MAIN_COLUMNS])
        lines.append(f"## Dataset: {dataset}")
        lines.append("")
        lines.append("### Abundance Metrics")
        lines.append("")
        lines.extend(PROFILE_PUBLIC_NOTE)
        if has_collection_rows:
            lines.append(COLLECTION_PROFILE_NOTE)
        lines.append("")
        lines.append("| " + " | ".join(header) + " |")
        lines.append("| " + " | ".join(["---"] * len(header)) + " |")
        for rec in sorted(rows, key=lambda r: (r.get("tool") or "", r.get("db_name") or "")):
            metrics = rec.get("metrics", {})
            row = [rec.get("tool") or "", rec.get("db_name") or ""]
            if has_collection_rows:
                row.append(_format_value(rec.get("sample_count")))
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

    header = BUILD_COLUMNS
    rows_by_key = _parse_build_readme_rows(readme_path)
    # Build outputs may contain very large DB directories; avoid a full recursive walk.
    for build_root in [root, *_supplemental_build_roots(root)]:
        for meta_path in build_root.glob("*/*/meta.json"):
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
            if tool in BUILD_README_EXCLUDED_TOOLS:
                continue
            db_name = meta.get("db_name") or meta_path.parent.name
            db_size = _resolve_build_db_size(meta_path, meta)
            rows_by_key[(str(tool), str(db_name))] = [
                str(tool),
                str(db_name),
                _format_value(meta.get("elapsed_seconds")),
                _format_value(resource.get("max_rss_kb")),
                db_size,
                _format_value(meta.get("started_at")),
                _format_value(meta.get("finished_at")),
            ]

    # Backfill DB size for preserved rows that don't have current meta.json.
    for (tool, db_name), row in rows_by_key.items():
        if len(row) < len(header):
            row.extend([""] * (len(header) - len(row)))
        if row[4]:
            continue
        row[4] = _resolve_build_db_size_from_run_dir(root, tool, db_name)

    lines.append("| " + " | ".join(header) + " |")
    lines.append("| " + " | ".join(["---"] * len(header)) + " |")
    for _, row in sorted(rows_by_key.items(), key=lambda r: (r[0][0], r[0][1])):
        if len(row) < len(header):
            row = row + [""] * (len(header) - len(row))
        elif len(row) > len(header):
            row = row[: len(header)]
        lines.append("| " + " | ".join(row) + " |")

    readme_path.write_text("\n".join(lines) + "\n")
