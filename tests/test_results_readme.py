from __future__ import annotations

import json
from pathlib import Path

from chimera_bench.core.results_readme import write_builds_readme, write_classify_readme, write_profile_readme


def test_write_profile_readme_overwrites_and_drops_stale_rows(tmp_path: Path):
    runs_root = tmp_path / "runs"
    profile_root = tmp_path / "profile"
    runs_root.mkdir()
    profile_root.mkdir()

    existing = "\n".join(
        [
            "# Profile Results",
            "",
            "Auto-generated. Do not edit.",
            "",
            "## Dataset: ds1",
            "",
            "### Abundance Metrics (UNK)",
            "",
            "| Tool | DB | Elapsed (s) | Max RSS (GB) | Completeness (species) | Purity (species) | L1 Norm (species) | Completeness (genus) | Purity (genus) | L1 Norm (genus) | Weighted UniFrac |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
            # ganon row should be preserved even if no run exists on disk
            "| ganon | cami_refseq | 1 | 2 | 1 | 1 | 0.5 | 1 | 1 | 0.5 | 0.25 |",
            # sylph row will be updated from on-disk run
            "| sylph | cami_refseq | 9 | 9 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 |",
            "",
        ]
    )
    (profile_root / "README.md").write_text(existing + "\n")

    run_dir = runs_root / "sylph" / "ds1"
    (run_dir / "logs").mkdir(parents=True)
    (run_dir / "outputs").mkdir(parents=True)
    meta = {
        "exp": "sylph",
        "tool": "sylph",
        "dataset": "ds1",
        "db": "/tmp/cami_refseq",
        "db_name": "cami_refseq",
        "elapsed_seconds": 10.0,
        "resource": {"max_rss_kb": 1024 * 1024},
        "outputs": {},
    }
    (run_dir / "meta.json").write_text(json.dumps(meta))
    metrics = {
        "completeness_species": 0.5,
        "purity_species": 0.25,
        "l1_norm_species": 0.75,
        "completeness_genus": 0.6,
        "purity_genus": 0.3,
        "l1_norm_genus": 0.4,
        "weighted_unifrac": 0.123456,
    }
    (run_dir / "metrics.json").write_text(json.dumps(metrics))

    write_profile_readme(profile_root, runs_root)

    text = (profile_root / "README.md").read_text()
    assert "| ganon | cami_refseq |" not in text
    assert "### Abundance Metrics (UNK)" not in text
    assert "本表按 OPAL core 计算 profile 结果" in text
    assert "只统计工具原生输出的 profile 文件" in text
    assert "部分工具的历史结果仍需刷新" in text
    assert "| Tool | DB | Elapsed (s) | Max RSS (GB) | Completeness (species) |" in text
    assert "| sylph | cami_refseq | 10 | 1 | 0.5 | 0.25 | 0.75 | 0.6 | 0.3 | 0.4 | 0.123456 |" in text


def test_write_classify_readme_clears_when_no_per_read_runs(tmp_path: Path):
    runs_root = tmp_path / "runs"
    runs_root.mkdir()

    existing = "\n".join(
        [
            "# Classify Results",
            "",
            "Auto-generated. Do not edit.",
            "",
            "## Dataset: ds1",
            "",
            "### Per-read Metrics (UNK)",
            "",
            "| Tool | DB | Elapsed (s) | Max RSS (GB) | Precision (species, UNK) | Recall (species, UNK) | F1 (species, UNK) | Precision (genus, UNK) | Recall (genus, UNK) | F1 (genus, UNK) |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
            "| ganon | cami_refseq | 1 | 2 | 0.1 | 0.2 | 0.3 | 0.1 | 0.2 | 0.3 |",
            "",
        ]
    )
    (runs_root / "README.md").write_text(existing + "\n")

    # Create a run that has no per-read metrics (e.g. profile-only), so classify README should stay.
    run_dir = runs_root / "sylph" / "ds1"
    run_dir.mkdir(parents=True)
    (run_dir / "meta.json").write_text(json.dumps({"tool": "sylph", "dataset": "ds1"}))
    (run_dir / "metrics.json").write_text(json.dumps({"completeness_species": 1.0}))

    write_classify_readme(runs_root)

    text = (runs_root / "README.md").read_text()
    assert "| ganon | cami_refseq |" not in text
    assert "## Dataset:" not in text


def test_results_readmes_display_ganon_as_ganon2(tmp_path: Path):
    runs_root = tmp_path / "runs"
    profile_root = tmp_path / "profile"
    runs_root.mkdir()
    profile_root.mkdir()

    run_dir = runs_root / "ganon" / "ds1"
    (run_dir / "logs").mkdir(parents=True)
    (run_dir / "outputs").mkdir(parents=True)
    (run_dir / "meta.json").write_text(
        json.dumps(
            {
                "exp": "ganon",
                "tool": "ganon",
                "dataset": "ds1",
                "db": "/tmp/cami_refseq",
                "db_name": "cami_refseq",
                "elapsed_seconds": 10.0,
                "resource": {"max_rss_kb": 1024 * 1024},
                "outputs": {},
            }
        )
    )
    (run_dir / "metrics.json").write_text(
        json.dumps(
            {
                "per_read_precision_species": 0.1,
                "per_read_recall_species": 0.2,
                "per_read_f1_species": 0.3,
                "per_read_precision_genus": 0.4,
                "per_read_recall_genus": 0.5,
                "per_read_f1_genus": 0.6,
                "completeness_species": 0.5,
                "purity_species": 0.25,
                "l1_norm_species": 0.75,
                "completeness_genus": 0.6,
                "purity_genus": 0.3,
                "l1_norm_genus": 0.4,
                "weighted_unifrac": 0.125,
            }
        )
    )

    write_classify_readme(runs_root)
    classify_text = (runs_root / "README.md").read_text()
    assert "本表按 `species/genus` 层级判断分类是否正确" in classify_text
    assert "原始真值若为 `strain/subspecies/isolate`" in classify_text
    assert "Truth Mapped Rate / Pred Mapped Rate 会同时展示" in classify_text
    assert "部分工具的历史结果仍需刷新" in classify_text
    assert "| ganon2 | cami_refseq |" in classify_text
    assert "| ganon | cami_refseq |" not in classify_text

    write_profile_readme(profile_root, runs_root)
    profile_text = (profile_root / "README.md").read_text()
    assert "Weighted UniFrac 是基于整棵 taxonomy tree 的全局距离" in profile_text
    assert "只统计工具原生输出的 profile 文件" in profile_text
    assert "部分工具的历史结果仍需刷新" in profile_text
    assert "| ganon2 | cami_refseq |" in profile_text
    assert "| ganon | cami_refseq |" not in profile_text


def test_results_readmes_prefer_public_rank_aligned_metrics_when_available(tmp_path: Path):
    runs_root = tmp_path / "runs"
    profile_root = tmp_path / "profile"
    runs_root.mkdir()
    profile_root.mkdir()

    run_dir = runs_root / "chimera" / "ds1"
    run_dir.mkdir(parents=True)
    (run_dir / "meta.json").write_text(
        json.dumps(
            {
                "exp": "chimera",
                "tool": "chimera",
                "dataset": "ds1",
                "db": "/tmp/cami_refseq",
                "db_name": "cami_refseq",
                "elapsed_seconds": 10.0,
                "resource": {"max_rss_kb": 1024 * 1024},
                "outputs": {},
                "return_code": 0,
            }
        )
    )
    (run_dir / "metrics.json").write_text(
        json.dumps(
            {
                "per_read_precision_species": 0.11,
                "per_read_recall_species": 0.12,
                "per_read_f1_species": 0.13,
                "per_read_precision_genus": 0.14,
                "per_read_recall_genus": 0.15,
                "per_read_f1_genus": 0.16,
                "exact_per_read_precision_species": 0.21,
                "exact_per_read_recall_species": 0.22,
                "exact_per_read_f1_species": 0.23,
                "exact_per_read_truth_mapped_rate_species": 0.61,
                "exact_per_read_pred_mapped_rate_species": 0.62,
                "exact_per_read_precision_genus": 0.24,
                "exact_per_read_recall_genus": 0.25,
                "exact_per_read_f1_genus": 0.26,
                "exact_per_read_truth_mapped_rate_genus": 0.71,
                "exact_per_read_pred_mapped_rate_genus": 0.72,
                "completeness_species": 0.31,
                "purity_species": 0.32,
                "l1_norm_species": 0.33,
                "completeness_genus": 0.37,
                "purity_genus": 0.38,
                "l1_norm_genus": 0.39,
                "weighted_unifrac": 0.41,
            }
        )
    )

    write_classify_readme(runs_root)
    classify_text = (runs_root / "README.md").read_text()
    assert "| chimera | cami_refseq | 10 | 1 | 0.61 | 0.62 | 0.21 | 0.22 | 0.23 | 0.71 | 0.72 | 0.24 | 0.25 | 0.26 |" in classify_text
    assert "| chimera | cami_refseq | 10 | 1 | 0.11 | 0.12 | 0.13 |" not in classify_text

    write_profile_readme(profile_root, runs_root)
    profile_text = (profile_root / "README.md").read_text()
    assert "| chimera | cami_refseq | 10 | 1 | 0.31 | 0.32 | 0.33 | 0.37 | 0.38 | 0.39 | 0.41 |" in profile_text


def test_write_builds_readme_preserves_existing_rows_and_adds_successful_meta(tmp_path: Path):
    builds_root = tmp_path / "builds"
    builds_root.mkdir()

    existing = "\n".join(
        [
            "# Build Results",
            "",
            "Auto-generated. Do not edit.",
            "",
            "| Tool | DB Name | Elapsed Seconds | Max RSS (KB) | Started At | Finished At |",
            "| --- | --- | --- | --- | --- | --- |",
            "| ganon | cami_refseq | 1 | 2 | a | b |",
            "| sylph | cami_refseq | 3 | 4 | c | d |",
            "",
        ]
    )
    (builds_root / "README.md").write_text(existing + "\n")

    failed_dir = builds_root / "taxor" / "fail_db"
    failed_dir.mkdir(parents=True)
    (failed_dir / "meta.json").write_text(
        json.dumps(
            {
                "tool": "taxor",
                "db_name": "fail_db",
                "return_code": 1,
                "elapsed_seconds": 999.0,
                "resource": {"max_rss_kb": 999},
                "started_at": "x",
                "finished_at": "y",
            }
        )
    )

    success_dir = builds_root / "taxor" / "cami_refseq"
    success_dir.mkdir(parents=True)
    db_dir = success_dir / "DB"
    db_dir.mkdir(parents=True)
    (db_dir / "db.bin").write_bytes(b"x" * 2048)
    (success_dir / "meta.json").write_text(
        json.dumps(
            {
                "tool": "taxor",
                "db_name": "cami_refseq",
                "return_code": 0,
                "elapsed_seconds": 10.0,
                "resource": {"max_rss_kb": 1024},
                "started_at": "2026-01-01T00:00:00+00:00",
                "finished_at": "2026-01-01T00:00:10+00:00",
            }
        )
    )

    write_builds_readme(builds_root)

    text = (builds_root / "README.md").read_text()
    assert "| Tool | DB Name | Elapsed Seconds | Max RSS (KB) | DB Size | Started At | Finished At |" in text
    assert "| ganon2 | cami_refseq | 1 | 2 |  | a | b |" in text
    assert "| ganon | cami_refseq |" not in text
    assert "| sylph | cami_refseq | 3 | 4 |  | c | d |" in text
    assert "| taxor | cami_refseq | 10 | 1024 | 2.00 KiB | 2026-01-01T00:00:00+00:00 | 2026-01-01T00:00:10+00:00 |" in text
    assert "fail_db" not in text


def test_write_builds_readme_drops_bracken_rows(tmp_path: Path):
    builds_root = tmp_path / "builds"
    builds_root.mkdir()
    (builds_root / "README.md").write_text(
        "\n".join(
            [
                "# Build Results",
                "",
                "Auto-generated. Do not edit.",
                "",
                "| Tool | DB Name | Elapsed Seconds | Max RSS (KB) | DB Size | Started At | Finished At |",
                "| --- | --- | --- | --- | --- | --- | --- |",
                "| bracken | cami_refseq | 1 | 2 | 3 GiB | a | b |",
                "| sylph | cami_refseq | 3 | 4 | 5 GiB | c | d |",
                "",
            ]
        )
        + "\n"
    )

    write_builds_readme(builds_root)

    text = (builds_root / "README.md").read_text()
    assert "| bracken | cami_refseq |" not in text
    assert "| sylph | cami_refseq | 3 | 4 | 5 GiB | c | d |" in text
