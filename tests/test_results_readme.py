from __future__ import annotations

import json
from pathlib import Path

from chimera_bench.core.metrics import METRIC_VERSION
from chimera_bench.core import results_readme as results_readme_module
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
        "profile_metric_version": METRIC_VERSION,
    }
    (run_dir / "metrics.json").write_text(json.dumps(metrics))

    old_run_dir = runs_root / "centrifuger" / "legacy"
    old_run_dir.mkdir(parents=True)
    (old_run_dir / "meta.json").write_text(
        json.dumps(
            {
                "tool": "centrifuger",
                "dataset": "legacy",
                "db_name": "cami_refseq",
                "elapsed_seconds": 1.0,
                "resource": {"max_rss_kb": 1024},
                "return_code": 0,
            }
        )
    )
    (old_run_dir / "metrics.json").write_text(
        json.dumps(
            {
                "completeness_species": 1.0,
                "purity_species": 1.0,
                "l1_norm_species": 0.0,
                "weighted_unifrac": 0.0,
                "profile_metric_version": "old-profile-version",
            }
        )
    )

    write_profile_readme(profile_root, runs_root)

    text = (profile_root / "README.md").read_text()
    assert "| ganon | cami_refseq |" not in text
    assert "### Abundance Metrics (UNK)" not in text
    assert "本表按 OPAL core 计算 profile 结果" in text
    assert "只统计工具原生输出的 profile 文件" in text
    assert "部分工具的历史结果仍需刷新" in text
    assert "只包含当前 profile 评估版本的结果" in text
    assert "legacy" not in text
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
                "profile_metric_version": METRIC_VERSION,
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
                "profile_metric_version": METRIC_VERSION,
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


def test_results_readmes_group_dataset_collection_samples(tmp_path: Path):
    runs_root = tmp_path / "runs"
    profile_root = tmp_path / "profile"
    runs_root.mkdir()
    profile_root.mkdir()

    samples = [
        (
            "S1",
            10.0,
            1024 * 1024,
            {
                "exact_per_read_truth_mapped_rate_species": 0.2,
                "exact_per_read_pred_mapped_rate_species": 0.4,
                "exact_per_read_precision_species": 0.5,
                "exact_per_read_recall_species": 0.6,
                "exact_per_read_f1_species": 0.7,
                "exact_per_read_truth_mapped_rate_genus": 0.8,
                "exact_per_read_pred_mapped_rate_genus": 1.0,
                "exact_per_read_precision_genus": 0.3,
                "exact_per_read_recall_genus": 0.4,
                "exact_per_read_f1_genus": 0.5,
                "completeness_species": 0.25,
                "purity_species": 0.1,
                "l1_norm_species": 1.5,
                "completeness_genus": 0.5,
                "purity_genus": 0.2,
                "l1_norm_genus": 1.0,
                "weighted_unifrac": 0.9,
                "profile_metric_version": METRIC_VERSION,
            },
        ),
        (
            "S2",
            20.0,
            3 * 1024 * 1024,
            {
                "exact_per_read_truth_mapped_rate_species": 0.6,
                "exact_per_read_pred_mapped_rate_species": 0.8,
                "exact_per_read_precision_species": 0.7,
                "exact_per_read_recall_species": 0.8,
                "exact_per_read_f1_species": 0.9,
                "exact_per_read_truth_mapped_rate_genus": 0.4,
                "exact_per_read_pred_mapped_rate_genus": 0.6,
                "exact_per_read_precision_genus": 0.5,
                "exact_per_read_recall_genus": 0.6,
                "exact_per_read_f1_genus": 0.7,
                "completeness_species": 0.75,
                "purity_species": 0.3,
                "l1_norm_species": 0.5,
                "completeness_genus": 0.7,
                "purity_genus": 0.4,
                "l1_norm_genus": 0.2,
                "weighted_unifrac": 0.5,
                "profile_metric_version": METRIC_VERSION,
            },
        ),
    ]

    for sample_id, elapsed_seconds, max_rss_kb, metrics in samples:
        run_dir = runs_root / "kraken2" / f"real.{sample_id}"
        run_dir.mkdir(parents=True)
        (run_dir / "meta.json").write_text(
            json.dumps(
                {
                    "exp": "kraken2",
                    "tool": "kraken2",
                    "dataset": f"real.{sample_id}",
                    "dataset_collection": "real",
                    "sample_id": sample_id,
                    "db": "/tmp/cami_refseq",
                    "db_name": "cami_refseq",
                    "elapsed_seconds": elapsed_seconds,
                    "resource": {"max_rss_kb": max_rss_kb},
                    "return_code": 0,
                }
            )
        )
        (run_dir / "metrics.json").write_text(json.dumps(metrics))

    write_classify_readme(runs_root)
    classify_text = (runs_root / "README.md").read_text()
    assert "## Dataset: real\n" in classify_text
    assert "## Dataset: real.S1" not in classify_text
    assert "## Dataset: real.S2" not in classify_text
    assert "| Tool | DB | Samples | Elapsed (s) | Max RSS (GB) |" in classify_text
    assert "集合聚合（collection aggregate）行的 `Samples` 为已完成 sample 数" in classify_text
    assert "| kraken2 | cami_refseq | 2 | 30 | 3 | 0.4 | 0.6 | 0.6 | 0.7 | 0.8 | 0.6 | 0.8 | 0.4 | 0.5 | 0.6 |" in classify_text

    write_profile_readme(profile_root, runs_root)
    profile_text = (profile_root / "README.md").read_text()
    assert "## Dataset: real\n" in profile_text
    assert "## Dataset: real.S1" not in profile_text
    assert "## Dataset: real.S2" not in profile_text
    assert "| Tool | DB | Samples | Elapsed (s) | Max RSS (GB) |" in profile_text
    assert "profile 指标为 sample 算术平均" in profile_text
    assert "| kraken2 | cami_refseq | 2 | 30 | 3 | 0.5 | 0.2 | 1 | 0.6 | 0.3 | 0.6 | 0.7 |" in profile_text


def test_results_readmes_use_collection_display_dataset(tmp_path: Path):
    runs_root = tmp_path / "runs"
    runs_root.mkdir()

    run_dir = runs_root / "ganon" / "short-batches.batch00"
    run_dir.mkdir(parents=True)
    (run_dir / "meta.json").write_text(
        json.dumps(
            {
                "exp": "ganon",
                "tool": "ganon",
                "dataset": "short-batches.batch00",
                "dataset_collection": "short-batches",
                "display_dataset": "short",
                "sample_id": "batch00",
                "db": "/tmp/cami_refseq",
                "db_name": "cami_refseq",
                "elapsed_seconds": 10.0,
                "resource": {"max_rss_kb": 1024 * 1024},
                "return_code": 0,
            }
        )
    )
    (run_dir / "metrics.json").write_text(
        json.dumps(
            {
                "exact_per_read_truth_mapped_rate_species": 1.0,
                "exact_per_read_pred_mapped_rate_species": 1.0,
                "exact_per_read_precision_species": 0.5,
                "exact_per_read_recall_species": 0.5,
                "exact_per_read_f1_species": 0.5,
            }
        )
    )

    write_classify_readme(runs_root)
    text = (runs_root / "README.md").read_text()

    assert "## Dataset: short\n" in text
    assert "## Dataset: short-batches" not in text
    assert "| ganon2 | cami_refseq | 1 |" in text


def test_main_result_readmes_prefer_supplemental_rerun_over_existing_rows(
    tmp_path: Path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        results_readme_module,
        "SUPPLEMENTAL_MAIN_RESULT_RUNS",
        [Path("results/reruns/current")],
    )

    main_run = Path("results/classify/chimera/ds1")
    main_run.mkdir(parents=True)
    (main_run / "meta.json").write_text(
        json.dumps(
            {
                "tool": "chimera",
                "dataset": "ds1",
                "db_name": "cami_refseq",
                "elapsed_seconds": 100.0,
                "resource": {"max_rss_kb": 1024 * 1024},
                "return_code": 0,
            }
        )
    )
    (main_run / "metrics.json").write_text(
        json.dumps(
            {
                "exact_per_read_truth_mapped_rate_species": 0.1,
                "exact_per_read_pred_mapped_rate_species": 0.2,
                "exact_per_read_precision_species": 0.3,
                "exact_per_read_recall_species": 0.4,
                "exact_per_read_f1_species": 0.5,
                "exact_per_read_truth_mapped_rate_genus": 0.11,
                "exact_per_read_pred_mapped_rate_genus": 0.21,
                "exact_per_read_precision_genus": 0.31,
                "exact_per_read_recall_genus": 0.41,
                "exact_per_read_f1_genus": 0.51,
                "completeness_species": 0.1,
                "purity_species": 0.2,
                "l1_norm_species": 0.3,
                "completeness_genus": 0.4,
                "purity_genus": 0.5,
                "l1_norm_genus": 0.6,
                "weighted_unifrac": 0.7,
                "profile_metric_version": METRIC_VERSION,
            }
        )
    )

    supplemental_run = Path("results/reruns/current/classify/chimera/ds1")
    supplemental_run.mkdir(parents=True)
    (supplemental_run / "meta.json").write_text(
        json.dumps(
            {
                "tool": "chimera",
                "dataset": "ds1",
                "db_name": "cami_refseq",
                "elapsed_seconds": 10.0,
                "resource": {"max_rss_kb": 2 * 1024 * 1024},
                "return_code": 0,
            }
        )
    )
    (supplemental_run / "metrics.json").write_text(
        json.dumps(
            {
                "exact_per_read_truth_mapped_rate_species": 0.6,
                "exact_per_read_pred_mapped_rate_species": 0.7,
                "exact_per_read_precision_species": 0.8,
                "exact_per_read_recall_species": 0.9,
                "exact_per_read_f1_species": 1.0,
                "exact_per_read_truth_mapped_rate_genus": 0.61,
                "exact_per_read_pred_mapped_rate_genus": 0.71,
                "exact_per_read_precision_genus": 0.81,
                "exact_per_read_recall_genus": 0.91,
                "exact_per_read_f1_genus": 0.99,
                "completeness_species": 0.6,
                "purity_species": 0.7,
                "l1_norm_species": 0.8,
                "completeness_genus": 0.61,
                "purity_genus": 0.71,
                "l1_norm_genus": 0.81,
                "weighted_unifrac": 0.91,
                "profile_metric_version": METRIC_VERSION,
            }
        )
    )

    profile_root = Path("results/profile")
    profile_root.mkdir(parents=True)

    write_classify_readme(Path("results/classify"))
    classify_text = Path("results/classify/README.md").read_text()
    assert classify_text.count("| chimera | cami_refseq |") == 1
    assert "| chimera | cami_refseq | 10 | 2 | 0.6 | 0.7 | 0.8 | 0.9 | 1 | 0.61 | 0.71 | 0.81 | 0.91 | 0.99 |" in classify_text
    assert "| chimera | cami_refseq | 100 | 1 | 0.1 |" not in classify_text

    write_profile_readme(profile_root, Path("results/classify"))
    profile_text = Path("results/profile/README.md").read_text()
    assert profile_text.count("| chimera | cami_refseq |") == 1
    assert "| chimera | cami_refseq | 10 | 2 | 0.6 | 0.7 | 0.8 | 0.61 | 0.71 | 0.81 | 0.91 |" in profile_text
    assert "| chimera | cami_refseq | 100 | 1 | 0.1 |" not in profile_text

    main_build = Path("results/builds/chimera/cami_refseq")
    main_build.mkdir(parents=True)
    (main_build / "meta.json").write_text(
        json.dumps(
            {
                "tool": "chimera",
                "db_name": "cami_refseq",
                "elapsed_seconds": 100.0,
                "resource": {"max_rss_kb": 1024},
                "started_at": "old-start",
                "finished_at": "old-finish",
                "return_code": 0,
            }
        )
    )
    supplemental_build = Path("results/reruns/current/builds/chimera/cami_refseq")
    supplemental_build.mkdir(parents=True)
    (supplemental_build / "meta.json").write_text(
        json.dumps(
            {
                "tool": "chimera",
                "db_name": "cami_refseq",
                "elapsed_seconds": 10.0,
                "resource": {"max_rss_kb": 2048},
                "started_at": "new-start",
                "finished_at": "new-finish",
                "return_code": 0,
            }
        )
    )

    write_builds_readme(Path("results/builds"))
    builds_text = Path("results/builds/README.md").read_text()
    assert builds_text.count("| chimera | cami_refseq |") == 1
    assert "| chimera | cami_refseq | 10 | 2048 |  | new-start | new-finish |" in builds_text
    assert "| chimera | cami_refseq | 100 | 1024 |" not in builds_text


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
