from __future__ import annotations

import json
from pathlib import Path

from chimera_bench.core.results_readme import write_classify_readme, write_profile_readme


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
            "| Tool | DB | Elapsed (s) | Max RSS (GB) | Presence Precision (species, UNK) | Presence Recall (species, UNK) | Presence F1 (species, UNK) | L1 (species, UNK) | TV (species, UNK) | Bray-Curtis (species, UNK) | Presence Precision (genus, UNK) | Presence Recall (genus, UNK) | Presence F1 (genus, UNK) | L1 (genus, UNK) | TV (genus, UNK) | Bray-Curtis (genus, UNK) |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
            # ganon row should be preserved even if no run exists on disk
            "| ganon | cami_refseq | 1 | 2 | 1 | 1 | 1 | 0.5 | 0.25 | 0.25 | 1 | 1 | 1 | 0.5 | 0.25 | 0.25 |",
            # sylph row will be updated from on-disk run
            "| sylph | cami_refseq | 9 | 9 | 0.1 | 0.2 | 0.3 | 1 | 0.005 | 0.005 | 0.1 | 0.2 | 0.3 | 1 | 0.005 | 0.005 |",
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
        "presence_precision_species": 0.5,
        "presence_recall_species": 0.25,
        "presence_f1_species": 0.3333333,
        "abundance_l1_species": 20.0,
        "abundance_tv_species": 0.1,
        "abundance_bc_species": 0.1,
        "presence_precision_genus": 0.6,
        "presence_recall_genus": 0.3,
        "presence_f1_genus": 0.4,
        "abundance_l1_genus": 40.0,
        "abundance_tv_genus": 0.2,
        "abundance_bc_genus": 0.2,
    }
    (run_dir / "metrics.json").write_text(json.dumps(metrics))

    write_profile_readme(profile_root, runs_root)

    text = (profile_root / "README.md").read_text()
    assert "| ganon | cami_refseq |" not in text
    assert "### Abundance Metrics (UNK)" not in text
    assert "| Tool | DB | Elapsed (s) | Max RSS (GB) | Presence Precision (species) |" in text
    assert "| sylph | cami_refseq | 10 | 1 | 0.5 | 0.25 | 0.333333 | 20 | 0.1 | 0.1 |" in text


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
    (run_dir / "metrics.json").write_text(json.dumps({"presence_precision_species": 1.0}))

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
                "presence_precision_species": 0.5,
                "presence_recall_species": 0.25,
                "presence_f1_species": 0.3333333,
                "abundance_l1_species": 20.0,
                "abundance_tv_species": 0.1,
                "abundance_bc_species": 0.1,
                "presence_precision_genus": 0.6,
                "presence_recall_genus": 0.3,
                "presence_f1_genus": 0.4,
                "abundance_l1_genus": 40.0,
                "abundance_tv_genus": 0.2,
                "abundance_bc_genus": 0.2,
            }
        )
    )

    write_classify_readme(runs_root)
    classify_text = (runs_root / "README.md").read_text()
    assert "| ganon2 | cami_refseq |" in classify_text
    assert "| ganon | cami_refseq |" not in classify_text

    write_profile_readme(profile_root, runs_root)
    profile_text = (profile_root / "README.md").read_text()
    assert "| ganon2 | cami_refseq |" in profile_text
    assert "| ganon | cami_refseq |" not in profile_text
