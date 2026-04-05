from __future__ import annotations

import json
from pathlib import Path

from chimera_bench.core.results_readme import write_classify_readme, write_profile_readme


def test_classify_and_profile_readmes_skip_failed_runs(tmp_path: Path):
    runs_root = tmp_path / "runs"
    profile_root = tmp_path / "profile"
    failed_dir = runs_root / "chimera" / "ds1"
    failed_dir.mkdir(parents=True)
    profile_root.mkdir()

    (failed_dir / "meta.json").write_text(
        json.dumps(
            {
                "exp": "chimera",
                "tool": "chimera",
                "dataset": "ds1",
                "db_name": "cami_refseq",
                "return_code": 137,
                "elapsed_seconds": 10.0,
                "resource": {"max_rss_kb": 1024},
            }
        )
    )
    (failed_dir / "metrics.json").write_text(
        json.dumps(
            {
                "per_read_precision_species": 0.9,
                "per_read_recall_species": 0.8,
                "per_read_f1_species": 0.85,
                "presence_precision_species": 0.7,
                "presence_recall_species": 0.6,
                "presence_f1_species": 0.65,
            }
        )
    )

    write_classify_readme(runs_root)
    write_profile_readme(profile_root, runs_root)

    classify_text = (runs_root / "README.md").read_text()
    profile_text = (profile_root / "README.md").read_text()
    assert "## Dataset:" not in classify_text
    assert "## Dataset:" not in profile_text
    assert "| chimera | cami_refseq |" not in classify_text
    assert "| chimera | cami_refseq |" not in profile_text
