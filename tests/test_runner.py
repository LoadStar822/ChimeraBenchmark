from __future__ import annotations

import json
from pathlib import Path

from chimera_bench.core.runner import Runner


class DummyTool:
    name = "chimera"
    output_basename = "ChimeraClassify"

    def __init__(self, steps):
        self._steps = steps

    def build_steps(self, **_kwargs):
        return self._steps


def test_runner_ignores_failed_step_outputs_and_clears_stale_metrics(tmp_path: Path):
    runs_root = tmp_path / "runs"
    runner = Runner(runs_root)

    run_dir = runs_root / "chimera" / "ds1"
    outputs_dir = run_dir / "outputs"
    outputs_dir.mkdir(parents=True)
    stale_tsv = outputs_dir / "ChimeraClassify.tsv"
    stale_tsv.write_text("read1\t123\n")
    (run_dir / "metrics.json").write_text(json.dumps({"per_read_precision_species": 0.99}))

    tool = DummyTool(
        [
            {
                "name": "classify",
                "cmd": ["fake-classify"],
                "outputs": {"classify_tsv": str(stale_tsv)},
            }
        ]
    )

    def executor(cmd, **_kwargs):
        assert cmd == ["fake-classify"]
        return 137

    result = runner.run(exp={"name": "chimera", "db": "/tmp/db"}, dataset={"name": "ds1"}, tool=tool, executor=executor)

    assert result["meta"]["return_code"] == 137
    assert result["meta"]["outputs"] == {}
    metrics = json.loads((run_dir / "metrics.json").read_text())
    assert metrics == {}
    readme_text = (runs_root / "README.md").read_text()
    assert "## Dataset:" not in readme_text


def test_runner_keeps_only_successful_step_outputs(tmp_path: Path):
    runs_root = tmp_path / "runs"
    runner = Runner(runs_root)

    run_dir = runs_root / "chimera" / "ds2"
    outputs_dir = run_dir / "outputs"
    outputs_dir.mkdir(parents=True)
    classify_tsv = outputs_dir / "ChimeraClassify.tsv"
    profile_tsv = outputs_dir / "ChimeraClassify_abundance.tsv"

    tool = DummyTool(
        [
            {
                "name": "classify",
                "cmd": ["fake-classify"],
                "outputs": {"classify_tsv": str(classify_tsv)},
            },
            {
                "name": "profile",
                "cmd": ["fake-profile"],
                "outputs": {"chimera_profile_tsv": str(profile_tsv)},
            },
        ]
    )

    def executor(cmd, **_kwargs):
        if cmd == ["fake-classify"]:
            classify_tsv.write_text("read1\t123\nread2\tunclassified\n")
            return 0
        assert cmd == ["fake-profile"]
        return 1

    result = runner.run(exp={"name": "chimera", "db": "/tmp/db"}, dataset={"name": "ds2"}, tool=tool, executor=executor)

    assert result["meta"]["return_code"] == 1
    assert result["meta"]["outputs"] == {"classify_tsv": str(classify_tsv)}
    metrics = json.loads((run_dir / "metrics.json").read_text())
    assert metrics["total_reads"] == 2
    assert metrics["classified_reads"] == 1
