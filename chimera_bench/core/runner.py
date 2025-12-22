from __future__ import annotations

import json
import time
from pathlib import Path

from .evaluator import summarize_classify_tsv
from ..io.layout import ensure_run_dirs, make_run_id


class Runner:
    def __init__(self, runs_root: Path) -> None:
        self.runs_root = runs_root

    def run(self, *, exp: dict, dataset: dict, tool, executor) -> dict:
        exp_name = exp.get("name", "exp")
        dataset_name = dataset.get("name", "dataset")
        run_id = make_run_id(exp_name, tool.name, dataset_name)
        run_dir = ensure_run_dirs(self.runs_root, exp_name, tool.name, dataset_name, run_id)
        out_prefix = str(run_dir / "outputs" / "ChimeraClassify")

        cmd, outputs = tool.build_cmd(dataset=dataset, exp=exp, out_prefix=out_prefix)
        stdout_path = run_dir / "logs" / "stdout.log"
        stderr_path = run_dir / "logs" / "stderr.log"

        start = time.time()
        rc = executor(cmd, cwd=run_dir, stdout_path=stdout_path, stderr_path=stderr_path)
        elapsed = time.time() - start

        meta = {
            "exp": exp_name,
            "dataset": dataset_name,
            "tool": tool.name,
            "cmd": cmd,
            "return_code": rc,
            "elapsed_seconds": elapsed,
            "outputs": outputs,
        }
        (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))

        metrics = {}
        classify_path = Path(outputs["classify_tsv"])
        if classify_path.exists():
            metrics = summarize_classify_tsv(classify_path)
        (run_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))

        return {"run_dir": str(run_dir), "metrics": metrics, "meta": meta}
