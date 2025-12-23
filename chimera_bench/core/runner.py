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
        basename = getattr(tool, "output_basename", tool.name)
        out_prefix = str(run_dir / "outputs" / basename)

        steps = None
        build_steps = getattr(tool, "build_steps", None)
        if callable(build_steps):
            steps = build_steps(dataset=dataset, exp=exp, out_prefix=out_prefix)
        if steps is None:
            cmd, outputs = tool.build_cmd(dataset=dataset, exp=exp, out_prefix=out_prefix)
            steps = [{"name": "run", "cmd": cmd, "outputs": outputs}]

        outputs_all = {}
        step_records = []
        total_start = time.time()
        for idx, step in enumerate(steps):
            name = step.get("name") or f"step{idx + 1}"
            stdout_path = run_dir / "logs" / f"{name}.stdout.log"
            stderr_path = run_dir / "logs" / f"{name}.stderr.log"
            start = time.time()
            rc = executor(step["cmd"], cwd=run_dir, stdout_path=stdout_path, stderr_path=stderr_path)
            elapsed = time.time() - start
            step_records.append(
                {
                    "name": name,
                    "cmd": step["cmd"],
                    "return_code": rc,
                    "elapsed_seconds": elapsed,
                    "stdout": str(stdout_path),
                    "stderr": str(stderr_path),
                }
            )
            outputs_all.update(step.get("outputs", {}))
            if rc != 0:
                break

        total_elapsed = time.time() - total_start
        meta = {
            "exp": exp_name,
            "dataset": dataset_name,
            "tool": tool.name,
            "steps": step_records,
            "return_code": step_records[-1]["return_code"] if step_records else None,
            "elapsed_seconds": total_elapsed,
            "outputs": outputs_all,
        }
        (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))

        metrics = {}
        classify_path_str = outputs_all.get("classify_tsv")
        if classify_path_str:
            classify_path = Path(classify_path_str)
            if classify_path.exists():
                metrics = summarize_classify_tsv(classify_path)
        (run_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))

        return {"run_dir": str(run_dir), "metrics": metrics, "meta": meta}
