from __future__ import annotations

import json
import time
from pathlib import Path

from .evaluator import summarize_classify_tsv, summarize_ganon_tre
from .metrics import evaluate_with_truth
from .results_readme import write_classify_readme, write_profile_readme
from .resources import aggregate_resources, parse_time_log
from ..io.layout import ensure_profile_dirs, ensure_run_dirs


class Runner:
    def __init__(self, runs_root: Path, profile_root: Path | None = None) -> None:
        self.runs_root = runs_root
        self.profile_root = profile_root

    def run(self, *, exp: dict, dataset: dict, tool, executor) -> dict:
        exp_name = exp.get("name", "exp")
        dataset_name = dataset.get("name", "dataset")
        run_dir = ensure_run_dirs(self.runs_root, exp_name, tool.name, dataset_name)
        basename = getattr(tool, "output_basename", tool.name)
        out_prefix = str((run_dir / "outputs" / basename).resolve())
        profile_dir = None
        profile_out_prefix = None
        if self.profile_root is not None:
            profile_dir = ensure_profile_dirs(self.profile_root, exp_name, tool.name, dataset_name)
            profile_out_prefix = str((profile_dir / "outputs" / f"{basename}_abundance").resolve())

        steps = None
        build_steps = getattr(tool, "build_steps", None)
        if callable(build_steps):
            try:
                steps = build_steps(
                    dataset=dataset,
                    exp=exp,
                    out_prefix=out_prefix,
                    profile_dir=profile_dir,
                    profile_out_prefix=profile_out_prefix,
                )
            except TypeError:
                steps = build_steps(dataset=dataset, exp=exp, out_prefix=out_prefix)
        if steps is None:
            try:
                cmd, outputs = tool.build_cmd(
                    dataset=dataset,
                    exp=exp,
                    out_prefix=out_prefix,
                    profile_dir=profile_dir,
                    profile_out_prefix=profile_out_prefix,
                )
            except TypeError:
                cmd, outputs = tool.build_cmd(dataset=dataset, exp=exp, out_prefix=out_prefix)
            steps = [{"name": "run", "cmd": cmd, "outputs": outputs}]

        outputs_all = {}
        step_records = []
        total_start = time.time()
        for idx, step in enumerate(steps):
            name = step.get("name") or f"step{idx + 1}"
            stdout_path = run_dir / "logs" / f"{name}.stdout.log"
            stderr_path = run_dir / "logs" / f"{name}.stderr.log"
            resource_path = run_dir / "logs" / f"{name}.time.log"
            start = time.time()
            rc = executor(
                step["cmd"],
                cwd=run_dir,
                stdout_path=stdout_path,
                stderr_path=stderr_path,
                resource_path=resource_path,
            )
            elapsed = time.time() - start
            resource = parse_time_log(resource_path)
            step_records.append(
                {
                    "name": name,
                    "cmd": step["cmd"],
                    "return_code": rc,
                    "elapsed_seconds": elapsed,
                    "stdout": str(stdout_path),
                    "stderr": str(stderr_path),
                    "resource_log": str(resource_path),
                    "resource": resource,
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
            "profile_dir": str(profile_dir) if profile_dir else None,
            "steps": step_records,
            "return_code": step_records[-1]["return_code"] if step_records else None,
            "elapsed_seconds": total_elapsed,
            "resource": aggregate_resources(step_records),
            "outputs": outputs_all,
        }
        (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))

        metrics = {}
        classify_path_str = outputs_all.get("classify_tsv")
        if classify_path_str:
            classify_path = Path(classify_path_str)
            if classify_path.exists():
                metrics = summarize_classify_tsv(classify_path)
        else:
            tre_path_str = outputs_all.get("report_reads_tre") or outputs_all.get("reads_tre")
            if tre_path_str:
                tre_path = Path(tre_path_str)
                if tre_path.exists():
                    metrics = summarize_ganon_tre(tre_path)

        truth_metrics = evaluate_with_truth(exp, dataset, outputs_all)
        if truth_metrics:
            metrics.update(truth_metrics)
        (run_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))

        write_classify_readme(self.runs_root)
        if self.profile_root is not None:
            write_profile_readme(self.profile_root, self.runs_root)

        return {"run_dir": str(run_dir), "metrics": metrics, "meta": meta}
