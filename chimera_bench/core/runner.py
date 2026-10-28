from __future__ import annotations

import json
import time
from datetime import datetime
from pathlib import Path

from .evaluator import summarize_classify_tsv, summarize_ganon_tre
from .metrics import evaluate_with_truth
from .resources import aggregate_resources, parse_time_log
from ..io.layout import claim_experiment_root, configuration_digest, ensure_run_dirs, require_fresh_run_dir, run_dir_path, write_json_within


def build_run_metrics(exp: dict, dataset: dict, outputs: dict) -> dict:
    metrics = {}
    classify_path_str = outputs.get("classify_tsv")
    if classify_path_str:
        classify_path = Path(classify_path_str)
        if classify_path.exists():
            metrics = summarize_classify_tsv(classify_path)
    else:
        tre_path_str = outputs.get("report_reads_tre") or outputs.get("reads_tre")
        if tre_path_str:
            tre_path = Path(tre_path_str)
            if tre_path.exists():
                metrics = summarize_ganon_tre(tre_path)

    truth_metrics = evaluate_with_truth(exp, dataset, outputs)
    if truth_metrics:
        metrics.update(truth_metrics)
    metrics["evaluation"] = {"config": dict(exp), "dataset": dict(dataset)}
    return metrics


class Runner:
    """Executes one experiment inside an explicit experiment root.

    The root is supplied by the caller (``--runs``) and never has a global
    default. Everything a run produces -- logs, native tool output and the
    tool's own profile output -- lives inside that run's directory. Runs are
    never overwritten: a target directory that already holds artifacts is
    rejected before the executor is invoked.
    """

    def __init__(self, runs_root: Path) -> None:
        self.runs_root = runs_root

    def run_dir_for(self, *, tool_name: str, dataset_name: str) -> Path:
        return run_dir_path(self.runs_root, tool_name, dataset_name)

    def require_fresh(self, *, tool_name: str, dataset_name: str) -> None:
        require_fresh_run_dir(self.run_dir_for(tool_name=tool_name, dataset_name=dataset_name))

    def run(self, *, exp: dict, dataset: dict, tool, executor) -> dict:
        exp_name = exp.get("name", "exp")
        dataset_name = dataset.get("name", "dataset")
        run_dir = self.run_dir_for(tool_name=tool.name, dataset_name=dataset_name)
        require_fresh_run_dir(run_dir)
        claim_experiment_root(
            self.runs_root,
            kind="runs",
            identity={
                "experiment": exp_name,
                "tool": tool.name,
                "config_sha256": configuration_digest(exp, getattr(tool, "config", None)),
            },
        )
        run_dir = ensure_run_dirs(self.runs_root, tool.name, dataset_name)

        basename = getattr(tool, "output_basename", tool.name)
        outputs_dir = run_dir / "outputs"
        out_prefix = str((outputs_dir / basename).resolve())
        # Profile output of the same run is co-located with its other outputs.
        profile_out_prefix = f"{out_prefix}_abundance"

        steps = None
        build_steps = getattr(tool, "build_steps", None)
        if callable(build_steps):
            try:
                steps = build_steps(
                    dataset=dataset,
                    exp=exp,
                    out_prefix=out_prefix,
                    profile_dir=str(outputs_dir.resolve()),
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
                    profile_dir=str(outputs_dir.resolve()),
                    profile_out_prefix=profile_out_prefix,
                )
            except TypeError:
                cmd, outputs = tool.build_cmd(dataset=dataset, exp=exp, out_prefix=out_prefix)
            steps = [{"name": "run", "cmd": cmd, "outputs": outputs}]

        outputs_all = {}
        step_records = []
        started_at = datetime.now().astimezone()
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
            if rc == 0:
                outputs_all.update(step.get("outputs", {}))
            if rc != 0:
                break

        total_elapsed = time.time() - total_start
        finished_at = datetime.now().astimezone()
        db_path = exp.get("db") or exp.get("db_prefix")
        db_name = None
        if db_path:
            db_name = Path(db_path).name
        succeeded = bool(step_records) and all(
            record["return_code"] == 0 for record in step_records
        )
        meta = {
            "exp": exp_name,
            "dataset": dataset_name,
            "dataset_collection": dataset.get("dataset_collection"),
            "display_dataset": dataset.get("display_dataset"),
            "sample_id": dataset.get("sample_id"),
            "sample_ids": dataset.get("sample_ids"),
            "batch_id": dataset.get("batch_id"),
            "tool": tool.name,
            "db": db_path,
            "db_name": db_name,
            "status": "evaluating" if succeeded else "failed",
            "started_at": started_at.isoformat(timespec="seconds"),
            "finished_at": finished_at.isoformat(timespec="seconds"),
            "steps": step_records,
            "return_code": step_records[-1]["return_code"] if step_records else None,
            "elapsed_seconds": total_elapsed,
            "resource": aggregate_resources(step_records),
            "outputs": outputs_all,
            "config": dict(exp),
            "dataset_config": dict(dataset),
            "tool_config": dict(getattr(tool, "config", None) or {}),
        }
        write_json_within(self.runs_root, run_dir / "meta.json", meta)

        metrics = None
        if succeeded:
            try:
                metrics = build_run_metrics(exp, dataset, outputs_all)
                write_json_within(self.runs_root, run_dir / "metrics.json", metrics)
            except Exception as exc:
                meta["status"] = "evaluation_failed"
                meta["evaluation_error"] = str(exc)
                write_json_within(self.runs_root, run_dir / "meta.json", meta)
                raise
            meta["status"] = "succeeded"
            write_json_within(self.runs_root, run_dir / "meta.json", meta)

        return {"run_dir": str(run_dir), "metrics": metrics, "meta": meta}
