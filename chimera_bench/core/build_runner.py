from __future__ import annotations

import json
import time
from datetime import datetime
from pathlib import Path

from .resources import aggregate_resources, parse_time_log
from ..io.layout import build_dir_path, claim_experiment_root, configuration_digest, ensure_build_dirs, require_fresh_run_dir


class BuildRunner:
    """Runs one database build inside an explicit build experiment root.

    The root is supplied by the caller (``--runs``). A build directory that
    already holds artifacts is rejected before the executor is invoked, and
    nothing outside the build directory (in particular no global README) is
    written by a build.
    """

    def __init__(self, runs_root: Path) -> None:
        self.runs_root = runs_root

    @staticmethod
    def resolve_db_name(build: dict) -> str:
        db_name = build.get("db_name")
        if db_name:
            return str(db_name)
        db_prefix = build.get("db_prefix") or build.get("db")
        if db_prefix:
            return Path(db_prefix).name
        return str(build.get("name", "build"))

    def build_dir_for(self, *, tool_name: str, build: dict) -> Path:
        return build_dir_path(self.runs_root, tool_name, self.resolve_db_name(build))

    def run(self, *, build: dict, tool, executor) -> dict:
        build_name = build.get("name", "build")
        db_prefix = build.get("db_prefix") or build.get("db")
        db_name = self.resolve_db_name(build)
        run_dir = self.build_dir_for(tool_name=tool.name, build=build)
        require_fresh_run_dir(run_dir)
        claim_experiment_root(
            self.runs_root,
            kind="builds",
            identity={
                "experiment": build_name,
                "tool": tool.name,
                "config_sha256": configuration_digest(build, getattr(tool, "config", None)),
            },
        )
        run_dir = ensure_build_dirs(self.runs_root, tool.name, db_name)
        execution_build = dict(build)
        db_prefix = str((run_dir / "DB" / db_name).resolve())
        execution_build["db_prefix"] = db_prefix
        if "db" in execution_build:
            execution_build["db"] = db_prefix

        build_steps = getattr(tool, "build_db_steps", None)
        if not callable(build_steps):
            raise ValueError(f"tool {tool.name} does not support build")
        steps = build_steps(build=execution_build, out_dir=str(run_dir.resolve()))

        outputs_all = {}
        step_records = []
        started_at = datetime.now().astimezone()
        total_start = time.time()
        for idx, step in enumerate(steps):
            name = step.get("name") or f"step{idx + 1}"
            stdout_path = run_dir / "logs" / f"{name}.stdout.log"
            stderr_path = run_dir / "logs" / f"{name}.stderr.log"
            resource_path = run_dir / "DB" / f"{name}.time.log"
            resource_path.parent.mkdir(parents=True, exist_ok=True)
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
        succeeded = bool(step_records) and all(
            record["return_code"] == 0 for record in step_records
        )
        meta = {
            "build": build_name,
            "tool": tool.name,
            "db_name": db_name,
            "db_prefix": db_prefix,
            "status": "succeeded" if succeeded else "failed",
            "started_at": started_at.isoformat(timespec="seconds"),
            "finished_at": finished_at.isoformat(timespec="seconds"),
            "steps": step_records,
            "return_code": step_records[-1]["return_code"] if step_records else None,
            "elapsed_seconds": total_elapsed,
            "resource": aggregate_resources(step_records),
            "outputs": outputs_all,
            "config": execution_build,
            "requested_config": dict(build),
            "tool_config": dict(getattr(tool, "config", None) or {}),
        }
        (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))

        return {"run_dir": str(run_dir), "meta": meta}
