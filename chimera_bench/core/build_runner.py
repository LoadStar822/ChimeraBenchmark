from __future__ import annotations

import json
import time
from pathlib import Path

from .resources import aggregate_resources, parse_time_log
from ..io.layout import ensure_run_dirs, make_run_id


class BuildRunner:
    def __init__(self, runs_root: Path) -> None:
        self.runs_root = runs_root

    def run(self, *, build: dict, tool, executor) -> dict:
        build_name = build.get("name", "build")
        run_id = make_run_id(build_name, tool.name, "build")
        run_dir = ensure_run_dirs(self.runs_root, build_name, tool.name, "build", run_id)

        build_steps = getattr(tool, "build_db_steps", None)
        if not callable(build_steps):
            raise ValueError(f"tool {tool.name} does not support build")
        steps = build_steps(build=build, out_dir=str(run_dir))

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
            "build": build_name,
            "tool": tool.name,
            "steps": step_records,
            "return_code": step_records[-1]["return_code"] if step_records else None,
            "elapsed_seconds": total_elapsed,
            "resource": aggregate_resources(step_records),
            "outputs": outputs_all,
        }
        (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))

        return {"run_dir": str(run_dir), "meta": meta}
