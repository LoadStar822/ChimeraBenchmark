from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess

from .config import load_yaml_dir
from .core.runner import Runner
from .core.build_runner import BuildRunner
from .core.reporter import write_summary
from .registry import TOOLS

DEFAULT_THREADS = 32


def _executor(cmd, cwd, stdout_path, stderr_path, resource_path):
    cwd = Path(cwd)
    cwd.mkdir(parents=True, exist_ok=True)
    stdout_path = Path(stdout_path)
    stderr_path = Path(stderr_path)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    timed_cmd = cmd
    if resource_path:
        resource_path = Path(resource_path).resolve()
        resource_path.parent.mkdir(parents=True, exist_ok=True)
        timed_cmd = ["/usr/bin/time", "-v", "-o", str(resource_path)] + cmd
    with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
        proc = subprocess.run(timed_cmd, cwd=cwd, stdout=out, stderr=err)
    return proc.returncode


def _resolve_datasets(exp: dict, datasets: dict, selected: list[str]) -> list[dict]:
    exp_datasets = exp.get("datasets")
    if exp_datasets is None:
        exp_dataset = exp.get("dataset")
        if not exp_dataset:
            raise ValueError("experiment must define dataset or datasets")
        exp_datasets = [exp_dataset]
    if selected:
        exp_datasets = [name for name in exp_datasets if name in set(selected)]
    resolved = []
    for name in exp_datasets:
        if name not in datasets:
            raise KeyError(f"dataset not found: {name}")
        data = dict(datasets[name])
        data["name"] = data.get("name", name)
        resolved.append(data)
    return resolved


def run_cmd(args) -> None:
    cfg_root = Path(args.config)
    datasets = load_yaml_dir(cfg_root / "datasets")
    exps = load_yaml_dir(cfg_root / "experiments")
    exp = dict(exps[args.exp])
    exp["name"] = exp.get("name", args.exp)
    exp["threads"] = DEFAULT_THREADS
    tool_name = exp.get("tool", "chimera")
    tool_cls = TOOLS.get(tool_name)
    tool_config = dict(exp.get("tool_config", {}))
    if tool_name == "chimera":
        tool_config.setdefault("bin", args.chimera_bin)
    if tool_name == "ganon":
        tool_config.setdefault("bin", args.ganon_bin)
        tool_config.setdefault("env", args.ganon_env)
    if tool_name == "sylph":
        tool_config.setdefault("bin", args.sylph_bin)
        tool_config.setdefault("env", args.sylph_env)

    runner = Runner(Path(args.runs), Path(args.profile) if args.profile else None)
    tool = tool_cls(tool_config)

    if args.dry_run:
        Path(args.runs).mkdir(parents=True, exist_ok=True)
        return

    selected = args.dataset or []
    for dataset in _resolve_datasets(exp, datasets, selected):
        runner.run(exp=exp, dataset=dataset, tool=tool, executor=_executor)


def report_cmd(args) -> None:
    runs_root = Path(args.runs) / args.exp
    run_records = []
    for meta_path in runs_root.rglob("meta.json"):
        meta = json.loads(meta_path.read_text())
        metrics_path = meta_path.parent / "metrics.json"
        metrics = json.loads(metrics_path.read_text()) if metrics_path.exists() else {}
        resource = meta.get("resource", {})
        if meta.get("elapsed_seconds") is not None:
            metrics["run_elapsed_seconds"] = meta.get("elapsed_seconds")
        for key, value in resource.items():
            metrics[f"resource_{key}"] = value
        run_records.append(
            {
                "exp": meta.get("exp"),
                "tool": meta.get("tool"),
                "dataset": meta.get("dataset"),
                "metrics": metrics,
            }
        )
    if args.dataset:
        allowed = set(args.dataset)
        run_records = [r for r in run_records if r.get("dataset") in allowed]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_summary(run_records, out)


def build_cmd(args) -> None:
    cfg_root = Path(args.config)
    builds = load_yaml_dir(cfg_root / "build")
    build = dict(builds[args.build])
    build["name"] = build.get("name", args.build)
    build["threads"] = DEFAULT_THREADS
    tool_name = build.get("tool", "ganon")
    tool_cls = TOOLS.get(tool_name)
    tool_config = dict(build.get("tool_config", {}))
    if tool_name == "ganon":
        tool_config.setdefault("bin", args.ganon_bin)
        tool_config.setdefault("env", args.ganon_env)
    if tool_name == "sylph":
        tool_config.setdefault("bin", args.sylph_bin)
        tool_config.setdefault("env", args.sylph_env)

    runner = BuildRunner(Path(args.runs))
    tool = tool_cls(tool_config)

    if args.dry_run:
        Path(args.runs).mkdir(parents=True, exist_ok=True)
        return

    runner.run(build=build, tool=tool, executor=_executor)


def main() -> None:
    p = argparse.ArgumentParser("chimera-bench")
    sub = p.add_subparsers(dest="cmd", required=True)

    run_p = sub.add_parser("run")
    run_p.add_argument("--exp", required=True)
    run_p.add_argument("--config", default="configs")
    run_p.add_argument("--runs", default="results/classify")
    run_p.add_argument("--profile", default="results/profile")
    run_p.add_argument("--chimera-bin", default="Chimera")
    run_p.add_argument("--ganon-bin", default="ganon")
    run_p.add_argument("--ganon-env", default="ganon")
    run_p.add_argument("--sylph-bin", default="sylph")
    run_p.add_argument("--sylph-env", default="sylph")
    run_p.add_argument("--dry-run", action="store_true")
    run_p.add_argument("--dataset", action="append", default=[])
    run_p.set_defaults(func=run_cmd)

    report_p = sub.add_parser("report")
    report_p.add_argument("--exp", required=True)
    report_p.add_argument("--runs", default="results/classify")
    report_p.add_argument("--out", default="results/reports/summary.tsv")
    report_p.add_argument("--dataset", action="append", default=[])
    report_p.set_defaults(func=report_cmd)

    build_p = sub.add_parser("build")
    build_p.add_argument("--build", required=True)
    build_p.add_argument("--config", default="configs")
    build_p.add_argument("--runs", default="results/builds")
    build_p.add_argument("--ganon-bin", default="ganon")
    build_p.add_argument("--ganon-env", default="ganon")
    build_p.add_argument("--sylph-bin", default="sylph")
    build_p.add_argument("--sylph-env", default="sylph")
    build_p.add_argument("--dry-run", action="store_true")
    build_p.set_defaults(func=build_cmd)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
