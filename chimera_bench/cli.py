from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess

from .config import load_yaml_dir
from .core.runner import Runner
from .core.reporter import write_summary
from .registry import TOOLS


def _executor(cmd, cwd, stdout_path, stderr_path):
    with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
        proc = subprocess.run(cmd, cwd=cwd, stdout=out, stderr=err)
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
    tool_name = exp.get("tool", "chimera")
    tool_cls = TOOLS.get(tool_name)
    tool_config = dict(exp.get("tool_config", {}))
    if tool_name == "chimera":
        tool_config.setdefault("bin", args.chimera_bin)
    if tool_name == "ganon":
        tool_config.setdefault("bin", args.ganon_bin)
        tool_config.setdefault("env", args.ganon_env)

    runner = Runner(Path(args.runs))
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
        run_records.append(
            {
                "exp": meta.get("exp"),
                "tool": meta.get("tool"),
                "dataset": meta.get("dataset"),
                "run_id": meta_path.parent.name,
                "metrics": metrics,
            }
        )
    if args.dataset:
        allowed = set(args.dataset)
        run_records = [r for r in run_records if r.get("dataset") in allowed]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_summary(run_records, out)


def main() -> None:
    p = argparse.ArgumentParser("chimera-bench")
    sub = p.add_subparsers(dest="cmd", required=True)

    run_p = sub.add_parser("run")
    run_p.add_argument("--exp", required=True)
    run_p.add_argument("--config", default="configs")
    run_p.add_argument("--runs", default="runs")
    run_p.add_argument("--chimera-bin", default="Chimera")
    run_p.add_argument("--ganon-bin", default="ganon")
    run_p.add_argument("--ganon-env", default="ganon")
    run_p.add_argument("--dry-run", action="store_true")
    run_p.add_argument("--dataset", action="append", default=[])
    run_p.set_defaults(func=run_cmd)

    report_p = sub.add_parser("report")
    report_p.add_argument("--exp", required=True)
    report_p.add_argument("--runs", default="runs")
    report_p.add_argument("--out", default="reports/summary.tsv")
    report_p.add_argument("--dataset", action="append", default=[])
    report_p.set_defaults(func=report_cmd)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
