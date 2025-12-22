from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess

from .config import load_yaml_dir
from .core.runner import Runner
from .core.reporter import write_summary
from .tools.chimera import ChimeraTool


def _executor(cmd, cwd, stdout_path, stderr_path):
    with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
        proc = subprocess.run(cmd, cwd=cwd, stdout=out, stderr=err)
    return proc.returncode


def run_cmd(args) -> None:
    cfg_root = Path(args.config)
    datasets = load_yaml_dir(cfg_root / "datasets")
    exps = load_yaml_dir(cfg_root / "experiments")
    exp = exps[args.exp]
    dataset = datasets[exp["dataset"]]
    exp["name"] = exp.get("name", args.exp)
    dataset["name"] = dataset.get("name", exp["dataset"])

    runner = Runner(Path(args.runs))
    tool = ChimeraTool({"bin": args.chimera_bin})

    if args.dry_run:
        Path(args.runs).mkdir(parents=True, exist_ok=True)
        return

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
    run_p.add_argument("--dry-run", action="store_true")
    run_p.set_defaults(func=run_cmd)

    report_p = sub.add_parser("report")
    report_p.add_argument("--exp", required=True)
    report_p.add_argument("--runs", default="runs")
    report_p.add_argument("--out", default="reports/summary.tsv")
    report_p.set_defaults(func=report_cmd)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
