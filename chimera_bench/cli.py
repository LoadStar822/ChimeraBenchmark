from __future__ import annotations

import argparse
import json
import resource
from pathlib import Path
import subprocess
import sys

from .catalog import write_catalog_outputs
from .config import expand_dataset_config, load_yaml_dir
from .dataset_prepare import prepare_dataset_inputs
from .paper_freeze import write_paper_tables
from .core.runner import Runner, build_run_metrics
from .core.build_runner import BuildRunner
from .core.reporter import write_summary
from .core.results_readme import write_classify_readme, write_profile_readme
from .registry import TOOLS

DEFAULT_THREADS = 32


def _make_executor(*, max_file_bytes: int | None = None):
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

        def _apply_limits():
            if max_file_bytes is not None:
                resource.setrlimit(resource.RLIMIT_FSIZE, (max_file_bytes, max_file_bytes))

        preexec_fn = _apply_limits if max_file_bytes is not None else None
        with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
            proc = subprocess.run(timed_cmd, cwd=cwd, stdout=out, stderr=err, preexec_fn=preexec_fn)
        return proc.returncode

    return _executor


def _resolve_datasets(exp: dict, datasets: dict, selected: list[str]) -> list[dict]:
    exp_datasets = exp.get("datasets")
    if exp_datasets is None:
        exp_dataset = exp.get("dataset")
        if not exp_dataset:
            raise ValueError("experiment must define dataset or datasets")
        exp_datasets = [exp_dataset]
    selected_set = set(selected)
    resolved = []
    for name in exp_datasets:
        if name not in datasets:
            raise KeyError(f"dataset not found: {name}")
        expanded = expand_dataset_config(name, datasets[name])
        if selected_set and name not in selected_set:
            expanded = [
                data
                for data in expanded
                if data.get("name") in selected_set
                or data.get("sample_id") in selected_set
                or data.get("dataset_collection") in selected_set
                or data.get("display_dataset") in selected_set
            ]
        resolved.extend(expanded)
    return resolved


def run_cmd(args) -> None:
    cfg_root = Path(args.config)
    datasets = load_yaml_dir(cfg_root / "datasets")
    exps = load_yaml_dir(cfg_root / "experiments")
    exp = dict(exps[args.exp])
    exp["name"] = exp.get("name", args.exp)
    exp.setdefault("threads", DEFAULT_THREADS)
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
    executor = _make_executor()

    selected = args.dataset or []
    resolved_datasets = _resolve_datasets(exp, datasets, selected)
    if args.dry_run:
        Path(args.runs).mkdir(parents=True, exist_ok=True)
        return

    failed_datasets = []
    for dataset in resolved_datasets:
        dataset = prepare_dataset_inputs(dataset)
        result = runner.run(exp=exp, dataset=dataset, tool=tool, executor=executor)
        meta = (result or {}).get("meta") if isinstance(result, dict) else None
        if isinstance(meta, dict) and meta.get("return_code") not in {None, 0}:
            failed_datasets.append((dataset.get("name", "dataset"), meta.get("return_code")))

    if failed_datasets:
        for dataset_name, return_code in failed_datasets:
            print(f"failed dataset: {dataset_name} return_code={return_code}", file=sys.stderr)
        raise SystemExit(1)


def catalog_cmd(args) -> None:
    write_catalog_outputs(
        config_root=Path(args.config),
        results_root=Path(args.results_root),
        resources_root=Path(getattr(args, "resources_root", "resources")),
        progress=True,
    )


def paper_freeze_cmd(args) -> None:
    counts = write_paper_tables(
        config_root=Path(args.config),
        results_root=Path(args.results_root),
    )
    print(json.dumps(counts, sort_keys=True))


def _collect_summary_records(runs_root: Path, exp_name: str, selected: list[str]) -> list[dict]:
    exp_root = runs_root / exp_name
    run_records = []
    for meta_path in exp_root.rglob("meta.json"):
        meta = json.loads(meta_path.read_text())
        if meta.get("return_code") not in {None, 0}:
            continue
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
                "dataset_collection": meta.get("dataset_collection"),
                "display_dataset": meta.get("display_dataset"),
                "sample_id": meta.get("sample_id"),
                "metrics": metrics,
            }
        )
    if selected:
        allowed = set(selected)
        run_records = [
            r
            for r in run_records
            if r.get("dataset") in allowed
            or r.get("dataset_collection") in allowed
            or r.get("display_dataset") in allowed
            or r.get("sample_id") in allowed
        ]
    return run_records


def report_cmd(args) -> None:
    run_records = _collect_summary_records(Path(args.runs), args.exp, args.dataset or [])
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_summary(run_records, out)


def recompute_cmd(args) -> None:
    cfg_root = Path(args.config)
    datasets = load_yaml_dir(cfg_root / "datasets")
    exps = load_yaml_dir(cfg_root / "experiments")
    exp = dict(exps[args.exp])
    exp["name"] = exp.get("name", args.exp)
    selected = args.dataset or []
    runs_root = Path(args.runs)
    exp_root = runs_root / args.exp

    for dataset in _resolve_datasets(exp, datasets, selected):
        dataset_name = dataset.get("name", "dataset")
        run_dir = exp_root / dataset_name
        meta_path = run_dir / "meta.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        if meta.get("return_code") not in {None, 0}:
            continue
        metrics = build_run_metrics(exp, dataset, meta.get("outputs") or {})
        (run_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))

    write_classify_readme(runs_root)
    if args.profile:
        write_profile_readme(Path(args.profile), runs_root)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_summary(_collect_summary_records(runs_root, args.exp, selected), out)


def build_cmd(args) -> None:
    cfg_root = Path(args.config)
    builds = load_yaml_dir(cfg_root / "build")
    build = dict(builds[args.build])
    build["name"] = build.get("name", args.build)
    build.setdefault("threads", DEFAULT_THREADS)
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
    executor = _make_executor()

    if args.dry_run:
        Path(args.runs).mkdir(parents=True, exist_ok=True)
        return

    result = runner.run(build=build, tool=tool, executor=executor)
    meta = (result or {}).get("meta") if isinstance(result, dict) else None
    if isinstance(meta, dict) and meta.get("return_code") not in {None, 0}:
        print(
            f"failed build: {build.get('name', args.build)} return_code={meta.get('return_code')}",
            file=sys.stderr,
        )
        raise SystemExit(1)


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
    report_p.add_argument("--out", default="resources/reports/summary.tsv")
    report_p.add_argument("--dataset", action="append", default=[])
    report_p.set_defaults(func=report_cmd)

    recompute_p = sub.add_parser("recompute")
    recompute_p.add_argument("--exp", required=True)
    recompute_p.add_argument("--config", default="configs")
    recompute_p.add_argument("--runs", default="results/classify")
    recompute_p.add_argument("--profile", default="results/profile")
    recompute_p.add_argument("--out", default="resources/reports/summary.tsv")
    recompute_p.add_argument("--dataset", action="append", default=[])
    recompute_p.set_defaults(func=recompute_cmd)

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

    catalog_p = sub.add_parser("catalog")
    catalog_p.add_argument("--config", default="configs")
    catalog_p.add_argument("--results-root", default="results")
    catalog_p.add_argument("--resources-root", default="resources")
    catalog_p.set_defaults(func=catalog_cmd)

    paper_p = sub.add_parser("paper-freeze")
    paper_p.add_argument("--config", default="configs")
    paper_p.add_argument("--results-root", default="results")
    paper_p.set_defaults(func=paper_freeze_cmd)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
