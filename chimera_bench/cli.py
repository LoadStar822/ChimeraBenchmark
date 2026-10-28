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
from .core.runner import Runner, build_run_metrics
from .core.build_runner import BuildRunner
from .core.reporter import write_summary
from .io.layout import (
    ExperimentRootError,
    ReusedRunDirectoryError,
    claim_experiment_root,
    configuration_digest,
    require_experiment_root,
    write_json_within,
)
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

    runs_root = Path(args.runs)
    runner = Runner(runs_root)
    tool = tool_cls(tool_config)
    executor = _make_executor()

    selected = args.dataset or []
    resolved_datasets = _resolve_datasets(exp, datasets, selected)
    if not resolved_datasets:
        raise SystemExit(f"no datasets selected for experiment {args.exp}")

    if args.dry_run:
        runs_root.mkdir(parents=True, exist_ok=True)
        return

    try:
        claim_experiment_root(
            runs_root,
            kind="runs",
            identity={
                "experiment": exp["name"],
                "tool": tool.name,
                "config_sha256": configuration_digest(exp, getattr(tool, "config", None)),
            },
        )
        for dataset in resolved_datasets:
            runner.require_fresh(
                tool_name=tool.name, dataset_name=dataset.get("name", "dataset")
            )
    except (ExperimentRootError, ReusedRunDirectoryError) as exc:
        raise SystemExit(str(exc)) from exc

    failed_datasets = []
    for dataset in resolved_datasets:
        dataset = prepare_dataset_inputs(dataset)
        result = runner.run(exp=exp, dataset=dataset, tool=tool, executor=executor)
        meta = (result or {}).get("meta") if isinstance(result, dict) else None
        if isinstance(meta, dict) and meta.get("status") != "succeeded":
            failed_datasets.append((dataset.get("name", "dataset"), meta.get("return_code")))

    if failed_datasets:
        for dataset_name, return_code in failed_datasets:
            print(f"failed dataset: {dataset_name} return_code={return_code}", file=sys.stderr)
        raise SystemExit(1)


def catalog_cmd(args) -> None:
    write_catalog_outputs(
        config_root=Path(args.config),
        resources_root=Path(getattr(args, "resources_root", "resources")),
        progress=True,
    )


def _run_completed(meta: dict) -> bool:
    return meta.get("status") == "succeeded" and meta.get("return_code") == 0


def _collect_summary_records(exp_root: Path, selected: list[str]) -> list[dict]:
    run_records = []
    for meta_path in sorted(exp_root.rglob("meta.json")):
        try:
            meta = json.loads(meta_path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(meta, dict) or not _run_completed(meta):
            continue
        metrics_path = meta_path.parent / "metrics.json"
        try:
            metrics = json.loads(metrics_path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(metrics, dict):
            continue
        resource = meta.get("resource") or {}
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


def _validated_export_path(exp_root: Path, out: Path) -> Path:
    resolved_root = exp_root.resolve()
    resolved_out = out.resolve()
    if not resolved_out.is_relative_to(resolved_root):
        raise SystemExit(
            f"report output must stay inside the experiment root {resolved_root}: {resolved_out}"
        )
    if resolved_out.parent != resolved_root or resolved_out.suffix.lower() != ".tsv":
        raise SystemExit(
            f"report output must be a top-level .tsv export, not a run artifact: {resolved_out}"
        )
    if resolved_out.is_dir():
        raise SystemExit(f"report output is a directory: {resolved_out}")
    return resolved_out


def report_cmd(args) -> None:
    exp_root = Path(args.runs)
    out = _validated_export_path(exp_root, Path(args.out))
    try:
        require_experiment_root(exp_root, kind="runs")
    except ExperimentRootError as exc:
        raise SystemExit(str(exc)) from exc
    run_records = _collect_summary_records(exp_root, args.dataset or [])
    if not run_records:
        raise SystemExit(f"no completed runs matched under {exp_root}")
    out.parent.mkdir(parents=True, exist_ok=True)
    write_summary(run_records, out)
    print(out)


def recompute_cmd(args) -> None:
    cfg_root = Path(args.config)
    exp_root = Path(args.runs)
    datasets = load_yaml_dir(cfg_root / "datasets")
    exps = load_yaml_dir(cfg_root / "experiments")
    exp = dict(exps[args.exp])
    exp["name"] = exp.get("name", args.exp)
    selected = args.dataset or []
    datasets_by_name = {
        dataset.get("name", "dataset"): dataset
        for dataset in _resolve_datasets(exp, datasets, selected)
    }

    try:
        require_experiment_root(
            exp_root,
            kind="runs",
            identity={"experiment": exp["name"]},
        )
    except ExperimentRootError as exc:
        raise SystemExit(str(exc)) from exc

    updated = 0
    skipped = {"unreadable": 0, "incomplete": 0, "other_experiment": 0, "unmatched_dataset": 0}
    for meta_path in sorted(exp_root.rglob("meta.json")):
        try:
            meta = json.loads(meta_path.read_text())
        except (OSError, json.JSONDecodeError):
            skipped["unreadable"] += 1
            continue
        if not isinstance(meta, dict) or not (
            _run_completed(meta)
            or (meta.get("status") == "evaluation_failed" and meta.get("return_code") == 0)
        ):
            skipped["incomplete"] += 1
            continue
        if meta.get("exp") != exp["name"]:
            skipped["other_experiment"] += 1
            continue
        dataset = datasets_by_name.get(str(meta.get("dataset")))
        if dataset is None:
            skipped["unmatched_dataset"] += 1
            continue
        metrics = build_run_metrics(exp, dataset, meta.get("outputs") or {})
        try:
            write_json_within(exp_root, meta_path.parent / "metrics.json", metrics)
            if meta.get("status") == "evaluation_failed":
                meta["status"] = "succeeded"
                meta.pop("evaluation_error", None)
                write_json_within(exp_root, meta_path, meta)
        except ExperimentRootError as exc:
            raise SystemExit(str(exc)) from exc
        updated += 1

    print(f"recomputed {updated} run(s) under {exp_root}")
    reasons = ", ".join(f"{name}={count}" for name, count in skipped.items() if count)
    if reasons:
        print(f"skipped runs: {reasons}", file=sys.stderr)
    elif not updated:
        print(f"no completed runs found under {exp_root}", file=sys.stderr)


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

    runs_root = Path(args.runs)
    runner = BuildRunner(runs_root)
    tool = tool_cls(tool_config)
    executor = _make_executor()

    if args.dry_run:
        runs_root.mkdir(parents=True, exist_ok=True)
        return

    try:
        claim_experiment_root(
            runs_root,
            kind="builds",
            identity={
                "experiment": build["name"],
                "tool": tool.name,
                "config_sha256": configuration_digest(build, getattr(tool, "config", None)),
            },
        )
        result = runner.run(build=build, tool=tool, executor=executor)
    except (ExperimentRootError, ReusedRunDirectoryError) as exc:
        raise SystemExit(str(exc)) from exc

    meta = (result or {}).get("meta") if isinstance(result, dict) else None
    if isinstance(meta, dict) and meta.get("status") != "succeeded":
        print(
            f"failed build: {build['name']} return_code={meta.get('return_code')}",
            file=sys.stderr,
        )
        raise SystemExit(1)


def main() -> None:
    p = argparse.ArgumentParser("chimera-bench")
    sub = p.add_subparsers(dest="cmd", required=True)

    run_p = sub.add_parser("run")
    run_p.add_argument("--exp", required=True)
    run_p.add_argument("--config", default="configs")
    run_p.add_argument(
        "--runs",
        required=True,
        help="experiment root directory; created and marked for this experiment, never shared",
    )
    run_p.add_argument("--chimera-bin", default="Chimera")
    run_p.add_argument("--ganon-bin", default="ganon")
    run_p.add_argument("--ganon-env", default="ganon")
    run_p.add_argument("--sylph-bin", default="sylph")
    run_p.add_argument("--sylph-env", default="sylph")
    run_p.add_argument("--dry-run", action="store_true")
    run_p.add_argument("--dataset", action="append", default=[])
    run_p.set_defaults(func=run_cmd)

    report_p = sub.add_parser("report")
    report_p.add_argument("--runs", required=True, help="experiment root to export from")
    report_p.add_argument("--out", required=True, help="TSV path inside the experiment root")
    report_p.add_argument("--dataset", action="append", default=[])
    report_p.set_defaults(func=report_cmd)

    recompute_p = sub.add_parser("recompute")
    recompute_p.add_argument("--exp", required=True)
    recompute_p.add_argument("--config", default="configs")
    recompute_p.add_argument(
        "--runs", required=True, help="experiment root holding the completed runs"
    )
    recompute_p.add_argument("--dataset", action="append", default=[])
    recompute_p.set_defaults(func=recompute_cmd)

    build_p = sub.add_parser("build")
    build_p.add_argument("--build", required=True)
    build_p.add_argument("--config", default="configs")
    build_p.add_argument(
        "--runs",
        required=True,
        help="build experiment root directory; created and marked for this build, never shared",
    )
    build_p.add_argument("--ganon-bin", default="ganon")
    build_p.add_argument("--ganon-env", default="ganon")
    build_p.add_argument("--sylph-bin", default="sylph")
    build_p.add_argument("--sylph-env", default="sylph")
    build_p.add_argument("--dry-run", action="store_true")
    build_p.set_defaults(func=build_cmd)

    catalog_p = sub.add_parser("catalog")
    catalog_p.add_argument("--config", default="configs")
    catalog_p.add_argument("--resources-root", default="resources")
    catalog_p.set_defaults(func=catalog_cmd)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
