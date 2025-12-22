from __future__ import annotations

from pathlib import Path


def write_summary(runs: list[dict], path: Path) -> None:
    if not runs:
        path.write_text("")
        return
    keys = ["exp", "tool", "dataset", "run_id"]
    metric_keys = sorted({k for r in runs for k in r.get("metrics", {}).keys()})
    header = keys + metric_keys
    lines = ["\t".join(header)]
    for r in runs:
        row = [r.get(k, "") for k in keys]
        metrics = r.get("metrics", {})
        row += [str(metrics.get(k, "")) for k in metric_keys]
        lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n")
