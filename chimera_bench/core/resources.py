from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable


def _parse_elapsed_to_seconds(value: str) -> float | None:
    value = value.strip()
    if not value:
        return None
    parts = value.split(":")
    try:
        if len(parts) == 3:
            hours = float(parts[0])
            minutes = float(parts[1])
            seconds = float(parts[2])
        elif len(parts) == 2:
            hours = 0.0
            minutes = float(parts[0])
            seconds = float(parts[1])
        else:
            hours = 0.0
            minutes = 0.0
            seconds = float(parts[0])
    except ValueError:
        return None
    return hours * 3600.0 + minutes * 60.0 + seconds


def parse_time_log(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    data: Dict[str, Any] = {}
    for raw in path.read_text().splitlines():
        if ": " not in raw:
            continue
        key, value = raw.split(": ", 1)
        key = key.strip()
        value = value.strip()
        if key == "User time (seconds)":
            try:
                data["user_time_seconds"] = float(value)
            except ValueError:
                continue
        elif key == "System time (seconds)":
            try:
                data["system_time_seconds"] = float(value)
            except ValueError:
                continue
        elif key == "Percent of CPU this job got":
            value = value.replace("%", "")
            try:
                data["cpu_percent"] = float(value)
            except ValueError:
                continue
        elif key == "Elapsed (wall clock) time (h:mm:ss or m:ss)":
            data["elapsed_wall"] = value
            parsed = _parse_elapsed_to_seconds(value)
            if parsed is not None:
                data["elapsed_wall_seconds"] = parsed
        elif key == "Maximum resident set size (kbytes)":
            try:
                data["max_rss_kb"] = int(value)
            except ValueError:
                continue
        elif key == "Exit status":
            try:
                data["exit_status"] = int(value)
            except ValueError:
                continue
    return data


def aggregate_resources(step_records: Iterable[dict]) -> Dict[str, Any]:
    max_rss = None
    user_total = 0.0
    sys_total = 0.0
    saw_user = False
    saw_sys = False
    for step in step_records:
        resource = step.get("resource") or {}
        rss = resource.get("max_rss_kb")
        if isinstance(rss, int):
            max_rss = rss if max_rss is None else max(max_rss, rss)
        user = resource.get("user_time_seconds")
        if isinstance(user, (int, float)):
            user_total += float(user)
            saw_user = True
        sys = resource.get("system_time_seconds")
        if isinstance(sys, (int, float)):
            sys_total += float(sys)
            saw_sys = True
    aggregated: Dict[str, Any] = {}
    if max_rss is not None:
        aggregated["max_rss_kb"] = max_rss
    if saw_user:
        aggregated["user_time_seconds"] = user_total
    if saw_sys:
        aggregated["system_time_seconds"] = sys_total
    return aggregated
