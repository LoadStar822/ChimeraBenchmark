from __future__ import annotations

from pathlib import Path


def ensure_run_dirs(root: Path, exp: str, tool: str, dataset: str) -> Path:
    run_dir = root / exp / tool / dataset
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    (run_dir / "outputs").mkdir(parents=True, exist_ok=True)
    return run_dir


def ensure_build_dirs(root: Path, tool: str, db_name: str) -> Path:
    run_dir = root / tool / db_name
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    (run_dir / "outputs").mkdir(parents=True, exist_ok=True)
    return run_dir


def ensure_profile_dirs(root: Path, exp: str, tool: str, dataset: str) -> Path:
    run_dir = root / exp / tool / dataset
    (run_dir / "outputs").mkdir(parents=True, exist_ok=True)
    return run_dir
