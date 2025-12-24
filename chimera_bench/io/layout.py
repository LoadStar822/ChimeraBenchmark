from __future__ import annotations

from datetime import datetime
import hashlib
from pathlib import Path


def make_run_id(exp: str, tool: str, dataset: str) -> str:
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    key = f"{exp}:{tool}:{dataset}:{ts}"
    short = hashlib.sha1(key.encode()).hexdigest()[:8]
    return f"{ts}-{short}"


def ensure_run_dirs(root: Path, exp: str, tool: str, dataset: str, run_id: str) -> Path:
    run_dir = root / exp / tool / dataset / run_id
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    (run_dir / "outputs").mkdir(parents=True, exist_ok=True)
    return run_dir


def ensure_build_dirs(root: Path, tool: str, db_name: str) -> Path:
    run_dir = root / tool / db_name
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    (run_dir / "outputs").mkdir(parents=True, exist_ok=True)
    return run_dir
