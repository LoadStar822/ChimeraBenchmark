from __future__ import annotations

import hashlib
import json
import os
import tempfile
from datetime import datetime
from pathlib import Path

EXPERIMENT_MARKER = "experiment.json"
MARKER_SCHEMA_VERSION = 1


class ReusedRunDirectoryError(RuntimeError):
    """Raised when a run directory already holds artifacts from a previous run."""


class ExperimentRootError(RuntimeError):
    """Raised when an experiment root is unmarked, foreign or otherwise unusable."""

def configuration_digest(config: dict, tool_config: dict | None = None) -> str:
    payload = json.dumps([config, tool_config or {}], sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _scoped_pair(root: Path, first: str, second: str) -> Path:
    for component in (first, second):
        if not isinstance(component, str) or not component or component in {".", ".."} or any(c in component for c in "/\\:"):
            raise ExperimentRootError(f"Invalid experiment path component: {component!r}")
    path = root / first / second
    if not path.resolve().is_relative_to(root.resolve()):
        raise ExperimentRootError(f"Experiment path escapes its root: {path}")
    return path



def write_json_within(root: Path, path: Path, payload: dict) -> None:
    if not path.resolve().is_relative_to(root.resolve()) or path.is_symlink():
        raise ExperimentRootError(f"JSON output escapes its experiment or is a symlink: {path}")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(payload, handle, indent=2)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def run_dir_path(root: Path, tool: str, dataset: str) -> Path:
    """Run directory for one tool/dataset pair under an experiment root."""
    return _scoped_pair(root, tool, dataset)


def build_dir_path(root: Path, tool: str, db_name: str) -> Path:
    """Build directory for one tool/database pair under a build experiment root."""
    return _scoped_pair(root, tool, db_name)


def ensure_run_dirs(root: Path, tool: str, dataset: str) -> Path:
    run_dir = run_dir_path(root, tool, dataset)
    try:
        run_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise ReusedRunDirectoryError(f"Run directory already exists: {run_dir}") from exc
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    (run_dir / "outputs").mkdir(parents=True, exist_ok=True)
    return run_dir


def ensure_build_dirs(root: Path, tool: str, db_name: str) -> Path:
    run_dir = build_dir_path(root, tool, db_name)
    try:
        run_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise ReusedRunDirectoryError(f"Build directory already exists: {run_dir}") from exc
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    (run_dir / "outputs").mkdir(parents=True, exist_ok=True)
    return run_dir


def require_fresh_run_dir(run_dir: Path) -> None:
    """Reject a run directory that already contains artifacts.

    Runs never overwrite previous results: a new root (``--runs``) must be
    selected instead. Directory creation reserves the run before execution.
    """
    if run_dir.exists() or run_dir.is_symlink():
        raise ReusedRunDirectoryError(
            f"run directory already exists: {run_dir}; select a new experiment root"
        )


def read_experiment_marker(root: Path) -> dict | None:
    marker_path = root / EXPERIMENT_MARKER
    if not marker_path.exists():
        return None
    try:
        data = json.loads(marker_path.read_text())
    except json.JSONDecodeError as exc:
        raise ExperimentRootError(f"unreadable experiment marker {marker_path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ExperimentRootError(f"invalid experiment marker {marker_path}: expected a JSON object")
    return data


def _require_matching_marker(root: Path, marker: dict, *, kind: str, identity: dict) -> None:
    if marker.get("schema_version") != MARKER_SCHEMA_VERSION:
        raise ExperimentRootError(f"Unsupported experiment marker schema: {root}")
    marker_kind = marker.get("kind")
    if marker_kind != kind:
        raise ExperimentRootError(
            f"experiment root {root} holds {marker_kind!r} results, not {kind!r}; "
            "select a root created by the matching command"
        )
    mismatched = [
        f"{key}: root has {marker.get(key)!r}, expected {value!r}"
        for key, value in identity.items()
        if marker.get(key) != value
    ]
    if mismatched:
        raise ExperimentRootError(
            f"experiment root {root} belongs to a different experiment ({'; '.join(mismatched)})"
        )


def claim_experiment_root(root: Path, *, kind: str, identity: dict) -> dict:
    """Validate an experiment root for execution, creating its marker when fresh.

    A root that already contains files without a marker is never adopted: that
    is how aggregate or formal result trees (e.g. a shared ``results`` tree)
    stay read-only for execution.
    """
    existing = read_experiment_marker(root)
    if existing is not None:
        _require_matching_marker(root, existing, kind=kind, identity=identity)
        return existing
    if root.exists() and any(root.iterdir()):
        raise ExperimentRootError(
            f"{root} already contains files but no {EXPERIMENT_MARKER} marker; refusing to "
            "write into an unmarked results tree, select a new experiment root"
        )
    marker = {
        "schema_version": MARKER_SCHEMA_VERSION,
        "kind": kind,
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        **identity,
    }
    root.mkdir(parents=True, exist_ok=True)
    try:
        with (root / EXPERIMENT_MARKER).open("x") as handle:
            handle.write(json.dumps(marker, indent=2))
    except FileExistsError:
        existing = read_experiment_marker(root)
        if existing is None:
            raise ExperimentRootError(f"Concurrent experiment root mutation: {root}")
        _require_matching_marker(root, existing, kind=kind, identity=identity)
        return existing
    return marker


def require_experiment_root(root: Path, *, kind: str, identity: dict | None = None) -> dict:
    """Validate that a read-only consumer targets a marked experiment root."""
    marker = read_experiment_marker(root)
    if marker is None:
        raise ExperimentRootError(
            f"{root} is not an experiment root marked with {EXPERIMENT_MARKER}; "
            "execution roots are created by run/build"
        )
    _require_matching_marker(root, marker, kind=kind, identity=identity or {})
    return marker
