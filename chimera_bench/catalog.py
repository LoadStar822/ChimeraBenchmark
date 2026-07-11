from __future__ import annotations

import csv
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from .config import load_yaml_dir
from .dataset_prepare import resolve_strain_madness_reads


PUBLIC_DATASET_ORDER = [
    "cami2-marine-long-sample0",
    "cami2-marine-long",
    "cami2-marine-short",
    "cami-strain-madness-long-sample0",
    "cami-strain-madness-short-sample0",
    "cami-strain-madness-long",
    "cami-strain-madness-short",
    "atcc-illumina",
    "atcc-hifi",
    "zymo-gridion-even",
    "zymo-gridion-log",
    "zymo-promethion-even",
    "zymo-promethion-log",
    "prjna637878-supported19",
    "prjna637878-supported19-single-read",
]
PUBLIC_BUILD_DB_ORDER = ["cami_refseq", "refseq_complete"]
DATASET_DISPLAY_NAMES = {
    "cami2-marine-long-sample0": "CAMI II Marine (long-read contigs, sample 0)",
    "cami2-marine-long": "CAMI II Marine (long-read contigs)",
    "cami2-marine-short": "CAMI II Marine (short-read contigs)",
    "cami-strain-madness-long-sample0": "CAMI II Strain Madness (long reads, sample 0)",
    "cami-strain-madness-short-sample0": "CAMI II Strain Madness (short reads, sample 0)",
    "cami-strain-madness-long": "CAMI II Strain Madness (long reads)",
    "cami-strain-madness-short": "CAMI II Strain Madness (short reads)",
    "atcc-illumina": "ATCC MSA-1003 (Illumina)",
    "atcc-hifi": "ATCC MSA-1003 (PacBio HiFi)",
    "zymo-gridion-even": "ZymoBIOMICS Even (GridION)",
    "zymo-gridion-log": "ZymoBIOMICS Log (GridION)",
    "zymo-promethion-even": "ZymoBIOMICS Even (PromethION)",
    "zymo-promethion-log": "ZymoBIOMICS Log (PromethION)",
    "prjna637878-supported19": "PRJNA637878 Supported Fecal Genomes",
    "prjna637878-supported19-single-read": "PRJNA637878 Supported Fecal Genomes (single-read)",
}
DB_DISPLAY_NAMES = {
    "cami_refseq": "CAMI RefSeq",
    "refseq_complete": "RefSeq Complete",
}
DB_SOURCE_NAMES = {
    "cami_refseq": "NCBI RefSeq genomic, CAMI snapshot",
    "refseq_complete": "NCBI RefSeq complete genomes",
}
MISSING = "—"
CATALOG_CACHE_VERSION = 2
_SEQKIT_CMD: list[str] | None | bool = False


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="strict")
    return path.open("r", encoding="utf-8", errors="strict")


def _file_signature(path: Path) -> str:
    stat = path.stat()
    return f"{stat.st_size}:{stat.st_mtime_ns}"


def _read_cache(cache_path: Path) -> dict[str, Any]:
    if not cache_path.exists():
        return {"version": CATALOG_CACHE_VERSION, "files": {}}
    try:
        payload = json.loads(cache_path.read_text())
    except json.JSONDecodeError:
        return {"version": CATALOG_CACHE_VERSION, "files": {}}
    if not isinstance(payload, dict):
        return {"version": CATALOG_CACHE_VERSION, "files": {}}
    if payload.get("version") != CATALOG_CACHE_VERSION:
        return {"version": CATALOG_CACHE_VERSION, "files": {}}
    if not isinstance(payload.get("files"), dict):
        payload["files"] = {}
    return payload


def _write_cache(cache_path: Path, payload: dict[str, Any]) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    payload["version"] = CATALOG_CACHE_VERSION
    cache_path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def _detect_sequence_format(path: Path) -> str:
    lower = path.name.lower()
    fasta_suffixes = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
    fastq_suffixes = (".fq", ".fastq", ".fq.gz", ".fastq.gz")
    if lower.endswith(fasta_suffixes):
        return "fasta"
    if lower.endswith(fastq_suffixes):
        return "fastq"
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                return "fasta"
            if line.startswith("@"):
                return "fastq"
            break
    raise ValueError(f"unsupported sequence file format: {path}")


def _parse_int(value: Any) -> int:
    return int(str(value).strip().replace(",", ""))


def _parse_float(value: Any) -> float | None:
    text = str(value).strip().replace(",", "")
    if not text or text in {"-", "NA", "na"}:
        return None
    return float(text)


def _format_decimal(value: float | None, digits: int = 1) -> str:
    if value is None:
        return MISSING
    return f"{value:.{digits}f}"


def _format_gb(size_bytes: int) -> str:
    return f"{size_bytes / 1_000_000_000:.1f}"


def _format_int(value: int | None) -> str:
    if value is None:
        return MISSING
    return f"{value:,}"


def _catalog_number_or_text(value: Any) -> int | str:
    if isinstance(value, int):
        return value
    text = str(value).strip()
    if text in {"", MISSING}:
        return MISSING
    return _parse_int(text)


def _n50_from_hist(length_counts: dict[int, int] | Counter[int], total_bases: int) -> int:
    if total_bases <= 0:
        return 0
    threshold = total_bases / 2
    running = 0
    for length in sorted(length_counts, reverse=True):
        running += int(length) * int(length_counts[length])
        if running >= threshold:
            return int(length)
    return 0


def _stats_complete(stats: dict[str, Any], *, require_histogram: bool) -> bool:
    required = {
        "format",
        "records",
        "total_bases",
        "min_len",
        "mean_len",
        "max_len",
        "n50",
        "gc_percent",
        "size_bytes",
    }
    if not required.issubset(stats):
        return False
    if stats.get("format") == "fastq" and "q30_percent" not in stats:
        return False
    if require_histogram and "length_counts" not in stats:
        return False
    return True


def _finalize_sequence_stats(
    *,
    fmt: str,
    records: int,
    total_bases: int,
    min_len: int,
    max_len: int,
    length_counts: Counter[int],
    gc_bases: int,
    q30_bases: int | None,
    size_bytes: int,
    keep_histogram: bool,
) -> dict[str, Any]:
    mean_len = (total_bases / records) if records else 0.0
    stats: dict[str, Any] = {
        "format": fmt,
        "records": records,
        "total_bases": total_bases,
        "min_len": min_len if records else 0,
        "mean_len": mean_len,
        "max_len": max_len if records else 0,
        "n50": _n50_from_hist(length_counts, total_bases),
        "gc_percent": (gc_bases * 100 / total_bases) if total_bases else 0.0,
        "size_bytes": size_bytes,
    }
    if fmt == "fastq":
        stats["q30_percent"] = (q30_bases or 0) * 100 / total_bases if total_bases else 0.0
    else:
        stats["q30_percent"] = None
    if keep_histogram:
        stats["length_counts"] = {str(length): count for length, count in length_counts.items()}
        stats["gc_bases"] = gc_bases
        if q30_bases is not None:
            stats["q30_bases"] = q30_bases
    return stats


def _scan_fastq(path: Path, *, keep_histogram: bool) -> dict[str, Any]:
    records = 0
    total_bases = 0
    min_len: int | None = None
    max_len = 0
    gc_bases = 0
    q30_bases = 0
    length_counts: Counter[int] = Counter()
    with _open_text(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n\r")
            plus = fh.readline()
            qual = fh.readline().rstrip("\n\r")
            if not seq or not plus or not qual:
                raise ValueError(f"FASTQ record truncated: {path}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"invalid FASTQ record: {path}")
            if len(seq) != len(qual):
                raise ValueError(f"FASTQ sequence/quality length mismatch: {path}")
            length = len(seq)
            records += 1
            total_bases += length
            min_len = length if min_len is None else min(min_len, length)
            max_len = max(max_len, length)
            length_counts[length] += 1
            for base in seq:
                if base in "GgCc":
                    gc_bases += 1
            for char in qual:
                if ord(char) - 33 >= 30:
                    q30_bases += 1
    return _finalize_sequence_stats(
        fmt="fastq",
        records=records,
        total_bases=total_bases,
        min_len=min_len or 0,
        max_len=max_len,
        length_counts=length_counts,
        gc_bases=gc_bases,
        q30_bases=q30_bases,
        size_bytes=path.stat().st_size,
        keep_histogram=keep_histogram,
    )


def _scan_fasta(path: Path, *, keep_histogram: bool) -> dict[str, Any]:
    records = 0
    total_bases = 0
    min_len: int | None = None
    max_len = 0
    current_bases = 0
    current_gc = 0
    gc_bases = 0
    length_counts: Counter[int] = Counter()
    saw_header = False

    def finish_record() -> None:
        nonlocal records, total_bases, min_len, max_len, current_bases, current_gc, gc_bases
        records += 1
        total_bases += current_bases
        gc_bases += current_gc
        min_len = current_bases if min_len is None else min(min_len, current_bases)
        max_len = max(max_len, current_bases)
        length_counts[current_bases] += 1

    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if saw_header:
                    finish_record()
                saw_header = True
                current_bases = 0
                current_gc = 0
                continue
            if not saw_header:
                raise ValueError(f"FASTA sequence before header: {path}")
            current_bases += len(line)
            for base in line:
                if base in "GgCc":
                    current_gc += 1
    if saw_header:
        finish_record()
    return _finalize_sequence_stats(
        fmt="fasta",
        records=records,
        total_bases=total_bases,
        min_len=min_len or 0,
        max_len=max_len,
        length_counts=length_counts,
        gc_bases=gc_bases,
        q30_bases=None,
        size_bytes=path.stat().st_size,
        keep_histogram=keep_histogram,
    )


def _seqkit_cmd() -> list[str] | None:
    global _SEQKIT_CMD
    if _SEQKIT_CMD is not False:
        return list(_SEQKIT_CMD) if _SEQKIT_CMD else None

    seqkit_bin = shutil.which("seqkit")
    if seqkit_bin:
        _SEQKIT_CMD = [seqkit_bin]
        return list(_SEQKIT_CMD)

    conda_bin = shutil.which("conda")
    if not conda_bin:
        _SEQKIT_CMD = None
        return None
    conda_root = Path(conda_bin).resolve().parent.parent
    env_seqkit = conda_root / "envs" / "seqkit" / "bin" / "seqkit"
    if env_seqkit.exists():
        _SEQKIT_CMD = [str(env_seqkit)]
        return list(_SEQKIT_CMD)

    probe = subprocess.run(
        [conda_bin, "run", "-n", "seqkit", "seqkit", "version"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if probe.returncode != 0:
        _SEQKIT_CMD = None
        return None

    _SEQKIT_CMD = [conda_bin, "run", "-n", "seqkit", "seqkit"]
    return list(_SEQKIT_CMD)


def _cache_hit(cache: dict[str, Any], path: Path, *, require_histogram: bool) -> tuple[str, dict[str, Any] | None]:
    key = str(path)
    signature = _file_signature(path)
    cached = cache.setdefault("files", {}).get(key)
    if isinstance(cached, dict) and cached.get("signature") == signature:
        stats = cached.get("stats")
        if isinstance(stats, dict) and _stats_complete(stats, require_histogram=require_histogram):
            return signature, dict(stats)
    return signature, None


def _store_cache(cache: dict[str, Any], path: Path, *, signature: str, stats: dict[str, Any]) -> None:
    cache.setdefault("files", {})[str(path)] = {
        "signature": signature,
        "stats": stats,
        "updated_at": _now_iso(),
    }


def _scan_with_seqkit(paths: list[Path]) -> dict[str, dict[str, Any]]:
    cmd = _seqkit_cmd()
    if not cmd or not paths:
        return {}
    proc = subprocess.run(
        cmd + ["stats", "-Ta", *[str(path) for path in paths]],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return {}

    stats_by_path: dict[str, dict[str, Any]] = {}
    reader = csv.DictReader(proc.stdout.splitlines(), delimiter="\t")
    for row in reader:
        file_path = row.get("file")
        if not file_path:
            continue
        path = Path(file_path)
        fmt = str(row.get("format", "")).lower()
        if fmt not in {"fastq", "fasta"}:
            continue
        total_bases = _parse_int(row["sum_len"])
        records = _parse_int(row["num_seqs"])
        stats_by_path[str(path)] = {
            "format": fmt,
            "records": records,
            "total_bases": total_bases,
            "min_len": _parse_int(row["min_len"]),
            "mean_len": float(_parse_float(row["avg_len"]) or 0.0),
            "max_len": _parse_int(row["max_len"]),
            "n50": _parse_int(row["N50"]),
            "gc_percent": float(_parse_float(row.get("GC(%)")) or 0.0),
            "q30_percent": float(_parse_float(row.get("Q30(%)")) or 0.0) if fmt == "fastq" else None,
            "size_bytes": path.stat().st_size,
        }
    return stats_by_path


def _stats_from_seqkit_row(row: dict[str, str], *, path: Path | None, size_bytes: int) -> dict[str, Any] | None:
    fmt = str(row.get("format", "")).lower()
    if fmt not in {"fastq", "fasta"}:
        return None
    stats = {
        "format": fmt,
        "records": _parse_int(row["num_seqs"]),
        "total_bases": _parse_int(row["sum_len"]),
        "min_len": _parse_int(row["min_len"]),
        "mean_len": float(_parse_float(row["avg_len"]) or 0.0),
        "max_len": _parse_int(row["max_len"]),
        "n50": _parse_int(row["N50"]),
        "gc_percent": float(_parse_float(row.get("GC(%)")) or 0.0),
        "q30_percent": float(_parse_float(row.get("Q30(%)")) or 0.0) if fmt == "fastq" else None,
        "size_bytes": size_bytes,
    }
    if path is not None:
        stats["path"] = str(path)
    return stats


def _scan_group_with_seqkit(paths: list[Path]) -> dict[str, Any] | None:
    cmd = _seqkit_cmd()
    if not cmd or not paths:
        return None
    if any(path.suffix == ".gz" for path in paths):
        return None
    formats = {_detect_sequence_format(path) for path in paths}
    if len(formats) != 1:
        return None
    cat_cmd = "cat " + " ".join(shlex.quote(str(path)) for path in paths)
    seqkit_cmd = " ".join(shlex.quote(part) for part in [*cmd, "stats", "-Ta", "-"])
    proc = subprocess.run(
        ["bash", "-lc", f"{cat_cmd} | {seqkit_cmd}"],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return None
    reader = csv.DictReader(proc.stdout.splitlines(), delimiter="\t")
    rows = list(reader)
    if len(rows) != 1:
        return None
    return _stats_from_seqkit_row(rows[0], path=None, size_bytes=sum(path.stat().st_size for path in paths))


def scan_sequence_files(paths: list[Path], *, cache: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    cache = cache if cache is not None else {"version": CATALOG_CACHE_VERSION, "files": {}}
    require_histogram = len(paths) > 1
    signatures: dict[str, str] = {}
    results: dict[str, dict[str, Any]] = {}
    pending: list[Path] = []

    for path in paths:
        signature, cached = _cache_hit(cache, path, require_histogram=require_histogram)
        signatures[str(path)] = signature
        if cached is not None:
            results[str(path)] = cached
        else:
            pending.append(path)

    seqkit_pending = pending if not require_histogram else []
    seqkit_results = _scan_with_seqkit(seqkit_pending)
    for path in pending:
        key = str(path)
        fmt = _detect_sequence_format(path)
        stats = seqkit_results.get(key)
        if stats is None:
            if fmt == "fastq":
                stats = _scan_fastq(path, keep_histogram=require_histogram)
            else:
                stats = _scan_fasta(path, keep_histogram=require_histogram)
        _store_cache(cache, path, signature=signatures[key], stats=stats)
        results[key] = dict(stats)

    return [dict(results[str(path)]) for path in paths]


def scan_sequence_group(paths: list[Path], *, cache: dict[str, Any] | None = None) -> dict[str, Any]:
    cache = cache if cache is not None else {"version": CATALOG_CACHE_VERSION, "files": {}}
    if len(paths) == 1:
        return scan_sequence_files(paths, cache=cache)[0]

    signatures = [_file_signature(path) for path in paths]
    group_key = "group:" + "|".join(f"{path}:{signature}" for path, signature in zip(paths, signatures))
    cached = cache.setdefault("files", {}).get(group_key)
    if isinstance(cached, dict):
        stats = cached.get("stats")
        if isinstance(stats, dict) and _stats_complete(stats, require_histogram=False):
            return dict(stats)

    stats = _scan_group_with_seqkit(paths)
    if stats is None:
        stats_list = scan_sequence_files(paths, cache=cache)
        stats = _merge_sequence_stats(stats_list)
    cache.setdefault("files", {})[group_key] = {
        "signature": "|".join(signatures),
        "stats": stats,
        "updated_at": _now_iso(),
    }
    return dict(stats)


def _db_name_from_path(value: str | Path | None) -> str:
    if not value:
        return ""
    name = Path(value).name
    for suffix in (".imcf", ".syldb", ".hixf"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _is_public_dataset(name: str, data: dict[str, Any]) -> bool:
    if name == "example":
        return False
    catalog = data.get("catalog")
    if isinstance(catalog, dict) and catalog.get("exclude"):
        return False
    return bool(data.get("group") or data.get("truth_dir") or data.get("truth_profile"))


def _is_public_build(name: str, data: dict[str, Any]) -> bool:
    if name == "ganon-local":
        return False
    return _db_name_from_path(data.get("db_prefix") or data.get("db")) != "local"


def _ordered_dataset_names(datasets: dict[str, dict[str, Any]]) -> list[str]:
    names = {name for name, data in datasets.items() if _is_public_dataset(name, data)}
    ordered = [name for name in PUBLIC_DATASET_ORDER if name in names]
    ordered.extend(sorted(names - set(ordered)))
    return ordered


def _ordered_build_db_names(builds: dict[str, dict[str, Any]]) -> list[str]:
    names = {
        _db_name_from_path(data.get("db_prefix") or data.get("db"))
        for name, data in builds.items()
        if _is_public_build(name, data)
    }
    names.discard("")
    ordered = [name for name in PUBLIC_BUILD_DB_ORDER if name in names]
    ordered.extend(sorted(names - set(ordered)))
    return ordered


def _dataset_input_paths(dataset: dict[str, Any]) -> list[Path]:
    samples = dataset.get("samples")
    if samples:
        paths: list[Path] = []
        for sample in samples:
            if sample.get("paired"):
                paths.extend(Path(path) for path in sample["paired"])
            elif sample.get("reads"):
                paths.extend(Path(path) for path in sample["reads"])
            else:
                raise ValueError(f"dataset sample has no inputs: {dataset.get('name')}:{sample.get('sample_id')}")
        return paths
    if dataset.get("paired"):
        return [Path(path) for path in dataset["paired"]]
    if dataset.get("reads"):
        return [Path(path) for path in dataset["reads"]]
    source = dataset.get("source") or {}
    if source.get("kind") == "strain_madness" and source.get("read_type") == "long":
        sample_ids = dataset.get("sample_ids") or source.get("sample_ids") or []
        return resolve_strain_madness_reads(Path(source["root"]), sample_ids)
    return []


def _dataset_sample_count(dataset: dict[str, Any]) -> int:
    if dataset.get("samples"):
        return len(dataset["samples"])
    if dataset.get("paired"):
        return 1
    source = dataset.get("source") or {}
    if source.get("kind") == "strain_madness":
        sample_ids = dataset.get("sample_ids") or source.get("sample_ids") or []
        return len(sample_ids)
    return len(dataset.get("reads") or [])


def _dataset_truth(dataset: dict[str, Any]) -> str:
    if dataset.get("truth_label"):
        return str(dataset["truth_label"])
    if dataset.get("truth_dir"):
        return "per-read mapping + profile"
    if dataset.get("truth_profile"):
        return "profile only"
    if dataset.get("truth_map"):
        return "per-read mapping only"
    return "unknown"


def _dataset_input_type(dataset: dict[str, Any], paths: list[Path]) -> str:
    if dataset.get("input_type"):
        return str(dataset["input_type"])
    samples = dataset.get("samples")
    if samples:
        if all(sample.get("paired") for sample in samples):
            return "paired FASTQ"
        if all(sample.get("reads") for sample in samples):
            formats = {_detect_sequence_format(Path(path)) for sample in samples for path in sample["reads"]}
            if formats == {"fastq"}:
                return "single FASTQ"
            if formats == {"fasta"}:
                return "contig FASTA"
        return "mixed FASTA/FASTQ"
    if dataset.get("paired"):
        return "paired FASTQ"
    formats = {_detect_sequence_format(path) for path in paths}
    if formats == {"fasta"}:
        return "contig FASTA"
    if formats == {"fastq"}:
        return "single FASTQ"
    return "mixed FASTA/FASTQ"


def _dataset_catalog_row_from_metadata(
    name: str,
    dataset: dict[str, Any],
    paths: list[Path],
) -> dict[str, Any] | None:
    catalog = dataset.get("catalog")
    if not catalog:
        return None
    if not isinstance(catalog, dict):
        raise ValueError(f"dataset catalog metadata must be an object: {name}")
    required = [
        "total_size_gb",
        "samples",
        "input_type",
        "reads_or_contigs",
        "base_pairs_bp",
        "mean_length_bp",
        "truth",
    ]
    missing = [key for key in required if key not in catalog]
    if missing:
        raise ValueError(f"dataset catalog metadata missing fields for {name}: {', '.join(missing)}")
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)
    return {
        "dataset": name,
        "dataset_name": DATASET_DISPLAY_NAMES.get(name, name),
        "total_size_gb": str(catalog["total_size_gb"]),
        "samples": int(catalog["samples"]),
        "input_type": str(catalog["input_type"]),
        "reads_or_contigs": _catalog_number_or_text(catalog["reads_or_contigs"]),
        "base_pairs_bp": _catalog_number_or_text(catalog["base_pairs_bp"]),
        "mean_length_bp": _catalog_number_or_text(catalog["mean_length_bp"]),
        "n50_bp": catalog.get("n50_bp", MISSING),
        "gc_percent": catalog.get("gc_percent", MISSING),
        "q30_percent": catalog.get("q30_percent", MISSING),
        "truth": str(catalog["truth"]),
    }


def _hist_from_stats(stats: dict[str, Any]) -> Counter[int]:
    raw = stats.get("length_counts")
    if not isinstance(raw, dict):
        raise ValueError("missing length histogram for multi-file dataset")
    return Counter({int(length): int(count) for length, count in raw.items()})


def _merge_sequence_stats(stats_list: list[dict[str, Any]]) -> dict[str, Any]:
    if not stats_list:
        raise ValueError("empty stats list")
    if len(stats_list) == 1:
        return dict(stats_list[0])

    records = sum(int(stats["records"]) for stats in stats_list)
    total_bases = sum(int(stats["total_bases"]) for stats in stats_list)
    total_bytes = sum(int(stats["size_bytes"]) for stats in stats_list)
    gc_bases = sum(int(stats.get("gc_bases", 0)) for stats in stats_list)
    has_fastq = any(stats.get("format") == "fastq" for stats in stats_list)
    q30_bases = sum(int(stats.get("q30_bases", 0)) for stats in stats_list) if has_fastq else None
    length_counts: Counter[int] = Counter()
    for stats in stats_list:
        length_counts.update(_hist_from_stats(stats))
    return {
        "records": records,
        "total_bases": total_bases,
        "size_bytes": total_bytes,
        "min_len": min(int(stats["min_len"]) for stats in stats_list),
        "mean_len": (total_bases / records) if records else 0.0,
        "max_len": max(int(stats["max_len"]) for stats in stats_list),
        "n50": _n50_from_hist(length_counts, total_bases),
        "gc_percent": (gc_bases * 100 / total_bases) if total_bases else 0.0,
        "q30_percent": (q30_bases * 100 / total_bases) if q30_bases is not None and total_bases else None,
    }


def collect_dataset_rows(
    *,
    config_root: Path,
    cache_path: Path,
    dataset_names: Iterable[str] | None = None,
    progress: bool = False,
) -> list[dict[str, Any]]:
    datasets = load_yaml_dir(config_root / "datasets")
    names = list(dataset_names) if dataset_names is not None else _ordered_dataset_names(datasets)
    cache = _read_cache(cache_path)
    rows: list[dict[str, Any]] = []

    for name in names:
        if name not in datasets:
            continue
        dataset = dict(datasets[name])
        if not _is_public_dataset(name, dataset):
            continue
        paths = _dataset_input_paths(dataset)
        if not paths:
            raise ValueError(f"dataset has no inputs: {name}")
        if progress:
            print(f"[catalog] dataset={name} start files={len(paths)}")
        metadata_row = _dataset_catalog_row_from_metadata(name, dataset, paths)
        if metadata_row is not None:
            rows.append(metadata_row)
            if progress:
                print(f"[catalog] dataset={name} done source=metadata")
            continue
        for path in paths:
            if not path.exists():
                raise FileNotFoundError(path)
        merged = scan_sequence_group(paths, cache=cache)
        row = {
            "dataset": name,
            "dataset_name": DATASET_DISPLAY_NAMES.get(name, name),
            "total_size_gb": _format_gb(int(merged["size_bytes"])),
            "samples": _dataset_sample_count(dataset),
            "input_type": _dataset_input_type(dataset, paths),
            "reads_or_contigs": int(merged["records"]),
            "base_pairs_bp": int(merged["total_bases"]),
            "mean_length_bp": int(round(float(merged["mean_len"]))),
            "n50_bp": int(merged["n50"]),
            "gc_percent": _format_decimal(float(merged["gc_percent"]), 2),
            "q30_percent": _format_decimal(merged.get("q30_percent"), 2),
            "truth": _dataset_truth(dataset),
        }
        rows.append(row)
        _write_cache(cache_path, cache)
        if progress:
            print(
                "[catalog] dataset="
                + name
                + f" done bytes={merged['size_bytes']} records={merged['records']} bases={merged['total_bases']}"
            )

    return rows


def _select_build_source(builds: list[dict[str, Any]]) -> tuple[Path, str]:
    for build in builds:
        build_cfg = build.get("build") or {}
        value = build_cfg.get("target_tsv") or build_cfg.get("input_tsv") or build_cfg.get("input_file")
        if value:
            return Path(value), str(build_cfg.get("taxonomy_dir") or "")
    raise ValueError("reference DB build has no target TSV")


def _target_entries(target_path: Path) -> list[tuple[Path, str]]:
    entries: list[tuple[Path, str]] = []
    with target_path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"invalid target row: {target_path}: {line}")
            path = Path(parts[0])
            if not path.is_absolute():
                path = (target_path.parent / path).resolve()
            entries.append((path, parts[1]))
    return entries


def _accession_from_path(path: Path) -> str:
    match = re.match(r"^(GC[AF]_\d+\.\d+)", path.name)
    if not match:
        raise ValueError(f"cannot parse assembly accession from path: {path}")
    return match.group(1)


def _assembly_summary_candidates(target_path: Path, entries: list[tuple[Path, str]]) -> list[Path]:
    roots: list[Path] = []
    for base in [target_path.parent, *(entries[0][0].parents if entries else [])]:
        roots.append(base)
    candidates: list[Path] = []
    seen: set[Path] = set()
    for root in roots:
        for path in sorted(root.glob("assembly_summary*.txt")):
            if path not in seen:
                seen.add(path)
                candidates.append(path)
    return candidates


def _load_assembly_rows(summary_path: Path, accessions: set[str]) -> dict[str, dict[str, int | str]]:
    rows: dict[str, dict[str, int | str]] = {}
    with summary_path.open("r", encoding="utf-8", errors="strict") as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 31:
                continue
            accession = parts[0]
            if accession not in accessions:
                continue
            rows[accession] = {
                "species_taxid": parts[6],
                "genome_size": _parse_int(parts[25]),
                "contig_count": _parse_int(parts[30]),
            }
            if len(rows) == len(accessions):
                break
    return rows


def _find_assembly_rows(target_path: Path, entries: list[tuple[Path, str]]) -> dict[str, dict[str, int | str]]:
    accessions = {_accession_from_path(path) for path, _taxid in entries}
    merged_rows: dict[str, dict[str, int | str]] = {}
    used_paths: list[Path] = []
    for candidate in _assembly_summary_candidates(target_path, entries):
        rows = _load_assembly_rows(candidate, accessions - set(merged_rows))
        if rows:
            merged_rows.update(rows)
            used_paths.append(candidate)
        if len(merged_rows) == len(accessions):
            return merged_rows
    missing = sorted(accessions - set(merged_rows))[:5]
    raise ValueError(f"assembly summary incomplete for {target_path}; used={used_paths}; missing={missing}")


def _reference_db_row(db_name: str, target_path: Path) -> dict[str, Any]:
    if not target_path.exists():
        raise FileNotFoundError(target_path)
    entries = _target_entries(target_path)
    assembly_rows = _find_assembly_rows(target_path, entries)
    total_size = 0
    total_sequences = 0
    base_pairs = 0
    species: set[str] = set()
    for path, _taxid in entries:
        if not path.exists():
            raise FileNotFoundError(path)
        accession = _accession_from_path(path)
        meta = assembly_rows[accession]
        total_size += path.stat().st_size
        total_sequences += int(meta["contig_count"])
        base_pairs += int(meta["genome_size"])
        species.add(str(meta["species_taxid"]))
    return {
        "db_name": db_name,
        "dataset_name": DB_DISPLAY_NAMES.get(db_name, db_name),
        "total_size_gb": _format_gb(total_size),
        "total_sequences": total_sequences,
        "base_pairs_bp": base_pairs,
        "assemblies": len(entries),
        "species_count": len(species),
        "source": DB_SOURCE_NAMES.get(db_name, "NCBI assembly summary"),
    }


def collect_build_rows(
    *,
    config_root: Path,
    cache_path: Path,
    build_db_names: Iterable[str] | None = None,
) -> list[dict[str, Any]]:
    del cache_path
    builds = load_yaml_dir(config_root / "build")
    names = list(build_db_names) if build_db_names is not None else _ordered_build_db_names(builds)
    grouped: dict[str, list[dict[str, Any]]] = {}
    for build_name, build in builds.items():
        if not _is_public_build(build_name, build):
            continue
        db_name = _db_name_from_path(build.get("db_prefix") or build.get("db"))
        grouped.setdefault(db_name, []).append(dict(build))

    rows: list[dict[str, Any]] = []
    for db_name in names:
        group = grouped.get(db_name)
        if not group:
            continue
        target_path, _taxonomy_source = _select_build_source(group)
        rows.append(_reference_db_row(db_name, target_path))
    return rows


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def _markdown_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]]) -> list[str]:
    header = [label for label, _ in columns]
    lines = [
        "| " + " | ".join(header) + " |",
        "| " + " | ".join(["---"] * len(header)) + " |",
    ]
    for row in rows:
        cells: list[str] = []
        for _label, key in columns:
            value = row.get(key, "")
            if isinstance(value, int):
                cells.append(_format_int(value))
            elif isinstance(value, str) and value.isdigit():
                cells.append(_format_int(int(value)))
            else:
                cells.append(str(value))
        lines.append("| " + " | ".join(cells) + " |")
    return lines


def _results_readme_tail() -> list[str]:
    return [
        "",
        "## 评估任务与口径",
        "",
        "本 benchmark 报告两个分类层级：`species` 和 `genus`。",
        "逐读段分类（per-read classification）评估每条 read 或 contig 的分类结果。",
        "丰度画像（abundance profiling）评估整个样本的物种组成估计。",
        "",
        "在 `species` 或 `genus` 层级计算逐读段结果时，先将真值分类编号（taxid）和预测 taxid 都提升到同一目标层级，再判断是否一致。",
        "这样可以公平处理原始真值为 strain、subspecies 或 isolate 的情况。",
        "",
        "## 数据集特定口径",
        "",
        "### PRJNA637878",
        "",
        "PRJNA637878 的部分测序运行（SRA run）包含未配对读段（singleton reads），表现为 `spots_with_mates` 小于 `spots`。",
        "`prjna637878-supported19-single-read` 将 R1 和 R2 作为独立单读段评估（single-read evaluation），并在生成 FASTQ 和真值表时为读段标识符（read identifier）增加 `__mate1` 或 `__mate2` 配对端标记（mate tag）。",
        "该设置不改变原始实验的 paired-end 来源；它只规定评估时所有工具使用同一批单读段输入，且不使用 mate-pair 信息。",
        "",
        "## 核心结果汇总",
        "",
        "跨数据集汇总由输入和评估口径一致、状态为 `complete` 的运行生成；对应运行信息记录在 `results/paper_run_manifest.tsv`。",
        "逐读段分类汇总覆盖 CAMI II Strain Madness long/short、CAMI II Marine long/short 和 `prjna637878-supported19-single-read`，比较 Chimera、Centrifuger 与 Kraken2。",
        "丰度画像结果独立汇总工具原生输出的 sample-level abundance estimates。",
        "",
        "机器可读结果位于 `results/builds/summary.tsv`、`results/classify/summary.tsv`、`results/classify/sample_metrics.tsv` 和 `results/profile/summary.tsv`。",
        "",
        "## 真实队列 clade-level 实验",
        "",
        "Fna C2 真实队列实验位于 `results/real/`。",
        "该实验在默认 `species/genus` benchmark 之外，评估原论文定义的近邻 clade 信号及其 read-level evidence。",
        "",
        "`results/real/README.md` 汇总三队列 Fna C2 分析的实验设计、指标、结果和配套机器表。",
        "",
        "Detection floor 先根据原论文报告的 C2 abundance 和当前输入深度换算 expected C2 reads，再将达到指定下限的 C2-positive 样本纳入阳性组。",
        "该口径用于区分“总体是否有 C2 信号”和“当 C2 信号达到可检测强度时，工具是否能把样本拉开”。",
        "",
        "## 指标说明",
        "",
        "### Per-read 指标",
        "",
        "逐读段指标衡量单条序列是否被分到正确的目标分类层级。",
        "",
        "- `Precision`：被工具分到某个目标层级的序列中，有多少是正确的。",
        "- `Recall`：真值中可评估的序列里，有多少被正确找回。",
        "- `F1`：precision 和 recall 的调和平均。",
        "- `Truth Mapped Rate`：真值中有多少序列能映射到目标层级。",
        "- `Pred Mapped Rate`：预测中有多少序列能映射到目标层级。",
        "",
        "当 mapped rate 明显低于 `1` 时，F1 代表可映射子集上的表现，不应单独解释为全数据集表现。",
        "",
        "### Profile 指标",
        "",
        "Profile 表使用 OPAL 核心指标（OPAL core metrics）：",
        "- `Completeness`：真值中出现的 taxa 有多少被预测到。",
        "- `Purity`：预测到的 taxa 中有多少确实存在于真值。",
        "- `L1 Norm`：真值和预测丰度分布之间的 L1 距离，范围为 `0..2`。",
        "- `Weighted UniFrac`：沿分类体系树（taxonomy tree）比较两边的累积丰度差异。",
        "",
        "Completeness 和 purity 越高越好；L1 Norm 和 Weighted UniFrac 越低越好。",
        "",
        "Profile 结果只使用工具原生输出的丰度文件。",
        "如果某工具没有原生 profile 输出，它不会出现在 profile 结果表中。",
        "",
        "### 字段说明",
        "",
        "- `Elapsed (s)`：运行耗时，单位为秒。",
        "- `Max RSS (GB)`：最大常驻内存，单位为 GB。",
        "- `Total Size (GB)`：输入文件或参考库文件的十进制 GB 大小。",
        "- `Q30 (%)`：FASTQ 碱基质量分数不低于 30 的比例。",
        "",
        "## 工具说明",
        "",
        "### Bracken",
        "",
        "Bracken 基于 Kraken2 输出进行丰度重估，不提供逐读段分类结果。",
        "因此 Bracken 只出现在 profile 结果表中。",
        "运行时间和内存按完整 `kraken2 -> bracken` 流程统计。",
        "",
        "Bracken 的数据库需要额外生成 `database100mers.kmer_distrib`。",
        "该步骤使用 Bracken 默认读长设置 `read_len=100`。",
        "",
        "### ganon2",
        "",
        "ganon2 在 CAMI II Strain Madness short reads 上使用固定批次（fixed-size batch）运行，当前为每批 2 个样本。",
        "原因是默认期望最大化（expectation maximization, EM）重分配在全量拼接输入上超过可用内存。",
        "该设置保留 ganon2 默认的多匹配处理方式，但 EM 的估计范围为每个 batch 内部。",
        "结果表中这一路径仍对应 CAMI II Strain Madness short-read dataset；`Samples` 表示完成的 batch 数。",
        "",
        "ganon2 不报告 PRJNA637878 single-read lane 结果。",
        "该数据集在默认 EM 重分配阶段峰值内存过高，在当前服务器上未能完成运行。",
        "失败样本的 Max RSS 约为 `495..596 GB`，进程在 reassigning reads 阶段被系统终止。",
        "",
        "### Centrifuger",
        "",
        "Centrifuger 的逐读段分类使用默认推断参数。",
        "Profile 结果由 `centrifuger-quant` 从 classify 输出生成，并使用 CAMI profile 格式导出。",
        "",
        "### Taxor",
        "",
        "Taxor 条目来自早期运行，与当前跨工具汇总的运行批次不同。",
        "",
        "## 软件版本",
        "",
        "- kraken2: 2.1.3",
        "- bracken: 2.9",
        "- centrifuger: 1.1.0-r291",
        "- ganon2: 2.1.0",
        "- sylph: 0.8.1",
        "- taxor: 0.1.3（SeqAn 3.4.0-rc.1）",
        "",
    ]


def _results_readme_header() -> list[str]:
    return [
        "# Benchmark 数据与结果说明",
        "",
        "本文件说明 ChimeraBenchmark 的数据、参考库、评估任务、指标和结果表。",
        "它面向希望复核数据来源、评估口径和结果的读者。",
        "",
        "## 目录",
        "",
        "- [数据与结果导航](#数据与结果导航)",
        "- [结果使用说明](#结果使用说明)",
        "- [Reference Database Summary](#reference-database-summary)",
        "- [Benchmark Dataset Summary](#benchmark-dataset-summary)",
        "- [评估任务与口径](#评估任务与口径)",
        "- [数据集特定口径](#数据集特定口径)",
        "- [核心结果汇总](#核心结果汇总)",
        "- [真实队列 clade-level 实验](#真实队列-clade-level-实验)",
        "- [指标说明](#指标说明)",
        "- [工具说明](#工具说明)",
        "- [软件版本](#软件版本)",
        "",
        "## 数据与结果导航",
        "",
        "- `results/builds/README.md`：数据库构建成本。",
        "- `results/classify/README.md`：逐读段分类准确性。",
        "- `results/profile/README.md`：样本丰度画像准确性。",
        "- `results/real/README.md`：真实队列 clade-level 实验结果。",
        "- `results/paper_run_manifest.tsv`：汇总结果所对应的运行、工具版本和完成状态。",
        "- `results/*/summary.tsv`：build、classify 和 profile 的机器可读汇总表。",
        "- `resources/reports/db_catalog.tsv`：Reference Database Summary 的 TSV 版本。",
        "- `resources/reports/dataset_catalog.tsv`：Benchmark Dataset Summary 的 TSV 版本。",
        "",
        "默认报告的分类层级为 `species` 和 `genus`。",
        "",
        "## 结果使用说明",
        "",
        "- 默认 species/genus benchmark 中所有工具使用相同线程数，按 `threads=32` 解释；Fna C2 真实队列实验统一使用 128 threads，并在其 README 中单独记录。",
        "- 除必要的输出格式选项外，分类和丰度推断使用工具默认参数。",
        "- 评估统一使用 NCBI 分类体系（NCBI taxonomy）快照 `resources/taxonomy/ncbi_20260408/`。",
        "- `ganon` 在结果表中显示为 `ganon2`，对应 ganon2 软件版本。",
        "- Taxor 条目来自较早的运行批次；对应版本和设置见各数据集结果。",
        "",
    ]


def write_results_readme(
    results_root: Path,
    *,
    build_rows: list[dict[str, Any]],
    dataset_rows: list[dict[str, Any]],
) -> None:
    lines = _results_readme_header() + [
        "## Reference Database Summary",
        "",
        "本表概括 benchmark 使用的参考库规模。",
        "它描述参考库本体，而不是某个工具构建出的数据库目录。",
        "统计来源为 reference target TSV 和 NCBI assembly summary 元数据。",
        "",
    ]
    lines.extend(
        _markdown_table(
            build_rows,
            [
                ("Dataset Name", "dataset_name"),
                ("Total Size (GB)", "total_size_gb"),
                ("Total Sequences", "total_sequences"),
                ("Base Pairs (bp)", "base_pairs_bp"),
                ("Assemblies", "assemblies"),
                ("Species Count", "species_count"),
                ("Source", "source"),
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Benchmark Dataset Summary",
            "",
            "本表概括 benchmark 使用的实际输入文件。",
            "统计来自 `configs/datasets/*.yaml`、输入序列文件扫描或数据集随附元数据。",
            "",
        ]
    )
    lines.extend(
        _markdown_table(
            dataset_rows,
            [
                ("Dataset Name", "dataset_name"),
                ("Total Size (GB)", "total_size_gb"),
                ("Samples", "samples"),
                ("Input Type", "input_type"),
                ("Reads / Contigs", "reads_or_contigs"),
                ("Base Pairs (bp)", "base_pairs_bp"),
                ("Mean Length (bp)", "mean_length_bp"),
                ("N50 (bp)", "n50_bp"),
                ("GC (%)", "gc_percent"),
                ("Q30 (%)", "q30_percent"),
                ("Truth", "truth"),
            ],
        )
    )
    lines.extend(_results_readme_tail())
    results_root.mkdir(parents=True, exist_ok=True)
    (results_root / "README.md").write_text("\n".join(lines))


def write_catalog_outputs(
    *,
    config_root: Path,
    results_root: Path,
    resources_root: Path,
    progress: bool = False,
) -> dict[str, list[dict[str, Any]]]:
    cache_path = resources_root / "cache" / "catalog_cache.json"
    dataset_rows = collect_dataset_rows(
        config_root=config_root,
        cache_path=cache_path,
        progress=progress,
    )
    build_rows = collect_build_rows(config_root=config_root, cache_path=cache_path)
    write_tsv(
        resources_root / "reports" / "dataset_catalog.tsv",
        dataset_rows,
        [
            "dataset_name",
            "total_size_gb",
            "samples",
            "input_type",
            "reads_or_contigs",
            "base_pairs_bp",
            "mean_length_bp",
            "n50_bp",
            "gc_percent",
            "q30_percent",
            "truth",
        ],
    )
    write_tsv(
        resources_root / "reports" / "db_catalog.tsv",
        build_rows,
        [
            "dataset_name",
            "total_size_gb",
            "total_sequences",
            "base_pairs_bp",
            "assemblies",
            "species_count",
            "source",
        ],
    )
    write_results_readme(results_root, build_rows=build_rows, dataset_rows=dataset_rows)
    return {"build_rows": build_rows, "dataset_rows": dataset_rows}
