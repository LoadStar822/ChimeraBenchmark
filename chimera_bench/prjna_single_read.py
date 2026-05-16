from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
from pathlib import Path
from typing import Any, Iterable

import yaml


DEFAULT_SOURCE_ROOT = Path(
    "/mnt/sda/CommonData/TaxonomyDataset/True/真实实验/群落结构线/当前主线/"
    "prjna637878_culture_supported_truth_20260411"
)
DEFAULT_SOURCE_CONFIG = Path("configs/datasets/prjna637878-supported19.yaml")
DEFAULT_OUT_ROOT = Path("resources/derived/prjna637878-supported19-single-read")
DEFAULT_AUDIT = Path("resources/reports/prjna637878_single_read_input_audit.tsv")


def _load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise ValueError(f"invalid YAML object: {path}")
    return data


def _load_manifest_by_sample(source_root: Path) -> dict[str, dict[str, str]]:
    path = source_root / "meta" / "PRJNA637878_supported19_sample_manifest.tsv"
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        return {row["sample_id"]: row for row in csv.DictReader(fh, delimiter="\t")}


def _load_runinfo(source_root: Path, biosample: str) -> dict[str, str]:
    if not biosample:
        return {}
    path = source_root / "meta" / f"{biosample}.runinfo.csv"
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    return rows[0] if rows else {}


def _mate_suffix(mate: str) -> bytes:
    if mate == "R1":
        return b"__mate1"
    if mate == "R2":
        return b"__mate2"
    raise ValueError(f"unsupported mate: {mate}")


def mate_marked_read_id(read_id: str, mate: str) -> str:
    base = read_id.strip().split()[0]
    if base.startswith("@"):
        base = base[1:]
    if base.endswith("/1") or base.endswith("/2"):
        base = base[:-2]
    if base.endswith("__mate1") or base.endswith("__mate2"):
        base = base.rsplit("__mate", 1)[0]
    return f"{base}{_mate_suffix(mate).decode('ascii')}"


def _mate_marked_header(header: bytes, mate: str) -> bytes:
    if not header.startswith(b"@"):
        raise ValueError(f"FASTQ header does not start with @: {header[:80]!r}")
    token = header[1:].strip().split(None, 1)[0]
    if token.endswith(b"/1") or token.endswith(b"/2"):
        token = token[:-2]
    if token.endswith(b"__mate1") or token.endswith(b"__mate2"):
        token = token.rsplit(b"__mate", 1)[0]
    return b"@" + token + _mate_suffix(mate) + b"\n"


def _gzip_ok(path: Path) -> bool:
    return subprocess.run(
        ["gzip", "-t", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    ).returncode == 0


def rewrite_fastq_with_mate(
    *,
    input_path: Path,
    output_path: Path,
    mate: str,
    pigz_threads: int,
    force: bool,
) -> dict[str, int]:
    if output_path.exists() and not force:
        if not _gzip_ok(output_path):
            raise ValueError(f"existing gzip failed validation: {output_path}")
        return {"records": -1, "bases": -1, "size_bytes": output_path.stat().st_size}

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_name(output_path.name + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()

    records = 0
    bases = 0
    with gzip.open(input_path, "rb") as in_fh, tmp_path.open("wb") as raw_out:
        proc = subprocess.Popen(
            ["pigz", "-c", "-p", str(pigz_threads)],
            stdin=subprocess.PIPE,
            stdout=raw_out,
        )
        assert proc.stdin is not None
        try:
            while True:
                header = in_fh.readline()
                if not header:
                    break
                seq = in_fh.readline()
                plus = in_fh.readline()
                qual = in_fh.readline()
                if not (seq and plus and qual):
                    raise ValueError(f"truncated FASTQ record in {input_path}")
                proc.stdin.write(_mate_marked_header(header, mate))
                proc.stdin.write(seq)
                proc.stdin.write(plus)
                proc.stdin.write(qual)
                records += 1
                bases += len(seq.rstrip(b"\r\n"))
        finally:
            proc.stdin.close()
            rc = proc.wait()
        if rc != 0:
            raise RuntimeError(f"pigz failed for {output_path} with rc={rc}")

    if not _gzip_ok(tmp_path):
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(f"generated gzip failed validation: {tmp_path}")
    tmp_path.replace(output_path)
    return {"records": records, "bases": bases, "size_bytes": output_path.stat().st_size}


def write_single_read_truth(*, species_truth: Path, output_path: Path, force: bool) -> int:
    if output_path.exists() and not force:
        return -1
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_name(output_path.name + ".tmp")
    rows = 0
    with species_truth.open("r", encoding="utf-8", newline="") as in_fh, tmp_path.open(
        "w", encoding="utf-8", newline=""
    ) as out_fh:
        reader = csv.DictReader(in_fh, delimiter="\t")
        required = {"read_id", "mate", "species_label"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{species_truth} missing columns: {', '.join(sorted(missing))}")
        writer = csv.DictWriter(
            out_fh,
            fieldnames=["read_id", "species_label", "source_read_id", "mate"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        seen: set[str] = set()
        for row in reader:
            mate = row["mate"].strip()
            if mate not in {"R1", "R2"}:
                continue
            read_id = mate_marked_read_id(row["read_id"], mate)
            if read_id in seen:
                raise ValueError(f"duplicate mate-level truth read_id: {read_id}")
            seen.add(read_id)
            species_label = row["species_label"].strip()
            if not species_label:
                continue
            writer.writerow(
                {
                    "read_id": read_id,
                    "species_label": species_label,
                    "source_read_id": row["read_id"].strip(),
                    "mate": mate,
                }
            )
            rows += 1
    tmp_path.replace(output_path)
    return rows


def _selected_samples(samples: Iterable[dict[str, Any]], selected: set[str]) -> list[dict[str, Any]]:
    out = []
    for sample in samples:
        sample_id = str(sample.get("sample_id") or "")
        if not selected or sample_id in selected:
            out.append(sample)
    return out


def build_single_read_assets(
    *,
    source_config: Path,
    source_root: Path,
    out_root: Path,
    audit_path: Path,
    selected: set[str],
    pigz_threads: int,
    force: bool,
    check_only: bool,
) -> list[dict[str, Any]]:
    cfg = _load_yaml(source_config)
    samples = _selected_samples(cfg.get("samples") or [], selected)
    if not samples:
        raise ValueError("no samples selected")

    manifest_by_sample = _load_manifest_by_sample(source_root)
    rows: list[dict[str, Any]] = []
    for sample in samples:
        sample_id = str(sample["sample_id"])
        paired = list(sample.get("paired") or [])
        if len(paired) != 2:
            raise ValueError(f"sample must have exactly two paired inputs: {sample_id}")

        source_r1 = Path(paired[0])
        source_r2 = Path(paired[1])
        species_truth = source_root / "truth" / f"{sample_id}.species_truth.tsv"
        out_r1 = out_root / "fastq" / f"{sample_id}_R1.single.fastq.gz"
        out_r2 = out_root / "fastq" / f"{sample_id}_R2.single.fastq.gz"
        out_truth = out_root / "truth" / f"{sample_id}.single_read_truth.tsv"
        for path in (source_r1, source_r2, species_truth):
            if not path.exists():
                raise FileNotFoundError(path)

        manifest_row = manifest_by_sample.get(sample_id, {})
        runinfo = _load_runinfo(source_root, manifest_row.get("biosample", ""))
        if check_only:
            r1_stats = {"records": 0, "bases": 0, "size_bytes": 0}
            r2_stats = {"records": 0, "bases": 0, "size_bytes": 0}
            truth_rows = 0
        else:
            r1_stats = rewrite_fastq_with_mate(
                input_path=source_r1,
                output_path=out_r1,
                mate="R1",
                pigz_threads=pigz_threads,
                force=force,
            )
            r2_stats = rewrite_fastq_with_mate(
                input_path=source_r2,
                output_path=out_r2,
                mate="R2",
                pigz_threads=pigz_threads,
                force=force,
            )
            truth_rows = write_single_read_truth(species_truth=species_truth, output_path=out_truth, force=force)

        rows.append(
            {
                "sample_id": sample_id,
                "source_r1": str(source_r1),
                "source_r2": str(source_r2),
                "out_r1": str(out_r1),
                "out_r2": str(out_r2),
                "truth": str(out_truth),
                "source_r1_records": r1_stats["records"],
                "source_r2_records": r2_stats["records"],
                "single_read_records": (
                    r1_stats["records"] + r2_stats["records"]
                    if r1_stats["records"] >= 0 and r2_stats["records"] >= 0
                    else -1
                ),
                "single_read_bases": (
                    r1_stats["bases"] + r2_stats["bases"]
                    if r1_stats["bases"] >= 0 and r2_stats["bases"] >= 0
                    else -1
                ),
                "truth_rows": truth_rows,
                "run": runinfo.get("Run", ""),
                "runinfo_spots": runinfo.get("spots", ""),
                "runinfo_spots_with_mates": runinfo.get("spots_with_mates", ""),
                "runinfo_unpaired_spots": (
                    int(runinfo["spots"]) - int(runinfo["spots_with_mates"])
                    if runinfo.get("spots") and runinfo.get("spots_with_mates")
                    else ""
                ),
            }
        )

    if not check_only:
        audit_path.parent.mkdir(parents=True, exist_ok=True)
        with audit_path.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        metadata = {
            "source_config": str(source_config),
            "source_root": str(source_root),
            "out_root": str(out_root),
            "audit_path": str(audit_path),
            "samples": [row["sample_id"] for row in rows],
        }
        out_root.mkdir(parents=True, exist_ok=True)
        (out_root / "manifest.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return rows


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Build PRJNA637878 supported19 single-read benchmark assets.")
    ap.add_argument("--source-config", type=Path, default=DEFAULT_SOURCE_CONFIG)
    ap.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    ap.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    ap.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    ap.add_argument("--samples", action="append", default=[])
    ap.add_argument("--pigz-threads", type=int, default=8)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--check-only", action="store_true")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    selected = {sample for group in args.samples for sample in group.split(",") if sample}
    rows = build_single_read_assets(
        source_config=args.source_config,
        source_root=args.source_root,
        out_root=args.out_root,
        audit_path=args.audit,
        selected=selected,
        pigz_threads=args.pigz_threads,
        force=args.force,
        check_only=args.check_only,
    )
    for row in rows:
        print(
            "\t".join(
                [
                    str(row["sample_id"]),
                    str(row["single_read_records"]),
                    str(row["single_read_bases"]),
                    str(row["truth_rows"]),
                    str(row["runinfo_unpaired_spots"]),
                ]
            )
        )


if __name__ == "__main__":
    main()
