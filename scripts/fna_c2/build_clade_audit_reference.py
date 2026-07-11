#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path
from typing import Iterable, TextIO


def open_maybe_gzip(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def parse_fasta(path: Path) -> Iterable[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with open_maybe_gzip(path) as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].strip()
                chunks = []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def safe_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", value).strip("_") or "na"


def build_clade_references(reference_manifest: Path, out_dir: Path) -> None:
    rows = read_tsv(reference_manifest)
    required = {"path", "taxid", "paper_clade_projection", "reference_name"}
    if not rows:
        raise ValueError(f"reference manifest is empty: {reference_manifest}")
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"reference manifest is missing columns: {sorted(missing)}")

    out_dir.mkdir(parents=True, exist_ok=True)
    reference_dir = out_dir / "clades"
    reference_dir.mkdir(parents=True, exist_ok=True)
    seen: dict[tuple[str, str], dict[str, object]] = {}
    records: list[dict[str, object]] = []
    duplicates: list[dict[str, object]] = []
    summary: dict[str, dict[str, int]] = defaultdict(
        lambda: {
            "input_records": 0,
            "input_bases": 0,
            "retained_records": 0,
            "retained_bases": 0,
            "duplicate_records": 0,
            "duplicate_bases": 0,
        }
    )

    with ExitStack() as stack:
        clade_handles: dict[str, TextIO] = {}
        for manifest_index, row in enumerate(rows, start=1):
            source_path = Path(row["path"])
            clade = row["paper_clade_projection"]
            taxid = row["taxid"]
            reference_name = row["reference_name"]
            if clade not in clade_handles:
                clade_handles[clade] = stack.enter_context(
                    (reference_dir / f"{safe_token(clade)}.fna").open("w")
                )
            output = clade_handles[clade]
            for sequence_index, (original_name, sequence) in enumerate(
                parse_fasta(source_path), start=1
            ):
                if not sequence:
                    raise ValueError(f"empty FASTA sequence in {source_path}: {original_name}")
                digest = hashlib.sha256(sequence.encode()).hexdigest()
                stats = summary[clade]
                stats["input_records"] += 1
                stats["input_bases"] += len(sequence)
                key = (clade, digest)
                if key in seen:
                    kept = seen[key]
                    stats["duplicate_records"] += 1
                    stats["duplicate_bases"] += len(sequence)
                    duplicates.append(
                        {
                            "clade": clade,
                            "sha256": digest,
                            "sequence_length": len(sequence),
                            "duplicate_source_path": str(source_path),
                            "duplicate_reference_name": reference_name,
                            "duplicate_original_name": original_name,
                            "kept_target_id": kept["target_id"],
                            "kept_source_path": kept["source_path"],
                            "kept_reference_name": kept["reference_name"],
                            "kept_original_name": kept["original_name"],
                        }
                    )
                    continue

                target_id = (
                    f"ref{manifest_index:05d}_{sequence_index:05d}"
                    f"|clade={safe_token(clade)}"
                    f"|taxid={safe_token(taxid)}"
                    f"|ref={safe_token(reference_name)}"
                    f"|orig={safe_token(original_name)}"
                )
                record = {
                    "target_id": target_id,
                    "clade": clade,
                    "taxid": taxid,
                    "reference_name": reference_name,
                    "original_name": original_name,
                    "source_path": str(source_path),
                    "sequence_length": len(sequence),
                    "sha256": digest,
                }
                seen[key] = record
                records.append(record)
                stats["retained_records"] += 1
                stats["retained_bases"] += len(sequence)
                output.write(f">{target_id}\n")
                for offset in range(0, len(sequence), 80):
                    output.write(sequence[offset : offset + 80] + "\n")

    summary_rows = [
        {"clade": clade, **values}
        for clade, values in sorted(summary.items())
    ]
    write_tsv(
        out_dir / "reference_records.tsv",
        records,
        [
            "target_id",
            "clade",
            "taxid",
            "reference_name",
            "original_name",
            "source_path",
            "sequence_length",
            "sha256",
        ],
    )
    write_tsv(
        out_dir / "reference_duplicates.tsv",
        duplicates,
        [
            "clade",
            "sha256",
            "sequence_length",
            "duplicate_source_path",
            "duplicate_reference_name",
            "duplicate_original_name",
            "kept_target_id",
            "kept_source_path",
            "kept_reference_name",
            "kept_original_name",
        ],
    )
    write_tsv(
        out_dir / "reference_summary.tsv",
        summary_rows,
        [
            "clade",
            "input_records",
            "input_bases",
            "retained_records",
            "retained_bases",
            "duplicate_records",
            "duplicate_bases",
        ],
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Build exact-deduplicated per-clade FASTA references for Fna read audits."
    )
    parser.add_argument("--reference-manifest", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    build_clade_references(args.reference_manifest, args.out_dir)
    print(f"wrote deduplicated clade references to {args.out_dir}")


if __name__ == "__main__":
    main()
