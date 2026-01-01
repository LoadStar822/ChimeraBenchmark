from __future__ import annotations

import argparse
import re
from pathlib import Path


def _accession_from_path(path_str: str) -> str | None:
    name = Path(path_str).name
    parts = name.split("_")
    if len(parts) >= 2 and parts[0] in {"GCF", "GCA"}:
        return f"{parts[0]}_{parts[1]}"
    match = re.search(r"(GCF|GCA)_\d+\.\d+", name)
    if match:
        return match.group(0)
    return None


def write_taxor_input(*, assembly_summary: Path, target_tsv: Path, out_path: Path) -> int:
    assembly_map: dict[str, str] = {}
    with assembly_summary.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 20:
                continue
            accession = parts[0]
            taxid = parts[6]
            if accession and taxid:
                assembly_map[accession] = taxid

    out_path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with target_tsv.open("r", encoding="utf-8", errors="ignore") as fin, out_path.open(
        "w",
        encoding="utf-8",
    ) as fout:
        for raw in fin:
            line = raw.strip()
            if not line:
                continue
            filepath = line.split("\t", 1)[0]
            accession = _accession_from_path(filepath)
            if not accession:
                continue
            taxid = assembly_map.get(accession)
            if not taxid:
                continue
            fout.write(f"{accession}\t{taxid}\t{filepath}\n")
            written += 1
    return written


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembly-summary", required=True)
    ap.add_argument("--target-tsv", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    written = write_taxor_input(
        assembly_summary=Path(args.assembly_summary),
        target_tsv=Path(args.target_tsv),
        out_path=Path(args.out),
    )
    print(f"Wrote {written} records to {args.out}")


if __name__ == "__main__":
    main()

