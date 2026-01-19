from __future__ import annotations

import argparse
from pathlib import Path


RANK_MAP = {
    "D": "domain",
    "P": "phylum",
    "C": "class",
    "O": "order",
    "F": "family",
    "G": "genus",
    "S": "species",
}


def _rank_name(raw: str) -> str:
    key = raw.strip()
    if not key:
        return "species"
    if key in RANK_MAP:
        return RANK_MAP[key]
    return key.lower()


def convert(*, bracken_tsv: Path, out_path: Path) -> None:
    if not bracken_tsv.exists():
        raise FileNotFoundError(bracken_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = None
    idx_taxid = idx_level = idx_est = None
    rows: list[tuple[int, str, float]] = []

    with bracken_tsv.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                lower = [p.strip().lower() for p in parts]
                if "taxonomy_id" in lower and "new_est_reads" in lower:
                    header = lower
                    idx_taxid = header.index("taxonomy_id")
                    idx_est = header.index("new_est_reads")
                    if "taxonomy_lvl" in header:
                        idx_level = header.index("taxonomy_lvl")
                    elif "taxonomy_level" in header:
                        idx_level = header.index("taxonomy_level")
                    continue
                header = []

            if header:
                if idx_taxid is None or idx_est is None:
                    continue
                if len(parts) <= max(idx_taxid, idx_est):
                    continue
                taxid_text = parts[idx_taxid].strip()
                level_text = (
                    parts[idx_level].strip() if idx_level is not None and len(parts) > idx_level else ""
                )
                est_text = parts[idx_est].strip()
            else:
                # Fallback for header-less files:
                # name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads
                if len(parts) < 6:
                    continue
                taxid_text = parts[1].strip()
                level_text = parts[2].strip()
                est_text = parts[5].strip()

            if not taxid_text.isdigit():
                continue
            taxid = int(taxid_text)
            if taxid <= 0:
                continue
            try:
                est = float(est_text)
            except ValueError:
                continue
            if est <= 0:
                continue
            rows.append((taxid, _rank_name(level_text), est))

    with out_path.open("w", encoding="utf-8") as out:
        out.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
        for taxid, rank, value in rows:
            out.write(f"{taxid}\t{rank}\t-\t-\t{value}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    convert(bracken_tsv=Path(args.input), out_path=Path(args.out))


if __name__ == "__main__":
    main()

