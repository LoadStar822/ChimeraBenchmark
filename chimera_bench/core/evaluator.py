from __future__ import annotations

from pathlib import Path


def summarize_classify_tsv(path: Path) -> dict:
    total = 0
    unclassified = 0
    taxids = set()
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            total += 1
            parts = line.split("\t")
            if len(parts) < 2:
                unclassified += 1
                continue
            tokens = [t for t in parts[1:] if t]
            if not tokens or tokens[0] == "unclassified":
                unclassified += 1
                continue
            first = tokens[0]
            if ":" in first:
                taxid = first.split(":", 1)[0]
                if taxid and taxid != "unclassified":
                    taxids.add(taxid)
            else:
                unclassified += 1
    return {
        "total_reads": total,
        "unclassified_reads": unclassified,
        "classified_reads": max(0, total - unclassified),
        "unique_taxids": len(taxids),
    }
