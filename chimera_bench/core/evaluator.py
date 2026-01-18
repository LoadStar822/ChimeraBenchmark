from __future__ import annotations

from pathlib import Path


def summarize_classify_tsv(path: Path) -> dict:
    total = 0
    unclassified = 0
    taxids: set[int] = set()
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#") or line.startswith("@"):
                continue
            total += 1
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            tokens = [t for t in parts[1:] if t]
            if not tokens:
                unclassified += 1
                continue
            first = tokens[0]
            taxid_str = first.split(":", 1)[0].split("|", 1)[0].strip()
            if taxid_str in {"unclassified", "-", "0", "NA"}:
                unclassified += 1
                continue
            try:
                taxid = int(taxid_str)
            except ValueError:
                unclassified += 1
                continue
            taxids.add(taxid)
    return {
        "total_reads": total,
        "unclassified_reads": unclassified,
        "classified_reads": max(0, total - unclassified),
        "unique_taxids": len(taxids),
    }


def summarize_ganon_tre(path: Path) -> dict:
    classified = None
    unclassified = None
    taxids = set()
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rank = parts[0]
            taxid = parts[1]
            try:
                count = int(float(parts[-2]))
            except (ValueError, IndexError):
                continue
            if rank == "unclassified":
                unclassified = count
            elif rank == "root":
                classified = count
            if (
                taxid
                and taxid not in {"-", "0"}
                and rank not in {"unclassified", "root"}
                and count > 0
            ):
                taxids.add(taxid)
    if classified is None and unclassified is None:
        return {}
    classified_reads = classified or 0
    unclassified_reads = unclassified or 0
    total_reads = classified_reads + unclassified_reads
    return {
        "total_reads": total_reads,
        "unclassified_reads": unclassified_reads,
        "classified_reads": classified_reads,
        "unique_taxids": len(taxids),
    }
