from __future__ import annotations

import argparse
import re
import sys
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


def load_nodes(path: Path) -> dict[str, tuple[str, str]]:
    nodes: dict[str, tuple[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            parts = [p.strip() for p in raw.split("|")]
            if len(parts) >= 3 and parts[0]:
                nodes[parts[0]] = (parts[1], parts[2])
    return nodes


def write_ganon_input(*, target_tsv: Path, nodes_path: Path, out_path: Path) -> dict[str, int]:
    """Write a ganon build-custom --input-file (file, target, node) from the library definition.

    Targets are assembly accessions and nodes are the library taxids, so `--level species`
    groups genomes by the species of their library taxid. Taxids missing from the
    taxonomy are skipped and counted.
    """
    nodes = load_nodes(nodes_path)
    counts = {"written": 0, "no_accession": 0, "taxid_not_in_taxonomy": 0}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with target_tsv.open("r", encoding="utf-8", errors="ignore") as fin, out_path.open("w", encoding="utf-8") as fout:
        for raw in fin:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            filepath = fields[0]
            taxid = fields[1].strip() if len(fields) > 1 else ""
            accession = _accession_from_path(filepath)
            if not accession:
                counts["no_accession"] += 1
                continue
            if taxid not in nodes:
                counts["taxid_not_in_taxonomy"] += 1
                print(f"skip {accession}: taxid {taxid} not in {nodes_path}", file=sys.stderr)
                continue
            fout.write(f"{filepath}\t{accession}\t{taxid}\n")
            counts["written"] += 1
    return counts


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-tsv", required=True)
    ap.add_argument("--nodes", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    counts = write_ganon_input(target_tsv=Path(args.target_tsv), nodes_path=Path(args.nodes), out_path=Path(args.out))
    print(
        f"Wrote {counts['written']} records to {args.out}; skipped {counts['no_accession']} without accession "
        f"and {counts['taxid_not_in_taxonomy']} with taxids missing from the taxonomy"
    )


if __name__ == "__main__":
    main()
