from __future__ import annotations

import argparse
import os
from pathlib import Path


def _load_nodes(path: Path) -> dict[int, int]:
    parent: dict[int, int] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 2:
                continue
            try:
                taxid = int(parts[0])
                parent_taxid = int(parts[1])
            except ValueError:
                continue
            parent[taxid] = parent_taxid
    return parent


def _load_scientific_names(path: Path) -> dict[int, str]:
    names: dict[int, str] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            parts = [p.strip() for p in raw.split("|")]
            if len(parts) < 4:
                continue
            if parts[3] != "scientific name":
                continue
            try:
                taxid = int(parts[0])
            except ValueError:
                continue
            name = parts[1]
            if name:
                names.setdefault(taxid, name)
    return names


def _build_lineage_fn(
    *,
    parent: dict[int, int],
    sci_name: dict[int, str],
) -> callable:
    cache: dict[int, tuple[str, str]] = {}

    def lineage(taxid: int) -> tuple[str, str]:
        cached = cache.get(taxid)
        if cached is not None:
            return cached

        seen: set[int] = set()
        ids: list[int] = []
        current = taxid
        while current not in seen:
            seen.add(current)
            ids.append(current)
            p = parent.get(current)
            if p is None or p == current:
                break
            current = p
        ids.reverse()

        names = [sci_name.get(t, str(t)) for t in ids]
        tax_str = "|".join(names) if names else "root"
        tax_id_str = "|".join(str(t) for t in ids) if ids else "1"
        cache[taxid] = (tax_str, tax_id_str)
        return cache[taxid]

    return lineage


def _inspect_header(path: Path) -> tuple[list[str] | None, dict[str, int]]:
    header_cols: list[str] | None = None
    col_index: dict[str, int] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            if not raw.startswith("#"):
                break
            line = raw.strip()
            if not line:
                continue
            cols = [c.lstrip("#") for c in line.split("\t")]
            header_cols = cols
            col_index = {c.lstrip("#"): i for i, c in enumerate(cols)}
            break
    return header_cols, col_index


def needs_tax_fix(path: Path) -> bool:
    header_cols, col_index = _inspect_header(path)
    if not header_cols:
        return True
    if "TAX_STR" not in col_index or "TAX_ID_STR" not in col_index or "TAXID" not in col_index:
        return True
    i_tax_str = col_index["TAX_STR"]
    i_tax_id_str = col_index["TAX_ID_STR"]

    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= max(i_tax_str, i_tax_id_str):
                return True
            return (not parts[i_tax_str]) or (not parts[i_tax_id_str])
    return True


def fix_taxor_search_file(
    *,
    search_file: Path,
    nodes_dmp: Path | None,
    names_dmp: Path | None,
    force: bool,
) -> None:
    if not force and not needs_tax_fix(search_file):
        return

    parent: dict[int, int] = {}
    sci_name: dict[int, str] = {}
    if nodes_dmp and nodes_dmp.exists():
        parent = _load_nodes(nodes_dmp)
    if names_dmp and names_dmp.exists():
        sci_name = _load_scientific_names(names_dmp)

    lineage = _build_lineage_fn(parent=parent, sci_name=sci_name)

    header_cols, col_index = _inspect_header(search_file)
    if not header_cols:
        raise ValueError(f"taxor search file missing header: {search_file}")

    required = ("ACCESSION", "REFERENCE_NAME", "TAXID", "TAX_STR", "TAX_ID_STR")
    for col in required:
        if col not in col_index:
            raise ValueError(f"taxor search header missing {col}: {search_file}")

    i_accession = col_index["ACCESSION"]
    i_ref_name = col_index["REFERENCE_NAME"]
    i_taxid = col_index["TAXID"]
    i_tax_str = col_index["TAX_STR"]
    i_tax_id_str = col_index["TAX_ID_STR"]
    pad_len = max(col_index.values()) + 1

    tmp_path = search_file.with_suffix(search_file.suffix + ".tmp")
    with search_file.open("r", encoding="utf-8", errors="ignore") as fin, tmp_path.open(
        "w", encoding="utf-8"
    ) as fout:
        for raw in fin:
            if raw.startswith("#"):
                fout.write(raw)
                continue
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < pad_len:
                parts.extend([""] * (pad_len - len(parts)))

            if not parts[i_ref_name] and parts[i_accession]:
                parts[i_ref_name] = parts[i_accession]

            taxid = None
            if parts[i_taxid].isdigit():
                taxid = int(parts[i_taxid])

            if not parts[i_tax_str] or not parts[i_tax_id_str]:
                if taxid is None:
                    tax_str, tax_id_str = ("root", "1")
                else:
                    tax_str, tax_id_str = lineage(taxid)
                    if not tax_str:
                        tax_str = str(taxid)
                    if not tax_id_str:
                        tax_id_str = str(taxid)
                if not parts[i_tax_str]:
                    parts[i_tax_str] = tax_str
                if not parts[i_tax_id_str]:
                    parts[i_tax_id_str] = tax_id_str

            fout.write("\t".join(parts) + "\n")

    os.replace(tmp_path, search_file)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--search-file", required=True)
    ap.add_argument("--nodes-dmp")
    ap.add_argument("--names-dmp")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    fix_taxor_search_file(
        search_file=Path(args.search_file),
        nodes_dmp=Path(args.nodes_dmp) if args.nodes_dmp else None,
        names_dmp=Path(args.names_dmp) if args.names_dmp else None,
        force=bool(args.force),
    )


if __name__ == "__main__":
    main()

