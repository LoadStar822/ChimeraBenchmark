from __future__ import annotations

import argparse
import json
import re
import shutil
import tarfile
from pathlib import Path
from typing import Iterable


ASSEMBLY_RE = re.compile(r"(GC[AF]_\d+\.\d+)")
INTERMEDIATE_RE = re.compile(r"\d+_(?:diffIdx|info)")
TAXONOMY_FILES = ("nodes.dmp", "names.dmp", "merged.dmp")


def _signatures(paths: Iterable[Path]) -> list[dict[str, int | str]]:
    signatures = []
    for path in paths:
        stat = path.stat()
        signatures.append({"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    return signatures


def prepare_query(*, inputs: list[Path], out_path: Path) -> None:
    if len(inputs) < 2:
        raise ValueError("query preparation requires at least two inputs")
    signatures = _signatures(inputs)
    manifest_path = Path(f"{out_path}.manifest.json")
    if out_path.exists() and manifest_path.exists():
        if json.loads(manifest_path.read_text()) == signatures:
            return

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = Path(f"{out_path}.tmp")
    try:
        with tmp_path.open("wb") as fout:
            for path in inputs:
                with path.open("rb") as fin:
                    shutil.copyfileobj(fin, fout, length=16 * 1024 * 1024)
        tmp_path.replace(out_path)
        manifest_path.write_text(json.dumps(signatures, indent=2))
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def prepare_db(
    *,
    target_tsv: Path,
    fasta_list: Path,
    assembly_map: Path,
    taxonomy_dir: Path,
    db_dir: Path,
) -> None:
    files: list[str] = []
    mappings: list[tuple[str, str]] = []
    seen: dict[str, str] = {}
    fasta_list.parent.mkdir(parents=True, exist_ok=True)
    alias_root = fasta_list.parent / "path_aliases"
    with target_tsv.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"invalid target TSV row {line_number}: {line}")
            fasta_path, taxid = parts[0], parts[1]
            path = Path(fasta_path)
            if not path.is_file():
                raise FileNotFoundError(path)
            assemblies = tuple(dict.fromkeys(ASSEMBLY_RE.findall(path.name)))
            if not assemblies:
                raise ValueError(f"assembly accession missing from target path: {path}")
            assembly = assemblies[0]
            previous = seen.get(assembly)
            if previous is not None and previous != taxid:
                raise ValueError(f"conflicting taxids for {assembly}: {previous}, {taxid}")
            if previous is not None:
                continue
            seen[assembly] = taxid
            mappings.append((assembly, taxid))
            listed_path = path.absolute()
            if len(assemblies) > 1:
                alias_root.mkdir(exist_ok=True)
                suffix = ".fna.gz" if path.name.endswith(".gz") else ".fna"
                listed_path = alias_root / f"{assembly}{suffix}"
                if listed_path.exists() or listed_path.is_symlink():
                    if listed_path.resolve() != path.resolve():
                        raise ValueError(f"conflicting FASTA paths for {assembly}")
                else:
                    listed_path.symlink_to(path.absolute())
            files.append(str(listed_path))

    if not files:
        raise ValueError(f"target TSV is empty: {target_tsv}")

    fasta_list.write_text("".join(f"{path}\n" for path in files))
    assembly_map.write_text("".join(f"{assembly}\t{taxid}\n" for assembly, taxid in mappings))

    taxonomy_out = db_dir / "taxonomy"
    taxonomy_out.mkdir(parents=True, exist_ok=True)
    for name in TAXONOMY_FILES:
        source = taxonomy_dir / name
        target = taxonomy_out / name
        if source.is_file():
            shutil.copy2(source, target)
            continue
        archive_path = taxonomy_dir / "taxdump.tar.gz"
        if name != "merged.dmp" or not archive_path.is_file():
            raise FileNotFoundError(source)
        with tarfile.open(archive_path, "r:gz") as archive:
            member = archive.extractfile(name)
            if member is None:
                raise FileNotFoundError(f"{name} in {archive_path}")
            with member, target.open("wb") as output:
                shutil.copyfileobj(member, output)


def cleanup_db(*, db_dir: Path) -> tuple[int, int]:
    removed = 0
    removed_bytes = 0
    for path in db_dir.iterdir():
        if path.is_file() and INTERMEDIATE_RE.fullmatch(path.name):
            removed_bytes += path.stat().st_size
            path.unlink()
            removed += 1
    return removed, removed_bytes


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    query_parser = subparsers.add_parser("query")
    query_parser.add_argument("--out", required=True)
    query_parser.add_argument("--inputs", nargs="+", required=True)

    db_parser = subparsers.add_parser("db")
    db_parser.add_argument("--target-tsv", required=True)
    db_parser.add_argument("--fasta-list", required=True)
    db_parser.add_argument("--assembly-map", required=True)
    db_parser.add_argument("--taxonomy-dir", required=True)
    db_parser.add_argument("--db-dir", required=True)
    cleanup_parser = subparsers.add_parser("cleanup")
    cleanup_parser.add_argument("--db-dir", required=True)

    args = parser.parse_args()
    if args.command == "query":
        prepare_query(inputs=[Path(path) for path in args.inputs], out_path=Path(args.out))
    elif args.command == "db":
        prepare_db(
            target_tsv=Path(args.target_tsv),
            fasta_list=Path(args.fasta_list),
            assembly_map=Path(args.assembly_map),
            taxonomy_dir=Path(args.taxonomy_dir),
            db_dir=Path(args.db_dir),
        )
    else:
        cleanup_db(db_dir=Path(args.db_dir))


if __name__ == "__main__":
    main()
