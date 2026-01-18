from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path


FASTA_EXTS = (
    ".fna.gz",
    ".fa.gz",
    ".fasta.gz",
    ".fna",
    ".fa",
    ".fasta",
)


def _strip_fasta_ext(name: str) -> str:
    for ext in FASTA_EXTS:
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def prepare_library(*, target_tsv: Path, out_dir: Path, force: bool) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    library_path = out_dir / "library.fna"
    tmp_library = out_dir / "library.fna.tmp"

    if library_path.exists() and not force:
        raise SystemExit(f"library file already exists: {library_path} (use --force to overwrite)")

    if tmp_library.exists():
        tmp_library.unlink()

    processed = 0
    skipped = 0

    with target_tsv.open("r", encoding="utf-8", errors="ignore") as f_in, tmp_library.open(
        "w", encoding="utf-8"
    ) as lib_out:
        for raw in f_in:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            genome_path = Path(parts[0].strip())
            taxid = parts[1].strip()
            if not genome_path.exists():
                skipped += 1
                raise FileNotFoundError(genome_path)

            basename = _strip_fasta_ext(genome_path.name)
            map_path = out_dir / f"prelim_map_{basename}.txt"

            with map_path.open("w", encoding="utf-8") as map_out, _open_text(genome_path) as fh:
                for fasta_line in fh:
                    if fasta_line.startswith(">"):
                        header = fasta_line[1:].rstrip("\n")
                        if not header:
                            continue
                        parts_h = header.split(None, 1)
                        seq_id = parts_h[0]
                        rest = (" " + parts_h[1]) if len(parts_h) > 1 else ""
                        if "|kraken:taxid|" not in seq_id:
                            seq_id = f"{seq_id}|kraken:taxid|{taxid}"
                        lib_out.write(f">{seq_id}{rest}\n")
                        map_out.write(f"TAXID\t{seq_id}\t{taxid}\n")
                    else:
                        lib_out.write(fasta_line)

            processed += 1
            if processed % 100 == 0:
                print(f"[kraken2] prepared {processed} genomes (skipped={skipped})")

    os.replace(tmp_library, library_path)
    print(
        f"kraken2 library prep done: processed={processed} written={processed - skipped} skipped={skipped}"
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-tsv", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    prepare_library(
        target_tsv=Path(args.target_tsv),
        out_dir=Path(args.out_dir),
        force=bool(args.force),
    )


if __name__ == "__main__":
    main()

