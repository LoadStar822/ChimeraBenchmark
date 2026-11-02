from __future__ import annotations

import argparse
from pathlib import Path


QUANT_COLUMNS = ["taxID", "taxRank", "genomeSize", "numReads", "numUniqueReads", "abundance"]


def convert_centrifuger_output(*, input_path: Path, out_path: Path) -> None:
    """Convert Centrifuger stdout TSV to ChimeraBenchmark classify TSV.

    Centrifuger output columns (from README):
      readID seqID taxID score 2ndBestScore hitLength queryLength numMatches
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    with input_path.open("r", encoding="utf-8", errors="ignore") as fin, out_path.open(
        "w", encoding="utf-8"
    ) as fout:
        for raw in fin:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            read_id = parts[0].strip()
            taxid = parts[2].strip()
            if not read_id:
                continue
            lower_id = read_id.lower()
            if lower_id == "readid" and taxid.lower() == "taxid":
                continue
            if read_id in seen:
                continue
            seen.add(read_id)
            if not taxid or taxid in {"0", "-", "unclassified"}:
                fout.write(f"{read_id}\tunclassified\n")
            else:
                fout.write(f"{read_id}\t{taxid}\n")


def convert_centrifuger_quant(*, quant_path: Path, out_path: Path) -> None:
    """Convert centrifuger-quant output (--output-format 0) to a species CAMI profile.

    Abundances are the species-level read counts (numReads), in which a read
    with k best hits counts 1/k towards each.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with quant_path.open("r", encoding="utf-8", errors="strict") as fin, out_path.open(
        "w", encoding="utf-8"
    ) as fout:
        header = fin.readline().rstrip("\n").split("\t")
        if header[-len(QUANT_COLUMNS):] != QUANT_COLUMNS:
            raise ValueError(f"unexpected centrifuger-quant header in {quant_path}: {header}")
        fout.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
        for raw in fin:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) <= len(QUANT_COLUMNS):
                continue
            # Columns are taken from the end, after the taxon name.
            taxid, rank, count = parts[-6], parts[-5], int(parts[-3])
            if rank == "species" and taxid.isdigit() and taxid != "0" and count > 0:
                fout.write(f"{taxid}\tspecies\t-\t-\t{count}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input")
    ap.add_argument("--out")
    ap.add_argument("--quant")
    ap.add_argument("--profile-out")
    args = ap.parse_args()
    if bool(args.input) != bool(args.out):
        ap.error("--input and --out must be provided together")
    if bool(args.quant) != bool(args.profile_out):
        ap.error("--quant and --profile-out must be provided together")
    if not (args.input or args.quant):
        ap.error("provide --input/--out or --quant/--profile-out")

    if args.input:
        convert_centrifuger_output(input_path=Path(args.input), out_path=Path(args.out))
    if args.quant:
        convert_centrifuger_quant(quant_path=Path(args.quant), out_path=Path(args.profile_out))


if __name__ == "__main__":
    main()

