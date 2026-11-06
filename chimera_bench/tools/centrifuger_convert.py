from __future__ import annotations

import argparse
import math
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


def convert_centrifuger_quant(
    *, quant_path: Path, out_path: Path, profile_method: str = "native-sequence"
) -> None:
    """Convert centrifuger-quant output (--output-format 0) to a species CAMI profile.

    The default estimates sequence abundance as native EM abundance times the
    reported species genome size, normalized over species. This conversion is
    not an inversion of the EM model or an estimate of its posterior read counts.
    The read-counts option retains the previous numReads-based comparison.
    """
    if profile_method not in {"native-sequence", "read-counts"}:
        raise ValueError(f"unknown profile method: {profile_method}")
    values = {}
    with quant_path.open("r", encoding="utf-8", errors="strict") as fin:
        header = fin.readline().rstrip("\n").split("\t")
        if header[-len(QUANT_COLUMNS):] != QUANT_COLUMNS:
            raise ValueError(f"unexpected centrifuger-quant header in {quant_path}: {header}")
        for lineno, raw in enumerate(fin, 2):
            if not raw.strip():
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < len(QUANT_COLUMNS):
                raise ValueError(f"incomplete quant row at {quant_path}:{lineno}")
            # Columns are taken from the end, after the taxon name.
            taxid, rank = parts[-6], parts[-5]
            if rank != "species":
                continue
            if not taxid.isdigit() or int(taxid) <= 0:
                raise ValueError(f"invalid species taxid at {quant_path}:{lineno}")
            if taxid in values:
                raise ValueError(f"duplicate species {taxid} in {quant_path}")
            if profile_method == "native-sequence":
                abundance = float(parts[-1])
                length = float(parts[-4])
                if not math.isfinite(abundance) or abundance < 0:
                    raise ValueError(f"invalid abundance at {quant_path}:{lineno}")
                if not math.isfinite(length) or length < 0 or (abundance > 0 and length == 0):
                    raise ValueError(f"invalid genome size at {quant_path}:{lineno}")
                value = abundance * length
            else:
                value = float(parts[-3])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"invalid profile weight at {quant_path}:{lineno}")
            values[taxid] = value
    total = sum(values.values())
    if not math.isfinite(total):
        raise ValueError(f"nonfinite profile total in {quant_path}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fout:
        fout.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
        for taxid, value in values.items():
            if value > 0:
                fout.write(f"{taxid}\tspecies\t-\t-\t{100 * value / total:.17g}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input")
    ap.add_argument("--out")
    ap.add_argument("--quant")
    ap.add_argument("--profile-out")
    ap.add_argument("--profile-method", choices=("native-sequence", "read-counts"), default="native-sequence")
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
        convert_centrifuger_quant(quant_path=Path(args.quant), out_path=Path(args.profile_out), profile_method=args.profile_method)


if __name__ == "__main__":
    main()

