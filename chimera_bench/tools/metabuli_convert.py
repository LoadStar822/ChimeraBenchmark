from __future__ import annotations

import argparse
from pathlib import Path


def convert_metabuli_output(*, input_path: Path, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with input_path.open("r", encoding="utf-8", errors="strict") as fin, out_path.open(
        "w", encoding="utf-8"
    ) as fout:
        for raw in fin:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            status, read_id, taxid = parts[:3]
            if not read_id:
                continue
            if status != "1" or taxid in {"0", "-", "unclassified", ""}:
                fout.write(f"{read_id}\tunclassified\n")
            else:
                fout.write(f"{read_id}\t{taxid}\n")


def convert_metabuli_report(*, report_path: Path, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("r", encoding="utf-8", errors="strict") as fin, out_path.open(
        "w", encoding="utf-8"
    ) as fout:
        fout.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
        for raw in fin:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5 or parts[3].strip().lower() != "species":
                continue
            taxid = parts[4].strip()
            try:
                count = int(parts[1])
            except ValueError:
                continue
            if taxid.isdigit() and taxid != "0" and count > 0:
                fout.write(f"{taxid}\tspecies\t-\t-\t{count}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--report")
    parser.add_argument("--profile-out")
    args = parser.parse_args()
    convert_metabuli_output(input_path=Path(args.input), out_path=Path(args.out))
    if bool(args.report) != bool(args.profile_out):
        parser.error("--report and --profile-out must be provided together")
    if args.report:
        convert_metabuli_report(report_path=Path(args.report), out_path=Path(args.profile_out))


if __name__ == "__main__":
    main()
