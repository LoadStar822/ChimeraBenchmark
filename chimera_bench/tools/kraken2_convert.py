from __future__ import annotations

import argparse
from pathlib import Path


def convert_kraken2_output(*, input_path: Path, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
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
            status = parts[0].strip()
            read_id = parts[1].strip()
            taxid = parts[2].strip()
            if not read_id:
                continue
            if status == "U" or taxid in {"0", "-", "unclassified", ""}:
                fout.write(f"{read_id}\tunclassified\n")
            else:
                fout.write(f"{read_id}\t{taxid}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    convert_kraken2_output(input_path=Path(args.input), out_path=Path(args.out))


if __name__ == "__main__":
    main()

