from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

# ganon exits non-zero when a classification assigns no reads: the EM reassignment
# finds an empty .all file and `ganon report` has nothing to report. The .rep file
# still records the read counts, so such a run yields empty outputs (no taxa)
# instead of a failed run. Any other failure is passed through unchanged.


def classified_reads(rep_path: Path) -> int | None:
    if not rep_path.exists():
        return None
    with rep_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if parts[0] == "#total_classified" and len(parts) > 1:
                try:
                    return int(parts[1])
                except ValueError:
                    return None
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rep", required=True, help=".rep file written by the classification")
    ap.add_argument("--empty-output", action="append", default=[], help="output to create empty when no read is classified")
    ap.add_argument("cmd", nargs=argparse.REMAINDER, help="ganon command, after --")
    args = ap.parse_args()
    cmd = args.cmd[1:] if args.cmd[:1] == ["--"] else args.cmd
    if not cmd:
        ap.error("missing ganon command")

    rc = subprocess.call(cmd)
    if rc != 0 and classified_reads(Path(args.rep)) == 0:
        for out in args.empty_output:
            path = Path(out)
            if not path.exists():
                path.write_text("", encoding="utf-8")
        print(f"ganon classified no reads (see {args.rep}); empty outputs: {' '.join(args.empty_output)}", file=sys.stderr)
        return 0
    return rc


if __name__ == "__main__":
    sys.exit(main())
