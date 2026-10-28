from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any, Dict, List

LONG_READ_PLATFORMS = {"nanopore", "pacbio-hifi", "long-read"}
CONTIG_PLATFORMS = {"contigs", "short-read-contigs"}

# Streams the input files, in order, with bare "+" separator lines into the MetaPhlAn
# command. Arguments: the number of input files, the files, then the command.
BARE_PLUS_PIPE = r"""set -o pipefail
n=$1; shift
files=("${@:1:n}"); shift "$n"
for f in "${files[@]}"; do gzip -cdf -- "$f" || exit 1; done |
  awk 'NR % 4 == 3 && /^\+/ {print "+"; next} {print}' | "$@"
"""


class MetaPhlAnTool:
    """MetaPhlAn 4 profiling with its own marker database.

    Reads go through bowtie2, or minimap2 with `--long_reads` for long-read platforms.
    `-t rel_ab_w_read_stats` adds the estimated reads (bases for long reads) per clade,
    the sequence-abundance unit used by the main profile table. Paired files are given
    as one comma-separated input, MetaPhlAn's usual way of profiling paired reads
    (its -1/-2 options are for --subsampling_paired only).

    MetaPhlAn parses FASTQ with Biopython, which rejects a "+" line repeating a title
    other than the "@" line's (the PRJNA637878 single-read files keep the original SRA
    title there). Such inputs reach MetaPhlAn on stdin with bare "+" lines; sequences
    and qualities are unchanged.
    """

    name = "metaphlan"
    output_basename = "metaphlan"

    def __init__(self, config: Dict[str, Any] | None = None) -> None:
        self.config = config or {}

    def _base_cmd(self, *, stdin: bool = False) -> List[str]:
        bin_path = self.config.get("bin", "metaphlan")
        env_prefix = self.config.get("env_prefix")
        # conda run hands its stdin to the command only when it does not capture output.
        run = ["conda", "run", "--no-capture-output"] if stdin else ["conda", "run"]
        if env_prefix:
            # An environment unpacked outside conda's envs directory is addressed by path.
            return run + ["-p", str(env_prefix), bin_path]
        return run + ["-n", self.config.get("env", "metaphlan"), bin_path]

    @staticmethod
    def _input_type(path: str) -> str:
        name = path[:-3] if path.endswith(".gz") else path
        return "fasta" if name.endswith((".fa", ".fna", ".fasta")) else "fastq"

    @staticmethod
    def _plus_title_differs(path: str) -> bool:
        with (gzip.open if path.endswith(".gz") else open)(path, "rt") as fh:
            header, _, plus = fh.readline(), fh.readline(), fh.readline()
        title = plus[1:].rstrip()
        return header.startswith("@") and bool(title) and title != header[1:].rstrip()

    def build_steps(
        self,
        *,
        dataset: Dict[str, Any],
        exp: Dict[str, Any],
        out_prefix: str,
        profile_dir: str | None = None,
        profile_out_prefix: str | None = None,
    ):
        db_dir = exp.get("db_dir") or exp.get("db")
        index = exp.get("index")
        if not db_dir or not index:
            raise ValueError("metaphlan requires db_dir and index in experiment config")
        platform = str(dataset.get("platform", ""))
        if platform in CONTIG_PLATFORMS:
            raise ValueError(f"metaphlan profiles reads; dataset {dataset.get('name')} is {platform}")
        threads = str(exp.get("threads", 32))
        tool_args = list(exp.get("tool_args", []))

        if "reads" in dataset:
            inputs = list(dataset["reads"])
        elif "paired" in dataset:
            inputs = list(dataset["paired"])
        else:
            raise ValueError("dataset must define reads or paired")

        prefix = profile_out_prefix or out_prefix
        profile_path = f"{prefix}_profile.txt"
        mapout_path = f"{prefix}.mapout.bz2"
        input_type = self._input_type(inputs[0])
        if input_type == "fastq" and any(self._plus_title_differs(path) for path in inputs):
            cmd = ["bash", "-c", BARE_PLUS_PIPE, "bash", str(len(inputs)), *inputs]
            cmd += self._base_cmd(stdin=True)
        else:
            cmd = self._base_cmd() + [",".join(inputs)]
        cmd += [
            "--input_type",
            input_type,
            "--db_dir",
            str(db_dir),
            "--index",
            str(index),
            "--offline",
            "--nproc",
            threads,
            "-t",
            "rel_ab_w_read_stats",
            "--mapout",
            mapout_path,
            "--tmp_dir",
            str(Path(prefix).parent),
            "--sample_id",
            str(dataset.get("name", "sample")),
            "-o",
            profile_path,
        ]
        if platform in LONG_READ_PLATFORMS:
            cmd.append("--long_reads")
        cmd += tool_args

        return [
            {
                "name": "profile",
                "cmd": cmd,
                "outputs": {"metaphlan_profile": profile_path, "metaphlan_mapout": mapout_path},
            }
        ]
