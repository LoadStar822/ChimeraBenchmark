from __future__ import annotations

from pathlib import Path
import shlex
from typing import Any, Dict, List


class TaxorTool:
    name = "taxor"
    output_basename = "taxor"

    def __init__(self, config: Dict[str, Any] | None = None) -> None:
        self.config = config or {}

    def _resolve_nodes_dmp(self, exp: Dict[str, Any]) -> str | None:
        for key in (
            "taxonomy_nodes_dmp",
            "nodes_dmp",
            "taxonomy_nodes",
            "coverage_nodes_dmp",
        ):
            value = exp.get(key)
            if value:
                return str(value)
        tax_dir = exp.get("coverage_taxonomy_dir")
        if tax_dir:
            return str(Path(tax_dir) / "nodes.dmp")
        return None

    def _resolve_names_dmp(self, exp: Dict[str, Any], nodes_dmp: str | None) -> str | None:
        for key in (
            "taxonomy_names_dmp",
            "names_dmp",
            "taxonomy_names",
            "coverage_names_dmp",
        ):
            value = exp.get(key)
            if value:
                return str(value)
        tax_dir = exp.get("coverage_taxonomy_dir")
        if tax_dir:
            return str(Path(tax_dir) / "names.dmp")
        if nodes_dmp:
            return str(Path(nodes_dmp).with_name("names.dmp"))
        return None

    def _base_cmd(self) -> List[str]:
        env = self.config.get("env", "taxor")
        bin_path = self.config.get("bin", "taxor")
        return ["conda", "run", "-n", env, bin_path]

    def build_steps(
        self,
        *,
        dataset: Dict[str, Any],
        exp: Dict[str, Any],
        out_prefix: str,
        profile_dir: str | None = None,
        profile_out_prefix: str | None = None,
    ):
        db_prefix = exp.get("db") or exp.get("db_prefix")
        if not db_prefix:
            raise ValueError("taxor requires db prefix in experiment config")

        threads_raw = exp.get("threads", 1)
        try:
            threads = int(threads_raw)
        except (TypeError, ValueError):
            threads = 1
        threads = max(1, min(threads, 32))

        tool_args = list(exp.get("tool_args", []))

        index_path = Path(db_prefix)
        if index_path.suffix != ".hixf":
            index_path = Path(f"{db_prefix}.hixf")

        search_out = f"{out_prefix}_search.tsv"

        if "reads" in dataset:
            query_inputs = list(dataset["reads"])
        elif "paired" in dataset:
            query_inputs = list(dataset["paired"])
        else:
            raise ValueError("dataset must define reads or paired")

        if not query_inputs:
            raise ValueError("taxor requires at least one query file")

        is_paired = "paired" in dataset
        if is_paired and len(query_inputs) != 2:
            raise ValueError("taxor paired dataset must provide exactly two read files")

        search_cmd: List[str]
        if len(query_inputs) == 1:
            search_cmd = self._base_cmd() + [
                "search",
                "--index-file",
                str(index_path),
                "--query-file",
                str(query_inputs[0]),
                "--output-file",
                search_out,
                "--threads",
                str(threads),
            ] + tool_args
        else:
            first = str(query_inputs[0])
            fasta_exts = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
            ext = ".fasta" if first.endswith(fasta_exts) else ".fq"
            query_path = str(Path(out_prefix).with_name(f"{Path(out_prefix).name}_query{ext}"))
            search_list = self._base_cmd() + [
                "search",
                "--index-file",
                str(index_path),
                "--query-file",
                query_path,
                "--output-file",
                search_out,
                "--threads",
                str(threads),
            ] + tool_args
            query_q = shlex.quote(query_path)

            if is_paired:
                python_script = r"""
import gzip
import sys
from pathlib import Path

path = Path(sys.argv[1])
suffix = sys.argv[2]

def open_text(p: Path):
    if str(p).endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="ignore")
    return p.open("r", encoding="utf-8", errors="ignore")

def rewrite_header(line: str, prefix: str) -> str:
    if not line.startswith(prefix):
        return line
    content = line[len(prefix):].rstrip("\n")
    if not content:
        return line
    parts = content.split(None, 1)
    seq_id = parts[0] if parts else ""
    rest = (" " + parts[1]) if len(parts) > 1 else ""
    if not seq_id:
        return line
    return f"{prefix}{seq_id}{suffix}{rest}\n"

with open_text(path) as fh:
    first = fh.readline()
    if not first:
        raise SystemExit(0)
    if first.startswith("@"):
        # FASTQ: header, seq, plus, qual
        header = first
        while header:
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not seq or not plus or not qual:
                break
            sys.stdout.write(rewrite_header(header, "@"))
            sys.stdout.write(seq)
            sys.stdout.write(rewrite_header(plus, "+"))
            sys.stdout.write(qual)
            header = fh.readline()
    elif first.startswith(">"):
        # FASTA: rewrite only header lines
        line = first
        while line:
            if line.startswith(">"):
                sys.stdout.write(rewrite_header(line, ">"))
            else:
                sys.stdout.write(line)
            line = fh.readline()
    else:
        sys.stdout.write(first)
        for line in fh:
            sys.stdout.write(line)
"""
                r1 = shlex.quote(str(query_inputs[0]))
                r2 = shlex.quote(str(query_inputs[1]))
                script = "\n".join(
                    [
                        "set -euo pipefail",
                        f"rm -f {query_q}",
                        f"python -u - {r1} /1 <<'PY' > {query_q}\n{python_script}\nPY",
                        f"python -u - {r2} /2 <<'PY' >> {query_q}\n{python_script}\nPY",
                        "cleanup() {",
                        f"  rm -f {query_q} || true",
                        "}",
                        "trap cleanup EXIT",
                        shlex.join(search_list),
                    ]
                )
            else:
                script = "\n".join(
                    [
                        "set -euo pipefail",
                        f"rm -f {query_q}",
                        f"cat {' '.join(shlex.quote(str(p)) for p in query_inputs)} > {query_q}",
                        "cleanup() {",
                        f"  rm -f {query_q} || true",
                        "}",
                        "trap cleanup EXIT",
                        shlex.join(search_list),
                    ]
                )
            search_cmd = ["bash", "-lc", script]

        cami_report_file = f"{profile_out_prefix}_genomic.tsv" if profile_out_prefix else f"{out_prefix}_genomic.tsv"
        seq_abundance_file = f"{profile_out_prefix}.tsv" if profile_out_prefix else f"{out_prefix}_abundance.tsv"
        binning_file = f"{out_prefix}_binning.tsv"
        sample_id = dataset.get("name") or "sample"

        profile_cmd = self._base_cmd() + [
            "profile",
            "--search-file",
            search_out,
            "--cami-report-file",
            cami_report_file,
            "--seq-abundance-file",
            seq_abundance_file,
            "--binning-file",
            binning_file,
            "--sample-id",
            str(sample_id),
            "--threads",
            str(threads),
        ] + tool_args

        nodes_dmp = self._resolve_nodes_dmp(exp)
        names_dmp = self._resolve_names_dmp(exp, nodes_dmp)
        fix_script = str(Path(__file__).resolve().with_name("taxor_fix_search.py"))
        fix_cmd = ["python", fix_script, "--search-file", search_out]
        if nodes_dmp:
            fix_cmd += ["--nodes-dmp", str(nodes_dmp)]
        if names_dmp:
            fix_cmd += ["--names-dmp", str(names_dmp)]

        return [
            {
                "name": "search",
                "cmd": search_cmd,
                "outputs": {"search_file": search_out},
            },
            {
                "name": "fix_search",
                "cmd": fix_cmd,
                "outputs": {"search_file": search_out},
            },
            {
                "name": "profile",
                "cmd": profile_cmd,
                "outputs": {
                    "cami_profile_tsv": seq_abundance_file,
                    "classify_tsv": binning_file,
                },
            },
        ]

    def build_db_steps(self, *, build: Dict[str, Any], out_dir: str):
        db_prefix = build.get("db_prefix") or build.get("db")
        if not db_prefix:
            raise ValueError("taxor build requires db_prefix")

        threads_raw = build.get("threads", 1)
        try:
            threads = int(threads_raw)
        except (TypeError, ValueError):
            threads = 1
        threads = max(1, min(threads, 32))

        build_cfg = build.get("build", {})
        input_sequence_dir = build_cfg.get("input_sequence_dir")
        if not input_sequence_dir:
            raise ValueError("taxor build requires build.input_sequence_dir")

        output_filename = build_cfg.get("output_filename") or f"{db_prefix}.hixf"

        steps = []

        db_dir = str(Path(db_prefix).parent)
        if db_dir and db_dir != ".":
            steps.append({"name": "mkdir_db", "cmd": ["mkdir", "-p", db_dir]})

        input_file = build_cfg.get("input_file")
        if input_file:
            input_file_path = str(input_file)
        else:
            assembly_summary = build_cfg.get("assembly_summary")
            target_tsv = build_cfg.get("target_tsv")
            if not assembly_summary or not target_tsv:
                raise ValueError("taxor build requires build.input_file or (build.assembly_summary and build.target_tsv)")
            input_file_path = str(Path(out_dir) / "outputs" / "taxor_input.tsv")
            prep_script = str(Path(__file__).resolve().with_name("taxor_prep.py"))
            steps.append(
                {
                    "name": "prepare_input",
                    "cmd": [
                        "python",
                        prep_script,
                        "--assembly-summary",
                        str(assembly_summary),
                        "--target-tsv",
                        str(target_tsv),
                        "--out",
                        input_file_path,
                    ],
                    "outputs": {"input_file": input_file_path},
                }
            )

        cmd = self._base_cmd() + [
            "build",
            "--input-file",
            input_file_path,
            "--input-sequence-dir",
            str(input_sequence_dir),
            "--output-filename",
            str(output_filename),
            "--threads",
            str(threads),
        ]

        kmer_size = build_cfg.get("kmer_size")
        if kmer_size is not None:
            cmd += ["--kmer-size", str(kmer_size)]
        syncmer_size = build_cfg.get("syncmer_size")
        if syncmer_size is not None:
            cmd += ["--syncmer-size", str(syncmer_size)]
        if build_cfg.get("use_syncmer"):
            cmd.append("--use-syncmer")

        extra_args = list(build_cfg.get("args", []))
        cmd += extra_args

        steps.append(
            {
                "name": "build_db",
                "cmd": cmd,
                "outputs": {
                    "db_prefix": db_prefix,
                    "db_file": str(output_filename),
                    "input_file": input_file_path,
                },
            }
        )
        return steps
