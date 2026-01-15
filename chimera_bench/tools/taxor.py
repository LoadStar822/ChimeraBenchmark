from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List


class TaxorTool:
    name = "taxor"
    output_basename = "taxor"

    def __init__(self, config: Dict[str, Any] | None = None) -> None:
        self.config = config or {}

    def _base_cmd(self) -> List[str]:
        env = self.config.get("env", "taxor")
        bin_path = self.config.get("bin", "taxor")
        return ["conda", "run", "-n", env, bin_path]

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
