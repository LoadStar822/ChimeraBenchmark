from __future__ import annotations

from pathlib import Path
import shlex
from typing import Any, Dict, List


class Kraken2Tool:
    name = "kraken2"
    output_basename = "kraken2"

    def __init__(self, config: Dict[str, Any] | None = None) -> None:
        self.config = config or {}

    def _base_cmd(self) -> List[str]:
        env = self.config.get("env", "kraken")
        bin_path = self.config.get("bin", "kraken2")
        return ["conda", "run", "-n", env, bin_path]

    def _build_cmd(self) -> List[str]:
        env = self.config.get("env", "kraken")
        bin_path = self.config.get("build_bin", "kraken2-build")
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
        db_dir = exp.get("db") or exp.get("db_prefix")
        if not db_dir:
            raise ValueError("kraken2 requires db directory in experiment config (exp.db)")
        threads = str(exp.get("threads", 32))
        tool_args = list(exp.get("tool_args", []))

        out_file = f"{out_prefix}.out"
        report_file = f"{out_prefix}.report"
        classify_tsv = f"{out_prefix}_classify.tsv"

        classify_cmd = self._base_cmd() + [
            "--db",
            str(db_dir),
            "--threads",
            threads,
            "--output",
            out_file,
            "--report",
            report_file,
        ]

        if "reads" in dataset:
            classify_cmd += [*dataset["reads"]]
        elif "paired" in dataset:
            paired = dataset["paired"]
            if len(paired) != 2:
                raise ValueError("kraken2 paired dataset must provide exactly two read files")
            classify_cmd += ["--paired", paired[0], paired[1]]
        else:
            raise ValueError("dataset must define reads or paired")

        classify_cmd += tool_args

        convert_script = str(Path(__file__).resolve().with_name("kraken2_convert.py"))
        convert_cmd = ["python", convert_script, "--input", out_file, "--out", classify_tsv]

        return [
            {
                "name": "classify",
                "cmd": classify_cmd,
                "outputs": {"kraken2_out": out_file, "kraken2_report": report_file},
            },
            {
                "name": "convert",
                "cmd": convert_cmd,
                "outputs": {"classify_tsv": classify_tsv},
            },
        ]

    def build_db_steps(self, *, build: Dict[str, Any], out_dir: str):
        db_prefix = build.get("db_prefix") or build.get("db")
        if not db_prefix:
            raise ValueError("kraken2 build requires db_prefix")

        threads = str(build.get("threads", 32))
        build_cfg = build.get("build", {})
        target_tsv = build_cfg.get("target_tsv")
        if not target_tsv:
            raise ValueError("kraken2 build requires build.target_tsv")

        tax_dir = build_cfg.get("taxonomy_dir")
        nodes = build_cfg.get("taxonomy_nodes_dmp") or (
            str(Path(tax_dir) / "nodes.dmp") if tax_dir else None
        )
        names = build_cfg.get("taxonomy_names_dmp") or (
            str(Path(tax_dir) / "names.dmp") if tax_dir else None
        )
        merged = build_cfg.get("taxonomy_merged_dmp") or (
            str(Path(tax_dir) / "merged.dmp") if tax_dir else None
        )
        delnodes = build_cfg.get("taxonomy_delnodes_dmp") or (
            str(Path(tax_dir) / "delnodes.dmp") if tax_dir else None
        )
        if not (nodes and names and merged and delnodes):
            raise ValueError("kraken2 build requires taxonomy_dir or explicit taxonomy_*.dmp paths")

        taxonomy_tax = build_cfg.get("taxonomy_tax")

        out_dir_path = Path(out_dir).resolve()
        db_dir = Path(db_prefix)
        if not db_dir.is_absolute():
            db_dir = out_dir_path / db_dir

        prep_script = str(Path(__file__).resolve().with_name("kraken2_prep_library.py"))

        mkdir_cmd = [
            "bash",
            "-lc",
            f"mkdir -p {shlex.quote(str(db_dir / 'taxonomy'))} {shlex.quote(str(db_dir / 'library' / 'added'))}",
        ]

        copy_tax_cmd = [
            "bash",
            "-lc",
            "\n".join(
                [
                    "set -euo pipefail",
                    f"cp -f {shlex.quote(str(nodes))} {shlex.quote(str(db_dir / 'taxonomy' / 'nodes.dmp'))}",
                    f"cp -f {shlex.quote(str(names))} {shlex.quote(str(db_dir / 'taxonomy' / 'names.dmp'))}",
                    f"cp -f {shlex.quote(str(merged))} {shlex.quote(str(db_dir / 'taxonomy' / 'merged.dmp'))}",
                    f"cp -f {shlex.quote(str(delnodes))} {shlex.quote(str(db_dir / 'taxonomy' / 'delnodes.dmp'))}",
                ]
            ),
        ]

        prepare_cmd = [
            "python",
            prep_script,
            "--target-tsv",
            str(target_tsv),
            "--out-dir",
            str(db_dir / "library" / "added"),
        ]

        build_cmd = self._build_cmd() + [
            "--db",
            str(db_dir),
            "--build",
            "--threads",
            threads,
        ]

        steps = [
            {"name": "mkdir_db", "cmd": mkdir_cmd},
            {"name": "copy_taxonomy", "cmd": copy_tax_cmd},
            {
                "name": "prepare_library",
                "cmd": prepare_cmd,
                "outputs": {"target_tsv": str(target_tsv)},
            },
            {
                "name": "build_db",
                "cmd": build_cmd,
                "outputs": {"db_prefix": str(db_dir)},
            },
        ]

        if taxonomy_tax:
            tax_source = Path(taxonomy_tax)
            tax_target = Path(f"{db_dir}.tax")
            steps.append(
                {
                    "name": "copy_taxonomy_tax",
                    "cmd": [
                        "bash",
                        "-lc",
                        f"cp -f {shlex.quote(str(tax_source))} {shlex.quote(str(tax_target))}",
                    ],
                    "outputs": {"taxonomy": str(tax_target)},
                }
            )

        return steps

