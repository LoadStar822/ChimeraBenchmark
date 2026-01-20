from __future__ import annotations

from pathlib import Path
import shlex
from typing import Any, Dict, List


class CentrifugerTool:
    name = "centrifuger"
    output_basename = "centrifuger"

    def __init__(self, config: Dict[str, Any] | None = None) -> None:
        self.config = config or {}

    def _bin(self) -> str:
        return self.config.get("bin", "centrifuger")

    def _build_bin(self) -> str:
        build_bin = self.config.get("build_bin")
        if build_bin:
            return str(build_bin)
        return str(Path(self._bin()).with_name("centrifuger-build"))

    def _quant_bin(self) -> str:
        quant_bin = self.config.get("quant_bin")
        if quant_bin:
            return str(quant_bin)
        return str(Path(self._bin()).with_name("centrifuger-quant"))

    @staticmethod
    def _resolve_prefix(path: str | Path, out_dir: Path) -> Path:
        prefix = Path(path)
        if not prefix.is_absolute():
            prefix = out_dir / prefix
        return prefix

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
            raise ValueError("centrifuger requires db prefix in experiment config (exp.db)")

        threads = str(exp.get("threads", 32))
        tool_args = list(exp.get("tool_args", []))

        raw_tsv = f"{out_prefix}.tsv"
        classify_tsv = f"{out_prefix}_classify.tsv"

        classify_parts: List[str] = [self._bin(), "-x", str(db_prefix), "-t", threads]
        if "paired" in dataset:
            paired = dataset["paired"]
            if len(paired) != 2:
                raise ValueError("centrifuger paired dataset must provide exactly two read files")
            classify_parts += ["-1", paired[0], "-2", paired[1]]
        elif "reads" in dataset:
            reads = list(dataset["reads"])
            if not reads:
                raise ValueError("centrifuger reads dataset must include at least one file")
            for path in reads:
                classify_parts += ["-u", str(path)]
        else:
            raise ValueError("dataset must define reads or paired")

        classify_parts += tool_args

        classify_script = "\n".join(
            [
                "set -euo pipefail",
                f"{shlex.join(classify_parts)} > {shlex.quote(raw_tsv)}",
            ]
        )
        classify_cmd = ["bash", "-lc", classify_script]

        convert_script = str(Path(__file__).resolve().with_name("centrifuger_convert.py"))
        convert_cmd = ["python", convert_script, "--input", raw_tsv, "--out", classify_tsv]

        steps = [
            {
                "name": "classify",
                "cmd": classify_cmd,
                "outputs": {"centrifuger_tsv": raw_tsv},
            },
            {"name": "convert", "cmd": convert_cmd, "outputs": {"classify_tsv": classify_tsv}},
        ]

        # Only run quantification (profile) when profile output dir is enabled.
        if profile_out_prefix:
            cami_profile = f"{profile_out_prefix}_cami_profile.tsv"
            quant_parts = [
                self._quant_bin(),
                "-x",
                str(db_prefix),
                "-c",
                raw_tsv,
                "--output-format",
                "2",
            ]
            quant_script = "\n".join(
                [
                    "set -euo pipefail",
                    f"{shlex.join(quant_parts)} > {shlex.quote(cami_profile)}",
                ]
            )
            steps.append(
                {
                    "name": "quant",
                    "cmd": ["bash", "-lc", quant_script],
                    "outputs": {"cami_profile_tsv": cami_profile},
                }
            )

        return steps

    def build_db_steps(self, *, build: Dict[str, Any], out_dir: str):
        db_prefix = build.get("db_prefix") or build.get("db")
        if not db_prefix:
            raise ValueError("centrifuger build requires db_prefix")

        threads = str(build.get("threads", 32))
        build_cfg = build.get("build", {})
        target_tsv = build_cfg.get("target_tsv")
        if not target_tsv:
            raise ValueError("centrifuger build requires build.target_tsv")

        tax_dir = build_cfg.get("taxonomy_dir")
        nodes = build_cfg.get("taxonomy_nodes_dmp") or (
            str(Path(tax_dir) / "nodes.dmp") if tax_dir else None
        )
        names = build_cfg.get("taxonomy_names_dmp") or (
            str(Path(tax_dir) / "names.dmp") if tax_dir else None
        )
        if not (nodes and names):
            raise ValueError("centrifuger build requires taxonomy_dir or explicit taxonomy_nodes_dmp/names_dmp")

        out_dir_path = Path(out_dir).resolve()
        db_prefix_path = self._resolve_prefix(db_prefix, out_dir_path)

        mkdir_cmd = ["mkdir", "-p", str(db_prefix_path.parent)]

        cmd: List[str] = [
            self._build_bin(),
            "-l",
            str(target_tsv),
            "--taxonomy-tree",
            str(nodes),
            "--name-table",
            str(names),
            "-o",
            str(db_prefix_path),
            "-t",
            threads,
        ]

        conversion_table = build_cfg.get("conversion_table") or build_cfg.get("seqid2taxid_map")
        if conversion_table:
            cmd += ["--conversion-table", str(conversion_table)]

        build_mem = build_cfg.get("build_mem")
        if build_mem:
            cmd += ["--build-mem", str(build_mem)]

        cmd += list(build_cfg.get("args", []))

        return [
            {"name": "mkdir_db", "cmd": mkdir_cmd},
            {
                "name": "build_db",
                "cmd": cmd,
                "outputs": {"db_prefix": str(db_prefix_path), "target_tsv": str(target_tsv)},
            },
        ]

