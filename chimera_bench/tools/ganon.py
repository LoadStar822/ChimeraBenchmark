from __future__ import annotations

from typing import Any, Dict, List


class GanonTool:
    name = "ganon"
    output_basename = "ganon"

    def __init__(self, config: Dict[str, Any] | None = None) -> None:
        self.config = config or {}

    def _base_cmd(self) -> List[str]:
        env = self.config.get("env", "ganon")
        bin_path = self.config.get("bin", "ganon")
        return ["conda", "run", "-n", env, bin_path]

    def build_steps(self, *, dataset: Dict[str, Any], exp: Dict[str, Any], out_prefix: str):
        db_prefix = exp.get("db") or exp.get("db_prefix")
        if not db_prefix:
            raise ValueError("ganon requires db prefix in experiment config")
        threads = str(exp.get("threads", 192))
        tool_args = list(exp.get("tool_args", []))

        classify_cmd = self._base_cmd() + [
            "classify",
            "--db-prefix",
            db_prefix,
            "--output-prefix",
            out_prefix,
            "--threads",
            threads,
            "--skip-report",
        ]
        if "reads" in dataset:
            classify_cmd += ["--single-reads", *dataset["reads"]]
        elif "paired" in dataset:
            classify_cmd += ["--paired-reads", *dataset["paired"]]
        else:
            raise ValueError("dataset must define reads or paired")
        classify_cmd += tool_args

        rep_path = f"{out_prefix}.rep"
        reads_prefix = f"{out_prefix}_reads"
        abundance_prefix = f"{out_prefix}_abundance"

        report_reads_cmd = self._base_cmd() + [
            "report",
            "--input",
            rep_path,
            "--output-prefix",
            reads_prefix,
            "--db-prefix",
            db_prefix,
            "--output-format",
            "tsv",
            "--report-type",
            "reads",
        ]

        report_abundance_cmd = self._base_cmd() + [
            "report",
            "--input",
            rep_path,
            "--output-prefix",
            abundance_prefix,
            "--db-prefix",
            db_prefix,
            "--output-format",
            "tsv",
            "--report-type",
            "abundance",
        ]

        return [
            {
                "name": "classify",
                "cmd": classify_cmd,
                "outputs": {"rep": rep_path},
            },
            {
                "name": "report_reads",
                "cmd": report_reads_cmd,
                "outputs": {"report_reads_tre": f"{reads_prefix}.tre"},
            },
            {
                "name": "report_abundance",
                "cmd": report_abundance_cmd,
                "outputs": {"report_abundance_tre": f"{abundance_prefix}.tre"},
            },
        ]
