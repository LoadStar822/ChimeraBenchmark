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
            "--output-one",
            "--output-all",
            "--output-unclassified",
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
        one_path = f"{out_prefix}.one"
        all_path = f"{out_prefix}.all"
        unc_path = f"{out_prefix}.unc"
        reads_prefix = f"{out_prefix}_reads"
        abundance_prefix = profile_out_prefix or f"{out_prefix}_abundance"

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
                "outputs": {
                    "rep": rep_path,
                    "classify_one": one_path,
                    "classify_all": all_path,
                    "classify_unc": unc_path,
                },
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

    def build_db_steps(self, *, build: Dict[str, Any], out_dir: str):
        db_prefix = build.get("db_prefix") or build.get("db")
        if not db_prefix:
            raise ValueError("ganon build requires db_prefix")
        threads = str(build.get("threads", 192))
        build_cfg = build.get("build", {})
        mode = build_cfg.get("mode", "custom")
        extra_args = list(build_cfg.get("args", []))

        if mode == "custom":
            cmd = self._base_cmd() + [
                "build-custom",
                "--db-prefix",
                db_prefix,
                "--threads",
                threads,
            ]
            inputs = build_cfg.get("input") or build_cfg.get("inputs")
            input_file = build_cfg.get("input_file")
            if inputs and input_file:
                raise ValueError("ganon build-custom expects input or input_file, not both")
            if inputs:
                if isinstance(inputs, str):
                    inputs = [inputs]
                cmd += ["--input", *inputs]
            elif input_file:
                cmd += ["--input-file", input_file]
            else:
                raise ValueError("ganon build-custom requires input or input_file")
            input_target = build_cfg.get("input_target")
            if input_target:
                cmd += ["--input-target", input_target]
            level = build_cfg.get("level")
            if level:
                cmd += ["--level", level]
            taxonomy_files = build_cfg.get("taxonomy_files")
            if taxonomy_files:
                if isinstance(taxonomy_files, str):
                    taxonomy_files = [taxonomy_files]
                cmd += ["--taxonomy-files", *taxonomy_files]
            cmd += extra_args
        elif mode == "build":
            cmd = self._base_cmd() + [
                "build",
                "--db-prefix",
                db_prefix,
                "--threads",
                threads,
            ]
            source = build_cfg.get("source")
            if source:
                cmd += ["--source", source]
            if build_cfg.get("complete_genomes"):
                cmd += ["--complete-genomes"]
            organism_group = build_cfg.get("organism_group")
            if organism_group:
                cmd += ["--organism-group", organism_group]
            taxid = build_cfg.get("taxid")
            if taxid:
                cmd += ["--taxid", str(taxid)]
            taxonomy = build_cfg.get("taxonomy")
            if taxonomy:
                cmd += ["--taxonomy", taxonomy]
            cmd += extra_args
        else:
            raise ValueError(f"ganon build mode not supported: {mode}")

        outputs = {"db_prefix": db_prefix}
        return [
            {
                "name": "build_db",
                "cmd": cmd,
                "outputs": outputs,
            }
        ]
