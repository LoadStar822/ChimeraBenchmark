from __future__ import annotations

import sys
from pathlib import Path
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
        if self.config.get("separate_binning"):
            return self._separate_binning_steps(dataset=dataset, db_prefix=db_prefix, threads=threads, tool_args=tool_args, out_prefix=out_prefix)
        output_one = bool(self.config.get("output_one", False))
        output_all = self.config.get("output_all", False)
        output_unclassified = self.config.get("output_unclassified", True)
        skip_report = self.config.get("skip_report", True)
        multiple_matches = self.config.get("multiple_matches")

        classify_cmd = self._base_cmd() + [
            "classify",
            "--db-prefix",
            db_prefix,
            "--output-prefix",
            out_prefix,
            "--threads",
            threads,
        ]
        if output_one:
            classify_cmd.append("--output-one")
        if output_all:
            classify_cmd.append("--output-all")
        if output_unclassified:
            classify_cmd.append("--output-unclassified")
        if skip_report:
            classify_cmd.append("--skip-report")
        if multiple_matches:
            classify_cmd += ["--multiple-matches", str(multiple_matches)]
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
                    **({"classify_one": one_path} if output_one else {}),
                    **({"classify_all": all_path} if output_all else {}),
                    **({"classify_unc": unc_path} if output_unclassified else {}),
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

    def _read_args(self, dataset: Dict[str, Any]) -> List[str]:
        if "reads" in dataset:
            return ["--single-reads", *dataset["reads"]]
        if "paired" in dataset:
            return ["--paired-reads", *dataset["paired"]]
        raise ValueError("dataset must define reads or paired")

    def _report_cmd(self, *, db_prefix: str, rep_path: str, output_prefix: str, report_type: str) -> List[str]:
        return self._base_cmd() + [
            "report",
            "--input",
            rep_path,
            "--output-prefix",
            output_prefix,
            "--db-prefix",
            db_prefix,
            "--output-format",
            "tsv",
            "--report-type",
            report_type,
        ]

    def _separate_binning_steps(self, *, dataset: Dict[str, Any], db_prefix: str, threads: str, tool_args: List[str], out_prefix: str):
        """Per-read assignments from ganon's binning preset, the profile from a default run.

        ganon documents `--binning` for read binning and its defaults for profiling, so the
        two tasks use two classifications of the same reads. The profile is reported as
        sequence abundance (`--report-type reads`), the unit shared by all tools in the
        main profile table; the default taxonomic-abundance report is kept alongside.
        Each step runs through ganon_run.py, so a classification that assigns no reads
        yields empty outputs instead of ganon's non-zero exit.
        """
        read_args = self._read_args(dataset)
        binning_prefix = f"{out_prefix}_binning"
        binning_cmd = self._base_cmd() + [
            "classify",
            "--db-prefix",
            db_prefix,
            "--output-prefix",
            binning_prefix,
            "--threads",
            threads,
            "--binning",
            "--output-one",
            "--output-unclassified",
            "--skip-report",
            *read_args,
            *tool_args,
        ]
        profile_cmd = self._base_cmd() + [
            "classify",
            "--db-prefix",
            db_prefix,
            "--output-prefix",
            out_prefix,
            "--threads",
            threads,
            *read_args,
            *tool_args,
        ]
        rep_path = f"{out_prefix}.rep"
        binning_rep = f"{binning_prefix}.rep"
        reads_prefix = f"{out_prefix}_reads"
        report_cmd = self._report_cmd(db_prefix=db_prefix, rep_path=rep_path, output_prefix=reads_prefix, report_type="reads")
        return [
            {
                "name": "classify_binning",
                "cmd": self._tolerate_no_classification(binning_cmd, rep=binning_rep, empty_outputs=[f"{binning_prefix}.one"]),
                "outputs": {
                    "classify_one": f"{binning_prefix}.one",
                    "classify_unc": f"{binning_prefix}.unc",
                    "binning_rep": binning_rep,
                },
            },
            {
                "name": "classify_profile",
                "cmd": self._tolerate_no_classification(profile_cmd, rep=rep_path, empty_outputs=[f"{out_prefix}.tre"]),
                "outputs": {"rep": rep_path, "report_abundance_tre": f"{out_prefix}.tre"},
            },
            {
                "name": "report_reads",
                "cmd": self._tolerate_no_classification(report_cmd, rep=rep_path, empty_outputs=[f"{reads_prefix}.tre"]),
                "outputs": {"report_reads_tre": f"{reads_prefix}.tre"},
            },
        ]

    @staticmethod
    def _tolerate_no_classification(cmd: List[str], *, rep: str, empty_outputs: List[str]) -> List[str]:
        wrapper = [sys.executable, str(Path(__file__).resolve().with_name("ganon_run.py")), "--rep", rep]
        for out in empty_outputs:
            wrapper += ["--empty-output", out]
        return wrapper + ["--", *cmd]

    def build_db_steps(self, *, build: Dict[str, Any], out_dir: str):
        db_prefix = build.get("db_prefix") or build.get("db")
        if not db_prefix:
            raise ValueError("ganon build requires db_prefix")
        threads = str(build.get("threads", 192))
        build_cfg = build.get("build", {})
        mode = build_cfg.get("mode", "custom")
        extra_args = list(build_cfg.get("args", []))

        steps = []
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
            target_tsv = build_cfg.get("target_tsv")
            if target_tsv and not inputs and not input_file:
                taxonomy_files = build_cfg.get("taxonomy_files") or []
                if isinstance(taxonomy_files, str):
                    taxonomy_files = [taxonomy_files]
                nodes = next((f for f in taxonomy_files if Path(f).name == "nodes.dmp"), None)
                if not nodes:
                    raise ValueError("ganon build from target_tsv requires nodes.dmp in taxonomy_files")
                input_file = str(Path(out_dir) / "outputs" / "ganon_input.tsv")
                steps.append(
                    {
                        "name": "prepare_input",
                        "cmd": [
                            sys.executable,
                            str(Path(__file__).resolve().with_name("ganon_prep.py")),
                            "--target-tsv",
                            str(target_tsv),
                            "--nodes",
                            str(nodes),
                            "--out",
                            input_file,
                        ],
                        "outputs": {"input_file": input_file},
                    }
                )
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
        steps.append(
            {
                "name": "build_db",
                "cmd": cmd,
                "outputs": outputs,
            }
        )
        return steps
