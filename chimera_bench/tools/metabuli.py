from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List


class MetabuliTool:
    name = "metabuli"
    output_basename = "metabuli"

    def __init__(self, config: Dict[str, Any] | None = None) -> None:
        self.config = config or {}

    def _base_cmd(self) -> List[str]:
        env = self.config.get("env", "metabuli")
        bin_path = self.config.get("bin", "metabuli")
        return ["conda", "run", "-n", env, bin_path]

    @staticmethod
    def _sequence_mode(dataset: Dict[str, Any]) -> str:
        platform = str(dataset.get("platform") or dataset.get("input_type") or "").lower()
        long_markers = ("contig", "long", "nanopore", "pacbio", "hifi")
        return "3" if any(marker in platform for marker in long_markers) else "1"

    @staticmethod
    def _query_suffix(path: str) -> str:
        name = Path(path).name.lower()
        for suffix in (".fastq.gz", ".fq.gz", ".fasta.gz", ".fna.gz", ".fa.gz", ".fastq", ".fq", ".fasta", ".fna", ".fa"):
            if name.endswith(suffix):
                return suffix
        raise ValueError(f"unsupported Metabuli query file: {path}")

    def build_steps(
        self,
        *,
        dataset: Dict[str, Any],
        exp: Dict[str, Any],
        out_prefix: str,
        profile_dir: str | None = None,
        profile_out_prefix: str | None = None,
    ):
        del profile_dir
        db_dir = exp.get("db") or exp.get("db_prefix")
        if not db_dir:
            raise ValueError("metabuli requires db directory in experiment config")

        threads = str(exp.get("threads", 32))
        tool_args = list(exp.get("tool_args", []))
        out_prefix_path = Path(out_prefix).resolve()
        output_dir = out_prefix_path.parent
        job_id = out_prefix_path.name
        raw_tsv = output_dir / f"{job_id}_classifications.tsv"
        report_tsv = output_dir / f"{job_id}_report.tsv"
        classify_tsv = Path(f"{out_prefix}_classify.tsv").resolve()
        steps = []

        if "paired" in dataset:
            paired = list(dataset["paired"])
            if len(paired) != 2:
                raise ValueError("metabuli paired dataset must provide exactly two read files")
            query_args = paired
            mode_args: List[str] = []
        elif "reads" in dataset:
            reads = list(dataset["reads"])
            if not reads:
                raise ValueError("metabuli reads dataset must provide at least one file")
            query_path = reads[0]
            if len(reads) > 1:
                cache_root = Path(self.config.get("query_cache_dir", "builds/datasets/metabuli")).resolve()
                query_path = str(cache_root / str(dataset.get("name", "dataset")) / f"queries{self._query_suffix(reads[0])}")
                prep_script = str(Path(__file__).resolve().with_name("metabuli_prep.py"))
                steps.append(
                    {
                        "name": "prepare_query",
                        "cmd": ["python", prep_script, "query", "--out", query_path, "--inputs", *reads],
                        "outputs": {"prepared_query": query_path},
                    }
                )
            query_args = [query_path]
            mode_args = ["--seq-mode", self._sequence_mode(dataset)]
        else:
            raise ValueError("dataset must define reads or paired")

        classify_cmd = self._base_cmd() + [
            "classify",
            *mode_args,
            *query_args,
            str(db_dir),
            str(output_dir),
            job_id,
            "--threads",
            threads,
            *tool_args,
        ]
        steps.append(
            {
                "name": "classify",
                "cmd": classify_cmd,
                "outputs": {
                    "metabuli_classifications_tsv": str(raw_tsv),
                    "metabuli_report_tsv": str(report_tsv),
                },
            }
        )

        convert_script = str(Path(__file__).resolve().with_name("metabuli_convert.py"))
        convert_cmd = ["python", convert_script, "--input", str(raw_tsv), "--out", str(classify_tsv)]
        convert_outputs = {"classify_tsv": str(classify_tsv)}
        if profile_out_prefix:
            profile_tsv = Path(f"{profile_out_prefix}_cami_profile.tsv").resolve()
            convert_cmd += ["--report", str(report_tsv), "--profile-out", str(profile_tsv)]
            convert_outputs["cami_profile_tsv"] = str(profile_tsv)
        steps.append({"name": "convert", "cmd": convert_cmd, "outputs": convert_outputs})
        return steps

    def build_db_steps(self, *, build: Dict[str, Any], out_dir: str):
        db_prefix = build.get("db_prefix") or build.get("db")
        if not db_prefix:
            raise ValueError("metabuli build requires db_prefix")

        build_cfg = build.get("build", {})
        target_tsv = build_cfg.get("target_tsv")
        taxonomy_dir = build_cfg.get("taxonomy_dir")
        if not target_tsv:
            raise ValueError("metabuli build requires build.target_tsv")
        if not taxonomy_dir:
            raise ValueError("metabuli build requires build.taxonomy_dir")

        out_dir_path = Path(out_dir).resolve()
        db_dir = Path(db_prefix)
        if not db_dir.is_absolute():
            db_dir = out_dir_path / db_dir
        inputs_dir = db_dir / "build_inputs"
        fasta_list = inputs_dir / "fasta.list"
        assembly_map = inputs_dir / "assembly_taxid.map"
        prep_script = str(Path(__file__).resolve().with_name("metabuli_prep.py"))

        prep_cmd = [
            "python",
            prep_script,
            "db",
            "--target-tsv",
            str(target_tsv),
            "--fasta-list",
            str(fasta_list),
            "--assembly-map",
            str(assembly_map),
            "--taxonomy-dir",
            str(taxonomy_dir),
            "--db-dir",
            str(db_dir),
        ]
        build_cmd = self._base_cmd() + [
            "build",
            "--gtdb",
            "1",
            str(db_dir),
            str(fasta_list),
            str(assembly_map),
            "--threads",
            str(build.get("threads", 32)),
            *list(build_cfg.get("args", [])),
        ]
        cleanup_cmd = ["python", prep_script, "cleanup", "--db-dir", str(db_dir)]
        return [
            {
                "name": "prepare_db",
                "cmd": prep_cmd,
                "outputs": {"target_tsv": str(target_tsv), "fasta_list": str(fasta_list)},
            },
            {
                "name": "build_db",
                "cmd": build_cmd,
                "outputs": {"db_prefix": str(db_dir)},
            },
            {
                "name": "cleanup_db",
                "cmd": cleanup_cmd,
            },
        ]
