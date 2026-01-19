from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List


class BrackenTool:
    name = "bracken"
    output_basename = "bracken"

    def __init__(self, config: Dict[str, Any] | None = None) -> None:
        self.config = config or {}

    @staticmethod
    def _repo_root() -> Path:
        return Path(__file__).resolve().parents[2]

    def _env(self) -> str:
        return self.config.get("env", "kraken")

    def _base_cmd(self) -> List[str]:
        return ["conda", "run", "-n", self._env(), self.config.get("bin", "bracken")]

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
            raise ValueError("bracken requires db directory in experiment config (exp.db)")

        dataset_name = dataset.get("name") or ""
        if not dataset_name:
            raise ValueError("dataset missing name")

        kraken2_runs = exp.get("kraken2_runs") or "results/classify/kraken2"
        report_name = exp.get("kraken2_report_name") or "kraken2.report"
        report_path = Path(kraken2_runs) / dataset_name / "outputs" / report_name
        if not report_path.is_absolute():
            report_path = (self._repo_root() / report_path).resolve()
        if not report_path.exists():
            raise FileNotFoundError(report_path)

        out_prefix_path = Path(out_prefix).resolve()
        if profile_out_prefix is not None:
            out_prefix_path = Path(profile_out_prefix).resolve()

        bracken_out = str(out_prefix_path) + ".tsv"
        cami_profile = str(out_prefix_path) + "_cami_profile.tsv"

        # Use Bracken defaults (fairness): read_len=100, level=S, threshold=10.
        bracken_cmd = self._base_cmd() + [
            "-d",
            str(db_dir),
            "-i",
            str(report_path),
            "-o",
            bracken_out,
        ]

        convert_script = str(Path(__file__).resolve().with_name("bracken_to_cami.py"))
        convert_cmd = ["python", convert_script, "--input", bracken_out, "--out", cami_profile]

        return [
            {"name": "bracken", "cmd": bracken_cmd, "outputs": {"bracken_tsv": bracken_out}},
            {"name": "convert", "cmd": convert_cmd, "outputs": {"cami_profile_tsv": cami_profile}},
        ]

    def build_db_steps(self, *, build: Dict[str, Any], out_dir: str):
        """
        Build Bracken kmer distribution files inside an existing Kraken2 DB.

        Bracken requires `database100mers.kmer_distrib` inside the Kraken2 DB directory.
        We intentionally keep Bracken defaults (k=35, read_len=100) for fairness.
        """

        db_prefix = build.get("db_prefix") or build.get("db")
        if not db_prefix:
            raise ValueError("bracken build requires db_prefix (kraken2 DB directory)")

        threads = str(build.get("threads", 32))
        build_cfg = build.get("build", {})
        target_tsv = build_cfg.get("target_tsv")
        cleanup_library = bool(build_cfg.get("cleanup_library_fna", True))
        if not target_tsv:
            raise ValueError("bracken build requires build.target_tsv (CAMIRefseq target.tsv)")

        db_dir = Path(db_prefix).resolve()
        if not db_dir.exists():
            raise FileNotFoundError(db_dir)

        distrib_path = db_dir / "database100mers.kmer_distrib"
        database_kraken = db_dir / "database.kraken"
        database_kraken_tmp = db_dir / "database.kraken.tmp"
        database_100mers = db_dir / "database100mers.kraken"

        out_added = db_dir / "library" / "added"
        library_fna = out_added / "library.fna"
        prep_script = str(Path(__file__).resolve().with_name("kraken2_prep_library.py"))

        prepare_cmd = [
            "bash",
            "-lc",
            "\n".join(
                [
                    "set -euo pipefail",
                    f"if [ -s {distrib_path!s} ]; then echo '[bracken] kmer_distrib exists, skip prepare'; exit 0; fi",
                    f"if [ -s {library_fna!s} ]; then echo '[bracken] library.fna exists, skip prepare'; exit 0; fi",
                    f"python {prep_script} --target-tsv {target_tsv} --out-dir {out_added!s} --force",
                ]
            ),
        ]

        build_database_cmd = [
            "bash",
            "-lc",
            "\n".join(
                [
                    "set -euo pipefail",
                    f"if [ -s {distrib_path!s} ]; then echo '[bracken] kmer_distrib exists, skip database.kraken'; exit 0; fi",
                    f"if [ -s {database_kraken!s} ]; then echo '[bracken] database.kraken exists, skip'; exit 0; fi",
                    f"if [ ! -s {library_fna!s} ]; then echo '[bracken] missing library.fna, run prepare step first' >&2; exit 2; fi",
                    f"rm -f {database_kraken_tmp!s}",
                    " ".join(
                        [
                            "conda",
                            "run",
                            "-n",
                            self._env(),
                            "kraken2",
                            "--db",
                            str(db_dir),
                            "--threads",
                            threads,
                            str(library_fna),
                            ">",
                            str(database_kraken_tmp),
                        ]
                    ),
                    f"mv -f {database_kraken_tmp!s} {database_kraken!s}",
                ]
            ),
        ]

        kmer2read_cmd = [
            "bash",
            "-lc",
            "\n".join(
                [
                    "set -euo pipefail",
                    f"if [ -s {distrib_path!s} ]; then echo '[bracken] kmer_distrib exists, skip kmer2read_distr'; exit 0; fi",
                    f"if [ ! -s {database_kraken!s} ]; then echo '[bracken] missing database.kraken, run build_database step first' >&2; exit 2; fi",
                    " ".join(
                        [
                            "conda",
                            "run",
                            "-n",
                            self._env(),
                            "kmer2read_distr",
                            "--seqid2taxid",
                            str(db_dir / "seqid2taxid.map"),
                            "--taxonomy",
                            str(db_dir / "taxonomy"),
                            "--kraken",
                            str(database_kraken),
                            "--output",
                            str(database_100mers),
                            "-k",
                            "35",
                            "-l",
                            "100",
                            "-t",
                            threads,
                        ]
                    ),
                ]
            ),
        ]

        distrib_cmd = [
            "bash",
            "-lc",
            "\n".join(
                [
                    "set -euo pipefail",
                    f"if [ -s {distrib_path!s} ]; then echo '[bracken] kmer_distrib exists, skip generate'; exit 0; fi",
                    f"if [ ! -s {database_100mers!s} ]; then echo '[bracken] missing database100mers.kraken, run kmer2read_distr step first' >&2; exit 2; fi",
                    " ".join(
                        [
                            "conda",
                            "run",
                            "-n",
                            self._env(),
                            "generate_kmer_distribution.py",
                            "-i",
                            str(database_100mers),
                            "-o",
                            str(distrib_path),
                        ]
                    ),
                ]
            ),
        ]

        steps = [
            {"name": "prepare_library", "cmd": prepare_cmd},
            {"name": "build_database", "cmd": build_database_cmd},
            {"name": "kmer2read_distr", "cmd": kmer2read_cmd},
            {"name": "generate_distrib", "cmd": distrib_cmd},
        ]

        if cleanup_library:
            cleanup_cmd = [
                "bash",
                "-lc",
                "\n".join(
                    [
                        "set -euo pipefail",
                        f"if [ -s {distrib_path!s} ] && [ -f {library_fna!s} ]; then rm -f {library_fna!s}; fi",
                        f"if [ -s {distrib_path!s} ] && [ -f {database_100mers!s} ]; then rm -f {database_100mers!s}; fi",
                    ]
                ),
            ]
            steps.append({"name": "cleanup_library", "cmd": cleanup_cmd})

        return steps

