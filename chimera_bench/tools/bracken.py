from __future__ import annotations

from pathlib import Path
import shlex
from typing import Any, Dict, List


class BrackenTool:
    name = "bracken"
    output_basename = "bracken"

    def __init__(self, config: Dict[str, Any] | None = None) -> None:
        self.config = config or {}

    def _env(self) -> str:
        return self.config.get("env", "kraken")

    def _base_cmd(self) -> List[str]:
        return ["conda", "run", "-n", self._env(), self.config.get("bin", "bracken")]

    def _kraken2_cmd(self) -> List[str]:
        kraken2_bin = self.config.get("kraken2_bin", "kraken2")
        kraken2_env = self.config.get("kraken2_env", self._env())
        return ["conda", "run", "-n", kraken2_env, kraken2_bin]

    @staticmethod
    def _resolve_db_dir(db_prefix: str, out_dir: str) -> Path:
        db_path = Path(db_prefix)
        if db_path.is_absolute():
            return db_path.resolve()
        return (Path(out_dir).resolve() / db_path).resolve()

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

        out_prefix_path = Path(out_prefix).resolve()
        if profile_out_prefix is not None:
            out_prefix_path = Path(profile_out_prefix).resolve()

        run_prefix_path = Path(out_prefix).resolve()
        kraken2_out = str(run_prefix_path.with_name(f"{run_prefix_path.name}_kraken2.out"))
        kraken2_report = str(run_prefix_path.with_name(f"{run_prefix_path.name}_kraken2.report"))
        bracken_out = str(out_prefix_path) + ".tsv"
        cami_profile = str(out_prefix_path) + "_cami_profile.tsv"
        threads = str(exp.get("threads", 32))

        classify_cmd = self._kraken2_cmd() + [
            "--db",
            str(db_dir),
            "--threads",
            threads,
            "--output",
            kraken2_out,
            "--report",
            kraken2_report,
        ]
        if "reads" in dataset:
            classify_cmd += [str(path) for path in dataset["reads"]]
        elif "paired" in dataset:
            paired = list(dataset["paired"])
            if len(paired) != 2:
                raise ValueError("bracken paired dataset must provide exactly two read files")
            classify_cmd += ["--paired", str(paired[0]), str(paired[1])]
        else:
            raise ValueError("dataset must define reads or paired")

        # Use Bracken defaults (fairness): read_len=100, level=S, threshold=10.
        bracken_cmd = self._base_cmd() + [
            "-d",
            str(db_dir),
            "-i",
            kraken2_report,
            "-o",
            bracken_out,
        ]

        convert_script = str(Path(__file__).resolve().with_name("bracken_to_cami.py"))
        convert_cmd = ["python", convert_script, "--input", bracken_out, "--out", cami_profile]

        return [
            {
                "name": "kraken2_classify",
                "cmd": classify_cmd,
                "outputs": {"kraken2_out": kraken2_out, "kraken2_report": kraken2_report},
            },
            {"name": "bracken", "cmd": bracken_cmd, "outputs": {"bracken_tsv": bracken_out}},
            {"name": "convert", "cmd": convert_cmd, "outputs": {"cami_profile_tsv": cami_profile}},
        ]

    def build_db_steps(self, *, build: Dict[str, Any], out_dir: str):
        """
        Build Bracken kmer distribution files in a dedicated Bracken DB directory.

        Bracken requires `database100mers.kmer_distrib` inside the DB directory.
        If `source_db_prefix` is provided, required Kraken2 index artifacts are linked
        from source DB into Bracken DB to avoid duplicating large files.
        We intentionally keep Bracken defaults (k=35, read_len=100) for fairness.
        """

        db_prefix = build.get("db_prefix") or build.get("db")
        if not db_prefix:
            raise ValueError("bracken build requires db_prefix")

        threads = str(build.get("threads", 32))
        build_cfg = build.get("build", {})
        target_tsv = build_cfg.get("target_tsv")
        source_db_prefix = build_cfg.get("source_db_prefix") or build.get("source_db_prefix")
        cleanup_library = bool(build_cfg.get("cleanup_library_fna", True))
        if not target_tsv:
            raise ValueError("bracken build requires build.target_tsv (CAMIRefseq target.tsv)")

        db_dir = self._resolve_db_dir(str(db_prefix), out_dir)
        source_db_dir: Path | None = None
        if source_db_prefix:
            source_db_dir = self._resolve_db_dir(str(source_db_prefix), out_dir)
            if not source_db_dir.exists():
                raise FileNotFoundError(source_db_dir)

        distrib_path = db_dir / "database100mers.kmer_distrib"
        database_kraken = db_dir / "database.kraken"
        database_kraken_tmp = db_dir / "database.kraken.tmp"
        database_100mers = db_dir / "database100mers.kraken"

        out_added = db_dir / "library" / "added"
        library_fna = out_added / "library.fna"
        prep_script = str(Path(__file__).resolve().with_name("kraken2_prep_library.py"))

        steps = []
        if source_db_dir and source_db_dir != db_dir:
            q_db_dir = shlex.quote(str(db_dir))
            q_src_dir = shlex.quote(str(source_db_dir))
            required_files = ["hash.k2d", "opts.k2d", "taxo.k2d", "seqid2taxid.map"]
            stage_lines = [
                "set -euo pipefail",
                f"mkdir -p {q_db_dir}",
            ]
            for filename in required_files:
                src = source_db_dir / filename
                dst = db_dir / filename
                q_src = shlex.quote(str(src))
                q_dst = shlex.quote(str(dst))
                stage_lines.extend(
                    [
                        f"if [ ! -e {q_dst} ]; then",
                        f"  if [ ! -e {q_src} ]; then echo '[bracken] missing source file: {src}' >&2; exit 2; fi",
                        f"  ln {q_src} {q_dst}",
                        "fi",
                    ]
                )
            src_taxonomy = source_db_dir / "taxonomy"
            dst_taxonomy = db_dir / "taxonomy"
            q_src_tax = shlex.quote(str(src_taxonomy))
            q_dst_tax = shlex.quote(str(dst_taxonomy))
            stage_lines.extend(
                [
                    f"if [ ! -e {q_dst_tax} ]; then",
                    f"  if [ ! -d {q_src_tax} ]; then echo '[bracken] missing source taxonomy dir: {src_taxonomy}' >&2; exit 2; fi",
                    f"  ln -s {q_src_tax} {q_dst_tax}",
                    "fi",
                    # Reuse existing generated files when they already exist in source DB.
                    f"if [ ! -e {shlex.quote(str(database_kraken))} ] && [ -s {shlex.quote(str(source_db_dir / 'database.kraken'))} ]; then ln {shlex.quote(str(source_db_dir / 'database.kraken'))} {shlex.quote(str(database_kraken))}; fi",
                    f"if [ ! -e {shlex.quote(str(distrib_path))} ] && [ -s {shlex.quote(str(source_db_dir / 'database100mers.kmer_distrib'))} ]; then ln {shlex.quote(str(source_db_dir / 'database100mers.kmer_distrib'))} {shlex.quote(str(distrib_path))}; fi",
                ]
            )
            stage_cmd = ["bash", "-lc", "\n".join(stage_lines)]
            steps.append({"name": "stage_kraken2_db", "cmd": stage_cmd})
        else:
            ensure_db_cmd = ["bash", "-lc", f"set -euo pipefail\nmkdir -p {shlex.quote(str(db_dir))}"]
            steps.append({"name": "ensure_db_dir", "cmd": ensure_db_cmd})

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

        steps.extend(
            [
                {"name": "prepare_library", "cmd": prepare_cmd},
                {"name": "build_database", "cmd": build_database_cmd},
                {"name": "kmer2read_distr", "cmd": kmer2read_cmd},
                {"name": "generate_distrib", "cmd": distrib_cmd},
            ]
        )

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
