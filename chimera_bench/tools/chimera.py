from __future__ import annotations

from pathlib import Path
import shlex
from typing import Any, Dict, List, Tuple


class ChimeraTool:
    name = "chimera"
    output_basename = "ChimeraClassify"

    def __init__(self, config: Dict[str, Any] | None = None) -> None:
        self.config = config or {}

    def _bin(self) -> str:
        return self.config.get("bin", "Chimera")

    def _python(self) -> str:
        return self.config.get("python", "python")

    def _profile_script(self) -> str:
        script = self.config.get("profile_script") or self.config.get("chimera_py")
        if script:
            return str(script)

        bin_path = Path(self._bin())
        if bin_path.is_absolute():
            repo_root = bin_path.parent.parent
            candidate = repo_root / "chimera.py"
            if candidate.exists():
                return str(candidate)

        # fallback: rely on PATH/relative resolution
        return "chimera.py"

    @staticmethod
    def _resolve_prefix(path: str | Path, out_dir: Path) -> Path:
        prefix = Path(path)
        if not prefix.is_absolute():
            prefix = out_dir / prefix
        return prefix

    @staticmethod
    def _resolve_imcf(path: str | Path) -> Path:
        db_path = Path(path)
        if db_path.suffix != ".imcf":
            candidate = db_path.with_suffix(".imcf")
            if candidate.exists():
                return candidate
        return db_path

    def build_cmd(
        self,
        *,
        dataset: Dict[str, Any],
        exp: Dict[str, Any],
        out_prefix: str,
        profile_dir: str | None = None,
        profile_out_prefix: str | None = None,
    ) -> Tuple[List[str], Dict[str, str]]:
        bin_path = self._bin()
        db = exp["db"]
        db_path = self._resolve_imcf(db)
        threads = str(exp.get("threads", 192))
        tool_args = list(exp.get("tool_args", []))

        cmd = [bin_path, "classify", "-d", str(db_path), "-o", out_prefix, "-t", threads]
        if "reads" in dataset:
            cmd += ["-i", *dataset["reads"]]
        elif "paired" in dataset:
            cmd += ["-p", *dataset["paired"]]
        else:
            raise ValueError("dataset must define reads or paired")
        cmd += tool_args

        outputs = {"classify_tsv": f"{out_prefix}.tsv"}
        return cmd, outputs

    def build_steps(
        self,
        *,
        dataset: Dict[str, Any],
        exp: Dict[str, Any],
        out_prefix: str,
        profile_dir: str | None = None,
        profile_out_prefix: str | None = None,
    ):
        classify_cmd, outputs = self.build_cmd(
            dataset=dataset,
            exp=exp,
            out_prefix=out_prefix,
            profile_dir=profile_dir,
            profile_out_prefix=profile_out_prefix,
        )

        steps = [{"name": "classify", "cmd": classify_cmd, "outputs": outputs}]

        # Chimera's abundance profiling is implemented in the Python CLI `chimera.py profile`,
        # which consumes one or more `ChimeraClassify.tsv` files and outputs a species-level TSV.
        if profile_out_prefix:
            classify_tsv = f"{out_prefix}.tsv"
            profile_base = str(Path(profile_out_prefix).resolve())
            profile_tsv = f"{profile_base}.tsv"
            profile_cmd = [
                self._python(),
                self._profile_script(),
                "profile",
                "-i",
                classify_tsv,
                "-o",
                profile_base,
            ]
            steps.append(
                {
                    "name": "profile",
                    "cmd": profile_cmd,
                    "outputs": {"chimera_profile_tsv": profile_tsv},
                }
            )

        return steps

    def build_db_steps(self, *, build: Dict[str, Any], out_dir: str):
        db_prefix = build.get("db_prefix") or build.get("db")
        if not db_prefix:
            raise ValueError("chimera build requires db_prefix")

        threads = str(build.get("threads", 32))
        build_cfg = build.get("build", {})
        input_tsv = build_cfg.get("input_tsv") or build_cfg.get("target_tsv") or build_cfg.get("input")
        if not input_tsv:
            raise ValueError("chimera build requires build.input_tsv (or build.target_tsv)")

        out_dir_path = Path(out_dir).resolve()
        db_prefix_path = self._resolve_prefix(db_prefix, out_dir_path)
        out_base = db_prefix_path.with_suffix("") if db_prefix_path.suffix == ".imcf" else db_prefix_path

        mkdir_cmd = ["mkdir", "-p", str(out_base.parent)]

        cmd: List[str] = [
            self._bin(),
            "build",
            "-i",
            str(input_tsv),
            "-o",
            str(out_base),
            "-t",
            threads,
        ] + list(build_cfg.get("args", []))

        steps = [
            {"name": "mkdir_db", "cmd": mkdir_cmd},
            {
                "name": "build_db",
                "cmd": cmd,
                "outputs": {"db_prefix": str(out_base), "db_file": str(out_base.with_suffix(".imcf"))},
            },
        ]

        taxonomy_tax = build_cfg.get("taxonomy_tax") or build_cfg.get("taxonomy") or build_cfg.get("tax_path")
        if taxonomy_tax:
            tax_source = Path(taxonomy_tax)
            tax_target = Path(f"{out_base}.tax")
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
