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

    @staticmethod
    def _resolve_prefix(path: str | Path, out_dir: Path) -> Path:
        prefix = Path(path)
        if not prefix.is_absolute():
            prefix = out_dir / prefix
        return prefix

    @staticmethod
    def _resolve_database_path(path: str | Path) -> Path:
        db_path = Path(path)
        if db_path.exists():
            return db_path
        if db_path.suffix != ".imcf":
            candidate = db_path.with_suffix(".imcf")
            if candidate.exists():
                return candidate
        return db_path

    @staticmethod
    def _native_profile_tsv_path(out_prefix: str) -> str:
        classify_tsv = Path(f"{out_prefix}.tsv")
        if classify_tsv.name == "ChimeraClassify.tsv":
            return str(classify_tsv.with_name("ChimeraProfile.tsv"))
        return str(classify_tsv.with_suffix(".profile.tsv"))

    @staticmethod
    def _cami_profile_tsv_path(out_prefix: str) -> str:
        classify_tsv = Path(f"{out_prefix}.tsv")
        if classify_tsv.name == "ChimeraClassify.tsv":
            return str(classify_tsv.with_name("ChimeraProfile.cami.tsv"))
        return str(classify_tsv.with_suffix(".profile.cami.tsv"))

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
        db_path = self._resolve_database_path(db)
        threads = str(exp.get("threads", 192))
        tool_args = list(exp.get("tool_args", []))
        if "--profile-cami" not in tool_args:
            tool_args.append("--profile-cami")

        cmd = [bin_path, "classify", "-d", str(db_path), "-o", out_prefix, "-t", threads]
        if "reads" in dataset:
            cmd += ["-i", *dataset["reads"]]
        elif "paired" in dataset:
            cmd += ["-p", *dataset["paired"]]
        else:
            raise ValueError("dataset must define reads or paired")
        cmd += tool_args

        outputs = {
            "classify_tsv": f"{out_prefix}.tsv",
            "chimera_profile_tsv": self._native_profile_tsv_path(out_prefix),
            "cami_profile_tsv": self._cami_profile_tsv_path(out_prefix),
        }
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
        taxonomy_dir = build_cfg.get("taxonomy_dir")
        if taxonomy_dir:
            cmd += ["--taxonomy-dir", str(taxonomy_dir)]

        steps = [
            {"name": "mkdir_db", "cmd": mkdir_cmd},
            {
                "name": "build_db",
                "cmd": cmd,
                "outputs": {
                    "db_prefix": str(out_base),
                    "db_dir": str(out_base),
                    "db_file": str(out_base / "core.imcf"),
                },
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
