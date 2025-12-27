from __future__ import annotations

from pathlib import Path
import shlex
from typing import Any, Dict, List


class SylphTool:
    name = "sylph"
    output_basename = "sylph"

    def __init__(self, config: Dict[str, Any] | None = None) -> None:
        self.config = config or {}

    def _base_cmd(self) -> List[str]:
        env = self.config.get("env", "sylph")
        bin_path = self.config.get("bin", "sylph")
        return ["conda", "run", "-n", env, bin_path]

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
            raise ValueError("sylph requires db prefix in experiment config")
        threads = str(exp.get("threads", 192))
        tool_args = list(exp.get("tool_args", []))

        db_path = Path(db_prefix)
        if db_path.suffix != ".syldb":
            db_path = Path(f"{db_prefix}.syldb")

        if profile_out_prefix:
            out_tsv = f"{profile_out_prefix}.tsv"
        else:
            out_tsv = f"{out_prefix}_profile.tsv"

        profile_cmd = self._base_cmd() + [
            "profile",
            "--output-file",
            out_tsv,
            "-t",
            threads,
            str(db_path),
        ]

        if "reads" in dataset:
            profile_cmd += ["-r", *dataset["reads"]]
        elif "paired" in dataset:
            paired = dataset["paired"]
            if len(paired) < 2:
                raise ValueError("sylph paired dataset must provide two read files")
            profile_cmd += ["-1", paired[0], "-2", paired[1]]
        else:
            raise ValueError("dataset must define reads or paired")

        profile_cmd += tool_args

        return [
            {
                "name": "profile",
                "cmd": profile_cmd,
                "outputs": {"profile_tsv": out_tsv},
            }
        ]

    def build_db_steps(self, *, build: Dict[str, Any], out_dir: str):
        db_prefix = build.get("db_prefix") or build.get("db")
        if not db_prefix:
            raise ValueError("sylph build requires db_prefix")
        threads = str(build.get("threads", 192))
        build_cfg = build.get("build", {})
        mode = build_cfg.get("mode", "sketch")
        if mode != "sketch":
            raise ValueError(f"sylph build mode not supported: {mode}")

        out_dir_path = Path(out_dir).resolve()
        db_prefix_path = self._resolve_prefix(db_prefix, out_dir_path)
        db_dir = db_prefix_path.parent
        list_path = out_dir_path / "genomes.list"

        genomes_list = build_cfg.get("genomes_list")
        genomes_dir = build_cfg.get("genomes_dir")
        target_tsv = build_cfg.get("target_tsv")
        target_col = int(build_cfg.get("target_col", 0))
        allow_missing = bool(build_cfg.get("allow_missing", False))
        genome_glob = build_cfg.get("genome_glob")
        extra_args = list(build_cfg.get("args", []))
        taxonomy_source = build_cfg.get("taxonomy_source")

        sources = [bool(genomes_list), bool(genomes_dir), bool(target_tsv)]
        if sum(sources) > 1:
            raise ValueError("sylph build expects exactly one of genomes_list, genomes_dir, target_tsv")

        if genomes_list:
            list_path = Path(genomes_list)
            prepare_cmd = [
                "bash",
                "-lc",
                f"mkdir -p {shlex.quote(str(db_dir))}",
            ]
        elif target_tsv:
            target_tsv = Path(target_tsv)
            script = f"""
from pathlib import Path
target_tsv = Path({str(target_tsv)!r})
list_path = Path({str(list_path)!r})
col = {target_col}
allow_missing = {allow_missing}
paths = []
missing = []
with target_tsv.open('r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split('\\t')
        if len(parts) <= col:
            raise SystemExit(f'Invalid line in {{target_tsv}}: expected >= {{col+1}} columns')
        p = Path(parts[col])
        paths.append(p)
        if not p.exists():
            missing.append(p)
if missing and not allow_missing:
    raise SystemExit(f'Missing {{len(missing)}} genome files from target list (example: {{missing[0]}})')
list_path.parent.mkdir(parents=True, exist_ok=True)
with list_path.open('w') as out:
    wrote = 0
    for p in paths:
        if p.exists():
            out.write(str(p) + '\\n')
            wrote += 1
        elif allow_missing:
            continue
print(f'target entries: {{len(paths)}}; missing: {{len(missing)}}; wrote: {{wrote}}')
"""
            prepare_cmd = [
                "bash",
                "-lc",
                f"mkdir -p {shlex.quote(str(db_dir))} && python - <<'PY'\n{script}\nPY",
            ]
        elif genomes_dir:
            genomes_dir = Path(genomes_dir)
            if not genome_glob:
                patterns = ["*.fna.gz", "*.fa.gz", "*.fasta.gz", "*.fna", "*.fa", "*.fasta"]
            else:
                patterns = [genome_glob] if isinstance(genome_glob, str) else list(genome_glob)
            script = f"""
from pathlib import Path
genomes_dir = Path({genomes_dir!r})
list_path = Path({str(list_path)!r})
patterns = {patterns!r}
files = []
for pat in patterns:
    files.extend(sorted(genomes_dir.glob(pat)))
if genomes_dir.is_file():
    files = [genomes_dir]
if not files:
    raise SystemExit(f'No genome files found in {{genomes_dir}} with patterns {{patterns}}')
list_path.parent.mkdir(parents=True, exist_ok=True)
with list_path.open('w') as f:
    for p in files:
        f.write(str(p) + '\\n')
"""
            prepare_cmd = [
                "bash",
                "-lc",
                f"mkdir -p {shlex.quote(str(db_dir))} && python - <<'PY'\n{script}\nPY",
            ]
        else:
            raise ValueError("sylph build requires genomes_list, genomes_dir, or target_tsv")

        sketch_cmd = self._base_cmd() + [
            "sketch",
            "--gl",
            str(list_path),
            "-o",
            str(db_prefix_path),
            "-t",
            threads,
        ] + extra_args

        steps = [
            {
                "name": "prepare_genomes",
                "cmd": prepare_cmd,
                "outputs": {"genomes_list": str(list_path)},
            },
            {
                "name": "build_db",
                "cmd": sketch_cmd,
                "outputs": {
                    "db_prefix": str(db_prefix_path),
                    "syldb": f"{db_prefix_path}.syldb",
                },
            },
        ]

        if taxonomy_source:
            taxonomy_source = Path(taxonomy_source)
            taxonomy_target = Path(f"{db_prefix_path}.tax")
            steps.append(
                {
                    "name": "copy_taxonomy",
                    "cmd": [
                        "bash",
                        "-lc",
                        f"cp -f {shlex.quote(str(taxonomy_source))} {shlex.quote(str(taxonomy_target))}",
                    ],
                    "outputs": {"taxonomy": str(taxonomy_target)},
                }
            )

        return steps
