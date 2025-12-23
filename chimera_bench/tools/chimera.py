from __future__ import annotations

from typing import Any, Dict, List, Tuple


class ChimeraTool:
    name = "chimera"
    output_basename = "ChimeraClassify"

    def __init__(self, config: Dict[str, Any] | None = None) -> None:
        self.config = config or {}

    def build_cmd(
        self, *, dataset: Dict[str, Any], exp: Dict[str, Any], out_prefix: str
    ) -> Tuple[List[str], Dict[str, str]]:
        bin_path = self.config.get("bin", "Chimera")
        db = exp["db"]
        threads = str(exp.get("threads", 192))
        tool_args = list(exp.get("tool_args", []))

        cmd = [bin_path, "classify", "-d", db, "-o", out_prefix, "-t", threads]
        if "reads" in dataset:
            cmd += ["-i", *dataset["reads"]]
        elif "paired" in dataset:
            cmd += ["-p", *dataset["paired"]]
        else:
            raise ValueError("dataset must define reads or paired")
        cmd += tool_args

        outputs = {"classify_tsv": f"{out_prefix}.tsv"}
        return cmd, outputs
