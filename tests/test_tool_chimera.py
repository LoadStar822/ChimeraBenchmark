from __future__ import annotations

from pathlib import Path

from chimera_bench.tools.chimera import ChimeraTool


def test_chimera_profile_uses_classify_evidence_tsv(tmp_path: Path):
    tool = ChimeraTool({"bin": "/repo/build/Chimera", "python": "python", "profile_script": "/repo/chimera.py"})
    out_prefix = str(tmp_path / "run" / "outputs" / "ChimeraClassify")
    profile_prefix = str(tmp_path / "profile" / "outputs" / "ChimeraClassify_abundance")

    steps = tool.build_steps(
        dataset={"reads": ["/reads/sample.fq"]},
        exp={"db": "/db/cami_refseq", "threads": 4},
        out_prefix=out_prefix,
        profile_out_prefix=profile_prefix,
    )

    assert steps[0]["outputs"]["classify_tsv"] == f"{out_prefix}.tsv"
    assert steps[0]["outputs"]["chimera_evidence_tsv"] == str(
        tmp_path / "run" / "outputs" / "ChimeraEvidence.tsv"
    )
    assert steps[1]["cmd"] == [
        "python",
        "/repo/chimera.py",
        "profile",
        "-i",
        str(tmp_path / "run" / "outputs" / "ChimeraEvidence.tsv"),
        "-o",
        str(Path(profile_prefix).resolve()),
    ]
    assert steps[1]["outputs"]["chimera_profile_tsv"] == f"{Path(profile_prefix).resolve()}.tsv"
