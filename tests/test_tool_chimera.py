from __future__ import annotations

from chimera_bench.tools.chimera import ChimeraTool


def test_chimera_classify_uses_native_cami_profile(tmp_path):
    tool = ChimeraTool({"bin": "/repo/build/Chimera"})
    out_prefix = str(tmp_path / "run" / "outputs" / "ChimeraClassify")
    profile_prefix = str(tmp_path / "profile" / "outputs" / "ChimeraClassify_abundance")

    steps = tool.build_steps(
        dataset={"reads": ["/reads/sample.fq"]},
        exp={"db": "/db/cami_refseq", "threads": 4},
        out_prefix=out_prefix,
        profile_out_prefix=profile_prefix,
    )

    assert len(steps) == 1
    assert steps[0]["cmd"] == [
        "/repo/build/Chimera",
        "classify",
        "-d",
        "/db/cami_refseq",
        "-o",
        out_prefix,
        "-t",
        "4",
        "-i",
        "/reads/sample.fq",
        "--profile-cami",
    ]
    assert steps[0]["outputs"]["classify_tsv"] == f"{out_prefix}.tsv"
    assert steps[0]["outputs"]["chimera_profile_tsv"] == str(
        tmp_path / "run" / "outputs" / "ChimeraProfile.tsv"
    )
    assert steps[0]["outputs"]["cami_profile_tsv"] == str(
        tmp_path / "run" / "outputs" / "ChimeraProfile.cami.tsv"
    )


def test_chimera_build_outputs_database_directory(tmp_path):
    tool = ChimeraTool({"bin": "/repo/build/Chimera"})
    target = tmp_path / "target.tsv"
    target.write_text("/ref/genome.fa\t562\n")

    steps = tool.build_db_steps(
        build={
            "db_prefix": "DB/cami_refseq",
            "threads": 8,
            "build": {"input_tsv": str(target), "taxonomy_dir": "/taxdump"},
        },
        out_dir=str(tmp_path / "builds" / "chimera" / "cami_refseq"),
    )

    build_step = steps[1]
    assert build_step["cmd"] == [
        "/repo/build/Chimera",
        "build",
        "-i",
        str(target),
        "-o",
        str(tmp_path / "builds" / "chimera" / "cami_refseq" / "DB" / "cami_refseq"),
        "-t",
        "8",
        "--taxonomy-dir",
        "/taxdump",
    ]
    assert build_step["outputs"] == {
        "db_prefix": str(tmp_path / "builds" / "chimera" / "cami_refseq" / "DB" / "cami_refseq"),
        "db_dir": str(tmp_path / "builds" / "chimera" / "cami_refseq" / "DB" / "cami_refseq"),
        "db_file": str(tmp_path / "builds" / "chimera" / "cami_refseq" / "DB" / "cami_refseq" / "core.imcf"),
    }
