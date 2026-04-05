from chimera_bench.tools.bracken import BrackenTool


def test_bracken_build_db_steps_stage_source_db(tmp_path):
    tool = BrackenTool({"env": "kraken", "bin": "bracken"})
    source_db = tmp_path / "kraken2" / "cami_refseq"
    source_db.mkdir(parents=True)
    steps = tool.build_db_steps(
        build={
            "db_prefix": "DB/cami_refseq",
            "threads": 32,
            "build": {
                "target_tsv": "/data/target.tsv",
                "source_db_prefix": str(source_db),
                "cleanup_library_fna": True,
            },
        },
        out_dir="/tmp/out",
    )

    assert steps[0]["name"] == "stage_kraken2_db"
    stage_cmd = steps[0]["cmd"]
    assert stage_cmd[0] == "bash"
    assert stage_cmd[1] == "-lc"
    assert f"ln {source_db}/hash.k2d /tmp/out/DB/cami_refseq/hash.k2d" in stage_cmd[2]
    assert f"ln -s {source_db}/taxonomy /tmp/out/DB/cami_refseq/taxonomy" in stage_cmd[2]
    assert "database100mers.kmer_distrib" in stage_cmd[2]

    step_names = [s["name"] for s in steps]
    assert step_names[1:5] == ["prepare_library", "build_database", "kmer2read_distr", "generate_distrib"]
    assert step_names[-1] == "cleanup_library"


def test_bracken_build_db_steps_without_source_db():
    tool = BrackenTool({"env": "kraken", "bin": "bracken"})
    steps = tool.build_db_steps(
        build={
            "db_prefix": "DB/cami_refseq",
            "threads": 32,
            "build": {
                "target_tsv": "/data/target.tsv",
                "cleanup_library_fna": False,
            },
        },
        out_dir="/tmp/out",
    )

    assert steps[0]["name"] == "ensure_db_dir"
    ensure_cmd = steps[0]["cmd"]
    assert ensure_cmd[0] == "bash"
    assert ensure_cmd[1] == "-lc"
    assert "mkdir -p /tmp/out/DB/cami_refseq" in ensure_cmd[2]
    assert steps[-1]["name"] == "generate_distrib"
