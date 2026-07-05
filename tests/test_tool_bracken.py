from chimera_bench.tools.bracken import BrackenTool


def test_bracken_build_steps_run_end_to_end_pipeline():
    tool = BrackenTool({"env": "kraken", "bin": "bracken"})
    steps = tool.build_steps(
        dataset={"reads": ["reads.fastq"]},
        exp={
            "db": "/db/cami_refseq",
            "threads": 32,
        },
        out_prefix="/tmp/runs/bracken/outputs/bracken",
        profile_out_prefix="/tmp/profile/bracken/outputs/bracken_abundance",
    )

    assert [step["name"] for step in steps] == ["kraken2_classify", "bracken", "convert"]

    classify_cmd = steps[0]["cmd"]
    assert classify_cmd[:6] == ["conda", "run", "-n", "kraken", "kraken2", "--db"]
    assert "/db/cami_refseq" in classify_cmd
    assert "--report" in classify_cmd
    report_path = classify_cmd[classify_cmd.index("--report") + 1]
    assert report_path.endswith("_kraken2.report")

    bracken_cmd = steps[1]["cmd"]
    assert bracken_cmd[:5] == ["conda", "run", "-n", "kraken", "bracken"]
    assert bracken_cmd[bracken_cmd.index("-i") + 1] == report_path
    assert "results/classify/kraken2" not in " ".join(bracken_cmd)


def test_bracken_build_steps_passes_optional_kraken2_args():
    tool = BrackenTool({"env": "kraken", "bin": "bracken"})
    steps = tool.build_steps(
        dataset={"reads": ["reads.fastq"]},
        exp={
            "db": "/db/cami_refseq",
            "threads": 32,
            "kraken2_tool_args": ["--memory-mapping"],
        },
        out_prefix="/tmp/runs/bracken/outputs/bracken",
        profile_out_prefix="/tmp/profile/bracken/outputs/bracken_abundance",
    )

    classify_cmd = steps[0]["cmd"]
    assert classify_cmd[-1] == "--memory-mapping"


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
