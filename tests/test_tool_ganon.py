from chimera_bench.tools.ganon import GanonTool


def test_ganon_build_steps_single_reads():
    tool = GanonTool({"env": "ganon", "bin": "ganon"})
    steps = tool.build_steps(
        dataset={"reads": ["/r1.fq"]},
        exp={"db": "/db/prefix", "threads": 8, "tool_args": ["--keep-zeros"]},
        out_prefix="/tmp/out/ganon",
    )
    assert steps[0]["cmd"][:6] == ["conda", "run", "-n", "ganon", "ganon", "classify"]
    assert "--db-prefix" in steps[0]["cmd"]
    assert "/db/prefix" in steps[0]["cmd"]
    assert "--single-reads" in steps[0]["cmd"]
    assert "/r1.fq" in steps[0]["cmd"]
    assert "--output-prefix" in steps[0]["cmd"]
    assert "/tmp/out/ganon" in steps[0]["cmd"]
    assert "--threads" in steps[0]["cmd"]
    assert "8" in steps[0]["cmd"]
    assert "--skip-report" in steps[0]["cmd"]
    assert "--keep-zeros" in steps[0]["cmd"]

    assert steps[1]["cmd"][4:6] == ["ganon", "report"]
    assert "--report-type" in steps[1]["cmd"]
    assert "reads" in steps[1]["cmd"]

    assert steps[2]["cmd"][4:6] == ["ganon", "report"]
    assert "--report-type" in steps[2]["cmd"]
    assert "abundance" in steps[2]["cmd"]

    assert steps[1]["outputs"]["report_reads_tre"].endswith("_reads.tre")
    assert steps[2]["outputs"]["report_abundance_tre"].endswith("_abundance.tre")


def test_ganon_build_db_custom():
    tool = GanonTool({"env": "ganon", "bin": "ganon"})
    steps = tool.build_db_steps(
        build={
            "db_prefix": "/db/prefix",
            "threads": 4,
            "build": {
                "mode": "custom",
                "input": ["/ref1.fa", "/ref2.fa"],
                "input_target": "sequence",
                "level": "species",
                "taxonomy_files": ["/tax/nodes.dmp", "/tax/names.dmp"],
                "args": ["--min-length", "100"],
            },
        },
        out_dir="/tmp/out",
    )
    cmd = steps[0]["cmd"]
    assert cmd[:6] == ["conda", "run", "-n", "ganon", "ganon", "build-custom"]
    assert "--db-prefix" in cmd
    assert "/db/prefix" in cmd
    assert "--threads" in cmd
    assert "4" in cmd
    assert "--input" in cmd
    assert "/ref1.fa" in cmd
    assert "/ref2.fa" in cmd
    assert "--input-target" in cmd
    assert "sequence" in cmd
    assert "--level" in cmd
    assert "species" in cmd
    assert "--taxonomy-files" in cmd
    assert "/tax/nodes.dmp" in cmd
    assert "/tax/names.dmp" in cmd
    assert "--min-length" in cmd
    assert "100" in cmd
    assert steps[0]["outputs"]["db_prefix"] == "/db/prefix"


def test_ganon_output_one_auto_disabled_without_truth_mapping():
    tool = GanonTool({"env": "ganon", "bin": "ganon"})
    steps = tool.build_steps(
        dataset={"reads": ["/r1.fq"]},
        exp={"db": "/db/prefix", "threads": 8, "tool_args": []},
        out_prefix="/tmp/out/ganon",
    )
    assert "--output-one" not in steps[0]["cmd"]
    assert "classify_one" not in steps[0]["outputs"]


def test_ganon_output_one_auto_enabled_with_truth_dir():
    tool = GanonTool({"env": "ganon", "bin": "ganon"})
    steps = tool.build_steps(
        dataset={"reads": ["/r1.fq"], "truth_dir": "/truth"},
        exp={"db": "/db/prefix", "threads": 8, "tool_args": []},
        out_prefix="/tmp/out/ganon",
    )
    assert "--output-one" in steps[0]["cmd"]
    assert steps[0]["outputs"]["classify_one"].endswith(".one")
