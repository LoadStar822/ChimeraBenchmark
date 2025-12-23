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
