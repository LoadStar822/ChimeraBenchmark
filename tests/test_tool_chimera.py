from chimera_bench.tools.chimera import ChimeraTool


def test_chimera_build_cmd_single_reads():
    tool = ChimeraTool({"bin": "Chimera"})
    cmd, outputs = tool.build_cmd(
        dataset={"reads": ["/r1.fq"]},
        exp={"db": "/db/imcf", "threads": 16, "tool_args": ["--em"]},
        out_prefix="/tmp/out/ChimeraClassify",
    )
    assert cmd[:4] == ["Chimera", "classify", "-d", "/db/imcf"]
    assert "-i" in cmd and "/r1.fq" in cmd
    assert "-o" in cmd and "/tmp/out/ChimeraClassify" in cmd
    assert "-t" in cmd and "16" in cmd
    assert "--em" in cmd
    assert outputs["classify_tsv"].endswith(".tsv")
