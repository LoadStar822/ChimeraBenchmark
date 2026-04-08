from chimera_bench.tools.ganon import GanonTool


def test_ganon_build_steps_do_not_enable_output_one_from_truth():
    tool = GanonTool({"env": "ganon", "bin": "ganon"})
    steps = tool.build_steps(
        dataset={"reads": ["reads.fq"], "truth_dir": "/truth/mapping"},
        exp={"db": "/db/cami_refseq", "threads": 32},
        out_prefix="/tmp/runs/ganon/outputs/ganon",
        profile_out_prefix="/tmp/profile/ganon/outputs/ganon_abundance",
    )

    classify_cmd = steps[0]["cmd"]
    assert "--output-one" not in classify_cmd
    assert "--multiple-matches" not in classify_cmd


def test_ganon_build_steps_use_explicit_output_flags_only():
    tool = GanonTool(
        {
            "env": "ganon",
            "bin": "ganon",
            "output_one": True,
            "output_unclassified": True,
            "skip_report": True,
        }
    )
    steps = tool.build_steps(
        dataset={"reads": ["reads.fq"]},
        exp={"db": "/db/cami_refseq", "threads": 32},
        out_prefix="/tmp/runs/ganon/outputs/ganon",
        profile_out_prefix="/tmp/profile/ganon/outputs/ganon_abundance",
    )

    classify_cmd = steps[0]["cmd"]
    assert "--output-one" in classify_cmd
    assert "--output-unclassified" in classify_cmd
    assert "--skip-report" in classify_cmd
    assert "--multiple-matches" not in classify_cmd
