from chimera_bench.tools.taxor import TaxorTool


def test_taxor_build_db_steps_with_prep():
    tool = TaxorTool({"env": "taxor", "bin": "taxor"})
    steps = tool.build_db_steps(
        build={
            "db_prefix": "DB/cami_refseq",
            "threads": 128,
            "build": {
                "assembly_summary": "/data/assembly_summary.txt",
                "target_tsv": "/data/target.tsv",
                "input_sequence_dir": "/data/seqdir",
                "kmer_size": 22,
                "syncmer_size": 12,
                "use_syncmer": True,
            },
        },
        out_dir="/tmp/out",
    )

    assert [s["name"] for s in steps] == ["mkdir_db", "prepare_input", "build_db"]

    assert steps[0]["cmd"] == ["mkdir", "-p", "DB"]

    prep_cmd = steps[1]["cmd"]
    assert prep_cmd[0] == "python"
    assert prep_cmd[1].endswith("chimera_bench/tools/taxor_prep.py")
    assert "--assembly-summary" in prep_cmd
    assert "/data/assembly_summary.txt" in prep_cmd
    assert "--target-tsv" in prep_cmd
    assert "/data/target.tsv" in prep_cmd
    assert "--out" in prep_cmd
    out_path = prep_cmd[prep_cmd.index("--out") + 1]
    assert out_path.startswith("/tmp/out/")
    assert out_path.endswith("/outputs/taxor_input.tsv")

    cmd = steps[2]["cmd"]
    assert cmd[:6] == ["conda", "run", "-n", "taxor", "taxor", "build"]
    assert "--input-file" in cmd
    assert out_path in cmd
    assert "--input-sequence-dir" in cmd
    assert "/data/seqdir" in cmd
    assert "--output-filename" in cmd
    assert "DB/cami_refseq.hixf" in cmd
    assert "--threads" in cmd
    assert "32" in cmd
    assert "--kmer-size" in cmd
    assert "22" in cmd
    assert "--syncmer-size" in cmd
    assert "12" in cmd
    assert "--use-syncmer" in cmd
