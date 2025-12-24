from chimera_bench.core.evaluator import summarize_classify_tsv, summarize_ganon_tre


def test_summarize_classify_tsv_counts_unclassified(tmp_path):
    p = tmp_path / "out.tsv"
    p.write_text("""r1\t123:1\n""" """r2\tunclassified\n""")

    metrics = summarize_classify_tsv(p)

    assert metrics["total_reads"] == 2
    assert metrics["unclassified_reads"] == 1
    assert metrics["classified_reads"] == 1
    assert metrics["unique_taxids"] == 1


def test_summarize_ganon_tre_counts(tmp_path):
    p = tmp_path / "reads.tre"
    p.write_text(
        "\n".join(
            [
                "unclassified\t-\t-\tunclassified\t0\t0\t0\t3\t10.0",
                "root\t1\t1\troot\t0\t0\t7\t7\t90.0",
                "species\t123\t1|2|123\tFoo\t0\t0\t4\t4\t40.0",
                "species\t456\t1|2|456\tBar\t0\t0\t3\t3\t30.0",
            ]
        )
        + "\n"
    )

    metrics = summarize_ganon_tre(p)

    assert metrics["total_reads"] == 10
    assert metrics["unclassified_reads"] == 3
    assert metrics["classified_reads"] == 7
    assert metrics["unique_taxids"] == 2
