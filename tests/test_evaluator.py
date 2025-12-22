from chimera_bench.core.evaluator import summarize_classify_tsv


def test_summarize_classify_tsv_counts_unclassified(tmp_path):
    p = tmp_path / "out.tsv"
    p.write_text("""r1\t123:1\n""" """r2\tunclassified\n""")

    metrics = summarize_classify_tsv(p)

    assert metrics["total_reads"] == 2
    assert metrics["unclassified_reads"] == 1
    assert metrics["classified_reads"] == 1
    assert metrics["unique_taxids"] == 1
