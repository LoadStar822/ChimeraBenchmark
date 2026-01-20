from pathlib import Path

from chimera_bench.tools.centrifuger_convert import convert_centrifuger_output


def test_convert_centrifuger_output_dedup_and_unclassified(tmp_path: Path):
    raw = tmp_path / "raw.tsv"
    raw.write_text(
        "\n".join(
            [
                "readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches",
                "r1\tMT019531.1\t2697049\t4225\t0\t80\t80\t1",
                "r1\tALT\t123\t1\t0\t10\t80\t1",
                "r2\tX\t0\t0\t0\t0\t80\t0",
                "r3\tY\t-\t0\t0\t0\t80\t0",
                "",
            ]
        )
    )
    out = tmp_path / "classify.tsv"
    convert_centrifuger_output(input_path=raw, out_path=out)
    assert out.read_text() == "r1\t2697049\nr2\tunclassified\nr3\tunclassified\n"

