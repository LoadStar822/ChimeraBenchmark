from pathlib import Path

from chimera_bench.core.reporter import write_summary


def test_write_summary(tmp_path: Path):
    runs = [
        {
            "exp": "e",
            "tool": "chimera",
            "dataset": "d",
            "metrics": {"total_reads": 2},
        },
    ]
    out = tmp_path / "summary.tsv"
    write_summary(runs, out)
    text = out.read_text()
    assert "total_reads" in text
    assert "2" in text
