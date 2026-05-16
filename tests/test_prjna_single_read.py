from __future__ import annotations

import gzip
from pathlib import Path

import yaml

from chimera_bench.prjna_single_read import build_single_read_assets, mate_marked_read_id


def _write_gzip(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write(text)


def test_mate_marked_read_id_replaces_existing_suffix() -> None:
    assert mate_marked_read_id("@SRR1.7/2 comment", "R1") == "SRR1.7__mate1"
    assert mate_marked_read_id("SRR1.7__mate1", "R2") == "SRR1.7__mate2"


def test_build_single_read_assets_rewrites_fastq_and_truth(tmp_path: Path) -> None:
    source_root = tmp_path / "source"
    fastq = source_root / "fastq"
    truth = source_root / "truth"
    meta = source_root / "meta"
    meta.mkdir(parents=True)

    r1 = fastq / "S1_R1.fastq.gz"
    r2 = fastq / "S1_R2.fastq.gz"
    _write_gzip(r1, "@SRR1.1 extra\nACGT\n+\nIIII\n")
    _write_gzip(r2, "@SRR1.1\nTGCA\n+\nJJJJ\n")
    truth.mkdir(parents=True)
    (truth / "S1.species_truth.tsv").write_text(
        "\n".join(
            [
                "read_id\tmate\tspecies_label\tbest_AS\tbest_mapq\tn_top_hits",
                "SRR1.1\tR1\tSpecies A\t1\t60\t1",
                "SRR1.1\tR2\tSpecies A\t1\t60\t1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (meta / "PRJNA637878_supported19_sample_manifest.tsv").write_text(
        "sample_id\tbiosample\nS1\tSAMN1\n", encoding="utf-8"
    )
    (meta / "SAMN1.runinfo.csv").write_text(
        "Run,spots,spots_with_mates\nSRR1,2,1\n", encoding="utf-8"
    )

    cfg_path = tmp_path / "dataset.yaml"
    cfg_path.write_text(
        yaml.safe_dump(
            {
                "name": "source",
                "samples": [
                    {
                        "sample_id": "S1",
                        "paired": [str(r1), str(r2)],
                    }
                ],
            },
            sort_keys=False,
            allow_unicode=True,
        ),
        encoding="utf-8",
    )

    out_root = tmp_path / "derived"
    audit = tmp_path / "audit.tsv"
    rows = build_single_read_assets(
        source_config=cfg_path,
        source_root=source_root,
        out_root=out_root,
        audit_path=audit,
        selected=set(),
        pigz_threads=1,
        force=True,
        check_only=False,
    )

    assert rows[0]["single_read_records"] == 2
    with gzip.open(out_root / "fastq" / "S1_R1.single.fastq.gz", "rt", encoding="utf-8") as fh:
        assert fh.readline().strip() == "@SRR1.1__mate1"
    with gzip.open(out_root / "fastq" / "S1_R2.single.fastq.gz", "rt", encoding="utf-8") as fh:
        assert fh.readline().strip() == "@SRR1.1__mate2"

    single_truth = (out_root / "truth" / "S1.single_read_truth.tsv").read_text(encoding="utf-8")
    assert "SRR1.1__mate1\tSpecies A\tSRR1.1\tR1" in single_truth
    assert "SRR1.1__mate2\tSpecies A\tSRR1.1\tR2" in single_truth
    assert "runinfo_unpaired_spots" in audit.read_text(encoding="utf-8")
