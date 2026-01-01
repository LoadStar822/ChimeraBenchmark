from pathlib import Path

from chimera_bench.tools.taxor_prep import write_taxor_input


def test_write_taxor_input_joins_target_and_assembly_summary(tmp_path: Path):
    assembly = tmp_path / "assembly_summary_refseq.txt"
    cols = ["x"] * 20
    cols[0] = "GCF_000001.1"
    cols[6] = "123"
    cols[19] = "https://example.org/genomes/all/GCF/000/001/GCF_000001.1_ASM1"
    assembly.write_text("# header\n" + "\t".join(cols) + "\n", encoding="utf-8")

    target = tmp_path / "target.tsv"
    target.write_text(
        "/data/GCF_000001.1_ASM1_genomic.fna.gz\t100\n"
        "/data/GCF_000002.1_ASM2_genomic.fna.gz\t200\n",
        encoding="utf-8",
    )

    out = tmp_path / "taxor_input.tsv"
    written = write_taxor_input(assembly_summary=assembly, target_tsv=target, out_path=out)

    assert written == 1
    assert out.read_text(encoding="utf-8") == "GCF_000001.1\t123\t/data/GCF_000001.1_ASM1_genomic.fna.gz\n"

