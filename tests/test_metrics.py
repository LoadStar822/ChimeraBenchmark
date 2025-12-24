from math import isclose
from pathlib import Path

from chimera_bench.core.metrics import evaluate_with_truth


def test_evaluate_with_truth_per_read_and_abundance(tmp_path: Path):
    tax = tmp_path / "test.tax"
    tax.write_text(
        "\n".join(
            [
                "1\t1\tno rank\troot\t0",
                "2\t1\tgenus\tGenusA\t0",
                "3\t2\tspecies\tSpeciesA\t0",
                "4\t1\tgenus\tGenusB\t0",
                "5\t4\tspecies\tSpeciesB\t0",
            ]
        )
        + "\n"
    )

    mapping = tmp_path / "mapping.tsv"
    mapping.write_text(
        "#anonymous_contig_id\tgenome_id\ttax_id\tcontig_id\tnumber_reads\tstart_position\tend_position\n"
        "c1\tOtu1\t3\tX\t10\t0\t0\n"
        "c2\tOtu2\t5\tY\t5\t0\t0\n"
    )

    classify = tmp_path / "pred.tsv"
    classify.write_text("c1\t3:1\n" "c2\tunclassified\n")

    tre = tmp_path / "pred.tre"
    tre.write_text(
        "\n".join(
            [
                "species\t3\t1|2|3\tSpeciesA\t0\t0\t12\t12\t80.0",
                "species\t5\t1|4|5\tSpeciesB\t0\t0\t3\t3\t20.0",
                "genus\t2\t1|2\tGenusA\t0\t0\t12\t12\t80.0",
                "genus\t4\t1|4\tGenusB\t0\t0\t3\t3\t20.0",
            ]
        )
        + "\n"
    )

    exp = {"taxonomy": str(tax)}
    dataset = {"truth_map": str(mapping)}
    outputs = {"classify_tsv": str(classify), "report_abundance_tre": str(tre)}

    metrics = evaluate_with_truth(exp, dataset, outputs)

    assert isclose(metrics["per_read_classified_rate"], 0.5, rel_tol=1e-6)
    assert isclose(metrics["per_read_precision_species"], 1.0, rel_tol=1e-6)
    assert isclose(metrics["per_read_recall_species"], 0.5, rel_tol=1e-6)

    assert isclose(metrics["abundance_l1_species"], 0.2666666667, rel_tol=1e-6)
    assert isclose(metrics["abundance_tv_species"], 0.1333333333, rel_tol=1e-6)
    assert isclose(metrics["abundance_bc_species"], 0.1333333333, rel_tol=1e-6)
    assert isclose(metrics["presence_precision_species"], 1.0, rel_tol=1e-6)
    assert isclose(metrics["presence_recall_species"], 1.0, rel_tol=1e-6)
