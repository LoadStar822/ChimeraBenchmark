from math import isclose
from pathlib import Path

from chimera_bench.core.metrics import evaluate_with_truth, parse_ganon_one, load_taxonomy


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

    assert isclose(metrics["abundance_l1_species"], 26.666666667, rel_tol=1e-6)
    assert isclose(metrics["abundance_tv_species"], 0.1333333333, rel_tol=1e-6)
    assert isclose(metrics["abundance_bc_species"], 0.1333333333, rel_tol=1e-6)
    assert isclose(metrics["presence_precision_species"], 1.0, rel_tol=1e-6)
    assert isclose(metrics["presence_recall_species"], 1.0, rel_tol=1e-6)


def test_parse_ganon_one_uses_file_mapping(tmp_path: Path):
    tax = tmp_path / "test.tax"
    tax.write_text(
        "\n".join(
            [
                "1\t1\tno rank\troot\t0",
                "2\t1\tgenus\tGenusA\t0",
                "3\t2\tspecies\tSpeciesA\t0",
                "GCF_000001.1\t3\tfile\tGCF_000001.1\t123",
            ]
        )
        + "\n"
    )
    taxonomy, file_map, _ = load_taxonomy(tax)

    one = tmp_path / "pred.one"
    one.write_text("r1\tGCF_000001.1\t42\n")

    preds = parse_ganon_one(one, file_map)

    assert preds["r1"] == 3


def test_evaluate_with_truth_profile_only(tmp_path: Path):
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

    profile = tmp_path / "truth.tsv"
    profile.write_text("species\tpercent\nSpeciesA\t80\nSpeciesB\t20\n")

    tre = tmp_path / "pred.tre"
    tre.write_text(
        "\n".join(
            [
                "species\t3\t1|2|3\tSpeciesA\t0\t0\t12\t12\t80.0",
                "species\t5\t1|4|5\tSpeciesB\t0\t0\t3\t3\t20.0",
            ]
        )
        + "\n"
    )

    exp = {"taxonomy": str(tax)}
    dataset = {"truth_profile": str(profile)}
    outputs = {"report_abundance_tre": str(tre)}

    metrics = evaluate_with_truth(exp, dataset, outputs)

    assert isclose(metrics["presence_precision_species"], 1.0, rel_tol=1e-6)
    assert metrics.get("per_read_precision_species") is None


def test_descendant_per_read_penalizes_ancestor(tmp_path: Path):
    tax = tmp_path / "test.tax"
    tax.write_text(
        "\n".join(
            [
                "1\t1\tno rank\troot\t0",
                "2\t1\tgenus\tGenusA\t0",
                "3\t2\tspecies\tSpeciesA\t0",
            ]
        )
        + "\n"
    )

    mapping = tmp_path / "mapping.tsv"
    mapping.write_text(
        "#anonymous_contig_id\tgenome_id\ttax_id\tcontig_id\tnumber_reads\tstart_position\tend_position\n"
        "c1\tOtu1\t3\tX\t10\t0\t0\n"
        "c2\tOtu2\t3\tY\t5\t0\t0\n"
    )

    classify = tmp_path / "pred.tsv"
    classify.write_text("c1\t3:1\n" "c2\t2:1\n")

    exp = {"taxonomy": str(tax)}
    dataset = {"truth_map": str(mapping)}
    outputs = {"classify_tsv": str(classify)}

    metrics = evaluate_with_truth(exp, dataset, outputs)

    assert isclose(metrics["per_read_precision_species"], 0.5, rel_tol=1e-6)
    assert isclose(metrics["per_read_recall_species"], 0.5, rel_tol=1e-6)
    assert isclose(metrics["exact_per_read_precision_species"], 0.5, rel_tol=1e-6)
    assert isclose(metrics["exact_per_read_recall_species"], 0.5, rel_tol=1e-6)


def test_truth_profile_uses_names_dmp_mapping(tmp_path: Path):
    tax = tmp_path / "test.tax"
    tax.write_text(
        "\n".join(
            [
                "1\t1\tno rank\troot\t0",
                "10\t1\tgenus\tGenusA\t0",
                "11\t10\tspecies\tSpeciesA\t0",
            ]
        )
        + "\n"
    )

    nodes = tmp_path / "nodes.dmp"
    nodes.write_text(
        "\n".join(
            [
                "1 | 1 | no rank |",
                "10 | 1 | genus |",
                "11 | 10 | species |",
                "20 | 1 | genus |",
                "21 | 20 | species |",
                "30 | 1 | genus |",
                "31 | 30 | species |",
            ]
        )
        + "\n"
    )

    names = tmp_path / "names.dmp"
    names.write_text(
        "\n".join(
            [
                "11 | SpeciesA | | scientific name |",
                "21 | SpeciesB | | scientific name |",
                "31 | SpeciesC | | scientific name |",
            ]
        )
        + "\n"
    )

    truth = tmp_path / "truth.tsv"
    truth.write_text("species\tpercent\nSpeciesA\t33.3333\nSpeciesB\t33.3333\nSpeciesC\t33.3333\n")

    tre = tmp_path / "pred.tre"
    tre.write_text("species\t11\t1|10|11\tSpeciesA\t0\t0\t1\t1\t100.0\n")

    exp = {"taxonomy": str(tax), "nodes_dmp": str(nodes), "names_dmp": str(names)}
    dataset = {"truth_profile": str(truth)}
    outputs = {"report_abundance_tre": str(tre)}

    metrics = evaluate_with_truth(exp, dataset, outputs)

    assert isclose(metrics["presence_recall_species"], 1.0 / 3.0, rel_tol=1e-6)


def test_evaluate_with_truth_accepts_legacy_sylph_profile_key(tmp_path: Path):
    tax = tmp_path / "test.tax"
    tax.write_text(
        "\n".join(
            [
                "1\t1\tno rank\troot\t0",
                "10\t1\tgenus\tGenusA\t0",
                "11\t10\tspecies\tSpeciesA\t0",
                "20\t1\tgenus\tGenusB\t0",
                "21\t20\tspecies\tSpeciesB\t0",
                "GenomeA\t11\tfile\tGenomeA\t0",
                "GenomeB\t21\tfile\tGenomeB\t0",
            ]
        )
        + "\n"
    )

    truth = tmp_path / "truth.tsv"
    truth.write_text("species\tpercent\nSpeciesA\t80\nSpeciesB\t20\n")

    profile = tmp_path / "sylph_profile.tsv"
    profile.write_text("genome_file\ttaxonomic_abundance\nGenomeA\t80\nGenomeB\t20\n")

    exp = {"taxonomy": str(tax)}
    dataset = {"truth_profile": str(truth)}
    outputs = {"sylph_profile_tsv": str(profile)}

    metrics = evaluate_with_truth(exp, dataset, outputs)

    assert isclose(metrics["presence_precision_species"], 1.0, rel_tol=1e-6)


def test_abundance_metrics_penalize_unmapped_taxa_and_no_unk_metrics(tmp_path: Path):
    tax = tmp_path / "test.tax"
    tax.write_text(
        "\n".join(
            [
                "1\t1\tno rank\troot\t0",
                "2\t1\tgenus\tGenusA\t0",
                "3\t2\tspecies\tSpeciesA\t0",
            ]
        )
        + "\n"
    )

    truth = tmp_path / "truth.tsv"
    truth.write_text("species\tpercent\nSpeciesA\t100\n")

    # Include an unknown taxid (999) in the prediction. Bench-style descendant-aware
    # abundance metrics should treat it as an extra FP taxon (not drop it).
    tre = tmp_path / "pred.tre"
    tre.write_text(
        "\n".join(
            [
                "species\t3\t1|2|3\tSpeciesA\t0\t0\t50\t50\t50.0",
                "species\t999\t-\tUnknown\t0\t0\t50\t50\t50.0",
            ]
        )
        + "\n"
    )

    exp = {"taxonomy": str(tax)}
    dataset = {"truth_profile": str(truth)}
    outputs = {"report_abundance_tre": str(tre)}

    metrics = evaluate_with_truth(exp, dataset, outputs)

    assert isclose(metrics["presence_precision_species"], 0.5, rel_tol=1e-6)
    assert isclose(metrics["presence_recall_species"], 1.0, rel_tol=1e-6)
    assert isclose(metrics["abundance_l1_species"], 100.0, rel_tol=1e-6)
    assert all("_unk" not in key for key in metrics.keys())
