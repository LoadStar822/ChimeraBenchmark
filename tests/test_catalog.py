from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import pytest

from chimera_bench.catalog import collect_build_rows, collect_dataset_rows, scan_sequence_files
from chimera_bench.cli import catalog_cmd


def _write_assembly_summary(path: Path, rows: list[tuple[str, str, int, int]]) -> None:
    path.write_text(
        "\n".join(
            [
                "#assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tasm_submitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\texcluded_from_refseq\trelation_to_type_material\tasm_not_live_date\tassembly_type\tgroup\tgenome_size\tgenome_size_ungapped\tgc_percent\treplicon_count\tscaffold_count\tcontig_count\tannotation_provider\tannotation_name\tannotation_date\ttotal_gene_count\tprotein_coding_gene_count\tnon_coding_gene_count\tpubmed_id",
                *[
                    "\t".join(
                        [
                            accession,
                            "PRJNA0",
                            "SAMN0",
                            "na",
                            "na",
                            species_taxid,
                            species_taxid,
                            "Species",
                            "na",
                            "na",
                            "latest",
                            "Complete Genome",
                            "Major",
                            "Full",
                            "2026-01-01",
                            "ASM",
                            "Submitter",
                            "na",
                            "na",
                            "ftp",
                            "na",
                            "na",
                            "na",
                            "haploid",
                            "bacteria",
                            str(genome_size),
                            str(genome_size),
                            "50.0",
                            "1",
                            str(contig_count),
                            str(contig_count),
                            "NCBI",
                            "PGAP",
                            "2026-01-01",
                            "0",
                            "0",
                            "0",
                            "na",
                        ]
                    )
                    for accession, species_taxid, genome_size, contig_count in rows
                ],
            ]
        )
        + "\n"
    )


def test_scan_sequence_files_supports_paper_metrics(tmp_path: Path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@r1\nACGT\n+\nIIII\n@r2\nAA\n+\n!!\n")
    fasta = tmp_path / "contigs.fasta"
    fasta.write_text(">c1\nACGT\n>c2\nAA\nTT\n")

    stats = scan_sequence_files([fastq, fasta])

    assert stats[0]["format"] == "fastq"
    assert stats[0]["records"] == 2
    assert stats[0]["total_bases"] == 6
    assert stats[0]["mean_len"] == 3
    assert stats[0]["n50"] == 4
    assert stats[0]["gc_percent"] == pytest.approx(33.333333)
    assert stats[0]["q30_percent"] == pytest.approx(66.666667)
    assert stats[1]["format"] == "fasta"
    assert stats[1]["records"] == 2
    assert stats[1]["total_bases"] == 8
    assert stats[1]["mean_len"] == 4
    assert stats[1]["n50"] == 4
    assert stats[1]["gc_percent"] == pytest.approx(25.0)
    assert stats[1]["q30_percent"] is None


def test_collect_dataset_rows_outputs_paper_fields(tmp_path: Path):
    config_root = tmp_path / "configs"
    datasets_root = config_root / "datasets"
    datasets_root.mkdir(parents=True)

    r1 = tmp_path / "illumina_R1.fastq"
    r2 = tmp_path / "illumina_R2.fastq"
    r1.write_text("@a\nAAAA\n+\nIIII\n")
    r2.write_text("@a\nTTTT\n+\nIIII\n")
    truth_profile = tmp_path / "truth.tsv"
    truth_profile.write_text("species\tpercent\nSpeciesA\t100\n")
    (datasets_root / "atcc-illumina.yaml").write_text(
        "\n".join(
            [
                "name: atcc-illumina",
                "group: atcc",
                "platform: illumina",
                "paired:",
                f"  - {r1}",
                f"  - {r2}",
                f"truth_profile: {truth_profile}",
            ]
        )
        + "\n"
    )

    rows = collect_dataset_rows(
        config_root=config_root,
        cache_path=tmp_path / "results" / ".cache" / "catalog_cache.json",
    )
    row = rows[0]

    assert row["dataset_name"] == "ATCC MSA-1003 (Illumina)"
    assert row["samples"] == 1
    assert row["input_type"] == "paired FASTQ"
    assert row["reads_or_contigs"] == 2
    assert row["base_pairs_bp"] == 8
    assert row["mean_length_bp"] == 4
    assert row["n50_bp"] == 4
    assert row["gc_percent"] == "0.00"
    assert row["q30_percent"] == "100.00"
    assert row["truth"] == "profile only"
    assert "supports_classify" not in row
    assert "input_summary" not in row


def test_collect_dataset_rows_uses_collection_catalog_metadata(tmp_path: Path, monkeypatch):
    config_root = tmp_path / "configs"
    datasets_root = config_root / "datasets"
    datasets_root.mkdir(parents=True)

    r1 = tmp_path / "sample_R1.fastq.gz"
    r2 = tmp_path / "sample_R2.fastq.gz"
    r1.write_bytes(b"not-scanned")
    r2.write_bytes(b"not-scanned")
    truth = tmp_path / "truth.tsv"
    profile = tmp_path / "profile.tsv"
    truth.write_text("read_id\tspecies_label\nr1\tSpeciesA\n")
    profile.write_text("species_label\ttruth_abundance\nSpeciesA\t1.0\n")
    (datasets_root / "real.yaml").write_text(
        "\n".join(
            [
                "name: real",
                "group: real",
                "truth_label: local per-read + profile",
                "truth_map_format: species_label",
                "catalog:",
                "  total_size_gb: '1.2'",
                "  samples: 1",
                "  input_type: paired FASTQ",
                "  reads_or_contigs: 2",
                "  base_pairs_bp: 250",
                "  mean_length_bp: 125",
                "  truth: local per-read + profile",
                "samples:",
                "  - sample_id: S1",
                "    paired:",
                f"      - {r1}",
                f"      - {r2}",
                f"    truth_map: {truth}",
                f"    truth_profile: {profile}",
            ]
        )
        + "\n"
    )

    def fail_scan(*_args, **_kwargs):
        raise AssertionError("catalog metadata dataset should not scan sequence files")

    monkeypatch.setattr("chimera_bench.catalog.scan_sequence_group", fail_scan)

    rows = collect_dataset_rows(
        config_root=config_root,
        cache_path=tmp_path / "results" / ".cache" / "catalog_cache.json",
    )

    assert rows == [
        {
            "dataset": "real",
            "dataset_name": "real",
            "total_size_gb": "1.2",
            "samples": 1,
            "input_type": "paired FASTQ",
            "reads_or_contigs": 2,
            "base_pairs_bp": 250,
            "mean_length_bp": 125,
            "n50_bp": "—",
            "gc_percent": "—",
            "q30_percent": "—",
            "truth": "local per-read + profile",
        }
    ]


def test_collect_dataset_rows_skips_catalog_excluded_execution_variants(tmp_path: Path):
    config_root = tmp_path / "configs"
    datasets_root = config_root / "datasets"
    datasets_root.mkdir(parents=True)

    r1 = tmp_path / "batch_R1.fastq"
    r2 = tmp_path / "batch_R2.fastq"
    r1.write_text("@r/1\nA\n+\nI\n")
    r2.write_text("@r/2\nT\n+\nI\n")
    (datasets_root / "batch_variant.yaml").write_text(
        "\n".join(
            [
                "name: batch-variant",
                "group: strain-madness",
                "catalog:",
                "  exclude: true",
                "samples:",
                "  - sample_id: batch00",
                "    paired:",
                f"      - {r1}",
                f"      - {r2}",
            ]
        )
        + "\n"
    )

    rows = collect_dataset_rows(
        config_root=config_root,
        cache_path=tmp_path / "results" / ".cache" / "catalog_cache.json",
        dataset_names=["batch-variant"],
    )

    assert rows == []


def test_collect_dataset_rows_resolves_strain_source_and_text_catalog_fields(tmp_path: Path, monkeypatch):
    config_root = tmp_path / "configs"
    datasets_root = config_root / "datasets"
    datasets_root.mkdir(parents=True)

    source_root = tmp_path / "strain" / "long_read"
    reads = source_root / "strmgCAMI2_sample_0" / "long_read" / "run0" / "reads" / "anonymous_reads.fq"
    reads.parent.mkdir(parents=True)
    reads.write_text("@S0R0\nACGT\n+\nIIII\n")
    (datasets_root / "strain.yaml").write_text(
        "\n".join(
            [
                "name: strain",
                "group: strain-madness",
                "truth_label: per-read mapping; profile from read abundance",
                "sample_ids: [0]",
                "catalog:",
                "  total_size_gb: '4.0'",
                "  samples: 1",
                "  input_type: single FASTQ",
                "  reads_or_contigs: 1",
                "  base_pairs_bp: '—'",
                "  mean_length_bp: '—'",
                "  truth: per-read mapping; profile from read abundance",
                "source:",
                "  kind: strain_madness",
                "  read_type: long",
                f"  root: {source_root}",
                f"truth_dir: {source_root}",
            ]
        )
        + "\n"
    )

    def fail_scan(*_args, **_kwargs):
        raise AssertionError("catalog metadata dataset should not scan sequence files")

    monkeypatch.setattr("chimera_bench.catalog.scan_sequence_group", fail_scan)

    rows = collect_dataset_rows(
        config_root=config_root,
        cache_path=tmp_path / "results" / ".cache" / "catalog_cache.json",
    )

    assert rows[0]["reads_or_contigs"] == 1
    assert rows[0]["base_pairs_bp"] == "—"
    assert rows[0]["mean_length_bp"] == "—"
    assert rows[0]["truth"] == "per-read mapping; profile from read abundance"


def test_collect_build_rows_uses_reference_metadata(tmp_path: Path):
    config_root = tmp_path / "configs"
    builds_root = config_root / "build"
    builds_root.mkdir(parents=True)
    ref_root = tmp_path / "ref"
    ref_root.mkdir()

    ref1 = ref_root / "GCF_000000001.1_ASM1_genomic.fna.gz"
    ref2 = ref_root / "GCF_000000002.1_ASM2_genomic.fna.gz"
    ref1.write_text("aaa")
    ref2.write_text("bbbb")
    target = ref_root / "target.tsv"
    target.write_text(f"{ref1}\t111\n{ref2}\t222\n")
    _write_assembly_summary(
        ref_root / "assembly_summary.txt",
        [
            ("GCF_000000001.1", "10", 1000, 2),
            ("GCF_000000002.1", "10", 2000, 3),
        ],
    )

    (builds_root / "ganon_cami_refseq.yaml").write_text(
        "\n".join(
            [
                "name: ganon-cami-refseq",
                "tool: ganon",
                "db_prefix: DB/cami_refseq",
                "build:",
                f"  target_tsv: {target}",
                f"  taxonomy_dir: {tmp_path / 'tax'}",
            ]
        )
        + "\n"
    )

    rows = collect_build_rows(
        config_root=config_root,
        cache_path=tmp_path / "results" / ".cache" / "catalog_cache.json",
    )
    row = rows[0]

    assert row["dataset_name"] == "CAMI RefSeq"
    assert row["total_size_gb"] == "0.0"
    assert row["total_sequences"] == 5
    assert row["base_pairs_bp"] == 3000
    assert row["assemblies"] == 2
    assert row["species_count"] == 1
    assert "builder_tools" not in row


def test_collect_build_rows_fails_on_missing_assembly_metadata(tmp_path: Path):
    config_root = tmp_path / "configs"
    builds_root = config_root / "build"
    builds_root.mkdir(parents=True)
    ref_root = tmp_path / "ref"
    ref_root.mkdir()
    ref1 = ref_root / "GCF_000000001.1_ASM1_genomic.fna.gz"
    ref1.write_text("aaa")
    target = ref_root / "target.tsv"
    target.write_text(f"{ref1}\t111\n")
    _write_assembly_summary(ref_root / "assembly_summary.txt", [])
    (builds_root / "ganon_cami_refseq.yaml").write_text(
        "\n".join(
            [
                "name: ganon-cami-refseq",
                "tool: ganon",
                "db_prefix: DB/cami_refseq",
                "build:",
                f"  target_tsv: {target}",
            ]
        )
        + "\n"
    )

    with pytest.raises(ValueError, match="assembly summary incomplete"):
        collect_build_rows(
            config_root=config_root,
            cache_path=tmp_path / "results" / ".cache" / "catalog_cache.json",
        )


def test_catalog_cmd_writes_paper_reports_and_results_readme(tmp_path: Path):
    config_root = tmp_path / "configs"
    datasets_root = config_root / "datasets"
    builds_root = config_root / "build"
    datasets_root.mkdir(parents=True)
    builds_root.mkdir(parents=True)

    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n")
    truth_profile = tmp_path / "truth.tsv"
    truth_profile.write_text("species\tpercent\nSpeciesA\t100\n")
    (datasets_root / "atcc-hifi.yaml").write_text(
        "\n".join(
            [
                "name: atcc-hifi",
                "group: atcc",
                "platform: pacbio-hifi",
                "reads:",
                f"  - {reads}",
                f"truth_profile: {truth_profile}",
            ]
        )
        + "\n"
    )

    ref = tmp_path / "GCF_000000001.1_ASM1_genomic.fna.gz"
    ref.write_text("aaa")
    target = tmp_path / "target.tsv"
    target.write_text(f"{ref}\t1\n")
    _write_assembly_summary(tmp_path / "assembly_summary.txt", [("GCF_000000001.1", "1", 100, 1)])
    (builds_root / "ganon_cami_refseq.yaml").write_text(
        "\n".join(
            [
                "name: ganon-cami-refseq",
                "tool: ganon",
                "db_prefix: DB/cami_refseq",
                "build:",
                f"  target_tsv: {target}",
            ]
        )
        + "\n"
    )

    results_root = tmp_path / "results"
    resources_root = tmp_path / "resources"
    catalog_cmd(
        Namespace(
            config=str(config_root),
            results_root=str(results_root),
            resources_root=str(resources_root),
        )
    )

    db_catalog = (resources_root / "reports" / "db_catalog.tsv").read_text()
    dataset_catalog = (resources_root / "reports" / "dataset_catalog.tsv").read_text()
    results_readme = (results_root / "README.md").read_text()

    assert "dataset_name\ttotal_size_gb\ttotal_sequences" in db_catalog
    assert "CAMI RefSeq" in db_catalog
    assert "dataset_name\ttotal_size_gb\tsamples\tinput_type" in dataset_catalog
    assert "ATCC MSA-1003 (PacBio HiFi)" in dataset_catalog
    assert "## 目录" in results_readme
    assert "## 数据与结果导航" in results_readme
    assert "## 结果使用说明" in results_readme
    assert "## Reference Database Summary" in results_readme
    assert "## Benchmark Dataset Summary" in results_readme
    assert "## 评估任务与口径" in results_readme
    assert "## 指标说明" in results_readme
    assert "## 工具说明" in results_readme
    assert "Builder Tools" not in results_readme
    assert "Public Experiment Tools" not in results_readme
    assert "Input Summary" not in results_readme
    assert "GiB" not in results_readme
    assert "### Bracken" in results_readme
    assert "合同" not in results_readme
    assert "落盘" not in results_readme
    assert "机器可核对" not in results_readme
    assert not (results_root / "reports").exists()
    assert not (results_root / ".cache").exists()
