from __future__ import annotations

import hashlib
from pathlib import Path

from chimera_bench.real_fna_tables import (
    build_audit_rows,
    build_reference_rows,
    build_sample_manifest_rows,
    build_signal_rows,
)


def test_build_sample_manifest_joins_paper_metadata_and_input_depth() -> None:
    source_rows = [
        {
            "sample": "2015_Yu_SAMPLE01_c2Positive_head3000k",
            "cohort": "2015 Yu et al",
            "condition": "disease",
            "role": "paper_c2_positive",
            "run_accession": "ERR000001",
            "paper_fna_c1_pct": "0.1",
            "paper_fna_c2_pct": "0.2",
            "paper_fna_total_pct": "0.3",
        }
    ]
    metadata_rows = [
        {
            "cohort": "2015 Yu et al",
            "paper_sample": "SAMPLE01",
            "condition": "disease",
            "age": 55.0,
            "bmi": 24.5,
            "sex": "male",
            "source_reads": 12_000_000,
            "published_fna_c1_percent": 0.1,
            "published_fna_c2_percent": 0.2,
        }
    ]
    signal_rows = [
        {
            "tool": "chimera_default",
            "sample": "2015_Yu_SAMPLE01_c2Positive_head3000k",
            "total_or_detected": "3000000",
            "expected_c2_reads_at_depth": "6000",
        }
    ]

    rows = build_sample_manifest_rows(metadata_rows, source_rows, signal_rows)

    assert rows == [
        {
            "sample": "2015_Yu_SAMPLE01_c2Positive_head3000k",
            "paper_sample_id": "SAMPLE01",
            "run_accession": "ERR000001",
            "cohort": "2015 Yu et al",
            "condition": "disease",
            "role": "paper_c2_positive",
            "age": 55.0,
            "sex": "male",
            "bmi": 24.5,
            "source_reads": 12000000,
            "input_reads_used": 3000000,
            "input_selection": "R1_head3m_or_all_available_if_short",
            "paper_fna_c1_percent": 0.1,
            "paper_fna_c2_percent": 0.2,
            "paper_fna_total_percent": 0.3,
            "expected_c2_reads_at_depth": 6000.0,
            "input_qc": "pass",
        }
    ]


def test_build_signal_and_audit_rows_use_plot_ready_units() -> None:
    manifest = {
        "sample-a": {
            "sample": "sample-a",
            "cohort": "YachidaS_2019",
            "condition": "control",
            "role": "paper_zero_fna",
            "input_reads_used": 2_000_000,
            "paper_fna_c2_percent": 0.0,
        }
    }
    signal_rows = [
        {
            "tool": "kraken2_safe_lf01",
            "sample": "sample-a",
            "unit": "reads_per_million",
            "Fna_C2_signal": "4.5",
            "Fna_C1_signal": "1.0",
            "non_C1C2_signal": "2.0",
            "other_or_unmapped_signal": "999992.5",
            "expected_c2_reads_at_depth": "0",
        }
    ]
    audit_rows = [
        {
            "sample": "sample-a",
            "selected_reads": "20",
            "strict_c2_best_reads": "4",
            "c2_tied_best_reads": "3",
            "c2_supported_reads": "7",
            "non_c2_best_reads": "5",
            "unmapped_reads": "8",
            "strict_c2_best_rate": "0.2",
            "c2_tied_best_rate": "0.15",
            "c2_supported_rate": "0.35",
            "non_c2_best_rate": "0.25",
            "unmapped_rate": "0.4",
        }
    ]

    signals = build_signal_rows(signal_rows, manifest)
    audits = build_audit_rows("Kraken2_LF01", audit_rows, manifest)

    assert signals[0]["tool"] == "Kraken2_LF01"
    assert signals[0]["fna_c2_signal"] == 4.5
    assert audits[0]["candidate_reads_per_million"] == 10.0
    assert audits[0]["strict_c2_best_reads_per_million"] == 2.0
    assert audits[0]["c2_supported_reads_per_million"] == 3.5


def test_build_reference_rows_removes_local_path_and_hashes_fasta(tmp_path: Path) -> None:
    fasta = tmp_path / "GCF_000001.1_genomic.fna"
    fasta.write_text(">contig-a\nACGT\n>contig-b\nAAA\n")

    rows = build_reference_rows(
        [
            {
                "path": str(fasta),
                "taxid": "851",
                "paper_clade_projection": "non_C1C2",
                "reference_name": fasta.name,
            }
        ]
    )

    assert rows == [
        {
            "reference_name": "GCF_000001.1_genomic.fna",
            "accession": "GCF_000001.1",
            "taxid": 851,
            "paper_clade": "non_C1C2",
            "source_kind": "RefSeq",
            "sequence_records": 2,
            "sequence_bases": 7,
            "sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
        }
    ]
