from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import pytest

from chimera_bench import cli as cli_mod
from chimera_bench.cli import paper_freeze_cmd
from chimera_bench.paper_freeze import (
    _dataset_sample_counts,
    _real_manifest_rows,
    _validate_complete_matrix,
    build_classify_rows,
    build_classify_sample_rows,
    build_profile_rows,
)


def _record(
    *,
    dataset: str,
    tool: str,
    metric_prefix: str = "per_read",
    collection: str | None = None,
    sample_id: str | None = None,
) -> dict:
    metrics = {
        "run_elapsed_seconds": 12.5,
        "resource_max_rss_gb": 3.25,
    }
    if metric_prefix == "per_read":
        metrics.update(
            {
                "per_read_truth_mapped_rate_species": 0.9,
                "per_read_pred_mapped_rate_species": 0.8,
                "per_read_precision_species": 0.7,
                "per_read_recall_species": 0.6,
                "per_read_f1_species": 0.65,
                "per_read_truth_mapped_rate_genus": 0.95,
                "per_read_pred_mapped_rate_genus": 0.85,
                "per_read_precision_genus": 0.8,
                "per_read_recall_genus": 0.75,
                "per_read_f1_genus": 0.77,
            }
        )
    else:
        metrics.update(
            {
                "completeness_species": 0.8,
                "purity_species": 0.7,
                "l1_norm_species": 0.4,
                "completeness_genus": 0.9,
                "purity_genus": 0.8,
                "l1_norm_genus": 0.3,
                "weighted_unifrac": 0.2,
                "profile_metric_version": "opal-core-v1",
            }
        )
    return {
        "tool": tool,
        "dataset": dataset,
        "dataset_collection": collection,
        "display_dataset": collection,
        "sample_id": sample_id,
        "db_name": "cami_refseq",
        "metrics": metrics,
    }


def test_classify_paper_rows_exclude_ganon_sample0_and_legacy_prjna() -> None:
    records = [
        _record(dataset="cami2-marine-long", tool=tool)
        for tool in ("chimera", "centrifuger", "kraken2", "ganon2")
    ]
    records.extend(
        [
            _record(dataset="cami2-marine-long-sample0", tool="chimera"),
            _record(dataset="prjna637878-supported19", tool="chimera"),
        ]
    )

    rows = build_classify_rows(records, {"cami2-marine-long": 10})

    assert [(row["dataset"], row["tool"]) for row in rows] == [
        ("cami2-marine-long", "centrifuger"),
        ("cami2-marine-long", "chimera"),
        ("cami2-marine-long", "kraken2"),
    ]
    assert {row["samples"] for row in rows} == {10}
    assert {row["paper_scope"] for row in rows} == {"main"}


def test_profile_paper_rows_are_supplementary_and_exclude_ganon() -> None:
    records = [
        _record(dataset="atcc-hifi", tool=tool, metric_prefix="profile")
        for tool in ("bracken", "centrifuger", "chimera", "sylph", "ganon2")
    ]
    records.append(
        _record(dataset="cami2-marine-long-sample0", tool="chimera", metric_prefix="profile")
    )

    rows = build_profile_rows(records, {"atcc-hifi": 1})

    assert [row["tool"] for row in rows] == ["bracken", "centrifuger", "chimera", "sylph"]
    assert {row["paper_scope"] for row in rows} == {"supplementary"}


def test_classify_sample_rows_only_export_complete_prjna_single_read_lane() -> None:
    records = [
        _record(
            dataset=f"prjna637878-supported19-single-read.{sample}",
            collection="prjna637878-supported19-single-read",
            sample_id=sample,
            tool=tool,
        )
        for tool in ("centrifuger", "chimera", "kraken2", "ganon2")
        for sample in ("sample-a", "sample-b")
    ]
    records.append(
        _record(
            dataset="prjna637878-supported19.sample-a",
            collection="prjna637878-supported19",
            sample_id="sample-a",
            tool="chimera",
        )
    )

    rows = build_classify_sample_rows(records)

    assert len(rows) == 6
    assert {row["dataset"] for row in rows} == {"prjna637878-supported19-single-read"}
    assert {row["tool"] for row in rows} == {"centrifuger", "chimera", "kraken2"}
    assert {row["sample_id"] for row in rows} == {"sample-a", "sample-b"}


def test_paper_freeze_command_uses_requested_roots(tmp_path: Path, monkeypatch) -> None:
    calls = []
    monkeypatch.setattr(
        cli_mod,
        "write_paper_tables",
        lambda **kwargs: calls.append(kwargs),
    )

    paper_freeze_cmd(
        Namespace(
            config=str(tmp_path / "configs"),
            results_root=str(tmp_path / "results"),
        )
    )

    assert calls == [
        {
            "config_root": tmp_path / "configs",
            "results_root": tmp_path / "results",
        }
    ]


def test_dataset_sample_counts_use_multi_file_dataset_units(tmp_path: Path) -> None:
    datasets = tmp_path / "configs" / "datasets"
    datasets.mkdir(parents=True)
    (datasets / "marine.yaml").write_text(
        "\n".join(
            [
                "name: marine",
                "group: cami2-marine",
                "reads:",
                "  - sample0.fasta",
                "  - sample1.fasta",
            ]
        )
        + "\n"
    )

    assert _dataset_sample_counts(tmp_path / "configs")["marine"] == 2


def test_real_manifest_records_sample_signal_and_read_audit_layers(tmp_path: Path) -> None:
    real = tmp_path / "results" / "real" / "fna_c2_crc3_head3m"
    real.mkdir(parents=True)
    (real / "sample_level_signals.tsv").write_text(
        "sample\ttool\n"
        "sample-a\tChimera\n"
        "sample-a\tCentrifuger\n"
        "sample-a\tKraken2_LF01\n"
        "sample-a\tsylph\n"
        "sample-b\tChimera\n"
        "sample-b\tCentrifuger\n"
        "sample-b\tKraken2_LF01\n"
        "sample-b\tsylph\n"
    )
    (real / "read_audit_sample_metrics.tsv").write_text(
        "sample\ttool\n"
        "sample-a\tChimera\n"
        "sample-a\tCentrifuger\n"
        "sample-a\tKraken2_LF01\n"
        "sample-b\tChimera\n"
        "sample-b\tCentrifuger\n"
        "sample-b\tKraken2_LF01\n"
    )

    rows = _real_manifest_rows(tmp_path / "results")

    assert len(rows) == 7
    assert {row["samples"] for row in rows} == {2}
    assert sum(row["task"] == "real_read_audit" for row in rows) == 3
    assert next(row for row in rows if row["tool"] == "Chimera")["tool_version"] == "1.6.3"
    assert next(row for row in rows if row["tool"] == "Kraken2_LF01")["tool_version"] == "2.1.3"


def test_real_manifest_rejects_equal_sized_but_different_sample_sets(tmp_path: Path) -> None:
    real = tmp_path / "results" / "real" / "fna_c2_crc3_head3m"
    real.mkdir(parents=True)
    (real / "sample_level_signals.tsv").write_text(
        "sample\ttool\n"
        "signal-only\tChimera\n"
        "signal-only\tCentrifuger\n"
        "signal-only\tKraken2_LF01\n"
        "signal-only\tsylph\n"
    )
    (real / "read_audit_sample_metrics.tsv").write_text(
        "sample\ttool\n"
        "audit-only\tChimera\n"
        "audit-only\tCentrifuger\n"
        "audit-only\tKraken2_LF01\n"
    )

    with pytest.raises(ValueError, match="sample sets differ"):
        _real_manifest_rows(tmp_path / "results")


def test_real_manifest_rejects_incomplete_per_tool_coverage(tmp_path: Path) -> None:
    real = tmp_path / "results" / "real" / "fna_c2_crc3_head3m"
    real.mkdir(parents=True)
    (real / "sample_level_signals.tsv").write_text(
        "sample\ttool\n"
        "sample-a\tChimera\n"
        "sample-a\tCentrifuger\n"
        "sample-a\tKraken2_LF01\n"
        "sample-a\tsylph\n"
        "sample-b\tCentrifuger\n"
        "sample-b\tKraken2_LF01\n"
        "sample-b\tsylph\n"
    )
    (real / "read_audit_sample_metrics.tsv").write_text(
        "sample\ttool\n"
        "sample-a\tChimera\n"
        "sample-a\tCentrifuger\n"
        "sample-a\tKraken2_LF01\n"
        "sample-b\tChimera\n"
        "sample-b\tCentrifuger\n"
        "sample-b\tKraken2_LF01\n"
    )

    with pytest.raises(ValueError, match="sample-level signal.*coverage"):
        _real_manifest_rows(tmp_path / "results")


def test_classify_rows_reject_partial_collection_marked_as_complete() -> None:
    records = [
        _record(
            dataset="sample-a",
            collection="cami-strain-madness-long",
            sample_id="sample-a",
            tool="chimera",
        )
    ]

    with pytest.raises(ValueError, match="expected 100 completed samples, found 1"):
        build_classify_rows(records, {"cami-strain-madness-long": 100})


def test_complete_matrix_rejects_missing_tool_dataset_pair() -> None:
    rows = [
        {"dataset": "dataset-a", "tool": "chimera"},
        {"dataset": "dataset-a", "tool": "centrifuger"},
    ]

    with pytest.raises(ValueError, match="missing=.*kraken2"):
        _validate_complete_matrix(
            rows,
            datasets=("dataset-a",),
            tools=("chimera", "centrifuger", "kraken2"),
            label="classify",
        )
