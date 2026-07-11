from __future__ import annotations

import sys
from pathlib import Path

import pytest


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts/fna_c2"
sys.path.insert(0, str(SCRIPTS))

from analyze_crc_association import (  # noqa: E402
    load_combined_sample_tables,
    load_frozen_metadata,
    match_metadata_rows,
)
from analyze_c2_specificity import main as specificity_main  # noqa: E402
from analyze_crc_sample_signal import read_signal_tables  # noqa: E402
from bootstrap_read_audit_auc import build_rows as build_bootstrap_rows  # noqa: E402
from make_detection_floor_tables import main as sample_summary_main  # noqa: E402


REAL_RESULTS = Path(__file__).resolve().parents[1] / "results/real/fna_c2_crc3_head3m"


def test_load_frozen_metadata_uses_public_sample_manifest(tmp_path: Path) -> None:
    manifest = tmp_path / "sample_manifest.tsv"
    manifest.write_text(
        "sample\tpaper_sample_id\tcohort\tcondition\tage\tbmi\tsex\tsource_reads"
        "\tpaper_fna_c1_percent\tpaper_fna_c2_percent\n"
        "sample-a\tPAPER-A\t2015 Yu et al\tdisease\t55\t24.5\tmale\t12000000\t0.1\t0.2\n"
    )

    rows = load_frozen_metadata(manifest)

    assert rows == [
        {
            "sample": "sample-a",
            "cohort": "2015 Yu et al",
            "paper_sample": "PAPER-A",
            "condition": "disease",
            "age": 55.0,
            "bmi": 24.5,
            "sex": "male",
            "source_reads": 12_000_000.0,
            "published_fna_c1_percent": 0.1,
            "published_fna_c2_percent": 0.2,
        }
    ]


def test_frozen_metadata_matches_exact_sample_not_name_fragment() -> None:
    metadata = [
        {
            "sample": "sample-a",
            "paper_sample": "PAPER-A",
            "cohort": "2015 Yu et al",
            "condition": "disease",
        }
    ]
    samples = [
        {
            "sample": "prefix_PAPER-A_suffix",
            "cohort": "2015 Yu et al",
            "condition": "disease",
        },
        {
            "sample": "sample-a",
            "cohort": "2015 Yu et al",
            "condition": "disease",
        },
    ]

    matched = match_metadata_rows(metadata, samples)

    assert len(matched) == 1
    assert matched[0]["sample"] == "sample-a"


def test_load_combined_sample_tables_groups_tools_and_rejects_duplicates(tmp_path: Path) -> None:
    table = tmp_path / "audit.tsv"
    table.write_text(
        "sample\ttool\tcohort\tcondition\n"
        "sample-a\tChimera\t2015 Yu et al\tdisease\n"
        "sample-a\tCentrifuger\t2015 Yu et al\tdisease\n"
    )

    grouped = load_combined_sample_tables([table])

    assert {tool: len(rows) for tool, rows in grouped.items()} == {
        "Chimera": 1,
        "Centrifuger": 1,
    }

    with pytest.raises(ValueError, match="duplicate tool/sample"):
        load_combined_sample_tables([table, table])


def test_read_signal_tables_accepts_frozen_normalized_columns(tmp_path: Path) -> None:
    table = tmp_path / "signals.tsv"
    table.write_text(
        "sample\ttool\tcohort\tcondition\tsignal_unit\tfna_c2_signal\n"
        "sample-a\tChimera\t2015 Yu et al\tdisease\treads_per_million\t12.5\n"
    )

    rows = read_signal_tables([table])

    assert rows == [
        {
            "sample": "sample-a",
            "tool": "Chimera",
            "cohort": "2015 Yu et al",
            "condition": "disease",
            "signal_unit": "reads_per_million",
            "fna_c2_signal": "12.5",
            "unit": "reads_per_million",
            "Fna_C2_signal": "12.5",
        }
    ]


def test_bootstrap_read_audit_builds_complete_three_tool_table() -> None:
    cohorts = ("2015 Yu et al", "2018 Wirbel et al", "YachidaS_2019")
    manifest = []
    audits = []
    for cohort_index, cohort in enumerate(cohorts):
        for role, suffix, expected in (
            ("paper_c2_positive", "positive", 6000.0),
            ("paper_zero_fna", "negative", 0.0),
        ):
            sample = f"sample-{cohort_index}-{suffix}"
            manifest.append(
                {
                    "sample": sample,
                    "cohort": cohort,
                    "role": role,
                    "expected_c2_reads_at_depth": str(expected),
                }
            )
            for tool_index, tool in enumerate(("Chimera", "Centrifuger", "Kraken2_LF01")):
                base = 10.0 + tool_index if role == "paper_c2_positive" else float(tool_index)
                audits.append(
                    {
                        "sample": sample,
                        "tool": tool,
                        "cohort": cohort,
                        "role": role,
                        "candidate_reads_per_million": str(base + 2.0),
                        "c2_supported_reads_per_million": str(base + 1.0),
                        "strict_c2_best_reads_per_million": str(base),
                    }
                )

    rows = build_bootstrap_rows(manifest, audits, iterations=20, seed=7)

    assert len(rows) == 45
    assert {row["tool"] for row in rows} == {"Chimera", "Centrifuger", "Kraken2_LF01"}
    assert all(row["bootstrap_iterations"] == 20 for row in rows)
    assert all(
        row["chimera_minus_tool_auc"] is not None
        for row in rows
        if row["tool"] != "Chimera"
    )


def test_sample_summary_script_reproduces_frozen_tables(tmp_path: Path) -> None:
    overall = tmp_path / "sample_level_overall_metrics.tsv"
    floors = tmp_path / "sample_level_detection_floor_auc.tsv"

    sample_summary_main(
        [
            "--sample-manifest",
            str(REAL_RESULTS / "sample_manifest.tsv"),
            "--signal-table",
            str(REAL_RESULTS / "sample_level_signals.tsv"),
            "--overall-out",
            str(overall),
            "--floor-out",
            str(floors),
        ]
    )

    assert overall.read_bytes() == (REAL_RESULTS / overall.name).read_bytes()
    assert floors.read_bytes() == (REAL_RESULTS / floors.name).read_bytes()


def test_specificity_script_reproduces_frozen_table(tmp_path: Path) -> None:
    specificity_main(
        [
            "--sample-manifest",
            str(REAL_RESULTS / "sample_manifest.tsv"),
            "--audit-table",
            str(REAL_RESULTS / "read_audit_sample_metrics.tsv"),
            "--signal-table",
            str(REAL_RESULTS / "sample_level_signals.tsv"),
            "--out-dir",
            str(tmp_path),
        ]
    )

    generated = tmp_path / "c2_specificity_partial_association.tsv"
    assert generated.read_bytes() == (REAL_RESULTS / generated.name).read_bytes()
