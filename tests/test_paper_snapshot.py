from __future__ import annotations

import bisect
import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FNA = RESULTS / "real/fna_c2_crc3_head3m"


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _auc(positive: list[float], negative: list[float]) -> float:
    ordered = sorted(negative)
    wins = 0.0
    for value in positive:
        lower = bisect.bisect_left(ordered, value)
        upper = bisect.bisect_right(ordered, value)
        wins += lower + 0.5 * (upper - lower)
    return wins / (len(positive) * len(negative))


def test_paper_manifest_contains_only_formal_classify_matrix() -> None:
    rows = _read(RESULTS / "paper_run_manifest.tsv")
    classify = [row for row in rows if row["task"] == "classify"]

    assert len(classify) == 15
    assert {row["tool"] for row in classify} == {"centrifuger", "chimera", "kraken2"}
    assert all(row["status"] == "complete" and row["paper_scope"] == "main" for row in classify)
    assert not any("sample0" in row["dataset"] for row in classify)
    assert not any(row["dataset"] == "prjna637878-supported19" for row in classify)
    assert not any(row["tool"] == "ganon2" for row in rows)

    prjna = [
        row for row in classify if row["dataset"] == "prjna637878-supported19-single-read"
    ]
    assert len(prjna) == 3
    assert {row["samples"] for row in prjna} == {"19"}


def test_paper_classify_tables_are_complete_and_plot_ready() -> None:
    summary = _read(RESULTS / "classify/summary.tsv")
    samples = _read(RESULTS / "classify/sample_metrics.tsv")

    assert len(summary) == 15
    assert all(row["per_read_f1_species"] and row["per_read_f1_genus"] for row in summary)
    assert Counter(row["tool"] for row in samples) == Counter(
        {"centrifuger": 19, "chimera": 19, "kraken2": 19}
    )
    assert len({(row["sample_id"], row["tool"]) for row in samples}) == 57


def test_fna_plot_tables_cover_all_samples_tools_and_references() -> None:
    manifest = _read(FNA / "sample_manifest.tsv")
    signals = _read(FNA / "sample_level_signals.tsv")
    audits = _read(FNA / "read_audit_sample_metrics.tsv")
    references = _read(FNA / "reference_manifest.tsv")

    sample_names = {row["sample"] for row in manifest}
    assert len(manifest) == len(sample_names) == 760
    assert Counter(row["role"] for row in manifest) == Counter(
        {"paper_zero_fna": 619, "paper_c2_positive": 137, "paper_fna_positive_c2_zero": 4}
    )
    assert all(row["input_qc"] == "pass" for row in manifest)
    assert sum(int(row["input_reads_used"]) < 3_000_000 for row in manifest) == 2

    assert len(signals) == 3040
    assert Counter(row["tool"] for row in signals) == Counter(
        {"Chimera": 760, "Centrifuger": 760, "Kraken2_LF01": 760, "sylph": 760}
    )
    assert {row["sample"] for row in signals} == sample_names

    assert len(audits) == 2280
    assert Counter(row["tool"] for row in audits) == Counter(
        {"Chimera": 760, "Centrifuger": 760, "Kraken2_LF01": 760}
    )
    for row in audits:
        assert (
            int(row["strict_c2_best_reads"])
            + int(row["c2_tied_best_reads"])
            + int(row["non_c2_best_reads"])
            + int(row["unmapped_reads"])
            == int(row["candidate_reads"])
        )

    assert len(references) == 68
    assert sum(int(row["sequence_records"]) for row in references) == 1196
    assert sum(int(row["sequence_bases"]) for row in references) == 507_495_341
    assert all(len(row["sha256"]) == 64 for row in references)

    for path in (
        FNA / "sample_manifest.tsv",
        FNA / "sample_level_signals.tsv",
        FNA / "read_audit_sample_metrics.tsv",
        FNA / "reference_manifest.tsv",
    ):
        assert "/mnt/" not in path.read_text()


def test_fna_long_table_reproduces_primary_strict_auc() -> None:
    rows = [
        row
        for row in _read(FNA / "read_audit_sample_metrics.tsv")
        if row["tool"] == "Chimera"
    ]
    positive = [
        float(row["strict_c2_best_reads_per_million"])
        for row in rows
        if row["role"] == "paper_c2_positive"
    ]
    negative = [
        float(row["strict_c2_best_reads_per_million"])
        for row in rows
        if row["role"] == "paper_zero_fna"
    ]

    assert abs(_auc(positive, negative) - 0.880753) < 5e-7


def test_public_readmes_use_reader_facing_scientific_language() -> None:
    readmes = [
        ROOT / "README.md",
        RESULTS / "README.md",
        RESULTS / "classify/README.md",
        RESULTS / "profile/README.md",
        RESULTS / "real/README.md",
    ]
    forbidden = (
        "Paper-ready snapshot",
        "paper-ready",
        "论文正文和作图只使用",
        "论文主表只引用",
        "论文补充表只引用",
        "不承担 Chimera 的主要性能结论",
        "这不是缺失运行",
        "retained for traceability",
        "formal cross-method superiority test",
        "### Main Claim",
    )

    for path in readmes:
        text = path.read_text()
        for phrase in forbidden:
            assert phrase not in text, f"{path}: internal wording remains: {phrase}"

    real_readme = (RESULTS / "real/README.md").read_text()
    assert "### Summary" in real_readme
    assert "sylph 输出样本层面的丰度估计，不提供逐 read assignment" in real_readme
