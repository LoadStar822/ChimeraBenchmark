from __future__ import annotations

import json
from argparse import Namespace
from pathlib import Path

import pytest

from chimera_bench import cli as cli_mod
from chimera_bench.cli import recompute_cmd, run_cmd
from chimera_bench.dataset_prepare import prepare_dataset_inputs


def test_recompute_rewrites_metrics_and_skips_failed_runs(tmp_path: Path):
    config_root = tmp_path / "configs"
    (config_root / "datasets").mkdir(parents=True)
    (config_root / "experiments").mkdir(parents=True)

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

    truth = tmp_path / "truth.tsv"
    truth.write_text("species\tpercent\nSpeciesA\t80\nSpeciesB\t20\n")

    mapping = tmp_path / "mapping.tsv"
    mapping.write_text(
        "#anonymous_contig_id\tgenome_id\ttax_id\tcontig_id\tnumber_reads\tstart_position\tend_position\n"
        "c1\tOtu1\t3\tX\t10\t0\t0\n"
        "c2\tOtu2\t5\tY\t5\t0\t0\n"
    )

    classify = tmp_path / "pred.tsv"
    classify.write_text("c1\t3:1\n" "c2\tunclassified\n")

    pred = tmp_path / "pred.tre"
    pred.write_text(
        "\n".join(
            [
                "species\t3\t1|2|3\tSpeciesA\t0\t0\t8\t8\t80.0",
                "species\t5\t1|4|5\tSpeciesB\t0\t0\t2\t2\t20.0",
            ]
        )
        + "\n"
    )

    (config_root / "datasets" / "ds1.yaml").write_text(
        "\n".join(
            [
                "truth_map: " + str(mapping),
                "truth_profile: " + str(truth),
            ]
        )
        + "\n"
    )
    (config_root / "datasets" / "ds_fail.yaml").write_text(
        "\n".join(
            [
                "truth_map: " + str(mapping),
                "truth_profile: " + str(truth),
            ]
        )
        + "\n"
    )
    (config_root / "experiments" / "chimera.yaml").write_text(
        "\n".join(
            [
                "tool: chimera",
                "taxonomy: " + str(tax),
                "datasets:",
                "  - ds1",
                "  - ds_fail",
            ]
        )
        + "\n"
    )

    runs_root = tmp_path / "results" / "classify"
    profile_root = tmp_path / "results" / "profile"
    profile_root.mkdir(parents=True)

    success_dir = runs_root / "chimera" / "ds1"
    success_dir.mkdir(parents=True)
    (success_dir / "meta.json").write_text(
        json.dumps(
            {
                "exp": "chimera",
                "tool": "chimera",
                "dataset": "ds1",
                "db_name": "cami_refseq",
                "return_code": 0,
                "elapsed_seconds": 12.0,
                "resource": {"max_rss_kb": 2048},
                "outputs": {"classify_tsv": str(classify), "report_abundance_tre": str(pred)},
            }
        )
    )

    failed_dir = runs_root / "chimera" / "ds_fail"
    failed_dir.mkdir(parents=True)
    (failed_dir / "meta.json").write_text(
        json.dumps(
            {
                "exp": "chimera",
                "tool": "chimera",
                "dataset": "ds_fail",
                "db_name": "cami_refseq",
                "return_code": 137,
                "elapsed_seconds": 99.0,
                "resource": {"max_rss_kb": 4096},
                "outputs": {"report_abundance_tre": str(pred)},
            }
        )
    )
    (failed_dir / "metrics.json").write_text(json.dumps({"presence_precision_species": 0.9}))

    out = tmp_path / "results" / "reports" / "chimera_summary.tsv"
    args = Namespace(
        exp="chimera",
        config=str(config_root),
        runs=str(runs_root),
        profile=str(profile_root),
        out=str(out),
        dataset=[],
    )

    recompute_cmd(args)

    metrics = json.loads((success_dir / "metrics.json").read_text())
    assert metrics["exact_per_read_precision_species"] == 1.0
    assert metrics["exact_per_read_recall_species"] == 0.5
    assert metrics["exact_per_read_truth_mapped_rate_species"] == 1.0
    assert metrics["completeness_species"] == 1.0
    assert metrics["purity_species"] == 1.0
    assert metrics["l1_norm_species"] == 0.0
    assert metrics["weighted_unifrac"] == 0.0

    classify_readme = (runs_root / "README.md").read_text()
    assert "Truth Mapped Rate (species)" in classify_readme
    assert "Pred Mapped Rate (species)" in classify_readme
    assert "ds_fail" not in classify_readme

    profile_readme = (profile_root / "README.md").read_text()
    assert "Weighted UniFrac" in profile_readme
    assert "ds_fail" not in profile_readme

    summary = out.read_text()
    assert "weighted_unifrac" in summary.splitlines()[0]
    assert "ds1" in summary
    assert "ds_fail" not in summary


def test_run_cmd_does_not_skip_large_taxor_dataset(tmp_path: Path, monkeypatch):
    config_root = tmp_path / "configs"
    (config_root / "datasets").mkdir(parents=True)
    (config_root / "experiments").mkdir(parents=True)

    reads = tmp_path / "reads.fastq"
    with reads.open("wb") as fh:
        fh.seek(100 * 1000 * 1000 * 1000)
        fh.write(b"\0")

    (config_root / "datasets" / "ds1.yaml").write_text(
        "\n".join(
            [
                "reads:",
                "  - " + str(reads),
            ]
        )
        + "\n"
    )
    (config_root / "experiments" / "taxor.yaml").write_text(
        "\n".join(
            [
                "tool: taxor",
                "db: /db/cami_refseq",
                "datasets:",
                "  - ds1",
            ]
        )
        + "\n"
    )

    calls = []

    class FakeTool:
        name = "taxor"

        def __init__(self, _config):
            pass

    class FakeRunner:
        def __init__(self, _runs_root, _profile_root=None):
            pass

        def run(self, *, exp, dataset, tool, executor):
            calls.append((exp["name"], dataset["name"], tool.name, callable(executor)))

    monkeypatch.setattr(cli_mod, "Runner", FakeRunner)
    monkeypatch.setattr(cli_mod.TOOLS, "get", lambda _name: FakeTool)

    args = Namespace(
        exp="taxor",
        config=str(config_root),
        runs=str(tmp_path / "results" / "classify"),
        profile=str(tmp_path / "results" / "profile"),
        chimera_bin="Chimera",
        ganon_bin="ganon",
        ganon_env="ganon",
        sylph_bin="sylph",
        sylph_env="sylph",
        dry_run=False,
        dataset=[],
    )

    run_cmd(args)

    assert calls == [("taxor", "ds1", "taxor", True)]


def test_run_cmd_expands_dataset_collection_samples(tmp_path: Path, monkeypatch):
    config_root = tmp_path / "configs"
    (config_root / "datasets").mkdir(parents=True)
    (config_root / "experiments").mkdir(parents=True)

    for name in ("a_R1.fq", "a_R2.fq", "b_R1.fq", "b_R2.fq"):
        (tmp_path / name).write_text("@r\nA\n+\nI\n")
    (config_root / "datasets" / "collection.yaml").write_text(
        "\n".join(
            [
                "name: collection",
                "group: real",
                "samples:",
                "  - sample_id: sampleA",
                "    paired:",
                f"      - {tmp_path / 'a_R1.fq'}",
                f"      - {tmp_path / 'a_R2.fq'}",
                "  - sample_id: sampleB",
                "    paired:",
                f"      - {tmp_path / 'b_R1.fq'}",
                f"      - {tmp_path / 'b_R2.fq'}",
            ]
        )
        + "\n"
    )
    (config_root / "experiments" / "chimera.yaml").write_text(
        "\n".join(
            [
                "tool: chimera",
                "db: /db/cami_refseq",
                "datasets:",
                "  - collection",
            ]
        )
        + "\n"
    )

    calls = []

    class FakeTool:
        name = "chimera"

        def __init__(self, _config):
            pass

    class FakeRunner:
        def __init__(self, _runs_root, _profile_root=None):
            pass

        def run(self, *, exp, dataset, tool, executor):
            calls.append(
                (
                    dataset["name"],
                    dataset["dataset_collection"],
                    dataset["sample_id"],
                    tuple(Path(p).name for p in dataset["paired"]),
                )
            )

    monkeypatch.setattr(cli_mod, "Runner", FakeRunner)
    monkeypatch.setattr(cli_mod.TOOLS, "get", lambda _name: FakeTool)

    args = Namespace(
        exp="chimera",
        config=str(config_root),
        runs=str(tmp_path / "results" / "classify"),
        profile=str(tmp_path / "results" / "profile"),
        chimera_bin="Chimera",
        ganon_bin="ganon",
        ganon_env="ganon",
        sylph_bin="sylph",
        sylph_env="sylph",
        dry_run=False,
        dataset=["collection"],
    )

    run_cmd(args)

    assert calls == [
        ("collection.sampleA", "collection", "sampleA", ("a_R1.fq", "a_R2.fq")),
        ("collection.sampleB", "collection", "sampleB", ("b_R1.fq", "b_R2.fq")),
    ]


def test_run_cmd_selects_dataset_by_display_dataset(tmp_path: Path, monkeypatch):
    config_root = tmp_path / "configs"
    (config_root / "datasets").mkdir(parents=True)
    (config_root / "experiments").mkdir(parents=True)

    for name in ("b0_R1.fq", "b0_R2.fq", "b1_R1.fq", "b1_R2.fq"):
        (tmp_path / name).write_text("@r\nA\n+\nI\n")
    (config_root / "datasets" / "batches.yaml").write_text(
        "\n".join(
            [
                "name: ds-ganon-batches",
                "display_dataset: ds",
                "samples:",
                "  - sample_id: batch00",
                "    paired:",
                f"      - {tmp_path / 'b0_R1.fq'}",
                f"      - {tmp_path / 'b0_R2.fq'}",
                "  - sample_id: batch01",
                "    paired:",
                f"      - {tmp_path / 'b1_R1.fq'}",
                f"      - {tmp_path / 'b1_R2.fq'}",
            ]
        )
        + "\n"
    )
    (config_root / "experiments" / "ganon.yaml").write_text(
        "\n".join(
            [
                "tool: ganon",
                "db: /db/cami_refseq",
                "datasets:",
                "  - ds-ganon-batches",
            ]
        )
        + "\n"
    )

    calls = []

    class FakeTool:
        name = "ganon"

        def __init__(self, _config):
            pass

    class FakeRunner:
        def __init__(self, _runs_root, _profile_root=None):
            pass

        def run(self, *, exp, dataset, tool, executor):
            calls.append((dataset["name"], dataset["display_dataset"], dataset["sample_id"]))
            return {"meta": {"return_code": 0}}

    monkeypatch.setattr(cli_mod, "Runner", FakeRunner)
    monkeypatch.setattr(cli_mod.TOOLS, "get", lambda _name: FakeTool)

    args = Namespace(
        exp="ganon",
        config=str(config_root),
        runs=str(tmp_path / "results" / "classify"),
        profile=str(tmp_path / "results" / "profile"),
        chimera_bin="Chimera",
        ganon_bin="ganon",
        ganon_env="ganon",
        sylph_bin="sylph",
        sylph_env="sylph",
        dry_run=False,
        dataset=["ds"],
    )

    run_cmd(args)

    assert calls == [
        ("ds-ganon-batches.batch00", "ds", "batch00"),
        ("ds-ganon-batches.batch01", "ds", "batch01"),
    ]


def test_run_cmd_returns_nonzero_after_dataset_failures(tmp_path: Path, monkeypatch):
    config_root = tmp_path / "configs"
    (config_root / "datasets").mkdir(parents=True)
    (config_root / "experiments").mkdir(parents=True)

    for name in ("ok.fq", "bad.fq"):
        (tmp_path / name).write_text("@r\nA\n+\nI\n")
    (config_root / "datasets" / "ok.yaml").write_text(f"reads:\n  - {tmp_path / 'ok.fq'}\n")
    (config_root / "datasets" / "bad.yaml").write_text(f"reads:\n  - {tmp_path / 'bad.fq'}\n")
    (config_root / "experiments" / "ganon.yaml").write_text(
        "\n".join(
            [
                "tool: ganon",
                "db: /db/cami_refseq",
                "datasets:",
                "  - ok",
                "  - bad",
            ]
        )
        + "\n"
    )

    calls = []

    class FakeTool:
        name = "ganon"

        def __init__(self, _config):
            pass

    class FakeRunner:
        def __init__(self, _runs_root, _profile_root=None):
            pass

        def run(self, *, exp, dataset, tool, executor):
            calls.append(dataset["name"])
            return {"meta": {"return_code": 137 if dataset["name"] == "bad" else 0}}

    monkeypatch.setattr(cli_mod, "Runner", FakeRunner)
    monkeypatch.setattr(cli_mod.TOOLS, "get", lambda _name: FakeTool)

    args = Namespace(
        exp="ganon",
        config=str(config_root),
        runs=str(tmp_path / "results" / "classify"),
        profile=str(tmp_path / "results" / "profile"),
        chimera_bin="Chimera",
        ganon_bin="ganon",
        ganon_env="ganon",
        sylph_bin="sylph",
        sylph_env="sylph",
        dry_run=False,
        dataset=[],
    )

    with pytest.raises(SystemExit) as exc:
        run_cmd(args)

    assert exc.value.code == 1
    assert calls == ["ok", "bad"]


def test_run_cmd_prepares_strain_madness_long_source(tmp_path: Path, monkeypatch):
    config_root = tmp_path / "configs"
    (config_root / "datasets").mkdir(parents=True)
    (config_root / "experiments").mkdir(parents=True)

    source_root = tmp_path / "strain" / "long_read"
    reads = source_root / "strmgCAMI2_sample_0" / "long_read" / "run0" / "reads" / "anonymous_reads.fq"
    reads.parent.mkdir(parents=True)
    reads.write_text("@S0R0\nACGT\n+\nIIII\n")

    (config_root / "datasets" / "strain.yaml").write_text(
        "\n".join(
            [
                "name: strain",
                "sample_ids: [0]",
                "source:",
                "  kind: strain_madness",
                "  read_type: long",
                f"  root: {source_root}",
                f"truth_dir: {source_root}",
            ]
        )
        + "\n"
    )
    (config_root / "experiments" / "chimera.yaml").write_text(
        "\n".join(
            [
                "tool: chimera",
                "db: /db/cami_refseq",
                "datasets:",
                "  - strain",
            ]
        )
        + "\n"
    )

    calls = []

    class FakeTool:
        name = "chimera"

        def __init__(self, _config):
            pass

    class FakeRunner:
        def __init__(self, _runs_root, _profile_root=None):
            pass

        def run(self, *, exp, dataset, tool, executor):
            calls.append((dataset["name"], dataset["sample_ids"], dataset["reads"]))

    monkeypatch.setattr(cli_mod, "Runner", FakeRunner)
    monkeypatch.setattr(cli_mod.TOOLS, "get", lambda _name: FakeTool)

    args = Namespace(
        exp="chimera",
        config=str(config_root),
        runs=str(tmp_path / "results" / "classify"),
        profile=str(tmp_path / "results" / "profile"),
        chimera_bin="Chimera",
        ganon_bin="ganon",
        ganon_env="ganon",
        sylph_bin="sylph",
        sylph_env="sylph",
        dry_run=False,
        dataset=[],
    )

    run_cmd(args)

    assert calls == [("strain", [0], [str(reads)])]


def test_prepare_strain_madness_short_source_writes_paired_fastqs(tmp_path: Path):
    source_root = tmp_path / "strain" / "short_read"
    interleaved = source_root / "strmgCAMI2_sample_0" / "short_read" / "run0" / "reads" / "anonymous_reads.fq"
    interleaved.parent.mkdir(parents=True)
    interleaved.write_text(
        "@S0R0/1\nACGT\n+\nIIII\n"
        "@S0R0/2\nTGCA\n+\nJJJJ\n"
    )

    r1 = tmp_path / "prepared" / "sample_0_R1.fq"
    r2 = tmp_path / "prepared" / "sample_0_R2.fq"
    dataset = {
        "name": "strain-short",
        "sample_ids": [0],
        "paired": [str(r1), str(r2)],
        "source": {
            "kind": "strain_madness",
            "read_type": "short",
            "root": str(source_root),
            "sample_ids": [0],
        },
    }

    prepared = prepare_dataset_inputs(dataset)

    assert prepared["paired"] == [str(r1), str(r2)]
    assert r1.read_text() == "@S0R0/1\nACGT\n+\nIIII\n"
    assert r2.read_text() == "@S0R0/2\nTGCA\n+\nJJJJ\n"
    manifest = json.loads((tmp_path / "prepared" / "sample_0_R1.fq.manifest.json").read_text())
    assert manifest["inputs"] == [str(interleaved)]
    assert manifest["records"] == 2
