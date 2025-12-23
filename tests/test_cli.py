from pathlib import Path

from chimera_bench.cli import main


def test_cli_run_dry(tmp_path, monkeypatch):
    cfg = tmp_path / "configs"
    (cfg / "datasets").mkdir(parents=True)
    (cfg / "experiments").mkdir(parents=True)
    (cfg / "datasets" / "d.yaml").write_text("name: d\nreads: [/r.fq]\n")
    (cfg / "experiments" / "e.yaml").write_text(
        "name: e\ntool: chimera\ndatasets: [d]\ndb: /db\n"
    )

    out_root = tmp_path / "runs"
    argv = [
        "chimera-bench",
        "run",
        "--exp",
        "e",
        "--config",
        str(cfg),
        "--runs",
        str(out_root),
        "--dry-run",
    ]
    monkeypatch.setattr("sys.argv", argv)

    main()
    assert out_root.exists()
