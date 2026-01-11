from pathlib import Path

from chimera_bench.cli import DEFAULT_THREADS, _executor, main


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


def test_cli_build_dry(tmp_path, monkeypatch):
    cfg = tmp_path / "configs"
    (cfg / "build").mkdir(parents=True)
    (cfg / "build" / "b.yaml").write_text(
        "\n".join(
            [
                "name: b",
                "tool: ganon",
                "db_prefix: /db/prefix",
                "threads: 1",
                "build:",
                "  mode: custom",
                "  input: [/ref.fa]",
            ]
        )
        + "\n"
    )

    out_root = tmp_path / "builds"
    argv = [
        "chimera-bench",
        "build",
        "--build",
        "b",
        "--config",
        str(cfg),
        "--runs",
        str(out_root),
        "--dry-run",
    ]
    monkeypatch.setattr("sys.argv", argv)

    main()
    assert out_root.exists()


def test_default_threads_is_32():
    assert DEFAULT_THREADS == 32


def test_executor_creates_parent_dirs_for_logs(tmp_path: Path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()

    stdout_path = tmp_path / "nested" / "logs" / "stdout.log"
    stderr_path = tmp_path / "nested" / "logs" / "stderr.log"
    resource_path = tmp_path / "nested" / "logs" / "time.log"

    rc = _executor(
        ["bash", "-lc", "echo hello"],
        cwd=run_dir,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        resource_path=resource_path,
    )

    assert rc == 0
    assert stdout_path.exists()
    assert "hello" in stdout_path.read_text(encoding="utf-8", errors="ignore")
    assert stderr_path.exists()
    assert resource_path.exists()
