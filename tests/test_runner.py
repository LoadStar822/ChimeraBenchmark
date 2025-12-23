from pathlib import Path

from chimera_bench.core.runner import Runner


class FakeTool:
    name = "chimera"
    output_basename = "ChimeraClassify"

    def build_cmd(self, *, dataset, exp, out_prefix):
        return ["echo", "ok"], {"classify_tsv": f"{out_prefix}.tsv"}


def test_runner_creates_run_dir_and_meta(tmp_path: Path):
    runner = Runner(tmp_path / "runs")
    dataset = {"reads": ["/r1.fq"]}
    exp = {"name": "exp1", "db": "/db", "threads": 1}

    def exec_stub(cmd, cwd, stdout_path, stderr_path, resource_path):
        Path(stdout_path).write_text("ok")
        Path(stderr_path).write_text("")
        Path(resource_path).write_text("User time (seconds): 0.10\nSystem time (seconds): 0.05\n")
        out_tsv = Path(cwd) / "outputs" / "ChimeraClassify.tsv"
        out_tsv.write_text("r1\t123:1\n")
        return 0

    result = runner.run(exp=exp, dataset=dataset, tool=FakeTool(), executor=exec_stub)

    meta = Path(result["run_dir"]) / "meta.json"
    assert meta.exists()


def test_runner_supports_build_steps(tmp_path: Path):
    runner = Runner(tmp_path / "runs")
    dataset = {"reads": ["/r1.fq"]}
    exp = {"name": "exp2", "db": "/db", "threads": 1}
    calls = []

    class StepTool:
        name = "ganon"

        def build_steps(self, *, dataset, exp, out_prefix):
            return [
                {"name": "classify", "cmd": ["echo", "ok"], "outputs": {"rep": f"{out_prefix}.rep"}},
                {
                    "name": "report",
                    "cmd": ["echo", "ok"],
                    "outputs": {"report_reads_tre": f"{out_prefix}_reads.tre"},
                },
            ]

    def exec_stub(cmd, cwd, stdout_path, stderr_path, resource_path):
        calls.append((cmd, stdout_path, stderr_path))
        Path(stdout_path).write_text("ok")
        Path(stderr_path).write_text("")
        Path(resource_path).write_text("")
        return 0

    result = runner.run(exp=exp, dataset=dataset, tool=StepTool(), executor=exec_stub)

    assert len(calls) == 2
    run_dir = Path(result["run_dir"])
    assert (run_dir / "logs" / "classify.stdout.log").exists()
    assert (run_dir / "logs" / "report.stdout.log").exists()
