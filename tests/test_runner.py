from pathlib import Path

from chimera_bench.core.runner import Runner


class FakeTool:
    name = "chimera"

    def build_cmd(self, *, dataset, exp, out_prefix):
        return ["echo", "ok"], {"classify_tsv": f"{out_prefix}.tsv"}


def test_runner_creates_run_dir_and_meta(tmp_path: Path):
    runner = Runner(tmp_path / "runs")
    dataset = {"reads": ["/r1.fq"]}
    exp = {"name": "exp1", "db": "/db", "threads": 1}

    def exec_stub(cmd, cwd, stdout_path, stderr_path):
        Path(stdout_path).write_text("ok")
        Path(stderr_path).write_text("")
        out_tsv = Path(cwd) / "outputs" / "ChimeraClassify.tsv"
        out_tsv.write_text("r1\t123:1\n")
        return 0

    result = runner.run(exp=exp, dataset=dataset, tool=FakeTool(), executor=exec_stub)

    meta = Path(result["run_dir"]) / "meta.json"
    assert meta.exists()
