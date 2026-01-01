from pathlib import Path
import json

from chimera_bench.core.build_runner import BuildRunner


class FakeBuildTool:
    name = "ganon"

    def build_db_steps(self, *, build, out_dir):
        return [{"name": "build_db", "cmd": ["echo", "ok"], "outputs": {"db_prefix": build.get("db_prefix")}}]


def test_build_runner_uses_tool_db_layout(tmp_path: Path):
    runner = BuildRunner(tmp_path / "builds")
    build = {"name": "cami", "db_prefix": "DB/cami_refseq"}

    def exec_stub(cmd, cwd, stdout_path, stderr_path, resource_path):
        Path(stdout_path).write_text("ok")
        Path(stderr_path).write_text("")
        Path(resource_path).write_text("")
        return 0

    result = runner.run(build=build, tool=FakeBuildTool(), executor=exec_stub)
    run_dir = Path(result["run_dir"])
    assert run_dir.parts[-2] == "ganon"
    assert run_dir.parts[-1] == "cami_refseq"


def test_build_runner_writes_time_log_under_db_dir(tmp_path: Path):
    runner = BuildRunner(tmp_path / "builds")
    build = {"name": "cami", "db_prefix": "DB/cami_refseq"}

    seen = {}

    def exec_stub(cmd, cwd, stdout_path, stderr_path, resource_path):
        seen["resource_path"] = resource_path
        Path(stdout_path).write_text("ok")
        Path(stderr_path).write_text("")
        Path(resource_path).parent.mkdir(parents=True, exist_ok=True)
        Path(resource_path).write_text("Maximum resident set size (kbytes): 123\nExit status: 0\n")
        return 0

    result = runner.run(build=build, tool=FakeBuildTool(), executor=exec_stub)
    run_dir = Path(result["run_dir"])
    resource_path = Path(seen["resource_path"])
    assert resource_path.parent == run_dir / "DB"

    meta = json.loads((run_dir / "meta.json").read_text())
    assert meta["resource"]["max_rss_kb"] == 123
