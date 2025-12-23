from pathlib import Path

from chimera_bench.core.resources import aggregate_resources, parse_time_log


def test_parse_time_log(tmp_path: Path):
    log = tmp_path / "time.log"
    log.write_text(
        "\n".join(
            [
                "User time (seconds): 1.23",
                "System time (seconds): 0.45",
                "Percent of CPU this job got: 99%",
                "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.70",
                "Maximum resident set size (kbytes): 123456",
                "Exit status: 0",
            ]
        )
        + "\n"
    )
    parsed = parse_time_log(log)
    assert parsed["user_time_seconds"] == 1.23
    assert parsed["system_time_seconds"] == 0.45
    assert parsed["cpu_percent"] == 99.0
    assert parsed["elapsed_wall_seconds"] == 1.70
    assert parsed["max_rss_kb"] == 123456
    assert parsed["exit_status"] == 0


def test_aggregate_resources():
    steps = [
        {"resource": {"max_rss_kb": 1000, "user_time_seconds": 1.0, "system_time_seconds": 0.5}},
        {"resource": {"max_rss_kb": 2000, "user_time_seconds": 2.5, "system_time_seconds": 0.25}},
        {"resource": {}},
    ]
    agg = aggregate_resources(steps)
    assert agg["max_rss_kb"] == 2000
    assert agg["user_time_seconds"] == 3.5
    assert agg["system_time_seconds"] == 0.75
