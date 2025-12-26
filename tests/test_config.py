from pathlib import Path
import yaml

from chimera_bench import config


def test_load_yaml_dir_uses_name_or_filename(tmp_path: Path):
    d = tmp_path / "datasets"
    d.mkdir()
    (d / "cami.yaml").write_text(
        """
name: cami-test
reads:
  - /tmp/a.fq
"""
    )
    (d / "atcc.yaml").write_text(
        """
reads:
  - /tmp/b.fq
"""
    )

    out = config.load_yaml_dir(d)

    assert "cami-test" in out
    assert "atcc" in out
    assert out["cami-test"]["reads"] == ["/tmp/a.fq"]
