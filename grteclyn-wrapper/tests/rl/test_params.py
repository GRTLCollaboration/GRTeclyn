from __future__ import annotations

from pathlib import Path

from grteclyn_wrapper.rl.params import read_rl_num_lumps


def test_read_rl_num_lumps_from_params(tmp_path: Path) -> None:
    params = tmp_path / "params.txt"
    params.write_text(
        "rl_num_lumps = 3\n# other\nrl_enabled = 1\n",
        encoding="utf-8",
    )
    assert read_rl_num_lumps(params) == 3


def test_read_rl_num_lumps_defaults_to_one(tmp_path: Path) -> None:
    params = tmp_path / "params.txt"
    params.write_text("rl_enabled = 1\n", encoding="utf-8")
    assert read_rl_num_lumps(params) == 1
