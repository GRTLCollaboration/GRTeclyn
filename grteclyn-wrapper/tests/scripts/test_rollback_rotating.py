"""Tests for rotating-wormhole rollback prefix support."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path


def _rollback_script() -> Path:
    return (
        Path(__file__).resolve().parents[2]
        / "scripts"
        / "wormhole"
        / "rollback"
    )


def test_rollback_recognizes_rotating_wormhole_prefixes(tmp_path: Path) -> None:
    data_dir = tmp_path / "output"
    data_dir.mkdir()
    (data_dir / "RotatingWormholeChk00040").mkdir()
    (data_dir / "RotatingWormholeChk00040" / "Header").write_text(
        "0\n1\n0.2\n", encoding="utf-8"
    )
    (data_dir / "RotatingWormholeChk00080").mkdir()
    (data_dir / "RotatingWormholePlt00050").mkdir()
    (data_dir / "RotatingWormholePlt00090").mkdir()
    (data_dir / "data").mkdir()
    (data_dir / "data" / "constraint_norms.dat").write_text(
        "0.0 0.1\n" * 100, encoding="utf-8"
    )
    (data_dir / "data" / "collapse_diagnostics.dat").write_text(
        "0.0 0.1\n" * 100, encoding="utf-8"
    )
    (data_dir / "small_data").mkdir()
    (data_dir / "small_data" / "consume_state.json").write_text(
        '{"RotatingWormholePlt00050": {}, "RotatingWormholePlt00090": {}}\n',
        encoding="utf-8",
    )

    completed = subprocess.run(
        [
            str(_rollback_script()),
            "--step",
            "40",
            "--data",
            str(data_dir),
            "--dry-run",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    assert "Removed checkpoint dirs (> 40):   1" in completed.stdout
    assert "Removed plotfile dirs (> 40):     2" in completed.stdout
    assert "Trimmed consume_state.json keys:      2" in completed.stdout


def test_consume_state_pattern_includes_rotating_prefix() -> None:
    pattern = re.compile(
        r"(?:SupportedWormholePlt|RotatingWormholePlt|WormholePlt|plt)(\d+)"
    )
    assert pattern.fullmatch("RotatingWormholePlt00040") is not None
    assert int(pattern.fullmatch("RotatingWormholePlt00040").group(1)) == 40
