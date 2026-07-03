"""Unit tests for the directional Psi4 metric reader."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from grteclyn_wrapper.metrics.diagnostics.psi4 import read_psi4_metrics
from grteclyn_wrapper.metrics.types.psi4 import Psi4Metrics


def _write_psi4_directional(tmp_path: Path, rows: list[tuple[float, float, float, float]]) -> Path:
    path = tmp_path / "psi4_directional.dat"
    lines = ["# time  P_total  P_z_beam  beam_ratio"]
    for t, p_total, p_z, ratio in rows:
        lines.append(f"{t:.6e}  {p_total:.6e}  {p_z:.6e}  {ratio:.6e}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_read_psi4_metrics_basic(tmp_path: Path) -> None:
    rows = [
        (0.0, 1.0, 0.1, 0.1),
        (1.0, 2.0, 0.5, 0.25),
        (2.0, 3.0, 1.5, 0.5),
    ]
    path = _write_psi4_directional(tmp_path, rows)
    metrics = read_psi4_metrics(tmp_path)
    assert metrics is not None
    assert metrics.n_samples == 3
    np.testing.assert_allclose(metrics.mean_total_power, 2.0)
    np.testing.assert_allclose(metrics.mean_z_beam_power, 0.7, rtol=1e-3)
    np.testing.assert_allclose(metrics.peak_total_power, 3.0)
    np.testing.assert_allclose(metrics.peak_z_beam_power, 1.5)
    np.testing.assert_allclose(metrics.final_beam_ratio, 0.5)


def test_read_psi4_metrics_missing_file(tmp_path: Path) -> None:
    assert read_psi4_metrics(tmp_path) is None


def test_read_psi4_metrics_empty_file(tmp_path: Path) -> None:
    (tmp_path / "psi4_directional.dat").write_text("# header only\n", encoding="utf-8")
    assert read_psi4_metrics(tmp_path) is None


def test_score_signal(tmp_path: Path) -> None:
    rows = [
        (0.0, 1.0, 0.5, 0.5),
        (1.0, 1.0, 0.5, 0.5),
    ]
    _write_psi4_directional(tmp_path, rows)
    metrics = read_psi4_metrics(tmp_path)
    assert metrics is not None
    np.testing.assert_allclose(metrics.score_signal, 1.5)
