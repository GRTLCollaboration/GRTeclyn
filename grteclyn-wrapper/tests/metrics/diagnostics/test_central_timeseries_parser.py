"""Tests for central_timeseries.dat parsing."""

from __future__ import annotations

from pathlib import Path

import pytest

from grteclyn_wrapper.metrics.diagnostics.central import read_central_field_metrics
from grteclyn_wrapper.metrics.types.central import CentralFieldMetrics


def _write_central(path: Path, rows: list[str]) -> None:
    header = "# time  rho_req  lapse  scalar_activity  phi_re  phi_im"
    path.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")


def test_parse_valid_central_timeseries(tmp_path: Path) -> None:
    path = tmp_path / "central_timeseries.dat"
    _write_central(
        path,
        [
            "0.0  1.0e-4  0.99  0.01  0.0  0.0",
            "1.0  2.0e-3  0.95  0.02  0.0  0.0",
            "2.0  5.0e-3  0.90  0.03  0.0  0.0",
            "3.0  1.0e-3  0.92  0.01  0.0  0.0",
        ],
    )
    metrics = read_central_field_metrics(path)
    assert metrics is not None
    assert metrics.n_frames == 4
    assert metrics.peak_rho_req_at_origin == pytest.approx(5.0e-3)
    assert metrics.peak_rho_req_time == pytest.approx(2.0)
    assert metrics.initial_rho_req_at_origin == pytest.approx(1.0e-4)
    assert metrics.min_lapse_at_origin == pytest.approx(0.90)
    assert 0.0 <= metrics.wave_chromaticity <= 1.0


def test_missing_file_returns_none(tmp_path: Path) -> None:
    assert read_central_field_metrics(tmp_path / "missing.dat") is None


def test_empty_file_returns_none(tmp_path: Path) -> None:
    path = tmp_path / "central_timeseries.dat"
    path.write_text("# time  rho_req  lapse  scalar_activity  phi_re  phi_im\n", encoding="utf-8")
    assert read_central_field_metrics(path) is None


def test_single_row_edge_case(tmp_path: Path) -> None:
    path = tmp_path / "central_timeseries.dat"
    _write_central(path, ["1.5  3.0e-3  0.88  0.05  0.0  0.0"])
    metrics = read_central_field_metrics(path)
    assert metrics is not None
    assert metrics.peak_rho_req_at_origin == pytest.approx(3.0e-3)
    assert metrics.wave_chromaticity == 0.0


def test_legacy_file_without_geometric_columns(tmp_path: Path) -> None:
    """Legacy 6-column files parse with no geometric data (chi_drop=0)."""
    path = tmp_path / "central_timeseries.dat"
    _write_central(
        path,
        [
            "0.0  1.0e-4  0.99  0.01  0.0  0.0",
            "1.0  2.0e-3  0.95  0.02  0.0  0.0",
        ],
    )
    metrics = read_central_field_metrics(path)
    assert metrics is not None
    assert metrics.has_geometric_data is False
    assert metrics.chi_drop == 0.0
    assert metrics.peak_abs_weyl4 == 0.0
    assert metrics.peak_abs_K == 0.0


def test_parse_geometric_columns(tmp_path: Path) -> None:
    """chi, K, weyl4 (cols 10,11,12) are parsed into geometric properties."""
    path = tmp_path / "central_timeseries.dat"
    header = (
        "# time  rho_req  lapse  scalar_activity  phi_re  phi_im  "
        "noether_charge  phase_coherence  ham_abs  mom_abs  chi  K  weyl4"
    )
    rows = [
        # t  rho  lapse act  phi_re phi_im noeth phase ham  mom  chi   K     weyl4
        "0.0  1e-4  0.99  0.01  0 0  0 0  0 0  1.00  0.0   0.0",
        "1.0  2e-3  0.95  0.02  0 0  0 0  0 0  0.70  0.5   3e-3",
        "2.0  5e-3  0.80  0.03  0 0  0 0  0 0  0.30  1.2  -1.5e-2",
        "3.0  1e-3  0.85  0.01  0 0  0 0  0 0  0.45  0.8   5e-3",
    ]
    path.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
    metrics = read_central_field_metrics(path)
    assert metrics is not None
    assert metrics.has_geometric_data is True
    # min chi = 0.30 -> chi_drop = 0.70
    assert metrics.min_chi_at_origin == pytest.approx(0.30)
    assert metrics.chi_drop == pytest.approx(0.70)
    # peak |K| = 1.2
    assert metrics.peak_abs_K == pytest.approx(1.2)
    # peak |weyl4| = 1.5e-2, occurs at t=2.0
    assert metrics.peak_abs_weyl4 == pytest.approx(1.5e-2)
    assert metrics.weyl4_peak_time == pytest.approx(2.0)


def test_focusing_efficiency_uses_baseline() -> None:
    metrics = CentralFieldMetrics(
        n_frames=2,
        t=(0.0, 1.0),
        rho_req=(1.0e-8, 1.0e-2),
        lapse=(1.0, 0.9),
        scalar_activity=(0.01, 0.02),
        peak_rho_req_at_origin=1.0e-2,
        peak_rho_req_time=1.0,
        initial_rho_req_at_origin=1.0e-8,
        min_lapse_at_origin=0.9,
        wave_chromaticity=0.0,
        initial_rho_baseline=1.0e-2,
    )
    assert metrics.focusing_efficiency == pytest.approx(1.0)
