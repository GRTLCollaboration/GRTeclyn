"""Tests for central_radial_profile.dat parsing."""

from __future__ import annotations

from pathlib import Path

from grteclyn_wrapper.metrics.diagnostics.central_radial import read_central_radial_profile


def test_read_radial_profile_derived_metrics(tmp_path: Path) -> None:
    path = tmp_path / "central_radial_profile.dat"
    path.write_text(
        "\n".join(
            [
                "# time=0.0000000000000000e+00",
                "# r  rho_req  lapse  scalar_activity",
                "5.0000000000000000e-01  1.0000000000000000e-02  9.9000000000000000e-01  1.0000000000000000e-02",
                "3.0000000000000000e+00  5.0000000000000000e-03  9.8000000000000000e-01  5.0000000000000000e-03",
                "# time=2.0000000000000000e+00",
                "# r  rho_req  lapse  scalar_activity",
                "5.0000000000000000e-01  2.0000000000000000e-02  8.0000000000000000e-01  2.0000000000000000e-02",
                "3.0000000000000000e+00  4.0000000000000000e-03  9.5000000000000000e-01  4.0000000000000000e-03",
            ]
        ),
        encoding="utf-8",
    )
    metrics = read_central_radial_profile(path, dx_finest=0.25, ring_radius=3.0)
    assert metrics is not None
    assert metrics.peak_radius == 0.5
    assert metrics.compression_ratio > 1.0
    assert metrics.initial_rho_at_ring == 0.005
