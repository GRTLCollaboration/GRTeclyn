"""Tests for central metrics collector wiring."""

from __future__ import annotations

from pathlib import Path

from grteclyn_wrapper.metrics import read_episode_metrics
from grteclyn_wrapper.metrics.catalog import list_metrics


def _write_minimal_episode(root: Path, *, central_rows: str | None) -> None:
    data = root / "data"
    small = root / "small_data"
    data.mkdir(parents=True)
    small.mkdir(parents=True)
    with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as fh:
        fh.write("0.0  1.0  0.98  0.02  0  0  0  5.0  0.05  5.0  0  0  0  0\n")
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as fh:
        fh.write("0.0  1e-4  1e-4  0  0  0\n")
    if central_rows is not None:
        header = "# time  rho_req  lapse  scalar_activity  phi_re  phi_im"
        (small / "central_timeseries.dat").write_text(
            header + "\n" + central_rows + "\n",
            encoding="utf-8",
        )


def test_read_episode_metrics_loads_central_when_present(tmp_path: Path) -> None:
    ep = tmp_path / "ep"
    _write_minimal_episode(
        ep,
        central_rows="0.0  1.0e-4  0.99  0.01  0.0  0.0\n1.0  2.0e-3  0.95  0.02  0.0  0.0",
    )
    metrics = read_episode_metrics(ep)
    assert metrics.central is not None
    assert metrics.central.peak_rho_req_at_origin > 0.0


def test_read_episode_metrics_central_none_when_absent(tmp_path: Path) -> None:
    ep = tmp_path / "ep"
    _write_minimal_episode(ep, central_rows=None)
    metrics = read_episode_metrics(ep)
    assert metrics.central is None


def test_catalog_lists_central_group() -> None:
    groups = {spec.group for spec in list_metrics()}
    assert "central" in groups
