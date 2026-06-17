"""Tests for 4D evolving null-geodesic FTL probe."""

from __future__ import annotations

import json

import numpy as np

from grteclyn_wrapper.metrics.probes import warpfactory as wf
from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import (
    EvolvingGeodesicFtlReport,
    compute_evolving_geodesic_ftl,
    compute_evolving_geodesic_ftl_best_direction,
    compute_evolving_geodesic_ftl_from_analytic,
    compute_evolving_geodesic_ftl_from_plotfiles,
    integrate_null_ray_on_field,
    write_evolving_geodesic_json,
)
from grteclyn_wrapper.metrics.probes.ftl.geodesic import (
    integrate_null_ray,
    partial_inverse_metric,
)
from grteclyn_wrapper.metrics.probes.ftl.metric_field import (
    EvolvingMetricField,
    StaticMetricField,
    evolving_field_from_analytic_stack,
)


def _flat_grid(n: int = 17):
    g = np.zeros((n, n, n, 4, 4))
    g[..., 0, 0] = -1.0
    for i in range(1, 4):
        g[..., i, i] = 1.0
    spacing = (0.5, 0.5, 0.5)
    origin = np.array([-4.0, -4.0, -4.0])
    return g, origin, spacing


def test_evolving_matches_frozen_on_stationary_stack():
    g, origin, spacing = _flat_grid()
    dg_inv = partial_inverse_metric(g, spacing)
    frozen = integrate_null_ray(
        g, dg_inv, origin, spacing, x_start=-3.5, x_end=3.5, y0=0.0, z0=0.0
    )
    stack = np.stack([g, g, g], axis=0)
    field = EvolvingMetricField(
        g_stack=stack,
        times=np.array([0.0, 1.0, 2.0]),
        origin=origin,
        spacing=(1.0, *spacing),
    )
    report = compute_evolving_geodesic_ftl(field, t_emit=0.0, n_rays=3)
    assert frozen.reached
    assert report.n_reached > 0
    rel = abs(report.f_geo - (frozen.t_flat - frozen.t_coord) / frozen.t_flat) / frozen.t_flat
    assert rel < 0.05


def test_integrate_null_ray_on_field_flat_reaches_detector():
    g, origin, spacing = _flat_grid()
    field = StaticMetricField(g=g, origin=origin, spatial_spacing=spacing)
    res = integrate_null_ray_on_field(
        field, x_start=-3.5, x_end=3.5, y0=0.0, z0=0.0, t0=0.0
    )
    assert res.reached
    assert res.t_coord is not None
    assert abs(res.t_coord - res.t_flat) / res.t_flat < 0.05


def test_alcubierre_evolving_finds_shortcut():
    g, spacing = wf.alcubierre_metric(
        velocity=0.5, bubble_radius=2.0, sigma=2.0, half_width=4.0, n_space=20, dt=0.2
    )
    report = compute_evolving_geodesic_ftl_from_analytic(g, spacing, n_rays=5)
    assert report.f_geo > 0.1
    assert report.n_reached == report.n_rays


def _swap_spatial_xy(
    g: np.ndarray, spacing: tuple[float, float, float, float]
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    """Permute an x-aligned metric grid to y-aligned (swap spatial x <-> y)."""
    perm4 = [0, 2, 1, 3]
    g_y = g.transpose(0, 2, 1, 3, 4, 5)
    g_y = g_y[..., perm4, :][..., :, perm4]
    dt, dx, dy, dz = spacing
    return g_y, (dt, dy, dx, dz)


def test_multi_direction_probe_finds_y_aligned_shortcut():
    g, spacing = wf.alcubierre_metric(
        velocity=0.5, bubble_radius=2.0, sigma=2.0, half_width=4.0, n_space=20, dt=0.2
    )
    g_y, spacing_y = _swap_spatial_xy(g, spacing)
    field = evolving_field_from_analytic_stack(g_y, spacing_y)

    x_report = compute_evolving_geodesic_ftl(field, n_rays=5, axis=0)
    y_report = compute_evolving_geodesic_ftl(field, n_rays=5, axis=1)
    assert x_report.f_geo < 0.05
    assert y_report.f_geo > 0.1
    assert y_report.n_reached == y_report.n_rays

    best = compute_evolving_geodesic_ftl_best_direction(
        field, directions=("x", "y", "z"), n_rays=5
    )
    assert best.f_geo > 0.1
    assert any("best_direction=y" in n for n in best.notes)


def test_too_few_plotfiles_returns_none():
    assert compute_evolving_geodesic_ftl_from_plotfiles(["a", "b"]) is None


def test_write_evolving_geodesic_json_handles_numpy_scalars(tmp_path):
    report = EvolvingGeodesicFtlReport(
        f_geo=np.float64(0.014),
        f_geo_frozen_peak=np.float64(0.058),
        t_emit=np.float64(0.0),
        t_arrival=np.float64(14.2),
        t_flat=np.float64(14.4),
        n_rays=5,
        n_reached=5,
        max_h_drift=np.float64(0.002),
        h_quality_ok=np.bool_(True),
        max_h_rel_drift=np.float64(0.002),
        notes=("numpy scalars",),
    )
    path = tmp_path / "evolving_geodesic.json"
    write_evolving_geodesic_json(path, report)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["f_geo"] == 0.014
    assert payload["h_quality_ok"] is True


def test_patch_ftl_timeseries_evolving_columns(tmp_path):
    from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import patch_ftl_timeseries_evolving

    path = tmp_path / "ftl_timeseries.dat"
    path.write_text(
        "# time  f_op  f_geo  geo_trustworthy  max_local_speed  "
        "superluminal_fraction  max_shift  structure_coherence  reachable  "
        "n_rays  n_reached  max_h_rel_drift\n"
        "0.0  0.1  0.05  1  1.1  0.2  0.0  nan  1  5  5  0.001\n",
        encoding="utf-8",
    )
    patch_ftl_timeseries_evolving(path, f_geo_evol=0.12, f_geo_evol_ok=True)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert "f_geo_evol" in lines[0]
    parts = lines[1].split()
    assert len(parts) == 14
    assert float(parts[12]) == 0.12
    assert parts[13] == "1"
