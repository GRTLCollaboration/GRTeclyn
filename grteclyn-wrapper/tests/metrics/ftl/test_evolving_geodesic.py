"""Tests for 4D evolving null-geodesic FTL probe."""

from __future__ import annotations

import json

import numpy as np

from grteclyn_wrapper.metrics.probes import warpfactory as wf
from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import (
    EvolvingGeodesicFtlReport,
    _emission_times,
    compute_evolving_geodesic_ftl,
    compute_evolving_geodesic_ftl_best_direction,
    compute_evolving_geodesic_ftl_emission_sweep,
    compute_evolving_geodesic_ftl_from_analytic,
    compute_evolving_geodesic_ftl_from_plotfiles,
    integrate_null_ray_on_field,
    read_evolving_geodesic_json,
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
    # allow_frozen_tail: this stack is genuinely stationary, so completing the
    # flight through the clamped final slice is exact (the strict default would
    # correctly reject a real run's stack this short).
    report = compute_evolving_geodesic_ftl(
        field, t_emit=0.0, n_rays=3, allow_frozen_tail=True
    )
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

    # Analytic stack spans only ±0.4 in t; lenient tail is the point here.
    x_report = compute_evolving_geodesic_ftl(
        field, n_rays=5, axis=0, allow_frozen_tail=True
    )
    y_report = compute_evolving_geodesic_ftl(
        field, n_rays=5, axis=1, allow_frozen_tail=True
    )
    assert x_report.f_geo < 0.05
    assert y_report.f_geo > 0.1
    assert y_report.n_reached == y_report.n_rays

    best = compute_evolving_geodesic_ftl_best_direction(
        field, directions=("x", "y", "z"), n_rays=5, allow_frozen_tail=True
    )
    assert best.f_geo > 0.1
    assert any("best_direction=y" in n for n in best.notes)


def test_emission_times_disabled_and_static_return_none():
    g, origin, spacing = _flat_grid()
    stack = np.stack([g, g, g, g, g], axis=0)
    field = EvolvingMetricField(
        g_stack=stack,
        times=np.array([0.0, 1.0, 2.0, 3.0, 4.0]),
        origin=origin,
        spacing=(1.0, *spacing),
    )
    # Disabled by interval<=0 or a single launch.
    assert _emission_times(field, interval=0.0, max_emissions=4) is None
    assert _emission_times(field, interval=1.0, max_emissions=1) is None
    # Non-evolving fields never sweep.
    static = StaticMetricField(g=g, origin=origin, spatial_spacing=spacing)
    assert _emission_times(static, interval=1.0, max_emissions=4) is None
    # Enabled: launches at 0,1,2 (capped at last slice time 4.0).
    assert _emission_times(field, interval=1.0, max_emissions=3) == [0.0, 1.0, 2.0]
    # max_emissions beyond the slice window is clamped to available times.
    assert _emission_times(field, interval=2.0, max_emissions=9) == [0.0, 2.0, 4.0]


def test_emission_sweep_maps_launch_times_and_reports_peak():
    g, origin, spacing = _flat_grid()
    stack = np.stack([g] * 5, axis=0)
    field = EvolvingMetricField(
        g_stack=stack,
        times=np.array([0.0, 1.0, 2.0, 3.0, 4.0]),
        origin=origin,
        spacing=(1.0, *spacing),
    )
    emit_times = _emission_times(field, interval=1.0, max_emissions=3)
    report = compute_evolving_geodesic_ftl_emission_sweep(
        field, directions=("x",), emit_times=emit_times, n_rays=3
    )
    # One sweep entry per launch time, in order.
    assert len(report.emit_sweep) == 3
    assert [row[0] for row in report.emit_sweep] == [0.0, 1.0, 2.0]
    # Headline f_geo is the peak over launch times.
    assert report.f_geo == max(row[1] for row in report.emit_sweep)
    assert any(n.startswith("emit_sweep:") for n in report.notes)
    assert any("peak f_geo=" in n for n in report.notes)


def test_emission_sweep_captures_time_dependent_peak():
    """On a moving Alcubierre bubble the best launch time need not be t=0."""
    g, spacing = wf.alcubierre_metric(
        velocity=0.5, bubble_radius=2.0, sigma=2.0, half_width=4.0, n_space=20, dt=0.2
    )
    field = evolving_field_from_analytic_stack(g, spacing)
    emit_times = _emission_times(field, interval=float(field.times[1] - field.times[0]),
                                 max_emissions=3)
    if emit_times is None:  # single-slice stack -> nothing to sweep
        return
    report = compute_evolving_geodesic_ftl_emission_sweep(
        field, directions=("x",), emit_times=emit_times, n_rays=5
    )
    assert len(report.emit_sweep) == len(emit_times)
    # The reported f_geo equals the best entry in the map.
    assert report.f_geo == max(row[1] for row in report.emit_sweep)


def test_emission_sweep_json_roundtrip(tmp_path):
    report = EvolvingGeodesicFtlReport(
        f_geo=0.12,
        f_geo_frozen_peak=0.2,
        t_emit=4.0,
        t_arrival=18.0,
        t_flat=14.0,
        n_rays=5,
        n_reached=5,
        max_h_drift=0.001,
        h_quality_ok=True,
        max_h_rel_drift=0.001,
        notes=("emit_sweep: t=0.00->f=0.050(n5), t=4.00->f=0.120(n5)",),
        emit_sweep=((0.0, 0.05, 5), (4.0, 0.12, 5)),
    )
    path = tmp_path / "evolving_geodesic.json"
    write_evolving_geodesic_json(path, report)
    loaded = read_evolving_geodesic_json(path)
    assert loaded is not None
    assert loaded.emit_sweep == ((0.0, 0.05, 5), (4.0, 0.12, 5))
    assert loaded.f_geo == 0.12


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
