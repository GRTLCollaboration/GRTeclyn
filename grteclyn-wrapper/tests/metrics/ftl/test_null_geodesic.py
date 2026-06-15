"""Tests for gauge-invariant null geodesic ray-tracing."""

from __future__ import annotations

import numpy as np

from grteclyn_wrapper.metrics.probes.ftl import geodesic as geo
from grteclyn_wrapper.metrics.probes.ftl.geodesic import (
    GEO_REFINE_N,
    GeodesicFtlReport,
    build_metric_3d_from_plotfile,
    compute_geodesic_ftl_from_plotfile,
    future_null_cov,
    integrate_null_ray,
    inverse_metric_4d,
    null_hamiltonian,
    partial_inverse_metric,
    project_null,
)


def _shifted_metric_grid(n: int, beta_x: float):
    """Constant ADM metric with a pure +x shift ``beta_x`` (lapse 1, flat gamma)."""
    g = np.zeros((n, n, n, 4, 4))
    g[..., 0, 0] = beta_x * beta_x - 1.0  # beta_sq - alpha^2
    g[..., 0, 1] = g[..., 1, 0] = beta_x  # beta lowered with gamma = I
    for i in range(1, 4):
        g[..., i, i] = 1.0
    return g


def test_flat_metric_null_projection():
    g = np.zeros((4, 4, 4, 4, 4))
    g[..., 0, 0] = -1.0
    for i in range(1, 4):
        g[..., i, i] = 1.0
    ginv = inverse_metric_4d(g)
    k = project_null(g[0, 0, 0], ginv[0, 0, 0], g[0, 0, 0] @ np.array([1.0, 1.0, 0.0, 0.0]))
    assert abs(null_hamiltonian(ginv[0, 0, 0], k)) < 1e-10
    assert (ginv[0, 0, 0] @ k)[0] > 0.0


def test_integrate_null_ray_flat_space():
    n = 17
    g = np.zeros((n, n, n, 4, 4))
    g[..., 0, 0] = -1.0
    for i in range(1, 4):
        g[..., i, i] = 1.0
    spacing = (0.5, 0.5, 0.5)
    origin = np.array([-4.0, -4.0, -4.0])
    dg_inv = partial_inverse_metric(g, spacing)
    res = integrate_null_ray(
        g,
        dg_inv,
        origin,
        spacing,
        x_start=-3.5,
        x_end=3.5,
        y0=0.0,
        z0=0.0,
    )
    assert res.reached
    assert res.t_coord is not None
    assert abs(res.t_coord - res.t_flat) / res.t_flat < 0.05


def test_future_null_cov_propagates_forward_under_shift():
    """A +x-directed null ray must keep moving in +x even with a strong shift.

    Regression for the initialisation bug: ``k^mu=(1,1,0,0)`` + a generic null
    projection could select the backward root, sending the ray off the -x edge.
    """
    beta_x = 0.5
    g = _shifted_metric_grid(3, beta_x)[0, 0, 0]
    ginv = inverse_metric_4d(g)
    k = future_null_cov(g, np.array([1.0, 0.0, 0.0]))
    assert abs(null_hamiltonian(ginv, k)) < 1e-10  # genuinely null
    dx = ginv @ k
    assert dx[0] > 0.0  # future-directed
    assert dx[1] > 0.0  # propagates toward +x detector, not backward


def test_integrate_null_ray_reaches_detector_under_shift():
    """The full integrator reaches the +x detector through a shifted region."""
    n = 17
    g = _shifted_metric_grid(n, 0.5)
    spacing = (0.5, 0.5, 0.5)
    origin = np.array([-4.0, -4.0, -4.0])
    dg_inv = partial_inverse_metric(g, spacing)
    res = integrate_null_ray(
        g, dg_inv, origin, spacing, x_start=-3.5, x_end=3.5, y0=0.0, z0=0.0
    )
    assert res.reached
    assert res.max_h_rel < 1e-6  # constant metric: stays exactly null


def _stub_report(*, f_geo: float, h_ok: float, rel_drift: float) -> GeodesicFtlReport:
    return GeodesicFtlReport(
        f_geo=f_geo,
        t_min=1.0,
        t_flat=2.0,
        n_rays=5,
        n_reached=5,
        max_h_drift=1.0e-3,
        h_quality_ok=h_ok,
        max_h_rel_drift=rel_drift,
    )


def test_reliability_reprobe_certifies_shortcut_rejected_at_qd_resolution(monkeypatch):
    """A real shortcut that fails the H gate only because of QD-resolution
    discretization (Alcubierre: 65^3 rel-H 2.2e-2 FAIL, 129^3 5e-3 PASS) must be
    re-traced at >96^3 and certified, not silently zeroed."""
    calls: list[int] = []

    def fake(plotfile, *, n, half_width, n_rays, h_tol):
        calls.append(n)
        if n <= 65:
            return _stub_report(f_geo=0.31, h_ok=False, rel_drift=2.2e-2)
        return _stub_report(f_geo=0.31, h_ok=True, rel_drift=5.0e-3)

    monkeypatch.setattr(geo, "_geodesic_report_at_resolution", fake)
    report = compute_geodesic_ftl_from_plotfile("dummy", n=65)

    assert report is not None
    assert report.h_quality_ok  # certified after the re-probe
    assert calls == [65, GEO_REFINE_N]  # escalated exactly once, to >96
    assert GEO_REFINE_N > 96
    assert any("reliability re-probe" in note for note in report.notes)


def test_no_reprobe_when_base_resolution_is_already_reliable(monkeypatch):
    """A shortcut that already passes the gate at QD resolution is returned as-is;
    the expensive re-probe must not fire."""
    calls: list[int] = []

    def fake(plotfile, *, n, half_width, n_rays, h_tol):
        calls.append(n)
        return _stub_report(f_geo=0.05, h_ok=True, rel_drift=4.0e-3)

    monkeypatch.setattr(geo, "_geodesic_report_at_resolution", fake)
    report = compute_geodesic_ftl_from_plotfile("dummy", n=65)

    assert report is not None and report.h_quality_ok
    assert calls == [65]  # no escalation


def test_no_reprobe_without_a_shortcut(monkeypatch):
    """No coordinate shortcut (f_geo below floor) means there is nothing to
    certify, so a failed gate must not trigger a re-probe."""
    calls: list[int] = []

    def fake(plotfile, *, n, half_width, n_rays, h_tol):
        calls.append(n)
        return _stub_report(f_geo=0.0, h_ok=False, rel_drift=3.0e-2)

    monkeypatch.setattr(geo, "_geodesic_report_at_resolution", fake)
    report = compute_geodesic_ftl_from_plotfile("dummy", n=65)

    assert report is not None and not report.h_quality_ok
    assert calls == [65]  # no escalation: nothing to refine
