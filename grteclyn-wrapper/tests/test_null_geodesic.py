"""Tests for gauge-invariant null geodesic ray-tracing."""

from __future__ import annotations

import numpy as np

from grteclyn_wrapper.metrics.null_geodesic import (
    build_metric_3d_from_plotfile,
    integrate_null_ray,
    inverse_metric_4d,
    null_hamiltonian,
    partial_inverse_metric,
    project_null,
)


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
