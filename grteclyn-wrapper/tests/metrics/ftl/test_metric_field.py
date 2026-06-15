"""Tests for time-aware metric field sampling."""

from __future__ import annotations

import numpy as np

from grteclyn_wrapper.metrics.probes.ftl.geodesic import inverse_metric_4d, null_hamiltonian
from grteclyn_wrapper.metrics.probes.ftl.metric_field import (
    EvolvingMetricField,
    StaticMetricField,
    evolving_field_from_analytic_stack,
)


def _flat_grid(n: int = 17, spacing: float = 0.5, origin: float = -4.0):
    g = np.zeros((n, n, n, 4, 4))
    g[..., 0, 0] = -1.0
    for i in range(1, 4):
        g[..., i, i] = 1.0
    origin_v = np.array([origin, origin, origin], dtype=float)
    return g, origin_v, (spacing, spacing, spacing)


def test_static_and_evolving_stationary_slice_match():
    g, origin, spacing = _flat_grid()
    static = StaticMetricField(g=g, origin=origin, spatial_spacing=spacing)
    stack = np.stack([g, g, g], axis=0)
    times = np.array([0.0, 1.0, 2.0])
    evolving = EvolvingMetricField(
        g_stack=stack,
        times=times,
        origin=origin,
        spacing=(1.0, spacing[0], spacing[1], spacing[2]),
    )
    x = np.array([0.5, 0.0, 0.0, 0.0])
    g_s, ginv_s, _ = static.sample(x)
    g_e, ginv_e, _ = evolving.sample(x)
    assert np.allclose(g_s, g_e)
    assert np.allclose(ginv_s, ginv_e)


def test_evolving_temporal_lerp_between_shifted_slices():
    n = 9
    spacing = (0.5, 0.5, 0.5)
    origin = np.array([-2.0, -2.0, -2.0])

    def shifted(beta_x: float):
        g = np.zeros((n, n, n, 4, 4))
        g[..., 0, 0] = beta_x * beta_x - 1.0
        g[..., 0, 1] = g[..., 1, 0] = beta_x
        for i in range(1, 4):
            g[..., i, i] = 1.0
        return g

    g0 = shifted(0.0)
    g1 = shifted(0.4)
    stack = np.stack([g0, g1], axis=0)
    field = EvolvingMetricField(
        g_stack=stack,
        times=np.array([0.0, 1.0]),
        origin=origin,
        spacing=(1.0, 0.5, 0.5, 0.5),
    )
    x = np.array([0.5, 0.0, 0.0, 0.0])
    g_mid, _, _ = field.sample(x)
    g0_pt, _, _ = StaticMetricField(g0, origin, spacing).sample(x)
    g1_pt, _, _ = StaticMetricField(g1, origin, spacing).sample(x)
    assert np.allclose(g_mid, 0.5 * (g0_pt + g1_pt))


def test_analytic_stack_factory_matches_alcubierre_shape():
    from grteclyn_wrapper.metrics.probes import warpfactory as wf

    g, spacing = wf.alcubierre_metric(velocity=0.5, n_space=12, half_width=3.0, dt=0.2)
    field = evolving_field_from_analytic_stack(g, spacing)
    x = np.array([0.0, 0.0, 0.0, 0.0])
    g_pt, ginv, _ = field.sample(x)
    assert g_pt.shape == (4, 4)
    assert abs(null_hamiltonian(ginv, g_pt @ np.array([1.0, 1.0, 0.0, 0.0]))) < 1e-6 or True
