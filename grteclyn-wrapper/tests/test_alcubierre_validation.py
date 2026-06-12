"""Metric-first positive control: the FTL probes on a textbook Alcubierre warp.

These tests prescribe the analytic Alcubierre 3+1 metric directly (no constraint
solve, no evolution) and assert that the *same* probes the QD campaign uses

  * detect the gauge-invariant superluminal shortcut (``f_geo`` > 0), and
  * flag the exotic (NEC-violating) matter the metric requires (min NEC < 0),

while a flat (``v_s = 0``) control stays dark on both.  If these ever fail, a
negative QD result can no longer be trusted -- the probe itself is broken.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

from grteclyn_wrapper.metrics.probes.ftl.geodesic import (
    integrate_null_ray,
    partial_inverse_metric,
)
from grteclyn_wrapper.metrics.probes.warpfactory import stress_energy


def _alcubierre_metric(
    *, n: int, half_width: float, v_s: float, radius: float, sigma: float
) -> tuple[NDArray, NDArray, tuple[float, float, float]]:
    """Covariant 4-metric for an Alcubierre bubble at the origin moving along +x.

    lapse=1, gamma=I, beta^x = -v_s f(r_s); matches the ADM->4-metric convention
    of ``build_metric_3d_from_plotfile`` so the probe sees a plotfile-equivalent.
    """
    axis = np.linspace(-half_width, half_width, n)
    spacing = (axis[1] - axis[0],) * 3
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    r = np.sqrt(x * x + y * y + z * z)
    f = (np.tanh(sigma * (r + radius)) - np.tanh(sigma * (r - radius))) / (
        2.0 * math.tanh(sigma * radius)
    )
    beta_x = -v_s * f

    g = np.zeros((n, n, n, 4, 4), dtype=float)
    g[..., 0, 0] = beta_x * beta_x - 1.0
    g[..., 0, 1] = g[..., 1, 0] = beta_x
    for i in range(1, 4):
        g[..., i, i] = 1.0
    origin = np.array([-half_width, -half_width, -half_width], dtype=float)
    return g, origin, spacing


def _f_geo(g: NDArray, origin: NDArray, spacing: tuple[float, float, float]) -> float:
    dg_inv = partial_inverse_metric(g, spacing)
    shape = g.shape[:3]
    cy = origin[1] + 0.5 * (shape[1] - 1) * spacing[1]
    cz = origin[2] + 0.5 * (shape[2] - 1) * spacing[2]
    x_start = origin[0] + 0.05 * (shape[0] - 1) * spacing[0]
    x_end = origin[0] + 0.95 * (shape[0] - 1) * spacing[0]
    t_flat = x_end - x_start
    res = integrate_null_ray(
        g, dg_inv, origin, spacing, x_start=x_start, x_end=x_end, y0=cy, z0=cz
    )
    assert res.reached and res.t_coord is not None
    return max(0.0, (t_flat - res.t_coord) / t_flat)


def _min_nec(g: NDArray, spacing: tuple[float, float, float]) -> float:
    g_t = np.repeat(g[np.newaxis, ...], 5, axis=0)
    T = stress_energy(g_t, (1.0, *spacing))[2]
    c = (slice(3, -3),) * 3
    Tc = T[c]
    worst = math.inf
    for ang in np.linspace(0.0, math.pi, 8):
        l = np.array([1.0, math.cos(ang), math.sin(ang), 0.0])
        worst = min(worst, float(np.einsum("...ab,a,b->...", Tc, l, l).min()))
    return worst


# Soft wall (sigma=1) so a modest grid resolves it -- keeps the test fast while
# the shortcut magnitude and NEC sign are already converged (see the convergence
# study in scripts/validation/alcubierre_metric_validation.py).
_KW = dict(n=49, half_width=10.0, radius=4.0, sigma=1.0)


def test_geodesic_probe_detects_alcubierre_shortcut() -> None:
    g, o, s = _alcubierre_metric(v_s=2.0, **_KW)
    f_geo = _f_geo(g, o, s)
    # v_s=2 bubble => ~30%+ coordinate-time shortcut; converged value ~0.33.
    assert f_geo > 0.1


def test_energy_condition_probe_flags_alcubierre_exotic_matter() -> None:
    g, _origin, s = _alcubierre_metric(v_s=2.0, **_KW)
    assert _min_nec(g, s) < 0.0


def test_flat_control_stays_dark() -> None:
    g, o, s = _alcubierre_metric(v_s=0.0, **_KW)
    assert _f_geo(g, o, s) < 1.0e-3
    assert abs(_min_nec(g, s)) < 1.0e-9
