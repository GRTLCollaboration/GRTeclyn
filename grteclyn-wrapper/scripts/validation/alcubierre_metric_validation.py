#!/usr/bin/env python3
"""Validate the FTL metrics on an analytic Alcubierre warp metric.

This is a *metric-first* positive control: we prescribe the textbook Alcubierre
3+1 metric directly (no GRTresna constraint solve, no GRTeclyn evolution) and run
the exact probes the QD campaign uses --

  * ``integrate_null_ray`` -> the gauge-invariant geodesic shortcut ``f_geo``
  * ``stress_energy`` -> the matter the metric *requires* (T = G/8pi), and its
    null-energy-condition slack.

Purpose: prove that (a) our geodesic probe genuinely registers a real
superluminal shortcut when one exists, and (b) our energy-condition probe
correctly flags the exotic (NEC-violating) matter Alcubierre demands.  If the
probe lights up here and stays dark on physical-scalar candidates, the negative
QD results are trustworthy; if it stays dark *here*, the probe has a bug.

Alcubierre (Class. Quantum Grav. 11 (1994) L73), bubble centred at the origin,
moving along +x with speed ``v_s``::

    alpha = 1,   gamma_ij = delta_ij,   beta^x = -v_s f(r_s)

    f(r_s) = [tanh(sigma (r_s + R)) - tanh(sigma (r_s - R))] / (2 tanh(sigma R))

with ``r_s = sqrt(x^2 + y^2 + z^2)``.  A +x null ray then has coordinate speed
``1 + v_s f`` -- superluminal inside the bubble wall -- so its crossing time is
shorter than flat space (``f_geo > 0``).
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from numpy.typing import NDArray

from grteclyn_wrapper.metrics.probes.ftl.geodesic import (
    integrate_null_ray,
    partial_inverse_metric,
)
from grteclyn_wrapper.metrics.probes.warpfactory import stress_energy


def alcubierre_shape(r: NDArray, *, radius: float, sigma: float) -> NDArray:
    """Alcubierre top-hat shape function ``f(r_s)`` in [0, 1]."""
    return (
        np.tanh(sigma * (r + radius)) - np.tanh(sigma * (r - radius))
    ) / (2.0 * math.tanh(sigma * radius))


def build_alcubierre_metric(
    *,
    n: int,
    half_width: float,
    v_s: float,
    radius: float,
    sigma: float,
) -> tuple[NDArray, NDArray, tuple[float, float, float]]:
    """Covariant 4-metric ``g_{mu nu}(x,y,z)`` on an ``(n,n,n,4,4)`` grid.

    Matches the ADM->4-metric convention used by ``build_metric_3d_from_plotfile``
    (lapse 1, flat spatial metric, pure shift), so the geodesic probe sees exactly
    what it would off a plotfile.
    """
    axis = np.linspace(-half_width, half_width, n)
    spacing = (axis[1] - axis[0],) * 3
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    r = np.sqrt(x * x + y * y + z * z)

    f = alcubierre_shape(r, radius=radius, sigma=sigma)
    beta_x = -v_s * f  # contravariant == covariant (gamma = I)

    g = np.zeros((n, n, n, 4, 4), dtype=float)
    g[..., 0, 0] = beta_x * beta_x - 1.0  # beta_i beta^i - alpha^2
    g[..., 0, 1] = g[..., 1, 0] = beta_x  # g_tx = beta_x
    for i in range(1, 4):
        g[..., i, i] = 1.0

    origin = np.array([-half_width, -half_width, -half_width], dtype=float)
    return g, origin, spacing


def geodesic_ftl(
    g: NDArray,
    origin: NDArray,
    spacing: tuple[float, float, float],
    *,
    n_rays: int = 5,
) -> dict:
    """Run the campaign's null-ray fan on an in-memory metric (mirrors
    ``compute_geodesic_ftl_from_plotfile`` but skips the yt plotfile load)."""
    dg_inv = partial_inverse_metric(g, spacing)
    shape = g.shape[:3]
    cy = origin[1] + 0.5 * (shape[1] - 1) * spacing[1]
    cz = origin[2] + 0.5 * (shape[2] - 1) * spacing[2]
    x_start = origin[0] + 0.05 * (shape[0] - 1) * spacing[0]
    x_end = origin[0] + 0.95 * (shape[0] - 1) * spacing[0]
    t_flat = x_end - x_start

    offsets = np.linspace(-0.1, 0.1, max(1, n_rays)) * (shape[2] - 1) * spacing[2]
    reached = []
    max_h_rel = 0.0
    for dz in offsets:
        res = integrate_null_ray(
            g, dg_inv, origin, spacing,
            x_start=x_start, x_end=x_end, y0=cy, z0=cz + float(dz),
        )
        max_h_rel = max(max_h_rel, res.max_h_rel)
        if res.reached and res.t_coord is not None:
            reached.append(res.t_coord)

    if not reached:
        return {"f_geo": 0.0, "t_min": None, "t_flat": t_flat,
                "n_reached": 0, "n_rays": len(offsets), "max_h_rel": max_h_rel}
    t_min = min(reached)
    return {
        "f_geo": max(0.0, (t_flat - t_min) / t_flat),
        "t_min": t_min, "t_flat": t_flat,
        "n_reached": len(reached), "n_rays": len(offsets),
        "max_h_rel": max_h_rel,
    }


def min_nec(g: NDArray, spacing: tuple[float, float, float]) -> float:
    """Minimum null-energy-condition contraction ``T_{mu nu} l^mu l^nu`` over a
    fan of spatial null directions, dropping the noisy boundary layer.

    Negative => the metric requires exotic (NEC-violating) matter.
    """
    # stress_energy differentiates along time too, so feed >=5 identical slices
    # (static metric => d_t g = 0) and read the central one.
    g_t = np.repeat(g[np.newaxis, ...], 5, axis=0)
    T = stress_energy(g_t, (1.0, *spacing))[2]
    c = (slice(3, -3),) * 3  # crop finite-difference boundary
    Tc = T[c]
    worst = math.inf
    # Null vectors l^mu = (1, n_hat) are null in flat-background coordinates; good
    # enough to expose the sign of the NEC integrand for this diagnostic.
    for ang in np.linspace(0.0, math.pi, 8):
        n_hat = np.array([math.cos(ang), math.sin(ang), 0.0])
        l = np.array([1.0, *n_hat])
        val = np.einsum("...ab,a,b->...", Tc, l, l)
        worst = min(worst, float(val.min()))
    return worst


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=65, help="grid points per axis")
    ap.add_argument("--half-width", type=float, default=10.0)
    ap.add_argument("--v-s", type=float, default=2.0, help="bubble speed (>1 = superluminal)")
    ap.add_argument("--radius", type=float, default=4.0, help="bubble radius R")
    ap.add_argument("--sigma", type=float, default=2.0, help="wall steepness")
    args = ap.parse_args()

    print("=" * 70)
    print("Alcubierre metric-first validation of the FTL probes")
    print("=" * 70)
    print(f"grid {args.n}^3  half_width={args.half_width}  "
          f"R={args.radius}  sigma={args.sigma}")
    print()

    # Flat control (v_s = 0): probe MUST report no shortcut and no exotic matter.
    g0, o0, s0 = build_alcubierre_metric(
        n=args.n, half_width=args.half_width, v_s=0.0,
        radius=args.radius, sigma=args.sigma,
    )
    flat = geodesic_ftl(g0, o0, s0)
    print(f"[control v_s=0.0] f_geo={flat['f_geo']:.4f}  "
          f"reached={flat['n_reached']}/{flat['n_rays']}  "
          f"rel_H_drift={flat['max_h_rel']:.2e}  min_NEC={min_nec(g0, s0):+.3e}")

    # Superluminal bubble: probe SHOULD report a shortcut AND exotic matter.
    g, o, s = build_alcubierre_metric(
        n=args.n, half_width=args.half_width, v_s=args.v_s,
        radius=args.radius, sigma=args.sigma,
    )
    warp = geodesic_ftl(g, o, s)
    nec = min_nec(g, s)
    print(f"[warp    v_s={args.v_s:.1f}] f_geo={warp['f_geo']:.4f}  "
          f"reached={warp['n_reached']}/{warp['n_rays']}  "
          f"rel_H_drift={warp['max_h_rel']:.2e}  min_NEC={nec:+.3e}")
    print()

    h_ok = warp["max_h_rel"] <= 1.0e-2
    print("VERDICT")
    print(f"  geodesic probe detects shortcut : {warp['f_geo'] > 0.01}  "
          f"(f_geo={warp['f_geo']:.3f})")
    print(f"  null-ray bundle reliable (H gate): {h_ok}  "
          f"(rel drift {warp['max_h_rel']:.2e} <= 1e-2)")
    print(f"  energy-condition probe flags NEC : {nec < 0.0}  "
          f"(min_NEC={nec:+.3e})")
    print(f"  flat control stays dark          : "
          f"{flat['f_geo'] < 0.01 and abs(min_nec(g0, s0)) < 1e-6}")


if __name__ == "__main__":
    main()
