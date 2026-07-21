#!/usr/bin/env python3
"""Referee-response diagnostics for candidate 146 (RM run).

Two products, both built from the *same* cached metric stack and the *same*
ray-fan geometry used by the production evolving-geodesic probe:

1. Endpoint gauge table (referee critique A): lapse alpha, shift beta^i, and the
   spatial metric gamma_ij evaluated at the fixed emission point A and receiver
   point B over the emission window t in [0, 15]. Establishes whether A and B
   behave as (nearly) stationary, asymptotically flat observers so that the
   coordinate-time advantage equals a proper-time advantage.

2. Along-ray null-constraint profile (referee critique D): the normalized null
   Hamiltonian k_mu k^mu / (k^t)^2 sampled every integration step of the
   headline ray (t_emit = 4), demonstrating that the ray stays orders of
   magnitude tighter than the 1e-2 reliability gate.

Run:
    cd grteclyn-wrapper
    PYTHONPATH=src python ../research/neuralspacetime/analysis/endpoint_gauge_and_null_constraint.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from grteclyn_wrapper.metrics.probes.ftl.metric_stack_cache import (
    evolving_field_from_metric_stack_cache,
    metric_stack_dir,
)
from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import _ray_fan_geometry
from grteclyn_wrapper.core.site_paths import grteclyn_root
from grteclyn_wrapper.metrics.probes.ftl.geodesic import (
    future_null_cov,
    null_hamiltonian,
    project_null,
    _hamiltonian_rhs,
)

REPO = grteclyn_root()
RUN = REPO / "runs/grtresna_promote/bcma_rm_L128_N256_t30_hq_eval000146"
OUT = REPO / "research/neuralspacetime/article/data"

T_EMIT = 4.0
AXIS = 0
EMIT_WINDOW = np.arange(0.0, 15.001, 1.0)


def adm_from_g(g: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Return (alpha, beta^i, gamma_ij) from a 4-metric g_{mu nu}.

    g_00 = -alpha^2 + beta_k beta^k, g_0i = beta_i, g_ij = gamma_ij.
    """
    gamma = g[1:, 1:]
    beta_low = g[0, 1:]
    gamma_inv = np.linalg.inv(gamma)
    beta_up = gamma_inv @ beta_low
    alpha2 = float(beta_up @ beta_low - g[0, 0])
    alpha = float(np.sqrt(alpha2)) if alpha2 > 0 else float("nan")
    return alpha, beta_up, gamma


def endpoint_gauge_table(field) -> str:
    prop_start, prop_end, y0, z0, t_flat, _ = _ray_fan_geometry(
        field, axis=AXIS, n_rays=5
    )
    A = (prop_start, y0, z0)
    B = (prop_end, y0, z0)
    print(f"# t_flat = {t_flat:.6f}")
    print(f"# A = (x={A[0]:.4f}, y={A[1]:.4f}, z={A[2]:.4f})")
    print(f"# B = (x={B[0]:.4f}, y={B[1]:.4f}, z={B[2]:.4f})")

    rows = ["t alpha_A beta_x_A beta_mag_A alpha_B beta_x_B beta_mag_B"]
    for t in EMIT_WINDOW:
        out = []
        for pt in (A, B):
            x = np.array([t, pt[0], pt[1], pt[2]], dtype=float)
            g = field.sample_g(x)
            alpha, beta, _ = adm_from_g(g)
            out.append((alpha, float(beta[0]), float(np.linalg.norm(beta))))
        rows.append(
            f"{t:.2f} "
            f"{out[0][0]:.6f} {out[0][1]:+.6e} {out[0][2]:.6e} "
            f"{out[1][0]:.6f} {out[1][1]:+.6e} {out[1][2]:.6e}"
        )
    return "\n".join(rows) + "\n"


def proper_time_correction(field) -> dict:
    """Receiver-clock (proper-time) advantage in the actual curved spacetime.

    A static receiver at B measures proper time dtau = sqrt(-g_00(B)) dt. The
    warp signal arrives at coordinate t_arr; the matched flat-rate signal would
    arrive at t_emit + t_flat. Expressing both arrivals in the receiver's own
    proper time removes any lapse/shift gauge ambiguity in the timing claim.
    """
    prop_start, prop_end, y0, z0, t_flat, _ = _ray_fan_geometry(
        field, axis=AXIS, n_rays=5
    )
    B = (prop_end, y0, z0)

    # From the RM evolving_geodesic.json headline ray.
    t_arr_curved = 15.48955004403862

    def dtau_dt(t: float) -> float:
        g = field.sample_g(np.array([t, B[0], B[1], B[2]], dtype=float))
        val = -float(g[0, 0])
        return float(np.sqrt(val)) if val > 0 else float("nan")

    def proper_interval(t0: float, t1: float, n: int = 4000) -> float:
        ts = np.linspace(t0, t1, n)
        vals = np.array([dtau_dt(t) for t in ts])
        return float(np.trapz(vals, ts))

    t_arr_flat = T_EMIT + t_flat
    tau_curved = proper_interval(T_EMIT, t_arr_curved)
    tau_flat = proper_interval(T_EMIT, t_arr_flat)
    f_proper = (tau_flat - tau_curved) / tau_flat if tau_flat > 0 else 0.0
    return {
        "t_flat": t_flat,
        "t_arr_curved": t_arr_curved,
        "t_arr_flat": t_arr_flat,
        "tau_curved": tau_curved,
        "tau_flat": tau_flat,
        "f_geo_coord": (t_flat - (t_arr_curved - T_EMIT)) / t_flat,
        "f_geo_proper": f_proper,
    }


def along_ray_null_constraint(field) -> tuple[str, dict]:
    """Re-trace the central headline ray, logging k_mu k^mu / (k^t)^2 per step."""
    prop_start, prop_end, y0, z0, _, _ = _ray_fan_geometry(field, axis=AXIS, n_rays=5)
    prop_idx = AXIS + 1
    x = np.array([T_EMIT, prop_start, y0, z0], dtype=float)
    n_hat = np.zeros(3)
    n_hat[AXIS] = 1.0
    g, ginv, _ = field.sample(x)
    k = future_null_cov(g, n_hat)

    ds = 0.05
    h_tol = 1.0e-6
    rows = ["s x_pos t kk_over_kt2 abs_H"]
    s = 0.0
    max_norm = 0.0
    for _ in range(50_000):
        if x[prop_idx] >= prop_end:
            break
        g, ginv, dg = field.sample(x)
        H = null_hamiltonian(ginv, k)          # 0.5 g^{ab} k_a k_b
        kk = 2.0 * H                            # g^{ab} k_a k_b = k_mu k^mu
        kt = float((ginv @ k)[0])               # k^t = g^{t nu} k_nu
        norm = abs(kk) / (kt * kt) if kt != 0 else 0.0
        max_norm = max(max_norm, norm)
        rows.append(
            f"{s:.4f} {x[prop_idx]:.5f} {x[0]:.5f} {norm:.6e} {abs(H):.6e}"
        )
        if abs(H) > h_tol:
            k = project_null(g, ginv, k, dx_ref=(ginv @ k)[1:])

        def rhs(xp, kp):
            _g, ginvp, dgp = field.sample(xp)
            return _hamiltonian_rhs(ginvp, dgp, kp)

        k1x, k1k = rhs(x, k)
        k2x, k2k = rhs(x + 0.5 * ds * k1x, k + 0.5 * ds * k1k)
        k3x, k3k = rhs(x + 0.5 * ds * k2x, k + 0.5 * ds * k2k)
        k4x, k4k = rhs(x + ds * k3x, k + ds * k3k)
        x = x + (ds / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
        k = k + (ds / 6.0) * (k1k + 2 * k2k + 2 * k3k + k4k)
        s += ds

    return "\n".join(rows) + "\n", {"max_kk_over_kt2": max_norm}


def main() -> None:
    cache = metric_stack_dir(RUN / "small_data")
    field = evolving_field_from_metric_stack_cache(cache)
    if field is None:
        raise SystemExit(f"could not build metric field from {cache}")

    OUT.mkdir(parents=True, exist_ok=True)

    gauge = endpoint_gauge_table(field)
    (OUT / "endpoint_gauge_rm.txt").write_text(gauge, encoding="utf-8")
    print(gauge)

    ray, summary = along_ray_null_constraint(field)
    (OUT / "null_constraint_ray_rm.txt").write_text(ray, encoding="utf-8")
    print(f"# max k_mu k^mu / (k^t)^2 along headline ray = "
          f"{summary['max_kk_over_kt2']:.6e}")

    pt = proper_time_correction(field)
    print("# proper-time correction (receiver clock at B):")
    for key, val in pt.items():
        print(f"#   {key} = {val:.6f}")


if __name__ == "__main__":
    main()
