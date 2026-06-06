"""Gauge-invariant null geodesic ray-tracing on evolved GRTeclyn plotfiles.

Integrates null geodesics in the Hamiltonian (position-momentum) formulation::

    dx^mu/dλ = g^{mu nu} k_nu
    dk_mu/dλ = -1/2 (partial_mu g^{ab}) k_a k_b

This avoids explicit Christoffel symbols (second derivatives of the metric),
which are noisy on AMR-derived grids.  Only first derivatives of the inverse
metric are required.

The primary output is ``f_geo``: the fractional proper-time / coordinate-time
shortcut relative to flat-space propagation between asymptotic emitter and
detector planes at ``x = ±L_boundary``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from .warpfactory import _d_dx, stress_energy


@dataclass(frozen=True)
class NullRayResult:
    """Outcome of tracing one null ray."""

    reached: bool
    t_coord: float | None
    t_flat: float
    max_h_drift: float
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class GeodesicFtlReport:
    """Gauge-invariant null-geodesic travel-time diagnostic."""

    f_geo: float
    t_min: float | None
    t_flat: float
    n_rays: int
    n_reached: int
    max_h_drift: float
    h_quality_ok: bool
    notes: tuple[str, ...] = ()


def build_metric_3d_from_plotfile(
    plotfile: str | Path,
    *,
    n: int = 65,
    half_width: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], tuple[float, float, float]]:
    """Sample a static 4-metric ``g_{mu nu}(x,y,z)`` from one plotfile.

    Returns ``(g, origin, spacing)`` where ``g`` has shape ``(nx,ny,nz,4,4)``,
    ``origin`` is the grid lower corner ``(x0,y0,z0)``, and ``spacing`` is
    ``(dx,dy,dz)``.
    """
    import yt  # local import: heavy optional dependency

    ds = yt.load(str(plotfile))
    center = np.array([float(c.to_value()) for c in ds.domain_center])
    dleft = np.array([float(c.to_value()) for c in ds.domain_left_edge])
    dright = np.array([float(c.to_value()) for c in ds.domain_right_edge])

    if half_width is None:
        half_width = 0.45 * float(ds.domain_width[0].to_value())
    left = np.maximum(center - half_width, dleft)
    right = np.minimum(center + half_width, dright)
    dims = (n, n, n)
    ag = ds.arbitrary_grid(left, right, dims)

    def field(name: str) -> NDArray[np.float64]:
        try:
            return np.asarray(ag["boxlib", name], dtype=float)
        except Exception:  # noqa: BLE001
            return np.asarray(ag[name], dtype=float)

    chi = np.clip(field("chi"), 1.0e-10, None)
    inv_chi = 1.0 / chi
    alpha = np.clip(field("lapse"), 1.0e-10, None)
    beta_up = np.stack([field("shift1"), field("shift2"), field("shift3")], axis=-1)

    gamma = np.zeros(chi.shape + (3, 3), dtype=float)
    gamma[..., 0, 0] = field("h11") * inv_chi
    gamma[..., 0, 1] = gamma[..., 1, 0] = field("h12") * inv_chi
    gamma[..., 0, 2] = gamma[..., 2, 0] = field("h13") * inv_chi
    gamma[..., 1, 1] = field("h22") * inv_chi
    gamma[..., 1, 2] = gamma[..., 2, 1] = field("h23") * inv_chi
    gamma[..., 2, 2] = field("h33") * inv_chi

    beta_low = np.einsum("...ij,...j->...i", gamma, beta_up)
    beta_sq = np.einsum("...i,...i->...", beta_low, beta_up)

    g = np.zeros(chi.shape + (4, 4), dtype=float)
    g[..., 0, 0] = beta_sq - alpha * alpha
    g[..., 0, 1:] = beta_low
    g[..., 1:, 0] = beta_low
    g[..., 1:, 1:] = gamma

    spacing = tuple((right - left) / (np.array(dims) - 1))
    return g, left.astype(float), spacing


def inverse_metric_4d(g: NDArray[np.float64]) -> NDArray[np.float64]:
    """Invert a batch of 4x4 metrics.  Shape ``(..., 4, 4)`` -> same."""
    return np.linalg.inv(g)


def partial_inverse_metric(
    g: NDArray[np.float64],
    spacing: Sequence[float],
) -> NDArray[np.float64]:
    """``partial_mu g^{ab}`` for a static spatial grid.

    Returns array of shape ``(..., 4, 4, 4)`` where the last index is ``mu``.
    Time derivatives are zero (single Cauchy slice).
    """
    ginv = inverse_metric_4d(g)
    out = np.zeros(g.shape[:-2] + (4, 4, 4), dtype=float)
    for mu in range(1, 4):
        out[..., :, :, mu] = _d_dx(ginv, float(spacing[mu - 1]), axis=mu - 1)
    return out


def _clamp_index(
    x: NDArray[np.float64],
    origin: NDArray[np.float64],
    spacing: Sequence[float],
    shape: Sequence[int],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Continuous grid coordinates and fractional cell position."""
    rel = (x[1:] - origin[:3]) / np.asarray(spacing, dtype=float)
    n = np.array(shape, dtype=float)
    rel_clamped = np.clip(rel, 0.0, n - 1.001)
    i0 = np.floor(rel_clamped).astype(int)
    frac = rel_clamped - i0
    return i0, frac


def _trilinear(
    field: NDArray[np.float64],
    i0: NDArray[np.int64],
    frac: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Trilinear interpolation on a 3D field with trailing dimensions."""
    i0 = np.clip(i0, 0, np.array(field.shape[:3]) - 2)
    fx, fy, fz = frac[0], frac[1], frac[2]
    c000 = field[i0[0], i0[1], i0[2]]
    c100 = field[i0[0] + 1, i0[1], i0[2]]
    c010 = field[i0[0], i0[1] + 1, i0[2]]
    c110 = field[i0[0] + 1, i0[1] + 1, i0[2]]
    c001 = field[i0[0], i0[1], i0[2] + 1]
    c101 = field[i0[0] + 1, i0[1], i0[2] + 1]
    c011 = field[i0[0], i0[1] + 1, i0[2] + 1]
    c111 = field[i0[0] + 1, i0[1] + 1, i0[2] + 1]
    c00 = c000 * (1 - fx) + c100 * fx
    c10 = c010 * (1 - fx) + c110 * fx
    c01 = c001 * (1 - fx) + c101 * fx
    c11 = c011 * (1 - fx) + c111 * fx
    c0 = c00 * (1 - fy) + c10 * fy
    c1 = c01 * (1 - fy) + c11 * fy
    return c0 * (1 - fz) + c1 * fz


def _sample_metric(
    g: NDArray[np.float64],
    dg_inv: NDArray[np.float64],
    x: NDArray[np.float64],
    origin: NDArray[np.float64],
    spacing: Sequence[float],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    i0, frac = _clamp_index(x, origin, spacing, g.shape[:3])
    ginv = inverse_metric_4d(g)
    g_pt = _trilinear(g, i0, frac)
    ginv_pt = _trilinear(ginv, i0, frac)
    dg_pt = _trilinear(dg_inv, i0, frac)
    return g_pt, ginv_pt, dg_pt


def null_hamiltonian(ginv: NDArray[np.float64], k_cov: NDArray[np.float64]) -> float:
    """``H = 1/2 g^{mu nu} k_mu k_nu`` (should be 0 for null rays)."""
    return 0.5 * float(k_cov @ ginv @ k_cov)


def project_null(
    g: NDArray[np.float64],
    ginv: NDArray[np.float64],
    k_cov: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Rescale spatial covariant components to restore ``H = 0``."""
    k = k_cov.copy()
    g_spatial = g[1:, 1:]
    ginv_spatial = ginv[1:, 1:]
    k_sp = k[1:]
    a = float(k_sp @ ginv_spatial @ k_sp)
    b = 2.0 * float(np.dot(ginv[0, 1:], k_sp)) * k[0]
    c = float(ginv[0, 0]) * k[0] * k[0]
    disc = b * b - 4.0 * a * c
    if a <= 0.0 or disc < 0.0:
        return k
    # Choose root that keeps k^0 > 0 (future directed).
    k0_new = (-b + math.sqrt(disc)) / (2.0 * a)
    if k0_new <= 0.0:
        k0_new = (-b - math.sqrt(disc)) / (2.0 * a)
    k[0] = k0_new
    return k


def _hamiltonian_rhs(
    ginv: NDArray[np.float64],
    dg_inv: NDArray[np.float64],
    k_cov: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    dx = ginv @ k_cov
    dk = np.zeros(4, dtype=float)
    for mu in range(4):
        dk[mu] = -0.5 * float(k_cov @ dg_inv[..., mu] @ k_cov)
    return dx, dk


def integrate_null_ray(
    g: NDArray[np.float64],
    dg_inv: NDArray[np.float64],
    origin: NDArray[np.float64],
    spacing: Sequence[float],
    *,
    x_start: float,
    x_end: float,
    y0: float,
    z0: float,
    t0: float = 0.0,
    max_steps: int = 50_000,
    ds_init: float = 0.05,
    h_tol: float = 1.0e-6,
) -> NullRayResult:
    """Trace one null ray from ``x_start`` to ``x_end`` at fixed ``(y,z)``."""
    x = np.array([t0, x_start, y0, z0], dtype=float)
    t_flat = abs(x_end - x_start)

    g_pt, ginv0, _ = _sample_metric(g, dg_inv, x, origin, spacing)
    k_up = np.array([1.0, 1.0, 0.0, 0.0])
    k = g_pt @ k_up
    k = project_null(g_pt, ginv0, k)

    max_h = 0.0
    lam = 0.0
    ds = ds_init

    for _ in range(max_steps):
        if x[1] >= x_end:
            t_coord = abs(float(x[0] - t0))
            return NullRayResult(
                reached=True,
                t_coord=t_coord,
                t_flat=t_flat,
                max_h_drift=max_h,
            )

        g_pt, ginv_pt, dg_pt = _sample_metric(g, dg_inv, x, origin, spacing)
        h = abs(null_hamiltonian(ginv_pt, k))
        max_h = max(max_h, h)
        if h > h_tol:
            k = project_null(g_pt, ginv_pt, k)

        def rhs(xp: NDArray[np.float64], kp: NDArray[np.float64]) -> tuple[NDArray, NDArray]:
            gp, ginvp, dgp = _sample_metric(g, dg_inv, xp, origin, spacing)
            return _hamiltonian_rhs(ginvp, dgp, kp)

        k1x, k1k = rhs(x, k)
        k2x, k2k = rhs(x + 0.5 * ds * k1x, k + 0.5 * ds * k1k)
        k3x, k3k = rhs(x + 0.5 * ds * k2x, k + 0.5 * ds * k2k)
        k4x, k4k = rhs(x + ds * k3x, k + ds * k3k)
        x = x + (ds / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
        k = k + (ds / 6.0) * (k1k + 2 * k2k + 2 * k3k + k4k)
        lam += ds

        if x[1] < origin[0] or x[1] > origin[0] + (g.shape[0] - 1) * spacing[0]:
            return NullRayResult(
                reached=False,
                t_coord=None,
                t_flat=t_flat,
                max_h_drift=max_h,
                notes=("ray left grid",),
            )

    return NullRayResult(
        reached=False,
        t_coord=None,
        t_flat=t_flat,
        max_h_drift=max_h,
        notes=("max_steps exceeded",),
    )


def compute_geodesic_ftl_from_plotfile(
    plotfile: str | Path,
    *,
    n: int = 65,
    half_width: float | None = None,
    n_rays: int = 5,
    h_tol: float = 1.0e-6,
) -> GeodesicFtlReport | None:
    """Run a fan of null rays across the x-z midplane and return ``f_geo``."""
    try:
        g, origin, spacing = build_metric_3d_from_plotfile(
            plotfile, n=n, half_width=half_width
        )
    except Exception:
        return None

    dg_inv = partial_inverse_metric(g, spacing)
    shape = g.shape[:3]
    cy = origin[1] + 0.5 * (shape[1] - 1) * spacing[1]
    cz = origin[2] + 0.5 * (shape[2] - 1) * spacing[2]

    x_start = origin[0] + 0.05 * (shape[0] - 1) * spacing[0]
    x_end = origin[0] + 0.95 * (shape[0] - 1) * spacing[0]
    t_flat = x_end - x_start

    offsets = np.linspace(-0.1, 0.1, max(1, n_rays)) * (shape[2] - 1) * spacing[2]
    results: list[NullRayResult] = []
    for dz in offsets:
        res = integrate_null_ray(
            g,
            dg_inv,
            origin,
            spacing,
            x_start=x_start,
            x_end=x_end,
            y0=cy,
            z0=cz + float(dz),
            h_tol=h_tol,
        )
        results.append(res)

    reached = [r for r in results if r.reached and r.t_coord is not None]
    if not reached:
        return GeodesicFtlReport(
            f_geo=0.0,
            t_min=None,
            t_flat=t_flat,
            n_rays=len(results),
            n_reached=0,
            max_h_drift=max((r.max_h_drift for r in results), default=0.0),
            h_quality_ok=False,
            notes=("no rays reached detector",),
        )

    t_min = min(r.t_coord for r in reached if r.t_coord is not None)
    max_h = max(r.max_h_drift for r in results)
    h_ok = max_h <= 10.0 * h_tol
    f_geo = max(0.0, (t_flat - t_min) / t_flat) if t_flat > 0 else 0.0

    notes: list[str] = []
    if not h_ok:
        notes.append(f"null constraint drift high (max H={max_h:.2e})")

    return GeodesicFtlReport(
        f_geo=f_geo,
        t_min=t_min,
        t_flat=t_flat,
        n_rays=len(results),
        n_reached=len(reached),
        max_h_drift=max_h,
        h_quality_ok=h_ok,
        notes=tuple(notes),
    )


def integrate_null_energy_along_ray(
    g: NDArray[np.float64],
    dg_inv: NDArray[np.float64],
    T: NDArray[np.float64],
    origin: NDArray[np.float64],
    spacing: Sequence[float],
    *,
    x_start: float,
    x_end: float,
    y0: float,
    z0: float,
    sampling_width: float = 1.0,
) -> float | None:
    """Integrate ``T_{mu nu} k^mu k^nu`` along a null ray (QEI trajectory check)."""
    x = np.array([0.0, x_start, y0, z0], dtype=float)
    _, ginv0, _ = _sample_metric(g, dg_inv, x, origin, spacing)
    k = project_null(_sample_metric(g, dg_inv, x, origin, spacing)[0], ginv0, np.array([-1.0, 1.0, 0.0, 0.0]))

    integral = 0.0
    ds = 0.05
    prev_lam = 0.0

    for _ in range(50_000):
        if x[1] >= x_end:
            return integral

        g_pt, ginv_pt, dg_pt = _sample_metric(g, dg_inv, x, origin, spacing)
        T_pt = _trilinear(T, *_clamp_index(x, origin, spacing, g.shape[:3]))
        ginv_pt = ginv_pt
        k_contra = ginv_pt @ k
        integrand = float(k @ T_pt @ k_contra)
        integral += integrand * ds * sampling_width

        dx, dk = _hamiltonian_rhs(ginv_pt, dg_pt, k)
        x = x + ds * dx
        k = k + ds * dk
        prev_lam += ds

        if x[1] < origin[1]:
            return None

    return None


def compute_trajectory_qei_from_plotfile(
    plotfile: str | Path,
    *,
    n: int = 65,
    half_width: float | None = None,
    tau0: float = 1.0,
    ford_roman_c: float = 1.0,
) -> float | None:
    """Trajectory-based QEI integral vs Ford-Roman ``-C/tau^4`` bound."""
    try:
        g, origin, spacing = build_metric_3d_from_plotfile(
            plotfile, n=n, half_width=half_width
        )
        g_stack = g[np.newaxis, ...]
        T = stress_energy(g_stack, (1.0, *spacing))[0]
    except Exception:
        return None

    dg_inv = partial_inverse_metric(g, spacing)
    shape = g.shape[:3]
    cy = origin[1] + 0.5 * (shape[1] - 1) * spacing[1]
    cz = origin[2] + 0.5 * (shape[2] - 1) * spacing[2]
    x_start = origin[0] + 0.05 * (shape[0] - 1) * spacing[0]
    x_end = origin[0] + 0.95 * (shape[0] - 1) * spacing[0]

    val = integrate_null_energy_along_ray(
        g, dg_inv, T, origin, spacing,
        x_start=x_start, x_end=x_end, y0=cy, z0=cz,
        sampling_width=tau0,
    )
    if val is None:
        return None
    bound = -ford_roman_c / (tau0 ** 4)
    return float(val - bound)
