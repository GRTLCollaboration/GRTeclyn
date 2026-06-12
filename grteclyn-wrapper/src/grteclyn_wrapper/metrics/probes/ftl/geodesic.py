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
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from ..warpfactory import _d_dx, stress_energy


# Reliability gate on the null constraint ``H = ½ g^{μν} k_μ k_ν`` (exactly 0 for
# a null ray).  The drift is graded *relative* to the natural magnitude of the
# terms that build ``H`` -- ``½(|g_tt (k^t)² | + |k_i γ^{ij} k_j|)`` -- because
# the metric is only C0 (trilinear) interpolated, so the absolute per-step drift
# has an interpolation-set floor (~1e-3) that no step size can beat.  A
# convergence study on the v7 elites (step 0.05→0.0125, grid 65→129) showed
# ``f_geo`` stable to ~5% while this relative drift held at 1e-3-2e-3, so 1e-2 is
# a conservative "integration stayed on the null cone" bar that the old absolute
# 1e-5 bound could never satisfy.
H_REL_TOL: float = 1.0e-2


# Resolution for the reliability re-probe.  At QD resolution (65^3) the relative
# H-drift gate is dominated by a C0 (trilinear) interpolation discretization
# artifact, not real null-cone error: an Alcubierre convergence study showed the
# drift halving per refinement (65^3: 2.2e-2 FAIL, 97^3: 1.1e-2 borderline,
# 129^3: 5e-3 PASS) while ``f_geo`` stayed at ~0.315.  A genuine sharp-walled
# shortcut would therefore be silently zeroed at QD resolution.  129^3 is the
# first grid that reliably certifies even a sharp (sigma=2) wall, so when the
# base probe finds a shortcut but fails the gate we re-trace the rays here before
# discarding a possibly-real warp.
GEO_REFINE_N: int = 129

# Minimum coordinate shortcut worth paying for a higher-resolution re-probe
# (matches the scoring floor ``GEO_FTL_FLOOR``): below this there is nothing to
# certify, so the expensive re-trace is skipped.
GEO_REFINE_FLOOR: float = 1.0e-3


@dataclass(frozen=True)
class NullRayResult:
    """Outcome of tracing one null ray."""

    reached: bool
    t_coord: float | None
    t_flat: float
    max_h_drift: float
    max_h_rel: float = 0.0
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
    max_h_rel_drift: float = 0.0
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
    ginv: NDArray[np.float64],
    dg_inv: NDArray[np.float64],
    x: NDArray[np.float64],
    origin: NDArray[np.float64],
    spacing: Sequence[float],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    i0, frac = _clamp_index(x, origin, spacing, g.shape[:3])
    g_pt = _trilinear(g, i0, frac)
    ginv_pt = _trilinear(ginv, i0, frac)
    dg_pt = _trilinear(dg_inv, i0, frac)
    return g_pt, ginv_pt, dg_pt


def null_hamiltonian(ginv: NDArray[np.float64], k_cov: NDArray[np.float64]) -> float:
    """``H = 1/2 g^{mu nu} k_mu k_nu`` (should be 0 for null rays)."""
    return 0.5 * float(k_cov @ ginv @ k_cov)


def _null_relative_drift(ginv: NDArray[np.float64], k_cov: NDArray[np.float64]) -> float:
    """``|H|`` normalised by the magnitude of its time and spatial parts.

    ``H`` is the near-cancellation of a negative time term and a positive spatial
    term, each ``O(|k|^2)``.  Dividing by their summed magnitude yields a
    scale-free "how far off the null cone" measure that is meaningful regardless
    of the ray's momentum normalisation.
    """
    h = abs(0.5 * float(k_cov @ ginv @ k_cov))
    t_term = abs(float(ginv[0, 0]) * k_cov[0] * k_cov[0])
    sp_term = abs(float(k_cov[1:] @ ginv[1:, 1:] @ k_cov[1:]))
    scale = 0.5 * (t_term + sp_term)
    return h / scale if scale > 1.0e-30 else 0.0


def future_null_cov(
    g_pt: NDArray[np.float64],
    n_hat: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Covariant momentum of a *future-directed* null ray whose contravariant
    spatial velocity points along ``n_hat`` (a 3-vector).

    Solving on the contravariant side guarantees the ray's coordinate velocity
    ``dx^i/dλ = (g^{-1} k)^i`` is parallel to ``n_hat`` and ``dt/dλ > 0``.  The
    previous initialisation (``k^μ = (1,1,0,0)`` then a generic null projection)
    did not pin the propagation direction: the projection could select the
    *backward* null root, sending the ray straight off the −x boundary.  That is
    why almost no rays reached the detector and the gauge-invariant gate never
    passed.

    With ``k^μ = (k^t, n_hat)`` the null condition ``g_{μν} k^μ k^ν = 0`` is a
    quadratic in ``k^t`` with ``a = g_tt < 0`` and ``c = n·γ·n > 0``, so the
    roots have opposite sign and the unique positive (future-directed) root is
    taken.
    """
    a = float(g_pt[0, 0])
    b = 2.0 * float(g_pt[0, 1:] @ n_hat)
    c = float(n_hat @ g_pt[1:, 1:] @ n_hat)
    disc = b * b - 4.0 * a * c
    if a == 0.0 or disc < 0.0:
        k_up = np.array([1.0, *n_hat], dtype=float)
    else:
        sq = math.sqrt(disc)
        roots = ((-b + sq) / (2.0 * a), (-b - sq) / (2.0 * a))
        k_t = max(roots)  # opposite-sign roots => max() is the positive one
        if k_t <= 0.0:
            k_t = 1.0
        k_up = np.array([k_t, *n_hat], dtype=float)
    return g_pt @ k_up


def project_null(
    g: NDArray[np.float64],
    ginv: NDArray[np.float64],
    k_cov: NDArray[np.float64],
    dx_ref: NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Rescale spatial covariant components to restore ``H = 0``.

    When ``dx_ref`` (the current coordinate velocity direction) is supplied, the
    null root whose resulting spatial velocity stays aligned with it is
    preferred, so an in-flight re-projection cannot silently reverse the ray.
    """
    k = k_cov.copy()
    ginv_spatial = ginv[1:, 1:]
    k_sp = k[1:]
    a = float(k_sp @ ginv_spatial @ k_sp)
    b = 2.0 * float(np.dot(ginv[0, 1:], k_sp)) * k[0]
    c = float(ginv[0, 0]) * k[0] * k[0]
    disc = b * b - 4.0 * a * c
    if a <= 0.0 or disc < 0.0:
        return k
    # The quadratic solves for a scale factor on k_i with fixed k_0.
    roots = (
        (-b + math.sqrt(disc)) / (2.0 * a),
        (-b - math.sqrt(disc)) / (2.0 * a),
    )
    candidates: list[NDArray[np.float64]] = []
    for scale in roots:
        kk = k.copy()
        kk[1:] = scale * k_sp
        if np.isfinite(scale) and (ginv @ kk)[0] > 0.0:
            candidates.append(kk)
    if dx_ref is not None:
        aligned = [kk for kk in candidates if np.dot((ginv @ kk)[1:], dx_ref) > 0.0]
        if aligned:
            candidates = aligned
    if candidates:
        # Prefer the smallest positive correction to preserve ray direction.
        k = min(candidates, key=lambda kk: abs(np.linalg.norm(kk[1:]) - np.linalg.norm(k_sp)))
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
    ginv = inverse_metric_4d(g)

    g_pt, ginv0, _ = _sample_metric(g, ginv, dg_inv, x, origin, spacing)
    # Launch a future-directed null ray whose coordinate velocity points toward
    # the detector (+x); see ``future_null_cov``.
    n_hat = np.array([1.0, 0.0, 0.0])
    k = future_null_cov(g_pt, n_hat)

    max_h = 0.0
    max_h_rel = 0.0
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
                max_h_rel=max_h_rel,
            )

        g_pt, ginv_pt, dg_pt = _sample_metric(g, ginv, dg_inv, x, origin, spacing)
        h = abs(null_hamiltonian(ginv_pt, k))
        max_h = max(max_h, h)
        max_h_rel = max(max_h_rel, _null_relative_drift(ginv_pt, k))
        if h > h_tol:
            k = project_null(g_pt, ginv_pt, k, dx_ref=(ginv_pt @ k)[1:])

        def rhs(xp: NDArray[np.float64], kp: NDArray[np.float64]) -> tuple[NDArray, NDArray]:
            gp, ginvp, dgp = _sample_metric(g, ginv, dg_inv, xp, origin, spacing)
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
                max_h_rel=max_h_rel,
                notes=("ray left grid",),
            )

    return NullRayResult(
        reached=False,
        t_coord=None,
        t_flat=t_flat,
        max_h_drift=max_h,
        max_h_rel=max_h_rel,
        notes=("max_steps exceeded",),
    )


def _geodesic_report_at_resolution(
    plotfile: str | Path,
    *,
    n: int,
    half_width: float | None,
    n_rays: int,
    h_tol: float,
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
            max_h_rel_drift=max((r.max_h_rel for r in results), default=0.0),
            notes=("no rays reached detector",),
        )

    t_min = min(r.t_coord for r in reached if r.t_coord is not None)
    max_h = max(r.max_h_drift for r in results)
    max_h_rel = max(r.max_h_rel for r in results)
    # Trustworthy when the rays stayed on the null cone in a *relative* sense; the
    # absolute drift floor is set by C0 metric interpolation (see ``H_REL_TOL``).
    h_ok = max_h_rel <= H_REL_TOL
    f_geo = max(0.0, (t_flat - t_min) / t_flat) if t_flat > 0 else 0.0

    notes: list[str] = []
    if not h_ok:
        notes.append(
            f"null constraint drift high (rel H={max_h_rel:.2e}, abs H={max_h:.2e})"
        )

    return GeodesicFtlReport(
        f_geo=f_geo,
        t_min=t_min,
        t_flat=t_flat,
        n_rays=len(results),
        n_reached=len(reached),
        max_h_drift=max_h,
        h_quality_ok=h_ok,
        max_h_rel_drift=max_h_rel,
        notes=tuple(notes),
    )


def compute_geodesic_ftl_from_plotfile(
    plotfile: str | Path,
    *,
    n: int = 65,
    half_width: float | None = None,
    n_rays: int = 5,
    h_tol: float = 1.0e-6,
    refine_n: int | None = GEO_REFINE_N,
) -> GeodesicFtlReport | None:
    """Gauge-invariant null-geodesic shortcut with a reliability re-probe.

    Traces the ray fan at the (cheap) base resolution ``n``.  If that finds a
    coordinate shortcut (``f_geo > GEO_REFINE_FLOOR``) but the reliability gate
    fails -- the classic C0-interpolation discretization artifact, where the
    relative H drift is set by the grid spacing, not real null-cone error -- the
    fan is re-traced once at ``refine_n`` (> base ``n``, default 129^3) before the
    shortcut is discarded.  This stops a genuine sharp-walled warp from being
    silently zeroed at QD resolution while leaving already-reliable or
    no-shortcut measurements untouched (the re-probe never fires for them).
    """
    base = _geodesic_report_at_resolution(
        plotfile, n=n, half_width=half_width, n_rays=n_rays, h_tol=h_tol
    )
    if base is None:
        return None

    needs_refine = (
        refine_n is not None
        and refine_n > n
        and math.isfinite(base.f_geo)
        and base.f_geo > GEO_REFINE_FLOOR
        and not base.h_quality_ok
    )
    if not needs_refine:
        return base

    refined = _geodesic_report_at_resolution(
        plotfile, n=refine_n, half_width=half_width, n_rays=n_rays, h_tol=h_tol
    )
    if refined is None:
        return base

    note = (
        f"reliability re-probe at n={refine_n} "
        f"(base n={n} f_geo={base.f_geo:.3e}, rel-H={base.max_h_rel_drift:.2e} "
        f"-> refined f_geo={refined.f_geo:.3e}, rel-H={refined.max_h_rel_drift:.2e}, "
        f"h_quality_ok={refined.h_quality_ok})"
    )
    return replace(refined, notes=refined.notes + (note,))


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
    ginv = inverse_metric_4d(g)
    g0, ginv0, _ = _sample_metric(g, ginv, dg_inv, x, origin, spacing)
    k = project_null(g0, ginv0, np.array([-1.0, 1.0, 0.0, 0.0]))

    integral = 0.0
    ds = 0.05
    prev_lam = 0.0

    for _ in range(50_000):
        if x[1] >= x_end:
            return integral

        g_pt, ginv_pt, dg_pt = _sample_metric(g, ginv, dg_inv, x, origin, spacing)
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
