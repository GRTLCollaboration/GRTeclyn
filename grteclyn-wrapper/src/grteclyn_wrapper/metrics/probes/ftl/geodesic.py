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

from ..warpfactory import stress_energy
from .geodesic_interp import (
    clamp_index as _clamp_index,
    inverse_metric_4d,
    partial_inverse_metric,
    trilinear as _trilinear,
)
from .metric_field import StaticMetricField


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


# Capture detection thresholds.  A ray that wanders into a collapsed-lapse or
# deep-conformal region (a puncture throat / apparent horizon) has fallen into
# the hole; it is *captured*, not "failed to reach".  Counting these separately
# keeps a black hole on the ray path from silently degrading ``n_reached`` and
# poisoning the trust statistics (PuncturePlan §0.4).  Thresholds are extreme so
# healthy warped regions never misfire: a lapse this small or a conformal factor
# this deep only occurs at a genuine throat/horizon.
CAPTURE_LAPSE_MIN: float = 0.1
CAPTURE_CHI_MIN: float = 1.0e-3


def lapse_from_ginv(ginv_pt: NDArray[np.float64]) -> float:
    """Lapse ``alpha = 1/sqrt(-g^{tt})`` recovered from the inverse metric."""
    gtt = float(ginv_pt[0, 0])
    if gtt >= 0.0:
        return 0.0
    return 1.0 / math.sqrt(-gtt)


def chi_from_g(g_pt: NDArray[np.float64]) -> float:
    """Conformal factor ``chi = det(gamma)^{-1/3}`` (BSSN ``det(h_ij)=1``)."""
    try:
        det = float(np.linalg.det(g_pt[1:, 1:]))
    except np.linalg.LinAlgError:
        return 0.0
    if not math.isfinite(det) or det <= 0.0:
        return 0.0
    return det ** (-1.0 / 3.0)


def ray_is_captured(
    g_pt: NDArray[np.float64],
    ginv_pt: NDArray[np.float64],
    *,
    chi_min: float = CAPTURE_CHI_MIN,
    lapse_min: float = CAPTURE_LAPSE_MIN,
) -> bool:
    """True when the sampled point lies inside a puncture throat / horizon."""
    return lapse_from_ginv(ginv_pt) < lapse_min or chi_from_g(g_pt) < chi_min


@dataclass(frozen=True)
class NullRayResult:
    """Outcome of tracing one null ray."""

    reached: bool
    t_coord: float | None
    t_flat: float
    max_h_drift: float
    max_h_rel: float = 0.0
    notes: tuple[str, ...] = ()
    # True when the ray fell into a puncture throat / horizon (excluded from
    # f_geo, counted separately from plain non-arrivals).
    captured: bool = False


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
    n_captured: int = 0


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


def build_metric_3d_from_gridinit(
    gridinit_path: str | Path,
    *,
    n: int | None = None,
    half_width: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], tuple[float, float, float]]:
    """Sample a static 4-metric ``g_{mu nu}(x,y,z)`` from a ``.gridinit`` file.

    Returns ``(g, origin, spacing)`` where ``g`` has shape ``(nx,ny,nz,4,4)``,
    ``origin`` is the coordinate of the ``(0,0,0)`` sample (first cell centre),
    and ``spacing`` is ``(dx,dy,dz)``.

    The on-disk gridinit layout is ``(nz, ny, nx, n_comp)``; this helper
    transposes to the ``(nx, ny, nz)`` convention used by the geodesic probes.
    Optional ``n`` / ``half_width`` subsample a centred window (nearest-neighbour
    index selection) for cheaper probes on large grids.
    """
    from ....grtresna.io.gridinit import read_gridinit

    grid = read_gridinit(gridinit_path)
    # grid.data: (nz, ny, nx, n_comp) → component lookup
    name_to_idx = {name: i for i, name in enumerate(grid.comp_names)}

    def comp(name: str) -> NDArray[np.float64]:
        # Return (nx, ny, nz)
        arr = np.asarray(grid.data[..., name_to_idx[name]], dtype=float)
        return np.transpose(arr, (2, 1, 0))

    chi = np.clip(comp("chi"), 1.0e-10, None)
    inv_chi = 1.0 / chi
    alpha = np.clip(comp("lapse"), 1.0e-10, None)
    beta_up = np.stack([comp("shift1"), comp("shift2"), comp("shift3")], axis=-1)

    gamma = np.zeros(chi.shape + (3, 3), dtype=float)
    gamma[..., 0, 0] = comp("h11") * inv_chi
    gamma[..., 0, 1] = gamma[..., 1, 0] = comp("h12") * inv_chi
    gamma[..., 0, 2] = gamma[..., 2, 0] = comp("h13") * inv_chi
    gamma[..., 1, 1] = comp("h22") * inv_chi
    gamma[..., 1, 2] = gamma[..., 2, 1] = comp("h23") * inv_chi
    gamma[..., 2, 2] = comp("h33") * inv_chi

    beta_low = np.einsum("...ij,...j->...i", gamma, beta_up)
    beta_sq = np.einsum("...i,...i->...", beta_low, beta_up)

    g = np.zeros(chi.shape + (4, 4), dtype=float)
    g[..., 0, 0] = beta_sq - alpha * alpha
    g[..., 0, 1:] = beta_low
    g[..., 1:, 0] = beta_low
    g[..., 1:, 1:] = gamma

    # Cell-centre coordinates of the full grid.
    dx_xyz = np.asarray(grid.dx_xyz, dtype=float)
    origin_full = np.asarray(grid.origin, dtype=float) + 0.5 * dx_xyz
    nx, ny, nz = g.shape[:3]

    if n is None and half_width is None:
        return g, origin_full, (float(dx_xyz[0]), float(dx_xyz[1]), float(dx_xyz[2]))

    # Subsample a centred window.
    center = origin_full + 0.5 * (np.array([nx, ny, nz]) - 1.0) * dx_xyz
    if half_width is None:
        half_width = 0.45 * float(nx * dx_xyz[0])
    if n is None:
        n = min(nx, ny, nz)

    # Index range covering [center - half_width, center + half_width].
    def _slice(axis: int) -> slice:
        i0 = int(np.floor(((center[axis] - half_width) - origin_full[axis]) / dx_xyz[axis]))
        i1 = int(np.ceil(((center[axis] + half_width) - origin_full[axis]) / dx_xyz[axis]))
        i0 = max(0, i0)
        i1 = min([nx, ny, nz][axis] - 1, i1)
        return slice(i0, i1 + 1)

    sx, sy, sz = _slice(0), _slice(1), _slice(2)
    g_win = g[sx, sy, sz]
    origin_win = np.array(
        [
            origin_full[0] + sx.start * dx_xyz[0],
            origin_full[1] + sy.start * dx_xyz[1],
            origin_full[2] + sz.start * dx_xyz[2],
        ],
        dtype=float,
    )
    # Optional downsample to roughly n^3.
    step = max(1, int(np.floor(max(g_win.shape[:3]) / max(n, 1))))
    g_ds = g_win[::step, ::step, ::step]
    spacing = (
        float(dx_xyz[0] * step),
        float(dx_xyz[1] * step),
        float(dx_xyz[2] * step),
    )
    return g_ds, origin_win, spacing


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
    axis: int = 0,
    max_steps: int = 50_000,
    ds_init: float = 0.05,
    h_tol: float = 1.0e-6,
) -> NullRayResult:
    """Trace one null ray from ``x_start`` to ``x_end`` at fixed ``(y,z)``."""
    from .evolving_geodesic import integrate_null_ray_on_field

    field = StaticMetricField(g=g, origin=origin, spatial_spacing=tuple(spacing))
    return integrate_null_ray_on_field(
        field,
        x_start=x_start,
        x_end=x_end,
        y0=y0,
        z0=z0,
        t0=t0,
        axis=axis,
        max_steps=max_steps,
        ds_init=ds_init,
        h_tol=h_tol,
    )


_AXIS_LABELS = ("x", "y", "z")


def propagation_endpoints(
    origin: NDArray[np.float64],
    spacing: Sequence[float],
    shape: Sequence[int],
    *,
    axis: int = 0,
    half_length: float | None = None,
    margin_frac: float = 0.05,
) -> tuple[float, float]:
    """Return ``(prop_start, prop_end)`` along ``axis``.

    Default behaviour matches the historical 5%–95% box span.  When
    ``half_length`` is set, endpoints are placed symmetrically about the
    domain centre at ``± half_length``, clamped inside the margin band so
    short-support geometries are not diluted by long flat path segments.
    """
    lo = float(origin[axis]) + margin_frac * (shape[axis] - 1) * spacing[axis]
    hi = float(origin[axis]) + (1.0 - margin_frac) * (shape[axis] - 1) * spacing[axis]
    if half_length is None:
        return lo, hi
    center = float(origin[axis]) + 0.5 * (shape[axis] - 1) * spacing[axis]
    half = max(float(half_length), 2.0 * float(spacing[axis]))
    start = max(lo, center - half)
    end = min(hi, center + half)
    if end <= start:
        return lo, hi
    return start, end


def geodesic_report_from_metric_g(
    g: NDArray[np.float64],
    origin: NDArray[np.float64],
    spacing: Sequence[float],
    *,
    n_rays: int,
    h_tol: float,
    axis: int = 0,
    half_length: float | None = None,
) -> GeodesicFtlReport:
    """Run a fan of null rays on a pre-sampled static 4-metric grid.

    ``axis`` selects the propagation direction (0=x, 1=y, 2=z).  Transverse
    ray offsets are applied along the higher of the two remaining spatial axes.
    Optional ``half_length`` localises the emitter/detector about the domain
    centre (see :func:`propagation_endpoints`).
    """
    dg_inv = partial_inverse_metric(g, spacing)
    shape = g.shape[:3]
    transverse = [i for i in range(3) if i != axis]
    fix_spatial, fan_spatial = transverse  # ascending order

    prop_start, prop_end = propagation_endpoints(
        origin, spacing, shape, axis=axis, half_length=half_length
    )
    t_flat = prop_end - prop_start

    y0_center = origin[fix_spatial] + 0.5 * (shape[fix_spatial] - 1) * spacing[fix_spatial]
    z0_center = origin[fan_spatial] + 0.5 * (shape[fan_spatial] - 1) * spacing[fan_spatial]
    offsets = np.linspace(-0.1, 0.1, max(1, n_rays)) * (shape[fan_spatial] - 1) * spacing[fan_spatial]

    results: list[NullRayResult] = []
    for dz in offsets:
        res = integrate_null_ray(
            g,
            dg_inv,
            origin,
            spacing,
            x_start=prop_start,
            x_end=prop_end,
            y0=y0_center,
            z0=z0_center + float(dz),
            axis=axis,
            h_tol=h_tol,
        )
        results.append(res)

    reached = [r for r in results if r.reached and r.t_coord is not None]
    n_captured = sum(1 for r in results if r.captured)
    if not reached:
        note = f"no rays reached detector (axis={_AXIS_LABELS[axis]})"
        if n_captured:
            note += f"; {n_captured}/{len(results)} captured by puncture/horizon"
        return GeodesicFtlReport(
            f_geo=0.0,
            t_min=None,
            t_flat=t_flat,
            n_rays=len(results),
            n_reached=0,
            max_h_drift=max((r.max_h_drift for r in results), default=0.0),
            h_quality_ok=False,
            max_h_rel_drift=max((r.max_h_rel for r in results), default=0.0),
            notes=(note,),
            n_captured=n_captured,
        )

    t_min = min(r.t_coord for r in reached if r.t_coord is not None)
    max_h = max(r.max_h_drift for r in results)
    max_h_rel = max(r.max_h_rel for r in results)
    h_ok = max_h_rel <= H_REL_TOL
    f_geo = max(0.0, (t_flat - t_min) / t_flat) if t_flat > 0 else 0.0

    notes: list[str] = []
    if not h_ok:
        notes.append(
            f"null constraint drift high (rel H={max_h_rel:.2e}, abs H={max_h:.2e})"
        )
    if n_captured:
        notes.append(
            f"{n_captured}/{len(results)} rays captured by puncture/horizon "
            "(excluded from f_geo)"
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
        n_captured=n_captured,
    )


def _rays_complete(report: GeodesicFtlReport) -> bool:
    """True when every non-captured ray reached the detector."""
    return report.n_reached > 0 and report.n_reached == report.n_rays - report.n_captured


def _frozen_report_score(report: GeodesicFtlReport) -> float:
    """Comparable quality metric for choosing best direction.

    Rank primarily by raw ``f_geo`` so a high-drift warp axis is not discarded
    in favour of a flat transverse axis.  A tiny bonus prefers h-quality when
    shortcuts are otherwise equal.
    """
    if not _rays_complete(report):
        return -1.0
    return float(report.f_geo) + (1.0e-6 if report.h_quality_ok else 0.0)


def geodesic_report_best_direction(
    g: NDArray[np.float64],
    origin: NDArray[np.float64],
    spacing: Sequence[float],
    *,
    n_rays: int,
    h_tol: float,
    directions: Sequence[str] = ("x",),
    half_length: float | None = None,
) -> GeodesicFtlReport:
    """Scan principal axes and return the report with the best ``f_geo``."""
    axis_map = {label: idx for idx, label in enumerate(_AXIS_LABELS)}
    axes = [axis_map[d] for d in directions if d in axis_map]
    if not axes:
        axes = [0]

    reports = [
        geodesic_report_from_metric_g(
            g,
            origin,
            spacing,
            n_rays=n_rays,
            h_tol=h_tol,
            axis=ax,
            half_length=half_length,
        )
        for ax in axes
    ]

    best = reports[0]
    best_score = _frozen_report_score(best)
    best_axis_label = _AXIS_LABELS[axes[0]]
    for ax, rep in zip(axes, reports):
        score = _frozen_report_score(rep)
        if score > best_score:
            best_score = score
            best = rep
            best_axis_label = _AXIS_LABELS[ax]

    if len(axes) > 1:
        extra_notes = best.notes + (f"best_direction={best_axis_label}",)
        return replace(best, notes=extra_notes)
    return best


def _geodesic_report_at_resolution(
    plotfile: str | Path,
    *,
    n: int,
    half_width: float | None,
    n_rays: int,
    h_tol: float,
    directions: Sequence[str] = ("x",),
) -> GeodesicFtlReport | None:
    """Run a fan of null rays and return ``f_geo``, scanning ``directions``."""
    try:
        g, origin, spacing = build_metric_3d_from_plotfile(
            plotfile, n=n, half_width=half_width
        )
    except Exception:
        return None
    if len(directions) > 1:
        return geodesic_report_best_direction(
            g, origin, spacing, n_rays=n_rays, h_tol=h_tol, directions=directions
        )
    axis = _AXIS_LABELS.index(directions[0]) if directions[0] in _AXIS_LABELS else 0
    return geodesic_report_from_metric_g(
        g, origin, spacing, n_rays=n_rays, h_tol=h_tol, axis=axis
    )


def compute_geodesic_ftl_from_plotfile(
    plotfile: str | Path,
    *,
    n: int = 65,
    half_width: float | None = None,
    n_rays: int = 5,
    h_tol: float = 1.0e-6,
    refine_n: int | None = GEO_REFINE_N,
    directions: Sequence[str] = ("x",),
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
        plotfile, n=n, half_width=half_width, n_rays=n_rays, h_tol=h_tol,
        directions=directions,
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
        plotfile, n=refine_n, half_width=half_width, n_rays=n_rays, h_tol=h_tol,
        directions=directions,
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
