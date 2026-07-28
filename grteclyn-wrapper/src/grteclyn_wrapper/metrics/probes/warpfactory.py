"""Warp Factory-style energy-condition evaluation from a 4-metric.

This is a numpy port of the core methodology of Helmerich et al.,
"Analyzing Warp Drive Spacetimes with Warp Factory" (arXiv:2404.03095):

1. Specify the covariant 4-metric ``g_{mu nu}`` on a Cartesian ``(t,x,y,z)``
   grid.
2. Compute the Einstein tensor ``G_{mu nu}`` numerically (finite-difference
   Christoffel -> Ricci -> Einstein) and hence the stress-energy tensor
   ``T_{mu nu} = G_{mu nu} / 8 pi`` (geometrized units, G=c=1).
3. Build a local orthonormal observer frame from the 3+1 (ADM) split of the
   metric, sample a sphere of null and timelike observers, and contract the
   stress-energy tensor with each to evaluate the pointwise energy conditions
   (NEC, WEC, SEC, DEC), taking the minimum over observers at every point.

Unlike the recipe ``physical_metrics`` proxies (which only see the
constraint-derived energy density ``rho`` at ``t=0``), this evaluator needs a
*full* 4-metric with genuine time dependence, so it applies to:

* the analytic seed metrics (Alcubierre, Minkowski, ...) built here, and
* evolved ``GRTeclyn`` plotfile data reassembled into ``g_{mu nu}`` (future).

For a warp bubble moving at ``v>0`` it reproduces the canonical result that the
NEC/WEC are violated (negative energy density seen by some observer), while
flat Minkowski satisfies every condition to numerical precision.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

EIGHT_PI = 8.0 * math.pi
_ETA_TIME = -1.0  # signature (-,+,+,+)


@dataclass(frozen=True)
class EnergyConditionReport:
    """Energy-condition diagnostics over a grid.

    Two complementary verification pathways are reported (cf. Le, arXiv:2602.18023):

    * **Observer sampling** (``*_min``): the minimum of ``T_{mu nu} v^mu v^nu``
      over a discrete sphere of null/timelike observers -- a lower bound on the
      true worst case (Warp Factory style, Helmerich et al. arXiv:2404.03095).
    * **Hawking--Ellis eigenvalue test** (``*_slack_min``): at Type I
      stress-energy points the energy conditions reduce to *exact*, observer-
      independent algebraic inequalities in the eigenvalues ``(-rho, p_i)`` of
      the mixed tensor ``T^a_b``.  This is authoritative where it applies.

    A *negative* value signals a pointwise violation (exotic matter).
    """

    # Observer-sampling pathway (lower bounds).
    nec_min: float
    wec_min: float
    sec_min: float
    dec_min: float
    nec_violation_fraction: float
    wec_violation_fraction: float
    sec_violation_fraction: float
    dec_violation_fraction: float
    rho_min: float
    rho_eulerian_min: float
    n_observers: int
    n_points: int
    s_energy_conditions: float
    # Hawking--Ellis eigenvalue pathway (exact at Type I points).
    type_I_fraction: float = 1.0
    nec_slack_min: float | None = None
    wec_slack_min: float | None = None
    sec_slack_min: float | None = None
    dec_slack_min: float | None = None
    rho_invariant_min: float | None = None
    notes: tuple[str, ...] = ()


def fibonacci_sphere(n: int) -> NDArray[np.float64]:
    """``n`` approximately-uniform unit vectors on S^2 (shape ``(n, 3)``)."""
    if n < 1:
        raise ValueError("n must be >= 1")
    if n == 1:
        return np.array([[1.0, 0.0, 0.0]])
    i = np.arange(n, dtype=np.float64)
    phi = math.pi * (3.0 - math.sqrt(5.0))  # golden angle
    z = 1.0 - 2.0 * i / (n - 1)
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, 1.0))
    theta = phi * i
    return np.stack([r * np.cos(theta), r * np.sin(theta), z], axis=1)


def _bounded_reward(value: float, scale: float) -> float:
    if not math.isfinite(value) or value < 0:
        return 0.0
    return 1.0 / (1.0 + value / scale)


# --------------------------------------------------------------------------
# Curvature: covariant 4-metric grid -> stress-energy tensor.
# --------------------------------------------------------------------------

def _d_dx(field: NDArray, dx: float, axis: int) -> NDArray:
    """Fourth-order central derivative along ``axis`` (2nd order at the edges).

    Warp Factory and Lousto (arXiv:gr-qc/0503001) both stress that high-order
    convergence is essential; second-order differences inject enough noise into
    the curvature chain to spuriously complexify the stress-energy eigenvalues.
    Interior points use the 5-point stencil
    ``(-f_{i+2}+8 f_{i+1}-8 f_{i-1}+f_{i-2})/(12 dx)``; the two boundary layers
    (which callers crop) fall back to ``np.gradient``.
    """
    n = field.shape[axis]
    out = np.gradient(field, dx, axis=axis, edge_order=2)
    if n >= 5:
        def sl(start: int, stop: int) -> tuple:
            idx = [slice(None)] * field.ndim
            idx[axis] = slice(start, stop)
            return tuple(idx)

        interior = (
            -field[sl(4, n)] + 8.0 * field[sl(3, n - 1)]
            - 8.0 * field[sl(1, n - 3)] + field[sl(0, n - 4)]
        ) / (12.0 * dx)
        out[sl(2, n - 2)] = interior
    return out


def _partial(field: NDArray, spacing: Sequence[float]) -> NDArray:
    """Stack of partial derivatives along the 4 grid axes.

    ``field`` has shape ``S + tail`` with ``S=(Nt,Nx,Ny,Nz)``.  Returns an
    array of shape ``S + tail + (4,)`` whose last axis indexes the derivative
    direction ``mu``.
    """
    grads = [_d_dx(field, spacing[mu], axis=mu) for mu in range(4)]
    return np.stack(grads, axis=-1)


def einstein_tensor(g: NDArray, spacing: Sequence[float]) -> NDArray:
    """Einstein tensor ``G_{mu nu}`` from a covariant 4-metric grid.

    ``g`` has shape ``(Nt,Nx,Ny,Nz,4,4)``; ``spacing`` is ``(dt,dx,dy,dz)``.
    Uses second-order central differences (lower order at the boundary, which
    callers should crop -- see :func:`evaluate_four_metric`).
    """
    ginv = np.linalg.inv(g)

    # dg[..., i, j, mu] = d_mu g_ij
    dg = _partial(g, spacing)

    # Lowered-index Christoffel: Gamma_{d b c} = 1/2 (d_b g_dc + d_c g_db - d_d g_bc)
    d_b_g_dc = np.einsum("...dcb->...dbc", dg)   # d_b g_{dc}
    d_c_g_db = np.einsum("...dbc->...dbc", dg)   # d_c g_{db}
    d_d_g_bc = np.einsum("...bcd->...dbc", dg)   # d_d g_{bc}
    gamma_low = 0.5 * (d_b_g_dc + d_c_g_db - d_d_g_bc)

    # Gamma^a_{bc}
    gamma = np.einsum("...ad,...dbc->...abc", ginv, gamma_low)

    # dGamma[..., a, b, c, mu] = d_mu Gamma^a_{bc}
    dgamma = _partial(gamma, spacing)

    # Ricci R_{s n} = d_r Gamma^r_{n s} - d_n Gamma^r_{r s}
    #                 + Gamma^r_{r l} Gamma^l_{n s} - Gamma^r_{n l} Gamma^l_{r s}
    term1 = np.einsum("...rnsr->...sn", dgamma)  # d_r Gamma^r_{ns}
    term2 = np.einsum("...rrsn->...sn", dgamma)  # d_n Gamma^r_{rs}
    term3 = np.einsum("...rrl,...lns->...sn", gamma, gamma)
    term4 = np.einsum("...rnl,...lrs->...sn", gamma, gamma)
    ricci = term1 - term2 + term3 - term4

    ricci_scalar = np.einsum("...sn,...sn->...", ginv, ricci)
    g_einstein = ricci - 0.5 * ricci_scalar[..., None, None] * g
    return g_einstein


def stress_energy(g: NDArray, spacing: Sequence[float]) -> NDArray:
    """Stress-energy tensor ``T_{mu nu} = G_{mu nu} / (8 pi)``."""
    return einstein_tensor(g, spacing) / EIGHT_PI


# --------------------------------------------------------------------------
# 3+1 (ADM) split and local orthonormal observer frame.
# --------------------------------------------------------------------------

def decompose_adm(g: NDArray) -> tuple[NDArray, NDArray, NDArray]:
    """ADM lapse/shift/spatial-metric from a covariant 4-metric.

    Returns ``(alpha, beta_up, gamma)`` with shapes ``S``, ``S+(3,)`` and
    ``S+(3,3)`` respectively.
    """
    gamma = g[..., 1:, 1:]
    beta_low = g[..., 0, 1:]
    gamma_inv = np.linalg.inv(gamma)
    beta_up = np.einsum("...ij,...j->...i", gamma_inv, beta_low)
    beta_sq = np.einsum("...i,...i->...", beta_low, beta_up)
    alpha = np.sqrt(np.clip(beta_sq - g[..., 0, 0], 1.0e-12, None))
    return alpha, beta_up, gamma


def _spatial_triad(gamma: NDArray) -> NDArray:
    """Orthonormal spatial triad ``e_(a)^i`` (shape ``S+(3,3)``, column=a).

    ``gamma_ij e_(a)^i e_(b)^j = delta_ab`` via a Cholesky factorization.
    """
    L = np.linalg.cholesky(gamma)
    # E = L^{-T} satisfies E^T gamma E = I.
    eye = np.broadcast_to(np.eye(3), gamma.shape).copy()
    L_inv = np.linalg.solve(L, eye)
    return np.einsum("...ai->...ia", L_inv)  # transpose -> columns are e_(a)


def hawking_ellis_slacks(
    T: NDArray,
    g: NDArray,
    *,
    real_tol: float = 1.0e-2,
    vacuum_rel: float = 1.0e-2,
) -> dict[str, NDArray]:
    """Exact Type I energy-condition slacks from the eigenvalues of ``T^a_b``.

    Following Le (arXiv:2602.18023, Sec. 2.2), at a Hawking--Ellis Type I point
    the mixed stress-energy tensor ``T^a_b = g^{ac} T_{cb}`` has four real
    eigenvalues ``(-rho, p_1, p_2, p_3)`` with one timelike eigenvector.  The
    energy conditions then reduce to algebraic inequalities whose *slack* (the
    signed distance to the violation boundary) is observer-independent::

        NEC: min_i (rho + p_i)
        WEC: min(rho, min_i (rho + p_i))
        SEC: min( min_i (rho + p_i), rho + sum_i p_i )
        DEC: min_i (rho - |p_i|)

    Returns per-point arrays plus a boolean ``type_I`` mask.  At non-Type-I
    points the slacks are ``nan`` (the algebraic test does not apply there).
    """
    ginv = np.linalg.inv(g)
    T_mixed = np.einsum("...ac,...cb->...ab", ginv, T)

    w, V = np.linalg.eig(T_mixed)  # w: S+(4,), V columns are eigenvectors
    re = w.real
    im = w.imag

    # Tolerances scale with the *global* stress-energy magnitude: in a warp
    # spacetime the matter is concentrated on the bubble wall, and the
    # near-flat bulk (|T| << peak) is effectively vacuum -- its eigenvectors are
    # numerically arbitrary, so classifying it by causal character is
    # meaningless (cf. Le, arXiv:2602.18023, Sec. 2.2).
    point_scale = np.max(np.abs(re), axis=-1)
    global_scale = float(point_scale.max()) if point_scale.size else 0.0
    vac_floor = vacuum_rel * global_scale if global_scale > 0 else 1.0e-30
    near_vacuum = point_scale < vac_floor

    denom = np.maximum(point_scale, vac_floor)
    all_real = np.all(np.abs(im) < real_tol * denom[..., None], axis=-1)

    # Causal character of each eigenvector: g_ab v^a v^b (sign is scale-free).
    gnorm = np.einsum("...ab,...ak,...bk->...k", g.astype(complex), V, V).real
    n_timelike = np.sum(gnorm < 0.0, axis=-1)

    t_idx = np.argmin(gnorm, axis=-1)
    idx = np.arange(4)
    t_onehot = idx == t_idx[..., None]
    spacelike = ~t_onehot

    rho = -np.sum(np.where(t_onehot, re, 0.0), axis=-1)
    rpp = np.where(spacelike, rho[..., None] + re, np.inf)
    nec = rpp.min(axis=-1)
    wec = np.minimum(rho, nec)
    sum_p = np.sum(np.where(spacelike, re, 0.0), axis=-1)
    sec = np.minimum(nec, rho + sum_p)
    dec = np.where(spacelike, rho[..., None] - np.abs(re), np.inf).min(axis=-1)

    type_I = all_real & (n_timelike >= 1)
    type_I = type_I | near_vacuum

    # Near-vacuum points trivially satisfy every condition (slack 0).
    for arr in (rho, nec, wec, sec, dec):
        arr[near_vacuum] = 0.0

    mask = type_I
    def masked(arr: NDArray) -> NDArray:
        out = np.where(mask, arr, np.nan)
        return out

    return {
        "type_I": type_I,
        "rho": masked(rho),
        "nec": masked(nec),
        "wec": masked(wec),
        "sec": masked(sec),
        "dec": masked(dec),
    }


def energy_conditions(
    T: NDArray,
    g: NDArray,
    *,
    n_directions: int = 50,
    n_speeds: int = 4,
    max_speed: float = 0.9,
    tol: float = 1.0e-8,
    score_scale: float = 1.0e-3,
) -> EnergyConditionReport:
    """Evaluate NEC/WEC/SEC/DEC over a sphere of observers at each point.

    ``T`` and ``g`` are covariant tensors of shape ``S+(4,4)``.  Observers are
    built from the local ADM frame, so contractions are exact regardless of
    coordinates.
    """
    S = g.shape[:-2]
    n_points = int(np.prod(S)) if S else 1

    ginv = np.linalg.inv(g)
    alpha, beta_up, gamma = decompose_adm(g)
    triad = _spatial_triad(gamma)  # S+(3,3), column a -> e_(a)^i

    # Eulerian normal n^mu = (1/alpha)(1, -beta^i)
    n_up = np.empty(S + (4,))
    n_up[..., 0] = 1.0 / alpha
    n_up[..., 1:] = -beta_up / alpha[..., None]

    # Spatial tetrad 4-vectors s_(a)^mu = (0, e_(a)^i)
    s_up = np.zeros(S + (3, 4))  # [..., a, mu]
    s_up[..., 1:] = np.einsum("...ia->...ai", triad)

    directions = fibonacci_sphere(n_directions)  # (D,3)
    D = directions.shape[0]

    # Spatial direction 4-vectors d^mu = d^a s_(a)^mu  -> shape S+(D,4)
    d_up = np.einsum("da,...am->...dm", directions, s_up)

    T_scalar = np.einsum("...mn,...mn->...", ginv, T)  # trace T

    def quad(vec: NDArray) -> NDArray:
        """Contract T_{mu nu} v^mu v^nu over trailing observer axis."""
        return np.einsum("...mn,...km,...kn->...k", T, vec, vec)

    # --- NEC: null vectors k^mu = n^mu + d^mu ---
    k_up = n_up[..., None, :] + d_up  # S+(D,4)
    nec = quad(k_up)  # S+(D,)
    nec_min_pt = nec.min(axis=-1)

    # --- timelike observers u^mu = gamma_L (n^mu + speed * d_hat^mu) ---
    speeds = np.linspace(0.0, max_speed, n_speeds + 1)[1:] if n_speeds > 0 else np.array([0.0])
    wec_min_pt = np.full(S, np.inf)
    sec_min_pt = np.full(S, np.inf)
    dec_min_pt = np.full(S, np.inf)
    # Eulerian (at-rest) energy density rho = T_{mn} n^m n^n
    rho_eulerian = np.einsum("...mn,...m,...n->...", T, n_up, n_up)

    for speed in speeds:
        lorentz = 1.0 / math.sqrt(max(1.0 - speed * speed, 1.0e-12))
        u_up = lorentz * (n_up[..., None, :] + speed * d_up)  # S+(D,4)

        wec = quad(u_up)
        wec_min_pt = np.minimum(wec_min_pt, wec.min(axis=-1))

        # SEC: (T_{mn} - 1/2 T g_{mn}) u^m u^n
        sec = wec - 0.5 * T_scalar[..., None] * np.einsum(
            "...mn,...km,...kn->...k", g, u_up, u_up
        )
        sec_min_pt = np.minimum(sec_min_pt, sec.min(axis=-1))

        # DEC: Upsilon^mu = -T^mu_nu u^nu must be causal & future pointing.
        T_mixed = np.einsum("...ma,...an->...mn", ginv, T)  # T^m_n
        ups = -np.einsum("...mn,...kn->...km", T_mixed, u_up)  # S+(D,4)
        ups_norm = np.einsum("...mn,...km,...kn->...k", g, ups, ups)
        xi = -np.sign(ups_norm) * np.sqrt(np.abs(ups_norm))
        dec_min_pt = np.minimum(dec_min_pt, xi.min(axis=-1))

    def vfrac(arr: NDArray) -> float:
        return float(np.mean(arr < -tol))

    nec_min = float(nec_min_pt.min())
    wec_min = float(wec_min_pt.min())
    sec_min = float(sec_min_pt.min())
    dec_min = float(dec_min_pt.min())

    # Exact Hawking-Ellis eigenvalue pathway (authoritative at Type I points).
    he = hawking_ellis_slacks(T, g)
    type_I_fraction = float(np.mean(he["type_I"]))

    def nanmin(arr: NDArray) -> float | None:
        vals = arr[np.isfinite(arr)]
        return float(vals.min()) if vals.size else None

    nec_slack_min = nanmin(he["nec"])
    wec_slack_min = nanmin(he["wec"])
    sec_slack_min = nanmin(he["sec"])
    dec_slack_min = nanmin(he["dec"])
    rho_invariant_min = nanmin(he["rho"])

    notes: list[str] = []
    if nec_slack_min is not None and nec_slack_min < -tol:
        notes.append(
            f"NEC violated exactly (Type I slack min={nec_slack_min:.3e}); exotic matter required"
        )
    elif nec_min < -tol:
        notes.append(f"NEC violated (observer min={nec_min:.3e})")
    if wec_slack_min is not None and wec_slack_min < -tol:
        notes.append(f"WEC violated exactly (Type I slack min={wec_slack_min:.3e})")
    if type_I_fraction < 1.0:
        notes.append(
            f"{100.0 * (1.0 - type_I_fraction):.1f}% of points are non-Type-I; "
            "eigenvalue slacks omit those (observer bounds still apply)"
        )

    # Hybrid margin (cf. Le, arXiv:2602.18023): the most-violating of the exact
    # Type I slacks and the observer-sampled lower bounds, so a violation seen by
    # either pathway is never falsely certified as clean.
    candidates = [nec_min, wec_min, sec_min]
    candidates += [v for v in (nec_slack_min, wec_slack_min, sec_slack_min) if v is not None]
    worst = min(candidates)
    s_energy = _bounded_reward(max(0.0, -worst), score_scale)

    return EnergyConditionReport(
        nec_min=nec_min,
        wec_min=wec_min,
        sec_min=sec_min,
        dec_min=dec_min,
        nec_violation_fraction=vfrac(nec_min_pt),
        wec_violation_fraction=vfrac(wec_min_pt),
        sec_violation_fraction=vfrac(sec_min_pt),
        dec_violation_fraction=vfrac(dec_min_pt),
        rho_min=float(nec_min),
        rho_eulerian_min=float(rho_eulerian.min()),
        n_observers=int(D * len(speeds)),
        n_points=n_points,
        s_energy_conditions=s_energy,
        type_I_fraction=type_I_fraction,
        nec_slack_min=nec_slack_min,
        wec_slack_min=wec_slack_min,
        sec_slack_min=sec_slack_min,
        dec_slack_min=dec_slack_min,
        rho_invariant_min=rho_invariant_min,
        notes=tuple(notes),
    )


@dataclass(frozen=True)
class ExoticEnergyBudget:
    """Negative-energy ("exotic matter") content of a 4-metric, on one measure.

    All energies are proper-volume integrals of the Eulerian energy density
    ``rho = T_{mu nu} n^mu n^nu`` (``n`` = future unit normal to the ``t=const``
    slice) over the central time slice, in geometric units (G=c=1)::

        total_negative_energy = int_{rho<0} rho sqrt(det gamma) d^3x   (<= 0)

    ``exotic_energy`` (its magnitude) is the headline "how much exotic matter"
    number.  Because both an analytic Alcubierre bubble and an evolved candidate
    (via :func:`build_four_metric_stack_from_plotfiles`) reach this through the
    same :func:`stress_energy` + :func:`decompose_adm`, the measure is identical
    on both sides and the comparison is apples-to-apples.
    """

    total_negative_energy: float
    total_positive_energy: float
    min_rho: float
    max_rho: float
    negative_fraction: float
    n_points: int

    @property
    def exotic_energy(self) -> float:
        """Magnitude of the negative-energy integral (>= 0)."""
        return abs(self.total_negative_energy)

    @property
    def net_energy(self) -> float:
        return self.total_negative_energy + self.total_positive_energy


def eulerian_energy_density(g: NDArray, spacing: Sequence[float]) -> NDArray:
    """Eulerian energy density ``rho = T_{mu nu} n^mu n^nu`` on every grid point.

    ``n^mu = (1, -beta^i) / alpha`` is the future-pointing unit normal to the
    ``t=const`` slices, taken from the ADM split of ``g`` (:func:`decompose_adm`).
    Returns an array of shape ``g.shape[:-2]`` (the spacetime grid).
    """
    T = stress_energy(g, spacing)
    alpha, beta_up, _gamma = decompose_adm(g)
    n = np.zeros(g.shape[:-1], dtype=float)  # grid + (4,)
    n[..., 0] = 1.0 / alpha
    n[..., 1:] = -beta_up / alpha[..., None]
    return np.einsum("...ab,...a,...b->...", T, n, n)


def exotic_energy_budget(
    g: NDArray,
    spacing: Sequence[float],
    *,
    crop: int = 4,
    time_index: int | None = None,
) -> ExoticEnergyBudget:
    """Negative-energy ("exotic matter") budget of a 4-metric grid.

    Mirrors :func:`evaluate_four_metric`: evaluates the central time slice (or
    ``time_index`` if given) and drops ``crop`` finite-difference boundary cells
    on every spatial axis.  The Eulerian energy density is integrated against the
    proper volume element ``sqrt(det gamma) d^3x``.
    """
    rho_all = eulerian_energy_density(g, spacing)
    _alpha, _beta, gamma = decompose_adm(g)
    sqrt_gamma = np.sqrt(np.clip(np.linalg.det(gamma), 0.0, None))

    nt = g.shape[0]
    tc = nt // 2 if time_index is None else int(time_index)
    spatial = tuple(slice(crop, dim - crop) for dim in g.shape[1:4])
    sl = (tc,) + spatial
    rho = rho_all[sl]
    sg = sqrt_gamma[sl]
    if rho.size == 0:
        raise ValueError("grid too small after cropping; increase resolution or lower crop")

    spatial_spacing = np.asarray(spacing, dtype=float)[-3:]
    dV = float(np.prod(spatial_spacing))
    cell = sg * dV

    neg = rho < 0.0
    pos = rho > 0.0
    return ExoticEnergyBudget(
        total_negative_energy=float(np.sum(rho[neg] * cell[neg])),
        total_positive_energy=float(np.sum(rho[pos] * cell[pos])),
        min_rho=float(rho.min()),
        max_rho=float(rho.max()),
        negative_fraction=float(neg.mean()),
        n_points=int(rho.size),
    )


def evaluate_four_metric(
    g: NDArray,
    spacing: Sequence[float],
    *,
    crop: int = 4,
    **ec_kwargs,
) -> EnergyConditionReport:
    """Full pipeline: 4-metric grid -> T_{mu nu} -> energy conditions.

    The boundary ``crop`` cells along every axis are discarded (finite-
    difference stencils are inaccurate there), and only the central time slice
    is evaluated.
    """
    T = stress_energy(g, spacing)

    nt = g.shape[0]
    tc = nt // 2
    sl = (slice(tc, tc + 1),) + tuple(
        slice(crop, dim - crop) for dim in g.shape[1:4]
    )
    g_c = g[sl]
    T_c = T[sl]
    if min(g_c.shape[:4]) < 1:
        raise ValueError("grid too small after cropping; increase resolution")
    return energy_conditions(T_c, g_c, **ec_kwargs)


# --------------------------------------------------------------------------
# Analytic 4-metric builders on a (t,x,y,z) grid.
# --------------------------------------------------------------------------

def _grid(
    *,
    half_width: float,
    n_space: int,
    dt: float,
    n_time: int = 5,
) -> tuple[NDArray, NDArray, NDArray, NDArray, tuple[float, float, float, float]]:
    x = np.linspace(-half_width, half_width, n_space)
    dx = x[1] - x[0]
    t = (np.arange(n_time) - n_time // 2) * dt
    T, X, Y, Z = np.meshgrid(t, x, x, x, indexing="ij")
    return T, X, Y, Z, (dt, dx, dx, dx)


def minkowski_metric(half_width: float = 4.0, n_space: int = 16, dt: float = 0.25):
    """Flat metric grid; returns ``(g, spacing)``."""
    T, X, Y, Z, spacing = _grid(half_width=half_width, n_space=n_space, dt=dt)
    g = np.zeros(T.shape + (4, 4))
    g[..., 0, 0] = -1.0
    for i in range(1, 4):
        g[..., i, i] = 1.0
    return g, spacing


def _alcubierre_shape(rs: NDArray, radius: float, sigma: float) -> NDArray:
    norm = 2.0 * math.tanh(sigma * radius)
    return (np.tanh(sigma * (rs + radius)) - np.tanh(sigma * (rs - radius))) / norm


def alcubierre_metric(
    velocity: float = 0.5,
    bubble_radius: float = 2.0,
    sigma: float = 2.0,
    half_width: float = 4.0,
    n_space: int = 20,
    dt: float = 0.2,
    n_time: int = 5,
):
    """Alcubierre warp-bubble 4-metric grid; returns ``(g, spacing)``.

    ``ds^2 = -dt^2 + (dx - v f(r_s) dt)^2 + dy^2 + dz^2`` with the bubble centre
    moving as ``x_s(t) = v t`` (genuine time dependence so the Einstein tensor
    has non-trivial ``d_t`` terms).  ``n_time`` sets the stored time span
    (``±(n_time//2)·dt``); a stack shorter than the probe flight relies on the
    frozen-tail extrapolation of the evolving tracer.
    """
    T, X, Y, Z, spacing = _grid(
        half_width=half_width, n_space=n_space, dt=dt, n_time=n_time
    )
    xs = velocity * T
    rs = np.sqrt((X - xs) ** 2 + Y ** 2 + Z ** 2)
    f = _alcubierre_shape(rs, bubble_radius, sigma)
    vf = velocity * f

    g = np.zeros(T.shape + (4, 4))
    g[..., 0, 0] = -1.0 + vf ** 2
    g[..., 0, 1] = -vf
    g[..., 1, 0] = -vf
    g[..., 1, 1] = 1.0
    g[..., 2, 2] = 1.0
    g[..., 3, 3] = 1.0
    return g, spacing


# --------------------------------------------------------------------------
# Evolved GRTeclyn plotfiles -> effective stress-energy T^eff = G/8pi.
#
# A warp/portal whose exoticity is *geometric* (sourced by the shift/extrinsic
# curvature rather than a matter field) leaves the matter-sector energy
# conditions clean: the negative energy lives in T^eff_{mu nu}=G_{mu nu}/8 pi,
# which only the curvature of the evolved 4-metric reveals.  We reassemble the
# covariant 4-metric from a *stack* of consecutive plotfiles (so the Einstein
# tensor has genuine d_t terms) and run the same observer-robust + Hawking-Ellis
# energy-condition evaluator used for the analytic seeds.
# --------------------------------------------------------------------------

def build_four_metric_stack_from_plotfiles(
    plotfiles: Sequence[str],
    *,
    n_space: int = 32,
    half_width: float = 8.0,
) -> tuple[NDArray, tuple[float, float, float, float]]:
    """Assemble ``g_{mu nu}(t,x,y,z)`` from time-ordered GRTeclyn plotfiles.

    Each plotfile supplies one time slice of the ADM fields
    (``gamma_ij = h_ij/chi``, lapse ``alpha``, shift ``beta^i``) sampled on a
    common cubic box of ``n_space`` cells per axis centred on the domain.  The
    covariant 4-metric is built per slice via
    ``g_00 = beta_i beta^i - alpha^2``, ``g_0i = beta_i``, ``g_ij = gamma_ij``
    and stacked along the time axis.  Returns ``(g, (dt,dx,dy,dz))`` with
    ``g`` of shape ``(Nt, n_space, n_space, n_space, 4, 4)``.

    Requires at least three plotfiles (for the central time derivative).
    """
    import yt  # local import: heavy optional dependency

    paths = [str(p) for p in plotfiles]
    if len(paths) < 3:
        raise ValueError("need >= 3 consecutive plotfiles for d_t of the 4-metric")

    slices: list[NDArray] = []
    times: list[float] = []
    spacing_xyz: tuple[float, float, float] | None = None

    for path in paths:
        ds = yt.load(path)
        center = np.array([float(c.to_value()) for c in ds.domain_center])
        dleft = np.array([float(c.to_value()) for c in ds.domain_left_edge])
        dright = np.array([float(c.to_value()) for c in ds.domain_right_edge])
        left = np.maximum(center - half_width, dleft)
        right = np.minimum(center + half_width, dright)
        dims = (n_space, n_space, n_space)
        ag = ds.arbitrary_grid(left, right, dims)

        def field(name: str) -> NDArray:
            try:
                return np.asarray(ag["boxlib", name], dtype=float)
            except Exception:  # noqa: BLE001 - field-name fallback
                return np.asarray(ag[name], dtype=float)

        chi = np.clip(field("chi"), 1.0e-10, None)
        inv_chi = 1.0 / chi
        alpha = np.clip(field("lapse"), 1.0e-10, None)

        beta_up = np.stack(
            [field("shift1"), field("shift2"), field("shift3")], axis=-1
        )  # (n,n,n,3)
        gamma = np.zeros(chi.shape + (3, 3))
        gamma[..., 0, 0] = field("h11") * inv_chi
        gamma[..., 0, 1] = gamma[..., 1, 0] = field("h12") * inv_chi
        gamma[..., 0, 2] = gamma[..., 2, 0] = field("h13") * inv_chi
        gamma[..., 1, 1] = field("h22") * inv_chi
        gamma[..., 1, 2] = gamma[..., 2, 1] = field("h23") * inv_chi
        gamma[..., 2, 2] = field("h33") * inv_chi

        beta_low = np.einsum("...ij,...j->...i", gamma, beta_up)
        beta_sq = np.einsum("...i,...i->...", beta_low, beta_up)

        g = np.zeros(chi.shape + (4, 4))
        g[..., 0, 0] = beta_sq - alpha * alpha
        g[..., 0, 1:] = beta_low
        g[..., 1:, 0] = beta_low
        g[..., 1:, 1:] = gamma
        slices.append(g)
        times.append(float(ds.current_time))
        spacing_xyz = tuple((right - left) / (np.array(dims) - 1))

    g_stack = np.stack(slices, axis=0)  # (Nt, n, n, n, 4, 4)
    times_arr = np.asarray(times, dtype=float)
    dt = float(np.mean(np.diff(times_arr))) if len(times_arr) > 1 else 1.0
    if not math.isfinite(dt) or dt <= 0.0:
        dt = 1.0
    assert spacing_xyz is not None
    return g_stack, (dt, spacing_xyz[0], spacing_xyz[1], spacing_xyz[2])


def effective_energy_conditions_from_plotfiles(
    plotfiles: Sequence[str],
    *,
    n_space: int = 32,
    half_width: float = 8.0,
    crop: int = 4,
    **ec_kwargs,
) -> EnergyConditionReport:
    """Effective energy conditions ``T^eff=G/8pi`` of an evolved spacetime.

    Convenience wrapper: assemble the 4-metric stack from plotfiles and run the
    full Warp-Factory-style evaluator on the central time slice.  A negative
    ``nec_min`` / ``nec_slack_min`` certifies geometry-sourced exotic energy
    that the matter-sector conditions cannot see.
    """
    g, spacing = build_four_metric_stack_from_plotfiles(
        plotfiles, n_space=n_space, half_width=half_width
    )
    return evaluate_four_metric(g, spacing, crop=crop, **ec_kwargs)


# --------------------------------------------------------------------------
# Convergence self-test (cf. Lousto, arXiv:gr-qc/0503001).
# --------------------------------------------------------------------------

def _integrated_curvature(g: NDArray, spacing: Sequence[float], crop: int = 2) -> float:
    """Volume-integrated ``|G_{mu nu}|`` over the interior central time slice.

    A smooth scalar functional whose discretization error scales with the grid
    spacing; used to measure the observed order of convergence.
    """
    G = einstein_tensor(g, spacing)
    nt = g.shape[0]
    tc = nt // 2
    sl = (tc,) + tuple(slice(crop, dim - crop) for dim in g.shape[1:4])
    dV = spacing[1] * spacing[2] * spacing[3]
    return float(np.sum(np.abs(G[sl])) * dV)


def convergence_order(
    velocity: float = 0.5,
    resolutions: Sequence[int] = (16, 32, 64),
    *,
    half_width: float = 4.0,
    dt: float = 0.2,
) -> dict[str, object]:
    """Observed order of convergence of the finite-difference curvature chain.

    Builds the Alcubierre metric at successively doubled spatial resolutions,
    evaluates the smooth functional :func:`_integrated_curvature`, and estimates
    the convergence order by Richardson extrapolation::

        p = log2( |Q_c - Q_m| / |Q_m - Q_f| )

    for a coarse/medium/fine triple at fixed refinement ratio 2.  Reproduces the
    convergence-verification methodology of Lousto (arXiv:gr-qc/0503001).
    """
    res = list(resolutions)
    if len(res) < 3:
        raise ValueError("need at least three resolutions for a Richardson estimate")
    values: list[float] = []
    for n in res:
        g, spacing = alcubierre_metric(
            velocity=velocity, half_width=half_width, n_space=n, dt=dt
        )
        values.append(_integrated_curvature(g, spacing))

    orders: list[float] = []
    for i in range(len(values) - 2):
        d1 = abs(values[i] - values[i + 1])
        d2 = abs(values[i + 1] - values[i + 2])
        ratio = res[i + 2] / res[i + 1]  # typically 2
        if d2 > 0 and d1 > 0:
            orders.append(float(math.log(d1 / d2) / math.log(ratio)))
        else:
            orders.append(float("nan"))

    return {
        "resolutions": res,
        "functional": values,
        "observed_order": orders,
        "order_estimate": orders[-1] if orders else float("nan"),
    }
