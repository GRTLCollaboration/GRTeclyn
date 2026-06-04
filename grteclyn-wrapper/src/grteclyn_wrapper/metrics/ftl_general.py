"""General, mechanism-agnostic operational FTL diagnostics.

The existing :mod:`ftl_metrics` module measures an FTL shortcut along the
``x`` axis only and assumes the conformally-flat warp form
(``c = alpha*sqrt(chi) - beta_x``).  That is specific to axial warp drives.

This module instead asks the *general* question: given the full 3+1 data
``(alpha, beta^i, gamma_ij)`` of any candidate spacetime, can a causal signal
travel between two asymptotic endpoints faster than the flat-space baseline,
*regardless of the mechanism* (warp push, wormhole throat, portal, or some
geometry we have not named) or its orientation?

The two ingredients are:

1. ``coordinate_light_speed`` -- the exact one-way coordinate speed of light in
   an arbitrary direction, obtained from the null condition
   ``gamma_ij (u^i + beta^i)(u^j + beta^j) = alpha^2``.  No warp/conformal
   assumption: it takes the full (possibly off-diagonal) spatial metric and the
   full shift vector.

2. ``operational_ftl_on_grid`` -- a Dijkstra shortest-*coordinate-time* search
   over a field grid.  Edge cost is the coordinate time to cross the edge along
   its direction; the globally fastest causal path is compared to the flat
   baseline between the same endpoints.  Because the search is free to take any
   in-plane direction, it finds shortcuts of any mechanism/orientation, not
   just an assumed axial one.
"""

from __future__ import annotations

import heapq
import math
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

from ..initial_data.constrained_recipe import RecipeBasis


@dataclass(frozen=True)
class GeneralFtlReport:
    """Result of the general operational FTL search on a 2D field slice."""

    f_op: float  #: (t_flat - t_min) / t_flat, clipped to >= 0; > 0 => shortcut
    t_min: float | None  #: fastest causal coordinate-time from source to target
    t_flat: float  #: flat-space baseline coordinate-time (Euclidean distance)
    max_local_speed: float  #: max one-way coordinate light speed on the slice
    superluminal_fraction: float  #: fraction of cells with local speed > 1
    path_offaxis: bool  #: whether the fastest path leaves the straight axis
    reachable: bool  #: whether the target was reachable from the source
    notes: tuple[str, ...]
    max_shift: float = 0.0  #: max shift-vector magnitude on the sampled slice


def coordinate_light_speed(
    alpha: float,
    beta: Sequence[float],
    gamma: NDArray[np.float64],
    direction: Sequence[float],
) -> float:
    """Exact one-way coordinate speed of light along ``direction``.

    Solves the null condition ``gamma_ij (s n^i + beta^i)(s n^j + beta^j) =
    alpha^2`` for the forward coordinate speed ``s = |dx/dt|`` along the unit
    Euclidean direction ``n``.  Returns 0.0 if the direction is causally
    blocked (no positive real root), and ``inf`` is never returned.

    This is fully general: ``gamma`` may be any symmetric positive-definite
    3x3 (or 2x2) matrix and ``beta`` any shift vector, so a value > 1 flags a
    genuinely superluminal coordinate channel of any mechanism.
    """
    n = np.asarray(direction, dtype=float)
    norm = float(np.linalg.norm(n))
    if norm <= 0.0:
        return 0.0
    n = n / norm
    g = np.asarray(gamma, dtype=float)
    b = np.asarray(beta, dtype=float)

    a_coeff = float(n @ g @ n)
    b_coeff = float(n @ g @ b)
    c_coeff = float(b @ g @ b) - float(alpha) * float(alpha)
    if a_coeff <= 0.0:
        return 0.0
    disc = b_coeff * b_coeff - a_coeff * c_coeff
    if disc < 0.0:
        return 0.0
    s = (-b_coeff + math.sqrt(disc)) / a_coeff
    return s if s > 0.0 else 0.0


def _local_max_speed_grid(
    alpha: NDArray[np.float64],
    beta: NDArray[np.float64],
    gamma: NDArray[np.float64],
    spacing: tuple[float, float],
) -> NDArray[np.float64]:
    """Vectorised per-cell maximum forward coordinate light speed over a small
    set of probe directions (used only for the superluminal-fraction
    diagnostic; the path search uses the exact per-edge speed)."""
    di_phys, dj_phys = spacing
    probes = [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, -1.0)]
    a2 = alpha * alpha
    best = np.zeros_like(alpha)
    for d in probes:
        for sgn in (1.0, -1.0):
            vec = np.array([sgn * d[0] * di_phys, sgn * d[1] * dj_phys])
            norm = float(np.linalg.norm(vec))
            if norm <= 0.0:
                continue
            n = vec / norm
            # A = n.gamma.n, B = n.gamma.beta, C = beta.gamma.beta - alpha^2
            gn0 = gamma[:, :, 0, 0] * n[0] + gamma[:, :, 0, 1] * n[1]
            gn1 = gamma[:, :, 1, 0] * n[0] + gamma[:, :, 1, 1] * n[1]
            a_c = n[0] * gn0 + n[1] * gn1
            b_c = beta[:, :, 0] * gn0 + beta[:, :, 1] * gn1
            gb0 = gamma[:, :, 0, 0] * beta[:, :, 0] + gamma[:, :, 0, 1] * beta[:, :, 1]
            gb1 = gamma[:, :, 1, 0] * beta[:, :, 0] + gamma[:, :, 1, 1] * beta[:, :, 1]
            c_c = beta[:, :, 0] * gb0 + beta[:, :, 1] * gb1 - a2
            disc = b_c * b_c - a_c * c_c
            with np.errstate(invalid="ignore", divide="ignore"):
                s = (-b_c + np.sqrt(np.maximum(disc, 0.0))) / np.where(a_c > 0.0, a_c, np.nan)
            s = np.where(np.isfinite(s) & (s > 0.0), s, 0.0)
            best = np.maximum(best, s)
    return best


def _neighbor_offsets() -> list[tuple[int, int]]:
    return [
        (di, dj)
        for di in (-1, 0, 1)
        for dj in (-1, 0, 1)
        if not (di == 0 and dj == 0)
    ]


def operational_ftl_on_grid(
    alpha: NDArray[np.float64],
    beta: NDArray[np.float64],
    gamma: NDArray[np.float64],
    *,
    spacing: tuple[float, float],
    source: tuple[int, int],
    target: tuple[int, int],
) -> GeneralFtlReport:
    """Dijkstra shortest-coordinate-time search over a 2D field slice.

    Parameters
    ----------
    alpha : (Ni, Nj) lapse on the slice.
    beta : (Ni, Nj, 2) shift components in the slice plane.
    gamma : (Ni, Nj, 2, 2) spatial metric restricted to the slice plane.
    spacing : (dx_i, dx_j) coordinate spacing of the two slice axes.
    source, target : (i, j) grid indices of the asymptotic endpoints.

    The edge cost from cell ``p`` to neighbour ``q`` is the coordinate time
    ``|pq| / s`` where ``s`` is the local coordinate light speed along the edge
    direction (metric averaged over the edge).  The flat baseline is the
    Euclidean coordinate distance between the endpoints (asymptotic c = 1).
    """
    ni, nj = alpha.shape
    di_phys, dj_phys = spacing
    offsets = _neighbor_offsets()
    notes: list[str] = []

    # Local maximum coordinate speed (over a few probe directions) per cell,
    # for the superluminal-fraction diagnostic.  Vectorised over the grid.
    local_max = _local_max_speed_grid(alpha, beta, gamma, spacing)
    max_local_speed = float(local_max.max())
    superluminal_fraction = float(np.mean(local_max > 1.0))
    max_shift = float(np.linalg.norm(beta, axis=-1).max())

    dist = np.full((ni, nj), np.inf, dtype=float)
    prev = np.full((ni, nj, 2), -1, dtype=int)
    si, sj = source
    dist[si, sj] = 0.0
    heap: list[tuple[float, int, int]] = [(0.0, si, sj)]

    ti, tj = target
    while heap:
        d, i, j = heapq.heappop(heap)
        if d > dist[i, j]:
            continue
        if (i, j) == (ti, tj):
            break
        for di, dj in offsets:
            qi, qj = i + di, j + dj
            if not (0 <= qi < ni and 0 <= qj < nj):
                continue
            edge = np.array([di * di_phys, dj * dj_phys])
            length = float(np.linalg.norm(edge))
            # Average the metric over the edge endpoints.
            a_avg = 0.5 * (float(alpha[i, j]) + float(alpha[qi, qj]))
            b_avg = 0.5 * (beta[i, j] + beta[qi, qj])
            g_avg = 0.5 * (gamma[i, j] + gamma[qi, qj])
            speed = coordinate_light_speed(a_avg, b_avg, g_avg, edge)
            if speed <= 0.0:
                continue
            cost = length / speed
            nd = d + cost
            if nd < dist[qi, qj]:
                dist[qi, qj] = nd
                prev[qi, qj] = (i, j)
                heapq.heappush(heap, (nd, qi, qj))

    reachable = math.isfinite(dist[ti, tj])
    t_min = float(dist[ti, tj]) if reachable else None

    src = np.array([si * di_phys, sj * dj_phys])
    tgt = np.array([ti * di_phys, tj * dj_phys])
    t_flat = float(np.linalg.norm(tgt - src))

    f_op = 0.0
    if reachable and t_flat > 0.0 and t_min is not None:
        f_op = max(0.0, (t_flat - t_min) / t_flat)

    # Reconstruct the path to check whether the fastest route went off-axis.
    path_offaxis = False
    if reachable:
        axis = tgt - src
        axis_len = float(np.linalg.norm(axis))
        if axis_len > 0.0:
            axis_hat = axis / axis_len
            ci, cj = ti, tj
            max_perp = 0.0
            while (ci, cj) != (si, sj) and prev[ci, cj][0] >= 0:
                p = np.array([ci * di_phys, cj * dj_phys]) - src
                perp = p - float(p @ axis_hat) * axis_hat
                max_perp = max(max_perp, float(np.linalg.norm(perp)))
                ci, cj = int(prev[ci, cj][0]), int(prev[ci, cj][1])
            path_offaxis = max_perp > 2.0 * max(di_phys, dj_phys)
    else:
        notes.append("target unreachable: no forward causal path on slice")

    if max_local_speed > 1.0:
        notes.append(
            f"superluminal coordinate channel present (max c = {max_local_speed:.3f})"
        )

    return GeneralFtlReport(
        f_op=f_op,
        t_min=t_min,
        t_flat=t_flat,
        max_local_speed=max_local_speed,
        superluminal_fraction=superluminal_fraction,
        path_offaxis=path_offaxis,
        reachable=reachable,
        notes=tuple(notes),
        max_shift=max_shift,
    )


def _legendre_p(ell: int, mu: NDArray[np.float64]) -> NDArray[np.float64]:
    if ell == 0:
        return np.ones_like(mu)
    if ell == 1:
        return mu
    if ell == 2:
        return 0.5 * (3.0 * mu * mu - 1.0)
    if ell == 3:
        return 0.5 * (5.0 * mu**3 - 3.0 * mu)
    return np.zeros_like(mu)


def _angular_sum(
    overrides: Mapping[str, object],
    prefix: str,
    r: NDArray[np.float64],
    mu: NDArray[np.float64],
) -> NDArray[np.float64]:
    num = int(overrides.get(f"recipe_num_{prefix}_angular_modes", 0))
    delta = np.zeros_like(r)
    for n in range(num):
        ell = int(overrides.get(f"recipe_{prefix}_mode_ell_{n}", 0))
        amp = float(overrides.get(f"recipe_{prefix}_mode_amp_{n}", 0.0))
        rc = float(overrides.get(f"recipe_{prefix}_mode_rc_{n}", 0.0))
        rw = float(overrides.get(f"recipe_{prefix}_mode_rw_{n}", 1.0))
        radial = np.exp(-((r - rc) ** 2) / (2.0 * rw * rw))
        delta = delta + amp * radial * _legendre_p(ell, mu)
    return delta


def build_tslice_fields_xz(
    overrides: Mapping[str, object],
    *,
    L: float,
    n: int = 161,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], tuple[float, float]]:
    """Reconstruct t=0 ``(alpha, beta, gamma)`` on the x-z plane from the recipe.

    The shift axis (x) and the Legendre angular axis (z) both lie in this
    plane, so it is the natural slice on which to look for an axisymmetric
    shortcut.  ``gamma`` is conformally flat at t=0 (``gamma_ij = delta_ij /
    chi``); the general machinery is metric-agnostic and will use the full
    off-diagonal metric unchanged once it is fed from evolved plotfiles.
    """
    num_bases = int(overrides.get("recipe_num_bases", 4))
    basis = RecipeBasis(
        num_bases=num_bases,
        basis_width=float(overrides.get("recipe_basis_width", 1.0)),
        basis_radius_max=float(overrides.get("recipe_basis_radius_max", 8.0)),
    )

    def coeffs(prefix: str) -> list[float]:
        return [float(overrides.get(f"{prefix}_{k}", 0.0)) for k in range(num_bases)]

    axis = np.linspace(-L, L, n)
    xx, zz = np.meshgrid(axis, axis, indexing="ij")
    r = np.sqrt(xx * xx + zz * zz)
    mu = np.divide(zz, np.maximum(r, 1.0e-12))

    chi = basis.evaluate(r, float(overrides.get("recipe_chi_asymptotic", 1.0)), coeffs("recipe_chi_coeff"))
    chi = chi + _angular_sum(overrides, "chi", r, mu)
    chi = np.clip(chi, 1.0e-10, None)

    alpha = basis.evaluate(r, float(overrides.get("recipe_alpha_asymptotic", 1.0)), coeffs("recipe_alpha_coeff"))
    alpha = alpha + _angular_sum(overrides, "lapse", r, mu)
    alpha = np.clip(alpha, 1.0e-10, None)

    beta_x = basis.evaluate(r, float(overrides.get("recipe_beta_asymptotic", 0.0)), coeffs("recipe_beta_coeff"))
    beta_x = beta_x + _angular_sum(overrides, "beta", r, mu)

    beta = np.zeros((n, n, 2), dtype=float)
    beta[:, :, 0] = beta_x  # shift is along x

    gamma = np.zeros((n, n, 2, 2), dtype=float)
    inv_chi = 1.0 / chi
    gamma[:, :, 0, 0] = inv_chi
    gamma[:, :, 1, 1] = inv_chi

    spacing = (axis[1] - axis[0], axis[1] - axis[0])
    return alpha, beta, gamma, spacing


def compute_general_ftl(
    overrides: Mapping[str, object],
    *,
    L: float | None = None,
    n: int = 161,
) -> GeneralFtlReport:
    """Convenience: build the t=0 x-z slice and run the general FTL search
    between the two asymptotic endpoints on the x axis (z = 0)."""
    integration_L = L if L is not None else float(overrides.get("recipe_basis_radius_max", 8.0))
    alpha, beta, gamma, spacing = build_tslice_fields_xz(overrides, L=integration_L, n=n)
    mid = n // 2
    source = (1, mid)
    target = (n - 2, mid)
    return operational_ftl_on_grid(
        alpha, beta, gamma, spacing=spacing, source=source, target=target
    )


# --------------------------------------------------------------------------- #
# Evolved-spacetime FTL: run the same operational search on a plotfile slice.
# Requires the full metric in amr.plot_vars (chi h11..h33 lapse shift1..3).
# --------------------------------------------------------------------------- #

_PLOTFILE_RE = re.compile(r"(\d+)$")


def plotfile_is_complete(plotfile: str | Path, *, stable_seconds: float = 2.0) -> bool:
    """Return ``True`` when an AMReX plotfile directory is fully written.

    A plotfile that is still being flushed by the solver (the scoring step can
    race the final write) has either no top-level ``Header`` yet or files whose
    mtime is still advancing.  We require both: a present ``Header`` and no file
    under the directory modified within the last ``stable_seconds``.
    """
    p = Path(plotfile)
    if not (p / "Header").is_file():
        return False
    try:
        newest = max(
            (f.stat().st_mtime for f in p.rglob("*") if f.is_file()),
            default=0.0,
        )
    except OSError:
        return False
    return (time.time() - newest) >= stable_seconds


def wait_for_plotfile_complete(
    plotfile: str | Path,
    *,
    timeout: float = 30.0,
    stable_seconds: float = 2.0,
    poll: float = 1.0,
) -> bool:
    """Block until ``plotfile`` looks complete, or ``timeout`` elapses."""
    deadline = time.time() + max(0.0, timeout)
    while True:
        if plotfile_is_complete(plotfile, stable_seconds=stable_seconds):
            return True
        if time.time() >= deadline:
            return plotfile_is_complete(plotfile, stable_seconds=stable_seconds)
        time.sleep(poll)


def _scan_plotfiles(episode_dir: Path, *, complete_only: bool) -> dict[int, Path]:
    found: dict[int, Path] = {}
    for pattern in ("*Plt*", "plt*"):
        for p in episode_dir.rglob(pattern):
            if not p.is_dir():
                continue
            m = _PLOTFILE_RE.search(p.name)
            if not m:
                continue
            if complete_only and not plotfile_is_complete(p):
                continue
            found[int(m.group(1))] = p
    return found


def find_latest_plotfile(
    episode_dir: str | Path, *, complete_only: bool = True
) -> Path | None:
    """Return the highest-index AMReX plotfile directory under ``episode_dir``.

    Plotfiles are named ``<Prefix>Plt<NNNNN>`` (or ``plt<NNNNN>``); they may live
    directly in the episode or under an ``hdf5`` subdir.  Returns ``None`` when
    no plotfile is present (e.g. the streaming consumer deleted them).  By
    default only fully-written plotfiles are considered, so scoring never reads
    a half-flushed final plotfile (the evolved-FTL/effective-EC race)."""
    episode_dir = Path(episode_dir)
    if not episode_dir.exists():
        return None
    found = _scan_plotfiles(episode_dir, complete_only=complete_only)
    if not found:
        return None
    return found[max(found)]


def find_recent_plotfiles(
    episode_dir: str | Path, count: int = 5, *, complete_only: bool = True
) -> list[Path]:
    """Return up to ``count`` highest-index plotfiles, time-ordered ascending.

    Used to assemble a short time stack for finite-difference ``d_t`` of the
    evolved 4-metric (effective energy conditions).  Returns ``[]`` when fewer
    than two plotfiles survive.  Skips half-written plotfiles by default."""
    episode_dir = Path(episode_dir)
    if not episode_dir.exists():
        return []
    found = _scan_plotfiles(episode_dir, complete_only=complete_only)
    if not found:
        return []
    idx_sorted = sorted(found)[-count:]
    return [found[i] for i in idx_sorted]


def build_slice_fields_xz_from_plotfile(
    plotfile: str | Path,
    *,
    n: int = 129,
    L: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64],
           tuple[float, float], float]:
    """Sample ``(alpha, beta, gamma)`` on the x-z plane of an evolved plotfile.

    Uses the full evolved 3+1 fields (``gamma_ij = h_ij / chi``, lapse, shift),
    so the operational FTL search runs on the genuine evolved geometry rather
    than the t=0 reconstruction.  Returns the fields plus the slice spacing and
    the plotfile coordinate time.
    """
    import yt  # local import: heavy optional dependency

    ds = yt.load(str(plotfile))
    center = ds.domain_center
    # Per-axis half-extents: stay just inside the domain so the x endpoints sit
    # in the asymptotic region.  The x-z plane can be rectangular (e.g. a
    # half-z box), so width (x) and height (z) are sized independently.
    half_x = 0.49 * float(ds.domain_width[0].to_value())
    half_z = 0.49 * float(ds.domain_width[2].to_value())
    if L is not None:
        half_x = min(half_x, float(L))
        half_z = min(half_z, float(L))
    cy = float(center[1].to_value())
    slc = ds.slice(1, cy)  # axis 1 = y; image plane spans (x, z)
    frb = slc.to_frb(
        (2.0 * half_x, "code_length"),
        n,
        center=center,
        height=(2.0 * half_z, "code_length"),
    )

    def field(name: str) -> NDArray[np.float64]:
        try:
            arr = np.asarray(frb["boxlib", name], dtype=float)
        except Exception:  # noqa: BLE001 - field-name fallback
            arr = np.asarray(frb[name], dtype=float)
        # FRB is indexed [image_y=z, image_x=x]; transpose to [x, z].
        return arr.T

    chi = np.clip(field("chi"), 1.0e-10, None)
    alpha = np.clip(field("lapse"), 1.0e-10, None)
    inv_chi = 1.0 / chi

    beta = np.zeros((n, n, 2), dtype=float)
    beta[:, :, 0] = field("shift1")
    beta[:, :, 1] = field("shift3")

    gamma = np.zeros((n, n, 2, 2), dtype=float)
    gamma[:, :, 0, 0] = field("h11") * inv_chi
    gamma[:, :, 0, 1] = gamma[:, :, 1, 0] = field("h13") * inv_chi
    gamma[:, :, 1, 1] = field("h33") * inv_chi

    dx = 2.0 * half_x / (n - 1)
    dz = 2.0 * half_z / (n - 1)
    return alpha, beta, gamma, (dx, dz), float(ds.current_time)


def compute_general_ftl_from_plotfile(
    plotfile: str | Path,
    *,
    n: int = 129,
    L: float | None = None,
) -> GeneralFtlReport:
    """Run the operational FTL search on the evolved geometry of a plotfile."""
    alpha, beta, gamma, spacing, _t = build_slice_fields_xz_from_plotfile(
        plotfile, n=n, L=L
    )
    mid = n // 2
    return operational_ftl_on_grid(
        alpha, beta, gamma, spacing=spacing, source=(1, mid), target=(n - 2, mid)
    )
