"""Gauge-robust physical metrics computed from the t=0 recipe profiles.

These implement the proposed extensions of the paper (Sec. "Proposed
Extensions"):

* An *ANEC line proxy*: the energy density integrated along the FTL travel
  axis, a directional energy-condition indicator complementary to the
  Eulerian volume integral already used by ``score.py``.
* A *coordinate-invariant tidal/curvature proxy*: the spatial Ricci-scalar
  magnitude plus lapse acceleration-gradient along the axis, a proxy for the
  rider tidal force / structural-integrity bound.
* A *trapped-surface flag*: an interior local minimum of the areal radius,
  flagging a potential throat/horizon precursor at t=0.

All quantities are t=0 proxies evaluated from the smooth recipe basis; they
do not replace the full geodesic ANEC / electric-Weyl tetrad contraction,
which require the C++ curvature diagnostics (Phase 2).  They are clearly
labelled as proxies and feed the score as bounded rewards.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray

from ...initial_data.constrained_recipe import RecipeBasis, compute_ricci_scalar


@dataclass(frozen=True)
class PhysicalMetrics:
    """t=0 gauge-robust proxy metrics."""

    anec_line: float | None
    anec_min: float | None
    tidal_proxy: float | None
    ricci_max: float | None
    lapse_curvature_max: float | None
    has_trapped_proxy: bool
    s_anec: float | None
    s_tidal: float | None
    qei_spatial_proxy: float | None = None
    notes: tuple[str, ...] = ()


def _bounded_reward(value: float, scale: float) -> float:
    if not math.isfinite(value) or value < 0:
        return 0.0
    return 1.0 / (1.0 + value / scale)


def _coeff_list(overrides: Mapping[str, Any], prefix: str, num_bases: int) -> list[float]:
    return [float(overrides.get(f"{prefix}_{n}", 0.0)) for n in range(num_bases)]


def compute_physical_metrics(
    overrides: Mapping[str, Any],
    *,
    L: float | None = None,
    n_points: int = 1024,
    anec_scale: float = 0.1,
    tidal_scale: float = 1.0,
) -> PhysicalMetrics:
    """Compute t=0 ANEC-line, tidal/curvature, and trapped-surface proxies.

    Parameters
    ----------
    overrides
        Recipe params (must contain ``recipe_chi_coeff_*`` etc.).
    L
        Half-axis length for the line integrals.  Defaults to
        ``recipe_basis_radius_max``.
    anec_scale, tidal_scale
        Bounded-reward scales controlling how sharply each reward decays.
    """
    num_bases = int(overrides.get("recipe_num_bases", 4))
    basis_radius_max = float(overrides.get("recipe_basis_radius_max", 8.0))
    integration_L = float(L) if L is not None else basis_radius_max
    basis = RecipeBasis(
        num_bases=num_bases,
        basis_width=float(overrides.get("recipe_basis_width", 1.0)),
        basis_radius_max=basis_radius_max,
    )

    chi_asym = float(overrides.get("recipe_chi_asymptotic", 1.0))
    alpha_asym = float(overrides.get("recipe_alpha_asymptotic", 1.0))
    k_const = float(overrides.get("recipe_K_asymptotic", 0.0))
    chi_coeffs = _coeff_list(overrides, "recipe_chi_coeff", num_bases)
    alpha_coeffs = _coeff_list(overrides, "recipe_alpha_coeff", num_bases)

    notes: list[str] = []

    # Radial grid (avoid r=0 where the spherical Laplacian term diverges).
    r_min = max(integration_L / n_points, 1.0e-6)
    r = np.linspace(r_min, integration_L, n_points)

    chi = np.clip(basis.evaluate(r, chi_asym, chi_coeffs), 1.0e-10, None)
    alpha = np.clip(basis.evaluate(r, alpha_asym, alpha_coeffs), 1.0e-10, None)

    ricci = compute_ricci_scalar(r, chi)
    ricci = np.where(np.isfinite(ricci), ricci, 0.0)

    # Required (vacuum-Hamiltonian) energy density: A_ij = 0 in the recipe.
    rho_req = (ricci + (2.0 / 3.0) * k_const**2) / (16.0 * math.pi)

    # ANEC line proxy: energy density integrated along the (symmetric) axis.
    anec_line = 2.0 * float(np.trapezoid(rho_req, r))
    anec_min = float(np.min(rho_req))

    # Coordinate-invariant-ish tidal proxy: spatial curvature scale plus the
    # lapse acceleration-gradient (geodesic-deviation precursor).
    ricci_max = float(np.max(np.abs(ricci)))
    dr = r[1] - r[0]
    dalpha = np.gradient(alpha, dr, edge_order=2)
    d2alpha = np.gradient(dalpha, dr, edge_order=2)
    lapse_curvature_max = float(np.max(np.abs(d2alpha)))
    tidal_proxy = ricci_max + lapse_curvature_max

    # Trapped-surface proxy: interior local minimum of the areal radius.
    r_areal = r / np.sqrt(chi)
    has_trapped = False
    if r_areal.size >= 5:
        interior = r_areal[1:-1]
        left = r_areal[:-2]
        right = r_areal[2:]
        local_min = (interior <= left) & (interior <= right)
        r_asymp = float(r_areal[-1])
        if r_asymp > 1.0e-8:
            has_trapped = bool(np.any(local_min & (interior < 0.95 * r_asymp)))

    # Rewards: ANEC penalizes net-negative line energy (exotic matter along the
    # travel path); tidal penalizes large curvature.  Higher is better.
    s_anec = _bounded_reward(max(0.0, -anec_line), anec_scale)
    s_tidal = _bounded_reward(tidal_proxy, tidal_scale)

    rho_minus = np.minimum(rho_req, 0.0)
    if np.any(rho_minus < 0.0):
        mask = rho_minus < 0.0
        extent = float(np.max(r[mask]) - np.min(r[mask])) if mask.any() else 1.0
        extent = max(extent, 1.0e-3)
        qei_spatial = float(np.trapezoid(np.abs(rho_minus) * extent ** -4, r))
    else:
        qei_spatial = 0.0

    if anec_line < 0.0:
        notes.append("ANEC line integral negative (exotic matter along axis)")
    if has_trapped:
        notes.append("trapped-surface proxy: interior areal-radius minimum at t=0")

    return PhysicalMetrics(
        anec_line=anec_line,
        anec_min=anec_min,
        tidal_proxy=tidal_proxy,
        ricci_max=ricci_max,
        lapse_curvature_max=lapse_curvature_max,
        has_trapped_proxy=has_trapped,
        s_anec=s_anec,
        s_tidal=s_tidal,
        qei_spatial_proxy=qei_spatial,
        notes=tuple(notes),
    )
