"""Pre-flight constraint filter: evaluate the Hamiltonian constraint on a
coarse 1D radial grid before launching GRTeclyn, rejecting candidates
with massive violations to save GPU time.

The filter computes the vacuum Hamiltonian constraint residual
    H = R[chi] + (2/3) K^2 - A_ij A^ij  (with A_ij = 0 for the recipe)
and optionally accounts for the scalar field energy density.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

from .constrained_recipe import RecipeBasis, compute_ricci_scalar


@dataclass(frozen=True)
class PreflightResult:
    """Outcome of the pre-flight constraint check."""

    passed: bool
    hamiltonian_l2: float
    hamiltonian_max: float
    momentum_l2: float
    negative_chi_fraction: float
    reason: str


def _gradient_norm_sq(
    r: NDArray[np.float64],
    f: NDArray[np.float64],
    chi: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Compute |grad f|^2 = chi * (df/dr)^2 for conformally flat metric."""
    dr = r[1] - r[0]
    df = np.gradient(f, dr, edge_order=2)
    return chi * df * df


def preflight_check(
    overrides: Mapping[str, Any],
    *,
    ham_l2_threshold: float = 10.0,
    ham_max_threshold: float = 100.0,
    mom_l2_threshold: float = 10.0,
    negative_chi_threshold: float = 0.1,
    n_points: int = 512,
    phantom: bool = False,
) -> PreflightResult:
    """Evaluate the constraint residual for a proposed set of recipe
    coefficients without running GRTeclyn.

    Parameters
    ----------
    overrides : dict
        The proposed params overrides (must include recipe_* keys).
    ham_l2_threshold : float
        Maximum allowed L2 norm of the Hamiltonian constraint.
    ham_max_threshold : float
        Maximum allowed pointwise Hamiltonian constraint.
    mom_l2_threshold : float
        Maximum allowed L2 norm of the momentum constraint.
    negative_chi_threshold : float
        Maximum allowed fraction of grid points with chi < 0.
    n_points : int
        Number of radial grid points for the check.
    phantom : bool
        Whether the scalar field is phantom-coupled.

    Returns
    -------
    PreflightResult
    """
    num_bases = int(overrides.get("recipe_num_bases", 4))
    basis_width = float(overrides.get("recipe_basis_width", 1.0))
    basis_radius_max = float(overrides.get("recipe_basis_radius_max", 8.0))
    chi_asymptotic = float(overrides.get("recipe_chi_asymptotic", 1.0))
    K_asymptotic = float(overrides.get("recipe_K_asymptotic", 0.0))

    basis = RecipeBasis(
        num_bases=num_bases,
        basis_width=basis_width,
        basis_radius_max=basis_radius_max,
    )

    chi_coeffs = [
        float(overrides.get(f"recipe_chi_coeff_{n}", 0.0))
        for n in range(num_bases)
    ]
    K_coeffs = [
        float(overrides.get(f"recipe_K_coeff_{n}", 0.0))
        for n in range(num_bases)
    ]
    phi_coeffs = [
        float(overrides.get(f"recipe_phi_coeff_{n}", 0.0))
        for n in range(num_bases)
    ]
    Pi_coeffs = [
        float(overrides.get(f"recipe_Pi_coeff_{n}", 0.0))
        for n in range(num_bases)
    ]

    r_max = 2.0 * basis_radius_max
    r_min = r_max / n_points
    r = np.linspace(r_min, r_max, n_points)

    chi = basis.evaluate(r, chi_asymptotic, chi_coeffs)
    neg_chi_frac = float(np.mean(chi < 1.0e-10))
    chi = np.clip(chi, 1.0e-10, None)

    K = basis.evaluate(r, K_asymptotic, K_coeffs)
    phi = basis.evaluate(r, float(overrides.get("recipe_phi_asymptotic", 0.0)), phi_coeffs)
    Pi = basis.evaluate(r, float(overrides.get("recipe_Pi_asymptotic", 0.0)), Pi_coeffs)

    ricci = compute_ricci_scalar(r, chi)

    ham_vacuum = ricci + (2.0 / 3.0) * K * K

    sign = -1.0 if phantom else 1.0
    grad_phi_sq = _gradient_norm_sq(r, phi, chi)
    rho_matter = 0.5 * (Pi * Pi + grad_phi_sq)
    ham_residual = ham_vacuum - sign * 16.0 * math.pi * rho_matter

    dr = r[1] - r[0]
    ham_l2 = float(np.sqrt(np.mean(ham_residual**2)))
    ham_max = float(np.max(np.abs(ham_residual)))

    dK = np.gradient(K, dr, edge_order=2)
    mom_vacuum = -(2.0 / 3.0) * dK
    dPhi = np.gradient(phi, dr, edge_order=2)
    j_matter = -Pi * dPhi
    mom_residual = mom_vacuum - sign * 8.0 * math.pi * j_matter
    mom_l2 = float(np.sqrt(np.mean(mom_residual**2)))

    reasons = []
    if ham_l2 > ham_l2_threshold:
        reasons.append(f"ham_l2={ham_l2:.3g}>{ham_l2_threshold}")
    if ham_max > ham_max_threshold:
        reasons.append(f"ham_max={ham_max:.3g}>{ham_max_threshold}")
    if mom_l2 > mom_l2_threshold:
        reasons.append(f"mom_l2={mom_l2:.3g}>{mom_l2_threshold}")
    if neg_chi_frac > negative_chi_threshold:
        reasons.append(f"neg_chi={neg_chi_frac:.2%}>{negative_chi_threshold:.0%}")

    passed = len(reasons) == 0
    reason = "; ".join(reasons) if reasons else "ok"

    return PreflightResult(
        passed=passed,
        hamiltonian_l2=ham_l2,
        hamiltonian_max=ham_max,
        momentum_l2=mom_l2,
        negative_chi_fraction=neg_chi_frac,
        reason=reason,
    )
