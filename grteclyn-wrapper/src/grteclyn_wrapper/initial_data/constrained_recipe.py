"""Constrained initial-data recipe: derive phi from chi so that the
Hamiltonian constraint is approximately satisfied at t=0.

Physics
-------
For conformally-flat data (h_ij = delta_ij) with A_ij = 0, K = const,
Pi = 0, the momentum constraint is identically satisfied and the
Hamiltonian constraint reduces to:

    R[chi] + (2/3) K^2  =  s * 8 pi |grad phi|^2

where s = +1 for a normal scalar field and s = -1 for a phantom field.

The 3-Ricci scalar for gamma_ij = chi^{-1} delta_ij (with psi = chi^{-1/4})
is  R = -8 psi^{-5} Laplacian(psi)  in flat coordinates.

In spherical symmetry this gives a 1D ODE for dphi/dr which we integrate
to obtain a constraint-consistent phi(r), then fit it back to the Gaussian
basis used by RadialRecipeInitialData.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from numpy.typing import NDArray


@dataclass
class RecipeBasis:
    """Mirrors RadialRecipeInitialData::params_t Gaussian basis layout."""

    num_bases: int = 4
    basis_width: float = 1.0
    basis_radius_max: float = 8.0

    def nodes(self) -> NDArray[np.float64]:
        if self.num_bases <= 1:
            return np.array([0.0])
        return np.linspace(0.0, self.basis_radius_max, self.num_bases)

    def evaluate(
        self,
        r: NDArray[np.float64],
        asymptotic: float,
        coeffs: Sequence[float],
    ) -> NDArray[np.float64]:
        nodes = self.nodes()
        w2 = self.basis_width**2
        result = np.full_like(r, asymptotic, dtype=np.float64)
        for n in range(self.num_bases):
            dr = r - nodes[n]
            result += coeffs[n] * np.exp(-dr * dr / (2.0 * w2))
        return result


@dataclass
class ConstrainedResult:
    """Output of the constraint solver."""

    phi_profile: NDArray[np.float64]
    phi_coeffs: list[float]
    phi_asymptotic: float
    r: NDArray[np.float64]
    ricci_scalar: NDArray[np.float64]
    rho_required: NDArray[np.float64]
    integrand_negative_mask: NDArray[np.bool_]
    fit_residual: float


def _radial_grid(r_max: float, n_points: int) -> NDArray[np.float64]:
    """Uniform radial grid from a small r_min to r_max."""
    r_min = r_max / n_points
    return np.linspace(r_min, r_max, n_points)


def compute_ricci_scalar(
    r: NDArray[np.float64],
    chi: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Compute the 3-Ricci scalar R for gamma_ij = chi^{-1} delta_ij.

    Uses  R = -8 psi^{-5} Laplacian(psi)  where psi = chi^{-1/4},
    and the spherical Laplacian  Laplacian(f) = f'' + (2/r) f'.
    """
    psi = np.power(chi, -0.25)
    dr = r[1] - r[0]

    dpsi = np.gradient(psi, dr, edge_order=2)
    d2psi = np.gradient(dpsi, dr, edge_order=2)

    laplacian_psi = d2psi + (2.0 / r) * dpsi
    ricci = -8.0 * np.power(psi, -5.0) * laplacian_psi
    return ricci


def solve_phi_from_constraint(
    r: NDArray[np.float64],
    chi: NDArray[np.float64],
    ricci: NDArray[np.float64],
    K_const: float = 0.0,
    phantom: bool = False,
) -> tuple[NDArray[np.float64], NDArray[np.bool_]]:
    """Integrate the Hamiltonian constraint for dphi/dr, then integrate to
    get phi(r).

    Returns (phi, negative_mask) where negative_mask marks regions where
    the required (dphi/dr)^2 < 0 (unsatisfiable by the chosen field type).
    """
    sign = -1.0 if phantom else 1.0
    source = sign * chi * (ricci + (2.0 / 3.0) * K_const**2) / (8.0 * math.pi)

    negative_mask = source < 0.0
    dphi_dr_sq = np.clip(source, 0.0, None)
    dphi_dr = np.sqrt(dphi_dr_sq)

    dr = r[1] - r[0]
    phi = np.cumsum(dphi_dr) * dr
    phi -= phi[-1]

    return phi, negative_mask


def fit_gaussian_basis(
    r: NDArray[np.float64],
    profile: NDArray[np.float64],
    basis: RecipeBasis,
    asymptotic: float = 0.0,
) -> tuple[list[float], float]:
    """Least-squares fit of a radial profile to the Gaussian basis.

    Returns (coefficients, residual_rms).
    """
    nodes = basis.nodes()
    w2 = basis.basis_width**2
    n = basis.num_bases

    A = np.zeros((len(r), n))
    for k in range(n):
        dr = r - nodes[k]
        A[:, k] = np.exp(-dr * dr / (2.0 * w2))

    target = profile - asymptotic
    coeffs, residuals, _, _ = np.linalg.lstsq(A, target, rcond=None)
    fit = asymptotic + A @ coeffs
    rms = float(np.sqrt(np.mean((fit - profile) ** 2)))
    return [float(c) for c in coeffs], rms


def constrained_phi(
    basis: RecipeBasis,
    chi_asymptotic: float,
    chi_coeffs: Sequence[float],
    K_const: float = 0.0,
    phantom: bool = False,
    r_max: float | None = None,
    n_points: int = 2048,
) -> ConstrainedResult:
    """Main entry point: given chi basis coefficients, derive constraint-
    consistent phi coefficients.

    Parameters
    ----------
    basis : RecipeBasis
        Gaussian basis layout (shared between chi and phi).
    chi_asymptotic : float
        Asymptotic value of chi (should be 1.0 for asymptotically flat).
    chi_coeffs : sequence of float
        Gaussian basis coefficients for chi.
    K_const : float
        Spatially constant trace of extrinsic curvature.
    phantom : bool
        If True, use phantom scalar field coupling (negative kinetic term).
    r_max : float or None
        Outer radius of the integration grid.  Defaults to
        2 * basis_radius_max.
    n_points : int
        Number of radial grid points.

    Returns
    -------
    ConstrainedResult
        Contains the derived phi profile, fitted coefficients, diagnostics.
    """
    if r_max is None:
        r_max = 2.0 * basis.basis_radius_max
    r = _radial_grid(r_max, n_points)

    chi = basis.evaluate(r, chi_asymptotic, chi_coeffs)
    chi = np.clip(chi, 1.0e-10, None)

    ricci = compute_ricci_scalar(r, chi)
    phi, neg_mask = solve_phi_from_constraint(r, chi, ricci, K_const, phantom)

    sign = -1.0 if phantom else 1.0
    rho_required = sign * (ricci + (2.0 / 3.0) * K_const**2) / (16.0 * math.pi)

    phi_coeffs, residual = fit_gaussian_basis(r, phi, basis, asymptotic=0.0)

    return ConstrainedResult(
        phi_profile=phi,
        phi_coeffs=phi_coeffs,
        phi_asymptotic=0.0,
        r=r,
        ricci_scalar=ricci,
        rho_required=rho_required,
        integrand_negative_mask=neg_mask,
        fit_residual=residual,
    )


def constrained_overrides(
    overrides: dict[str, object],
    *,
    phantom: bool = False,
) -> dict[str, object]:
    """Given a dict of params overrides that includes chi recipe coefficients,
    derive constraint-consistent phi coefficients and inject them, while
    locking K to constant and Pi to zero.

    Modifies and returns *overrides* in place.
    """
    num_bases = int(overrides.get("recipe_num_bases", 4))
    basis_width = float(overrides.get("recipe_basis_width", 1.0))
    basis_radius_max = float(overrides.get("recipe_basis_radius_max", 8.0))
    chi_asymptotic = float(overrides.get("recipe_chi_asymptotic", 1.0))
    K_const = float(overrides.get("recipe_K_asymptotic", 0.0))

    basis = RecipeBasis(
        num_bases=num_bases,
        basis_width=basis_width,
        basis_radius_max=basis_radius_max,
    )
    chi_coeffs = [
        float(overrides.get(f"recipe_chi_coeff_{n}", 0.0))
        for n in range(num_bases)
    ]

    result = constrained_phi(
        basis=basis,
        chi_asymptotic=chi_asymptotic,
        chi_coeffs=chi_coeffs,
        K_const=K_const,
        phantom=phantom,
    )

    for n in range(num_bases):
        overrides[f"recipe_K_coeff_{n}"] = 0.0
        overrides[f"recipe_Pi_coeff_{n}"] = 0.0
        overrides[f"recipe_phi_coeff_{n}"] = result.phi_coeffs[n]

    overrides["recipe_phi_asymptotic"] = result.phi_asymptotic
    overrides["recipe_Pi_asymptotic"] = 0.0

    return overrides
