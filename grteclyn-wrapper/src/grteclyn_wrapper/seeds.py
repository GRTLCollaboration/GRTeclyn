"""Known-solution seeds for optimizer warm-start and regression benchmarks.

Each seed provides a set of RadialRecipe parameter overrides that represent
a known, analytically understood spacetime.  These serve two purposes:

1. **Regression benchmarks**: flat Minkowski should score perfectly on
   constraints; the EB wormhole should reproduce the known wormhole
   dynamics.
2. **Optimizer warm-start**: seeding CMA-ES with a known good geometry
   starts the search in a physically meaningful region.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

from .constrained_recipe import RecipeBasis, fit_gaussian_basis


@dataclass(frozen=True)
class Seed:
    """A named initial-data seed with RadialRecipe parameter overrides."""

    name: str
    description: str
    overrides: dict[str, Any]
    chi_vector: list[float]


def flat_minkowski(
    num_bases: int = 4,
    basis_width: float = 1.0,
    basis_radius_max: float = 8.0,
) -> Seed:
    """Flat Minkowski spacetime: chi=1, alpha=1, K=0, phi=0, Pi=0.

    This is the trivial solution.  All constraints are exactly satisfied,
    and the evolution should remain static to machine precision.
    """
    overrides: dict[str, Any] = {
        "recipe_num_bases": num_bases,
        "recipe_basis_width": basis_width,
        "recipe_basis_radius_max": basis_radius_max,
        "recipe_chi_asymptotic": 1.0,
        "recipe_alpha_asymptotic": 1.0,
        "recipe_K_asymptotic": 0.0,
        "recipe_phi_asymptotic": 0.0,
        "recipe_Pi_asymptotic": 0.0,
    }
    chi_vec = []
    for n in range(num_bases):
        overrides[f"recipe_chi_coeff_{n}"] = 0.0
        overrides[f"recipe_alpha_coeff_{n}"] = 0.0
        overrides[f"recipe_K_coeff_{n}"] = 0.0
        overrides[f"recipe_phi_coeff_{n}"] = 0.0
        overrides[f"recipe_Pi_coeff_{n}"] = 0.0
        chi_vec.append(0.0)

    return Seed(
        name="flat_minkowski",
        description="Trivial flat spacetime. All constraints exactly satisfied.",
        overrides=overrides,
        chi_vector=chi_vec,
    )


def ellis_bronnikov_wormhole(
    b0: float = 0.5,
    num_bases: int = 8,
    basis_width: float = 0.5,
    basis_radius_max: float = 8.0,
    n_fit_points: int = 2048,
) -> Seed:
    """Ellis-Bronnikov wormhole in isotropic coordinates.

    The exact chi profile is:
        chi(r) = (4*r^2 / (4*r^2 + b0^2))^2

    The exact scalar field profile (phantom) is:
        phi(r) = (1/sqrt(4*pi)) * arctan((r - b0^2/(4*r)) / b0)

    These are fitted to the Gaussian basis for use in RadialRecipe.
    K = 0, Pi = 0, alpha = 1 (exact static solution).
    """
    basis = RecipeBasis(
        num_bases=num_bases,
        basis_width=basis_width,
        basis_radius_max=basis_radius_max,
    )

    r_min = basis_radius_max / n_fit_points
    r = np.linspace(r_min, 2.0 * basis_radius_max, n_fit_points)

    chi_exact = (4.0 * r**2 / (4.0 * r**2 + b0**2)) ** 2
    chi_coeffs, chi_residual = fit_gaussian_basis(r, chi_exact, basis, asymptotic=1.0)

    phi_exact = (1.0 / math.sqrt(4.0 * math.pi)) * np.arctan(
        (r - b0**2 / (4.0 * r)) / b0
    )
    phi_asymptotic = float(phi_exact[-1])
    phi_coeffs, phi_residual = fit_gaussian_basis(
        r, phi_exact, basis, asymptotic=phi_asymptotic,
    )

    overrides: dict[str, Any] = {
        "recipe_num_bases": num_bases,
        "recipe_basis_width": basis_width,
        "recipe_basis_radius_max": basis_radius_max,
        "recipe_chi_asymptotic": 1.0,
        "recipe_alpha_asymptotic": 1.0,
        "recipe_K_asymptotic": 0.0,
        "recipe_phi_asymptotic": phi_asymptotic,
        "recipe_Pi_asymptotic": 0.0,
    }
    for n in range(num_bases):
        overrides[f"recipe_chi_coeff_{n}"] = chi_coeffs[n]
        overrides[f"recipe_alpha_coeff_{n}"] = 0.0
        overrides[f"recipe_K_coeff_{n}"] = 0.0
        overrides[f"recipe_phi_coeff_{n}"] = phi_coeffs[n]
        overrides[f"recipe_Pi_coeff_{n}"] = 0.0

    return Seed(
        name="ellis_bronnikov",
        description=(
            f"Ellis-Bronnikov wormhole (b0={b0}) in isotropic coordinates. "
            f"Gaussian fit residuals: chi={chi_residual:.2e}, phi={phi_residual:.2e}."
        ),
        overrides=overrides,
        chi_vector=chi_coeffs,
    )


def schwarzschild_puncture(
    M: float = 0.5,
    num_bases: int = 8,
    basis_width: float = 0.5,
    basis_radius_max: float = 8.0,
    n_fit_points: int = 2048,
) -> Seed:
    """Schwarzschild black hole in isotropic (puncture) coordinates.

    The exact chi profile is:
        psi(r) = 1 + M/(2*r)
        chi(r) = psi^{-4} = (1 + M/(2*r))^{-4}

    This is a vacuum solution (no scalar field), K=0, time-symmetric.
    """
    basis = RecipeBasis(
        num_bases=num_bases,
        basis_width=basis_width,
        basis_radius_max=basis_radius_max,
    )

    r_min = basis_radius_max / n_fit_points
    r = np.linspace(r_min, 2.0 * basis_radius_max, n_fit_points)

    psi = 1.0 + M / (2.0 * r)
    chi_exact = psi ** (-4.0)
    chi_coeffs, chi_residual = fit_gaussian_basis(r, chi_exact, basis, asymptotic=1.0)

    overrides: dict[str, Any] = {
        "recipe_num_bases": num_bases,
        "recipe_basis_width": basis_width,
        "recipe_basis_radius_max": basis_radius_max,
        "recipe_chi_asymptotic": 1.0,
        "recipe_alpha_asymptotic": 1.0,
        "recipe_K_asymptotic": 0.0,
        "recipe_phi_asymptotic": 0.0,
        "recipe_Pi_asymptotic": 0.0,
    }
    for n in range(num_bases):
        overrides[f"recipe_chi_coeff_{n}"] = chi_coeffs[n]
        overrides[f"recipe_alpha_coeff_{n}"] = 0.0
        overrides[f"recipe_K_coeff_{n}"] = 0.0
        overrides[f"recipe_phi_coeff_{n}"] = 0.0
        overrides[f"recipe_Pi_coeff_{n}"] = 0.0

    return Seed(
        name="schwarzschild_puncture",
        description=(
            f"Schwarzschild BH (M={M}) in isotropic puncture coords. "
            f"Vacuum solution, no scalar field. "
            f"Gaussian fit residual: chi={chi_residual:.2e}."
        ),
        overrides=overrides,
        chi_vector=chi_coeffs,
    )


ALL_SEEDS: dict[str, type[...] | None] = {
    "flat_minkowski": None,
    "ellis_bronnikov": None,
    "schwarzschild_puncture": None,
}


def get_seed(name: str, **kwargs: Any) -> Seed:
    """Get a named seed with optional parameter overrides."""
    factories = {
        "flat_minkowski": flat_minkowski,
        "ellis_bronnikov": ellis_bronnikov_wormhole,
        "schwarzschild_puncture": schwarzschild_puncture,
    }
    if name not in factories:
        known = ", ".join(sorted(factories))
        raise ValueError(f"Unknown seed {name!r}. Known seeds: {known}")
    return factories[name](**kwargs)


def list_seeds() -> list[str]:
    """Return names of all available seeds."""
    return ["flat_minkowski", "ellis_bronnikov", "schwarzschild_puncture"]
