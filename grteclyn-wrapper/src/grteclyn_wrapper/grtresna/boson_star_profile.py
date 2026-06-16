"""Radial equilibrium profiles for mini boson stars (lambda=0).

Phase 1 uses a Gaussian envelope with frequency omega ~ m (flat-space
mini-boson-star limit), tabulated on a radial grid and wrapped in a cubic
spline.  This matches the analytic profile painted by GRTresna BosonStarBH
(MatterParams Gaussian) and avoids the numerical stiffness of the full
self-gravitating ODE during pipeline bring-up.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline


@dataclass(frozen=True)
class BosonStarProfile:
    """Tabulated mini-boson-star radial data."""

    mass: float
    phi_c: float
    omega: float
    r: np.ndarray
    phi0: np.ndarray
    spline: CubicSpline
    adm_mass: float
    profile_width: float

    def eval_phi0(self, radius: float | np.ndarray) -> np.ndarray | float:
        r_arr = np.asarray(radius, dtype=float)
        scalar = bool(r_arr.ndim == 0)
        r_q = np.atleast_1d(r_arr)
        r_q = np.clip(r_q, self.r[0], self.r[-1])
        out = self.spline(r_q)
        return float(out[0]) if scalar else out

    def eval_dphi0_dr(self, radius: float | np.ndarray) -> np.ndarray | float:
        r_arr = np.asarray(radius, dtype=float)
        scalar = bool(r_arr.ndim == 0)
        r_q = np.atleast_1d(r_arr)
        r_q = np.clip(r_q, self.r[0], self.r[-1])
        out = self.spline.derivative(1)(r_q)
        return float(out[0]) if scalar else out


def _gaussian_phi0(r: np.ndarray, phi_c: float, width: float) -> np.ndarray:
    return phi_c * np.exp(-0.5 * (r / width) ** 2)


def _estimate_omega(mass: float, phi_c: float) -> float:
    """Mini-boson-star frequency (lambda=0).

    Phase 1 uses omega = m as the single source of truth shared with the
    GRTresna paint (BosonStarParams omega defaults to scalar_mass) and the
    post-solve Pi_im correction, so the constraint-solve ID and the evolved
    ID agree on the phase velocity.
    """
    del phi_c
    return mass


def solve_mini_boson_star(
    mass: float = 0.1,
    phi_c: float = 0.08,
    r_max: float = 80.0,
    n_points: int = 400,
    profile_width: float | None = None,
    omega_guess: float | None = None,
) -> BosonStarProfile:
    """Build a spline tabulation of a Gaussian mini-boson-star profile."""
    if mass <= 0.0:
        raise ValueError("mass must be positive")
    if phi_c <= 0.0:
        raise ValueError("phi_c must be positive")

    width = profile_width if profile_width is not None else max(8.0, 4.0 / mass)
    omega = omega_guess if omega_guess is not None else _estimate_omega(mass, phi_c)

    r = np.linspace(0.0, r_max, n_points)
    phi0 = _gaussian_phi0(r, phi_c, width)
    spline = CubicSpline(r, phi0, bc_type="natural")
    adm_mass = float(4.0 * math.pi * np.trapezoid(phi0 * phi0 * r * r, r))

    return BosonStarProfile(
        mass=mass,
        phi_c=phi_c,
        omega=float(omega),
        r=r,
        phi0=phi0,
        spline=spline,
        adm_mass=adm_mass,
        profile_width=width,
    )
