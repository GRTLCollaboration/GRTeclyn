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


def _sech_phi0(r: np.ndarray, phi_c: float, width: float) -> np.ndarray:
    """Bound boson-lump profile phi_c * sech(r / width).

    Regular at r=0 (sech(0)=1) and decays as 2 phi_c e^{-r/width} -- the correct
    EXPONENTIAL tail of a massive bound scalar, unlike a Gaussian which decays
    far too fast (e^{-r^2}) and squeezes the field into a non-equilibrium core.
    """
    return phi_c / np.cosh(r / width)


#: Lump-shape tags shared with GRTresna (BosonStarParams::lump_envelope) and the
#: GRTeclyn pump controller.  0 = Gaussian, 1 = tanh shell, 2 = sech bound lump.
PROFILE_GAUSSIAN = 0
PROFILE_TANH_SHELL = 1
PROFILE_SECH_BOUND = 2

# When omega >= mass the field is unbound (no exponential confinement); fall back
# to this finite decay scale so the width stays sane instead of diverging.
_UNBOUND_FALLBACK_KAPPA = 1.0 / 16.0


def bound_width(mass: float, omega: float) -> float:
    """Physical bound-state size R = 1 / sqrt(m^2 - omega^2).

    This is the decay length of a massive complex scalar oscillating at omega:
    the field cannot localize tighter than R without huge gradient energy, so a
    lump narrower than R is not a bound state and disperses.  Returns a large but
    finite width when omega >= mass (marginally/unbound).
    """
    kappa_sq = float(mass) * float(mass) - float(omega) * float(omega)
    kappa = math.sqrt(kappa_sq) if kappa_sq > 0.0 else _UNBOUND_FALLBACK_KAPPA
    return 1.0 / max(kappa, _UNBOUND_FALLBACK_KAPPA)


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
