"""Flat-space Q-ball radial ODE solver (Fix 3).

Stationary profile φ₀(r) solves::

    φ₀'' + (2/r) φ₀' = (m² − ω²) φ₀ − λ φ₀³ + μ φ₀⁵

with φ₀'(0) = 0 and φ₀(∞) → 0.  Outward integration with a bisection search on
φ₀(0) finds the positive ground-state soliton (smallest tail amplitude before
numerical blow-up).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from .qball_couplings import QBallCouplings


@dataclass(frozen=True)
class QBallRadialProfile:
    """Tabulated flat-space Q-ball profile φ₀(r)."""

    mass: float
    lam: float
    mu: float
    omega: float
    phi_c: float
    r: NDArray[np.float64]
    phi0: NDArray[np.float64]
    spline: CubicSpline

    def eval_phi0(self, radius: float | NDArray[np.float64]) -> NDArray[np.float64] | float:
        r_arr = np.asarray(radius, dtype=np.float64)
        scalar = bool(r_arr.ndim == 0)
        r_q = np.atleast_1d(r_arr)
        r_q = np.clip(r_q, float(self.r[0]), float(self.r[-1]))
        out = self.spline(r_q)
        return float(out[0]) if scalar else out


def _ode_rhs_outward(
    r: float,
    state: NDArray[np.float64],
    *,
    kappa_sq: float,
    lam: float,
    mu: float,
) -> list[float]:
    phi, dphi = float(state[0]), float(state[1])
    d2phi = kappa_sq * phi - lam * phi**3 + mu * phi**5
    if r > 1.0e-10:
        d2phi -= (2.0 / r) * dphi
    return [dphi, d2phi]


def _integrate_outward(
    phi_c: float,
    *,
    kappa_sq: float,
    lam: float,
    mu: float,
    r_max: float,
) -> tuple[bool, float, NDArray[np.float64], NDArray[np.float64]]:
    """Return (ok, tail_phi, r_grid, phi_grid)."""
    if phi_c <= 0.0:
        return False, math.inf, np.array([0.0]), np.array([0.0])

    def rhs(r, y):
        return _ode_rhs_outward(r, y, kappa_sq=kappa_sq, lam=lam, mu=mu)

    sol = solve_ivp(
        rhs,
        (1.0e-8, r_max),
        [phi_c, 0.0],
        method="RK45",
        max_step=min(0.02, r_max / 1000.0),
        rtol=1.0e-8,
        atol=1.0e-10,
        dense_output=True,
    )
    if not np.all(np.isfinite(sol.y)):
        return False, math.inf, np.array([0.0]), np.array([0.0])

    n_eval = max(200, int(float(sol.t[-1]) * 12))
    r_grid = np.linspace(0.0, float(sol.t[-1]), n_eval)
    phi_grid = np.asarray(sol.sol(r_grid)[0], dtype=np.float64)
    phi_grid[0] = phi_c
    if np.any(phi_grid < -1.0e-6) or np.max(np.abs(phi_grid)) > max(10.0 * phi_c, 1.0):
        return False, math.inf, r_grid, phi_grid
    tail = float(phi_grid[-1])
    if tail > max(0.5, 2.0 * phi_c):
        return False, math.inf, r_grid, phi_grid
    return True, tail, r_grid, phi_grid


def solve_qball_radial_profile(
    couplings: QBallCouplings,
    *,
    r_max: float | None = None,
) -> QBallRadialProfile:
    """Find the flat-space Q-ball profile for *couplings* via outward shooting."""
    couplings.validate()
    if couplings.lam <= 0.0 or couplings.mu <= 0.0:
        raise ValueError("Q-ball ODE requires lam > 0 and mu > 0")

    mass, lam, mu, omega = couplings.mass, couplings.lam, couplings.mu, couplings.omega
    kappa_sq = mass * mass - omega * omega
    if kappa_sq <= 0.0:
        raise ValueError("omega must be < mass for a bound Q-ball profile")

    width = couplings.bound_state_width
    r_max = float(r_max if r_max is not None else max(100.0, 16.0 * width))

    phi_lo = max(1.0e-4, 0.5 * couplings.core_amplitude)
    phi_hi = max(phi_lo * 1.01, 1.05 * couplings.core_amplitude)

    # Ensure upper bracket hits blow-up.
    for _ in range(20):
        ok, _, _, _ = _integrate_outward(
            phi_hi, kappa_sq=kappa_sq, lam=lam, mu=mu, r_max=r_max
        )
        if not ok:
            break
        phi_hi = min(phi_hi * 1.05, 3.0 * couplings.core_amplitude)
    else:
        raise RuntimeError("Q-ball ODE shoot failed: could not find unstable upper bracket")

    best_phi_c = phi_lo
    best_r: NDArray[np.float64] | None = None
    best_phi: NDArray[np.float64] | None = None
    for _ in range(48):
        mid = 0.5 * (phi_lo + phi_hi)
        ok, _, r_grid, phi_grid = _integrate_outward(
            mid, kappa_sq=kappa_sq, lam=lam, mu=mu, r_max=r_max
        )
        if ok:
            phi_lo = mid
            best_phi_c = mid
            best_r, best_phi = r_grid, phi_grid
        else:
            phi_hi = mid
        if phi_hi - phi_lo < 1.0e-4:
            break

    if best_r is None or best_phi is None:
        raise RuntimeError("Q-ball ODE shoot failed: no stable outward integration")

    spline = CubicSpline(best_r, best_phi, bc_type="natural")
    return QBallRadialProfile(
        mass=mass,
        lam=lam,
        mu=mu,
        omega=omega,
        phi_c=best_phi_c,
        r=best_r,
        phi0=best_phi,
        spline=spline,
    )


@lru_cache(maxsize=16)
def cached_qball_radial_profile(
    mass: float,
    lam: float,
    mu: float,
    omega: float,
) -> QBallRadialProfile:
    """Memoised profile lookup for grid painting hot paths."""
    return solve_qball_radial_profile(
        QBallCouplings(mass=mass, lam=lam, mu=mu, omega=omega, phi_core=0.0)
    )


def profile_for_lump(
    lump: dict,
    *,
    mass: float,
    lam: float,
    mu: float,
    omega: float,
) -> QBallRadialProfile:
    """Build or fetch the ODE profile for a lump dict (``profile == ODE``)."""
    amp = float(lump.get("amp", 0.0))
    base = cached_qball_radial_profile(mass, lam, mu, omega)
    if amp <= 0.0 or math.isclose(amp, base.phi_c, rel_tol=0.0, abs_tol=1.0e-12):
        return base
    scale = amp / base.phi_c
    r = base.r
    phi0 = base.phi0 * scale
    return QBallRadialProfile(
        mass=base.mass,
        lam=base.lam,
        mu=base.mu,
        omega=base.omega,
        phi_c=amp,
        r=r,
        phi0=phi0,
        spline=CubicSpline(r, phi0, bc_type="natural"),
    )
