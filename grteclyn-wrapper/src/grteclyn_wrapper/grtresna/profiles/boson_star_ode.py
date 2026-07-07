"""Self-gravitating boson star radial ODE solver (isotropic coordinates).

A stationary, spherically-symmetric boson star is a genuine equilibrium of the
*coupled* Einstein-Klein-Gordon system: gravity provides the binding, so no
external pump well is required and the configuration does not disperse when
evolved.

Metric ansatz -- **isotropic** form (spatially conformally flat), G=1::

    ds^2 = -alpha(r)^2 dt^2 + psi(r)^4 (dr^2 + r^2 dOmega^2)

This is exactly the slice GRTresna reconstructs (gamma_ij = psi^4 delta_ij on a
conformally-flat, maximally-sliced K=0 hypersurface), so ``phi0`` tabulated
against the isotropic radius ``r`` here is stationary initial data once GRTresna
solves the Lichnerowicz constraint for its own conformal factor.

Complex scalar with rotating phase::

    Phi(t, r) = phi0(r) * exp(-i * omega * t)

For a single star at rest the momentum density vanishes, so K_ij = 0 and the
geometry equations reduce to the Hamiltonian (Lichnerowicz) constraint plus the
maximal-slicing lapse equation.  With potential
``U = 1/2 m^2 phi^2 - 1/4 lam phi^4 + 1/6 mu phi^6`` and
``dU = m^2 phi - lam phi^3 + mu phi^5`` the coupled radial system is::

    psi''   = -(2/r) psi' - 2*pi * psi^5 * rho
    alpha'' = -(2/r + 2 psi'/psi) alpha' + 4*pi * psi^4 * alpha * (rho + S)
    phi0''  = -(2/r + alpha'/alpha + 2 psi'/psi) phi0'
              + psi^4 * (dU - omega^2/alpha^2 * phi0)

with the scalar energy density and stress trace (Phi = phi0 e^{-i omega t}
treated as two real components each carrying the standard 1/2 kinetic factor,
matching GRTeclyn's ``ComplexScalarField``)::

    rho   = 1/2 * omega^2/alpha^2 * phi0^2 + 1/2 * psi^-4 * phi0'^2 + U
    rho+S = 2 * omega^2/alpha^2 * phi0^2 - 2 * U

Regularity at the origin::

    phi0'(0) = 0,  psi'(0) = 0,  alpha'(0) = 0

and asymptotic flatness + localization::

    psi(inf) -> 1 + M/(2r),  alpha(inf) -> (1 - M/2r)/(1 + M/2r),  phi0(inf) -> 0

Solution strategy (single-parameter shooting + two exact rescalings)
--------------------------------------------------------------------
Integrate outward from a regular start with the gauge ``psi(0) = alpha(0) = 1``.
Only ``omega`` is the shooting eigenvalue: we bisect between the *overshoot*
branch (omega too high -> phi0 crosses zero and diverges) and the *undershoot*
branch (omega too low -> phi0 turns back up toward a nonzero asymptote).  This is
the same shooting philosophy as the flat-space ``qball_ode.py`` but on omega.

The integration gauge does not satisfy ``psi(inf)=1`` or ``alpha(inf)=1``, but
the system has two exact scaling symmetries that restore them without re-solving:

  1. ``(alpha -> c*alpha, omega -> c*omega)`` leaves the equations invariant
     (they depend on alpha only through omega^2/alpha^2 and the linear-in-alpha
     lapse equation).  Choosing ``c = 1/alpha_inf`` gives ``alpha_inf = 1`` and
     the physical frequency ``omega_phys = omega_int / alpha_inf < m``.
  2. ``(r -> s*r, psi -> psi/sqrt(s))`` with ``s = psi_inf^2`` leaves the
     equations invariant and gives ``psi_inf = 1``; the physical isotropic radius
     is ``r_phys = psi_inf^2 * r`` and ``psi_phys = psi / psi_inf``.  The proper
     (areal) radius ``psi^2 r`` is preserved by this map.

The ADM mass is read from the exterior conformal factor ``psi -> 1 + M/(2r)``.
For the mini boson star (``lam = mu = 0``) this reproduces the classic Kaup /
Ruffini-Bonazzola soliton with maximum ``M_max * m ~= 0.633`` (a coordinate-
independent check on the whole solver).
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

# State vector layout for the outward integration:
#   y[0] = phi0,  y[1] = dphi0/dr,  y[2] = psi,  y[3] = dpsi/dr,
#   y[4] = alpha, y[5] = dalpha/dr
_PHI = 0
_DPHI = 1
_PSI = 2
_DPSI = 3
_ALPHA = 4
_DALPHA = 5

# Geometric units with G = 1; the 2*pi / 4*pi factors are the standard
# York-Lichnerowicz normalisation for gamma_ij = psi^4 delta_ij, K = 0.
_TWO_PI = 2.0 * math.pi
_FOUR_PI = 4.0 * math.pi

# Tolerance used to decide that phi0 has "crossed zero" (overshoot) on a noisy
# integration.  Absolute, in field units.
_ZERO_CROSS_TOL = 1.0e-8


@dataclass(frozen=True)
class BosonStarODEProfile:
    """Tabulated self-gravitating boson-star radial data (isotropic radius)."""

    mass: float
    lam: float
    mu: float
    omega: float
    """Gravitational eigenvalue (single source of truth for paint + evolution)."""
    phi_c: float
    """Central field amplitude that produced this star."""
    r: NDArray[np.float64]
    """Isotropic radius -- the coordinate GRTresna reconstructs (gamma = psi^4 delta)."""
    phi0: NDArray[np.float64]
    psi: NDArray[np.float64]
    """Conformal factor psi(r) (spatial metric gamma_ij = psi^4 delta_ij)."""
    alpha: NDArray[np.float64]
    """Lapse alpha(r)."""
    adm_mass: float
    compactness: float
    """Rough compactness proxy 2*M / R_99 (R_99 = radius enclosing 99% of |phi|)."""
    spline_phi0: CubicSpline
    spline_psi: CubicSpline
    spline_alpha: CubicSpline

    def eval_phi0(self, radius: float | NDArray[np.float64]) -> NDArray[np.float64] | float:
        r_arr = np.asarray(radius, dtype=np.float64)
        scalar = bool(r_arr.ndim == 0)
        r_q = np.atleast_1d(r_arr)
        r_q = np.clip(r_q, float(self.r[0]), float(self.r[-1]))
        out = self.spline_phi0(r_q)
        return float(out[0]) if scalar else out

    def eval_dphi0_dr(self, radius: float | NDArray[np.float64]) -> NDArray[np.float64] | float:
        r_arr = np.asarray(radius, dtype=np.float64)
        scalar = bool(r_arr.ndim == 0)
        r_q = np.atleast_1d(r_arr)
        r_q = np.clip(r_q, float(self.r[0]), float(self.r[-1]))
        out = self.spline_phi0.derivative(1)(r_q)
        return float(out[0]) if scalar else out


def _potential(phi: float, *, mass: float, lam: float, mu: float) -> float:
    """U(phi) = 1/2 m^2 phi^2 - 1/4 lam phi^4 + 1/6 mu phi^6."""
    phi2 = phi * phi
    return 0.5 * mass * mass * phi2 - 0.25 * lam * phi2 * phi2 + (mu / 6.0) * phi2 * phi2 * phi2


def _d_potential(phi: float, *, mass: float, lam: float, mu: float) -> float:
    """dU/dphi = m^2 phi - lam phi^3 + mu phi^5."""
    phi2 = phi * phi
    return mass * mass * phi - lam * phi2 * phi + mu * phi2 * phi2 * phi


def _ode_rhs(
    r: float,
    state: NDArray[np.float64],
    *,
    mass: float,
    lam: float,
    mu: float,
    omega: float,
) -> list[float]:
    """RHS of the coupled Einstein-Klein-Gordon radial system (isotropic form)."""
    phi = float(state[_PHI])
    dphi = float(state[_DPHI])
    psi = float(state[_PSI])
    dpsi = float(state[_DPSI])
    alpha = float(state[_ALPHA])
    dalpha = float(state[_DALPHA])

    # Guard against integrator excursions producing non-positive metric
    # components (which would blow up the logarithmic terms).
    psi = max(psi, 1.0e-12)
    alpha = max(alpha, 1.0e-12)

    psi2 = psi * psi
    psi4 = psi2 * psi2
    inv_psi4 = 1.0 / psi4
    phi2 = phi * phi
    u_val = _potential(phi, mass=mass, lam=lam, mu=mu)
    du_val = _d_potential(phi, mass=mass, lam=lam, mu=mu)
    omega_over_alpha_sq = (omega / alpha) ** 2

    # Scalar energy density and stress trace (K_ij = 0 single star at rest).
    # NB: the two real components (phi1, phi2) of Phi = phi0 e^{-i omega t} each
    # carry the standard 1/2 kinetic factor, matching GRTeclyn's
    # ComplexScalarField::compute_emtensor (rho = 1/2 sum Pi_k^2 + 1/2 sum
    # |grad phi_k|^2 + V).  With sum Pi_k^2 = (omega/alpha)^2 phi0^2 this gives:
    #   rho    = 1/2 omega^2/alpha^2 phi^2 + 1/2 psi^-4 phi'^2 + U
    #   rho+S  = 2 omega^2/alpha^2 phi^2 - 2 U
    rho = 0.5 * omega_over_alpha_sq * phi2 + 0.5 * inv_psi4 * dphi * dphi + u_val
    rho_plus_s = 2.0 * omega_over_alpha_sq * phi2 - 2.0 * u_val

    if r > 1.0e-12:
        inv_r = 1.0 / r
        psi_log_deriv = dpsi / psi
        # Hamiltonian (Lichnerowicz) constraint: psi'' + (2/r) psi' = -2 pi psi^5 rho.
        d2psi = -2.0 * inv_r * dpsi - _TWO_PI * psi4 * psi * rho
        # Maximal-slicing lapse: alpha'' + (2/r + 2 psi'/psi) alpha' = 4 pi psi^4 alpha (rho+S).
        d2alpha = (
            -(2.0 * inv_r + 2.0 * psi_log_deriv) * dalpha
            + _FOUR_PI * psi4 * alpha * rho_plus_s
        )
        # Klein-Gordon: phi'' + (2/r + alpha'/alpha + 2 psi'/psi) phi' = psi^4 (dU - omega^2/alpha^2 phi).
        friction = 2.0 * inv_r + dalpha / alpha + 2.0 * psi_log_deriv
        d2phi = -friction * dphi + psi4 * (du_val - omega_over_alpha_sq * phi)
    else:
        # Regular origin: the 1/r friction terms vanish because the first
        # derivatives are zero there; the O(r^2) series is recovered within a
        # few integrator steps from r0 ~ 1e-4.
        d2psi = -_TWO_PI * psi4 * psi * rho / 3.0
        d2alpha = _FOUR_PI * psi4 * alpha * rho_plus_s / 3.0
        d2phi = psi4 * (du_val - omega_over_alpha_sq * phi) / 3.0

    return [dphi, d2phi, dpsi, d2psi, dalpha, d2alpha]


def _shoot(
    omega: float,
    phi_c: float,
    *,
    mass: float,
    lam: float,
    mu: float,
    r_max: float,
) -> tuple[str, object]:
    """Integrate the coupled system outward from a regular r~0 start.

    ``omega`` here is the *integration-frame* frequency with the gauge
    ``alpha(0) = psi(0) = 1``.  The physical frequency ``omega_phys =
    omega/alpha_inf < m`` is recovered by rescaling in the caller.

    Classifies the trajectory as:

      * ``"over"``  -- omega too high; the field is over-bound and phi0 descends
                      through zero (would go to -inf).
      * ``"under"`` -- omega too low; binding is too weak and phi0 grows away
                      (runs to +inf) or turns back up toward a nonzero asymptote.

    Returns ``(kind, sol)``.
    """
    if phi_c <= 0.0:
        raise ValueError("phi_c must be positive")
    if omega <= 0.0:
        raise ValueError(f"omega ({omega}) must be positive")

    r0 = 1.0e-4
    y0 = [phi_c, 0.0, 1.0, 0.0, 1.0, 0.0]

    def rhs(r, y):
        return _ode_rhs(r, y, mass=mass, lam=lam, mu=mu, omega=omega)

    # Overshoot: phi0 crosses zero (descending) -> over.
    def event_zero(r, y):
        return y[_PHI] + _ZERO_CROSS_TOL

    event_zero.terminal = True
    event_zero.direction = -1.0

    # Undershoot: phi0 grows past its central value -> ran away up (under).
    grow_cap = 1.5 * phi_c

    def event_grow(r, y):
        return y[_PHI] - grow_cap

    event_grow.terminal = True
    event_grow.direction = 1.0

    # Undershoot (near-separatrix): phi0' returns to >= 0 after a genuine
    # descent, meaning the field turns back up toward a nonzero asymptote.
    descended = 0.999 * phi_c

    def event_turn(r, y):
        return y[_DPHI] if (r > 1.0e-3 and y[_PHI] < descended) else -1.0

    event_turn.terminal = True
    event_turn.direction = 1.0

    sol = solve_ivp(
        rhs,
        (r0, r_max),
        y0,
        method="DOP853",
        max_step=max(1.0e-3, r_max / 4000.0),
        rtol=1.0e-9,
        atol=1.0e-11,
        dense_output=True,
        events=(event_zero, event_grow, event_turn),
    )
    fired_zero = sol.t_events[0].size > 0
    fired_grow = sol.t_events[1].size > 0
    fired_turn = sol.t_events[2].size > 0
    if fired_zero:
        return "over", sol
    if fired_grow or fired_turn:
        return "under", sol
    # Neither event fired within r_max: classify by the asymptotic state.
    phi_end = float(sol.y[_PHI, -1])
    dphi_end = float(sol.y[_DPHI, -1])
    if phi_end < 0.0:
        return "over", sol  # crossed zero without the event catching it
    if dphi_end > 0.0 and phi_end > 0.01 * phi_c:
        return "under", sol  # turned back up toward a nonzero asymptote
    if phi_end < 0.01 * phi_c:
        # Decayed to near zero -- near-eigenvalue.  Treat as over so bisection
        # keeps the tightest *under* trajectory as the kept solution.
        return "over", sol
    return "under", sol  # still positive, decreasing, but not yet decayed


def _asymptotic_constant(f_m: float, df_m: float, r_m: float) -> float:
    """Extract the r->inf constant of a field behaving as ``f = A + B/r``.

    From ``f(r) = A + B/r`` we have ``f' = -B/r^2``, so ``B = -r^2 f'`` and
    ``A = f + r f'``.  Used to read ``psi_inf`` and ``alpha_inf`` at a finite
    match radius where the metric is already near-Schwarzschild.
    """
    return f_m + r_m * df_m


def _build_profile(
    sol,
    *,
    omega_int: float,
    phi_c: float,
    mass: float,
    lam: float,
    mu: float,
    r_max: float,
) -> BosonStarODEProfile:
    """Sample the near-separatrix solution, rescale to the physical isotropic
    frame, and attach the analytic exterior.

    The ODE is integrated in the gauge ``psi(0)=alpha(0)=1``; two exact scaling
    symmetries then restore ``alpha_inf=1`` and ``psi_inf=1`` (see module
    docstring):

      * ``omega_phys = omega_int / alpha_inf`` (with ``alpha -> alpha/alpha_inf``)
      * ``r_phys = psi_inf^2 * r``, ``psi_phys = psi / psi_inf``.
    """
    r_end = float(sol.t[-1])
    n = max(600, int(r_end * 60))
    r_core = np.linspace(float(sol.t[0]), r_end, n)
    y = sol.sol(r_core)
    phi = np.asarray(y[_PHI], dtype=np.float64)
    dphi = np.asarray(y[_DPHI], dtype=np.float64)
    psi = np.asarray(y[_PSI], dtype=np.float64)
    dpsi = np.asarray(y[_DPSI], dtype=np.float64)
    alpha = np.asarray(y[_ALPHA], dtype=np.float64)
    dalpha = np.asarray(y[_DALPHA], dtype=np.float64)

    # Ensure the origin sample is exactly the boundary values.
    phi[0] = phi_c
    dphi[0] = 0.0
    psi[0] = 1.0
    dpsi[0] = 0.0
    alpha[0] = 1.0
    dalpha[0] = 0.0

    # --- Match radius: first point where the field has decayed to <5% of the
    # core while still descending.  Beyond this we attach the analytic tail so
    # the near-separatrix metric run-away never enters the tabulation.
    core = phi_c
    descended = (phi < 0.05 * core) & (dphi < 0.0) & (r_core > 1.0e-3)
    if np.any(descended):
        i_match = int(np.where(descended)[0][0])
    else:
        desc = (dphi < 0.0) & (r_core > 1.0e-3) & (phi > 0.0)
        i_match = int(np.where(desc)[0][-1]) if np.any(desc) else max(1, len(r_core) // 2)

    r_m = float(r_core[i_match])

    # --- Read the asymptotic constants at the match radius (psi, alpha ~ A+B/r).
    psi_inf = _asymptotic_constant(float(psi[i_match]), float(dpsi[i_match]), r_m)
    alpha_inf = _asymptotic_constant(float(alpha[i_match]), float(dalpha[i_match]), r_m)
    psi_inf = max(psi_inf, 1.0e-12)
    alpha_inf = max(alpha_inf, 1.0e-12)

    # --- Rescale to the physical (asymptotically flat) isotropic frame. --------
    # (1) alpha -> alpha/alpha_inf, omega -> omega/alpha_inf  (alpha_inf -> 1).
    # (2) r -> psi_inf^2 r, psi -> psi/psi_inf                (psi_inf  -> 1).
    omega_phys = omega_int / alpha_inf
    r_scale = psi_inf * psi_inf

    r_core_keep = r_core[: i_match + 1] * r_scale
    phi_core_keep = phi[: i_match + 1].copy()
    psi_core_keep = psi[: i_match + 1] / psi_inf
    alpha_core_keep = alpha[: i_match + 1] / alpha_inf

    r_m_phys = float(r_core_keep[-1])
    phi_m = float(phi_core_keep[-1])

    # --- ADM mass from the exterior conformal factor psi -> 1 + M/(2 r). -------
    psi_m = float(psi_core_keep[-1])
    adm_mass = float(max(2.0 * r_m_phys * (psi_m - 1.0), 0.0))

    kappa = math.sqrt(max(mass * mass - omega_phys * omega_phys, 1.0e-30))

    # Exterior: analytic Yukawa tail phi ~ e^{-kappa r}/r matched at r_m, plus the
    # isotropic Schwarzschild metric psi = 1 + M/2r, alpha = (1-M/2r)/(1+M/2r).
    tail_end = max(r_max, r_m_phys + 20.0 / max(kappa, 1.0e-6))
    n_tail = max(200, int((tail_end - r_m_phys) * 20))
    r_tail = np.linspace(r_m_phys, tail_end, n_tail)[1:]
    phi_tail = phi_m * (r_m_phys / r_tail) * np.exp(-kappa * (r_tail - r_m_phys))
    half_m_over_r = 0.5 * adm_mass / r_tail
    psi_tail = 1.0 + half_m_over_r
    alpha_tail = (1.0 - half_m_over_r) / (1.0 + half_m_over_r)

    r_full = np.concatenate([r_core_keep, r_tail])
    phi_full = np.concatenate([phi_core_keep, phi_tail])
    psi_full = np.concatenate([psi_core_keep, psi_tail])
    alpha_full = np.concatenate([alpha_core_keep, alpha_tail])

    # Prepend the exact origin so eval_phi0(0) is well defined.
    if r_full[0] > 0.0:
        r_full = np.concatenate([[0.0], r_full])
        phi_full = np.concatenate([[phi_c], phi_full])
        psi_full = np.concatenate([[float(psi_core_keep[0])], psi_full])
        alpha_full = np.concatenate([[float(alpha_core_keep[0])], alpha_full])

    # --- Compactness proxy: 2M / R_99. ---------------------------------------
    # R_99 is the isotropic radius enclosing 99% of the integrated |phi0|*r^2.
    integrand = np.abs(phi_full) * r_full ** 2
    cumulative = np.cumsum(integrand) * np.gradient(r_full)
    total = cumulative[-1] if cumulative.size else 1.0
    r99_idx = np.searchsorted(cumulative, 0.99 * total)
    r99 = float(r_full[min(r99_idx, r_full.size - 1)])
    compactness = (2.0 * adm_mass / r99) if r99 > 0.0 else 0.0

    spline_phi0 = CubicSpline(r_full, phi_full, bc_type="natural")
    spline_psi = CubicSpline(r_full, psi_full, bc_type="natural")
    spline_alpha = CubicSpline(r_full, alpha_full, bc_type="natural")

    return BosonStarODEProfile(
        mass=mass,
        lam=lam,
        mu=mu,
        omega=omega_phys,
        phi_c=phi_c,
        r=r_full,
        phi0=phi_full,
        psi=psi_full,
        alpha=alpha_full,
        adm_mass=adm_mass,
        compactness=compactness,
        spline_phi0=spline_phi0,
        spline_psi=spline_psi,
        spline_alpha=spline_alpha,
    )


def solve_selfgrav_boson_star(
    couplings: QBallCouplings,
    *,
    phi_c: float,
    r_max: float | None = None,
    omega_bracket: tuple[float, float] | None = None,
) -> BosonStarODEProfile:
    """Find the self-gravitating boson star profile for *couplings* at central
    amplitude *phi_c*, in isotropic coordinates.

    The frequency ``omega`` is the gravitational eigenvalue, found by bisection
    between the overshoot (omega too high -> phi0 crosses zero) and undershoot
    (omega too low -> phi0 turns back up) branches.

    For the mini boson star (``lam = mu = 0``) the result is the classic Kaup /
    Ruffini-Bonazzola soliton; the maximum ADM mass over the central-amplitude
    branch is ``M_max * m ~= 0.633`` (literature).  For the sextic family the
    same machinery handles a self-interacting star.
    """
    # NOTE: we do NOT call couplings.validate() -- its omega-band check applies
    # to the flat-space Q-ball where omega is an input.  For a self-gravitating
    # star omega is the eigenvalue we solve for, so only the potential couplings
    # (mass, lam, mu) are validated here.
    mass, lam, mu = couplings.mass, couplings.lam, couplings.mu
    if mass <= 0.0:
        raise ValueError("mass must be positive")
    if lam < 0.0 or mu < 0.0:
        raise ValueError("lam and mu must be non-negative")
    if lam > 0.0 and mu <= 0.0:
        raise ValueError("mu > 0 required when lam > 0 (sextic stabiliser)")
    if phi_c <= 0.0:
        raise ValueError("phi_c must be positive")

    r_max = float(r_max if r_max is not None else max(80.0, 30.0 / mass))

    def kind(omega: float) -> str:
        return _shoot(omega, phi_c, mass=mass, lam=lam, mu=mu, r_max=r_max)[0]

    # --- Bracket the integration-frame eigenvalue omega_int. ----------------
    # Gauge alpha(0)=psi(0)=1.  Ordering (same sign convention as the flat-space
    # solver):
    #   * omega_int too LOW  -> weak binding -> field grows away    -> "under".
    #   * omega_int too HIGH -> over-bound   -> field crosses zero  -> "over".
    if omega_bracket is not None:
        omega_lo, omega_hi = omega_bracket
    else:
        omega_lo = 0.5 * mass
        omega_hi = 1.5 * mass

    for _ in range(60):
        if kind(omega_lo) == "under":
            break
        omega_lo *= 0.8
        if omega_lo <= 1.0e-4 * mass:
            raise RuntimeError(
                "self-grav boson star: could not bracket the under (low-omega) side"
            )
    else:
        raise RuntimeError("self-grav boson star: could not bracket the under side")

    for _ in range(60):
        if kind(omega_hi) == "over":
            break
        omega_hi *= 1.25
        if omega_hi >= 1.0e4 * mass:
            raise RuntimeError(
                "self-grav boson star: could not bracket the over (high-omega) side"
            )
    else:
        raise RuntimeError("self-grav boson star: could not bracket the over side")

    # --- Bisect on omega_int to ~machine precision. -------------------------
    # Keep the tightest *under* trajectory (just below the eigenvalue): it tracks
    # the true soliton furthest down before the near-separatrix breaks.
    best_omega = omega_lo
    best_sol = None
    for _ in range(120):
        mid = 0.5 * (omega_lo + omega_hi)
        k, sol = _shoot(mid, phi_c, mass=mass, lam=lam, mu=mu, r_max=r_max)
        if k == "over":
            omega_hi = mid
        else:
            omega_lo = mid
            best_omega = mid
            best_sol = sol
        if omega_hi - omega_lo <= 1.0e-13 * max(1.0, mass):
            break

    if best_sol is None:
        _, best_sol = _shoot(omega_lo, phi_c, mass=mass, lam=lam, mu=mu, r_max=r_max)
        best_omega = omega_lo

    return _build_profile(
        best_sol,
        omega_int=best_omega,
        phi_c=phi_c,
        mass=mass,
        lam=lam,
        mu=mu,
        r_max=r_max,
    )


@lru_cache(maxsize=32)
def cached_selfgrav_profile(
    mass: float,
    lam: float,
    mu: float,
    phi_c: float,
) -> BosonStarODEProfile:
    """Memoised self-gravitating profile lookup for grid painting hot paths.

    ``omega`` is *not* a cache key because it is the eigenvalue we solve for;
    the cache is keyed on the physical inputs (couplings + central amplitude).
    """
    return solve_selfgrav_boson_star(
        QBallCouplings(mass=mass, lam=lam, mu=mu, omega=0.5 * mass, phi_core=phi_c),
        phi_c=phi_c,
    )


def profile_for_lump_selfgrav(
    lump: dict,
    *,
    mass: float,
    lam: float,
    mu: float,
) -> BosonStarODEProfile:
    """Build or fetch the self-gravitating ODE profile for a lump dict.

    The lump's ``amp`` (if present and nonzero) sets the central amplitude
    ``phi_c``; otherwise a default ``phi_c`` appropriate to the couplings is
    used (the flat-space equilibrium core for sextic potentials, or 0.1 for the
    mini star).
    """
    phi_c = float(lump.get("amp", 0.0))
    if phi_c <= 0.0:
        phi_c = default_selfgrav_phi_c(lam, mu)
    return cached_selfgrav_profile(mass, lam, mu, phi_c)


def default_selfgrav_phi_c(lam: float, mu: float) -> float:
    """Default central amplitude when a lump does not specify one.

    Mini star (``lam = mu = 0``) uses 0.1; the sextic family uses half the
    flat-space equilibrium core amplitude ``sqrt(3 lam / 4 mu)`` to stay on the
    stable (rising-mass) branch.
    """
    if lam > 0.0 and mu > 0.0:
        return 0.5 * math.sqrt(3.0 * lam / (4.0 * mu))
    return 0.1


def sech_width_for_profile(profile: BosonStarODEProfile) -> float:
    """Best sech ``1/cosh(r/w)`` width matching a self-gravitating profile.

    The closed-loop transport pump can only drive the field toward an analytic
    Gaussian/sech target (no tabulated target on the GPU), so we pick the sech
    width whose half-maximum coincides with the true profile's:
    ``phi0(r_half) = phi0(0)/2`` and ``sech(r_half/w) = 1/2`` give
    ``w = r_half / arccosh(2)``.  This matches the star's core far better than
    the tail-only length ``1/sqrt(m^2 - omega^2)`` (they differ by 10-25%).
    """
    core = float(profile.phi0[0])
    if core <= 0.0:
        raise ValueError("profile has non-positive central amplitude")
    r = np.asarray(profile.r, dtype=np.float64)
    phi = np.asarray(profile.phi0, dtype=np.float64)
    half = 0.5 * core
    below = np.where(phi <= half)[0]
    if below.size == 0:
        # Never reaches half the core within the table: fall back to R_99-ish.
        r_half = float(r[-1])
    else:
        i = int(below[0])
        if i == 0:
            r_half = float(r[0])
        else:
            # Linearly interpolate the crossing between r[i-1] and r[i].
            r0, r1 = float(r[i - 1]), float(r[i])
            p0, p1 = float(phi[i - 1]), float(phi[i])
            frac = 0.0 if p0 == p1 else (p0 - half) / (p0 - p1)
            r_half = r0 + frac * (r1 - r0)
    return float(r_half / math.acosh(2.0))


def selfgrav_pump_target_width(
    mass: float, lam: float, mu: float, phi_c: float
) -> float:
    """Sech pump-target width for the self-gravitating star at *phi_c*.

    Convenience wrapper used by the evolution wiring: solves (or fetches) the
    self-gravitating profile and returns the sech width that best matches it.
    """
    if phi_c <= 0.0:
        phi_c = default_selfgrav_phi_c(lam, mu)
    profile = cached_selfgrav_profile(mass, lam, mu, phi_c)
    return sech_width_for_profile(profile)
