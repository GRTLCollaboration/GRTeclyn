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


def _shoot(
    phi_c: float,
    *,
    kappa_sq: float,
    lam: float,
    mu: float,
    r_max: float,
):
    """Integrate φ₀ outward from φ(0)=φ_c and classify the trajectory.

    In the mechanical analogy (particle at position φ, "time" r, friction 2/r) the
    soliton is the knife-edge that comes to rest on the φ=0 hilltop as r→∞.  These
    sextic couplings also carry a competing false-vacuum minimum at φ_fv>0, so the
    two off-separatrix outcomes are:

      - **over** (φ_c too large): φ either descends through zero (→ −∞) *or*, for
        φ_c above the φ''(0)<0 window, climbs away to +∞.  Both mean "too much".
      - **under** (φ_c too small): φ descends to a positive minimum and turns back
        up, settling onto the false vacuum φ_fv instead of reaching 0.

    Returns ``(kind, sol)`` with ``kind ∈ {"over", "under"}``.
    """
    if phi_c <= 0.0:
        raise ValueError("phi_c must be positive")

    def rhs(r, y):
        return _ode_rhs_outward(r, y, kappa_sq=kappa_sq, lam=lam, mu=mu)

    grow_cap = max(2.0 * phi_c, 1.0e-6)

    def event_zero(r, y):  # φ crosses zero downward -> over (→ −∞)
        return y[0]

    event_zero.terminal = True
    event_zero.direction = -1.0

    def event_grow(r, y):  # φ climbs away above the core -> over (→ +∞)
        return y[0] - grow_cap

    event_grow.terminal = True
    event_grow.direction = 1.0

    descended = 0.999 * phi_c

    def event_turn(r, y):  # φ' returns to 0 after a genuine descent -> under
        # Only count a turning point once φ has actually dropped below φ_c; a
        # φ_c just above the window rises from r=0 and must NOT read as a turn.
        return y[1] if (r > 1.0e-6 and y[0] < descended) else -1.0

    event_turn.terminal = True
    event_turn.direction = 1.0

    sol = solve_ivp(
        rhs,
        (1.0e-8, r_max),
        [phi_c, 0.0],
        method="DOP853",
        max_step=max(1.0e-3, r_max / 4000.0),
        rtol=1.0e-10,
        atol=1.0e-12,
        dense_output=True,
        events=(event_zero, event_grow, event_turn),
    )
    fired_zero = sol.t_events[0].size > 0
    fired_grow = sol.t_events[1].size > 0
    fired_turn = sol.t_events[2].size > 0
    if fired_zero or fired_grow:
        return "over", sol
    if fired_turn:
        return "under", sol
    # Neither fired within r_max: asymptotically decaying (near-separatrix) or
    # settled flat -> treat as the under side (a small monotone tail).
    return "under", sol


def _build_localized_profile(
    phi_c: float,
    sol,
    *,
    kappa: float,
    r_max: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Truncate the near-separatrix shoot and attach the analytic Yukawa tail.

    The large-r linearisation φ'' + (2/r)φ' − κ²φ = 0 has the decaying solution
    φ ∝ e^{−κ r}/r, i.e. d(ln φ)/dr = −(κ + 1/r).  We match the integrated core
    to that tail at the radius where the measured log-derivative first reaches the
    asymptotic value, then continue analytically so φ₀(r) decays monotonically to
    zero instead of turning back up.
    """
    r_end = float(sol.t[-1])
    n = max(400, int(r_end * 40))
    r_core = np.linspace(0.0, r_end, n)
    y = sol.sol(r_core)
    phi = np.asarray(y[0], dtype=np.float64)
    dphi = np.asarray(y[1], dtype=np.float64)
    phi[0] = phi_c

    # Candidate match points: φ>0 and decreasing (φ'<0).  Match where the measured
    # log-derivative meets the Yukawa value −(κ + 1/r).
    with np.errstate(divide="ignore", invalid="ignore"):
        logderiv = dphi / phi
    target = -(kappa + 1.0 / np.maximum(r_core, 1.0e-9))
    valid = (phi > 0.0) & (dphi < 0.0) & (r_core > 0.0)
    gap = logderiv - target  # >0 in the core, → 0 approaching the tail
    idx = np.where(valid & (gap <= 0.0))[0]
    if idx.size:
        i_match = int(idx[0])
    elif np.any(valid):
        # Never quite reached the asymptote: match at the steepest-decay point.
        vi = np.where(valid)[0]
        i_match = int(vi[np.argmin(gap[vi])])
    else:
        i_match = max(1, int(np.argmin(phi)))

    r_m = float(r_core[i_match])
    phi_m = float(phi[i_match])
    r_core_keep = r_core[: i_match + 1]
    phi_core_keep = phi[: i_match + 1].copy()
    phi_core_keep[-1] = phi_m

    # Analytic Yukawa tail beyond the match radius, value-matched at r_m.
    tail_end = max(r_max, r_m + 20.0 / max(kappa, 1.0e-6))
    r_tail = np.linspace(r_m, tail_end, max(200, int((tail_end - r_m) * 20)))[1:]
    phi_tail = phi_m * (r_m / r_tail) * np.exp(-kappa * (r_tail - r_m))

    r_full = np.concatenate([r_core_keep, r_tail])
    phi_full = np.concatenate([phi_core_keep, phi_tail])
    return r_full, phi_full


def solve_qball_radial_profile(
    couplings: QBallCouplings,
    *,
    r_max: float | None = None,
) -> QBallRadialProfile:
    """Find the flat-space Q-ball profile for *couplings* via outward shooting.

    Uses event-based undershoot/overshoot classification, bisects φ(0) to ~machine
    precision, then truncates the near-separatrix solution at its asymptotic-match
    radius and attaches the analytic e^{−κ r}/r tail so the tabulated profile is a
    genuinely localized soliton (decays to ~0 within ``r_max``).
    """
    couplings.validate()
    if couplings.lam <= 0.0 or couplings.mu <= 0.0:
        raise ValueError("Q-ball ODE requires lam > 0 and mu > 0")

    mass, lam, mu, omega = couplings.mass, couplings.lam, couplings.mu, couplings.omega
    kappa_sq = mass * mass - omega * omega
    if kappa_sq <= 0.0:
        raise ValueError("omega must be < mass for a bound Q-ball profile")
    kappa = math.sqrt(kappa_sq)

    width = couplings.bound_state_width
    r_max = float(r_max if r_max is not None else max(60.0, 24.0 * width))

    # The soliton starts decreasing (φ''(0) < 0), which requires φ_c in the window
    # where λφ_c² − μφ_c⁴ > κ², i.e. φ_c² ∈ ((λ−√Δ)/2μ, (λ+√Δ)/2μ) with
    # Δ = λ² − 4μκ².  Outside it the linear/sextic terms make φ grow immediately.
    disc = lam * lam - 4.0 * mu * kappa_sq
    if disc <= 0.0:
        raise RuntimeError(
            "Q-ball ODE shoot: no φ''(0)<0 window (λ²<4μκ²); couplings admit no soliton"
        )
    sqrt_disc = math.sqrt(disc)
    phi_win_lo = math.sqrt((lam - sqrt_disc) / (2.0 * mu))
    phi_win_hi = math.sqrt((lam + sqrt_disc) / (2.0 * mu))

    def kind(phi_c: float) -> str:
        return _shoot(phi_c, kappa_sq=kappa_sq, lam=lam, mu=mu, r_max=r_max)[0]

    # The separatrix φ_c lies inside the φ''(0)<0 window (φ_win_lo, φ_win_hi):
    # below it the field turns back to the false vacuum ("under"); at/above
    # φ_win_hi it blows up ("over").  For thin-wall couplings φ_c sits a hair below
    # φ_win_hi, so bracket [under, just-below-φ_win_hi] directly rather than scan.
    phi_lo = phi_win_lo
    if kind(phi_lo) != "under":
        # Walk down until we find the false-vacuum (under) side.
        for _ in range(60):
            phi_lo *= 0.9
            if kind(phi_lo) == "under":
                break
        else:
            raise RuntimeError("Q-ball ODE shoot: could not bracket the under side")

    # At/above φ_win_hi the field blows up to +∞ (the "over" side).  Use a small
    # nudge past the window edge as the upper bracket.
    phi_hi = phi_win_hi * (1.0 + 1.0e-6)
    if kind(phi_hi) != "over":
        raise RuntimeError("Q-ball ODE shoot: could not bracket the over side")

    # Bisect to ~machine precision on φ_c (the knife-edge separatrix).  Keep the
    # tightest *under* trajectory: it tracks the true soliton down to its minimum.
    best_phi_c = phi_lo
    best_sol = None
    for _ in range(120):
        mid = 0.5 * (phi_lo + phi_hi)
        k, sol = _shoot(mid, kappa_sq=kappa_sq, lam=lam, mu=mu, r_max=r_max)
        if k == "under":
            phi_lo = mid
            best_phi_c = mid
            best_sol = sol
        else:
            phi_hi = mid
        if phi_hi - phi_lo <= 1.0e-15 * max(1.0, phi_hi):
            break

    if best_sol is None:
        # Degenerate: take the tightest undershoot we have.
        _, best_sol = _shoot(phi_lo, kappa_sq=kappa_sq, lam=lam, mu=mu, r_max=r_max)
        best_phi_c = phi_lo

    r_full, phi_full = _build_localized_profile(
        best_phi_c, best_sol, kappa=kappa, r_max=r_max
    )

    spline = CubicSpline(r_full, phi_full, bc_type="natural")
    return QBallRadialProfile(
        mass=mass,
        lam=lam,
        mu=mu,
        omega=omega,
        phi_c=best_phi_c,
        r=r_full,
        phi0=phi_full,
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
