"""Tests for the self-gravitating boson star ODE solver.

Validates:
  * Localization (the 2026-06-30 non-localization regression lesson).
  * Eigenvalue omega is in the bound-state range (0, m).
  * Zero nodes in the ground state.
  * ADM mass reproduces the Kaup limit M_max * m ~= 0.633 for the mini star.
  * Monotone decay (no unphysical turn-up).
  * Sextic family produces localized profiles too.
"""

from __future__ import annotations


import numpy as np
import pytest

from grteclyn_wrapper.grtresna.profiles.boson_star_ode import (
    BosonStarODEProfile,
    cached_selfgrav_profile,
    solve_selfgrav_boson_star,
)
from grteclyn_wrapper.grtresna.profiles.qball_couplings import QBallCouplings


@pytest.fixture(autouse=True)
def _clear_cache() -> None:
    cached_selfgrav_profile.cache_clear()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _radius_at_fraction(prof: BosonStarODEProfile, frac: float) -> float:
    """Smallest radius where phi0 has fallen below *frac* * phi0(0)."""
    core = float(prof.phi0[0])
    below = np.where(prof.phi0 / core < frac)[0]
    return float(prof.r[below[0]]) if below.size else float("inf")


def _count_nodes(prof: BosonStarODEProfile) -> int:
    """Count zero crossings of phi0 (excluding the tail where it's ~0)."""
    core = float(prof.phi0[0])
    threshold = 0.01 * core
    significant = prof.phi0[np.abs(prof.phi0) > threshold]
    if significant.size < 2:
        return 0
    sign_changes = np.diff(np.sign(significant))
    return int(np.sum(sign_changes != 0))


# ---------------------------------------------------------------------------
# Mini boson star (lambda = mu = 0)
# ---------------------------------------------------------------------------

def test_mini_star_eigenvalue_in_bound_range() -> None:
    """omega must lie in (0, m) for a bound state."""
    c = QBallCouplings.free_field(mass=1.0, omega=0.5)
    prof = solve_selfgrav_boson_star(c, phi_c=0.1)
    assert 0.0 < prof.omega < prof.mass


def test_mini_star_localized_and_monotone() -> None:
    """The profile must be a genuine localized soliton, not a flat-top slab.

    Regression guard for the 2026-06-30 bug where the flat-space Q-ball solver
    returned a near-constant box-spanning condensate.
    """
    c = QBallCouplings.free_field(mass=1.0, omega=0.5)
    prof = solve_selfgrav_boson_star(c, phi_c=0.1)
    core = float(prof.phi0[0])
    # Decays to <1% of core within the tabulated range.
    assert _radius_at_fraction(prof, 0.01) < float(prof.r[-1])
    # Actually reaches ~zero.
    assert prof.phi0.min() / core < 1.0e-2
    # Monotonically non-increasing (no unphysical turn-up).
    assert np.all(np.diff(prof.phi0) <= 1.0e-8)


def test_mini_star_zero_nodes() -> None:
    """Ground state must have no nodes (no sign changes in phi0)."""
    c = QBallCouplings.free_field(mass=1.0, omega=0.5)
    prof = solve_selfgrav_boson_star(c, phi_c=0.1)
    assert _count_nodes(prof) == 0


def test_mini_star_adm_mass_positive() -> None:
    """A self-gravitating star must have positive ADM mass."""
    c = QBallCouplings.free_field(mass=1.0, omega=0.5)
    prof = solve_selfgrav_boson_star(c, phi_c=0.1)
    assert prof.adm_mass > 0.0


def test_mini_star_kaup_limit() -> None:
    """The maximum ADM mass over the stable branch must reproduce M*m ~ 0.633.

    The Kaup limit (Kaup 1968, Ruffini & Bonazzola 1969) for a mini boson star
    (lambda=mu=0) is M_max * m ~= 0.633 in G=1 units.  We scan central
    amplitudes and check the peak mass is in the right ballpark.
    """
    mass = 1.0
    phi_c_values = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20]
    masses = []
    for phi_c in phi_c_values:
        c = QBallCouplings.free_field(mass=mass, omega=0.5 * mass)
        try:
            prof = solve_selfgrav_boson_star(c, phi_c=phi_c)
            if prof.adm_mass > 0.0 and prof.phi0.min() / prof.phi0[0] < 1e-2:
                masses.append(prof.adm_mass)
        except (RuntimeError, ValueError):
            pass  # phi_c too large -> no stable solution (past Kaup max)
    assert len(masses) >= 3, "expected at least 3 stable configurations in the scan"
    m_max = max(masses)
    # Literature: M_max * m ~ 0.633 (Kaup 1968).  The isotropic solver reproduces
    # this to ~2%; allow abs=0.1 for the coarse central-amplitude scan (the true
    # peak sits between the sampled phi_c values) and finite r_max truncation.
    assert m_max * mass == pytest.approx(0.633, abs=0.1)


def test_mini_star_metric_approaches_schwarzschild() -> None:
    """The conformal factor psi(r) and lapse alpha(r) must approach 1.

    Asymptotic flatness in isotropic coordinates: psi -> 1 and alpha -> 1 at
    large r (the two exact scaling rescalings pin both to unity).
    """
    c = QBallCouplings.free_field(mass=1.0, omega=0.5)
    prof = solve_selfgrav_boson_star(c, phi_c=0.1)
    assert abs(float(prof.psi[-1]) - 1.0) < 0.05
    assert abs(float(prof.alpha[-1]) - 1.0) < 0.05
    # The conformal factor is >= 1 everywhere (positive ADM mass) and peaks at
    # the centre.
    assert float(prof.psi[0]) > 1.0
    assert float(prof.psi.max()) == float(prof.psi[0])


def test_mini_star_isotropic_exterior_matches_schwarzschild() -> None:
    """Outside the star the metric must match the isotropic Schwarzschild form.

    In isotropic coordinates the vacuum exterior is
    ``psi = 1 + M/(2r)`` and ``alpha = (1 - M/2r)/(1 + M/2r)``.  We check both
    against the ADM mass read off the profile at a radius well outside the core.
    """
    c = QBallCouplings.free_field(mass=1.0, omega=0.5)
    prof = solve_selfgrav_boson_star(c, phi_c=0.08)
    m = prof.adm_mass
    # Sample at a radius where the field has decayed to ~0 (deep exterior).
    core = float(prof.phi0[0])
    ext = np.where((prof.phi0 / core < 1e-4) & (prof.r > 0.0))[0]
    assert ext.size > 0
    i = int(ext[len(ext) // 2])
    r = float(prof.r[i])
    half = 0.5 * m / r
    assert float(prof.psi[i]) == pytest.approx(1.0 + half, rel=2e-2)
    assert float(prof.alpha[i]) == pytest.approx((1.0 - half) / (1.0 + half), rel=2e-2)


def test_mini_star_areal_radius_is_monotone() -> None:
    """The areal radius R = psi^2 * r_iso must increase monotonically.

    A valid isotropic->areal coordinate map is strictly monotone; a non-monotone
    R would signal a corrupt conformal factor.
    """
    c = QBallCouplings.free_field(mass=1.0, omega=0.5)
    prof = solve_selfgrav_boson_star(c, phi_c=0.1)
    areal = prof.psi ** 2 * prof.r
    assert np.all(np.diff(areal) > -1e-9)


# ---------------------------------------------------------------------------
# Sextic family (lambda > 0, mu > 0)
# ---------------------------------------------------------------------------

def test_sextic_star_localized() -> None:
    """A sextic self-gravitating star must also be localized."""
    c = QBallCouplings(mass=1.0, lam=10.0, mu=100.0, omega=0.5, phi_core=0.0)
    # Use a modest central amplitude on the stable branch.
    prof = solve_selfgrav_boson_star(c, phi_c=0.1)
    core = float(prof.phi0[0])
    assert _radius_at_fraction(prof, 0.01) < float(prof.r[-1])
    assert prof.phi0.min() / core < 1.0e-2
    assert _count_nodes(prof) == 0
    assert 0.0 < prof.omega < prof.mass


def test_sextic_star_adm_mass_positive() -> None:
    c = QBallCouplings(mass=1.0, lam=10.0, mu=100.0, omega=0.5, phi_core=0.0)
    prof = solve_selfgrav_boson_star(c, phi_c=0.1)
    assert prof.adm_mass > 0.0


# ---------------------------------------------------------------------------
# Cache + lump helper
# ---------------------------------------------------------------------------

def test_cached_profile_returns_same_result() -> None:
    c = QBallCouplings.free_field(mass=1.0, omega=0.5)
    p1 = solve_selfgrav_boson_star(c, phi_c=0.08)
    p2 = cached_selfgrav_profile(1.0, 0.0, 0.0, 0.08)
    assert p1.omega == pytest.approx(p2.omega, rel=1e-6)
    assert p1.adm_mass == pytest.approx(p2.adm_mass, rel=1e-6)


def test_profile_for_lump_selfgrav_default_amplitude() -> None:
    from grteclyn_wrapper.grtresna.profiles.boson_star_ode import (
        profile_for_lump_selfgrav,
    )
    lump: dict = {}
    prof = profile_for_lump_selfgrav(lump, mass=1.0, lam=0.0, mu=0.0)
    assert prof.phi_c > 0.0
    assert 0.0 < prof.omega < 1.0


# ---------------------------------------------------------------------------
# Pump-target sech width fit
# ---------------------------------------------------------------------------

def test_sech_width_matches_half_maximum() -> None:
    """The fitted sech width must reproduce the profile's half-maximum radius.

    A sech ``1/cosh(r/w)`` has half-maximum at ``r = w * arccosh(2)``, so the
    fitted width must satisfy ``phi0(w * arccosh(2)) ~= phi0(0)/2``.
    """
    from grteclyn_wrapper.grtresna.profiles.boson_star_ode import (
        sech_width_for_profile,
    )
    c = QBallCouplings.free_field(mass=1.0, omega=0.5)
    prof = solve_selfgrav_boson_star(c, phi_c=0.1)
    w = sech_width_for_profile(prof)
    assert w > 0.0
    r_half = w * np.arccosh(2.0)
    assert float(prof.eval_phi0(r_half)) == pytest.approx(0.5 * float(prof.phi0[0]), rel=5e-2)


def test_selfgrav_pump_target_width_more_compact_than_tail_length() -> None:
    """The star's core is tighter than a tail-matched sech.

    The tail-only length ``1/sqrt(m^2 - omega^2)`` overestimates the core width;
    the profile-fitted width should be strictly smaller for a typical star.
    """
    from grteclyn_wrapper.grtresna.profiles.boson_star import bound_width
    from grteclyn_wrapper.grtresna.profiles.boson_star_ode import (
        selfgrav_pump_target_width,
    )
    prof = cached_selfgrav_profile(1.0, 0.0, 0.0, 0.1)
    fitted = selfgrav_pump_target_width(1.0, 0.0, 0.0, 0.1)
    tail = bound_width(1.0, prof.omega)
    assert 0.0 < fitted < tail
