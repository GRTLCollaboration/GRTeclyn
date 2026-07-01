"""Tests for flat-space Q-ball radial ODE solver."""

from __future__ import annotations

import numpy as np
import pytest

from grteclyn_wrapper.grtresna.qball_couplings import QBallCouplings
from grteclyn_wrapper.grtresna.qball_radial_ode import (
    cached_qball_radial_profile,
    solve_qball_radial_profile,
)


@pytest.fixture(autouse=True)
def _clear_profile_cache() -> None:
    cached_qball_radial_profile.cache_clear()


def test_standard_profile_peak_near_core_amplitude() -> None:
    c = QBallCouplings.standard()
    prof = solve_qball_radial_profile(c)
    assert prof.phi_c == pytest.approx(c.core_amplitude, rel=0.05)
    assert float(prof.eval_phi0(0.0)) == pytest.approx(prof.phi_c, rel=1e-3)
    assert float(prof.eval_phi0(float(prof.r[-1]))) < prof.phi_c


def test_stiff_profile_peak_is_smaller_than_standard() -> None:
    std = solve_qball_radial_profile(QBallCouplings.standard())
    stiff = solve_qball_radial_profile(QBallCouplings.stiff())
    assert stiff.phi_c == pytest.approx(QBallCouplings.stiff().core_amplitude, rel=0.05)
    assert stiff.phi_c < std.phi_c


def test_profile_scales_with_requested_lump_amp() -> None:
    from grteclyn_wrapper.grtresna.qball_radial_ode import profile_for_lump

    c = QBallCouplings.standard()
    base = solve_qball_radial_profile(c)
    lump = {
        "amp": 0.10,
        "profile": 3,
        "qball_mass": c.mass,
        "qball_lam": c.lam,
        "qball_mu": c.mu,
        "qball_omega": c.omega,
    }
    scaled = profile_for_lump(lump, mass=c.mass, lam=c.lam, mu=c.mu, omega=c.omega)
    assert float(scaled.eval_phi0(0.0)) == pytest.approx(0.10, rel=1e-3)
    assert float(scaled.eval_phi0(1.0)) == pytest.approx(
        0.10 / base.phi_c * float(base.eval_phi0(1.0)), rel=0.02
    )


def _radius_at_fraction(prof, frac: float) -> float:
    """Smallest radius where φ₀ has fallen below *frac*·φ₀(0)."""
    core = float(prof.phi0[0])
    below = np.where(prof.phi0 / core < frac)[0]
    return float(prof.r[below[0]]) if below.size else float("inf")


@pytest.mark.parametrize(
    "couplings",
    [
        QBallCouplings.standard(),
        QBallCouplings.stiff(),
        QBallCouplings.compact(),
    ],
)
def test_profile_is_localized_and_monotone(couplings) -> None:
    """The shoot must return a genuine bound soliton, not a flat-top condensate.

    Regression guard for the bug where the loose bisection returned a profile that
    plateaued near ~0.5·core out to r=100 (see NextSteps.md §3).
    """
    prof = solve_qball_radial_profile(couplings)
    core = float(prof.phi0[0])
    # Decays to <1% of core within the tabulated range.
    assert _radius_at_fraction(prof, 0.01) < float(prof.r[-1])
    # Actually reaches ~zero (no residual false-vacuum plateau).
    assert prof.phi0.min() / core < 1.0e-3
    # Monotonically non-increasing (no unphysical turn-up tail).
    assert np.all(np.diff(prof.phi0) <= 1.0e-9)


def test_compact_omega_gives_smaller_soliton_than_thinwall() -> None:
    """compact() (ω=0.8) must be box-fitting vs the thin-wall stiff() (ω=0.4)."""
    thin = solve_qball_radial_profile(QBallCouplings.stiff())
    comp = solve_qball_radial_profile(QBallCouplings.compact())
    assert _radius_at_fraction(comp, 0.01) < _radius_at_fraction(thin, 0.01)
    # Compact soliton should decay within roughly the lump-scale, not the whole box.
    assert _radius_at_fraction(comp, 0.01) < 12.0


def test_ode_profile_differs_from_sech() -> None:
    from grteclyn_wrapper.grtresna.boson_star_profile import bound_width

    c = QBallCouplings.standard()
    prof = solve_qball_radial_profile(c)
    r = np.linspace(0.5, 5.0, 20)
    width = bound_width(c.mass, c.omega)
    sech = prof.phi_c / np.cosh(r / width)
    ode = np.asarray(prof.eval_phi0(r), dtype=float)
    assert not np.allclose(ode, sech, rtol=0.05)
