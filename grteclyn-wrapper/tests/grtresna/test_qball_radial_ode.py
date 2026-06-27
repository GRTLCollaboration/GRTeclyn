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


def test_ode_profile_differs_from_sech() -> None:
    from grteclyn_wrapper.grtresna.boson_star_profile import bound_width

    c = QBallCouplings.standard()
    prof = solve_qball_radial_profile(c)
    r = np.linspace(0.5, 5.0, 20)
    width = bound_width(c.mass, c.omega)
    sech = prof.phi_c / np.cosh(r / width)
    ode = np.asarray(prof.eval_phi0(r), dtype=float)
    assert not np.allclose(ode, sech, rtol=0.05)
