"""Tests for boson_star_profile radial ODE solver."""

from __future__ import annotations

import numpy as np

from grteclyn_wrapper.grtresna.profiles.boson_star import solve_mini_boson_star


def test_profile_decays_at_infinity() -> None:
    profile = solve_mini_boson_star(mass=0.1, phi_c=0.05, r_max=60.0, n_points=200)
    assert profile.omega > 0.0
    assert profile.phi0[0] > 0.0
    assert profile.phi0[-1] / profile.phi0[0] < 0.5
    assert profile.adm_mass > 0.0


def test_spline_derivatives_continuous() -> None:
    profile = solve_mini_boson_star(mass=0.1, phi_c=0.05, r_max=60.0, n_points=200)
    r_test = profile.r[10:-10]
    d1 = profile.spline.derivative(1)(r_test)
    d2 = profile.spline.derivative(2)(r_test)
    assert np.all(np.isfinite(d1))
    assert np.all(np.isfinite(d2))
    # No large jumps between adjacent knot intervals (C2 spline).
    d1_left = profile.spline.derivative(1)(profile.r[1:-1])
    d1_right = profile.spline.derivative(1)(profile.r[1:-1] + 1.0e-8)
    assert np.max(np.abs(d1_left - d1_right)) < 1.0e-3


def test_eval_phi0_at_center() -> None:
    profile = solve_mini_boson_star(mass=0.1, phi_c=0.05, r_max=60.0, n_points=200)
    assert abs(profile.eval_phi0(profile.r[0]) - profile.phi0[0]) < 1.0e-10
