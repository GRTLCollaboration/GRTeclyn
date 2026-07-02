"""Tests for complex scalar EM tensor reference physics."""

from __future__ import annotations

import numpy as np
import pytest

from grteclyn_wrapper.grtresna.profiles.complex_scalar import (
    ComplexScalarState,
    compute_emtensor,
    conserved_charge_density,
    flat_minkowski_h,
    phase_velocity_pi2,
    potential_and_derivatives,
    random_test_state,
)


def test_potential_uses_modulus_squared() -> None:
    v1, _, _ = potential_and_derivatives(0.3, 0.4, mass=0.5, lam=0.0)
    v2, _, _ = potential_and_derivatives(0.5, 0.0, mass=0.5, lam=0.0)
    assert abs(v1 - v2) < 1.0e-12


def test_potential_derivatives_chain_rule() -> None:
    phi1, phi2, mass, lam = 0.2, 0.1, 0.6, 0.05
    eps = 1.0e-7
    v0, dv1, dv2 = potential_and_derivatives(phi1, phi2, mass, lam)
    vp1, _, _ = potential_and_derivatives(phi1 + eps, phi2, mass, lam)
    vp2, _, _ = potential_and_derivatives(phi1, phi2 + eps, mass, lam)
    assert abs((vp1 - v0) / eps - dv1) < 1.0e-4
    assert abs((vp2 - v0) / eps - dv2) < 1.0e-4


def test_emtensor_includes_kinetic_term(rng: np.random.Generator) -> None:
    """rho must carry 1/2 (Pi_re^2 + Pi_im^2): the boson-star phase energy."""
    state = random_test_state(rng)
    static = ComplexScalarState(
        phi1=state.phi1, phi2=state.phi2, pi1=0.0, pi2=0.0,
        chi=state.chi, h_uu=state.h_uu, h_ll=state.h_ll,
        dphi1=state.dphi1, dphi2=state.dphi2,
    )
    em = compute_emtensor(state, mass=0.1, lam=0.0)
    em_static = compute_emtensor(static, mass=0.1, lam=0.0)
    pi_sq = state.pi1 ** 2 + state.pi2 ** 2
    assert abs((em.rho - em_static.rho) - 0.5 * pi_sq) < 1.0e-12


def test_emtensor_real_scalar_limit(rng: np.random.Generator) -> None:
    h_uu, h_ll = flat_minkowski_h()
    phi1, pi1 = 0.25, -0.1
    dphi1 = rng.uniform(-0.1, 0.1, size=3)
    state = ComplexScalarState(
        phi1=phi1, phi2=0.0, pi1=pi1, pi2=0.0,
        chi=1.0, h_uu=h_uu, h_ll=h_ll,
        dphi1=dphi1, dphi2=np.zeros(3),
    )
    em = compute_emtensor(state, mass=0.2, lam=0.0)
    grad_sq = float(dphi1 @ dphi1)
    v = 0.5 * (0.2 * phi1) ** 2
    expected_rho = 0.5 * pi1 ** 2 + 0.5 * grad_sq + v
    assert abs(em.rho - expected_rho) < 1.0e-12


def test_emtensor_rho_closed_form(rng: np.random.Generator) -> None:
    """rho == 1/2 sum Pi^2 + 1/2 chi sum|grad phi|^2 + V(|Phi|^2)."""
    state = random_test_state(rng)
    mass, lam = 0.3, 0.02
    em = compute_emtensor(state, mass=mass, lam=lam)
    grad1 = float(state.dphi1 @ (state.chi * (state.h_uu @ state.dphi1)))
    grad2 = float(state.dphi2 @ (state.chi * (state.h_uu @ state.dphi2)))
    v, _, _ = potential_and_derivatives(state.phi1, state.phi2, mass, lam)
    expected = (
        0.5 * (state.pi1 ** 2 + state.pi2 ** 2)
        + 0.5 * (grad1 + grad2)
        + v
    )
    assert abs(em.rho - expected) < 1.0e-12


def test_conserved_charge_nonzero_with_phase_velocity() -> None:
    phi0, omega, alpha = 0.1, 0.08, 0.9
    pi2 = phase_velocity_pi2(phi0, omega, alpha)
    q = conserved_charge_density(phi0, 0.0, 0.0, pi2)
    assert q < 0.0  # wormhole convention: Pi2 = -omega phi0 / alpha
    assert abs(q) > 0.0


def test_conserved_charge_zero_without_phase_velocity() -> None:
    q = conserved_charge_density(0.1, 0.0, 0.0, 0.0)
    assert q == 0.0


# --- Phantom sign tests ---------------------------------------------------


def test_phantom_rho_negative(rng: np.random.Generator) -> None:
    """sign=-1 must flip rho: phantom boson star has negative energy density."""
    state = random_test_state(rng)
    mass = 0.2
    em_canonical = compute_emtensor(state, mass=mass, lam=0.0, sign=1.0)
    em_phantom = compute_emtensor(state, mass=mass, lam=0.0, sign=-1.0)
    assert abs(em_phantom.rho + em_canonical.rho) < 1.0e-12


def test_phantom_j_flipped(rng: np.random.Generator) -> None:
    """Momentum density reverses sign under phantom flip."""
    state = random_test_state(rng)
    em_can = compute_emtensor(state, mass=0.1, lam=0.0, sign=1.0)
    em_ph = compute_emtensor(state, mass=0.1, lam=0.0, sign=-1.0)
    np.testing.assert_allclose(em_ph.j, -em_can.j, atol=1.0e-12)


def test_phantom_S_flipped(rng: np.random.Generator) -> None:
    """Spatial stress tensor reverses sign under phantom flip."""
    state = random_test_state(rng)
    em_can = compute_emtensor(state, mass=0.1, lam=0.0, sign=1.0)
    em_ph = compute_emtensor(state, mass=0.1, lam=0.0, sign=-1.0)
    np.testing.assert_allclose(em_ph.S, -em_can.S, atol=1.0e-12)


def test_phantom_trS_flipped(rng: np.random.Generator) -> None:
    """trS reverses sign under phantom flip."""
    state = random_test_state(rng)
    em_can = compute_emtensor(state, mass=0.15, lam=0.01, sign=1.0)
    em_ph = compute_emtensor(state, mass=0.15, lam=0.01, sign=-1.0)
    assert abs(em_ph.trS + em_can.trS) < 1.0e-12


def test_phantom_sign_default_is_canonical(rng: np.random.Generator) -> None:
    """Omitting sign (default=1.0) must match explicit sign=1.0."""
    state = random_test_state(rng)
    em_default = compute_emtensor(state, mass=0.1, lam=0.0)
    em_explicit = compute_emtensor(state, mass=0.1, lam=0.0, sign=1.0)
    assert abs(em_default.rho - em_explicit.rho) < 1.0e-15
