"""Tests for QBallCouplings presets and derived physics."""

from __future__ import annotations

import math

import pytest

from grteclyn_wrapper.grtresna.profiles.qball_couplings import QBallCouplings


def test_standard_and_stiff_share_omega_min() -> None:
    std = QBallCouplings.standard()
    stiff = QBallCouplings.stiff()
    assert std.omega_min == pytest.approx(stiff.omega_min, rel=1e-3)
    assert std.omega_min == pytest.approx(0.316, rel=0.01)


def test_stiff_preserves_omega_min_and_scales_lambda() -> None:
    std = QBallCouplings.standard()
    stiff = QBallCouplings.stiff()
    assert stiff.lam == pytest.approx(4.0 * std.lam)
    assert stiff.omega_min == pytest.approx(std.omega_min, rel=1e-3)
    assert stiff.lam == pytest.approx(640.0)
    assert stiff.mu == pytest.approx(85333.0, rel=1e-3)


def test_stiff_equilibrium_core_is_half_standard() -> None:
    std = QBallCouplings.standard()
    stiff = QBallCouplings.stiff()
    assert std.core_amplitude == pytest.approx(0.15, rel=0.01)
    assert stiff.core_amplitude == pytest.approx(0.075, rel=0.01)
    assert std.scaled_preserving_band(4.0).mu == pytest.approx(85328.0, rel=1e-3)


def test_scaled_preserves_core_amplitude() -> None:
    std = QBallCouplings.standard()
    assert std.scaled(4.0).core_amplitude == pytest.approx(std.core_amplitude)


def test_to_search_overrides_keys() -> None:
    overrides = QBallCouplings.stiff().to_search_overrides()
    assert overrides == {
        "grtresna_scalar_mass": 1.0,
        "grtresna_scalar_lambda": 640.0,
        "grtresna_scalar_mu": 85333.0,
        "grtresna_bs_omega": 0.4,
    }


def test_validate_rejects_unbound_omega() -> None:
    bad = QBallCouplings(mass=1.0, lam=160.0, mu=5333.0, omega=0.2)
    with pytest.raises(ValueError, match="omega_min"):
        bad.validate()


def test_free_field_skips_qball_band_check() -> None:
    free = QBallCouplings.free_field(mass=1.0, omega=0.4)
    free.validate()
    assert free.lam == 0.0
    assert free.mu == 0.0


def test_bound_state_width_matches_boson_profile() -> None:
    c = QBallCouplings.standard()
    expected = 1.0 / math.sqrt(c.mass * c.mass - c.omega * c.omega)
    assert c.bound_state_width == pytest.approx(expected)


def test_merge_into_overrides() -> None:
    base = {"trajectory_mode": 1}
    merged = QBallCouplings.stiff().merge_into(base)
    assert merged["grtresna_scalar_mu"] == 85333.0
    assert merged["trajectory_mode"] == 1


def test_cap_well_depth_at_equilibrium_core() -> None:
    stiff = QBallCouplings.stiff()
    assert stiff.cap_well_depth(0.15) == pytest.approx(0.075, rel=0.01)
    assert stiff.cap_well_depth(0.05) == pytest.approx(0.05)


def test_seed_overrides_include_dispersion_flags() -> None:
    keys = QBallCouplings.stiff().seed_overrides(
        equilibrium_amplitude=True,
        ode_profile=True,
    )
    assert keys["grtresna_qball_equilibrium_amplitude"] == 1
    assert keys["grtresna_qball_ode_profile"] == 1
    assert keys["grtresna_scalar_lambda"] == 640.0


def test_seed_overrides_selfgrav_flag() -> None:
    keys = QBallCouplings.mini_selfgrav(mass=1.0).seed_overrides(selfgrav=True)
    assert keys["grtresna_bs_selfgrav"] == 1
    # Default (no selfgrav) must not emit the flag.
    assert "grtresna_bs_selfgrav" not in QBallCouplings.stiff().seed_overrides()


def test_with_equilibrium_paint_sets_phi_core() -> None:
    stiff = QBallCouplings.stiff().with_equilibrium_paint()
    assert stiff.phi_core == pytest.approx(stiff.core_amplitude)
