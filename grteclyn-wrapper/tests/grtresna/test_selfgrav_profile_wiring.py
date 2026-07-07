"""Tests for wiring the self-gravitating boson star through GRTresna params."""

from __future__ import annotations

from pathlib import Path

import pytest

from grteclyn_wrapper.grtresna.profiles.boson_star import (
    PROFILE_ODE_BOUND,
    PROFILE_SELFGRAV_BOUND,
    grtresna_lump_profile,
)
from grteclyn_wrapper.grtresna.profiles.boson_star_ode import cached_selfgrav_profile
from grteclyn_wrapper.grtresna.solver.config import GRTresnaConfig
from grteclyn_wrapper.grtresna.solver.params import write_grtresna_params


@pytest.fixture(autouse=True)
def _clear_cache() -> None:
    cached_selfgrav_profile.cache_clear()


def test_selfgrav_profile_maps_to_c_profile_3() -> None:
    """The Python selfgrav tag emits as C++ profile 3 (shared tabulated loader)."""
    assert grtresna_lump_profile(PROFILE_SELFGRAV_BOUND) == PROFILE_ODE_BOUND
    assert grtresna_lump_profile(PROFILE_ODE_BOUND) == PROFILE_ODE_BOUND


def _selfgrav_cfg() -> GRTresnaConfig:
    return GRTresnaConfig(
        matter_model="grtresna_complex_scalar",
        scalar_mass=1.0,
        scalar_lambda=0.0,
        scalar_mu=0.0,
        bs_phi_c=0.1,
        bs_omega=0.0,
        lumps=[
            {
                "amp": 0.1,
                "width": 5.0,
                "center": (0.0, 0.0, 0.0),
                "profile": PROFILE_SELFGRAV_BOUND,
            }
        ],
    )


def test_selfgrav_params_write_sets_eigenvalue(tmp_path: Path) -> None:
    """Writing params solves the eigenvalue and back-writes bs_omega + lump omega."""
    cfg = _selfgrav_cfg()
    params_path = tmp_path / "params.txt"
    write_grtresna_params(cfg, params_path)

    # Eigenvalue written back onto the config and the lump.
    assert 0.0 < cfg.bs_omega < cfg.scalar_mass
    assert cfg.lumps[0]["omega"] == pytest.approx(cfg.bs_omega)

    # A profile table was emitted and referenced from params.
    dat = tmp_path / "qball_profile.dat"
    assert dat.exists()
    text = params_path.read_text()
    assert "qball_profile_path" in text
    assert f"bs_omega = {cfg.bs_omega}" in text
    # profile is emitted as C++ profile 3.
    assert f"lump0_profile = {PROFILE_ODE_BOUND}" in text


def test_selfgrav_profile_table_is_localized(tmp_path: Path) -> None:
    """The emitted table must be a localized soliton (regression guard)."""
    cfg = _selfgrav_cfg()
    write_grtresna_params(cfg, tmp_path / "params.txt")
    dat = tmp_path / "qball_profile.dat"
    rows = [
        line.split()
        for line in dat.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    phi0 = [float(r[1]) for r in rows]
    core = phi0[0]
    assert core > 0.0
    # Decays to <1% of core by the end of the table.
    assert abs(phi0[-1]) / core < 1.0e-2


def test_selfgrav_profile_table_has_lapse_column(tmp_path: Path) -> None:
    """The self-grav table must emit a 3rd column alpha(r) for the stationary
    momentum Pi_im = -(omega/alpha) phi0.

    Regression guard for the two dominant dispersal bugs:
      * the GRTresna loader must be able to parse this file (comment header +
        3 columns), and
      * alpha must dip below 1 in the core (a boson star is stationary, not
        static) and rise back toward 1 in the tail.
    """
    cfg = _selfgrav_cfg()
    write_grtresna_params(cfg, tmp_path / "params.txt")
    dat = tmp_path / "qball_profile.dat"
    rows = [
        line.split()
        for line in dat.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    assert all(len(r) == 3 for r in rows), "self-grav table must have r phi0 alpha"
    alpha = [float(r[2]) for r in rows]
    assert 0.0 < alpha[0] < 1.0            # central lapse dips below 1
    assert alpha[-1] == pytest.approx(1.0, abs=0.05)  # asymptotically flat
    assert alpha[-1] > alpha[0]            # rises outward


def test_selfgrav_pump_target_uses_fitted_width() -> None:
    """Evolution overrides drive the pump toward the self-grav star's shape.

    The trap target is a sech (profile 2) whose width is fitted to the ODE star
    (more compact than the tail-only bound length), not the broad tail-matched
    sech used for plain Q-ball seeds.
    """
    from grteclyn_wrapper.grtresna.matter.wiring import (
        GRTresnaMatterMetadata,
        evolution_overrides_from_config,
    )
    from grteclyn_wrapper.grtresna.profiles.boson_star import PROFILE_SECH_BOUND, bound_width
    from grteclyn_wrapper.grtresna.profiles.boson_star_ode import (
        selfgrav_pump_target_width,
    )

    cfg = _selfgrav_cfg()
    cfg.bs_omega = 0.85  # populated by the eigenvalue solve in the real pipeline
    ov = evolution_overrides_from_config(cfg)
    assert ov["rl_pump_target_profile"] == PROFILE_SECH_BOUND
    expected = selfgrav_pump_target_width(1.0, 0.0, 0.0, 0.1)
    assert ov["rl_pump_target_width"] == pytest.approx(expected)
    # Tighter than the tail-only sech a plain Q-ball would use.
    assert ov["rl_pump_target_width"] < bound_width(1.0, cfg.bs_omega)

    # The self-grav marker survives metadata (de)serialisation.
    meta = GRTresnaMatterMetadata.from_config(cfg)
    assert meta.bs_selfgrav is True
    assert meta.to_evolution_overrides()["rl_pump_target_width"] == pytest.approx(expected)


def test_ode_qball_path_unchanged(tmp_path: Path) -> None:
    """A plain profile==3 Q-ball still emits its flat-space table (no regression)."""
    cfg = GRTresnaConfig(
        matter_model="grtresna_complex_scalar",
        scalar_mass=1.0,
        scalar_lambda=640.0,
        scalar_mu=85333.0,
        bs_omega=0.8,
        lumps=[
            {
                "amp": 0.075,
                "width": 1.667,
                "center": (0.0, 0.0, 0.0),
                "profile": PROFILE_ODE_BOUND,
            }
        ],
    )
    write_grtresna_params(cfg, tmp_path / "params.txt")
    dat = tmp_path / "qball_profile.dat"
    assert dat.exists()
    header = dat.read_text().splitlines()[0]
    assert "flat-space Q-ball" in header
    # bs_omega is NOT overwritten for the flat-space Q-ball path.
    assert cfg.bs_omega == 0.8
