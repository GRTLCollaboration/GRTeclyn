"""Integration tests for Q-ball ODE seed + equilibrium amplitude wiring."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

from grteclyn_wrapper.grtresna.boson_star_fields import (
    _lump_phi0_at_radius,
    paint_bicomplex_fields_on_grid,
)
from grteclyn_wrapper.grtresna.boson_star_profile import (
    PROFILE_ODE_BOUND,
    PROFILE_SECH_BOUND,
    grtresna_lump_profile,
)
from grteclyn_wrapper.grtresna.qball_couplings import QBallCouplings
from grteclyn_wrapper.grtresna.solver import _lump_lines, write_grtresna_params
from grteclyn_wrapper.search.optimize.config import build_grtresna_config


def _trajectory_boson_overrides(**extra: object) -> dict:
    """Minimal eval-122-like boson trajectory overrides."""
    base = {
        "grtresna_matter_model": "grtresna_bicomplex_scalar",
        "grtresna_matter_sector": "boson_star",
        "grtresna_scalar_mass": 1.0,
        "grtresna_scalar_lambda": 640.0,
        "grtresna_scalar_mu": 85333.0,
        "grtresna_bs_omega": 0.4,
        "grtresna_boost_lumps": 1,
        "grtresna_boost_v_max": 0.8,
        "trajectory_mode": 1,
        "trajectory_num_lumps": 2,
        "trajectory_lump0_R0": 5.0,
        "trajectory_lump0_omega_rot": -0.5,
        "trajectory_lump0_well_depth": 0.15,
        "trajectory_lump0_exotic": 0,
        "trajectory_lump1_R0": 4.0,
        "trajectory_lump1_omega_rot": -0.4,
        "trajectory_lump1_well_depth": 0.15,
        "trajectory_lump1_exotic": 0,
    }
    base.update(extra)
    return base


def test_build_grtresna_config_ode_and_equilibrium_lumps() -> None:
    overrides = _trajectory_boson_overrides(
        grtresna_qball_equilibrium_amplitude=1,
        grtresna_qball_ode_profile=1,
    )
    cfg = build_grtresna_config(overrides)
    assert len(cfg.lumps) == 2
    stiff = QBallCouplings.stiff()
    for lump in cfg.lumps:
        assert lump["profile"] == PROFILE_ODE_BOUND
        assert lump["amp"] == pytest.approx(stiff.core_amplitude, rel=0.01)
        assert lump["qball_lam"] == 640.0
        assert lump["qball_mu"] == 85333.0


def test_ode_lump_profile_differs_from_sech_at_same_amp() -> None:
    stiff = QBallCouplings.stiff()
    amp = stiff.core_amplitude
    width = stiff.bound_state_width
    r = np.linspace(0.5, 4.0, 32)
    ode_lump = {
        "amp": amp,
        "width": width,
        "profile": PROFILE_ODE_BOUND,
        "qball_mass": stiff.mass,
        "qball_lam": stiff.lam,
        "qball_mu": stiff.mu,
        "qball_omega": stiff.omega,
    }
    sech_lump = {"amp": amp, "width": width, "profile": PROFILE_SECH_BOUND}
    ode = _lump_phi0_at_radius(r, ode_lump)
    sech = _lump_phi0_at_radius(r, sech_lump)
    assert not np.allclose(ode, sech, rtol=0.05)


def test_paint_bicomplex_ode_lump_on_grid() -> None:
    stiff = QBallCouplings.stiff()
    amp = stiff.core_amplitude
    nx = ny = nz = 32
    dx = np.array([1.0, 1.0, 1.0])
    origin = np.array([0.0, 0.0, 0.0])
    center = (16.0, 16.0, 16.0)
    comp_names = ["chi"] * 25 + ["lapse"] + ["chi"] * 7
    data = np.zeros((nz, ny, nx, 33), dtype=np.float64)
    data[:, :, :, 25] = 1.0  # lapse

    lump = {
        "amp": amp,
        "width": stiff.bound_state_width,
        "center": center,
        "velocity": (0.0, 0.0, 0.0),
        "profile": PROFILE_ODE_BOUND,
        "exotic": 0,
        "qball_mass": stiff.mass,
        "qball_lam": stiff.lam,
        "qball_mu": stiff.mu,
        "qball_omega": stiff.omega,
    }
    painted, names = paint_bicomplex_fields_on_grid(
        data, comp_names, dx, origin, [lump], omega=stiff.omega
    )
    phi_idx = names.index("phi")
    phi = painted[:, :, :, phi_idx]
    peak = float(phi.max())
    assert peak == pytest.approx(amp, rel=0.05)
    assert peak < 0.10  # equilibrium stiff core, not 0.15 sech overshoot


def test_grtresna_params_map_ode_profile_to_sech() -> None:
    """GRTresna C++ only knows profiles 0–2; ODE (3) stays in-memory for gridinit."""
    overrides = _trajectory_boson_overrides(
        grtresna_qball_equilibrium_amplitude=1,
        grtresna_qball_ode_profile=1,
    )
    cfg = build_grtresna_config(overrides)
    assert all(lump["profile"] == PROFILE_ODE_BOUND for lump in cfg.lumps)

    lump_text = "\n".join(_lump_lines(cfg))
    assert "lump0_profile = 2" in lump_text
    assert "lump1_profile = 2" in lump_text
    assert "lump0_profile = 3" not in lump_text

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "params.txt"
        write_grtresna_params(cfg, path)
        text = path.read_text()
        assert "lump0_profile = 2" in text
        assert "lump1_profile = 2" in text

    assert grtresna_lump_profile(PROFILE_ODE_BOUND) == PROFILE_SECH_BOUND
    assert grtresna_lump_profile(PROFILE_SECH_BOUND) == PROFILE_SECH_BOUND
