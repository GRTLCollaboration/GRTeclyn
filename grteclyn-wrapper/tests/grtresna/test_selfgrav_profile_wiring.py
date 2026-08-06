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


def test_selfgrav_lump_repaints_star_profile_not_gaussian() -> None:
    """The evolution repaint must evaluate the dressed star table for profile 4.

    Regression: phi0_at_radius used to fall through to the Gaussian envelope
    for PROFILE_SELFGRAV_BOUND, painting evolution matter that contradicted
    the constraint solve's tabulated star.
    """
    import numpy as np

    from grteclyn_wrapper.grtresna.profiles.envelope import phi0_at_radius

    lump = {
        "amp": 0.1,
        "width": 5.0,
        "center": (0.0, 0.0, 0.0),
        "profile": PROFILE_SELFGRAV_BOUND,
        "qball_mass": 1.0,
        "qball_lam": 0.0,
        "qball_mu": 0.0,
    }
    star = cached_selfgrav_profile(1.0, 0.0, 0.0, 0.1)
    r = np.linspace(0.0, 8.0, 33)
    painted = phi0_at_radius(r, lump, raw_amp=True)
    assert np.allclose(painted, np.asarray(star.eval_phi0(r)), rtol=1.0e-10)
    gaussian = 0.1 * np.exp(-r * r / (2.0 * 5.0 * 5.0))
    assert not np.allclose(painted, gaussian, rtol=1.0e-2)


def test_selfgrav_at_omega_params_write_back_amp_and_eigenvalue(tmp_path: Path) -> None:
    """Sextic selfgrav lump with a target frequency: phi_c becomes the OUTPUT.

    The writer must back-write amp := solved phi_c (painters rescale the table
    by amp/phi_c, so anything else de-tunes the star) and bs_omega := the
    physical eigenvalue (== the request within solver tolerance).  The repaint
    dispatch must resolve the same star from the lump keys alone.
    """
    import numpy as np

    from grteclyn_wrapper.grtresna.profiles.boson_star_ode import (
        cached_selfgrav_at_omega,
    )
    from grteclyn_wrapper.grtresna.profiles.envelope import phi0_at_radius
    from grteclyn_wrapper.grtresna.solver.params import write_grtresna_params

    cfg = GRTresnaConfig(
        matter_model="grtresna_complex_scalar",
        scalar_mass=1.0,
        scalar_lambda=10240.0,
        scalar_mu=21845333.0,
        bs_phi_c=0.08,
        bs_omega=0.0,
        lumps=[
            {
                "amp": 0.08,
                "width": 1.2,
                "center": (0.0, 0.0, 0.0),
                "profile": PROFILE_SELFGRAV_BOUND,
                "qball_mass": 1.0,
                "qball_lam": 10240.0,
                "qball_mu": 21845333.0,
                "qball_omega": 0.55,
            }
        ],
    )
    write_grtresna_params(cfg, tmp_path / "params.txt")

    star = cached_selfgrav_at_omega(1.0, 10240.0, 21845333.0, 0.55)
    assert cfg.lumps[0]["amp"] == pytest.approx(star.phi_c, rel=0, abs=0)
    assert cfg.bs_omega == pytest.approx(star.omega)
    assert cfg.bs_omega == pytest.approx(0.55, rel=2.0e-4)
    # Repaint resolves the identical star from the lump dict alone.
    r = np.linspace(0.0, 10.0, 41)
    painted = phi0_at_radius(r, cfg.lumps[0], raw_amp=True)
    assert np.allclose(painted, np.asarray(star.eval_phi0(r)), rtol=1.0e-10)


def test_selfgrav_pair_gets_per_sector_tables(tmp_path: Path) -> None:
    """A (+,-) dressed pair must paint two DIFFERENT stars, one table each.

    The phantom star (repulsive self-gravity, gravity_sign=-1) is slightly
    puffed with negative ADM mass; sharing the canonical table would seed it
    off-equilibrium.  Each lump gets its own profile_path; the canonical table
    stays the global fallback.
    """
    import numpy as np

    from grteclyn_wrapper.grtresna.profiles.boson_star_ode import (
        cached_selfgrav_at_omega,
    )
    from grteclyn_wrapper.grtresna.profiles.envelope import phi0_at_radius
    from grteclyn_wrapper.grtresna.solver.params import write_grtresna_params

    def lump(x: float, exotic: int) -> dict:
        return {
            "amp": 0.08,
            "width": 1.2,
            "center": (x, 0.0, 0.0),
            "profile": PROFILE_SELFGRAV_BOUND,
            "exotic": exotic,
            "qball_mass": 1.0,
            "qball_lam": 10240.0,
            "qball_mu": 21845333.0,
            "qball_omega": 0.55,
        }

    cfg = GRTresnaConfig(
        matter_model="grtresna_bicomplex_scalar",
        scalar_mass=1.0,
        scalar_lambda=10240.0,
        scalar_mu=21845333.0,
        bs_phi_c=0.08,
        bs_omega=0.0,
        lumps=[lump(4.0, 0), lump(-4.0, 1)],
    )
    write_grtresna_params(cfg, tmp_path / "params.txt")

    canon = cached_selfgrav_at_omega(1.0, 10240.0, 21845333.0, 0.55, 1.0)
    phant = cached_selfgrav_at_omega(1.0, 10240.0, 21845333.0, 0.55, -1.0)
    assert canon.adm_mass > 0.0 > phant.adm_mass

    # Per-sector amp back-writes (the two stars differ).
    assert cfg.lumps[0]["amp"] == pytest.approx(canon.phi_c, abs=0)
    assert cfg.lumps[1]["amp"] == pytest.approx(phant.phi_c, abs=0)
    assert cfg.lumps[0]["amp"] != cfg.lumps[1]["amp"]
    # Both sectors share the requested frequency (global-omega paint stays valid).
    assert cfg.lumps[0]["omega"] == pytest.approx(0.55, rel=2.0e-4)
    assert cfg.lumps[1]["omega"] == pytest.approx(0.55, rel=2.0e-4)

    text = (tmp_path / "params.txt").read_text(encoding="utf-8")
    assert "lump0_profile_path = " in text
    assert "lump1_profile_path = " in text
    assert (tmp_path / "qball_profile.dat").exists()
    assert (tmp_path / "qball_profile_exotic.dat").exists()
    assert "qball_profile_exotic.dat" in cfg.lumps[1]["profile_path"]
    head = (tmp_path / "qball_profile_exotic.dat").read_text().splitlines()[0]
    assert "ADM_mass=-" in head

    # Repaint dispatch resolves each sector's own star.
    r = np.linspace(0.0, 10.0, 41)
    assert np.allclose(
        phi0_at_radius(r, cfg.lumps[1], raw_amp=True),
        np.asarray(phant.eval_phi0(r)), rtol=1.0e-10,
    )
    assert not np.allclose(
        phi0_at_radius(r, cfg.lumps[1], raw_amp=True),
        np.asarray(canon.eval_phi0(r)), rtol=1.0e-3,
    )


def test_selfgrav_pair_mixed_frequency_equal_mass(tmp_path: Path) -> None:
    """Equal-|ADM| (+,-) pair: the phantom lump requests its OWN frequency.

    Per-lump qball_omega must produce per-lump solves, per-lump bs_omega in
    params.txt (the C++ painter's per-lump U(1) phase velocity), and a repaint
    that resolves each lump's own star.  This is the wiring the equal-mass
    Bondi runaway cell depends on: phantom at omega=0.566 weighs the same as
    canonical at omega=0.550.
    """
    import numpy as np

    from grteclyn_wrapper.grtresna.profiles.boson_star_ode import (
        cached_selfgrav_at_omega,
    )
    from grteclyn_wrapper.grtresna.profiles.envelope import phi0_at_radius
    from grteclyn_wrapper.grtresna.solver.params import write_grtresna_params

    def lump(x: float, exotic: int, omega: float) -> dict:
        return {
            "amp": 0.08,
            "width": 1.2,
            "center": (x, 0.0, 0.0),
            "profile": PROFILE_SELFGRAV_BOUND,
            "exotic": exotic,
            "qball_mass": 1.0,
            "qball_lam": 10240.0,
            "qball_mu": 21845333.0,
            "qball_omega": omega,
        }

    cfg = GRTresnaConfig(
        matter_model="grtresna_bicomplex_scalar",
        scalar_mass=1.0,
        scalar_lambda=10240.0,
        scalar_mu=21845333.0,
        bs_phi_c=0.08,
        bs_omega=0.0,
        lumps=[lump(4.0, 0, 0.55), lump(-4.0, 1, 0.566)],
    )
    write_grtresna_params(cfg, tmp_path / "params.txt")

    canon = cached_selfgrav_at_omega(1.0, 10240.0, 21845333.0, 0.55, 1.0)
    phant = cached_selfgrav_at_omega(1.0, 10240.0, 21845333.0, 0.566, -1.0)

    # Each lump carries its own solved star: amp AND phase velocity.
    assert cfg.lumps[0]["amp"] == pytest.approx(canon.phi_c, abs=0)
    assert cfg.lumps[1]["amp"] == pytest.approx(phant.phi_c, abs=0)
    assert cfg.lumps[0]["bs_omega"] == pytest.approx(canon.omega)
    assert cfg.lumps[1]["bs_omega"] == pytest.approx(phant.omega)
    assert cfg.lumps[1]["bs_omega"] == pytest.approx(0.566, rel=2.0e-4)
    # Global omega stays the canonical star's (fallback for unset lumps).
    assert cfg.bs_omega == pytest.approx(canon.omega)

    # The C++ painter contract: per-lump bs_omega lines in params.txt.
    text = (tmp_path / "params.txt").read_text(encoding="utf-8")
    assert f"lump0_bs_omega = {cfg.lumps[0]['bs_omega']}" in text
    assert f"lump1_bs_omega = {cfg.lumps[1]['bs_omega']}" in text
    assert "lump1_profile_path = " in text
    assert (tmp_path / "qball_profile_exotic.dat").exists()
    head = (tmp_path / "qball_profile_exotic.dat").read_text().splitlines()[0]
    assert "omega_eigenvalue=0.566" in head

    # Repaint resolves the phantom's 0.566 star, not a 0.55 one.
    r = np.linspace(0.0, 10.0, 41)
    assert np.allclose(
        phi0_at_radius(r, cfg.lumps[1], raw_amp=True),
        np.asarray(phant.eval_phi0(r)), rtol=1.0e-10,
    )

    # The Python bicomplex repaint uses the per-lump phase velocity for Pi_im.
    from grteclyn_wrapper.grtresna.fields.boson_star import _boosted_lump_fields

    px = np.zeros((1, 1, 1)) - 4.0
    py = np.zeros((1, 1, 1))
    pz = np.zeros((1, 1, 1))
    inv_lapse = np.ones((1, 1, 1))
    _, _, _, pi_im = _boosted_lump_fields(
        cfg.lumps[1], px, py, pz, cfg.bs_omega, inv_lapse
    )
    expected = -phant.omega * float(
        np.asarray(phant.eval_phi0(np.array([0.0]))).ravel()[0]
    )
    assert float(pi_im[0, 0, 0]) == pytest.approx(expected, rel=1.0e-10)
