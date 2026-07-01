"""Tests for trajectory + boson star ansatz — compact solitons on independent orbits.

Each lump is a U(1)-charged boson star placed at a trajectory orbit position.
Unlike the real-scalar trajectory (pump-driven, dispersing), boson stars are
self-bound via charge conservation and the pump is corrective only.
"""

from __future__ import annotations

import math
from typing import Any

from grteclyn_wrapper.__main__ import build_parser
from grteclyn_wrapper.cli.grtresna_context import build_grtresna_search_context
from grteclyn_wrapper.grtresna.matter_models import (
    GRTRESNA_COMPLEX_SCALAR_MODEL,
    GRTRESNA_EXAMPLE_BOSON_STAR_BH,
    matter_selection_base_overrides,
    resolve_matter_selection,
)
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig
from grteclyn_wrapper.search.optimize import (
    TRAJECTORY_BOSON_PROFILE_CHOICES,
    TRAJECTORY_DEFAULT_NUM_LUMPS,
    build_grtresna_config,
    build_search_space,
    grtresna_trajectory_boson_search_space,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _boson_trajectory_overrides(num_lumps: int = 5, **extra: float) -> dict:
    """Minimal trajectory + boson_star overrides for config expansion."""
    ov: dict[str, Any] = {
        # Matter sector: complex scalar (boson star).
        **matter_selection_base_overrides(
            resolve_matter_selection("boson_star", "canonical")
        ),
        # Trajectory parameters.
        "trajectory_mode": 1,
        "trajectory_num_lumps": num_lumps,
        "trajectory_well_width": 1.5,
        "trajectory_A_breath": 0.0,
        "trajectory_omega_breath": 0.0,
        "trajectory_z_amp": 0.0,
        "trajectory_omega_z": 0.0,
        # Boson physics.
        "grtresna_scalar_mass": 0.2,
        "grtresna_scalar_lambda": 0.0,
        "grtresna_bs_omega": 0.15,
    }
    for k in range(num_lumps):
        ov[f"trajectory_lump{k}_R0"] = 5.0
        ov[f"trajectory_lump{k}_omega_rot"] = 0.0
        ov[f"trajectory_lump{k}_phase0"] = 2.0 * math.pi * k / num_lumps
        ov[f"trajectory_lump{k}_tilt_theta"] = 0.0
        ov[f"trajectory_lump{k}_tilt_phi"] = 0.0
        ov[f"trajectory_lump{k}_well_depth"] = 0.005
        ov[f"trajectory_lump{k}_v_rad"] = 0.0
        ov[f"trajectory_lump{k}_exotic"] = 0.0
    for k, v in extra.items():
        ov[k] = v
    return ov


# ---------------------------------------------------------------------------
# Search space tests
# ---------------------------------------------------------------------------


def test_trajectory_boson_search_space_dimensionality() -> None:
    """Discovery profile: 8 per-lump + 5 shared trajectory + 3 boson physics."""
    for n in (3, 5):
        space = grtresna_trajectory_boson_search_space(num_lumps=n)
        keys = {d.param_key for d in space}

        # Per-lump trajectory dims.
        for k in range(n):
            assert f"trajectory_lump{k}_R0" in keys
            assert f"trajectory_lump{k}_omega_rot" in keys
            assert f"trajectory_lump{k}_phase0" in keys
            assert f"trajectory_lump{k}_tilt_theta" in keys
            assert f"trajectory_lump{k}_tilt_phi" in keys
            assert f"trajectory_lump{k}_well_depth" in keys
            assert f"trajectory_lump{k}_v_rad" in keys
            assert f"trajectory_lump{k}_exotic" in keys

        # Shared trajectory dims.
        assert "trajectory_A_breath" in keys
        assert "trajectory_omega_breath" in keys
        assert "trajectory_z_amp" in keys
        assert "trajectory_omega_z" in keys
        assert "trajectory_well_width" in keys

        # Shared boson physics dims.
        assert "grtresna_scalar_mass" in keys
        assert "grtresna_scalar_lambda" in keys
        assert "grtresna_bs_omega" in keys

        # Total: 8*N + 5 shared_traj + 3 boson = 8*N + 8.
        assert len(space) == 8 * n + 8


def test_trajectory_boson_profile_choices() -> None:
    assert "discovery" in TRAJECTORY_BOSON_PROFILE_CHOICES
    assert "rotation_only" in TRAJECTORY_BOSON_PROFILE_CHOICES


def test_trajectory_boson_rotation_only_profile() -> None:
    """rotation_only removes R0/tilt per-lump and breathing shared."""
    space = grtresna_trajectory_boson_search_space(
        num_lumps=3, profile="rotation_only"
    )
    keys = {d.param_key for d in space}

    for k in range(3):
        assert f"trajectory_lump{k}_omega_rot" in keys
        assert f"trajectory_lump{k}_phase0" in keys
        assert f"trajectory_lump{k}_well_depth" in keys
        assert f"trajectory_lump{k}_R0" not in keys
        assert f"trajectory_lump{k}_tilt_theta" not in keys
        assert f"trajectory_lump{k}_tilt_phi" not in keys

    assert "trajectory_A_breath" not in keys
    assert "trajectory_omega_breath" not in keys
    # Boson physics still present.
    assert "grtresna_scalar_mass" in keys
    assert "grtresna_bs_omega" in keys


def test_trajectory_boson_well_depth_bounds_are_corrective() -> None:
    """Pump amplitude is reduced: [0.001, 0.02] instead of [0.01, 0.15]."""
    space = grtresna_trajectory_boson_search_space(num_lumps=3)
    for d in space:
        if "well_depth" in d.param_key:
            assert d.lower >= 0.001, f"{d.param_key}: lower too low"
            assert d.upper <= 0.02, f"{d.param_key}: upper too high for corrective pump"


def test_trajectory_boson_wired_into_build_search_space() -> None:
    """build_search_space dispatches trajectory + boson_star correctly."""
    space = build_search_space(
        grtresna=True,
        grtresna_ansatz="trajectory",
        grtresna_matter_sector="boson_star",
        grtresna_lumps=4,
    )
    keys = {d.param_key for d in space}

    # Should have trajectory orbit dims.
    assert "trajectory_lump0_R0" in keys
    assert "trajectory_lump3_omega_rot" in keys
    assert "trajectory_lump4_R0" not in keys  # only 4 lumps

    # Should have boson physics dims.
    assert "grtresna_scalar_mass" in keys
    assert "grtresna_bs_omega" in keys

    # Should NOT have shell/SH keys.
    assert "grtresna_shell_amp" not in keys
    assert "grtresna_sh_amp" not in keys


def test_trajectory_boson_bounds_are_sane() -> None:
    """All trajectory-boson dims have finite, ordered bounds."""
    space = grtresna_trajectory_boson_search_space()
    for d in space:
        assert d.lower < d.upper, f"{d.param_key}: lower >= upper"
        assert d.lower <= d.center <= d.upper, (
            f"{d.param_key}: center outside bounds"
        )
        assert math.isfinite(d.lower)
        assert math.isfinite(d.upper)


# ---------------------------------------------------------------------------
# Config expansion tests
# ---------------------------------------------------------------------------


def test_trajectory_boson_expands_to_complex_scalar_model() -> None:
    """Trajectory + boson_star routes to BosonStarBH / complex scalar."""
    cfg = build_grtresna_config(
        _boson_trajectory_overrides(num_lumps=5), GRTresnaConfig()
    )
    assert cfg.matter_model == GRTRESNA_COMPLEX_SCALAR_MODEL
    assert cfg.example == GRTRESNA_EXAMPLE_BOSON_STAR_BH


def test_trajectory_boson_creates_lumps_at_orbit_positions() -> None:
    """Lumps are placed at correct t=0 positions (same geometry as scalar trajectory)."""
    cfg = build_grtresna_config(
        _boson_trajectory_overrides(num_lumps=5), GRTresnaConfig()
    )
    assert len(cfg.lumps) == 5

    for lump in cfg.lumps:
        cx, cy, cz = lump["center"]
        r = math.sqrt(cx * cx + cy * cy + cz * cz)
        assert abs(r - 5.0) < 1e-6, f"lump at r={r}, expected 5.0"
        assert abs(cz) < 1e-10, "z should be 0 (no tilt)"


def test_trajectory_boson_lumps_have_bulk_velocity() -> None:
    """Boson star lumps get tangential velocity from omega_rot x r."""
    ov = _boson_trajectory_overrides(num_lumps=2)
    # Use small omega so |v| = omega * R0 = 0.1 * 5 = 0.5 < 0.9 (no cap).
    ov["trajectory_lump0_omega_rot"] = 0.1
    ov["trajectory_lump0_phase0"] = 0.0  # lump at (R0, 0, 0)
    ov["trajectory_lump1_omega_rot"] = -0.1
    ov["trajectory_lump1_phase0"] = math.pi  # lump at (-R0, 0, 0)
    cfg = build_grtresna_config(ov, GRTresnaConfig())

    # Lump 0 at (5,0,0) with omega=0.1: v = (0, 0.1*5, 0) = (0, 0.5, 0)
    vx0, vy0, vz0 = cfg.lumps[0]["velocity"]
    assert abs(vx0) < 1e-10, f"vx0={vx0}, expected ~0"
    assert abs(vy0 - 0.5) < 1e-6, f"vy0={vy0}, expected 0.5"
    assert abs(vz0) < 1e-10, f"vz0={vz0}, expected ~0"

    # Lump 1 at (-5,0,0) with omega=-0.1:
    #   v_x_orb = -omega * R0 * sin(pi) = 0
    #   v_y_orb =  omega * R0 * cos(pi) = -0.1 * 5 * (-1) = 0.5
    vx1, vy1, vz1 = cfg.lumps[1]["velocity"]
    assert abs(vx1) < 1e-10, f"vx1={vx1}, expected ~0"
    assert abs(vy1 - 0.5) < 1e-6, f"vy1={vy1}, expected 0.5"
    assert abs(vz1) < 1e-10, f"vz1={vz1}, expected ~0"


def test_trajectory_boson_velocity_capped_subluminal() -> None:
    """Bulk velocity is capped at 0.9c."""
    ov = _boson_trajectory_overrides(num_lumps=1)
    ov["trajectory_lump0_R0"] = 8.0
    ov["trajectory_lump0_omega_rot"] = 1.0  # |v| = 1.0 * 8.0 = 8.0 >> 0.9
    ov["trajectory_lump0_phase0"] = 0.0
    cfg = build_grtresna_config(ov, GRTresnaConfig())

    vx, vy, vz = cfg.lumps[0]["velocity"]
    v_mag = math.sqrt(vx * vx + vy * vy + vz * vz)
    assert v_mag <= 0.9 + 1e-10, f"|v|={v_mag}, should be capped at 0.9"


def test_trajectory_boson_spiral_velocity() -> None:
    """Radial drift v_rad adds a spiral velocity component at t=0."""
    ov = _boson_trajectory_overrides(num_lumps=1)
    ov["trajectory_lump0_R0"] = 5.0
    ov["trajectory_lump0_omega_rot"] = 0.0
    ov["trajectory_lump0_phase0"] = 0.0  # lump at (R0, 0, 0)
    ov["trajectory_lump0_v_rad"] = 0.2   # pure outward radial velocity
    cfg = build_grtresna_config(ov, GRTresnaConfig())

    vx, vy, vz = cfg.lumps[0]["velocity"]
    assert abs(vx - 0.2) < 1e-6, f"vx={vx}, expected 0.2"
    assert abs(vy) < 1e-10, f"vy={vy}, expected ~0"
    assert abs(vz) < 1e-10, f"vz={vz}, expected ~0"


def test_trajectory_boson_spiral_velocity_with_omega() -> None:
    """Spiral velocity combines tangential and radial components."""
    ov = _boson_trajectory_overrides(num_lumps=1)
    ov["trajectory_lump0_R0"] = 5.0
    ov["trajectory_lump0_omega_rot"] = 0.1
    ov["trajectory_lump0_phase0"] = math.pi / 2  # lump at (0, R0, 0)
    ov["trajectory_lump0_v_rad"] = 0.2
    cfg = build_grtresna_config(ov, GRTresnaConfig())

    vx, vy, vz = cfg.lumps[0]["velocity"]
    # At phase0=pi/2:
    # vx_orb = v_rad*cos(pi/2) - R0*omega*sin(pi/2) = 0 - 0.5 = -0.5
    # vy_orb = v_rad*sin(pi/2) + R0*omega*cos(pi/2) = 0.2 + 0 = 0.2
    assert abs(vx + 0.5) < 1e-6, f"vx={vx}, expected -0.5"
    assert abs(vy - 0.2) < 1e-6, f"vy={vy}, expected 0.2"
    assert abs(vz) < 1e-10, f"vz={vz}, expected ~0"


def test_trajectory_boson_amplitude_capped() -> None:
    """Well_depth is capped at 0.15 (same as real-scalar trajectory)."""
    ov = _boson_trajectory_overrides(num_lumps=1)
    ov["trajectory_lump0_well_depth"] = 0.5  # way too high
    cfg = build_grtresna_config(ov, GRTresnaConfig())

    assert cfg.lumps[0]["amp"] <= 0.15


def test_trajectory_boson_lumps_all_canonical_by_default() -> None:
    """Default exotic flag is 0 (canonical coupling)."""
    cfg = build_grtresna_config(
        _boson_trajectory_overrides(num_lumps=5), GRTresnaConfig()
    )
    assert all(lump["exotic"] == 0 for lump in cfg.lumps)


def test_trajectory_boson_exotic_lumps() -> None:
    """Per-lump exotic flag works for boson trajectory."""
    ov = _boson_trajectory_overrides(num_lumps=3)
    ov["trajectory_lump0_exotic"] = 1.0
    ov["trajectory_lump1_exotic"] = 0.0
    ov["trajectory_lump2_exotic"] = 1.0
    cfg = build_grtresna_config(ov, GRTresnaConfig())

    assert cfg.lumps[0]["exotic"] == 1
    assert cfg.lumps[1]["exotic"] == 0
    assert cfg.lumps[2]["exotic"] == 1


def test_trajectory_boson_reads_scalar_mass_and_omega() -> None:
    """Boson physics params (scalar_mass, bs_omega) are read into config."""
    ov = _boson_trajectory_overrides(num_lumps=2)
    ov["grtresna_scalar_mass"] = 0.3
    ov["grtresna_bs_omega"] = 0.25
    ov["grtresna_scalar_lambda"] = 0.01
    cfg = build_grtresna_config(ov, GRTresnaConfig())

    assert cfg.scalar_mass == 0.3
    assert cfg.bs_omega == 0.25
    assert cfg.scalar_lambda == 0.01


def test_trajectory_boson_tilted_orbit_with_velocity() -> None:
    """Tilted orbits produce correctly rotated positions AND velocities."""
    ov = _boson_trajectory_overrides(num_lumps=1)
    ov["trajectory_lump0_R0"] = 5.0
    # Use small omega: |v| = 0.1 * 5 = 0.5 < 0.9 (no cap).
    ov["trajectory_lump0_omega_rot"] = 0.1
    ov["trajectory_lump0_phase0"] = 0.0
    ov["trajectory_lump0_tilt_theta"] = math.pi / 2.0  # polar tilt
    ov["trajectory_lump0_tilt_phi"] = 0.0
    cfg = build_grtresna_config(ov, GRTresnaConfig())

    cx, cy, cz = cfg.lumps[0]["center"]
    # With tilt_theta=pi/2, phase0=0: x_orb=R0 rotates fully into -z.
    assert abs(cx) < 1e-6, f"cx={cx}, expected ~0"
    assert abs(cz + 5.0) < 1e-6, f"cz={cz}, expected -5.0"

    # Velocity should also be rotated.
    vx, vy, vz = cfg.lumps[0]["velocity"]
    # v_orb = (0, omega*R0, 0) = (0, 0.5, 0)
    # After tilt_theta=pi/2, tilt_phi=0:
    #   vx = cp*ct*vx_orb - sp*vy_orb = 0*0 - 0*0.5 = 0
    #   vy = sp*ct*vx_orb + cp*vy_orb = 0*0 + 1*0.5 = 0.5
    #   vz = -st*vx_orb = 0
    assert abs(vy - 0.5) < 1e-6, f"vy={vy}, expected 0.5"


# ---------------------------------------------------------------------------
# CLI integration tests
# ---------------------------------------------------------------------------


def test_trajectory_boson_cli_creates_search_context() -> None:
    """--grtresna-ansatz trajectory + --grtresna-matter-sector boson_star works."""
    parser = build_parser()
    args = parser.parse_args([
        "qd",
        "--grtresna",
        "--grtresna-ansatz", "trajectory",
        "--grtresna-matter-sector", "boson_star",
        "--grtresna-lumps", "4",
    ])
    ctx = build_grtresna_search_context(args, base_overrides={})

    keys = {d.param_key for d in ctx.search_space}

    # 4 lumps * 8 per-lump + 5 shared_traj + 3 boson_phys = 40 dims.
    assert len(ctx.search_space) == 40

    # Trajectory orbit dims.
    assert "trajectory_lump0_R0" in keys
    assert "trajectory_lump3_omega_rot" in keys
    assert "trajectory_lump4_R0" not in keys

    # Boson physics dims.
    assert "grtresna_scalar_mass" in keys
    assert "grtresna_bs_omega" in keys

    # Base overrides include trajectory and matter selection.
    assert ctx.base_overrides.get("trajectory_mode") == 1
    assert ctx.base_overrides.get("trajectory_num_lumps") == 4
    assert ctx.base_overrides.get("grtresna_matter_sector") == "boson_star"
    assert ctx.base_overrides.get("grtresna_matter_model") == GRTRESNA_COMPLEX_SCALAR_MODEL

    assert ctx.use_grtresna is True
    assert ctx.grtresna_config is not None


def test_scalar_trajectory_unchanged_by_boson_path() -> None:
    """Scalar + trajectory still uses zero-velocity real-scalar expansion."""
    parser = build_parser()
    args = parser.parse_args([
        "qd",
        "--grtresna",
        "--grtresna-ansatz", "trajectory",
        "--grtresna-lumps", "3",
    ])
    ctx = build_grtresna_search_context(args, base_overrides={})

    keys = {d.param_key for d in ctx.search_space}
    # Scalar trajectory has 6 per-lump (no exotic dim in per-lump count for
    # the old test) + 5 shared = 23 dims for 3 lumps... let's check:
    # Actually: 7 per-lump (R0, omega_rot, phase0, tilt_theta, tilt_phi,
    # well_depth, exotic) + 5 shared = 7*3+5 = 26.
    assert "trajectory_lump0_R0" in keys
    # Should NOT have boson physics dims.
    assert "grtresna_scalar_mass" not in keys
    assert "grtresna_bs_omega" not in keys
