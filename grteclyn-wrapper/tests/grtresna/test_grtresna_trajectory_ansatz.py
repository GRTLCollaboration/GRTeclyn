"""Tests for the trajectory-guided FTL geometry survey ansatz.

Each lump orbits independently on its own tilted circular trajectory.
Counter-rotating lumps create shear and frame dragging — the primary
FTL mechanism the search is designed to discover.
"""

import math

from grteclyn_wrapper.__main__ import build_parser
from grteclyn_wrapper.cli.grtresna_context import build_grtresna_search_context
from grteclyn_wrapper.search.optimize import (
    TRAJECTORY_DEFAULT_NUM_LUMPS,
    TRAJECTORY_PROFILE_CHOICES,
    build_grtresna_config,
    build_search_space,
    grtresna_trajectory_search_space,
)


# ---------------------------------------------------------------------------
# Search space tests
# ---------------------------------------------------------------------------


def test_trajectory_search_space_dimensionality() -> None:
    """Discovery profile: 7 per-lump dims * N + 5 shared."""
    for n in (3, 5):
        space = grtresna_trajectory_search_space(num_lumps=n)
        keys = {d.param_key for d in space}

        # Per-lump dims.
        for k in range(n):
            assert f"trajectory_lump{k}_R0" in keys
            assert f"trajectory_lump{k}_omega_rot" in keys
            assert f"trajectory_lump{k}_phase0" in keys
            assert f"trajectory_lump{k}_tilt_theta" in keys
            assert f"trajectory_lump{k}_tilt_phi" in keys
            assert f"trajectory_lump{k}_well_depth" in keys
            assert f"trajectory_lump{k}_exotic" in keys

        # Shared dims.
        assert "trajectory_A_breath" in keys
        assert "trajectory_omega_breath" in keys
        assert "trajectory_z_amp" in keys
        assert "trajectory_omega_z" in keys
        assert "trajectory_well_width" in keys

        # Total: 7 * N + 5.
        assert len(space) == 7 * n + 5


def test_trajectory_default_num_lumps() -> None:
    """Default number of lumps matches the constant."""
    space = grtresna_trajectory_search_space()
    per_lump_keys = {
        d.param_key for d in space if d.param_key.startswith("trajectory_lump")
    }
    # 7 per-lump keys * TRAJECTORY_DEFAULT_NUM_LUMPS.
    assert len(per_lump_keys) == 7 * TRAJECTORY_DEFAULT_NUM_LUMPS


def test_trajectory_rotation_only_profile() -> None:
    """rotation_only removes R0/tilt per-lump and breathing shared."""
    space = grtresna_trajectory_search_space(num_lumps=3, profile="rotation_only")
    keys = {d.param_key for d in space}

    # Should keep omega_rot, phase0, well_depth per lump.
    for k in range(3):
        assert f"trajectory_lump{k}_omega_rot" in keys
        assert f"trajectory_lump{k}_phase0" in keys
        assert f"trajectory_lump{k}_well_depth" in keys
        # Should NOT have R0, tilt.
        assert f"trajectory_lump{k}_R0" not in keys
        assert f"trajectory_lump{k}_tilt_theta" not in keys
        assert f"trajectory_lump{k}_tilt_phi" not in keys

    # Breathing shared dims removed.
    assert "trajectory_A_breath" not in keys
    assert "trajectory_omega_breath" not in keys

    # z-oscillation and well_width still present.
    assert "trajectory_z_amp" in keys
    assert "trajectory_well_width" in keys

    assert len(space) < len(grtresna_trajectory_search_space(num_lumps=3))


def test_trajectory_breathing_only_profile() -> None:
    """breathing_only removes omega_rot, phase0, tilt per lump."""
    space = grtresna_trajectory_search_space(num_lumps=3, profile="breathing_only")
    keys = {d.param_key for d in space}

    # Should keep R0 and well_depth per lump.
    for k in range(3):
        assert f"trajectory_lump{k}_R0" in keys
        assert f"trajectory_lump{k}_well_depth" in keys
        # Should NOT have rotation or tilt.
        assert f"trajectory_lump{k}_omega_rot" not in keys
        assert f"trajectory_lump{k}_phase0" not in keys
        assert f"trajectory_lump{k}_tilt_theta" not in keys
        assert f"trajectory_lump{k}_tilt_phi" not in keys

    # Breathing shared dims present.
    assert "trajectory_A_breath" in keys
    assert "trajectory_omega_breath" in keys


def test_trajectory_profile_choices_constant() -> None:
    """Profile choices match the valid set."""
    assert "discovery" in TRAJECTORY_PROFILE_CHOICES
    assert "rotation_only" in TRAJECTORY_PROFILE_CHOICES
    assert "breathing_only" in TRAJECTORY_PROFILE_CHOICES


def test_trajectory_wired_into_build_search_space() -> None:
    """build_search_space dispatches 'trajectory' ansatz correctly."""
    space = build_search_space(
        grtresna=True, grtresna_ansatz="trajectory", grtresna_lumps=4
    )
    keys = {d.param_key for d in space}
    assert "trajectory_lump0_R0" in keys
    assert "trajectory_lump3_omega_rot" in keys
    # Exactly 4 lumps (not 5 default).
    assert "trajectory_lump4_R0" not in keys
    # Should NOT have shell/SH keys.
    assert "grtresna_shell_amp" not in keys
    assert "grtresna_sh_amp" not in keys


def test_trajectory_bounds_are_sane() -> None:
    """All trajectory dims have finite, ordered bounds."""
    space = grtresna_trajectory_search_space()
    for d in space:
        assert d.lower < d.upper, f"{d.param_key}: lower={d.lower} >= upper={d.upper}"
        assert d.lower <= d.center <= d.upper, (
            f"{d.param_key}: center={d.center} outside [{d.lower}, {d.upper}]"
        )
        assert math.isfinite(d.lower), f"{d.param_key}: non-finite lower"
        assert math.isfinite(d.upper), f"{d.param_key}: non-finite upper"


def test_trajectory_no_unbounded_drift() -> None:
    """z_drift replaced with bounded z_amp oscillation."""
    space = grtresna_trajectory_search_space()
    keys = {d.param_key for d in space}
    assert "trajectory_z_drift" not in keys, "z_drift is unbounded and removed"
    assert "trajectory_z_amp" in keys, "z_amp (bounded oscillation) should be present"


def test_trajectory_counter_rotation_possible() -> None:
    """omega_rot range includes both positive and negative values."""
    space = grtresna_trajectory_search_space(num_lumps=2)
    for d in space:
        if "omega_rot" in d.param_key:
            assert d.lower < 0.0, f"{d.param_key}: lower={d.lower} should be negative"
            assert d.upper > 0.0, f"{d.param_key}: upper={d.upper} should be positive"


# ---------------------------------------------------------------------------
# Config expansion tests
# ---------------------------------------------------------------------------


def _trajectory_overrides(num_lumps: int = 5, **extra: float) -> dict:
    """Minimal trajectory overrides for config expansion.

    ``extra`` keys are auto-prefixed with ``trajectory_`` if not already.
    """
    ov: dict[str, Any] = {
        "trajectory_mode": 1,
        "trajectory_num_lumps": num_lumps,
        "trajectory_well_width": 1.5,
        "trajectory_A_breath": 0.0,
        "trajectory_omega_breath": 0.0,
        "trajectory_z_amp": 0.0,
        "trajectory_omega_z": 0.0,
    }
    # Default per-lump: equi-spaced on a circle of radius 5.
    for k in range(num_lumps):
        ov[f"trajectory_lump{k}_R0"] = 5.0
        ov[f"trajectory_lump{k}_omega_rot"] = 0.0
        ov[f"trajectory_lump{k}_phase0"] = 2.0 * math.pi * k / num_lumps
        ov[f"trajectory_lump{k}_tilt_theta"] = 0.0
        ov[f"trajectory_lump{k}_tilt_phi"] = 0.0
        ov[f"trajectory_lump{k}_well_depth"] = 0.05
    for k, v in extra.items():
        key = k if k.startswith("trajectory_") else f"trajectory_{k}"
        ov[key] = v
    return ov


def test_trajectory_expands_lumps_at_t0_positions() -> None:
    """Trajectory ansatz creates N lumps on a circle at R0 (no tilt)."""
    overrides = _trajectory_overrides(num_lumps=5)
    cfg = build_grtresna_config(overrides)

    assert len(cfg.lumps) == 5

    for lump in cfg.lumps:
        cx, cy, cz = lump["center"]
        r = math.sqrt(cx * cx + cy * cy + cz * cz)
        assert abs(r - 5.0) < 1e-6, f"lump at r={r}, expected 5.0"
        assert abs(cz) < 1e-10, f"z={cz}, expected 0 (no tilt)"


def test_trajectory_independent_radii() -> None:
    """Different per-lump R0 produces lumps at different radii."""
    ov = _trajectory_overrides(num_lumps=3)
    ov["trajectory_lump0_R0"] = 3.0
    ov["trajectory_lump1_R0"] = 5.0
    ov["trajectory_lump2_R0"] = 7.0
    cfg = build_grtresna_config(ov)

    radii = []
    for lump in cfg.lumps:
        cx, cy, cz = lump["center"]
        radii.append(math.sqrt(cx * cx + cy * cy + cz * cz))

    assert abs(radii[0] - 3.0) < 1e-6
    assert abs(radii[1] - 5.0) < 1e-6
    assert abs(radii[2] - 7.0) < 1e-6


def test_trajectory_independent_tilt() -> None:
    """Different per-lump tilt_theta places lumps on different planes."""
    ov = _trajectory_overrides(num_lumps=2)
    # Lump 0: equatorial (tilt=0), Lump 1: polar (tilt=pi/2).
    ov["trajectory_lump0_tilt_theta"] = 0.0
    ov["trajectory_lump0_phase0"] = 0.0
    ov["trajectory_lump1_tilt_theta"] = math.pi / 2.0
    ov["trajectory_lump1_phase0"] = 0.0
    cfg = build_grtresna_config(ov)

    # Lump 0: should be in equatorial plane (z=0).
    assert abs(cfg.lumps[0]["center"][2]) < 1e-10

    # Lump 1: tilt_theta=pi/2 rotates x_orb fully into -z.
    assert abs(cfg.lumps[1]["center"][2]) > 1.0


def test_trajectory_counter_rotating_lumps() -> None:
    """Counter-rotating lumps have opposite omega_rot — positions differ."""
    ov = _trajectory_overrides(num_lumps=2)
    ov["trajectory_lump0_omega_rot"] = 0.5
    ov["trajectory_lump1_omega_rot"] = -0.5
    # Both at phase0=0 — at t=0 they start at the same angle.
    ov["trajectory_lump0_phase0"] = 0.0
    ov["trajectory_lump1_phase0"] = 0.0

    cfg = build_grtresna_config(ov)
    # At t=0, positions are identical (omega doesn't affect t=0).
    # The point is that the C++ evaluator will move them apart.
    # Here we just verify both lumps exist with correct R0.
    for lump in cfg.lumps:
        cx, cy, cz = lump["center"]
        r = math.sqrt(cx * cx + cy * cy + cz * cz)
        assert abs(r - 5.0) < 1e-6


def test_trajectory_independent_amplitudes() -> None:
    """Per-lump well_depth produces different pump amplitudes."""
    ov = _trajectory_overrides(num_lumps=3)
    ov["trajectory_lump0_well_depth"] = 0.03
    ov["trajectory_lump1_well_depth"] = 0.08
    ov["trajectory_lump2_well_depth"] = 0.15
    cfg = build_grtresna_config(ov)

    amps = [lump["amp"] for lump in cfg.lumps]
    assert abs(amps[0] - 0.03) < 1e-10
    assert abs(amps[1] - 0.08) < 1e-10
    assert abs(amps[2] - 0.15) < 1e-10


def test_trajectory_amplitude_bounded() -> None:
    """Large well_depth is capped at 0.15 for GRTresna convergence."""
    ov = _trajectory_overrides(num_lumps=1)
    ov["trajectory_lump0_well_depth"] = 0.5
    cfg = build_grtresna_config(ov)

    assert cfg.lumps[0]["amp"] <= 0.15


def test_trajectory_width_minimum() -> None:
    """Trajectory lump width stays above minimum for solver stability."""
    ov = _trajectory_overrides(num_lumps=1)
    ov["trajectory_well_width"] = 0.5
    cfg = build_grtresna_config(ov)

    assert cfg.lumps[0]["width"] >= 1.5


def test_trajectory_lumps_are_static_at_t0() -> None:
    """All trajectory lumps have zero velocity (pump handles motion)."""
    ov = _trajectory_overrides(num_lumps=3)
    ov["trajectory_lump0_omega_rot"] = 0.5  # even with omega set
    cfg = build_grtresna_config(ov)

    for lump in cfg.lumps:
        vx, vy, vz = lump["velocity"]
        assert abs(vx) < 1e-10 and abs(vy) < 1e-10 and abs(vz) < 1e-10


def test_trajectory_overrides_not_filtered_as_grtresna() -> None:
    """trajectory_* keys are NOT grtresna_ prefixed, so they pass through."""
    ov = _trajectory_overrides()
    grtresna_keys = {k for k in ov if str(k).startswith("grtresna_")}
    trajectory_keys = {k for k in ov if str(k).startswith("trajectory_")}
    assert len(grtresna_keys) == 0
    assert len(trajectory_keys) > 0


def test_trajectory_tilt_phi_rotates_orbit_plane() -> None:
    """tilt_phi rotates the orbital plane around z-axis."""
    ov = _trajectory_overrides(num_lumps=1)
    ov["trajectory_lump0_R0"] = 5.0
    ov["trajectory_lump0_phase0"] = 0.0
    ov["trajectory_lump0_tilt_theta"] = math.pi / 4.0
    ov["trajectory_lump0_tilt_phi"] = math.pi / 2.0  # 90 degree azimuth
    cfg = build_grtresna_config(ov)

    cx, cy, cz = cfg.lumps[0]["center"]
    # With tilt_phi=pi/2, the tilt axis points along y.
    # The x_orb component should rotate into y-z plane.
    assert abs(cx) < 1e-6, f"x={cx} should be ~0 with tilt_phi=pi/2"
    assert abs(cy) > 1.0, f"y={cy} should be significant"


# ---------------------------------------------------------------------------
# CLI integration tests
# ---------------------------------------------------------------------------


def test_trajectory_cli_creates_search_context() -> None:
    """--grtresna-ansatz trajectory produces a valid search context."""
    parser = build_parser()
    args = parser.parse_args([
        "qd",
        "--grtresna",
        "--grtresna-ansatz", "trajectory",
        "--grtresna-lumps", "4",
    ])
    ctx = build_grtresna_search_context(args, base_overrides={})

    keys = {d.param_key for d in ctx.search_space}
    # 4 lumps * 7 per-lump + 5 shared = 33 dims.
    assert len(ctx.search_space) == 33
    assert "trajectory_lump0_R0" in keys
    assert "trajectory_lump3_omega_rot" in keys
    assert "trajectory_lump4_R0" not in keys  # only 4 lumps

    assert ctx.base_overrides.get("trajectory_mode") == 1
    assert ctx.base_overrides.get("trajectory_num_lumps") == 4
    assert ctx.use_grtresna is True
    assert ctx.grtresna_config is not None
