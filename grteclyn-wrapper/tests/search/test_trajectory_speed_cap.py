"""Tests for trajectory constraints: normalized speed fractions + retrograde.

The trajectory search space now treats ``trajectory_lump*_omega_rot`` and
``trajectory_lump*_v_rad`` as *normalized fractions* of ``trajectory_v_max``:

  * omega_rot_norm = -1..1  ->  tangential speed v_t = omega_rot_norm * v_max
  * v_rad_norm     = -1..1  ->  radial speed      v_r = v_rad_norm     * v_max

The combined normalized speed is clamped to the unit disk, then converted to
physical quantities:

  * omega_rot = v_t / R0
  * v_rad     = v_r

This gives the optimizer a smooth 2-D speed disk instead of a hard-cornered
box where every non-zero omega_rot hits the speed boundary.
"""

from __future__ import annotations

import math

from grteclyn_wrapper.search.optimize.candidates import (
    TRAJECTORY_R_MIN_DEFAULT,
    TRAJECTORY_STOP_TIME_DEFAULT,
    TRAJECTORY_V_MAX_DEFAULT,
    _clamp_trajectory_speed,
    _clamp_trajectory_spiral_radius,
    _enforce_retrograde,
    _min_spiral_v_rad,
    _vector_to_overrides,
)
from grteclyn_wrapper.search.optimize.dimension import SearchDimension


def _physical_speed(overrides: dict, lump_idx: int) -> tuple[float, float, float]:
    """Return (v_tan, v_rad, v_total) for a lump after conversion."""
    r0 = float(overrides.get(f"trajectory_lump{lump_idx}_R0", 0.0))
    omega = float(overrides.get(f"trajectory_lump{lump_idx}_omega_rot", 0.0))
    v_rad = float(overrides.get(f"trajectory_lump{lump_idx}_v_rad", 0.0))
    v_t = abs(omega) * r0
    return v_t, v_rad, math.sqrt(v_t * v_t + v_rad * v_rad)


# ---------------------------------------------------------------------------
# Normalized speed conversion
# ---------------------------------------------------------------------------


def test_clamp_converts_normalized_omega_to_physical() -> None:
    """omega_rot is now a normalized tangential-speed fraction."""
    ov = {
        "trajectory_lump0_R0": 4.0,
        "trajectory_lump0_omega_rot": 0.5,  # -> v_t = 0.5 * 0.3 = 0.15
    }
    _clamp_trajectory_speed(ov)
    v_t, v_rad, v_total = _physical_speed(ov, 0)
    assert math.isclose(v_t, 0.5 * TRAJECTORY_V_MAX_DEFAULT, rel_tol=1e-9)
    assert math.isclose(v_rad, 0.0, rel_tol=1e-9)
    assert math.isclose(v_total, 0.15, rel_tol=1e-9)
    # Physical omega_rot = v_t / R0.
    assert math.isclose(
        ov["trajectory_lump0_omega_rot"], 0.15 / 4.0, rel_tol=1e-9
    )


def test_clamp_caps_normalized_speed_to_unit_disk() -> None:
    """Combined normalized speed > 1 is scaled to the unit disk."""
    ov = {
        "trajectory_lump0_R0": 4.0,
        "trajectory_lump0_omega_rot": 0.8,
        "trajectory_lump0_v_rad": 0.8,
    }
    _clamp_trajectory_speed(ov)
    v_t, v_rad, v_total = _physical_speed(ov, 0)
    # norm_total = sqrt(0.8^2 + 0.8^2) = sqrt(1.28) ≈ 1.131
    # scale = 1 / norm_total ≈ 0.884
    # v_t = v_rad = 0.8 * scale * 0.3 ≈ 0.212
    assert math.isclose(v_total, TRAJECTORY_V_MAX_DEFAULT, rel_tol=1e-9)
    assert math.isclose(v_t, v_rad, rel_tol=1e-9)
    assert ov["trajectory_lump0_omega_rot"] > 0.0
    assert ov["trajectory_lump0_v_rad"] > 0.0


def test_clamp_respects_custom_v_max() -> None:
    ov = {
        "trajectory_lump0_R0": 5.0,
        "trajectory_lump0_omega_rot": 0.5,
        "trajectory_v_max": 0.1,
    }
    _clamp_trajectory_speed(ov)
    v_t, _, _ = _physical_speed(ov, 0)
    assert math.isclose(v_t, 0.05, rel_tol=1e-9)
    assert math.isclose(ov["trajectory_lump0_omega_rot"], 0.05 / 5.0, rel_tol=1e-9)


def test_clamp_disabled_when_v_max_non_positive() -> None:
    ov = {
        "trajectory_lump0_R0": 8.0,
        "trajectory_lump0_omega_rot": 0.5,
        "trajectory_v_max": 0.0,
    }
    _clamp_trajectory_speed(ov)
    # With v_max=0 the conversion is disabled; overrides stay as normalized.
    assert ov["trajectory_lump0_omega_rot"] == 0.5


def test_clamp_noop_for_non_trajectory_overrides() -> None:
    ov = {"grtresna_lump0_amp": 0.1, "grtresna_scalar_mass": 0.2}
    before = dict(ov)
    _clamp_trajectory_speed(ov)
    assert ov == before


def test_vector_to_overrides_applies_speed_cap() -> None:
    dims = [
        SearchDimension("trajectory_lump0_R0", 1.5, 8.0, 5.0),
        SearchDimension("trajectory_lump0_omega_rot", -1.0, 1.0, 0.0),
    ]
    overrides = _vector_to_overrides([8.0, 1.0], dims, {})
    v_t, _, _ = _physical_speed(overrides, 0)
    assert math.isclose(v_t, TRAJECTORY_V_MAX_DEFAULT, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# Retrograde enforcement
# ---------------------------------------------------------------------------


def test_retrograde_negates_prograde_when_enabled() -> None:
    ov = {
        "trajectory_lump0_omega_rot": 0.5,
        "trajectory_lump1_omega_rot": 0.3,
        "trajectory_retrograde_only": 1,
    }
    _enforce_retrograde(ov)
    assert ov["trajectory_lump0_omega_rot"] == -0.5
    assert ov["trajectory_lump1_omega_rot"] == -0.3


def test_retrograde_leaves_negative_untouched() -> None:
    ov = {
        "trajectory_lump0_omega_rot": -0.7,
        "trajectory_retrograde_only": 1,
    }
    _enforce_retrograde(ov)
    assert ov["trajectory_lump0_omega_rot"] == -0.7


def test_retrograde_leaves_zero_untouched() -> None:
    ov = {
        "trajectory_lump0_omega_rot": 0.0,
        "trajectory_retrograde_only": 1,
    }
    _enforce_retrograde(ov)
    assert ov["trajectory_lump0_omega_rot"] == 0.0


def test_retrograde_noop_when_disabled() -> None:
    ov = {
        "trajectory_lump0_omega_rot": 0.5,
        "trajectory_retrograde_only": 0,
    }
    _enforce_retrograde(ov)
    assert ov["trajectory_lump0_omega_rot"] == 0.5


def test_retrograde_noop_when_flag_absent() -> None:
    ov = {"trajectory_lump0_omega_rot": 0.5}
    _enforce_retrograde(ov)
    assert ov["trajectory_lump0_omega_rot"] == 0.5


def test_vector_to_overrides_applies_retrograde() -> None:
    dims = [
        SearchDimension("trajectory_lump0_R0", 1.5, 8.0, 5.0),
        SearchDimension("trajectory_lump0_omega_rot", -1.0, 1.0, 0.0),
    ]
    overrides = _vector_to_overrides(
        [4.0, 0.6], dims, {"trajectory_retrograde_only": 1}
    )
    assert overrides["trajectory_lump0_omega_rot"] <= 0.0
    v_t, _, _ = _physical_speed(overrides, 0)
    assert math.isclose(v_t, 0.6 * TRAJECTORY_V_MAX_DEFAULT, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# Spiral (radial drift) speed cap
# ---------------------------------------------------------------------------


def test_clamp_caps_combined_tangential_and_radial_speed() -> None:
    ov = {
        "trajectory_lump0_R0": 4.0,
        "trajectory_lump0_omega_rot": 0.8,
        "trajectory_lump0_v_rad": 0.8,
    }
    _clamp_trajectory_speed(ov)
    v_t, v_rad, v_total = _physical_speed(ov, 0)
    assert math.isclose(v_total, TRAJECTORY_V_MAX_DEFAULT, rel_tol=1e-9)
    # Both components scaled equally on the normalized disk.
    assert math.isclose(v_t, v_rad, rel_tol=1e-9)


def test_clamp_leaves_subluminal_spiral_untouched() -> None:
    ov = {
        "trajectory_lump0_R0": 4.0,
        "trajectory_lump0_omega_rot": 0.5,  # norm_total = 0.5 < 1
        "trajectory_lump0_v_rad": 0.0,
    }
    _clamp_trajectory_speed(ov)
    v_t, _, _ = _physical_speed(ov, 0)
    assert math.isclose(v_t, 0.5 * TRAJECTORY_V_MAX_DEFAULT, rel_tol=1e-9)
    assert math.isclose(
        ov["trajectory_lump0_omega_rot"],
        0.5 * TRAJECTORY_V_MAX_DEFAULT / 4.0,
        rel_tol=1e-9,
    )


def test_clamp_scales_omega_only_when_v_rad_missing() -> None:
    """When v_rad is absent it is treated as 0; only omega_rot is converted."""
    ov = {
        "trajectory_lump0_R0": 4.0,
        "trajectory_lump0_omega_rot": 0.5,
    }
    _clamp_trajectory_speed(ov)
    v_t, _, _ = _physical_speed(ov, 0)
    assert math.isclose(v_t, 0.5 * TRAJECTORY_V_MAX_DEFAULT, rel_tol=1e-9)
    assert "trajectory_lump0_v_rad" not in ov


# ---------------------------------------------------------------------------
# Spiral radius floor (inward drift vs stop_time)
# ---------------------------------------------------------------------------


def test_min_spiral_v_rad_worst_case_inward() -> None:
    # R0=1.5, stop=16, no breathing -> v_rad >= (0.1 - 1.5) / 16
    floor = _min_spiral_v_rad(
        1.5,
        a_breath=0.0,
        stop_time=16.0,
        r_min=TRAJECTORY_R_MIN_DEFAULT,
    )
    assert math.isclose(floor, -0.0875)


def test_clamp_spiral_radius_limits_inward_drift() -> None:
    ov = {
        "trajectory_lump0_R0": 1.5,
        "trajectory_lump0_v_rad": -0.3,
        "stop_time": TRAJECTORY_STOP_TIME_DEFAULT,
    }
    _clamp_trajectory_spiral_radius(ov)
    assert ov["trajectory_lump0_v_rad"] > -0.3
    assert ov["trajectory_lump0_v_rad"] >= _min_spiral_v_rad(
        1.5,
        a_breath=0.0,
        stop_time=TRAJECTORY_STOP_TIME_DEFAULT,
        r_min=TRAJECTORY_R_MIN_DEFAULT,
    )


def test_clamp_spiral_radius_leaves_outward_drift_untouched() -> None:
    ov = {
        "trajectory_lump0_R0": 1.5,
        "trajectory_lump0_v_rad": 0.2,
        "stop_time": TRAJECTORY_STOP_TIME_DEFAULT,
    }
    _clamp_trajectory_spiral_radius(ov)
    assert ov["trajectory_lump0_v_rad"] == 0.2


def test_vector_to_overrides_applies_spiral_radius_after_speed() -> None:
    """Speed conversion now happens before the spiral-radius floor."""
    dims = [
        SearchDimension("trajectory_lump0_R0", 1.5, 8.0, 1.5),
        SearchDimension("trajectory_lump0_omega_rot", -1.0, 1.0, -0.1),
        SearchDimension("trajectory_lump0_v_rad", -1.0, 1.0, -0.9),
    ]
    overrides = _vector_to_overrides([1.5, 0.0, -1.0], dims, {})
    # v_rad_norm = -1 -> v_rad_physical = -0.3, but inward floor is -0.0875
    assert overrides["trajectory_lump0_v_rad"] >= _min_spiral_v_rad(
        1.5,
        a_breath=0.0,
        stop_time=TRAJECTORY_STOP_TIME_DEFAULT,
        r_min=TRAJECTORY_R_MIN_DEFAULT,
    )
