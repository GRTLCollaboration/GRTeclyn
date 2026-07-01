"""Tests for trajectory constraints: sub-luminal speed cap + retrograde enforcement.

A trajectory lump is dragged along its orbit by a co-moving trap whose target
centre advances at v_t = R0 * |omega_rot| (geometric units, c = 1).  Superluminal
or strongly relativistic v_t cannot be followed by any soliton, so the decode
path clamps omega_rot per lump.

The optional retrograde constraint (``trajectory_retrograde_only=1``) negates
prograde omega_rot values — HQ validation showed counter-rotation is a
false-positive generator and all confirmed FTL configs are all-retrograde.
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


def test_clamp_caps_superluminal_eval122_lumps() -> None:
    # Eval-122 champion lumps: v_t up to ~6c.
    ov = {
        "trajectory_lump0_R0": 7.12,
        "trajectory_lump0_omega_rot": -0.851,
        "trajectory_lump1_R0": 3.86,
        "trajectory_lump1_omega_rot": -0.909,
    }
    _clamp_trajectory_speed(ov)
    for k, r0 in ((0, 7.12), (1, 3.86)):
        v_t = abs(ov[f"trajectory_lump{k}_omega_rot"]) * r0
        assert v_t <= TRAJECTORY_V_MAX_DEFAULT + 1e-9
    # Sign (orbit direction) is preserved.
    assert ov["trajectory_lump0_omega_rot"] < 0.0
    assert ov["trajectory_lump1_omega_rot"] < 0.0


def test_clamp_leaves_subluminal_lumps_untouched() -> None:
    ov = {"trajectory_lump0_R0": 4.0, "trajectory_lump0_omega_rot": 0.05}  # v_t = 0.2
    _clamp_trajectory_speed(ov)
    assert ov["trajectory_lump0_omega_rot"] == 0.05


def test_clamp_respects_custom_v_max() -> None:
    ov = {
        "trajectory_lump0_R0": 5.0,
        "trajectory_lump0_omega_rot": 0.5,  # v_t = 2.5
        "trajectory_v_max": 0.1,
    }
    _clamp_trajectory_speed(ov)
    assert math.isclose(abs(ov["trajectory_lump0_omega_rot"]) * 5.0, 0.1, rel_tol=1e-9)


def test_clamp_disabled_when_v_max_non_positive() -> None:
    ov = {
        "trajectory_lump0_R0": 8.0,
        "trajectory_lump0_omega_rot": 1.0,
        "trajectory_v_max": 0.0,
    }
    _clamp_trajectory_speed(ov)
    assert ov["trajectory_lump0_omega_rot"] == 1.0


def test_clamp_noop_for_non_trajectory_overrides() -> None:
    ov = {"grtresna_lump0_amp": 0.1, "grtresna_scalar_mass": 0.2}
    before = dict(ov)
    _clamp_trajectory_speed(ov)
    assert ov == before


def test_vector_to_overrides_applies_cap() -> None:
    dims = [
        SearchDimension("trajectory_lump0_R0", 1.5, 8.0, 5.0),
        SearchDimension("trajectory_lump0_omega_rot", -1.0, 1.0, 0.0),
    ]
    overrides = _vector_to_overrides([8.0, 1.0], dims, {})
    v_t = abs(overrides["trajectory_lump0_omega_rot"]) * overrides["trajectory_lump0_R0"]
    assert v_t <= TRAJECTORY_V_MAX_DEFAULT + 1e-9


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
    # Should be negated AND speed-capped.
    assert overrides["trajectory_lump0_omega_rot"] <= 0.0
    v_t = abs(overrides["trajectory_lump0_omega_rot"]) * overrides["trajectory_lump0_R0"]
    assert v_t <= TRAJECTORY_V_MAX_DEFAULT + 1e-9


# ---------------------------------------------------------------------------
# Spiral (radial drift) speed cap
# ---------------------------------------------------------------------------

def test_clamp_caps_combined_tangential_and_radial_speed() -> None:
    ov = {
        "trajectory_lump0_R0": 4.0,
        "trajectory_lump0_omega_rot": 0.1,  # v_t = 0.4
        "trajectory_lump0_v_rad": 0.4,      # total = 0.566 > 0.3
    }
    _clamp_trajectory_speed(ov)
    v_t = abs(ov["trajectory_lump0_omega_rot"]) * ov["trajectory_lump0_R0"]
    v_rad = abs(ov["trajectory_lump0_v_rad"])
    v_total = math.sqrt(v_t * v_t + v_rad * v_rad)
    assert v_total <= TRAJECTORY_V_MAX_DEFAULT + 1e-9
    # Both components should be scaled by the same factor.
    assert math.isclose(v_t / v_rad, 1.0, rel_tol=1e-9)


def test_clamp_leaves_subluminal_spiral_untouched() -> None:
    ov = {
        "trajectory_lump0_R0": 4.0,
        "trajectory_lump0_omega_rot": 0.05,  # v_t = 0.2
        "trajectory_lump0_v_rad": 0.05,      # total = 0.206 < 0.3
    }
    _clamp_trajectory_speed(ov)
    assert ov["trajectory_lump0_omega_rot"] == 0.05
    assert ov["trajectory_lump0_v_rad"] == 0.05


def test_clamp_scales_omega_only_when_v_rad_missing() -> None:
    """When v_rad is absent it is treated as 0; only omega_rot is scaled."""
    ov = {
        "trajectory_lump0_R0": 4.0,
        "trajectory_lump0_omega_rot": 0.1,  # v_t = 0.4
    }
    _clamp_trajectory_speed(ov)
    v_t = abs(ov["trajectory_lump0_omega_rot"]) * ov["trajectory_lump0_R0"]
    assert v_t <= TRAJECTORY_V_MAX_DEFAULT + 1e-9
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


def test_vector_to_overrides_applies_spiral_radius_before_speed() -> None:
    dims = [
        SearchDimension("trajectory_lump0_R0", 1.5, 8.0, 1.5),
        SearchDimension("trajectory_lump0_omega_rot", -1.0, 1.0, -0.1),
        SearchDimension("trajectory_lump0_v_rad", -0.3, 0.3, -0.3),
    ]
    overrides = _vector_to_overrides(
        [1.5, -0.1, -0.3],
        dims,
        {"stop_time": TRAJECTORY_STOP_TIME_DEFAULT},
    )
    assert overrides["trajectory_lump0_v_rad"] > -0.3
    v_t = abs(overrides["trajectory_lump0_omega_rot"]) * overrides["trajectory_lump0_R0"]
    v_rad = abs(overrides["trajectory_lump0_v_rad"])
    v_total = math.sqrt(v_t * v_t + v_rad * v_rad)
    assert v_total <= TRAJECTORY_V_MAX_DEFAULT + 1e-9
