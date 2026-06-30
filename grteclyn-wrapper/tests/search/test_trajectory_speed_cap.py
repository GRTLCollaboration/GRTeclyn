"""Tests for the sub-luminal / adiabatic trajectory-speed cap.

A trajectory lump is dragged along its orbit by a co-moving trap whose target
centre advances at v_t = R0 * |omega_rot| (geometric units, c = 1).  Superluminal
or strongly relativistic v_t cannot be followed by any soliton, so the decode
path clamps omega_rot per lump.
"""

from __future__ import annotations

import math

from grteclyn_wrapper.search.optimize.candidates import (
    TRAJECTORY_V_MAX_DEFAULT,
    _clamp_trajectory_speed,
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
