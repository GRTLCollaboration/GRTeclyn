import math

from grteclyn_wrapper.search.optimize import (
    build_grtresna_config,
    build_search_space,
    grtresna_ring_search_space,
)


def test_shell_ansatz_search_space_is_low_dimensional() -> None:
    space = build_search_space(grtresna=True, grtresna_lumps=5, grtresna_ansatz="shell")

    assert len(space) == 16
    assert {dim.param_key for dim in space} >= {
        "grtresna_shell_amp",
        "grtresna_shell_radius",
        "grtresna_shell_axis_theta",
        "grtresna_shell_axis_phi",
        "grtresna_shell_toroidal_velocity",
        "grtresna_shell_poloidal_velocity",
        "grtresna_shell_exotic_fraction",
    }


def test_shell_ansatz_stays_compact_relative_to_free() -> None:
    shell = build_search_space(grtresna=True, grtresna_lumps=5, grtresna_ansatz="shell")
    free = build_search_space(grtresna=True, grtresna_lumps=5, grtresna_ansatz="free")
    # The whole point: full-sphere coverage without the 55D free blow-up.
    assert len(shell) < len(free)
    assert len(free) == 55


def test_shell_ansatz_expands_to_lumps_on_full_sphere() -> None:
    overrides = {
        "grtresna_shell_lumps": 6,
        "grtresna_shell_amp": 0.15,
        "grtresna_shell_width": 3.0,
        "grtresna_shell_radius": 4.0,
        "grtresna_shell_thickness": 0.0,
        "grtresna_shell_axis_theta": 0.5 * math.pi,
        "grtresna_shell_axis_phi": 0.0,
        "grtresna_shell_toroidal_velocity": 0.4,
        "grtresna_shell_poloidal_velocity": 0.0,
        "grtresna_shell_radial_velocity": 0.0,
        "grtresna_shell_exotic_fraction": 0.5,
        "grtresna_shell_exotic_phase": 0.0,
        "grtresna_shell_mode": 1.0,
    }

    cfg = build_grtresna_config(overrides)

    assert len(cfg.lumps) == 6
    assert sum(lump["exotic"] for lump in cfg.lumps) == 3
    assert cfg.maximal_slicing

    # Centers should cover the full sphere, not lie on a single plane: at least
    # one lump must have a non-trivial component along EVERY axis (the planar
    # ring would pin one coordinate near zero for all lumps).
    for axis in range(3):
        spread = max(abs(lump["center"][axis]) for lump in cfg.lumps)
        assert spread > 0.5, f"lumps collapse along axis {axis}: shell is not 3D"

    # All centers sit (approximately) on the requested radius shell.
    for lump in cfg.lumps:
        cx, cy, cz = lump["center"]
        r = math.sqrt(cx * cx + cy * cy + cz * cz)
        assert abs(r - 4.0) < 1e-6


def test_shell_toroidal_current_carries_angular_momentum_about_axis() -> None:
    # Pure toroidal flow about the z-axis must give a net L_z (cross r x v),
    # i.e. the gravitomagnetic "motor" -- the physics the ring also has but now
    # reachable for an arbitrary axis.
    overrides = {
        "grtresna_shell_lumps": 8,
        "grtresna_shell_radius": 4.0,
        "grtresna_shell_thickness": 0.0,
        "grtresna_shell_axis_theta": 0.0,  # axis = +z
        "grtresna_shell_axis_phi": 0.0,
        "grtresna_shell_toroidal_velocity": 0.5,
        "grtresna_shell_poloidal_velocity": 0.0,
        "grtresna_shell_radial_velocity": 0.0,
        "grtresna_shell_exotic_fraction": 0.0,
    }
    cfg = build_grtresna_config(overrides)
    lz = 0.0
    for lump in cfg.lumps:
        cx, cy, _ = lump["center"]
        vx, vy, _ = lump["velocity"]
        lz += cx * vy - cy * vx
    assert lz > 0.5, "toroidal current produced no net angular momentum about the axis"


def test_shell_poloidal_current_is_not_reachable_by_planar_ring() -> None:
    # Poloidal (over-the-pole) flow puts vertical velocity on lumps that are
    # themselves off the equatorial plane -- a current topology the planar ring
    # cannot represent. Verify the expansion produces out-of-plane velocity.
    overrides = {
        "grtresna_shell_lumps": 8,
        "grtresna_shell_radius": 4.0,
        "grtresna_shell_axis_theta": 0.0,  # axis = +z
        "grtresna_shell_poloidal_velocity": 0.5,
        "grtresna_shell_toroidal_velocity": 0.0,
        "grtresna_shell_radial_velocity": 0.0,
        "grtresna_shell_exotic_fraction": 0.0,
    }
    cfg = build_grtresna_config(overrides)
    max_vz = max(abs(lump["velocity"][2]) for lump in cfg.lumps)
    assert max_vz > 0.1, "poloidal current produced no vertical (over-pole) flow"


def test_ring_ansatz_is_unchanged() -> None:
    # Regression guard: adding 'shell' must not perturb the established ring.
    assert len(grtresna_ring_search_space()) == 14


def test_shell_compact_profile_caps_width() -> None:
    compact = build_search_space(
        grtresna=True, grtresna_lumps=5, grtresna_ansatz="shell",
        grtresna_shell_profile="compact",
    )
    middle = build_search_space(
        grtresna=True, grtresna_lumps=5, grtresna_ansatz="shell",
        grtresna_shell_profile="middle",
    )
    compact_width = next(d for d in compact if d.param_key == "grtresna_shell_width")
    middle_width = next(d for d in middle if d.param_key == "grtresna_shell_width")

    assert compact_width.upper == 3.0
    assert middle_width.upper == 4.0
