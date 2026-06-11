import math

from grteclyn_wrapper.__main__ import build_parser
from grteclyn_wrapper.search.optimize import (
    build_grtresna_config,
    build_search_space,
    grtresna_ring_search_space,
)


def test_shell_ansatz_search_space_is_low_dimensional() -> None:
    space = build_search_space(grtresna=True, grtresna_lumps=5, grtresna_ansatz="shell")

    assert len(space) == 19
    assert {dim.param_key for dim in space} >= {
        "grtresna_shell_amp",
        "grtresna_shell_radius",
        "grtresna_shell_axis_theta",
        "grtresna_shell_axis_phi",
        "grtresna_shell_toroidal_velocity",
        "grtresna_shell_poloidal_velocity",
        "grtresna_shell_exotic_fraction",
        "grtresna_shift_seed",
        "grtresna_scalar_mass",
        "grtresna_shell_static",
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


def test_shift_seed_flows_into_grtresna_config() -> None:
    base_overrides = {
        "grtresna_shell_lumps": 5,
        "grtresna_shell_radius": 4.0,
        "grtresna_shell_toroidal_velocity": 0.4,
    }
    cfg_default = build_grtresna_config(dict(base_overrides))
    assert cfg_default.shift_seed == 0.0

    cfg_seeded = build_grtresna_config({**base_overrides, "grtresna_shift_seed": 0.3})
    assert cfg_seeded.shift_seed == 0.3


def test_scalar_mass_flows_into_grtresna_config() -> None:
    base_overrides = {
        "grtresna_shell_lumps": 5,
        "grtresna_shell_radius": 4.0,
        "grtresna_shell_toroidal_velocity": 0.4,
    }
    cfg_default = build_grtresna_config(dict(base_overrides))
    assert cfg_default.scalar_mass == 0.1  # GRTresnaConfig default, untouched

    cfg_heavy = build_grtresna_config({**base_overrides, "grtresna_scalar_mass": 1.2})
    assert cfg_heavy.scalar_mass == 1.2


def test_shell_fly_away_velocities_are_capped() -> None:
    space = build_search_space(grtresna=True, grtresna_lumps=5, grtresna_ansatz="shell")
    by_key = {dim.param_key: dim for dim in space}
    # Toroidal (warp motor) keeps its full range; radial/poloidal (net outflow,
    # the "fly away") are capped tighter so matter stays bound.
    assert by_key["grtresna_shell_radial_velocity"].upper == 0.3
    assert by_key["grtresna_shell_poloidal_velocity"].upper == 0.8
    assert by_key["grtresna_shell_toroidal_velocity"].upper == 1.2
    mass = by_key["grtresna_scalar_mass"]
    assert (mass.lower, mass.upper) == (0.3, 1.5)


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


def test_shell_static_toggle_zeroes_all_matter_currents() -> None:
    # With the static flag set, every lump must be momentum-free regardless of
    # the (otherwise active) toroidal/poloidal/radial velocity knobs.
    overrides = {
        "grtresna_shell_lumps": 6,
        "grtresna_shell_radius": 4.0,
        "grtresna_shell_axis_theta": 0.0,
        "grtresna_shell_toroidal_velocity": 0.8,
        "grtresna_shell_poloidal_velocity": 0.5,
        "grtresna_shell_radial_velocity": 0.3,
        "grtresna_shell_omega": 0.4,
        "grtresna_shell_static": 1.0,
    }
    cfg = build_grtresna_config(overrides)
    for lump in cfg.lumps:
        assert lump["velocity"] == (0.0, 0.0, 0.0)
        assert lump["omega"] == 0.0

    # Same knobs with the flag off keep the moving-matter behaviour.
    moving = build_grtresna_config({**overrides, "grtresna_shell_static": 0.0})
    assert max(abs(v) for lump in moving.lumps for v in lump["velocity"]) > 0.1


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


def test_qd_cli_accepts_grtresna_shell_channel_options() -> None:
    parser = build_parser()
    args = parser.parse_args([
        "--consume-plotfiles",
        "qd",
        "--descriptor-mode",
        "channel",
        "--objective-mode",
        "ftl_first",
        "--grtresna",
        "--grtresna-ansatz",
        "shell",
        "--grtresna-shell-profile",
        "compact",
        "--grtresna-full-z",
    ])

    assert args.command == "qd"
    assert args.descriptor_mode == "channel"
    assert args.objective_mode == "ftl_first"
    assert args.grtresna is True
    assert args.grtresna_ansatz == "shell"
    assert args.grtresna_shell_profile == "compact"
    assert args.grtresna_full_z is True
