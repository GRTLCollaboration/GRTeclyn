import argparse
import math

from grteclyn_wrapper.__main__ import build_parser
from grteclyn_wrapper.cli.grtresna_context import build_grtresna_search_context
from grteclyn_wrapper.search.optimize import (
    build_grtresna_config,
    build_search_space,
    grtresna_ring_search_space,
)


def test_shell_ansatz_search_space_is_low_dimensional() -> None:
    space = build_search_space(grtresna=True, grtresna_lumps=5, grtresna_ansatz="shell")

    assert len(space) == 23
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
        "grtresna_scalar_lambda",
        "grtresna_matter_layout",
        "grtresna_shell_static",
        "grtresna_shell_profile_fraction",
        "grtresna_shell_profile_phase",
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


def _shell_base_overrides() -> dict:
    return {
        "grtresna_shell_lumps": 8,
        "grtresna_shell_amp": 0.13,
        "grtresna_shell_width": 2.4,
        "grtresna_shell_radius": 4.0,
        "grtresna_shell_thickness": 0.5,
        "grtresna_shell_axis_theta": 0.5 * math.pi,
        "grtresna_shell_axis_phi": 0.0,
        "grtresna_shell_toroidal_velocity": 0.2,
        "grtresna_shell_poloidal_velocity": 0.0,
        "grtresna_shell_radial_velocity": 0.0,
        "grtresna_shell_exotic_fraction": 0.0,
    }


def test_scalar_lambda_flows_into_grtresna_config() -> None:
    base = _shell_base_overrides()
    cfg_default = build_grtresna_config(dict(base))
    assert cfg_default.scalar_lambda == 0.0

    cfg_lam = build_grtresna_config({**base, "grtresna_scalar_lambda": 0.05})
    assert cfg_lam.scalar_lambda == 0.05


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


def _axis_projection(centers: list[tuple[float, float, float]], axis_idx: int = 2) -> float:
    return max(c[axis_idx] for c in centers) - min(c[axis_idx] for c in centers)


def test_matter_layout_channel_is_elongated_along_axis() -> None:
    overrides = {**_shell_base_overrides(), "grtresna_matter_layout": 1.0}
    cfg = build_grtresna_config(overrides)
    centers = [tuple(lump["center"]) for lump in cfg.lumps]
    # axis_theta=pi/2, axis_phi=0 => polar axis is +x
    assert _axis_projection(centers, 0) > 5.0


def test_matter_layout_bipolar_splits_along_axis() -> None:
    overrides = {**_shell_base_overrides(), "grtresna_matter_layout": 2.0}
    cfg = build_grtresna_config(overrides)
    # axis_theta=pi/2, axis_phi=0 => polar axis is +x
    projections = [lump["center"][0] for lump in cfg.lumps]
    assert max(projections) > 0.5 and min(projections) < -0.5


def test_pin_dimension_removes_key_from_search_space() -> None:
    args = argparse.Namespace(
        nonspherical=False,
        grtresna=True,
        grtresna_lumps=5,
        grtresna_ansatz="shell",
        grtresna_shell_profile="compact",
        pin_dimension=["grtresna_matter_layout=2"],
        grtresna_full_z=False,
        grtresna_evolution_l_full=64.0,
        grtresna_evolution_n_full=64,
        grtresna_domain_l=128.0,
        grtresna_domain_nx=64,
        grtresna_domain_ny=64,
        grtresna_domain_nz=None,
        grtresna_gridinit_nx=64,
        grtresna_gridinit_ny=64,
        grtresna_gridinit_nz=64,
        grtresna_ranks=8,
        grtresna_iterations=50,
        grtresna_timeout=3600,
        grtresna_max_level=3,
        grtresna_refine_threshold=0.5,
        grtresna_regrid_radius=0.0,
        grtresna_coefficient_average_type="harmonic",
        grtresna_psi_relaxation=1.0,
        grtresna_psi_floor=-1.0,
        grtresna_jacobian_cap=-1.0,
        grtresna_keep_source=False,
        grtresna_max_ham_pct=5.0,
        grtresna_max_mom_pct=5.0,
        grtresna_solved_ftl_gate=False,
        grtresna_solved_ftl_min=0.0,
        grtresna_solved_ftl_max=1.0,
        grtresna_speed_mode=None,
        grtresna_speed_factor=1.0,
    )
    ctx = build_grtresna_search_context(args, {"grtresna_shell_lumps": 5})
    pinned_keys = {d.param_key for d in ctx.search_space}
    assert "grtresna_matter_layout" not in pinned_keys
    assert ctx.base_overrides["grtresna_matter_layout"] == 2.0
    cfg = build_grtresna_config(ctx.base_overrides)
    projections = [lump["center"][0] for lump in cfg.lumps]
    assert max(projections) > 0.5 and min(projections) < -0.5


def test_matter_layout_ring_is_planar() -> None:
    overrides = {**_shell_base_overrides(), "grtresna_matter_layout": 3.0}
    cfg = build_grtresna_config(overrides)
    centers = [tuple(lump["center"]) for lump in cfg.lumps]
    # Ring lies in the plane orthogonal to the polar axis (+x).
    x_spread = _axis_projection(centers, 0)
    yz_radii = [math.hypot(c[1], c[2]) for c in centers]
    assert x_spread < 1.5
    assert min(yz_radii) > 3.0
    assert max(yz_radii) < 5.0
    assert _axis_projection(centers, 1) > 2.0
    assert _axis_projection(centers, 2) > 2.0


def test_matter_layout_sphere_covers_all_axes() -> None:
    overrides = {**_shell_base_overrides(), "grtresna_matter_layout": 0.0}
    cfg = build_grtresna_config(overrides)
    for axis in range(3):
        spread = max(abs(lump["center"][axis]) for lump in cfg.lumps)
        assert spread > 0.5


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


def test_matter_layout_cloud_is_bounded_and_volumetric() -> None:
    # The cloud scatters lumps irregularly inside a ball -- unlike the sphere/
    # ring layouts it does NOT pin every lump to a single radius shell.
    overrides = {**_shell_base_overrides(), "grtresna_matter_layout": 4.0}
    cfg = build_grtresna_config(overrides)
    centers = [tuple(lump["center"]) for lump in cfg.lumps]
    span = 4.0 + 0.5  # radius + thickness from _shell_base_overrides
    radii = [math.sqrt(c[0] ** 2 + c[1] ** 2 + c[2] ** 2) for c in centers]
    # Bounded inside the ball.
    for r in radii:
        assert r <= span + 1e-9
    # Volumetric: the lumps span a range of radii (not a thin shell).
    assert max(radii) - min(radii) > 0.5
    # Still genuinely 3D (every axis populated).
    for axis in range(3):
        assert max(abs(c[axis]) for c in centers) > 0.3


def test_matter_layout_cloud_is_deterministic() -> None:
    overrides = {**_shell_base_overrides(), "grtresna_matter_layout": 4.0}
    a = build_grtresna_config(dict(overrides))
    b = build_grtresna_config(dict(overrides))
    assert [tuple(l["center"]) for l in a.lumps] == [tuple(l["center"]) for l in b.lumps]


def test_profile_fraction_assigns_top_hat_to_subset() -> None:
    base = {**_shell_base_overrides(), "grtresna_shell_lumps": 6}
    # Default (no profile fraction) keeps every lump Gaussian -- v13 behaviour.
    cfg0 = build_grtresna_config(dict(base))
    assert all(int(l.get("profile", 0)) == 0 for l in cfg0.lumps)
    # Half the lumps switch to the smoothed top-hat profile.
    cfg_half = build_grtresna_config({**base, "grtresna_shell_profile_fraction": 0.5})
    assert sum(int(l["profile"]) for l in cfg_half.lumps) == 3
    # All lumps switch.
    cfg_all = build_grtresna_config({**base, "grtresna_shell_profile_fraction": 1.0})
    assert all(int(l["profile"]) == 1 for l in cfg_all.lumps)


def test_top_hat_envelope_has_flatter_core_than_gaussian() -> None:
    from grteclyn_wrapper.grtresna.fields.lump import lump_phi_at

    width = 3.0
    gauss = {"amp": 0.2, "width": width, "center": (0.0, 0.0, 0.0), "mode": 0, "profile": 0}
    ball = {**gauss, "profile": 1}
    mid = (0.5 * width, 0.0, 0.0)
    g_ratio = lump_phi_at(gauss, mid) / lump_phi_at(gauss, (0.0, 0.0, 0.0))
    b_ratio = lump_phi_at(ball, mid) / lump_phi_at(ball, (0.0, 0.0, 0.0))
    # The top-hat holds near-peak density across the core; the Gaussian has
    # already decayed, so the top-hat is flatter (closer to 1) at mid-core.
    assert b_ratio > g_ratio


def test_solver_serializes_lump_profile() -> None:
    from grteclyn_wrapper.grtresna.solver import GRTresnaConfig, _lump_lines

    cfg = GRTresnaConfig()
    cfg.lumps = [
        {"amp": 0.1, "width": 3.0, "center": (0.0, 0.0, 0.0), "velocity": (0.0, 0.0, 0.0),
         "omega": 0.0, "mode": 0, "exotic": 0, "profile": 1},
        {"amp": 0.1, "width": 3.0, "center": (1.0, 0.0, 0.0), "velocity": (0.0, 0.0, 0.0),
         "omega": 0.0, "mode": 0, "exotic": 0, "profile": 0},
    ]
    text = "\n".join(_lump_lines(cfg))
    assert "lump0_profile = 1" in text
    assert "lump1_profile = 0" in text
