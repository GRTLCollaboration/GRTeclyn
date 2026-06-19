"""Tests for bosonic shell: shell geometry + complex scalar matter."""

from __future__ import annotations

import math

from grteclyn_wrapper.__main__ import build_parser
from grteclyn_wrapper.cli.grtresna_context import build_grtresna_search_context
from grteclyn_wrapper.grtresna.matter_models import (
    GRTRESNA_COMPLEX_SCALAR_MODEL,
    GRTRESNA_EXAMPLE_BOSON_STAR_BH,
    matter_selection_base_overrides,
    resolve_matter_selection,
)
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig, _lump_lines, write_grtresna_params
from grteclyn_wrapper.search.optimize import build_grtresna_config, build_search_space


def _boson_shell_overrides() -> dict:
    return {
        **matter_selection_base_overrides(
            resolve_matter_selection("boson_star", "canonical")
        ),
        "grtresna_shell_lumps": 5,
        "grtresna_shell_amp": 0.15,
        "grtresna_shell_width": 3.0,
        "grtresna_shell_radius": 4.0,
        "grtresna_shell_thickness": 0.5,
        "grtresna_shell_axis_theta": 0.5 * math.pi,
        "grtresna_shell_axis_phi": 0.0,
        "grtresna_shell_toroidal_velocity": 0.0,
        "grtresna_shell_static": 1.0,
        "grtresna_scalar_mass": 0.1,
        "grtresna_scalar_lambda": 0.0,
        "grtresna_bs_omega": 0.15,
    }


def test_boson_shell_search_space_has_geometry_not_exotic() -> None:
    space = build_search_space(
        grtresna=True,
        grtresna_ansatz="shell",
        grtresna_matter_sector="boson_star",
    )
    keys = {d.param_key for d in space}
    assert "grtresna_shell_amp" in keys
    assert "grtresna_shell_radius" in keys
    assert "grtresna_bs_omega" in keys
    assert "grtresna_shell_exotic_fraction" not in keys
    assert "grtresna_shell_exotic_phase" not in keys
    assert "grtresna_bs_phi_c" not in keys
    amp = next(d for d in space if d.param_key == "grtresna_shell_amp")
    mass = next(d for d in space if d.param_key == "grtresna_scalar_mass")
    assert amp.upper <= 0.12
    assert mass.upper <= 0.35
    assert "grtresna_shell_toroidal_velocity" not in keys


def test_boson_shell_static_zeros_lump_velocities() -> None:
    overrides = {**_boson_shell_overrides(), "grtresna_shell_static": 1.0}
    cfg = build_grtresna_config(overrides, GRTresnaConfig())
    for lump in cfg.lumps:
        assert lump["velocity"] == (0.0, 0.0, 0.0)
        assert lump["omega"] == 0.0


def test_boson_shell_params_include_lump_kinematics_keys(tmp_path) -> None:
    cfg = build_grtresna_config(_boson_shell_overrides(), GRTresnaConfig())
    path = tmp_path / "params.txt"
    write_grtresna_params(cfg, path)
    text = path.read_text(encoding="utf-8")
    assert "lump0_velocity" in text
    assert "lump0_omega" in text
    assert "lump0_mode" in text


def test_boson_shell_expands_to_canonical_lumps() -> None:
    cfg = build_grtresna_config(_boson_shell_overrides(), GRTresnaConfig())
    assert cfg.matter_model == GRTRESNA_COMPLEX_SCALAR_MODEL
    assert cfg.example == GRTRESNA_EXAMPLE_BOSON_STAR_BH
    assert len(cfg.lumps) == 5
    assert all(lump["exotic"] == 0 for lump in cfg.lumps)
    for lump in cfg.lumps:
        cx, cy, cz = lump["center"]
        r = math.sqrt(cx * cx + cy * cy + cz * cz)
        assert 3.0 < r < 5.5


def test_boson_shell_writes_num_lumps_params(tmp_path) -> None:
    cfg = build_grtresna_config(_boson_shell_overrides(), GRTresnaConfig())
    path = tmp_path / "params.txt"
    write_grtresna_params(cfg, path)
    text = path.read_text(encoding="utf-8")
    assert "num_lumps = 5" in text
    assert "lump0_amp" in text
    assert "bs_omega" in text
    assert "scalar_sign = 1" in text


def test_boson_shell_lump_lines_match_solver_keys() -> None:
    cfg = build_grtresna_config(_boson_shell_overrides(), GRTresnaConfig())
    lines = "\n".join(_lump_lines(cfg))
    assert "lump0_center" in lines
    assert "lump4_amp" in lines


def test_grtresna_context_injects_shell_lumps_for_boson() -> None:
    parser = build_parser()
    args = parser.parse_args(
        [
            "qd",
            "--grtresna",
            "--grtresna-ansatz",
            "shell",
            "--grtresna-matter-sector",
            "boson_star",
            "--grtresna-lumps",
            "6",
        ]
    )
    ctx = build_grtresna_search_context(args, {})
    assert ctx.base_overrides["grtresna_shell_lumps"] == 6
    assert ctx.base_overrides["grtresna_matter_sector"] == "boson_star"
