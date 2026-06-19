"""Tests for boson splash search space."""

from __future__ import annotations

from grteclyn_wrapper.grtresna.matter_models import GRTRESNA_COMPLEX_SCALAR_MODEL
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig
from grteclyn_wrapper.search.optimize import build_grtresna_config, build_search_space


def test_build_search_space_splash_ansatz() -> None:
    space = build_search_space(
        grtresna=True,
        grtresna_ansatz="splash",
        grtresna_matter_sector="boson_star",
    )
    keys = {d.param_key for d in space}
    assert len(space) == 7
    assert "grtresna_bs_omega" in keys


def test_splash_omega_unpinned() -> None:
    space = build_search_space(grtresna=True, grtresna_ansatz="splash")
    omega = next(d for d in space if d.param_key == "grtresna_bs_omega")
    assert omega.lower == 0.05
    assert omega.upper == 0.4
    assert omega.upper > omega.lower


def test_splash_profile_width_compact() -> None:
    space = build_search_space(grtresna=True, grtresna_ansatz="splash")
    width = next(d for d in space if d.param_key == "grtresna_bs_profile_width")
    assert width.lower == 2.0
    assert width.upper == 8.0
    assert width.initial == 3.5


def test_splash_phi_c_capped_for_convergence() -> None:
    space = build_search_space(grtresna=True, grtresna_ansatz="splash")
    phi_c = next(d for d in space if d.param_key == "grtresna_bs_phi_c")
    assert phi_c.upper == 0.12


def test_boson_star_omega_remains_pinned() -> None:
    space = build_search_space(grtresna=True, grtresna_ansatz="boson_star")
    omega = next(d for d in space if d.param_key == "grtresna_bs_omega")
    assert omega.lower == omega.upper == 0.0


def test_boson_star_phi_c_unchanged_for_ftl_campaigns() -> None:
    space = build_search_space(
        grtresna=True,
        grtresna_ansatz="boson_star",
        grtresna_matter_sector="boson_star",
    )
    phi_c = next(d for d in space if d.param_key == "grtresna_bs_phi_c")
    assert phi_c.upper == 0.15


def test_shell_ansatz_unaffected_by_splash() -> None:
    space = build_search_space(grtresna=True, grtresna_ansatz="shell")
    keys = {d.param_key for d in space}
    assert "grtresna_bs_omega" not in keys
    assert "grtresna_shell_amp" in keys


def test_build_grtresna_config_sets_bs_omega() -> None:
    cfg = build_grtresna_config(
        {
            "grtresna_matter_model": GRTRESNA_COMPLEX_SCALAR_MODEL,
            "grtresna_scalar_mass": 0.1,
            "grtresna_bs_phi_c": 0.08,
            "grtresna_bs_profile_width": 10.0,
            "grtresna_bs_omega": 0.35,
        },
        GRTresnaConfig(),
    )
    assert cfg.bs_omega == 0.35
