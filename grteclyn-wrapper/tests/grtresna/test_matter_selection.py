"""Tests for selectable matter sector × coupling."""

from __future__ import annotations

from grteclyn_wrapper.grtresna.matter_models import (
    GRTRESNA_COMPLEX_SCALAR_MODEL,
    GRTRESNA_INDEPENDENT_MATTER_MODEL,
    MATTER_COUPLING_EXOTIC,
    MATTER_SECTOR_BOSON_STAR,
    MATTER_SECTOR_SCALAR,
    matter_selection_base_overrides,
    resolve_matter_selection,
)
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig
from grteclyn_wrapper.search.optimize import build_grtresna_config, build_search_space


def test_matter_selection_four_modes() -> None:
    scalar = resolve_matter_selection("scalar", "canonical")
    assert scalar.model_id == GRTRESNA_INDEPENDENT_MATTER_MODEL
    assert not scalar.is_exotic

    scalar_ex = resolve_matter_selection("scalar", "exotic")
    assert scalar_ex.is_scalar and scalar_ex.is_exotic
    assert scalar_ex.model_id == GRTRESNA_INDEPENDENT_MATTER_MODEL

    boson = resolve_matter_selection("boson_star", "canonical")
    assert boson.model_id == GRTRESNA_COMPLEX_SCALAR_MODEL
    assert boson.scalar_sign == 1

    boson_ex = resolve_matter_selection("boson_star", "phantom")
    assert boson_ex.is_boson_star and boson_ex.is_exotic
    assert boson_ex.scalar_sign == -1


def test_scalar_and_boson_search_spaces_differ() -> None:
    shell = build_search_space(
        grtresna=True,
        grtresna_ansatz="shell",
        grtresna_matter_sector=MATTER_SECTOR_SCALAR,
    )
    boson = build_search_space(
        grtresna=True,
        grtresna_ansatz="shell",
        grtresna_matter_sector=MATTER_SECTOR_BOSON_STAR,
    )
    shell_keys = {d.param_key for d in shell}
    boson_keys = {d.param_key for d in boson}
    assert "grtresna_shell_amp" in shell_keys
    assert "grtresna_shell_amp" in boson_keys
    assert "grtresna_shell_exotic_fraction" in shell_keys
    assert "grtresna_shell_exotic_fraction" not in boson_keys
    assert "grtresna_bs_omega" in boson_keys
    assert "grtresna_bs_omega" not in shell_keys


def test_scalar_exotic_forces_all_lumps_exotic() -> None:
    overrides = {
        **matter_selection_base_overrides(
            resolve_matter_selection("scalar", "exotic")
        ),
        "grtresna_shell_lumps": 5,
        "grtresna_shell_amp": 0.13,
        "grtresna_shell_width": 2.4,
        "grtresna_shell_radius": 4.0,
        "grtresna_shell_thickness": 0.5,
        "grtresna_shell_toroidal_velocity": 0.2,
        "grtresna_shell_exotic_fraction": 0.0,
    }
    cfg = build_grtresna_config(overrides, GRTresnaConfig())
    assert cfg.matter_model == GRTRESNA_INDEPENDENT_MATTER_MODEL
    assert len(cfg.lumps) == 5
    assert all(lump["exotic"] == 1 for lump in cfg.lumps)
    assert cfg.maximal_slicing is True


def test_boson_exotic_base_overrides_set_sign() -> None:
    overrides = matter_selection_base_overrides(
        resolve_matter_selection("boson_star", "exotic")
    )
    assert overrides["grtresna_scalar_sign"] == -1.0
    cfg = build_grtresna_config(
        {**overrides, "grtresna_bs_phi_c": 0.08},
        GRTresnaConfig(),
    )
    assert cfg.scalar_sign == -1
    assert cfg.matter_model == GRTRESNA_COMPLEX_SCALAR_MODEL
