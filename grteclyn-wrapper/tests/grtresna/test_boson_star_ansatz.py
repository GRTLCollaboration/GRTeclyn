"""Tests for boson-star MAP-Elites / CMA-ES search integration."""

from __future__ import annotations

import json

from grteclyn_wrapper.grtresna.matter_models import (
    GRTRESNA_COMPLEX_SCALAR_MODEL,
    GRTRESNA_EXAMPLE_BOSON_STAR_BH,
    matter_model_for_ansatz,
)
from grteclyn_wrapper.grtresna.matter_wiring import (
    GRTresnaMatterMetadata,
    evolution_overrides_from_config,
    read_matter_metadata,
    write_matter_metadata,
)
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig
from grteclyn_wrapper.search.optimize import build_grtresna_config, build_search_space


def test_matter_model_for_boson_star_ansatz() -> None:
    assert matter_model_for_ansatz("boson_star") == GRTRESNA_COMPLEX_SCALAR_MODEL
    assert matter_model_for_ansatz("shell") == "grtresna_independent_scalars"


def test_boson_star_search_space_is_seven_dimensional() -> None:
    space = build_search_space(grtresna=True, grtresna_ansatz="boson_star")
    keys = {d.param_key for d in space}
    assert len(space) == 7
    assert "grtresna_bs_phi_c" in keys
    assert "grtresna_scalar_mass" in keys
    assert "grtresna_shell_amp" not in keys


def test_build_grtresna_config_boson_star() -> None:
    base = GRTresnaConfig(
        dphi=0.0,
        dpi=0.0,
        bh1_bare_mass=0.0,
        bh1_spin=(0.0, 0.0, 0.0),
    )
    cfg = build_grtresna_config(
        {
            "grtresna_matter_model": GRTRESNA_COMPLEX_SCALAR_MODEL,
            "grtresna_scalar_mass": 0.12,
            "grtresna_bs_phi_c": 0.06,
            "grtresna_bs_profile_width": 10.0,
            "grtresna_shift_seed": 0.2,
        },
        base,
    )
    assert cfg.matter_model == GRTRESNA_COMPLEX_SCALAR_MODEL
    assert cfg.example == GRTRESNA_EXAMPLE_BOSON_STAR_BH
    assert cfg.lumps == []
    assert cfg.scalar_mass == 0.12
    assert cfg.bs_phi_c == 0.06
    assert cfg.bs_profile_width == 10.0
    assert cfg.shift_seed == 0.2


def test_boson_star_evolution_overrides_match_solver() -> None:
    cfg = build_grtresna_config(
        {
            "grtresna_matter_model": GRTRESNA_COMPLEX_SCALAR_MODEL,
            "grtresna_scalar_mass": 0.1,
            "grtresna_scalar_lambda": 0.02,
            "grtresna_bs_phi_c": 0.08,
        },
    )
    overrides = evolution_overrides_from_config(cfg)
    assert overrides["recipe_matter_model"] == GRTRESNA_COMPLEX_SCALAR_MODEL
    assert overrides["recipe_scalar_mass"] == 0.1
    assert overrides["recipe_scalar_lambda"] == 0.02
    assert "phi2" in overrides["amr.plot_vars"]


def test_boson_star_matter_metadata_roundtrip(tmp_path) -> None:
    cfg = build_grtresna_config(
        {
            "grtresna_matter_model": GRTRESNA_COMPLEX_SCALAR_MODEL,
            "grtresna_scalar_mass": 0.1,
            "grtresna_bs_phi_c": 0.07,
            "grtresna_bs_profile_width": 9.0,
        },
    )
    path = tmp_path / "matter.json"
    write_matter_metadata(path, cfg)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["matter_model"] == GRTRESNA_COMPLEX_SCALAR_MODEL
    assert payload["bs_phi_c"] == 0.07
    meta = read_matter_metadata(path)
    assert isinstance(meta, GRTresnaMatterMetadata)
    assert meta.lump_count == 0
    assert meta.bs_profile_width == 9.0


def test_phantom_boson_star_enables_exotic_safe_solver() -> None:
    cfg = build_grtresna_config(
        {
            "grtresna_matter_model": GRTRESNA_COMPLEX_SCALAR_MODEL,
            "grtresna_scalar_sign": -1.0,
            "grtresna_bs_phi_c": 0.08,
        },
    )
    assert cfg.scalar_sign == -1
    assert cfg.maximal_slicing is True
    assert cfg.psi_relaxation == 0.6
