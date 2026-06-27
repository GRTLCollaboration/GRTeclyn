"""Tests for complex scalar matter wiring."""

from __future__ import annotations

import json

from grteclyn_wrapper.grtresna.boson_star_fields import rename_complex_scalar_components
from grteclyn_wrapper.grtresna.matter_wiring import (
    GRTRESNA_COMPLEX_SCALAR_MODEL,
    evolution_overrides_from_bicomplex_scalar,
    evolution_overrides_from_complex_scalar,
    evolution_overrides_from_metadata,
    plot_vars_for_bicomplex_scalar,
    plot_vars_for_complex_scalar,
    read_matter_metadata,
    write_matter_metadata,
)
from grteclyn_wrapper.grtresna.matter_models import GRTRESNA_BICOMPLEX_SCALAR_MODEL


def test_rename_complex_scalar_components() -> None:
    renamed = rename_complex_scalar_components(
        ["chi", "phi_re", "Pi_re", "phi_im", "Pi_im"]
    )
    assert renamed == ["chi", "phi", "Pi", "phi2", "Pi2"]


def test_evolution_overrides_for_complex_scalar() -> None:
    overrides = evolution_overrides_from_complex_scalar(mass=0.1, lam=0.0, mu=5333.0)
    assert overrides["recipe_matter_model"] == GRTRESNA_COMPLEX_SCALAR_MODEL
    assert overrides["recipe_scalar_mass"] == 0.1
    assert overrides["recipe_scalar_lambda"] == 0.0
    assert overrides["recipe_scalar_mu"] == 5333.0
    assert "phi_lump0" in overrides["amr.plot_vars"]


def test_read_matter_metadata_includes_scalar_lambda(tmp_path) -> None:
    payload = {
        "matter_model": "grtresna_independent_scalars",
        "num_scalar_fields": 1,
        "scalar_field_signs": [1],
        "scalar_mass": 0.1,
        "scalar_lambda": 0.02,
        "scalar_mu": 85333.0,
        "lump_count": 1,
    }
    path = tmp_path / "matter.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    meta = read_matter_metadata(path)
    assert meta.scalar_lambda == 0.02
    assert meta.scalar_mass == 0.1
    assert meta.scalar_mu == 85333.0


def test_bicomplex_metadata_round_trip_includes_scalar_mu(tmp_path) -> None:
    from grteclyn_wrapper.grtresna.solver import GRTresnaConfig

    cfg = GRTresnaConfig(
        matter_model=GRTRESNA_BICOMPLEX_SCALAR_MODEL,
        scalar_mass=1.0,
        scalar_lambda=640.0,
        scalar_mu=85333.0,
        bs_omega=0.4,
        lumps=[
            {"amp": 0.15, "width": 1.09, "center": (0.0, 0.0, 0.0), "exotic": 0},
            {"amp": 0.15, "width": 1.09, "center": (1.0, 0.0, 0.0), "exotic": 1},
        ],
    )
    path = write_matter_metadata(tmp_path / "matter.json", cfg)
    meta = read_matter_metadata(path)
    assert meta.scalar_mu == 85333.0
    assert meta.scalar_lambda == 640.0
    evo = evolution_overrides_from_metadata(meta)
    assert evo["recipe_scalar_mu"] == 85333.0
    assert evo["recipe_matter_model"] == GRTRESNA_BICOMPLEX_SCALAR_MODEL


def test_plot_vars_for_complex_scalar() -> None:
    names = plot_vars_for_complex_scalar()
    assert names == (
        "chi", "h11", "h12", "h13", "h22", "h23", "h33", "K",
        "lapse", "shift1", "shift2", "shift3",
        "phi", "Pi", "phi_lump0", "Pi_lump0",
    )


def test_plot_vars_for_bicomplex_scalar_includes_phantom_channels() -> None:
    names = plot_vars_for_bicomplex_scalar()
    # Canonical Im (lump0) plus phantom Re/Pi (lump1) and phantom Im (lump2).
    for nm in ("phi_lump0", "Pi_lump0", "phi_lump1", "Pi_lump1", "phi_lump2", "Pi_lump2"):
        assert nm in names


def test_evolution_overrides_bicomplex_uses_phantom_plot_vars() -> None:
    overrides = evolution_overrides_from_bicomplex_scalar(
        mass=0.15, lam=0.0, field_signs=(1, -1, 1)
    )
    assert overrides["recipe_matter_model"] == GRTRESNA_BICOMPLEX_SCALAR_MODEL
    assert overrides["recipe_scalar_field_signs"] == "1 -1 1"
    assert "phi_lump1" in overrides["amr.plot_vars"]


def test_evolution_overrides_enable_pd_controller() -> None:
    # bs_omega drives the closed-loop PD trap pump controller: the pump frequency
    # carries the boson phase velocity and the kp/kd gains switch on the feedback
    # controller (kp > 0).
    overrides = evolution_overrides_from_complex_scalar(
        mass=0.3, lam=0.0, bs_omega=0.25
    )
    assert overrides["trajectory_pump_frequency"] == 0.25
    assert overrides["rl_pump_kp"] > 0.0
    assert overrides["rl_pump_kd"] > 0.0


def test_evolution_overrides_bicomplex_enable_pd_controller() -> None:
    overrides = evolution_overrides_from_bicomplex_scalar(
        mass=0.3, lam=0.0, field_signs=(1, -1, 1), bs_omega=0.25
    )
    assert overrides["trajectory_pump_frequency"] == 0.25
    assert overrides["rl_pump_kp"] > 0.0
    assert overrides["rl_pump_kd"] > 0.0


def test_evolution_overrides_phantom_sign() -> None:
    overrides = evolution_overrides_from_complex_scalar(mass=0.1, lam=0.0, sign=-1.0)
    assert overrides["recipe_scalar_sign"] == -1.0


def test_evolution_overrides_canonical_omits_sign() -> None:
    overrides = evolution_overrides_from_complex_scalar(mass=0.1, lam=0.0, sign=1.0)
    assert "recipe_scalar_sign" not in overrides
