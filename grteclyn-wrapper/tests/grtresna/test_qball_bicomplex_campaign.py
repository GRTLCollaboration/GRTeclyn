"""Campaign wiring: qball_trajectory must select bicomplex matter."""

from __future__ import annotations

from pathlib import Path

from grteclyn_wrapper.grtresna.matter.models import GRTRESNA_BICOMPLEX_SCALAR_MODEL
from grteclyn_wrapper.grtresna.matter.wiring import evolution_overrides_from_config
from grteclyn_wrapper.search.optimize.config import build_grtresna_config

_RUN_SH = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "campaigns"
    / "qball_trajectory"
    / "run.sh"
)


def test_qball_run_sh_pins_bicomplex_matter_model() -> None:
    text = _RUN_SH.read_text(encoding="utf-8")
    assert "grtresna_matter_model=grtresna_bicomplex_scalar" in text
    assert "GRTRESNA_ALLOW_SIGN_MISMATCH" in text


def test_qball_style_overrides_build_bicomplex_with_matched_signs() -> None:
    overrides = {
        "grtresna_matter_model": GRTRESNA_BICOMPLEX_SCALAR_MODEL,
        "grtresna_matter_sector": "boson_star",
        "grtresna_scalar_mass": 1.0,
        "grtresna_scalar_lambda": 640.0,
        "grtresna_scalar_mu": 85333.0,
        "grtresna_bs_omega": 0.8,
        "grtresna_qball_ode_profile": 1,
        "grtresna_qball_equilibrium_amplitude": 1,
        "trajectory_mode": 1,
        "trajectory_num_lumps": 5,
        "trajectory_well_width": 1.667,
        "trajectory_retrograde_only": 1,
    }
    for k in range(5):
        overrides.update(
            {
                f"trajectory_lump{k}_R0": 3.0 + 0.5 * k,
                f"trajectory_lump{k}_omega_rot": -0.05,
                f"trajectory_lump{k}_phase0": 0.2 * k,
                f"trajectory_lump{k}_tilt_theta": 0.1,
                f"trajectory_lump{k}_tilt_phi": 0.0,
                f"trajectory_lump{k}_v_rad": 0.0,
                f"trajectory_lump{k}_well_depth": 0.15,
                f"trajectory_lump{k}_exotic": 0.9 if k else 0.2,
            }
        )
    cfg = build_grtresna_config(overrides)
    assert cfg.matter_model == GRTRESNA_BICOMPLEX_SCALAR_MODEL
    evo = evolution_overrides_from_config(cfg)
    assert evo["recipe_scalar_field_signs"] == "1 -1 -1 -1 -1"
    assert "rl_pump_kp" in evo
