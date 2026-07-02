"""Matter–geometry consistency helpers for the GRTresna ↔ GRTeclyn bridge."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from grteclyn_wrapper.grtresna.io import write_gridinit
from grteclyn_wrapper.grtresna.fields.lump import (
    lump_phi_at,
    lump_pi_at,
    lump_sign,
    paint_lump_fields_on_grid,
    shift_lump_centers_for_gridinit,
)
from grteclyn_wrapper.core.params import format_param_value, write_params
from grteclyn_wrapper.core.config import resolve_example
from grteclyn_wrapper.grtresna.matter.wiring import (
    GRTRESNA_INDEPENDENT_MATTER_MODEL,
    evolution_overrides_from_config,
    plot_vars_for_independent_scalars,
    write_matter_metadata,
)
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig
from grteclyn_wrapper.projection.postload_gate import (
    PostLoadGateConfig,
    evaluate_constraint_gate,
)
def _canonical_lump(**kwargs) -> dict:
    base = {
        "amp": 0.2,
        "width": 4.0,
        "center": (0.0, 0.0, 0.0),
        "velocity": (0.1, 0.0, 0.0),
        "omega": 0.0,
        "mode": 0,
        "exotic": 0,
    }
    base.update(kwargs)
    return base


def test_lump_sign_and_profiles_are_deterministic() -> None:
    lump = _canonical_lump()
    point = lump["center"]
    assert lump_sign(lump) == 1
    phi = lump_phi_at(lump, point)
    assert phi > 0.0
    pi = lump_pi_at(lump, (point[0] + 0.5, point[1], point[2]))
    assert abs(pi) > 0.0

    exotic = _canonical_lump(exotic=1)
    assert lump_sign(exotic) == -1


def test_paint_lump_fields_appends_named_channels() -> None:
    nz, ny, nx = 4, 4, 4
    data = np.zeros((nz, ny, nx, 2), dtype=np.float64)
    names = ["phi", "Pi"]
    dx = np.array([2.0, 2.0, 2.0])
    origin = np.array([0.0, 0.0, 0.0])
    lumps = [_canonical_lump(), _canonical_lump(exotic=1, amp=0.15)]

    painted, new_names = paint_lump_fields_on_grid(
        data, names, dx, origin, lumps,
    )

    assert painted.shape[-1] == 6
    assert new_names[-4:] == ["phi_lump0", "Pi_lump0", "phi_lump1", "Pi_lump1"]
    assert np.any(painted[:, :, :, 2] != 0.0)
    assert np.any(painted[:, :, :, 4] != 0.0)


def test_independent_vs_summed_kinetic_differs_on_overlap() -> None:
    """Freeze the legacy mismatch: summed |∇(φ1+φ2)|² ≠ Σ|∇φk|²."""
    lumps = [
        _canonical_lump(center=(-2.0, 0.0, 0.0)),
        _canonical_lump(center=(2.0, 0.0, 0.0), amp=0.18),
    ]
    point = (0.0, 0.0, 0.0)
    phi_sum = sum(lump_phi_at(lump, point) for lump in lumps)
    grad_sum_sq = 0.0
    for axis in range(3):
        lp = list(point)
        lm = list(point)
        eps = 1.0e-3
        lp[axis] += eps
        lm[axis] -= eps
        dphi = (
            sum(lump_phi_at(lump, lp) for lump in lumps)
            - sum(lump_phi_at(lump, lm) for lump in lumps)
        ) / (2.0 * eps)
        grad_sum_sq += dphi * dphi

    independent_grad_sq = 0.0
    for lump in lumps:
        for axis in range(3):
            lp = list(point)
            lm = list(point)
            eps = 1.0e-3
            lp[axis] += eps
            lm[axis] -= eps
            dphi = (lump_phi_at(lump, lp) - lump_phi_at(lump, lm)) / (2.0 * eps)
            independent_grad_sq += dphi * dphi

    assert independent_grad_sq > 0.0
    assert abs(grad_sum_sq - independent_grad_sq) > 1.0e-6


def test_matter_wiring_selects_independent_scalar_model() -> None:
    cfg = GRTresnaConfig(
        lumps=[_canonical_lump(), _canonical_lump(exotic=1, amp=0.12)],
        scalar_mass=0.1,
    )
    overrides = evolution_overrides_from_config(cfg)
    assert overrides["recipe_matter_model"] == GRTRESNA_INDEPENDENT_MATTER_MODEL
    assert overrides["recipe_num_scalar_fields"] == 2
    assert overrides["recipe_scalar_field_signs"] == "1 -1"
    assert overrides["recipe_scalar_mass"] == pytest.approx(0.1)
    assert overrides["recipe_scalar_lambda"] == pytest.approx(0.0)
    plot_vars = overrides["amr.plot_vars"]
    assert "phi_lump0" in plot_vars
    assert "Pi_lump1" in plot_vars


def test_plot_vars_render_unquoted_for_parmparse(tmp_path: Path) -> None:
    example = resolve_example("RadialRecipe")
    out = tmp_path / "params.txt"
    write_params(
        example.template,
        out,
        episode_dir=tmp_path / "episode",
        example=example,
        overrides={"amr.plot_vars": plot_vars_for_independent_scalars(2)},
    )
    line = next(l for l in out.read_text(encoding="utf-8").splitlines() if "amr.plot_vars" in l)
    assert 'amr.plot_vars = "' not in line
    assert "phi_lump0" in line
    assert "Pi_lump1" in line
    assert format_param_value("amr.plot_vars", plot_vars_for_independent_scalars(2)) == (
        "chi h11 h12 h13 h22 h23 h33 K lapse shift1 shift2 shift3 phi Pi "
        "phi_lump0 Pi_lump0 phi_lump1 Pi_lump1"
    )


def test_matter_metadata_round_trip(tmp_path: Path) -> None:
    cfg = GRTresnaConfig(lumps=[_canonical_lump()], scalar_mass=0.05, scalar_lambda=0.04)
    meta_path = write_matter_metadata(tmp_path / "episode.matter.json", cfg)
    payload = json.loads(meta_path.read_text(encoding="utf-8"))
    assert payload["matter_model"] == GRTRESNA_INDEPENDENT_MATTER_MODEL
    assert payload["num_scalar_fields"] == 1
    assert payload["scalar_lambda"] == pytest.approx(0.04)


def test_postload_gate_rejects_poor_constraints(tmp_path: Path) -> None:
    bad = tmp_path / "constraint_norms.dat"
    bad.write_text("0 1.0 0.5\n1 0.9 0.4\n", encoding="utf-8")
    result = evaluate_constraint_gate(
        bad,
        config=PostLoadGateConfig(max_hamiltonian_l2=1.0e-2, max_momentum_l2=1.0e-2),
    )
    assert not result.passed
    assert result.max_hamiltonian_l2 == pytest.approx(1.0)


def test_smoke_canonical_exotic_and_mixed_shell_wiring() -> None:
    canonical = evolution_overrides_from_config(
        GRTresnaConfig(lumps=[_canonical_lump()], scalar_mass=0.1)
    )
    assert canonical["recipe_scalar_field_signs"] == "1"

    exotic = evolution_overrides_from_config(
        GRTresnaConfig(lumps=[_canonical_lump(exotic=1)], scalar_mass=0.1)
    )
    assert exotic["recipe_scalar_field_signs"] == "-1"

    mixed = evolution_overrides_from_config(
        GRTresnaConfig(
            lumps=[
                _canonical_lump(center=(-2.0, 0.0, 0.0)),
                _canonical_lump(exotic=1, amp=0.12, center=(2.0, 0.0, 0.0)),
            ],
            scalar_mass=0.1,
        )
    )
    assert mixed["recipe_num_scalar_fields"] == 2
    assert mixed["recipe_scalar_field_signs"] == "1 -1"


def test_shift_lump_centers_aligns_shell_ansatz_with_grid_center() -> None:
    lumps = [_canonical_lump(center=(3.0, 1.0, 2.0))]
    shifted = shift_lump_centers_for_gridinit(
        lumps,
        grid_center=(32.0, 32.0, 0.0),
    )
    assert shifted[0]["center"] == pytest.approx((35.0, 33.0, 2.0))

    nz, ny, nx = 64, 64, 64
    dx = np.array([2.0, 2.0, 2.0])
    origin = np.array([-32.0, -32.0, -64.0])
    data = np.zeros((nz, ny, nx, 2), dtype=np.float64)
    names = ["chi", "lapse"]
    painted, new_names = paint_lump_fields_on_grid(
        data, names, dx, origin, shifted,
    )
    phi0 = painted[:, :, :, new_names.index("phi_lump0")]
    peak = np.unravel_index(int(np.argmax(phi0)), phi0.shape)
    assert all(26 <= coord <= 38 for coord in peak)


def test_gridinit_v2_round_trip_with_lump_channels(tmp_path: Path) -> None:
    nz, ny, nx = 3, 3, 3
    data = np.zeros((nz, ny, nx, 2), dtype=np.float64)
    names = ["phi", "Pi"]
    dx = np.array([1.0, 1.0, 1.0])
    origin = np.zeros(3)
    painted, names = paint_lump_fields_on_grid(
        data, names, dx, origin, [_canonical_lump(amp=0.1)],
    )
    path = tmp_path / "test.gridinit"
    write_gridinit(painted, names, dx, origin, path)
    assert path.exists()
    assert path.stat().st_size > 0
