"""Provenance guard: matter.json + params.txt must agree on signs."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from grteclyn_wrapper.grtresna.matter.sign_consistency import (
    SignMismatchError,
    check_id_evolution_sign_consistency,
    evolution_effective_signs,
)
from grteclyn_wrapper.grtresna.matter.wiring import read_matter_metadata
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig


def _assert_artifact_sign_consistency(episode_dir: Path) -> None:
    """Fail if initial_data.matter.json disagrees with params.txt evolution signs."""
    matter_path = episode_dir / "initial_data.matter.json"
    params_path = episode_dir / "params.txt"
    if not matter_path.is_file() or not params_path.is_file():
        pytest.skip("episode artifacts missing")

    meta = read_matter_metadata(matter_path)
    params: dict[str, str] = {}
    for line in params_path.read_text(encoding="utf-8").splitlines():
        body = line.split("#", 1)[0].strip()
        if "=" not in body:
            continue
        key, val = (part.strip() for part in body.split("=", 1))
        params[key] = val

    model = params.get("recipe_matter_model", meta.matter_model)
    if model == "grtresna_bicomplex_scalar":
        evo_signs = tuple(
            int(x) for x in params.get("recipe_scalar_field_signs", "").split()
        )
        assert evo_signs == meta.scalar_field_signs
        assert evo_signs  # must not be empty for exotic campaigns
        return

    if model == "grtresna_complex_scalar":
        # Reconstruct a cfg from metadata + any lump exotic flags in params.
        # Continuous trajectory_lump*_exotic in params are search knobs; the
        # evolution sign is global recipe_scalar_sign.
        sign = int(float(params.get("recipe_scalar_sign", meta.scalar_sign or 1)))
        # If params still list exotic lumps and model is single-complex with
        # sign=+1, that is the eval-118 failure mode.
        exotic_keys = [
            k for k in params if k.startswith("trajectory_lump") and k.endswith("_exotic")
        ]
        if exotic_keys and sign >= 0:
            rounded = [int(round(float(params[k]))) for k in sorted(exotic_keys)]
            if any(v == 1 for v in rounded):
                raise SignMismatchError(
                    f"params.txt has exotic lumps {rounded} but "
                    f"recipe_matter_model={model} with scalar_sign={sign}"
                )


def test_eval118_artifacts_are_flagged_as_mismatch() -> None:
    """Real corrupted eval-118 directory must fail the provenance guard."""
    episode = Path(
        "/home/jovyan/nachevsky/test/simulation/GRTeclyn/"
        "runs/grtresna_qd/qball_traj_spiral_v2/eval_000118"
    )
    if not episode.is_dir():
        pytest.skip("eval_000118 not present")
    with pytest.raises(SignMismatchError):
        _assert_artifact_sign_consistency(episode)


def test_synthetic_bicomplex_artifacts_pass(tmp_path: Path) -> None:
    matter = {
        "matter_model": "grtresna_bicomplex_scalar",
        "num_scalar_fields": 2,
        "scalar_field_signs": [1, -1],
        "scalar_mass": 1.0,
        "scalar_lambda": 640.0,
        "scalar_mu": 85333.0,
        "lump_count": 2,
        "scalar_sign": 1,
        "bs_omega": 0.8,
    }
    (tmp_path / "initial_data.matter.json").write_text(
        json.dumps(matter), encoding="utf-8"
    )
    (tmp_path / "params.txt").write_text(
        "recipe_matter_model = grtresna_bicomplex_scalar\n"
        "recipe_scalar_field_signs = 1 -1\n",
        encoding="utf-8",
    )
    _assert_artifact_sign_consistency(tmp_path)


def test_cfg_check_matches_guard_logic() -> None:
    bad = GRTresnaConfig(
        matter_model="grtresna_complex_scalar",
        scalar_sign=1,
        lumps=[
            {"amp": 0.1, "exotic": 0, "center": (0, 0, 0)},
            {"amp": 0.1, "exotic": 1, "center": (1, 0, 0)},
        ],
    )
    assert check_id_evolution_sign_consistency(bad).ok is False
    evo = evolution_effective_signs(
        "grtresna_complex_scalar", scalar_sign=1, n_lumps=2
    )
    assert evo == (1, 1)
