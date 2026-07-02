"""Tests for HQ replay matter-parameter restoration."""

from __future__ import annotations

import json
import sys
from pathlib import Path

# replay_eval lives under scripts/campaigns/hq/
_REPLAY_DIR = Path(__file__).resolve().parents[2] / "scripts" / "campaigns" / "hq"
if str(_REPLAY_DIR) not in sys.path:
    sys.path.insert(0, str(_REPLAY_DIR))

import replay_eval  # noqa: E402


def test_load_matter_replay_overrides_from_matter_json(tmp_path) -> None:
    eval_dir = tmp_path / "eval_000001"
    eval_dir.mkdir()
    payload = {
        "matter_model": "grtresna_bicomplex_scalar",
        "num_scalar_fields": 2,
        "scalar_field_signs": [1, -1],
        "scalar_mass": 1.0,
        "scalar_lambda": 640.0,
        "scalar_mu": 85333.0,
        "lump_count": 2,
        "bs_omega": 0.4,
    }
    (eval_dir / "initial_data.matter.json").write_text(
        json.dumps(payload), encoding="utf-8"
    )

    overrides = replay_eval._load_matter_replay_overrides(eval_dir)
    assert overrides["recipe_scalar_mu"] == 85333.0
    assert overrides["recipe_scalar_lambda"] == 640.0
    assert overrides["recipe_matter_model"] == "grtresna_bicomplex_scalar"


def test_load_matter_replay_overrides_from_params_txt(tmp_path) -> None:
    eval_dir = tmp_path / "eval_000002"
    eval_dir.mkdir()
    (eval_dir / "params.txt").write_text(
        "\n".join(
            [
                "recipe_matter_model = grtresna_bicomplex_scalar",
                "recipe_scalar_mass = 1.0",
                "recipe_scalar_lambda = 640.0",
                "recipe_scalar_mu = 85333.0",
                "recipe_num_scalar_fields = 2",
                "recipe_scalar_field_signs = 1 -1",
            ]
        ),
        encoding="utf-8",
    )

    overrides = replay_eval._load_matter_replay_overrides(eval_dir)
    assert overrides["recipe_scalar_mu"] == 85333.0
    assert overrides["recipe_scalar_lambda"] == 640.0


def test_apply_replay_consumer_env_boson_fields(monkeypatch) -> None:
    monkeypatch.delenv("GRTECLYN_FRAMES_FIELDS", raising=False)
    replay_eval._apply_replay_consumer_env(
        {"grtresna_matter_model": "grtresna_bicomplex_scalar"}
    )
    fields = replay_eval.os.environ["GRTECLYN_FRAMES_FIELDS"]
    assert "scalar_activity" in fields
    assert "phi_lump0" in fields
    assert "Weyl4_Re" in fields
    assert "chi" in fields
