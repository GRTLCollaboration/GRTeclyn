from __future__ import annotations

import json
from pathlib import Path

from grteclyn_wrapper.rl.seed import (
    load_elite_overrides,
    lump_seeds_from_eval,
    pump_geometry_from_overrides,
    rl_pump_params,
)


def test_pump_geometry_defaults() -> None:
    radius, width = pump_geometry_from_overrides({})
    assert radius == 5.0
    assert width == 1.5


def test_pump_geometry_from_overrides() -> None:
    radius, width = pump_geometry_from_overrides(
        {"recipe_basis_radius_max": 7.0, "recipe_basis_width": 2.0}
    )
    assert radius == 7.0
    assert width == 2.0


def _write_matter(eval_dir: Path, payload: dict) -> None:
    eval_dir.mkdir(parents=True, exist_ok=True)
    (eval_dir / "initial_data.matter.json").write_text(
        json.dumps(payload), encoding="utf-8"
    )


def test_lump_seeds_from_matter_centers(tmp_path: Path) -> None:
    _write_matter(
        tmp_path,
        {
            "lump_count": 3,
            "lump_centers": [[1.0, 0.0, 0.0], [-1.0, 2.0, 0.0], [0.0, 0.0, 3.0]],
        },
    )
    num_lumps, centers = lump_seeds_from_eval(tmp_path)
    assert num_lumps == 3
    assert centers[1] == (-1.0, 2.0, 0.0)


def test_lump_seeds_fallback_single_spotlight(tmp_path: Path) -> None:
    _write_matter(tmp_path, {"lump_count": 0})
    num_lumps, centers = lump_seeds_from_eval(tmp_path)
    assert num_lumps == 1
    assert centers == [(0.0, 0.0, 0.0)]


def test_lump_seeds_capped_at_max(tmp_path: Path) -> None:
    _write_matter(
        tmp_path,
        {"lump_centers": [[float(i), 0.0, 0.0] for i in range(8)]},
    )
    num_lumps, centers = lump_seeds_from_eval(tmp_path, max_lumps=5)
    assert num_lumps == 5
    assert len(centers) == 5


def test_rl_pump_params_emits_seed_strings(tmp_path: Path) -> None:
    _write_matter(
        tmp_path,
        {
            "lump_centers": [[1.5, 0.0, 0.0], [0.0, -2.5, 0.0]],
            "bs_profile_width": 2.0,
        },
    )
    params = rl_pump_params(tmp_path)
    assert params["rl_num_lumps"] == 2
    assert len(params["rl_lump_seed_x"].split()) == 2
    assert len(params["rl_lump_seed_y"].split()) == 2
    assert params["rl_pump_width"] == 2.0
    # x seeds carry the two lump x-offsets
    xs = [float(v) for v in params["rl_lump_seed_x"].split()]
    assert xs == [1.5, 0.0]
