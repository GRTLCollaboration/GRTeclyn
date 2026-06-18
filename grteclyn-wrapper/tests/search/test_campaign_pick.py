"""Tests for campaign trajectory pick helpers."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from grteclyn_wrapper.search.campaign_pick import (
    load_gpu_ok_records,
    pick_top_eval_id,
    pick_top_eval_ids,
)


def _write_trajectory(path: Path, rows: list[dict]) -> None:
    path.write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n",
        encoding="utf-8",
    )


def test_pick_top_eval_id_by_score(tmp_path: Path) -> None:
    traj = tmp_path / "trajectory.jsonl"
    _write_trajectory(
        traj,
        [
            {"eval": 1, "status": "gpu_ok", "score": 10.0},
            {"eval": 2, "status": "gpu_ok", "score": 50.0},
            {"eval": 3, "status": "gpu_ok", "score": 30.0},
        ],
    )
    assert pick_top_eval_id(traj) == 2
    assert pick_top_eval_ids(traj, top_k=2) == [2, 3]


def test_pick_top_eval_skips_non_gpu_ok(tmp_path: Path) -> None:
    traj = tmp_path / "trajectory.jsonl"
    _write_trajectory(
        traj,
        [
            {"eval": 1, "status": "grtresna_rejected", "score": 999.0},
            {"eval": 2, "status": "gpu_ok", "score": 5.0},
        ],
    )
    assert pick_top_eval_id(traj) == 2
    assert load_gpu_ok_records(traj)[0]["eval"] == 2


def test_pick_top_eval_rank_by_ftl_geo(tmp_path: Path) -> None:
    traj = tmp_path / "trajectory.jsonl"
    _write_trajectory(
        traj,
        [
            {
                "eval": 10,
                "status": "gpu_ok",
                "score": 100.0,
                "components": {"ftl_geo_evolving": 0.10},
            },
            {
                "eval": 20,
                "status": "gpu_ok",
                "score": 80.0,
                "components": {"ftl_geo_evolving": 0.25},
            },
        ],
    )
    assert pick_top_eval_id(traj, rank_by="ftl_geo_evolving") == 20


def test_pick_top_eval_empty_raises(tmp_path: Path) -> None:
    traj = tmp_path / "trajectory.jsonl"
    _write_trajectory(traj, [{"eval": 1, "status": "grtresna_rejected", "score": 1.0}])
    with pytest.raises(ValueError, match="no gpu_ok"):
        pick_top_eval_id(traj)


def test_pick_top_eval_includes_dry_run_when_requested(tmp_path: Path) -> None:
    traj = tmp_path / "trajectory.jsonl"
    _write_trajectory(
        traj,
        [
            {"eval": 1, "status": "dry_run", "dry_run": True, "score": 5.0},
            {"eval": 2, "status": "dry_run", "dry_run": True, "score": 10.0},
        ],
    )
    with pytest.raises(ValueError, match="no gpu_ok"):
        pick_top_eval_id(traj)
    assert pick_top_eval_id(traj, include_dry_run=True) == 2
