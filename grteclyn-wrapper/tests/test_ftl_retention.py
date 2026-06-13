from __future__ import annotations

import json
from pathlib import Path

from grteclyn_wrapper.search.ftl_retention import (
    FTL_CHAMPIONS_FILE,
    FTL_RETENTION_LOG,
    FtlChampionBoard,
    append_ftl_retention_events,
    compute_keep_eval_ids,
    save_ftl_champions,
)
from grteclyn_wrapper.search.qd_search.io import _prune_eval_dirs


def _gpu_ok_record(
    eval_id: int,
    *,
    score: float,
    f_geo_peak: float = 0.0,
    max_speed: float = 0.0,
    lifetime: float = 0.0,
    campaign: Path,
) -> dict:
    eval_dir = campaign / f"eval_{eval_id:06d}"
    eval_dir.mkdir(parents=True, exist_ok=True)
    return {
        "eval": eval_id,
        "status": "gpu_ok",
        "score": score,
        "episode": str(eval_dir),
        "descriptor_details": {
            "f_geo_peak": f_geo_peak,
            "t_at_f_geo_peak": 6.4,
            "f_op_peak": 0.0,
            "max_local_speed_peak": max_speed,
            "t_at_max_speed": 6.4 if max_speed > 1.0 else float("nan"),
            "superluminal_fraction_peak": 0.0,
            "ftl_lifetime": lifetime,
            "ftl_lifetime_fraction": lifetime,
            "n_frames": 7.0,
            "ftl_geo_timeavg": 0.01,
        },
        "components": {"ftl_geo_timeavg": 0.01},
    }


def test_champion_board_crowns_and_replaces() -> None:
    board = FtlChampionBoard()
    r1 = {
        "eval": 1,
        "status": "gpu_ok",
        "score": -10.0,
        "episode": "/tmp/e1",
        "descriptor_details": {"f_geo_peak": 0.02, "t_at_f_geo_peak": 6.4, "n_frames": 7.0},
        "components": {},
    }
    events = board.consider(r1)
    assert len(events) == 1
    assert events[0].metric == "f_geo_peak"
    assert board.champions["f_geo_peak"].eval_id == 1

    r2 = {
        "eval": 2,
        "status": "gpu_ok",
        "score": -20.0,
        "episode": "/tmp/e2",
        "descriptor_details": {"f_geo_peak": 0.05, "t_at_f_geo_peak": 9.6, "n_frames": 7.0},
        "components": {},
    }
    events = board.consider(r2)
    assert events[0].replaced_eval == 1
    assert board.champions["f_geo_peak"].eval_id == 2


def test_compute_keep_eval_ids_unions_score_and_champions(tmp_path: Path) -> None:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    records = [
        _gpu_ok_record(1, score=100.0, f_geo_peak=0.01, campaign=campaign),
        _gpu_ok_record(2, score=-50.0, f_geo_peak=0.05, campaign=campaign),
        _gpu_ok_record(3, score=-40.0, max_speed=1.12, campaign=campaign),
    ]
    board = FtlChampionBoard.rebuild(records)
    keep = compute_keep_eval_ids(
        records,
        keep_top_score=1,
        board=board,
        ftl_retention_enabled=True,
        protect_eval_ids=set(),
    )
    assert keep is not None
    assert 1 in keep
    assert 2 in keep
    assert 3 in keep


def test_prune_keeps_ftl_champion_outside_top_score(tmp_path: Path) -> None:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    records = [
        _gpu_ok_record(1, score=100.0, f_geo_peak=0.01, campaign=campaign),
        _gpu_ok_record(2, score=-50.0, f_geo_peak=0.05, campaign=campaign),
    ]
    board = FtlChampionBoard.rebuild(records)
    keep = compute_keep_eval_ids(
        records,
        keep_top_score=1,
        board=board,
        ftl_retention_enabled=True,
        protect_eval_ids=set(),
    )
    assert keep is not None
    deleted = _prune_eval_dirs(campaign, records, keep_eval_ids=keep)
    assert deleted == 0
    assert (campaign / "eval_000001").is_dir()
    assert (campaign / "eval_000002").is_dir()

    deleted = _prune_eval_dirs(campaign, records, keep_eval_ids={1})
    assert deleted == 1
    assert not (campaign / "eval_000002").is_dir()


def test_retention_log_and_champions_snapshot(tmp_path: Path) -> None:
    board = FtlChampionBoard()
    record = {
        "eval": 7,
        "status": "gpu_ok",
        "score": 3.0,
        "episode": "/tmp/e7",
        "descriptor_details": {"f_geo_peak": 0.04, "t_at_f_geo_peak": 6.4, "n_frames": 7.0},
        "components": {"ftl_geo_timeavg": 0.02},
    }
    events = board.consider(record)
    log_path = tmp_path / FTL_RETENTION_LOG
    champions_path = tmp_path / FTL_CHAMPIONS_FILE
    append_ftl_retention_events(log_path, events)
    save_ftl_champions(champions_path, board)

    lines = log_path.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) >= 1
    metrics = {json.loads(line)["metric"] for line in lines}
    assert "f_geo_peak" in metrics

    snap = json.loads(champions_path.read_text(encoding="utf-8"))
    assert snap["f_geo_peak"]["eval"] == 7
