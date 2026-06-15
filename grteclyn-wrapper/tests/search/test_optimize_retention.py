"""Smoke tests for CMA-ES (optimize) eval-dir retention + FTL champions.

Mirrors the QD retention behaviour: the optimize loop should keep only the
top-N-by-score eval dirs on disk plus one champion dir per FTL peak metric,
and write ftl_champions.json / ftl_retention.jsonl.
"""

from __future__ import annotations

import json
from pathlib import Path

from grteclyn_wrapper.search.ftl_retention import (
    FTL_CHAMPIONS_FILE,
    FTL_RETENTION_LOG,
    FtlChampionBoard,
)
from grteclyn_wrapper.search.optimize.driver import _apply_optimize_retention


def _gpu_ok_record(
    eval_id: int,
    *,
    score: float,
    f_geo_peak: float = 0.0,
    max_speed: float = 0.0,
    opt_dir: Path,
) -> dict:
    eval_dir = opt_dir / f"eval_{eval_id:06d}"
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
            "ftl_lifetime": 0.0,
            "ftl_lifetime_fraction": 0.0,
            "n_frames": 7.0,
            "ftl_geo_timeavg": 0.01,
        },
        "components": {"ftl_geo_timeavg": 0.01},
    }


def test_optimize_retention_keeps_top_and_ftl_champion(tmp_path: Path) -> None:
    opt_dir = tmp_path / "optimize_run"
    opt_dir.mkdir()
    # Prior generations (not protected): a top scorer, an FTL champion with a low
    # score, and a candidate that is neither.
    prior = [
        _gpu_ok_record(1, score=100.0, f_geo_peak=0.01, opt_dir=opt_dir),  # top score
        _gpu_ok_record(2, score=-50.0, f_geo_peak=0.09, opt_dir=opt_dir),  # f_geo champ
        _gpu_ok_record(3, score=10.0, f_geo_peak=0.0, opt_dir=opt_dir),    # prunable
    ]
    board = FtlChampionBoard.rebuild(prior)
    # Current generation crowns a new speed record -> exercises the audit log.
    current = [_gpu_ok_record(4, score=5.0, f_geo_peak=0.0, max_speed=1.5, opt_dir=opt_dir)]
    trajectory = [*prior, *current]

    _apply_optimize_retention(
        opt_dir=opt_dir,
        trajectory=trajectory,
        new_records=current,
        board=board,
        keep_top_eval_dirs=1,
        ftl_retention_enabled=True,
        champions_path=opt_dir / FTL_CHAMPIONS_FILE,
        retention_path=opt_dir / FTL_RETENTION_LOG,
    )

    # eval 1 = top score (kept), eval 2 = f_geo champion (kept despite low score),
    # eval 4 = current generation (protected), eval 3 = neither -> pruned.
    assert (opt_dir / "eval_000001").is_dir()
    assert (opt_dir / "eval_000002").is_dir()
    assert (opt_dir / "eval_000004").is_dir()
    assert not (opt_dir / "eval_000003").is_dir()

    champs = json.loads((opt_dir / FTL_CHAMPIONS_FILE).read_text(encoding="utf-8"))
    assert champs["f_geo_peak"]["eval"] == 2
    assert (opt_dir / FTL_RETENTION_LOG).exists()


def test_optimize_retention_top_only_prunes_ftl_peak(tmp_path: Path) -> None:
    """With ftl_retention off, only score decides: a high-FTL low-score dir is pruned."""
    opt_dir = tmp_path / "optimize_run"
    opt_dir.mkdir()
    prior = [
        _gpu_ok_record(1, score=100.0, f_geo_peak=0.01, opt_dir=opt_dir),
        _gpu_ok_record(2, score=-50.0, f_geo_peak=0.09, opt_dir=opt_dir),
    ]
    current = [_gpu_ok_record(3, score=5.0, f_geo_peak=0.0, opt_dir=opt_dir)]
    trajectory = [*prior, *current]

    _apply_optimize_retention(
        opt_dir=opt_dir,
        trajectory=trajectory,
        new_records=current,
        board=None,
        keep_top_eval_dirs=1,
        ftl_retention_enabled=False,
        champions_path=opt_dir / FTL_CHAMPIONS_FILE,
        retention_path=opt_dir / FTL_RETENTION_LOG,
    )

    # top-1 = eval 1, eval 3 protected (current gen); eval 2 has high f_geo but
    # retention is off, so it is pruned.
    assert (opt_dir / "eval_000001").is_dir()
    assert (opt_dir / "eval_000003").is_dir()
    assert not (opt_dir / "eval_000002").is_dir()
    assert not (opt_dir / FTL_CHAMPIONS_FILE).exists()


def test_optimize_retention_protects_current_generation(tmp_path: Path) -> None:
    """The just-evaluated generation is never pruned, even below the top-N cut."""
    opt_dir = tmp_path / "optimize_run"
    opt_dir.mkdir()
    old = _gpu_ok_record(1, score=100.0, opt_dir=opt_dir)
    current = [
        _gpu_ok_record(2, score=-99.0, opt_dir=opt_dir),
        _gpu_ok_record(3, score=-98.0, opt_dir=opt_dir),
    ]
    trajectory = [old, *current]

    _apply_optimize_retention(
        opt_dir=opt_dir,
        trajectory=trajectory,
        new_records=current,
        board=None,
        keep_top_eval_dirs=1,
        ftl_retention_enabled=False,
        champions_path=opt_dir / FTL_CHAMPIONS_FILE,
        retention_path=opt_dir / FTL_RETENTION_LOG,
    )

    # eval 1 is top-1; evals 2,3 are the protected current generation.
    assert (opt_dir / "eval_000001").is_dir()
    assert (opt_dir / "eval_000002").is_dir()
    assert (opt_dir / "eval_000003").is_dir()
