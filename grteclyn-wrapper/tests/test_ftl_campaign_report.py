from __future__ import annotations

import json
from pathlib import Path

from grteclyn_wrapper.metrics.diagnostics.ftl_timeseries import read_ftl_timeseries_metrics
from grteclyn_wrapper.search.ftl_campaign_report import (
    FtlSortKey,
    format_ftl_curve,
    load_campaign_ftl_summaries,
    load_eval_ftl_summary,
    rank_summaries,
    summary_from_timeseries,
)


def _write_timeseries(path: Path, rows: list[str]) -> None:
    header = (
        "# time  f_op  f_geo  geo_trustworthy  max_local_speed  superluminal_fraction  "
        "max_shift  structure_coherence  reachable  n_rays  n_reached  max_h_rel_drift"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")


def test_read_timeseries_and_summary_peaks(tmp_path: Path) -> None:
    ts_path = tmp_path / "small_data" / "ftl_timeseries.dat"
    _write_timeseries(
        ts_path,
        [
            "0.0  0.0  0.0  0  0.99  0.0  0.0  nan  1  0  0  0.0",
            "3.2  0.02  0.021  1  1.05  0.05  0.01  nan  1  5  5  1e-4",
            "6.4  0.05  0.033  1  1.10  0.55  0.04  nan  1  5  5  1e-4",
            "16.0  0.0  0.0  1  1.11  0.70  0.10  nan  1  5  5  5e-4",
        ],
    )
    ts = read_ftl_timeseries_metrics(ts_path)
    assert ts is not None
    summary = summary_from_timeseries(ts, eval_id=7, eval_dir=tmp_path)
    assert summary.f_geo_peak == 0.033
    assert summary.t_at_f_geo_peak == 6.4
    assert summary.f_geo_final == 0.0
    assert summary.max_local_speed_peak == 1.11
    assert summary.t_at_max_speed == 16.0
    assert summary.ftl_lifetime_fraction == 0.5


def test_mid_run_alpha_not_final_frame(tmp_path: Path) -> None:
    """Peak f_geo at t=6.4 must win even though final frame is zero."""
    eval_dir = tmp_path / "eval_000127"
    eval_dir.mkdir()
    _write_timeseries(
        eval_dir / "small_data" / "ftl_timeseries.dat",
        [
            "0.0  0.0  0.0  0  0.99  0.0  0.0  nan  1  0  0  0.0",
            "3.2  0.02  0.022  1  1.05  0.05  0.01  nan  1  5  5  1e-4",
            "6.4  0.05  0.033  1  1.10  0.55  0.04  nan  1  5  5  1e-4",
            "16.0  0.0  0.0  1  1.11  0.70  0.10  nan  1  5  5  5e-4",
        ],
    )
    summary = load_eval_ftl_summary(eval_dir)
    assert summary is not None
    assert summary.f_geo_peak == 0.033
    assert summary.t_at_f_geo_peak == 6.4
    assert summary.f_geo_final == 0.0

    curve = format_ftl_curve(summary)
    assert " *" in curve
    assert "6.4" in curve


def test_dat_preferred_and_reaggregates_wrong_score_json_peaks(tmp_path: Path) -> None:
    eval_dir = tmp_path / "eval_000010"
    eval_dir.mkdir()
    _write_timeseries(
        eval_dir / "small_data" / "ftl_timeseries.dat",
        ["6.4  0.04  0.05  1  1.06  0.1  0.02  nan  1  5  5  1e-4"],
    )
    (eval_dir / "score.json").write_text(
        json.dumps(
            {
                "score": {"total": 12.5, "components": {}},
                "metrics": {
                    "ftl_timeseries": {
                        "n_frames": 1,
                        "t": [6.4],
                        "f_op": [0.04],
                        "f_geo": [0.01],
                        "geo_trustworthy": [True],
                        "max_local_speed": [1.06],
                        "superluminal_fraction": [0.1],
                        "structure_coherence": [float("nan")],
                        "max_h_rel_drift": [1e-4],
                        "f_geo_peak": 0.01,
                        "t_at_f_geo_peak": 6.4,
                        "f_op_peak": 0.04,
                        "t_at_f_op_peak": 6.4,
                        "ftl_lifetime_fraction": 1.0,
                        "op_lifetime_fraction": 1.0,
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    summary = load_eval_ftl_summary(eval_dir)
    assert summary is not None
    assert summary.f_geo_peak == 0.05


def test_load_eval_from_score_json_when_no_dat(tmp_path: Path) -> None:
    eval_dir = tmp_path / "eval_000010"
    eval_dir.mkdir()
    (eval_dir / "score.json").write_text(
        json.dumps(
            {
                "score": {
                    "total": 12.5,
                    "components": {
                        "ftl_geo_timeavg": 0.01,
                        "operational_ftl_geodesic": 0.01,
                        "structural_persistence": 0.95,
                    },
                },
                "metrics": {
                    "ftl_timeseries": {
                        "n_frames": 2,
                        "t": [3.2, 6.4],
                        "f_op": [0.02, 0.0],
                        "f_geo": [0.01, 0.04],
                        "geo_trustworthy": [True, True],
                        "max_local_speed": [1.05, 1.10],
                        "superluminal_fraction": [0.05, 0.55],
                        "structure_coherence": [float("nan"), float("nan")],
                        "max_h_rel_drift": [1e-4, 1e-4],
                        "f_geo_peak": 0.01,
                        "t_at_f_geo_peak": 3.2,
                        "f_op_peak": 0.02,
                        "t_at_f_op_peak": 3.2,
                        "ftl_lifetime_fraction": 1.0,
                        "op_lifetime_fraction": 0.5,
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    summary = load_eval_ftl_summary(eval_dir, status="gpu_ok")
    assert summary is not None
    assert summary.eval_id == 10
    assert summary.score == 12.5
    assert summary.status == "gpu_ok"
    assert summary.structural_persistence == 0.95
    assert summary.f_geo_peak == 0.04
    assert summary.t_at_f_geo_peak == 6.4


def test_trajectory_fallback_when_eval_dir_pruned(tmp_path: Path) -> None:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    (campaign / "trajectory.jsonl").write_text(
        json.dumps(
            {
                "eval": 42,
                "status": "gpu_ok",
                "score": 9.0,
                "descriptor_details": {
                    "f_geo_peak": 0.04,
                    "t_at_f_geo_peak": 6.4,
                    "ftl_lifetime": 0.5,
                    "n_frames": 4,
                    "ftl_geo_timeavg": 0.02,
                },
                "components": {"structural_persistence": 0.9},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    summaries = load_campaign_ftl_summaries(campaign)
    assert len(summaries) == 1
    assert summaries[0].eval_id == 42
    assert summaries[0].f_geo_peak == 0.04
    assert summaries[0].t_at_f_geo_peak == 6.4
    assert summaries[0].has_frame_curve is False


def test_trajectory_fallback_includes_max_speed(tmp_path: Path) -> None:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    (campaign / "trajectory.jsonl").write_text(
        json.dumps(
            {
                "eval": 42,
                "status": "gpu_ok",
                "score": 9.0,
                "descriptor_details": {
                    "f_geo_peak": 0.04,
                    "t_at_f_geo_peak": 6.4,
                    "f_op_peak": 0.03,
                    "t_at_f_op_peak": 6.4,
                    "max_local_speed_peak": 1.12,
                    "t_at_max_speed": 16.0,
                    "superluminal_fraction_peak": 0.55,
                    "t_at_superluminal_peak": 9.6,
                    "ftl_lifetime": 0.5,
                    "n_frames": 4,
                    "ftl_geo_timeavg": 0.02,
                },
                "components": {"structural_persistence": 0.9, "ftl_geo_timeavg": 0.02},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    summaries = load_campaign_ftl_summaries(campaign)
    assert len(summaries) == 1
    assert summaries[0].max_local_speed_peak == 1.12
    assert summaries[0].f_op_peak == 0.03


def test_timeseries_speed_peak_aggregate(tmp_path: Path) -> None:
    ts_path = tmp_path / "small_data" / "ftl_timeseries.dat"
    _write_timeseries(
        ts_path,
        [
            "0.0  0.0  0.0  0  0.99  0.0  0.0  nan  1  0  0  0.0",
            "6.4  0.05  0.033  1  1.10  0.55  0.04  nan  1  5  5  1e-4",
            "16.0  0.0  0.0  1  1.11  0.70  0.10  nan  1  5  5  5e-4",
        ],
    )
    ts = read_ftl_timeseries_metrics(ts_path)
    assert ts is not None
    assert ts.max_local_speed_peak == 1.11
    assert ts.t_at_max_speed == 16.0
    assert ts.superluminal_fraction_peak == 0.70


def test_rank_campaign_by_max_speed(tmp_path: Path) -> None:
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    for eval_id, speed in ((1, 1.02), (2, 1.08)):
        eval_dir = campaign / f"eval_{eval_id:06d}"
        eval_dir.mkdir()
        _write_timeseries(
            eval_dir / "small_data" / "ftl_timeseries.dat",
            [f"3.2  0.0  0.0  0  {speed:.2f}  0.0  0.0  nan  1  0  0  0.0"],
        )
    (campaign / "trajectory.jsonl").write_text(
        '{"eval": 1, "status": "gpu_ok", "score": 1.0}\n'
        '{"eval": 2, "status": "gpu_ok", "score": 2.0}\n',
        encoding="utf-8",
    )
    summaries = load_campaign_ftl_summaries(campaign)
    ranked = rank_summaries(summaries, sort_by=FtlSortKey.MAX_SPEED)
    assert [s.eval_id for s in ranked] == [2, 1]
