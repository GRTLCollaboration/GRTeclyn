"""Summarize per-frame FTL metrics from completed QD / search campaign evals.

Peaks are computed by scanning **every** frame in ``small_data/ftl_timeseries.dat``
(or the per-frame arrays in ``score.json``), not from the final frame alone.
When eval dirs are pruned, falls back to peak fields stored in ``trajectory.jsonl``.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Iterable, Sequence

from ..metrics.diagnostics.ftl_timeseries import (
    read_ftl_timeseries_metrics,
    reaggregate_ftl_timeseries,
)
from ..metrics.types.diagnostics import FtlTimeSeriesMetrics
from .ftl_peak_metrics import from_trajectory_record as ftl_peaks_from_trajectory

_EVAL_DIR_RE = re.compile(r"^eval_(\d+)$")


class FtlSortKey(str, Enum):
    F_GEO_PEAK = "f_geo_peak"
    F_OP_PEAK = "f_op_peak"
    MAX_SPEED = "max_speed"
    SUPERLUMINAL = "superluminal"
    LIFETIME = "lifetime"
    TIMEAVG = "timeavg"
    SCORE = "score"


@dataclass(frozen=True)
class EvalFtlSummary:
    """Peak and end-state FTL observables for one eval episode."""

    eval_id: int
    eval_dir: Path
    status: str | None
    score: float | None
    n_frames: int
    f_geo_peak: float
    t_at_f_geo_peak: float | None
    f_geo_final: float
    f_op_peak: float
    t_at_f_op_peak: float | None
    f_op_final: float
    max_local_speed_peak: float
    t_at_max_speed: float | None
    superluminal_fraction_peak: float
    t_at_superluminal_peak: float | None
    ftl_lifetime_fraction: float
    op_lifetime_fraction: float
    ftl_geo_timeavg: float | None
    operational_ftl_geodesic: float | None
    structural_persistence: float | None
    has_frame_curve: bool = True


def _eval_id_from_dir(eval_dir: Path) -> int | None:
    match = _EVAL_DIR_RE.match(eval_dir.name)
    return int(match.group(1)) if match else None


def _peak_with_time(values: Sequence[float], times: Sequence[float]) -> tuple[float, float | None]:
    peak = float("-inf")
    t_peak: float | None = None
    for value, time in zip(values, times):
        if not (value == value):  # NaN
            continue
        if value > peak:
            peak = value
            t_peak = float(time)
    if peak == float("-inf"):
        return 0.0, None
    return peak, t_peak


def _final_value(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    last = float(values[-1])
    return last if last == last else 0.0


def summary_from_timeseries(
    ts: FtlTimeSeriesMetrics,
    *,
    eval_id: int,
    eval_dir: Path,
    status: str | None = None,
    score: float | None = None,
    ftl_geo_timeavg: float | None = None,
    operational_ftl_geodesic: float | None = None,
    structural_persistence: float | None = None,
) -> EvalFtlSummary:
    speed_peak, t_speed = _peak_with_time(ts.max_local_speed, ts.t)
    sup_peak, t_sup = _peak_with_time(ts.superluminal_fraction, ts.t)
    return EvalFtlSummary(
        eval_id=eval_id,
        eval_dir=eval_dir,
        status=status,
        score=score,
        n_frames=ts.n_frames,
        f_geo_peak=ts.f_geo_peak,
        t_at_f_geo_peak=ts.t_at_f_geo_peak,
        f_geo_final=_final_value(ts.f_geo),
        f_op_peak=ts.f_op_peak,
        t_at_f_op_peak=ts.t_at_f_op_peak,
        f_op_final=_final_value(ts.f_op),
        max_local_speed_peak=speed_peak,
        t_at_max_speed=t_speed,
        superluminal_fraction_peak=sup_peak,
        t_at_superluminal_peak=t_sup,
        ftl_lifetime_fraction=ts.ftl_lifetime_fraction,
        op_lifetime_fraction=ts.op_lifetime_fraction,
        ftl_geo_timeavg=ftl_geo_timeavg,
        operational_ftl_geodesic=operational_ftl_geodesic,
        structural_persistence=structural_persistence,
        has_frame_curve=True,
    )


def _load_score_payload(eval_dir: Path) -> dict[str, Any] | None:
    score_path = eval_dir / "score.json"
    if not score_path.is_file():
        return None
    try:
        return json.loads(score_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _timeseries_arrays_from_score_metrics(metrics: dict[str, Any]) -> FtlTimeSeriesMetrics | None:
    raw = metrics.get("ftl_timeseries")
    if not isinstance(raw, dict) or not raw.get("n_frames"):
        return None
    try:
        return FtlTimeSeriesMetrics(
            n_frames=int(raw["n_frames"]),
            t=tuple(float(x) for x in raw["t"]),
            f_op=tuple(float(x) for x in raw["f_op"]),
            f_geo=tuple(float(x) for x in raw["f_geo"]),
            geo_trustworthy=tuple(bool(x) for x in raw["geo_trustworthy"]),
            max_local_speed=tuple(float(x) for x in raw["max_local_speed"]),
            superluminal_fraction=tuple(float(x) for x in raw["superluminal_fraction"]),
            structure_coherence=tuple(
                float("nan") if x != x else float(x) for x in raw.get("structure_coherence", [])
            ),
            max_h_rel_drift=tuple(float(x) for x in raw["max_h_rel_drift"]),
            f_geo_peak=0.0,
            t_at_f_geo_peak=None,
            f_op_peak=0.0,
            t_at_f_op_peak=None,
            max_local_speed_peak=0.0,
            t_at_max_speed=None,
            superluminal_fraction_peak=0.0,
            t_at_superluminal_peak=None,
            ftl_lifetime_fraction=0.0,
            op_lifetime_fraction=0.0,
        )
    except (KeyError, TypeError, ValueError):
        return None


def load_eval_timeseries(eval_dir: Path) -> FtlTimeSeriesMetrics | None:
    """Load per-frame FTL arrays, preferring the raw ``.dat`` log on disk."""
    eval_dir = eval_dir.resolve()
    ts = read_ftl_timeseries_metrics(eval_dir / "small_data" / "ftl_timeseries.dat")
    if ts is not None:
        return ts

    payload = _load_score_payload(eval_dir)
    metrics = (payload or {}).get("metrics") or {}
    arrays = _timeseries_arrays_from_score_metrics(metrics)
    if arrays is None:
        return None
    return reaggregate_ftl_timeseries(arrays)


def load_eval_ftl_summary(
    eval_dir: Path,
    *,
    status: str | None = None,
) -> EvalFtlSummary | None:
    """Load one eval's time-resolved FTL summary by scanning all frames."""
    eval_dir = eval_dir.resolve()
    eval_id = _eval_id_from_dir(eval_dir)
    if eval_id is None:
        return None

    ts = load_eval_timeseries(eval_dir)
    if ts is None:
        return None

    payload = _load_score_payload(eval_dir)
    score_block = (payload or {}).get("score") or {}
    components = score_block.get("components") or {}
    score_total = score_block.get("total")
    return summary_from_timeseries(
        ts,
        eval_id=eval_id,
        eval_dir=eval_dir,
        status=status,
        score=float(score_total) if score_total is not None else None,
        ftl_geo_timeavg=_optional_float(components.get("ftl_geo_timeavg")),
        operational_ftl_geodesic=_optional_float(components.get("operational_ftl_geodesic")),
        structural_persistence=_optional_float(components.get("structural_persistence")),
    )


def _optional_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if out == out else None


def _iter_trajectory_records(campaign_dir: Path) -> Iterable[dict[str, Any]]:
    path = campaign_dir / "trajectory.jsonl"
    if not path.is_file():
        return
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(record, dict):
            yield record


def load_trajectory_status(campaign_dir: Path) -> dict[int, str]:
    status_by_eval: dict[int, str] = {}
    for record in _iter_trajectory_records(campaign_dir):
        eval_id = record.get("eval")
        status = record.get("status")
        if isinstance(eval_id, int) and isinstance(status, str):
            status_by_eval[eval_id] = status
    return status_by_eval


def summary_from_trajectory_record(
    record: dict[str, Any],
    *,
    campaign_dir: Path,
) -> EvalFtlSummary | None:
    """Build a summary from trajectory peak fields (no per-frame curve on disk)."""
    eval_id = record.get("eval")
    if not isinstance(eval_id, int):
        return None

    peaks = ftl_peaks_from_trajectory(record)
    if peaks is None:
        return None

    components = record.get("components") or {}
    episode = record.get("episode")
    if isinstance(episode, str) and episode:
        eval_dir = Path(episode)
    else:
        eval_dir = campaign_dir / f"eval_{eval_id:06d}"

    return EvalFtlSummary(
        eval_id=eval_id,
        eval_dir=eval_dir,
        status=record.get("status") if isinstance(record.get("status"), str) else None,
        score=_optional_float(record.get("score")),
        n_frames=peaks.n_frames,
        f_geo_peak=peaks.f_geo_peak,
        t_at_f_geo_peak=peaks.t_at_f_geo_peak,
        f_geo_final=0.0,
        f_op_peak=peaks.f_op_peak,
        t_at_f_op_peak=peaks.t_at_f_op_peak,
        f_op_final=0.0,
        max_local_speed_peak=peaks.max_local_speed_peak,
        t_at_max_speed=peaks.t_at_max_speed,
        superluminal_fraction_peak=peaks.superluminal_fraction_peak,
        t_at_superluminal_peak=peaks.t_at_superluminal_peak,
        ftl_lifetime_fraction=peaks.ftl_lifetime_fraction,
        op_lifetime_fraction=0.0,
        ftl_geo_timeavg=peaks.ftl_geo_timeavg or None,
        operational_ftl_geodesic=_optional_float(components.get("operational_ftl_geodesic")),
        structural_persistence=_optional_float(components.get("structural_persistence")),
        has_frame_curve=False,
    )


def load_campaign_ftl_summaries(
    campaign_dir: Path,
    *,
    include_trajectory: bool = True,
) -> list[EvalFtlSummary]:
    """Load FTL summaries for the campaign (trajectory + on-disk eval dirs)."""
    campaign_dir = campaign_dir.expanduser().resolve()
    status_by_eval = load_trajectory_status(campaign_dir)
    by_eval: dict[int, EvalFtlSummary] = {}

    if include_trajectory:
        for record in _iter_trajectory_records(campaign_dir):
            summary = summary_from_trajectory_record(record, campaign_dir=campaign_dir)
            if summary is not None:
                by_eval[summary.eval_id] = summary

    for eval_dir in sorted(campaign_dir.glob("eval_*")):
        if not eval_dir.is_dir():
            continue
        eval_id = _eval_id_from_dir(eval_dir)
        if eval_id is None:
            continue
        summary = load_eval_ftl_summary(
            eval_dir,
            status=status_by_eval.get(eval_id),
        )
        if summary is not None:
            by_eval[eval_id] = summary

    return [by_eval[k] for k in sorted(by_eval)]


def _sort_value(summary: EvalFtlSummary, key: FtlSortKey) -> float:
    if key is FtlSortKey.F_GEO_PEAK:
        return summary.f_geo_peak
    if key is FtlSortKey.F_OP_PEAK:
        return summary.f_op_peak
    if key is FtlSortKey.MAX_SPEED:
        return summary.max_local_speed_peak
    if key is FtlSortKey.SUPERLUMINAL:
        return summary.superluminal_fraction_peak
    if key is FtlSortKey.LIFETIME:
        return summary.ftl_lifetime_fraction
    if key is FtlSortKey.TIMEAVG:
        return summary.ftl_geo_timeavg or 0.0
    return summary.score if summary.score is not None else float("-inf")


def rank_summaries(
    summaries: Iterable[EvalFtlSummary],
    *,
    sort_by: FtlSortKey = FtlSortKey.F_GEO_PEAK,
    status: str | None = None,
    min_f_geo_peak: float = 0.0,
) -> list[EvalFtlSummary]:
    filtered = [
        s
        for s in summaries
        if s.f_geo_peak >= min_f_geo_peak
        and (status is None or s.status == status)
    ]
    return sorted(filtered, key=lambda s: _sort_value(s, sort_by), reverse=True)


def _pct(value: float) -> str:
    return f"{100.0 * value:.2f}%"


def _time(value: float | None) -> str:
    return f"t={value:.1f}" if value is not None else "t=?"


def _maybe_pct(value: float, *, available: bool) -> str:
    return _pct(value) if available else "n/a"


def format_summary_line(summary: EvalFtlSummary) -> str:
    score = f"{summary.score:>7.1f}" if summary.score is not None else "    n/a"
    status = summary.status or "?"
    curve = "curve" if summary.has_frame_curve else "traj"
    final = _maybe_pct(summary.f_geo_final, available=summary.has_frame_curve)
    speed = (
        f"{summary.max_local_speed_peak:.3f} @{_time(summary.t_at_max_speed):>8}"
        if summary.has_frame_curve
        else "n/a"
    )
    return (
        f"eval {summary.eval_id:>5}  score={score}  {status:<20}  [{curve}]  "
        f"f_geo peak={_pct(summary.f_geo_peak):>6} @{_time(summary.t_at_f_geo_peak):>8}  "
        f"final={final:>6}  "
        f"f_op peak={_pct(summary.f_op_peak):>6} @{_time(summary.t_at_f_op_peak):>8}  "
        f"speed max={speed}  "
        f"lifetime={summary.ftl_lifetime_fraction:>4.0%}  "
        f"timeavg={summary.ftl_geo_timeavg or 0.0:.4f}  "
        f"persist={summary.structural_persistence if summary.structural_persistence is not None else 0.0:.2f}"
    )


def format_campaign_report(
    summaries: Sequence[EvalFtlSummary],
    *,
    campaign_dir: Path,
    sort_by: FtlSortKey,
    top_n: int,
    total_in_campaign: int | None = None,
) -> str:
    with_curve = sum(1 for s in summaries if s.has_frame_curve)
    total = total_in_campaign if total_in_campaign is not None else len(summaries)
    lines = [
        f"campaign : {campaign_dir}",
        f"evals    : {len(summaries)} scored ({with_curve} with per-frame curve on disk)",
        f"total    : {total} in trajectory",
        f"sorted by: {sort_by.value}",
        "",
    ]
    if not summaries:
        lines.append("(no scored evals in trajectory or on-disk ftl_timeseries)")
        return "\n".join(lines)

    shown = list(summaries[:top_n])
    lines.extend(format_summary_line(s) for s in shown)
    if len(summaries) > top_n:
        lines.append(f"... {len(summaries) - top_n} more evals omitted")
    return "\n".join(lines)


def format_ftl_curve(summary: EvalFtlSummary) -> str:
    """Print every per-frame sample for one eval (shows mid-run alpha)."""
    ts = load_eval_timeseries(summary.eval_dir)
    if ts is None:
        return f"eval {summary.eval_id}: no per-frame ftl_timeseries on disk"

    lines = [
        f"eval {summary.eval_id}  n_frames={ts.n_frames}  "
        f"f_geo peak={_pct(ts.f_geo_peak)} @{_time(ts.t_at_f_geo_peak)}  "
        f"final={_pct(_final_value(ts.f_geo))}",
        "     t    f_geo  ok  f_op   speed  superlum",
    ]
    for i in range(ts.n_frames):
        trust = "Y" if ts.geo_trustworthy[i] else "N"
        marker = " *" if ts.t[i] == ts.t_at_f_geo_peak else ""
        lines.append(
            f"  {ts.t[i]:5.1f}  {_pct(ts.f_geo[i]):>6}  {trust}  "
            f"{_pct(ts.f_op[i]):>6}  {ts.max_local_speed[i]:5.3f}  "
            f"{ts.superluminal_fraction[i]:5.2f}{marker}"
        )
    return "\n".join(lines)
