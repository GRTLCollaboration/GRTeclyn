"""Prefix (in-flight) episode scoring from streaming diagnostics.

After each plotfile the consumer appends a row to ``ftl_timeseries.dat``.
This module re-reads the prefix of that stream plus collapse/constraint
series truncated at the same simulation time, runs ``score_episode``, and
appends one JSON line to ``small_data/score_timeseries.jsonl``.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from ..diagnostics.collapse import read_collapse_metrics
from ..diagnostics.ftl_timeseries import _aggregate_ftl_frames
from ..io.dat import numeric_rows
from ..probes.ftl.analytic import compute_ftl_metrics, load_overrides_from_episode
from ..probes.ftl.general import GeneralFtlReport
from ..probes.ftl.geodesic import H_REL_TOL, GeodesicFtlReport
from ..score import domain_half_width_for_episode, score_episode
from ..types.diagnostics import ConstraintMetrics, FtlPersistenceMetrics
from ..types.episode import EpisodeMetrics
from .context import build_episode_context


def _rows_up_to(path: Path, max_time: float, min_columns: int) -> list[list[float]]:
    return [row for row in numeric_rows(path, min_columns) if row[0] <= max_time + 1.0e-9]


def _write_prefix_dat(rows: list[list[float]], dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write("  ".join(f"{value:.10e}" for value in row) + "\n")


def read_prefix_collapse_metrics(path: Path, max_time: float):
    rows = _rows_up_to(path, max_time, 4)
    if not rows:
        return None
    tmp = path.parent / ".incremental_collapse_prefix.dat"
    _write_prefix_dat(rows, tmp)
    try:
        return read_collapse_metrics(tmp)
    finally:
        tmp.unlink(missing_ok=True)


def read_prefix_constraint_metrics(path: Path, max_time: float) -> ConstraintMetrics | None:
    rows = _rows_up_to(path, max_time, 3)
    if not rows:
        return None
    rho_rows = [row for row in rows if len(row) >= 6]
    min_rho_req = min((row[3] for row in rho_rows), default=None)
    max_rho_req = max((row[4] for row in rho_rows), default=None)
    max_int_neg = max((row[5] for row in rho_rows), default=None)
    final_peak_rho_req = rho_rows[-1][4] if rho_rows else None
    return ConstraintMetrics(
        final_time=rows[-1][0],
        max_hamiltonian_l2=max(row[1] for row in rows),
        max_momentum_l2=max(row[2] for row in rows),
        final_hamiltonian_l2=rows[-1][1],
        final_momentum_l2=rows[-1][2],
        min_rho_required=min_rho_req,
        max_rho_required=max_rho_req,
        integral_negative_rho=max_int_neg,
        final_peak_rho_required=final_peak_rho_req,
    )


def read_prefix_ftl_timeseries(path: Path, max_time: float):
    rows = _rows_up_to(path, max_time, 12)
    if not rows:
        return None
    return _aggregate_ftl_frames(
        t=tuple(float(r[0]) for r in rows),
        f_op=tuple(float(r[1]) for r in rows),
        f_geo=tuple(float(r[2]) for r in rows),
        geo_trustworthy=tuple(bool(round(r[3])) for r in rows),
        max_local_speed=tuple(float(r[4]) for r in rows),
        superluminal_fraction=tuple(float(r[5]) for r in rows),
        structure_coherence=tuple(float(r[7]) for r in rows),
        max_h_rel_drift=tuple(float(r[11]) for r in rows),
    )


def _general_ftl_from_timeseries_row(row: list[float]) -> GeneralFtlReport:
    coherence = float(row[7]) if math.isfinite(row[7]) else None
    return GeneralFtlReport(
        f_op=float(row[1]),
        t_min=None,
        t_flat=1.0,
        max_local_speed=float(row[4]),
        superluminal_fraction=float(row[5]),
        path_offaxis=False,
        reachable=bool(round(row[8])),
        notes=(),
        max_shift=float(row[6]),
        structure_coherence=coherence,
    )


def _geodesic_from_timeseries_row(row: list[float]) -> GeodesicFtlReport | None:
    n_rays = int(round(row[9]))
    n_reached = int(round(row[10]))
    if n_rays <= 0:
        return None
    max_h_rel = float(row[11])
    return GeodesicFtlReport(
        f_geo=float(row[2]),
        t_min=None,
        t_flat=1.0,
        n_rays=n_rays,
        n_reached=n_reached,
        max_h_drift=0.0,
        h_quality_ok=max_h_rel <= H_REL_TOL,
        max_h_rel_drift=max_h_rel,
        notes=(),
    )


def _ftl_persistence_from_rows(rows: list[list[float]]) -> FtlPersistenceMetrics | None:
    if len(rows) < 2:
        return None
    tail = rows[-5:]
    f_ops = [float(r[1]) for r in tail]
    speeds = [float(r[4]) for r in tail]
    shifts = [float(r[6]) for r in tail]
    return FtlPersistenceMetrics(
        n_samples=len(tail),
        f_op_min=float(min(f_ops)),
        f_op_median=float(np.median(f_ops)),
        f_op_last=float(f_ops[-1]),
        max_local_speed_min=float(min(speeds)),
        max_shift_max=float(max(shifts)),
    )


@dataclass
class _StaticFtlCache:
    ftl: Any = None
    general_ftl: Any = None
    general_ftl_solved: Any = None
    mechanism_descriptor: float | None = None


@dataclass
class IncrementalScoreWriter:
    """Append prefix scores to ``small_data/score_timeseries.jsonl``."""

    episode_dir: Path
    objective_mode: str = "weighted"
    target_stop_time: float | None = None
    ftl_L: float | None = None
    score_weights: Mapping[str, float] | None = None
    out_path: Path | None = None
    _static: _StaticFtlCache | None = field(default=None, init=False, repr=False)
    _last_t: float | None = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        self.episode_dir = Path(self.episode_dir)
        if self.out_path is None:
            self.out_path = self.episode_dir / "small_data" / "score_timeseries.jsonl"

    def _load_static(self) -> _StaticFtlCache:
        if self._static is not None:
            return self._static
        ctx = build_episode_context(self.episode_dir, ftl_L=self.ftl_L)
        overrides = load_overrides_from_episode(self.episode_dir)
        ftl = general_ftl = general_ftl_solved = mechanism_descriptor = None
        if overrides:
            try:
                ftl = compute_ftl_metrics(overrides, L=ctx.ftl_L)
            except Exception:
                pass
            try:
                from ..probes.ftl.general import compute_general_ftl

                general_ftl = compute_general_ftl(overrides, L=ctx.ftl_L, n=97)
            except Exception:
                pass
        if ctx.gridinit_path.is_file():
            try:
                from ..probes.ftl.solved import compute_solved_geometry_ftl
                from ...search.solved_ftl_gate import DEFAULT_SOLVED_FTL_GATE_POLICY

                solved = compute_solved_geometry_ftl(ctx.gridinit_path, L=ctx.ftl_L)
                if solved is not None and not DEFAULT_SOLVED_FTL_GATE_POLICY.is_degenerate(solved):
                    general_ftl_solved = solved.operational
                    mechanism_descriptor = solved.mechanisms.mechanism_descriptor
            except Exception:
                pass
        self._static = _StaticFtlCache(
            ftl=ftl,
            general_ftl=general_ftl,
            general_ftl_solved=general_ftl_solved,
            mechanism_descriptor=mechanism_descriptor,
        )
        return self._static

    def build_metrics_at(self, at_time: float) -> EpisodeMetrics | None:
        ctx = build_episode_context(self.episode_dir, ftl_L=self.ftl_L)
        ts_rows = _rows_up_to(ctx.ftl_timeseries_path, at_time, 12)
        if not ts_rows:
            return None

        ftl_timeseries = read_prefix_ftl_timeseries(ctx.ftl_timeseries_path, at_time)
        collapse = read_prefix_collapse_metrics(ctx.collapse_path, at_time)
        constraints = read_prefix_constraint_metrics(ctx.constraint_path, at_time)
        last_row = ts_rows[-1]
        general_ftl_evolved = _general_ftl_from_timeseries_row(last_row)
        geodesic_ftl = _geodesic_from_timeseries_row(last_row)
        ftl_persistence = _ftl_persistence_from_rows(ts_rows)
        static = self._load_static()

        return EpisodeMetrics(
            collapse=collapse,
            constraints=constraints,
            stability=None,
            comoving=None,
            ftl=static.ftl,
            general_ftl=static.general_ftl,
            general_ftl_solved=static.general_ftl_solved,
            mechanism_descriptor=static.mechanism_descriptor,
            general_ftl_evolved=general_ftl_evolved,
            geodesic_ftl=geodesic_ftl,
            ftl_persistence=ftl_persistence,
            ftl_timeseries=ftl_timeseries,
            termination_reason="incremental",
        )

    def append(self, at_time: float) -> dict[str, Any] | None:
        if self._last_t is not None and at_time <= self._last_t + 1.0e-12:
            return None
        metrics = self.build_metrics_at(at_time)
        if metrics is None:
            return None

        overrides = load_overrides_from_episode(self.episode_dir) or {}
        score = score_episode(
            metrics,
            target_stop_time=self.target_stop_time,
            weights=self.score_weights,
            objective_mode=self.objective_mode,
            domain_half_width=domain_half_width_for_episode(self.episode_dir, overrides),
        )
        ts = metrics.ftl_timeseries
        record: dict[str, Any] = {
            "t": at_time,
            "score": score.total,
            "objective_mode": self.objective_mode,
            "components": dict(score.components),
            "f_geo": float(ts.f_geo[-1]) if ts and ts.f_geo else 0.0,
            "f_geo_peak": float(ts.f_geo_peak) if ts else 0.0,
            "max_local_speed": float(ts.max_local_speed[-1]) if ts and ts.max_local_speed else 0.0,
            "horizon_penalty": score.components.get("horizon_penalty", 0.0),
            "operational_ftl_geodesic": score.components.get("operational_ftl_geodesic", 0.0),
        }
        self.out_path.parent.mkdir(parents=True, exist_ok=True)
        with self.out_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record) + "\n")
        self._last_t = at_time
        return record


def parse_score_weights_env(raw: str | None) -> dict[str, float] | None:
    if not raw:
        return None
    weights: dict[str, float] = {}
    for token in raw.split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        weights[key.strip()] = float(value.strip())
    return weights or None
