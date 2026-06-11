"""MAP-Elites quality-diversity search for the Spacetime Failure Atlas.

The stated deliverable is an atlas of how geometries behave/fail, which is a
quality-diversity problem rather than single-optimum search.  MAP-Elites keeps
the best candidate in each cell of a behavior-descriptor grid, illuminating the
whole constructibility map in one campaign while resisting premature
convergence.

Behavior descriptors (both in [0, 1]):
  * ``ftl_benefit``  -- operational FTL on constraint-satisfying initial data
    (``operational_ftl_solved``), else evolved/t=0 shortcuts.
  * ``mechanism``    -- how the shortcut is produced (0 shift-warp .. 1 portal);
    separates warp / throat / portal families in the archive.

Elite quality is the full episode score.  The archive directly supplies the
diverse training set for a future Level-4 agent.
"""

from __future__ import annotations

import json
import math
import shutil
import threading
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from ..core.config import ExampleConfig, ExecutableConfig, resolve_example
from ..core.episode import write_json
from ..core.evaluation import Evaluation, evaluate_overrides
from .optimize import DEFAULT_SEARCH_SPACE, SearchDimension
from .trajectory_log import (
    format_trajectory_line,
    infer_trajectory_status,
    trajectory_flags_from_evaluation,
)
from .validation_tiers import (
    DEFAULT_TIER_CONFIG,
    Tier,
    TierConfig,
    build_survivors,
    convergence_signals,
    evaluate_tiers,
    survivor_front,
)


def _path_closeness_from_report(report: Mapping[str, Any] | None) -> float:
    if not report or not bool(report.get("reachable", False)):
        return 0.0
    t_min = report.get("t_min")
    t_flat = float(report.get("t_flat") or 0.0)
    if t_min is None or not math.isfinite(float(t_min)) or t_flat <= 0.0:
        return 0.0
    t_min_f = float(t_min)
    if t_min_f <= t_flat:
        return 1.0
    excess = (t_min_f - t_flat) / t_flat
    return max(0.0, 1.0 - excess / 0.12)


def _best_ftl_report(metrics: Mapping[str, Any] | None) -> Mapping[str, Any] | None:
    metrics = metrics or {}
    return (
        metrics.get("general_ftl_evolved")
        or metrics.get("general_ftl_solved")
        or metrics.get("general_ftl")
    )


def _solved_ftl_report(metrics: Mapping[str, Any] | None) -> Mapping[str, Any] | None:
    """Report from the constructed (solved) initial data.

    The superluminal region built by GRTresna largely decays under evolution, so
    ``general_ftl_evolved`` reads ~0 for almost every candidate.  The solved data
    retains the discriminating signal (observed superluminal_fraction 0 - 0.30 and
    max_local_speed ~0.95 - 1.32), so behaviour descriptors read it instead.
    """
    metrics = metrics or {}
    return (
        metrics.get("general_ftl_solved")
        or metrics.get("general_ftl_evolved")
        or metrics.get("general_ftl")
    )


# speed_horizon descriptor axes.  ``x`` maps the fastest evolved coordinate light
# speed across the c=1 threshold (cone-tilt strength); ``y`` maps the minimum
# null expansion so a horizon-free slice (min_theta_plus > 0) lands above 0.5 and
# a trapped surface (<= 0) below it.  This explicitly fills the "fast but no
# trapped surface" niche without needing a non-zero shift (unlike ``channel``).
_SPEED_HORIZON_C_FLOOR = 0.9
_SPEED_HORIZON_C_TARGET = 1.3
_SPEED_HORIZON_THETA_SCALE = 0.5

# speed_super descriptor axes, both read from the solved (constructed) initial
# data where the signal lives (the evolved report collapses to ~0).  ``x`` is the
# cone-tilt strength calibrated to the observed solved range (max_local_speed
# ~0.95-1.32) so it spreads across the bins instead of saturating.  ``y`` is the
# superluminal_fraction (share of the slice whose local light speed exceeds c=1),
# rescaled to its observed solved ceiling (~0.30) so the realistic range fills the
# grid instead of clustering in bin 0 on a raw [0, 1] scale.  Together they
# separate a localized cone-tilt from a widespread superluminal region.
_SPEED_SUPER_C_FLOOR = 0.95
_SPEED_SUPER_C_TARGET = 1.20
# ``superluminal_fraction`` is now measured with a 0.05 margin above c=1 (see
# ``SUPERLUMINAL_MARGIN`` in the general FTL probe), so it tracks the genuinely
# superluminal *area* (cone-tilted lobes at c ~ 1.08-1.18) rather than the broad
# gauge-shift background (~1.03) that previously saturated it to 1.0.  The
# noise-free fraction spreads over ~0.04-0.17 for the observed candidates, so
# the descriptor target is lowered to keep that range across the grid instead of
# clustering in the top bin.
_SPEED_SUPER_FRACTION_TARGET = 0.15


def _speed_tilt_axis(
    report: Mapping[str, Any] | None,
    *,
    floor: float,
    target: float,
) -> float:
    if not report:
        return 0.0
    c = float(report.get("max_local_speed") or 0.0)
    if not math.isfinite(c) or c <= floor:
        return 0.0
    span = target - floor
    return float(np.clip((c - floor) / span, 0.0, 1.0))


def _speed_axis_from_report(report: Mapping[str, Any] | None) -> float:
    return _speed_tilt_axis(
        report, floor=_SPEED_HORIZON_C_FLOOR, target=_SPEED_HORIZON_C_TARGET
    )


def _superluminal_axis(
    report: Mapping[str, Any] | None,
    *,
    target: float = 1.0,
) -> float:
    if not report:
        return 0.0
    frac = report.get("superluminal_fraction")
    if frac is None or not math.isfinite(float(frac)):
        return 0.0
    span = target if target > 0.0 else 1.0
    return float(np.clip(float(frac) / span, 0.0, 1.0))


def _horizon_free_axis(metrics: Mapping[str, Any] | None) -> float:
    metrics = metrics or {}
    collapse = metrics.get("collapse") or {}
    theta = collapse.get("min_theta_plus")
    if theta is None or not math.isfinite(float(theta)):
        return 0.5  # unknown -> neutral cell, neither flagged safe nor trapped
    scaled = 0.5 + float(theta) / (2.0 * _SPEED_HORIZON_THETA_SCALE)
    return float(np.clip(scaled, 0.0, 1.0))


def _descriptor_details(
    components: Mapping[str, float],
    metrics: Mapping[str, Any] | None = None,
    *,
    mode: str = "legacy",
) -> dict[str, float]:
    if mode == "speed_horizon":
        report = _best_ftl_report(metrics)
        speed = _speed_axis_from_report(report)
        horizon_free = _horizon_free_axis(metrics)
        c = float(report.get("max_local_speed")) if report and report.get("max_local_speed") is not None else float("nan")
        collapse = (metrics or {}).get("collapse") or {}
        theta = float(collapse.get("min_theta_plus")) if collapse.get("min_theta_plus") is not None else float("nan")
        return {
            "x": speed,
            "y": horizon_free,
            "speed_tilt": speed,
            "horizon_free": horizon_free,
            "max_local_speed": c,
            "min_theta_plus": theta,
            "ftl_persistence": float(components.get("ftl_persistence", 0.0)),
            "operational_ftl": float(components.get("operational_ftl", 0.0)),
        }

    if mode == "speed_super":
        report = _solved_ftl_report(metrics)
        speed = _speed_tilt_axis(
            report, floor=_SPEED_SUPER_C_FLOOR, target=_SPEED_SUPER_C_TARGET
        )
        super_frac = _superluminal_axis(
            report, target=_SPEED_SUPER_FRACTION_TARGET
        )
        raw_frac = (
            float(report.get("superluminal_fraction"))
            if report and report.get("superluminal_fraction") is not None
            else float("nan")
        )
        c = float(report.get("max_local_speed")) if report and report.get("max_local_speed") is not None else float("nan")
        return {
            "x": speed,
            "y": super_frac,
            "speed_tilt": speed,
            "superluminal_fraction": super_frac,
            "superluminal_fraction_raw": raw_frac,
            "max_local_speed": c,
            "ftl_persistence": float(components.get("ftl_persistence", 0.0)),
            "operational_ftl": float(components.get("operational_ftl", 0.0)),
        }

    if mode == "channel":
        report = _best_ftl_report(metrics)
        path_closeness = _path_closeness_from_report(report)
        precursor = float(np.clip(components.get("ftl_precursor", 0.0), 0.0, 1.0))
        shift = float(np.clip(components.get("shift_drive", 0.0), 0.0, 1.0))
        mechanism_balance = float(math.sqrt(max(precursor, 0.0) * max(shift, 0.0)))
        t_min = float(report.get("t_min")) if report and report.get("t_min") is not None else float("nan")
        t_flat = float(report.get("t_flat")) if report and report.get("t_flat") is not None else float("nan")
        return {
            "x": path_closeness,
            "y": mechanism_balance,
            "path_closeness": path_closeness,
            "mechanism_balance": mechanism_balance,
            "shift_drive": shift,
            "ftl_precursor": precursor,
            "channel_progress": float(components.get("channel_progress", 0.0)),
            "operational_ftl": float(components.get("operational_ftl", 0.0)),
            "t_min": t_min,
            "t_flat": t_flat,
        }

    ftl_benefit = float(
        np.clip(
            max(
                components.get("operational_ftl_solved", 0.0),
                components.get("operational_ftl", 0.0),
                components.get("ftl_precursor", 0.0),
                components.get("ftl_shortcut", 0.0),
            ),
            0.0,
            1.0,
        )
    )
    mechanism = float(np.clip(components.get("mechanism_descriptor", 0.0), 0.0, 1.0))
    if mechanism <= 0.0:
        mechanism = float(np.clip(1.0 - components.get("anec_condition", 1.0), 0.0, 1.0))
    return {
        "x": ftl_benefit,
        "y": mechanism,
        "ftl_benefit": ftl_benefit,
        "mechanism": mechanism,
    }


def _descriptors(
    components: Mapping[str, float],
    metrics: Mapping[str, Any] | None = None,
    *,
    mode: str = "legacy",
) -> tuple[float, float]:
    details = _descriptor_details(components, metrics, mode=mode)
    return details["x"], details["y"]


def _bin_index(value: float, bins: int) -> int:
    return int(min(bins - 1, max(0, math.floor(value * bins))))


@dataclass
class Elite:
    cell: tuple[int, int]
    score: float
    descriptors: tuple[float, float]
    params: dict[str, float]
    episode: str | None
    tier: int = int(Tier.REJECTED)
    tier_name: str = "rejected"
    objectives: dict[str, float] = field(default_factory=dict)
    descriptor_details: dict[str, float] = field(default_factory=dict)


@dataclass
class QDArchive:
    bins: int
    cells: dict[tuple[int, int], Elite] = field(default_factory=dict)

    def insert(self, elite: Elite) -> bool:
        existing = self.cells.get(elite.cell)
        if existing is None or elite.score > existing.score:
            self.cells[elite.cell] = elite
            return True
        return False

    @property
    def coverage(self) -> float:
        return len(self.cells) / float(self.bins * self.bins)

    @property
    def best(self) -> Elite | None:
        if not self.cells:
            return None
        return max(self.cells.values(), key=lambda e: e.score)

    def tier_histogram(self) -> dict[str, int]:
        from .validation_tiers import TIER_NAMES

        hist: dict[str, int] = {name: 0 for name in TIER_NAMES.values()}
        for e in self.cells.values():
            hist[e.tier_name] = hist.get(e.tier_name, 0) + 1
        return hist

    def to_dict(self) -> dict[str, Any]:
        return {
            "bins": self.bins,
            "coverage": self.coverage,
            "num_elites": len(self.cells),
            "best_score": self.best.score if self.best else None,
            "tier_histogram": self.tier_histogram(),
            "cells": [
                {
                    "cell": list(e.cell),
                    "score": e.score,
                    "descriptors": list(e.descriptors),
                    "descriptor_details": e.descriptor_details,
                    "tier": e.tier,
                    "tier_name": e.tier_name,
                    "objectives": e.objectives,
                    "params": e.params,
                    "episode": e.episode,
                }
                for e in sorted(self.cells.values(), key=lambda e: (-e.tier, -e.score))
            ],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> QDArchive:
        archive = cls(bins=int(data["bins"]))
        for cell_data in data.get("cells", []):
            cell = tuple(cell_data["cell"])
            archive.cells[cell] = Elite(
                cell=cell,
                score=float(cell_data["score"]),
                descriptors=tuple(cell_data["descriptors"]),
                params=dict(cell_data["params"]),
                episode=cell_data.get("episode"),
                tier=int(cell_data.get("tier", int(Tier.REJECTED))),
                tier_name=str(cell_data.get("tier_name", "rejected")),
                objectives=dict(cell_data.get("objectives", {})),
                descriptor_details=dict(cell_data.get("descriptor_details", {})),
            )
        return archive


def _load_trajectory_records(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    records: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            records.append(json.loads(line))
    return records


def _prune_eval_dirs(
    qd_dir: Path,
    records: Sequence[Mapping[str, Any]],
    *,
    keep_top: int,
    protect_eval_ids: set[int] | None = None,
) -> int:
    """Delete completed eval dirs outside the top-scoring ``keep_top`` records.

    ``trajectory.jsonl`` remains the source of truth; this only reclaims bulky
    per-eval folders after their score/metrics have been recorded.  Current
    batch ids are protected so live consumers/writers are not raced.
    """
    if keep_top <= 0:
        return 0

    scored: list[tuple[float, int]] = []
    for rec in records:
        score = rec.get("score")
        eval_id = rec.get("eval")
        if not isinstance(score, (int, float)) or eval_id is None:
            continue
        try:
            scored.append((float(score), int(eval_id)))
        except (TypeError, ValueError):
            continue

    keep_ids = {eval_id for _score, eval_id in sorted(scored, reverse=True)[:keep_top]}
    protected = set(protect_eval_ids or set())
    scored_ids = {eval_id for _score, eval_id in scored}
    deleted = 0

    for eval_dir in qd_dir.glob("eval_*"):
        if not eval_dir.is_dir():
            continue
        suffix = eval_dir.name.removeprefix("eval_")
        if not suffix.isdigit():
            continue
        eval_id = int(suffix)
        if eval_id in keep_ids or eval_id in protected:
            continue
        if eval_id not in scored_ids:
            continue
        shutil.rmtree(eval_dir, ignore_errors=True)
        deleted += 1

    if deleted:
        print(
            f"[qd] pruned {deleted} eval dirs; kept top {keep_top} scored "
            f"+ protected {sorted(protected)}",
            flush=True,
        )
    return deleted


def _iterations_for_target_evals(
    *,
    target_evals: int,
    batch: int,
    completed_evals: int = 0,
    include_init: bool = True,
) -> int:
    """Return MAP-Elites batch count needed to reach ``target_evals``."""
    if target_evals <= 0:
        return 0
    if completed_evals > 0:
        remaining = target_evals - completed_evals
        if remaining <= 0:
            return 0
        return int(math.ceil(remaining / batch))
    init = batch if include_init else 0
    if target_evals <= init:
        return 0
    return int(math.ceil((target_evals - init) / batch))


# Fraction of each batch drawn by mutating an existing elite (vs. random
# exploration).  Raised from the original 0.7 so the search exploits the
# narrow feasible basin once it is seeded instead of repeatedly re-sampling
# pathological corners.
_ELITE_FRACTION = 0.85


def _reflect(value: float, lower: float, upper: float) -> float:
    """Fold ``value`` back into ``[lower, upper]`` by boundary reflection.

    Unlike hard clipping (which piles probability mass exactly on the bounds,
    where the GRTresna solve tends to diverge), reflection keeps the mutated
    value in the interior, preserving a smooth proposal density.
    """
    span = upper - lower
    if span <= 0.0:
        return lower
    t = (value - lower) % (2.0 * span)
    if t < 0.0:
        t += 2.0 * span
    if t > span:
        t = 2.0 * span - t
    return lower + t


def _feasible_bounds(
    dims: Sequence[SearchDimension],
    elites: Sequence[Elite],
) -> list[tuple[float, float]]:
    """Per-dimension (min, max) over the current elite params.

    Falls back to the full search bound for any dimension with no spread, so a
    single (or empty) elite set degrades gracefully to the full space.
    """
    bounds: list[tuple[float, float]] = []
    for d in dims:
        vals = [
            float(e.params[d.param_key])
            for e in elites
            if d.param_key in e.params
        ]
        if len(vals) >= 2 and max(vals) > min(vals):
            bounds.append((min(vals), max(vals)))
        else:
            bounds.append((d.lower, d.upper))
    return bounds


def _sample_random(dims: Sequence[SearchDimension], rng: np.random.Generator) -> list[float]:
    return [float(rng.uniform(d.lower, d.upper)) for d in dims]


def _sample_feasible_box(
    dims: Sequence[SearchDimension],
    bounds: Sequence[tuple[float, float]],
    rng: np.random.Generator,
) -> list[float]:
    """Uniform sample inside the bounding box of the feasible elites."""
    return [float(rng.uniform(lo, hi)) for (lo, hi) in bounds]


def _mutate_elite(
    elite: Elite,
    dims: Sequence[SearchDimension],
    rng: np.random.Generator,
    *,
    sigma: float = 0.15,
) -> list[float]:
    x = []
    for d in dims:
        base = float(elite.params.get(d.param_key, d.center))
        step = rng.normal(0.0, sigma * max(d.range, 1.0e-9))
        x.append(_reflect(base + step, d.lower, d.upper))
    return x


def _vector_to_overrides(
    x: Sequence[float],
    dims: Sequence[SearchDimension],
    base: Mapping[str, Any],
) -> dict[str, Any]:
    overrides = dict(base)
    for xi, d in zip(x, dims):
        overrides[d.param_key] = float(min(d.upper, max(d.lower, xi)))
    return overrides


def run_qd_search(
    *,
    runs_dir: Path,
    executable: ExecutableConfig | None = None,
    iterations: int = 10,
    batch_size: int | None = None,
    bins: int = 8,
    init_random: int | None = None,
    seed: int | None = None,
    base_overrides: Mapping[str, Any] | None = None,
    search_space: Sequence[SearchDimension] | None = None,
    template: Path | None = None,
    example: ExampleConfig | str = "RadialRecipe",
    name: str | None = None,
    dry_run: bool = False,
    constrained: bool = True,
    phantom: bool = True,
    use_preflight: bool = True,
    cuda_devices: str | None = None,
    gpu_ids: Sequence[int] | None = None,
    check_params: bool = True,
    score_weights: Mapping[str, float] | None = None,
    ftl_L: float | None = None,
    consume_plotfiles: bool = True,
    consumer_radii: Sequence[float] = (4.0, 8.0),
    consumer_keep_last: int = 1,
    tier_config: TierConfig = DEFAULT_TIER_CONFIG,
    survivor_min_tier: int = int(Tier.OPERATIONAL),
    descriptor_mode: str = "legacy",
    objective_mode: str = "weighted",
    grtresna: bool = False,
    grtresna_config: Any | None = None,
    grtresna_solved_ftl_gate: bool = False,
    solved_ftl_gate_config: Any | None = None,
    grtresna_convergence_config: Any | None = None,
    grtresna_postload_gate: bool = False,
    postload_gate_config: Any | None = None,
    resume: bool = False,
    target_evals: int | None = None,
    seed_overrides: Sequence[Mapping[str, Any]] | None = None,
    keep_top_eval_dirs: int = 0,
) -> QDArchive:
    example_cfg = example if isinstance(example, ExampleConfig) else resolve_example(example)
    dims = list(search_space or DEFAULT_SEARCH_SPACE)
    tpl = template or example_cfg.template
    base = dict(base_overrides or {})
    base.setdefault("N_full", 64)
    base.setdefault("max_level", 2)
    base.setdefault("stop_time", 2.0)
    base.setdefault("plot_interval", 10)
    base.setdefault("checkpoint_interval", -1)
    base.setdefault("dt_multiplier", 0.02)
    target_stop_time = float(base["stop_time"])

    rng = np.random.default_rng(seed)
    batch = batch_size or (len(gpu_ids) if gpu_ids else 4)
    n_init = init_random if init_random is not None else batch

    if name is None:
        if resume:
            raise ValueError("resume requires --name pointing at an existing qd_<timestamp> directory")
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        name = f"qd_{timestamp}"
    qd_dir = (runs_dir / name).expanduser().resolve()

    trajectory_path = qd_dir / "trajectory.jsonl"
    conv_history: list[dict[str, Any]] = []
    completed_evals = 0

    if resume:
        if not qd_dir.is_dir():
            raise FileNotFoundError(f"resume directory not found: {qd_dir}")
        archive_path = qd_dir / "archive.json"
        if not archive_path.exists():
            raise FileNotFoundError(f"resume requires archive.json in {qd_dir}")
        archive = QDArchive.from_dict(json.loads(archive_path.read_text(encoding="utf-8")))
        if archive.bins != bins:
            raise ValueError(
                f"resume bins mismatch: archive has {archive.bins}, requested {bins}"
            )
        trajectory = _load_trajectory_records(trajectory_path)
        if trajectory:
            completed_evals = max(int(rec["eval"]) for rec in trajectory if "eval" in rec)
        validation_path = qd_dir / "validation.json"
        if validation_path.exists():
            validation = json.loads(validation_path.read_text(encoding="utf-8"))
            conv_history = list(validation.get("history", []))
        metadata_path = qd_dir / "metadata.json"
        if metadata_path.exists() and seed is None:
            meta = json.loads(metadata_path.read_text(encoding="utf-8"))
            seed = meta.get("seed")
        if target_evals is not None:
            iterations = _iterations_for_target_evals(
                target_evals=target_evals,
                batch=batch,
                completed_evals=completed_evals,
            )
        print(
            f"[qd] resume {qd_dir.name}: completed={completed_evals} "
            f"elites={len(archive.cells)} planned_batches={iterations}"
            + (f" target_evals={target_evals}" if target_evals is not None else "")
        )
        if target_evals is not None and iterations == 0:
            print(f"[qd] target_evals={target_evals} already reached; nothing to do.")
            return archive
    else:
        qd_dir.mkdir(parents=True, exist_ok=False)
        archive = QDArchive(bins=bins)
        trajectory = []
        if target_evals is not None:
            iterations = _iterations_for_target_evals(
                target_evals=target_evals,
                batch=batch,
            )

    eval_counter = [completed_evals]
    counter_lock = threading.Lock()
    trajectory_lock = threading.Lock()

    metadata = {
        "mode": "qd_map_elites",
        "example": example_cfg.name,
        "iterations": iterations,
        "batch_size": batch,
        "bins": bins,
        "seed": seed,
        "gpu_ids": list(gpu_ids) if gpu_ids else None,
        "descriptor_mode": descriptor_mode,
        "base_overrides": base,
        "search_space": [
            {"key": d.param_key, "lower": d.lower, "upper": d.upper} for d in dims
        ],
        "grtresna": grtresna,
        "grtresna_solved_ftl_gate": grtresna_solved_ftl_gate,
        "keep_top_eval_dirs": keep_top_eval_dirs,
        "grtresna_max_ham_pct": getattr(
            grtresna_convergence_config, "max_ham_pct", 5.0,
        ),
        "grtresna_max_mom_pct": getattr(
            grtresna_convergence_config, "max_mom_pct", 5.0,
        ),
    }
    if target_evals is not None:
        metadata["target_evals"] = target_evals
    if resume:
        metadata["resumed_at"] = datetime.now(timezone.utc).isoformat()
        metadata["resume_from_evals"] = completed_evals
        metadata["resume_planned_batches"] = iterations
        existing_meta_path = qd_dir / "metadata.json"
        if existing_meta_path.exists():
            existing = json.loads(existing_meta_path.read_text(encoding="utf-8"))
            metadata["created_at"] = existing.get(
                "created_at", datetime.now(timezone.utc).isoformat()
            )
        else:
            metadata["created_at"] = datetime.now(timezone.utc).isoformat()
    else:
        metadata["created_at"] = datetime.now(timezone.utc).isoformat()
    write_json(qd_dir / "metadata.json", metadata)

    def _append_live_trajectory(record: Mapping[str, Any]) -> None:
        with trajectory_lock:
            with trajectory_path.open("a", encoding="utf-8") as fh:
                fh.write(format_trajectory_line(record))

    def _write_trajectory() -> None:
        with trajectory_lock:
            with trajectory_path.open("w", encoding="utf-8") as fh:
                for rec in trajectory:
                    fh.write(format_trajectory_line(rec))

    def _record_result(
        *,
        idx: int,
        x: list[float],
        res: Evaluation | None,
        insert_archive: bool,
    ) -> dict[str, Any]:
        overrides = _vector_to_overrides(x, dims, base)
        if res is None or res.preflight_rejected:
            record: dict[str, Any] = {
                "eval": idx,
                "preflight_rejected": True,
                "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
            }
            if res is not None:
                record.update(trajectory_flags_from_evaluation(res))
            record["status"] = infer_trajectory_status(record)
            return record

        descriptor_details = _descriptor_details(
            res.components, res.metrics, mode=descriptor_mode,
        )
        d1, d2 = descriptor_details["x"], descriptor_details["y"]
        cell = (_bin_index(d1, bins), _bin_index(d2, bins))
        params = {d.param_key: float(overrides[d.param_key]) for d in dims}
        assessment = evaluate_tiers(res.components, metrics=res.metrics, config=tier_config)

        flags = trajectory_flags_from_evaluation(res)
        status = infer_trajectory_status({**flags, "components": res.components})

        # Only GPU-evolved candidates carry meaningful behavior descriptors;
        # constraint-rejected / failed candidates all collapse onto the same
        # degenerate cell and would pollute coverage, so they are logged to the
        # trajectory but kept out of the archive grid.
        improved = None
        if insert_archive and status == "gpu_ok":
            elite = Elite(
                cell=cell,
                score=res.score,
                descriptors=(d1, d2),
                params=params,
                episode=res.episode_path,
                tier=assessment.tier,
                tier_name=assessment.tier_name,
                objectives=assessment.objectives,
                descriptor_details=descriptor_details,
            )
            improved = archive.insert(elite)

        record = {
            "eval": idx,
            "score": res.score,
            "cell": list(cell),
            "descriptors": [d1, d2],
            "descriptor_details": descriptor_details,
            "tier": assessment.tier,
            "tier_name": assessment.tier_name,
            "improved": improved,
            "episode": res.episode_path,
            "components": res.components,
            "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
        }
        record.update(flags)
        record["status"] = status
        return record

    def _evaluate_batch(vectors: list[list[float]]) -> list[tuple[int, list[float], Evaluation | None] | None]:
        results: list[tuple[int, list[float], Evaluation | None] | None] = [None] * len(vectors)

        def _eval_one(i: int, x: list[float], gpu: str | None) -> None:
            with counter_lock:
                eval_counter[0] += 1
                idx = eval_counter[0]
            overrides = _vector_to_overrides(x, dims, base)
            res = evaluate_overrides(
                overrides,
                out_dir=qd_dir,
                name=f"eval_{idx:06d}",
                example=example_cfg,
                template=tpl,
                executable=executable,
                constrained=constrained,
                phantom=phantom,
                use_preflight=use_preflight,
                cuda_devices=gpu if gpu is not None else cuda_devices,
                check_params=check_params,
                dry_run=dry_run,
                target_stop_time=target_stop_time,
                score_weights=score_weights,
                objective_mode=objective_mode,
                ftl_L=ftl_L,
                consume_plotfiles=consume_plotfiles,
                consumer_radii=consumer_radii,
                consumer_keep_last=consumer_keep_last,
                grtresna=grtresna,
                grtresna_base=grtresna_config,
                grtresna_solved_ftl_gate=grtresna_solved_ftl_gate,
                solved_ftl_gate_config=solved_ftl_gate_config,
                grtresna_convergence_config=grtresna_convergence_config,
                grtresna_postload_gate=grtresna_postload_gate,
                postload_gate_config=postload_gate_config,
            )
            results[i] = (idx, x, res)
            _append_live_trajectory(
                _record_result(idx=idx, x=x, res=res, insert_archive=False)
            )

        if gpu_ids and len(gpu_ids) > 1 and not dry_run:
            threads = []
            for i, x in enumerate(vectors):
                gpu = str(gpu_ids[i % len(gpu_ids)])
                t = threading.Thread(target=_eval_one, args=(i, x, gpu))
                t.start()
                threads.append(t)
            for t in threads:
                t.join()
        else:
            for i, x in enumerate(vectors):
                _eval_one(i, x, cuda_devices)
        return results

    def _ingest(results: list[tuple[int, list[float], Evaluation | None] | None]) -> None:
        protected_eval_ids: set[int] = set()
        for item in results:
            if item is None:
                continue
            idx, x, res = item
            protected_eval_ids.add(idx)
            trajectory.append(
                _record_result(idx=idx, x=x, res=res, insert_archive=True)
            )
        _write_trajectory()
        _prune_eval_dirs(
            qd_dir,
            trajectory,
            keep_top=int(keep_top_eval_dirs),
            protect_eval_ids=protected_eval_ids,
        )

    def _write_validation() -> dict[str, Any]:
        """Assess the archive's survivor front and archive-convergence signal."""
        items = [
            {
                "label": f"cell_{e.cell[0]}_{e.cell[1]}",
                "tier": e.tier,
                "score": e.score,
                "objectives": e.objectives,
                "episode": e.episode,
            }
            for e in archive.cells.values()
        ]
        survivors = build_survivors(items, min_tier=survivor_min_tier)
        front = survivor_front(survivors)
        best = archive.best
        snapshot = {
            "iteration": len(conv_history),
            "coverage": archive.coverage,
            "best_score": best.score if best else 0.0,
            "front_labels": [s.label for s in front],
        }
        conv_history.append(snapshot)
        signals = convergence_signals(conv_history)
        write_json(qd_dir / "validation.json", {
            "tier_histogram": archive.tier_histogram(),
            "survivor_min_tier": survivor_min_tier,
            "num_survivors": len(survivors),
            "survivor_front": [
                {
                    "label": s.label,
                    "tier": s.tier,
                    "tier_name": next(
                        (e.tier_name for e in archive.cells.values()
                         if f"cell_{e.cell[0]}_{e.cell[1]}" == s.label),
                        "",
                    ),
                    "score": s.score,
                    "objectives": s.objectives,
                    "episode": s.episode,
                }
                for s in front
            ],
            "convergence": signals,
            "history": conv_history,
        })
        return signals

    start_iter = len(conv_history)
    n_seed = 0 if resume else len(seed_overrides or [])
    total_target = (
        target_evals
        if target_evals is not None
        else completed_evals + (0 if resume else n_init) + n_seed + iterations * batch
    )
    print(
        f"[qd] MAP-Elites: {len(dims)}D, bins={bins}x{bins}, batch={batch}, "
        f"iters={iterations}, GPUs={gpu_ids or cuda_devices}, "
        f"evals={completed_evals}->{total_target}"
    )

    if not resume:
        # Warm-start: evaluate provided seed parameter sets (e.g. promoted
        # survivors) first so MAP-Elites mutates around known-good basins
        # instead of starting purely random.  Missing dims fall back to the
        # dimension centre (so search-space additions like grtresna_shift_seed
        # start neutral and are then explored by mutation).
        seed_vectors = [
            [float(s.get(d.param_key, d.center)) for d in dims]
            for s in (seed_overrides or [])
        ]
        init_vectors = seed_vectors + [
            _sample_random(dims, rng) for _ in range(n_init)
        ]
        _ingest(_evaluate_batch(init_vectors))
        write_json(qd_dir / "archive.json", archive.to_dict())
        _write_validation()

    for it in range(1, iterations + 1):
        vectors: list[list[float]] = []
        elites = list(archive.cells.values())
        feasible_bounds = _feasible_bounds(dims, elites) if elites else None
        for _ in range(batch):
            if elites and rng.random() < _ELITE_FRACTION:
                parent = elites[int(rng.integers(len(elites)))]
                vectors.append(_mutate_elite(parent, dims, rng))
            elif feasible_bounds is not None:
                # Bias exploration toward the box that already produced feasible
                # elites instead of the full (mostly pathological) search space.
                vectors.append(_sample_feasible_box(dims, feasible_bounds, rng))
            else:
                vectors.append(_sample_random(dims, rng))
        _ingest(_evaluate_batch(vectors))

        write_json(qd_dir / "archive.json", archive.to_dict())
        signals = _write_validation()
        best = archive.best
        n_front = len(conv_history[-1]["front_labels"])
        global_iter = start_iter + it
        print(
            f"[qd] iter {it}/{iterations} (total {global_iter}): "
            f"elites={len(archive.cells)} coverage={archive.coverage:.2f} "
            f"best={best.score if best else float('nan'):.4f} "
            f"front(T>={survivor_min_tier})={n_front} "
            f"converged={signals.get('converged')} evals={eval_counter[0]}/{total_target}"
        )
        if signals.get("converged"):
            print(f"[qd] archive converged at iter {it}: "
                  f"coverage_delta={signals.get('coverage_delta'):.3f} "
                  f"best_delta={signals.get('best_score_delta'):.4f} "
                  f"front_stable={signals.get('front_stable')}")

    write_json(qd_dir / "archive.json", archive.to_dict())
    final_hist = archive.tier_histogram()
    print(f"[qd] Done. elites={len(archive.cells)} coverage={archive.coverage:.2f} "
          f"tiers={final_hist} dir={qd_dir}")
    return archive
