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


def _descriptor_details(
    components: Mapping[str, float],
    metrics: Mapping[str, Any] | None = None,
    *,
    mode: str = "legacy",
) -> dict[str, float]:
    if mode == "channel":
        metrics = metrics or {}
        report = (
            metrics.get("general_ftl_evolved")
            or metrics.get("general_ftl_solved")
            or metrics.get("general_ftl")
        )
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


def _sample_random(dims: Sequence[SearchDimension], rng: np.random.Generator) -> list[float]:
    return [float(rng.uniform(d.lower, d.upper)) for d in dims]


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
        x.append(float(min(d.upper, max(d.lower, base + step))))
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
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        name = f"qd_{timestamp}"
    qd_dir = (runs_dir / name).expanduser().resolve()
    qd_dir.mkdir(parents=True, exist_ok=False)

    archive = QDArchive(bins=bins)
    eval_counter = [0]
    counter_lock = threading.Lock()
    trajectory_lock = threading.Lock()
    trajectory_path = qd_dir / "trajectory.jsonl"
    trajectory: list[dict[str, Any]] = []

    write_json(qd_dir / "metadata.json", {
        "created_at": datetime.now(timezone.utc).isoformat(),
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
        "grtresna_max_ham_pct": getattr(
            grtresna_convergence_config, "max_ham_pct", 5.0,
        ),
        "grtresna_max_mom_pct": getattr(
            grtresna_convergence_config, "max_mom_pct", 5.0,
        ),
    })

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
        improved = None
        if insert_archive:
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
        record.update(trajectory_flags_from_evaluation(res))
        record["status"] = infer_trajectory_status(record)
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
        for item in results:
            if item is None:
                continue
            idx, x, res = item
            trajectory.append(
                _record_result(idx=idx, x=x, res=res, insert_archive=True)
            )
        _write_trajectory()

    conv_history: list[dict[str, Any]] = []

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

    print(f"[qd] MAP-Elites: {len(dims)}D, bins={bins}x{bins}, batch={batch}, "
          f"iters={iterations}, GPUs={gpu_ids or cuda_devices}")

    # Initial random fill.
    init_vectors = [_sample_random(dims, rng) for _ in range(n_init)]
    _ingest(_evaluate_batch(init_vectors))
    write_json(qd_dir / "archive.json", archive.to_dict())
    _write_validation()

    for it in range(1, iterations + 1):
        vectors: list[list[float]] = []
        elites = list(archive.cells.values())
        for _ in range(batch):
            if elites and rng.random() < 0.7:
                parent = elites[int(rng.integers(len(elites)))]
                vectors.append(_mutate_elite(parent, dims, rng))
            else:
                vectors.append(_sample_random(dims, rng))
        _ingest(_evaluate_batch(vectors))

        write_json(qd_dir / "archive.json", archive.to_dict())
        signals = _write_validation()
        best = archive.best
        n_front = len(conv_history[-1]["front_labels"])
        print(f"[qd] iter {it}/{iterations}: elites={len(archive.cells)} "
              f"coverage={archive.coverage:.2f} best={best.score if best else float('nan'):.4f} "
              f"front(T>={survivor_min_tier})={n_front} "
              f"converged={signals.get('converged')} evals={eval_counter[0]}")
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
