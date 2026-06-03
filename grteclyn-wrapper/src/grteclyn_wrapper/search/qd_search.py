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
from .validation_tiers import (
    DEFAULT_TIER_CONFIG,
    Tier,
    TierConfig,
    build_survivors,
    convergence_signals,
    evaluate_tiers,
    survivor_front,
)


def _descriptors(components: Mapping[str, float]) -> tuple[float, float]:
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
    return ftl_benefit, mechanism


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
    tier_config: TierConfig = DEFAULT_TIER_CONFIG,
    survivor_min_tier: int = int(Tier.OPERATIONAL),
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
        "base_overrides": base,
        "search_space": [
            {"key": d.param_key, "lower": d.lower, "upper": d.upper} for d in dims
        ],
    })

    def _evaluate_batch(vectors: list[list[float]]) -> list[Evaluation | None]:
        results: list[Evaluation | None] = [None] * len(vectors)

        def _eval_one(i: int, x: list[float], gpu: str | None) -> None:
            eval_counter[0] += 1
            idx = eval_counter[0]
            overrides = _vector_to_overrides(x, dims, base)
            results[i] = evaluate_overrides(
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
                ftl_L=ftl_L,
                consume_plotfiles=consume_plotfiles,
                consumer_radii=consumer_radii,
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

    def _ingest(vectors: list[list[float]], results: list[Evaluation | None]) -> None:
        for x, res in zip(vectors, results):
            if res is None or res.preflight_rejected:
                trajectory.append({"preflight_rejected": True})
                continue
            d1, d2 = _descriptors(res.components)
            cell = (_bin_index(d1, bins), _bin_index(d2, bins))
            params = {d.param_key: float(_vector_to_overrides(x, dims, base)[d.param_key]) for d in dims}
            assessment = evaluate_tiers(res.components, metrics=res.metrics, config=tier_config)
            elite = Elite(
                cell=cell,
                score=res.score,
                descriptors=(d1, d2),
                params=params,
                episode=res.episode_path,
                tier=assessment.tier,
                tier_name=assessment.tier_name,
                objectives=assessment.objectives,
            )
            improved = archive.insert(elite)
            trajectory.append({
                "eval": eval_counter[0],
                "score": res.score,
                "cell": list(cell),
                "descriptors": [d1, d2],
                "tier": assessment.tier,
                "tier_name": assessment.tier_name,
                "improved": improved,
                "episode": res.episode_path,
            })

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
    _ingest(init_vectors, _evaluate_batch(init_vectors))

    for it in range(1, iterations + 1):
        vectors: list[list[float]] = []
        elites = list(archive.cells.values())
        for _ in range(batch):
            if elites and rng.random() < 0.7:
                parent = elites[int(rng.integers(len(elites)))]
                vectors.append(_mutate_elite(parent, dims, rng))
            else:
                vectors.append(_sample_random(dims, rng))
        _ingest(vectors, _evaluate_batch(vectors))

        write_json(qd_dir / "archive.json", archive.to_dict())
        with (qd_dir / "trajectory.jsonl").open("w", encoding="utf-8") as fh:
            for rec in trajectory:
                fh.write(json.dumps(rec, sort_keys=True) + "\n")
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
