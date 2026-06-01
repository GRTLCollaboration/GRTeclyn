"""MAP-Elites quality-diversity search for the Spacetime Failure Atlas.

The stated deliverable is an atlas of how geometries behave/fail, which is a
quality-diversity problem rather than single-optimum search.  MAP-Elites keeps
the best candidate in each cell of a behavior-descriptor grid, illuminating the
whole constructibility map in one campaign while resisting premature
convergence.

Behavior descriptors (both in [0, 1]):
  * ``ftl_benefit``  -- the log-amplified FTL shortcut (``ftl_shortcut``).
  * ``exoticity``    -- ``1 - anec_condition`` (how much exotic energy the
    geometry needs along the travel axis).

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


def _descriptors(components: Mapping[str, float]) -> tuple[float, float]:
    ftl_benefit = float(np.clip(components.get("ftl_shortcut", 0.0), 0.0, 1.0))
    exoticity = float(np.clip(1.0 - components.get("anec_condition", 1.0), 0.0, 1.0))
    return ftl_benefit, exoticity


def _bin_index(value: float, bins: int) -> int:
    return int(min(bins - 1, max(0, math.floor(value * bins))))


@dataclass
class Elite:
    cell: tuple[int, int]
    score: float
    descriptors: tuple[float, float]
    params: dict[str, float]
    episode: str | None


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

    def to_dict(self) -> dict[str, Any]:
        return {
            "bins": self.bins,
            "coverage": self.coverage,
            "num_elites": len(self.cells),
            "best_score": self.best.score if self.best else None,
            "cells": [
                {
                    "cell": list(e.cell),
                    "score": e.score,
                    "descriptors": list(e.descriptors),
                    "params": e.params,
                    "episode": e.episode,
                }
                for e in sorted(self.cells.values(), key=lambda e: -e.score)
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
            elite = Elite(cell=cell, score=res.score, descriptors=(d1, d2), params=params, episode=res.episode_path)
            improved = archive.insert(elite)
            trajectory.append({
                "eval": eval_counter[0],
                "score": res.score,
                "cell": list(cell),
                "descriptors": [d1, d2],
                "improved": improved,
                "episode": res.episode_path,
            })

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
        best = archive.best
        print(f"[qd] iter {it}/{iterations}: elites={len(archive.cells)} "
              f"coverage={archive.coverage:.2f} best={best.score if best else float('nan'):.4f} "
              f"evals={eval_counter[0]}")

    write_json(qd_dir / "archive.json", archive.to_dict())
    print(f"[qd] Done. elites={len(archive.cells)} coverage={archive.coverage:.2f} dir={qd_dir}")
    return archive
