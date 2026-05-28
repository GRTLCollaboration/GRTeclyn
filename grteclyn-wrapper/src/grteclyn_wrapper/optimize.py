"""CMA-ES optimization driver for the RadialRecipe metric search.

Replaces random sampling with gradient-free optimization over the
Gaussian basis coefficients.  Each evaluation runs a full GRTeclyn
episode and returns the negative score (CMA-ES minimizes).

Supports multi-GPU parallel evaluation within each generation.
"""

from __future__ import annotations

import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from .config import (
    DEFAULT_RADIAL_RECIPE_TEMPLATE,
    ExecutableConfig,
    ExampleConfig,
    resolve_example,
)
from .constrained_recipe import constrained_overrides
from .episode import Episode, create_episode, update_metadata, write_json
from .metrics import dataclass_to_dict, read_episode_metrics
from .params import write_params
from .preflight import preflight_check
from .runner import run_episode
from .score import Score, score_episode

try:
    import cma
except ImportError:
    cma = None  # type: ignore[assignment]

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]


@dataclass
class SearchDimension:
    """One axis in the optimizer's search space."""

    param_key: str
    lower: float
    upper: float
    initial: float | None = None

    @property
    def center(self) -> float:
        if self.initial is not None:
            return self.initial
        return 0.5 * (self.lower + self.upper)

    @property
    def range(self) -> float:
        return self.upper - self.lower


DEFAULT_SEARCH_SPACE: list[SearchDimension] = [
    SearchDimension("recipe_chi_coeff_0", -0.5, 0.1, -0.1),
    SearchDimension("recipe_chi_coeff_1", -0.3, 0.3, 0.0),
    SearchDimension("recipe_chi_coeff_2", -0.2, 0.2, 0.0),
    SearchDimension("recipe_chi_coeff_3", -0.2, 0.2, 0.0),
    SearchDimension("recipe_basis_width", 0.3, 3.0, 1.0),
]


@dataclass(frozen=True)
class OptimizeResult:
    best_params: dict[str, float]
    best_score: float
    best_episode: str
    generations: int
    evaluations: int
    trajectory: list[dict[str, Any]]


def _vector_to_overrides(
    x: Sequence[float],
    dims: Sequence[SearchDimension],
    base: Mapping[str, Any],
) -> dict[str, Any]:
    overrides = dict(base)
    for xi, dim in zip(x, dims):
        clamped = max(dim.lower, min(dim.upper, xi))
        overrides[dim.param_key] = clamped
    return overrides


def _assign_gpu(index: int, gpu_ids: Sequence[int]) -> str:
    """Round-robin GPU assignment for parallel evaluation."""
    return str(gpu_ids[index % len(gpu_ids)])


def _objective(
    x: Sequence[float],
    *,
    dims: Sequence[SearchDimension],
    base_overrides: Mapping[str, Any],
    opt_dir: Path,
    example: ExampleConfig,
    template: Path,
    executable: ExecutableConfig | None,
    eval_counter: list[int],
    constrained: bool,
    phantom: bool,
    use_preflight: bool,
    cuda_devices: str | None,
    check_params: bool,
    dry_run: bool,
    trajectory: list[dict[str, Any]],
    target_stop_time: float | None,
    score_weights: Mapping[str, float] | None,
    consume_plotfiles: bool = True,
    consumer_radii: Sequence[float] = (4.0, 8.0),
) -> float:
    """Evaluate one candidate.  Returns negative score (CMA-ES minimizes)."""
    eval_counter[0] += 1
    idx = eval_counter[0]

    overrides = _vector_to_overrides(x, dims, base_overrides)

    if constrained:
        constrained_overrides(overrides, phantom=phantom)

    if use_preflight:
        pf = preflight_check(overrides, phantom=phantom)
        if not pf.passed:
            record = {
                "eval": idx,
                "score": 0.0,
                "preflight_rejected": True,
                "reason": pf.reason,
                "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
            }
            trajectory.append(record)
            return 100.0

    episode = create_episode(
        opt_dir,
        name=f"eval_{idx:06d}",
        metadata={
            "mode": "optimize",
            "example": example.name,
            "eval_index": idx,
            "overrides": overrides,
        },
    )
    write_params(
        template, episode.params_path,
        episode_dir=episode.path, example=example, overrides=overrides,
    )

    exit_code: int | None = None
    if dry_run:
        update_metadata(episode, {"dry_run": True})
    else:
        if executable is None:
            raise ValueError("executable is required unless dry_run=True")
        try:
            result = run_episode(
                episode, executable,
                check_params=check_params,
                cuda_devices=cuda_devices,
                consume_plotfiles=consume_plotfiles,
                consumer_radii=consumer_radii,
                consumer_delete=True,
            )
            exit_code = result.returncode
        except Exception as exc:
            exit_code = 1
            update_metadata(episode, {
                "simulation_error": repr(exc),
                "simulation_exit_code": exit_code,
            })

    metrics = read_episode_metrics(episode.path)
    score = score_episode(
        metrics, target_stop_time=target_stop_time, weights=score_weights,
    )

    write_json(episode.score_path, {
        "score": dataclass_to_dict(score),
        "metrics": dataclass_to_dict(metrics),
    })

    record = {
        "eval": idx,
        "episode": str(episode.path),
        "exit_code": exit_code,
        "score": score.total,
        "components": score.components,
        "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
    }
    trajectory.append(record)
    return -score.total


def run_optimize(
    *,
    runs_dir: Path,
    executable: ExecutableConfig | None = None,
    max_generations: int = 50,
    population_size: int | None = None,
    sigma0: float = 0.3,
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
    x0: Sequence[float] | None = None,
    consume_plotfiles: bool = True,
    consumer_radii: Sequence[float] = (4.0, 8.0),
) -> OptimizeResult:
    """Run multi-GPU CMA-ES optimization loop.

    Parameters
    ----------
    gpu_ids : list of int, optional
        Available GPU indices for parallel evaluation. If None, uses
        cuda_devices for sequential mode. If provided, each member of
        the CMA-ES population is assigned a GPU round-robin.
    """
    if cma is None:
        raise ImportError(
            "CMA-ES optimization requires the 'cma' package. "
            "Install it with: pip install cma"
        )
    if np is None:
        raise ImportError("numpy is required for optimization.")

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
    base.setdefault("regrid_threshold", 0.01)
    max_lvl = int(base["max_level"])
    if "regrid_interval" not in base and max_lvl > 0:
        intervals = [16] * min(max_lvl, 2) + [8] * max(0, max_lvl - 2)
        base["regrid_interval"] = intervals

    target_stop_time = float(base["stop_time"])

    if name is None:
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        name = f"optimize_{timestamp}"
    opt_dir = (runs_dir / name).expanduser().resolve()
    opt_dir.mkdir(parents=True, exist_ok=False)

    if x0 is not None:
        initial = list(x0)
    else:
        initial = [d.center for d in dims]

    bounds = [[d.lower for d in dims], [d.upper for d in dims]]

    effective_popsize = population_size
    if effective_popsize is None and gpu_ids is not None:
        effective_popsize = len(gpu_ids)

    opts: dict[str, Any] = {
        "maxiter": max_generations,
        "bounds": bounds,
        "CMA_stds": [sigma0 * d.range for d in dims],
        "verbose": -1,
    }
    if effective_popsize is not None:
        opts["popsize"] = effective_popsize
    if seed is not None:
        opts["seed"] = seed

    eval_counter = [0]
    trajectory: list[dict[str, Any]] = []

    write_json(opt_dir / "metadata.json", {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "example": example_cfg.name,
        "max_generations": max_generations,
        "population_size": effective_popsize,
        "sigma0": sigma0,
        "seed": seed,
        "constrained": constrained,
        "phantom": phantom,
        "use_preflight": use_preflight,
        "dry_run": dry_run,
        "gpu_ids": list(gpu_ids) if gpu_ids else None,
        "consume_plotfiles": consume_plotfiles,
        "base_overrides": base,
        "search_space": [
            {"key": d.param_key, "lower": d.lower, "upper": d.upper, "initial": d.center}
            for d in dims
        ],
    })

    es = cma.CMAEvolutionStrategy(initial, sigma0, opts)

    best_score = -math.inf
    best_params: dict[str, float] = {}
    best_episode = ""
    gen = 0

    print(f"[optimize] Starting CMA-ES: {len(dims)}D, popsize={es.popsize}, "
          f"max_gen={max_generations}, GPUs={gpu_ids or cuda_devices}")

    while not es.stop():
        gen += 1
        solutions = es.ask()
        fitnesses = []

        if gpu_ids is not None and len(gpu_ids) > 1 and not dry_run:
            fitnesses = _evaluate_generation_parallel(
                solutions,
                dims=dims,
                base_overrides=base,
                opt_dir=opt_dir,
                example=example_cfg,
                template=tpl,
                executable=executable,
                eval_counter=eval_counter,
                constrained=constrained,
                phantom=phantom,
                use_preflight=use_preflight,
                gpu_ids=gpu_ids,
                check_params=check_params,
                trajectory=trajectory,
                target_stop_time=target_stop_time,
                score_weights=score_weights,
                consume_plotfiles=consume_plotfiles,
                consumer_radii=consumer_radii,
            )
        else:
            for sol in solutions:
                f = _objective(
                    sol,
                    dims=dims,
                    base_overrides=base,
                    opt_dir=opt_dir,
                    example=example_cfg,
                    template=tpl,
                    executable=executable,
                    eval_counter=eval_counter,
                    constrained=constrained,
                    phantom=phantom,
                    use_preflight=use_preflight,
                    cuda_devices=cuda_devices,
                    check_params=check_params,
                    dry_run=dry_run,
                    trajectory=trajectory,
                    target_stop_time=target_stop_time,
                    score_weights=score_weights,
                    consume_plotfiles=consume_plotfiles,
                    consumer_radii=consumer_radii,
                )
                fitnesses.append(f)

        es.tell(solutions, fitnesses)

        gen_scores = [rec.get("score", -math.inf) for rec in trajectory[-len(solutions):]]
        gen_best = max(gen_scores) if gen_scores else -math.inf
        gen_mean = sum(gen_scores) / len(gen_scores) if gen_scores else 0.0

        for rec in trajectory[-len(solutions):]:
            sc = rec.get("score", -math.inf)
            if sc > best_score:
                best_score = sc
                best_params = dict(rec.get("overrides", {}))
                best_episode = rec.get("episode", "")

        print(f"[optimize] gen {gen}/{max_generations}: "
              f"best={gen_best:.4f} mean={gen_mean:.4f} "
              f"all-time-best={best_score:.4f} evals={eval_counter[0]}")

        with (opt_dir / "trajectory.jsonl").open("w", encoding="utf-8") as fh:
            for rec in trajectory:
                fh.write(json.dumps(rec, sort_keys=True) + "\n")

    result = OptimizeResult(
        best_params=best_params,
        best_score=best_score,
        best_episode=best_episode,
        generations=gen,
        evaluations=eval_counter[0],
        trajectory=trajectory,
    )

    write_json(opt_dir / "result.json", {
        "best_params": result.best_params,
        "best_score": result.best_score,
        "best_episode": result.best_episode,
        "generations": result.generations,
        "evaluations": result.evaluations,
    })

    print(f"\n[optimize] Done. {gen} generations, {eval_counter[0]} evaluations.")
    print(f"[optimize] Best score: {best_score:.4f}")
    print(f"[optimize] Best episode: {best_episode}")
    print(f"[optimize] Results: {opt_dir}")

    return result


def _evaluate_generation_parallel(
    solutions: list,
    *,
    dims: Sequence[SearchDimension],
    base_overrides: Mapping[str, Any],
    opt_dir: Path,
    example: ExampleConfig,
    template: Path,
    executable: ExecutableConfig | None,
    eval_counter: list[int],
    constrained: bool,
    phantom: bool,
    use_preflight: bool,
    gpu_ids: Sequence[int],
    check_params: bool,
    trajectory: list[dict[str, Any]],
    target_stop_time: float | None,
    score_weights: Mapping[str, float] | None,
    consume_plotfiles: bool,
    consumer_radii: Sequence[float],
) -> list[float]:
    """Evaluate an entire CMA-ES generation in parallel across GPUs.

    Each solution is assigned a GPU round-robin, and all are launched
    concurrently via threads (GIL released during subprocess waits).
    """
    import threading

    fitnesses: list[float | None] = [None] * len(solutions)
    lock = threading.Lock()

    def _eval_one(idx_in_gen: int, sol) -> None:
        gpu = _assign_gpu(idx_in_gen, gpu_ids)
        f = _objective(
            sol,
            dims=dims,
            base_overrides=base_overrides,
            opt_dir=opt_dir,
            example=example,
            template=template,
            executable=executable,
            eval_counter=eval_counter,
            constrained=constrained,
            phantom=phantom,
            use_preflight=use_preflight,
            cuda_devices=gpu,
            check_params=check_params,
            dry_run=False,
            trajectory=trajectory,
            target_stop_time=target_stop_time,
            score_weights=score_weights,
            consume_plotfiles=consume_plotfiles,
            consumer_radii=consumer_radii,
        )
        with lock:
            fitnesses[idx_in_gen] = f

    threads = []
    for i, sol in enumerate(solutions):
        t = threading.Thread(target=_eval_one, args=(i, sol))
        t.start()
        threads.append(t)

    for t in threads:
        t.join()

    return [f if f is not None else 100.0 for f in fitnesses]
