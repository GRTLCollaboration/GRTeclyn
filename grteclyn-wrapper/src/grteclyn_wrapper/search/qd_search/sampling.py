"""Candidate sampling and mutation for MAP-Elites."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from ..optimize import SearchDimension
from ..pre_gpu.near_miss_pool import NearMissPool
from .archive import Elite, QDArchive


@dataclass(frozen=True)
class QdSamplingConfig:
    elite_fraction: float = 0.60
    shadow_fraction: float = 0.15
    near_miss_fraction: float = 0.15
    feasible_fraction: float = 0.05
    random_fraction: float = 0.05
    mutation_sigma: float = 0.15


def _reflect(value: float, lower: float, upper: float) -> float:
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


def _feasible_bounds_from_params(
    dims: Sequence[SearchDimension],
    param_sets: Sequence[Mapping[str, float]],
) -> list[tuple[float, float]] | None:
    if not param_sets:
        return None
    bounds: list[tuple[float, float]] = []
    for d in dims:
        vals = [
            float(params[d.param_key])
            for params in param_sets
            if d.param_key in params
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
    return [float(rng.uniform(lo, hi)) for (lo, hi) in bounds]


def _mutate_params(
    params: Mapping[str, float],
    dims: Sequence[SearchDimension],
    rng: np.random.Generator,
    *,
    sigma: float = 0.15,
) -> list[float]:
    x = []
    for d in dims:
        base = float(params.get(d.param_key, d.center))
        step = rng.normal(0.0, sigma * max(d.range, 1.0e-9))
        x.append(_reflect(base + step, d.lower, d.upper))
    return x


def _mutate_elite(
    elite: Elite,
    dims: Sequence[SearchDimension],
    rng: np.random.Generator,
    *,
    sigma: float = 0.15,
) -> list[float]:
    return _mutate_params(elite.params, dims, rng, sigma=sigma)


def sample_next_candidate(
    *,
    dims: Sequence[SearchDimension],
    rng: np.random.Generator,
    archive: QDArchive,
    pre_gpu_archive: QDArchive | None = None,
    near_miss_pool: NearMissPool | None = None,
    config: QdSamplingConfig | None = None,
) -> list[float]:
    """Sample the next MAP-Elites candidate from elites, shadow, near-miss, or random."""
    cfg = config or QdSamplingConfig()
    sources: list[tuple[str, float]] = []
    if archive.cells:
        sources.append(("elite", cfg.elite_fraction))
    if pre_gpu_archive is not None and pre_gpu_archive.cells:
        sources.append(("shadow", cfg.shadow_fraction))
    if near_miss_pool is not None and len(near_miss_pool) > 0:
        sources.append(("near_miss", cfg.near_miss_fraction))

    param_sets: list[Mapping[str, float]] = []
    if near_miss_pool is not None:
        param_sets.extend(near_miss_pool.param_sets())
    if pre_gpu_archive is not None:
        param_sets.extend(e.params for e in pre_gpu_archive.cells.values())
    feasible_bounds = _feasible_bounds_from_params(dims, param_sets)
    if feasible_bounds is None and archive.cells:
        feasible_bounds = _feasible_bounds(dims, list(archive.cells.values()))
    if feasible_bounds is not None:
        sources.append(("feasible", cfg.feasible_fraction))
    sources.append(("random", cfg.random_fraction))

    total = sum(weight for _name, weight in sources)
    roll = float(rng.random()) * total
    cumulative = 0.0
    choice = "random"
    for name, weight in sources:
        cumulative += weight
        if roll <= cumulative:
            choice = name
            break

    if choice == "elite":
        elites = list(archive.cells.values())
        parent = elites[int(rng.integers(len(elites)))]
        return _mutate_elite(parent, dims, rng, sigma=cfg.mutation_sigma)
    if choice == "shadow" and pre_gpu_archive is not None and pre_gpu_archive.cells:
        shadows = list(pre_gpu_archive.cells.values())
        parent = shadows[int(rng.integers(len(shadows)))]
        return _mutate_elite(parent, dims, rng, sigma=cfg.mutation_sigma)
    if choice == "near_miss" and near_miss_pool is not None:
        params = near_miss_pool.sample_params(rng)
        if params is not None:
            return _mutate_params(params, dims, rng, sigma=cfg.mutation_sigma)
    if choice == "feasible" and feasible_bounds is not None:
        return _sample_feasible_box(dims, feasible_bounds, rng)
    return _sample_random(dims, rng)
