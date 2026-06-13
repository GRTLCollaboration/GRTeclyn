"""Candidate sampling and mutation for MAP-Elites."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from ..optimize import SearchDimension
from .archive import Elite

_ELITE_FRACTION = 0.85


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


def _sample_random(dims: Sequence[SearchDimension], rng: np.random.Generator) -> list[float]:
    return [float(rng.uniform(d.lower, d.upper)) for d in dims]


def _sample_feasible_box(
    dims: Sequence[SearchDimension],
    bounds: Sequence[tuple[float, float]],
    rng: np.random.Generator,
) -> list[float]:
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
