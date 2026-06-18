"""Candidate vector utilities: warm-start, jitter, exotic injection."""

from __future__ import annotations

import json
import math
import random
from pathlib import Path
from typing import Any, Mapping, Sequence

from .config import _LUMP_KEY_RE
from .dimension import SearchDimension

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


def _clip_vector(x: Sequence[float], dims: Sequence[SearchDimension]) -> list[float]:
    """Clamp an optimizer vector to the configured search bounds."""
    return [
        max(dim.lower, min(dim.upper, float(xi)))
        for xi, dim in zip(x, dims)
    ]


def _overrides_to_vector(
    overrides: Mapping[str, Any],
    dims: Sequence[SearchDimension],
) -> list[float]:
    """Map a trajectory overrides dict back onto the current search vector."""
    vals: list[float] = []
    for dim in dims:
        raw = overrides.get(dim.param_key, dim.center)
        try:
            value = float(raw)
        except (TypeError, ValueError):
            value = dim.center
        vals.append(max(dim.lower, min(dim.upper, value)))
    return vals


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


def _load_warm_start_vectors(
    trajectories: Sequence[Path],
    dims: Sequence[SearchDimension],
    top_k: int,
    *,
    include_near_miss: bool = False,
    near_miss_k: int = 4,
) -> list[list[float]]:
    """Load top-scoring vectors from one or more prior trajectory.jsonl files."""
    all_records: list[dict[str, Any]] = []
    gpu_records: list[tuple[float, Mapping[str, Any]]] = []
    for path in trajectories:
        path = path.expanduser().resolve()
        if not path.exists():
            continue
        for rec in _load_trajectory_records(path):
            all_records.append(rec)
            if rec.get("surrogate_predicted"):
                continue
            status = rec.get("status")
            if isinstance(status, str) and status != "gpu_ok":
                continue
            score = rec.get("score")
            overrides = rec.get("overrides")
            if not isinstance(overrides, Mapping):
                continue
            try:
                score_f = float(score)
            except (TypeError, ValueError):
                continue
            if math.isfinite(score_f):
                gpu_records.append((score_f, overrides))
    gpu_records.sort(key=lambda item: item[0], reverse=True)
    vectors = [
        _overrides_to_vector(overrides, dims)
        for _, overrides in gpu_records[:top_k]
    ]
    if include_near_miss and near_miss_k > 0 and all_records:
        from ..pre_gpu.near_miss_pool import near_miss_vectors_from_trajectory

        near_vectors = near_miss_vectors_from_trajectory(
            all_records, dims, max_size=near_miss_k,
        )
        seen = {tuple(round(x, 8) for x in vec) for vec in vectors}
        for vec in near_vectors:
            key = tuple(round(x, 8) for x in vec)
            if key in seen:
                continue
            vectors.append(vec)
            seen.add(key)
    return vectors


def _random_vector(
    dims: Sequence[SearchDimension],
    rng: random.Random,
) -> list[float]:
    return [rng.uniform(dim.lower, dim.upper) for dim in dims]


def _jitter_vector(
    x: Sequence[float],
    dims: Sequence[SearchDimension],
    rng: random.Random,
    scale: float,
) -> list[float]:
    jittered = [
        float(xi) + rng.gauss(0.0, scale * dim.range)
        for xi, dim in zip(x, dims)
    ]
    return _clip_vector(jittered, dims)


def _force_exotic_template(
    x: Sequence[float],
    dims: Sequence[SearchDimension],
    rng: random.Random,
    template_index: int,
) -> list[float]:
    """Impose a categorical exotic/multi-lump pattern on a continuous vector."""
    vec = _clip_vector(x, dims)
    index = {dim.param_key: i for i, dim in enumerate(dims)}
    exotic_frac_key = next(
        (k for k in ("grtresna_ring_exotic_fraction", "grtresna_shell_exotic_fraction") if k in index),
        None,
    )
    if exotic_frac_key is not None:
        frac_i = index[exotic_frac_key]
        frac_dim = dims[frac_i]
        phase_key = exotic_frac_key.replace("_fraction", "_phase")
        phase_i = index.get(phase_key)
        patterns = (0.25, 0.40, 0.60, 1.0)
        vec[frac_i] = max(
            frac_dim.lower,
            min(frac_dim.upper, patterns[template_index % len(patterns)]),
        )
        if phase_i is not None:
            phase_dim = dims[phase_i]
            phase = (template_index % 8) * (2.0 * math.pi / 8.0)
            vec[phase_i] = max(phase_dim.lower, min(phase_dim.upper, phase))
        return vec

    lump_ids = sorted({
        int(m.group(1))
        for dim in dims
        if (m := _LUMP_KEY_RE.match(dim.param_key))
    })
    if not lump_ids:
        return vec

    n = len(lump_ids)
    if template_index % 4 == 0:
        exotic = {lump_ids[n // 2]}
    elif template_index % 4 == 1:
        exotic = {lump_ids[0], lump_ids[-1]}
    elif template_index % 4 == 2:
        exotic = {k for i, k in enumerate(lump_ids) if i % 2 == 0}
    else:
        exotic = set(lump_ids)

    radius = rng.uniform(2.0, 5.0)
    for pos, k in enumerate(lump_ids):
        angle = 2.0 * math.pi * pos / max(n, 1)
        tangential = 0.45 if k in exotic else 0.30
        assignments = {
            "amp": rng.uniform(0.12, 0.24),
            "width": rng.uniform(1.5, 3.5),
            "center_x": radius * math.cos(angle),
            "center_y": radius * math.sin(angle),
            "center_z": rng.uniform(-2.0, 2.0),
            "velocity_x": -tangential * math.sin(angle),
            "velocity_y": tangential * math.cos(angle),
            "velocity_z": rng.uniform(-0.15, 0.15),
            "omega": rng.uniform(-0.15, 0.15),
            "mode": 1.0 if pos % 2 == 0 else 2.0,
            "exotic": 1.0 if k in exotic else 0.0,
        }
        for suffix, value in assignments.items():
            key = f"grtresna_lump{k}_{suffix}"
            if key in index:
                i = index[key]
                dim = dims[i]
                vec[i] = max(dim.lower, min(dim.upper, value))
    return vec
