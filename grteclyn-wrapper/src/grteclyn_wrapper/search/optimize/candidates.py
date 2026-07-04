"""Candidate vector utilities: warm-start, jitter, exotic injection."""

from __future__ import annotations

import json
import math
import random
import re
from pathlib import Path
from typing import Any, Mapping, Sequence

from .config import _LUMP_KEY_RE
from .dimension import SearchDimension

# Hard cap on per-lump tangential speed v_t = R0 * |omega_rot| (geometric units,
# c = 1).  A trajectory lump is dragged along its orbit by a co-moving trap whose
# TARGET centre advances at v_t; if v_t is superluminal (or even mildly
# relativistic on a tight curve) no soliton can follow it adiabatically and the
# lump radiates itself apart (the "spotlight mismatch").  The historical search
# left R0 in [1.5, 8] and omega_rot in [-1, 1], so v_t reached ~8c -- physically
# unconfineable.  We clamp omega_rot per lump so v_t <= TRAJECTORY_V_MAX_DEFAULT
# unless a campaign overrides ``trajectory_v_max``.  Set ``trajectory_v_max <= 0``
# to disable the cap (NOT recommended).
TRAJECTORY_V_MAX_DEFAULT = 0.3
# Keep r(t) >= TRAJECTORY_R_MIN through stop_time (matches C++ TRAJECTORY_R_MIN).
TRAJECTORY_R_MIN_DEFAULT = 0.1
TRAJECTORY_STOP_TIME_DEFAULT = 16.0

_TRAJECTORY_OMEGA_RE = re.compile(r"^trajectory_lump(\d+)_omega_rot$")
_TRAJECTORY_V_RAD_RE = re.compile(r"^trajectory_lump(\d+)_v_rad$")


def _min_spiral_v_rad(
    r0: float,
    *,
    a_breath: float,
    stop_time: float,
    r_min: float,
) -> float:
    """Most-inward v_rad that keeps r(stop_time) >= r_min (worst-case breathing)."""
    if stop_time <= 0.0 or r0 <= 0.0:
        return 0.0
    return (r_min - r0 + abs(a_breath)) / stop_time


def _clamp_trajectory_spiral_radius(
    overrides: dict[str, Any],
    *,
    default_r_min: float = TRAJECTORY_R_MIN_DEFAULT,
    default_stop_time: float = TRAJECTORY_STOP_TIME_DEFAULT,
) -> dict[str, Any]:
    """Clamp inward ``v_rad`` so the spiral stays outside ``r_min`` until stop_time.

    Uses the same floor as C++ ``TrajectoryEvaluator`` (``TRAJECTORY_R_MIN``).
    Outward drift (``v_rad > 0``) is unchanged.  Mutates and returns *overrides*.
    """
    try:
        r_min = float(overrides.get("trajectory_r_min", default_r_min))
    except (TypeError, ValueError):
        r_min = default_r_min
    try:
        stop_time = float(overrides.get("stop_time", default_stop_time))
    except (TypeError, ValueError):
        stop_time = default_stop_time
    if not math.isfinite(r_min) or r_min <= 0.0 or stop_time <= 0.0:
        return overrides

    try:
        a_breath = float(overrides.get("trajectory_A_breath", 0.0))
    except (TypeError, ValueError):
        a_breath = 0.0
    try:
        v_rad_inward_floor = float(overrides.get("trajectory_v_rad_inward_floor", float("-inf")))
    except (TypeError, ValueError):
        v_rad_inward_floor = float("-inf")

    lump_indices: set[str] = set()
    for key in overrides:
        match = _TRAJECTORY_OMEGA_RE.match(str(key))
        if match is not None:
            lump_indices.add(match.group(1))
            continue
        match = _TRAJECTORY_V_RAD_RE.match(str(key))
        if match is not None:
            lump_indices.add(match.group(1))

    for lump_idx in lump_indices:
        v_rad_key = f"trajectory_lump{lump_idx}_v_rad"
        try:
            v_rad = float(overrides.get(v_rad_key, 0.0))
        except (TypeError, ValueError):
            v_rad = 0.0
        try:
            r0 = float(overrides.get(f"trajectory_lump{lump_idx}_R0", 0.0))
        except (TypeError, ValueError):
            r0 = 0.0
        if r0 <= 0.0:
            continue
        v_rad_floor = _min_spiral_v_rad(
            r0, a_breath=a_breath, stop_time=stop_time, r_min=r_min
        )
        if math.isfinite(v_rad_inward_floor):
            v_rad_floor = max(v_rad_floor, v_rad_inward_floor)
        if v_rad < v_rad_floor:
            overrides[v_rad_key] = v_rad_floor
    return overrides


def _clamp_trajectory_speed(
    overrides: dict[str, Any],
    *,
    default_v_max: float = TRAJECTORY_V_MAX_DEFAULT,
) -> dict[str, Any]:
    """Convert normalized trajectory speed fractions to physical omega_rot / v_rad.

    The search space now treats ``trajectory_lump*_omega_rot`` and
    ``trajectory_lump*_v_rad`` as *normalized fractions* of ``trajectory_v_max``:

      * omega_rot_norm = -1..1  ->  tangential speed v_t = omega_rot_norm * v_max
      * v_rad_norm     = -1..1  ->  radial speed      v_r = v_rad_norm     * v_max

    The combined normalized speed is clamped to the unit disk, then converted to
    physical quantities:

      * omega_rot = v_t / R0
      * v_rad     = v_r

    This turns the speed budget into a smooth 2-D disk for the optimizer, instead
    of a hard-cornered box where almost every non-zero omega_rot hits the speed
    boundary and every configuration ends up with exactly the same max speed.

    Mutates and returns *overrides*.  No-op when no ``trajectory_lump*_omega_rot``
    keys are present, so non-trajectory campaigns are unaffected.
    """
    try:
        v_max = float(overrides.get("trajectory_v_max", default_v_max))
    except (TypeError, ValueError):
        v_max = default_v_max
    if not math.isfinite(v_max) or v_max <= 0.0:
        return overrides

    for key in list(overrides.keys()):
        match = _TRAJECTORY_OMEGA_RE.match(str(key))
        if match is None:
            continue
        lump_idx = match.group(1)
        try:
            omega_norm = float(overrides[key])
        except (TypeError, ValueError):
            continue
        try:
            r0 = float(overrides.get(f"trajectory_lump{lump_idx}_R0", 0.0))
        except (TypeError, ValueError):
            r0 = 0.0
        if r0 <= 0.0:
            continue
        v_rad_key = f"trajectory_lump{lump_idx}_v_rad"
        try:
            v_rad_norm = float(overrides.get(v_rad_key, 0.0))
        except (TypeError, ValueError):
            v_rad_norm = 0.0

        # Clamp the combined normalized speed to the unit disk.
        norm_total = math.sqrt(omega_norm * omega_norm + v_rad_norm * v_rad_norm)
        if norm_total > 1.0:
            scale = 1.0 / norm_total
            omega_norm *= scale
            v_rad_norm *= scale

        # Convert normalized fractions to physical angular velocity / radial drift.
        overrides[key] = omega_norm * v_max / r0
        if v_rad_key in overrides:
            overrides[v_rad_key] = v_rad_norm * v_max
    return overrides


def _enforce_retrograde(
    overrides: dict[str, Any],
) -> dict[str, Any]:
    """Force all trajectory lumps retrograde (``omega_rot <= 0``).

    Mutates and returns *overrides*.  Controlled by
    ``trajectory_retrograde_only`` (default ``False``).  When enabled, any
    positive ``omega_rot`` is negated; exact zero is left as-is (a
    non-rotating lump is acceptable).  This removes 50% of the search space
    — HQ validation showed counter-rotation is a false-positive generator
    and all confirmed FTL configs are all-retrograde.
    """
    try:
        enabled = bool(int(overrides.get("trajectory_retrograde_only", 0)))
    except (TypeError, ValueError):
        enabled = False
    if not enabled:
        return overrides

    for key in list(overrides.keys()):
        match = _TRAJECTORY_OMEGA_RE.match(str(key))
        if match is None:
            continue
        try:
            omega_rot = float(overrides[key])
        except (TypeError, ValueError):
            continue
        if omega_rot > 0.0:
            overrides[key] = -omega_rot
    return overrides


def _vector_to_overrides(
    x: Sequence[float],
    dims: Sequence[SearchDimension],
    base: Mapping[str, Any],
) -> dict[str, Any]:
    overrides = dict(base)
    for xi, dim in zip(x, dims):
        clamped = max(dim.lower, min(dim.upper, xi))
        overrides[dim.param_key] = clamped
    # Convert normalized speed fractions to physical omega_rot / v_rad first.
    _clamp_trajectory_speed(overrides)
    # Keep spirals inside a physical radius band (operates on physical v_rad).
    _clamp_trajectory_spiral_radius(overrides)
    # Optionally force all-retrograde orbits (omega_rot <= 0).
    _enforce_retrograde(overrides)
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
