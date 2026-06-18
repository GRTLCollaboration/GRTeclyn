"""Pre-GPU behavior descriptors for shadow MAP-Elites and trajectory metrics."""

from __future__ import annotations

import math
import re
from typing import Any, Mapping

import numpy as np

from ..grtresna_convergence_gate import (
    DEFAULT_GRTRESNA_CONVERGENCE_CONFIG,
    GRTresnaConvergenceConfig,
)
from ..trajectory_log import infer_trajectory_status

_HAM_MOM_REASON_RE = re.compile(
    r"Ham=([\d.eE+-]+)%.*Mom=([\d.eE+-]+)%",
    re.DOTALL,
)

_SPEED_SUPER_C_FLOOR = 0.95
_SPEED_SUPER_C_TARGET = 1.20
_DEFAULT_SPAN = 100.0


def parse_convergence_from_reason(reason: str | None) -> dict[str, float] | None:
    if not reason:
        return None
    match = _HAM_MOM_REASON_RE.search(reason)
    if not match:
        return None
    try:
        return {
            "ham_pct": float(match.group(1)),
            "mom_pct": float(match.group(2)),
        }
    except (TypeError, ValueError):
        return None


def _convergence_from_record(record: Mapping[str, Any]) -> dict[str, float] | None:
    raw = record.get("grtresna_convergence")
    if isinstance(raw, Mapping):
        try:
            return {
                "ham_pct": float(raw["ham_pct"]),
                "mom_pct": float(raw["mom_pct"]),
            }
        except (KeyError, TypeError, ValueError):
            pass
    return parse_convergence_from_reason(str(record.get("reason") or ""))


def _axis_quality(pct: float | None, *, threshold: float, span: float = _DEFAULT_SPAN) -> float:
    if pct is None or not math.isfinite(float(pct)):
        return 0.0
    excess = max(0.0, float(pct) - threshold)
    return float(np.clip(1.0 - excess / span, 0.0, 1.0))


def _precursor_tilt(metrics: Mapping[str, Any] | None) -> float:
    if not metrics:
        return 0.0
    solved = metrics.get("solved_geometry_ftl")
    if not isinstance(solved, Mapping):
        return 0.0
    op = solved.get("operational")
    if not isinstance(op, Mapping):
        return 0.0
    c_raw = op.get("max_local_speed")
    if c_raw is None:
        return 0.0
    try:
        c = float(c_raw)
    except (TypeError, ValueError):
        return 0.0
    if not math.isfinite(c) or c <= _SPEED_SUPER_C_FLOOR:
        return 0.0
    span = _SPEED_SUPER_C_TARGET - _SPEED_SUPER_C_FLOOR
    return float(np.clip((c - _SPEED_SUPER_C_FLOOR) / span, 0.0, 1.0))


def _postload_margin(gate: Mapping[str, Any] | None) -> float:
    if not gate:
        return 0.0
    ham = gate.get("max_hamiltonian_l2")
    mom = gate.get("max_momentum_l2")
    margins: list[float] = []
    for value, limit in ((ham, 1.0e-2), (mom, 1.0e-2)):
        if value is None:
            continue
        try:
            v = float(value)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(v) or limit <= 0.0:
            continue
        margins.append(float(np.clip(1.0 - v / limit, 0.0, 1.0)))
    if not margins:
        return 0.0
    return float(min(margins))


def _score_fallback_quality(score: float | None) -> float:
    if score is None:
        return 0.5
    try:
        s = float(score)
    except (TypeError, ValueError):
        return 0.5
    if not math.isfinite(s):
        return 0.5
    return float(np.clip((s + 350.0) / 250.0, 0.0, 1.0))


def pre_gpu_descriptor_details(
    record: Mapping[str, Any],
    *,
    convergence_config: GRTresnaConvergenceConfig | None = None,
) -> dict[str, float]:
    """Status-aware pre-GPU descriptor axes in [0, 1]."""
    cfg = convergence_config or DEFAULT_GRTRESNA_CONVERGENCE_CONFIG
    status = record.get("status")
    if not isinstance(status, str):
        status = infer_trajectory_status(record)
    metrics = {
        k: record[k]
        for k in ("grtresna_convergence", "solved_geometry_ftl", "postload_gate")
        if k in record
    }
    convergence = _convergence_from_record(record)
    score = record.get("score")

    if status == "grtresna_rejected":
        ham = convergence.get("ham_pct") if convergence else None
        mom = convergence.get("mom_pct") if convergence else None
        x = _axis_quality(ham, threshold=cfg.max_ham_pct)
        y = _axis_quality(mom, threshold=cfg.max_mom_pct)
        if convergence is None:
            fallback = _score_fallback_quality(score)
            x = fallback
            y = 0.0
    elif status == "solved_ftl_rejected":
        x = 1.0
        y = _precursor_tilt(metrics)
    elif status == "postload_rejected":
        gate = metrics.get("postload_gate")
        gate_map = gate if isinstance(gate, Mapping) else None
        if convergence is not None:
            ham = convergence.get("ham_pct")
            mom = convergence.get("mom_pct")
            x = 0.5 * (
                _axis_quality(ham, threshold=cfg.max_ham_pct)
                + _axis_quality(mom, threshold=cfg.max_mom_pct)
            )
        else:
            x = _score_fallback_quality(score)
        y = _postload_margin(gate_map)
    else:
        x = _score_fallback_quality(score)
        y = 0.0

    return {
        "x": x,
        "y": y,
        "hamiltonian_quality": x if status == "grtresna_rejected" else float("nan"),
        "momentum_quality": y if status == "grtresna_rejected" else float("nan"),
        "convergence_quality": x,
        "precursor_tilt": y,
    }


def attach_pre_gpu_metrics_to_record(
    record: dict[str, Any],
    metrics: Mapping[str, Any] | None,
) -> None:
    """Copy structured pre-GPU metrics from an Evaluation onto a trajectory record."""
    if not metrics:
        return
    for key in ("grtresna_convergence", "solved_geometry_ftl", "postload_gate"):
        if key in metrics:
            record[key] = metrics[key]
