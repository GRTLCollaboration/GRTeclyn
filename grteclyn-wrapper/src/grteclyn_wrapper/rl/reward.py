from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class FenceResult:
    terminated: bool = False
    penalty: float = 0.0
    reason: str = ""


@dataclass
class TaxManState:
    dense_history: list[float] = field(default_factory=list)
    ftl_window: list[float] = field(default_factory=list)
    prev_level_ftl: float = 0.0
    horizon_penalized: bool = False


def compute_dense_reward(
    *,
    ftl_geo: float,
    l2_ham: float,
    min_lapse: float,
    state: TaxManState,
    window_k: int = 4,
    ftl_level_weight: float = 0.7,
    ftl_delta_weight: float = 0.3,
) -> float:
    state.ftl_window.append(ftl_geo)
    if len(state.ftl_window) > window_k:
        state.ftl_window.pop(0)
    level_ftl = sum(state.ftl_window) / max(len(state.ftl_window), 1)
    delta_ftl = level_ftl - state.prev_level_ftl
    state.prev_level_ftl = level_ftl
    ftl_term = ftl_level_weight * level_ftl + ftl_delta_weight * delta_ftl
    dense = (
        1000.0 * ftl_term
        - 500.0 * l2_ham
        - 50.0 * max(0.0, 0.2 - min_lapse)
    )
    state.dense_history.append(dense)
    return dense


def evaluate_fences(
    *,
    min_lapse: float,
    l2_ham: float,
    wec_violation_fraction: float | None,
    horizon_detected: bool,
    state: TaxManState,
) -> FenceResult:
    if wec_violation_fraction is not None and wec_violation_fraction > 0.15:
        return FenceResult(True, -5000.0, "wec")
    if l2_ham > 0.05:
        return FenceResult(True, 0.0, "l2_ham")
    if min_lapse < 0.05:
        return FenceResult(True, 0.0, "min_lapse")
    if horizon_detected and not state.horizon_penalized:
        state.horizon_penalized = True
        return FenceResult(True, -500.0, "horizon")
    return FenceResult(False, 0.0, "")


def baseline_dense_reward(*, l2_ham: float, min_lapse: float, max_abs_k: float, alpha: float = 0.0) -> float:
    return -500.0 * l2_ham - 50.0 * max(0.0, 0.2 - min_lapse) + alpha * max_abs_k


def frame_metrics_from_record(record: dict[str, Any]) -> dict[str, float | None]:
    return {
        "ftl_geo": float(record.get("f_geo", 0.0) or 0.0),
        "ftl_geo_evolving": float(record.get("ftl_geo_evolving", 0.0) or 0.0),
    }
