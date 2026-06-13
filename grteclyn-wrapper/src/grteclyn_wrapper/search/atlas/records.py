"""Episode classification and record building for the failure atlas."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from ...core.episode import Episode
from ...metrics import EpisodeMetrics, dataclass_to_dict
from ...metrics.score import Score
from .config import AtlasThresholds


def classify_episode(
    metrics: EpisodeMetrics,
    score: Score,
    *,
    exit_code: int | None,
    target_stop_time: float | None,
    thresholds: AtlasThresholds = AtlasThresholds(),
) -> list[str]:
    labels: list[str] = []

    if exit_code is not None and exit_code != 0:
        labels.append("solver_failed")

    if metrics.collapse is None and metrics.constraints is None:
        labels.append("missing_diagnostics")
        return labels

    constraints = metrics.constraints
    if constraints is not None:
        max_constraint = max(
            constraints.max_hamiltonian_l2 or 0.0,
            constraints.max_momentum_l2 or 0.0,
        )
        if max_constraint >= thresholds.constraint_blowup:
            labels.append("constraint_blowup")

    collapse = metrics.collapse
    if collapse is not None:
        if collapse.min_lapse is not None and collapse.min_lapse <= thresholds.lapse_collapse:
            labels.append("lapse_collapse")
        if collapse.max_horizon_radius is not None and collapse.max_horizon_radius > thresholds.horizon_radius:
            labels.append("horizon_formed")
        if abs(score.components.get("nontrivial_geometry", 0.0)) <= thresholds.trivial_geometry:
            labels.append("trivial_geometry")

    if (
        exit_code == 0
        and "missing_diagnostics" not in labels
        and (target_stop_time is None or score.components.get("survival", 0.0) >= 1.0)
    ):
        labels.append("completed")

    return labels or ["completed"]


def build_record(
    *,
    index: int,
    episode: Episode,
    overrides: Mapping[str, Any],
    exit_code: int | None,
    metrics: EpisodeMetrics,
    score: Score,
    labels: Sequence[str],
    target_stop_time: float | None,
) -> dict[str, Any]:
    collapse = metrics.collapse
    constraints = metrics.constraints
    final_time = None
    if collapse and collapse.final_time is not None:
        final_time = collapse.final_time
    elif constraints and constraints.final_time is not None:
        final_time = constraints.final_time

    return {
        "index": index,
        "episode": str(episode.path),
        "exit_code": exit_code,
        "overrides": dict(overrides),
        "metrics": dataclass_to_dict(metrics),
        "score": dataclass_to_dict(score),
        "labels": list(labels),
        "final_time": final_time,
        "target_stop_time": target_stop_time,
    }


def summarize_records(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    from collections import Counter

    label_counts: Counter[str] = Counter()
    best_record = None
    for record in records:
        label_counts.update(record.get("labels") or [])
        if best_record is None:
            best_record = record
            continue
        score = (record.get("score") or {}).get("total")
        best_score = (best_record.get("score") or {}).get("total")
        if score is not None and (best_score is None or score > best_score):
            best_record = record

    return {
        "count": len(records),
        "label_counts": dict(sorted(label_counts.items())),
        "best": {
            "episode": best_record.get("episode") if best_record else None,
            "score_total": (best_record.get("score") or {}).get("total") if best_record else None,
            "labels": best_record.get("labels") if best_record else [],
        },
    }


def write_score(episode: Episode, metrics: EpisodeMetrics, score: Score, labels: Sequence[str]) -> None:
    from ...core.episode import write_json

    write_json(
        episode.score_path,
        {
            "labels": list(labels),
            "score": dataclass_to_dict(score),
            "metrics": dataclass_to_dict(metrics),
        },
    )
