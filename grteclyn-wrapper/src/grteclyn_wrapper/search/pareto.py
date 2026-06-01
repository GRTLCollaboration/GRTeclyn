"""Multi-objective Pareto-front extraction over evaluated candidates.

The scalarized score ``S = sum_k w_k s_k`` forces hand-tuned trade-offs and is
the root cause of the flat-space-attractor weighting hacks.  This module reads
an optimizer ``trajectory.jsonl`` (or a QD ``archive.json``) and returns the
non-dominated set over the conflicting objectives, so the FTL-vs-exotic-energy
trade surface can be reported without committing to weights.

Objectives (all framed as *maximize*):
  * ``ftl``          = ftl_shortcut           (more shortcut is better)
  * ``anec``         = anec_condition          (less exotic energy is better)
  * ``stability``    = constraint_growth       (slower growth is better)
  * ``tidal``        = tidal_comfort           (lower curvature is better)
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

DEFAULT_OBJECTIVES: tuple[str, ...] = (
    "ftl_shortcut",
    "anec_condition",
    "constraint_growth",
    "tidal_comfort",
)


@dataclass(frozen=True)
class ParetoPoint:
    label: str
    objectives: dict[str, float]
    total: float
    episode: str | None


def dominates(a: Mapping[str, float], b: Mapping[str, float], keys: Sequence[str]) -> bool:
    """True if ``a`` Pareto-dominates ``b`` (maximization on every key)."""
    at_least_as_good = all(a.get(k, 0.0) >= b.get(k, 0.0) for k in keys)
    strictly_better = any(a.get(k, 0.0) > b.get(k, 0.0) for k in keys)
    return at_least_as_good and strictly_better


def pareto_front(
    points: Sequence[ParetoPoint],
    *,
    objectives: Sequence[str] = DEFAULT_OBJECTIVES,
) -> list[ParetoPoint]:
    front: list[ParetoPoint] = []
    for p in points:
        dominated = any(dominates(q.objectives, p.objectives, objectives) for q in points if q is not p)
        if not dominated:
            front.append(p)
    front.sort(key=lambda p: -p.total)
    return front


def _points_from_records(
    records: Iterable[Mapping],
    *,
    objectives: Sequence[str],
) -> list[ParetoPoint]:
    points: list[ParetoPoint] = []
    for i, rec in enumerate(records):
        components = rec.get("components") or {}
        if not components:
            continue
        obj = {k: float(components.get(k, 0.0)) for k in objectives}
        points.append(ParetoPoint(
            label=str(rec.get("eval", rec.get("episode", i))),
            objectives=obj,
            total=float(rec.get("score", 0.0)),
            episode=rec.get("episode"),
        ))
    return points


def load_trajectory_points(
    path: Path,
    *,
    objectives: Sequence[str] = DEFAULT_OBJECTIVES,
) -> list[ParetoPoint]:
    """Read points from an optimizer ``trajectory.jsonl``.

    Records lacking ``components`` (e.g. preflight-rejected) are skipped.  When
    the trajectory rows omit component breakdowns, the companion ``score.json``
    files are read from each episode directory if available.
    """
    records: list[dict] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        if not rec.get("components") and rec.get("episode"):
            score_json = Path(rec["episode"]) / "score.json"
            if score_json.exists():
                payload = json.loads(score_json.read_text(encoding="utf-8"))
                rec["components"] = payload.get("score", {}).get("components", {})
        records.append(rec)
    return _points_from_records(records, objectives=objectives)


def front_to_dict(front: Sequence[ParetoPoint]) -> dict:
    return {
        "num_points": len(front),
        "points": [
            {
                "label": p.label,
                "total": p.total,
                "objectives": p.objectives,
                "episode": p.episode,
            }
            for p in front
        ],
    }
