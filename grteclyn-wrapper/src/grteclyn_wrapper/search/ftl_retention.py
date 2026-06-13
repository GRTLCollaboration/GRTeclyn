"""FTL champion retention: keep eval dirs that hold peak-metric records."""

from __future__ import annotations

import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .ftl_peak_metrics import FtlPeakMetrics, GEO_FTL_FLOOR, from_trajectory_record

FTL_RETENTION_LOG = "ftl_retention.jsonl"
FTL_CHAMPIONS_FILE = "ftl_champions.json"


@dataclass(frozen=True)
class FtlMetricSpec:
    key: str
    extract: Callable[[FtlPeakMetrics], float]
    min_value: float
    time_field: str | None = None


@dataclass
class FtlChampion:
    eval_id: int
    value: float
    episode: str | None = None
    score: float | None = None
    t_at_peak: float | None = None


@dataclass(frozen=True)
class FtlRetentionEvent:
    event: str
    metric: str
    eval_id: int
    value: float
    score: float | None
    t_at_peak: float | None
    replaced_eval: int | None
    replaced_value: float | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "event": self.event,
            "metric": self.metric,
            "eval": self.eval_id,
            "value": self.value,
            "score": self.score,
            "t_at_peak": self.t_at_peak,
            "replaced_eval": self.replaced_eval,
            "replaced_value": self.replaced_value,
        }


DEFAULT_FTL_METRIC_SPECS: tuple[FtlMetricSpec, ...] = (
    FtlMetricSpec("f_geo_peak", lambda p: p.f_geo_peak, GEO_FTL_FLOOR, "t_at_f_geo_peak"),
    FtlMetricSpec("f_op_peak", lambda p: p.f_op_peak, GEO_FTL_FLOOR, "t_at_f_op_peak"),
    FtlMetricSpec("max_local_speed", lambda p: p.max_local_speed_peak, 1.0, "t_at_max_speed"),
    FtlMetricSpec(
        "superluminal_fraction",
        lambda p: p.superluminal_fraction_peak,
        0.0,
        "t_at_superluminal_peak",
    ),
    FtlMetricSpec("ftl_lifetime_fraction", lambda p: p.ftl_lifetime_fraction, 0.0),
    FtlMetricSpec("ftl_geo_timeavg", lambda p: p.ftl_geo_timeavg, 0.0),
)


@dataclass
class FtlChampionBoard:
    """One champion eval per FTL metric."""

    champions: dict[str, FtlChampion] = field(default_factory=dict)
    specs: tuple[FtlMetricSpec, ...] = DEFAULT_FTL_METRIC_SPECS

    def champion_eval_ids(self) -> set[int]:
        return {c.eval_id for c in self.champions.values()}

    def consider(self, record: Mapping[str, Any]) -> list[FtlRetentionEvent]:
        """Try to crown ``record`` on each metric; emit events for new records."""
        peaks = from_trajectory_record(record)
        if peaks is None:
            return []

        eval_id = record.get("eval")
        if not isinstance(eval_id, int):
            return []

        episode = record.get("episode")
        episode_str = str(episode) if isinstance(episode, str) else None
        score = record.get("score")
        score_f = float(score) if isinstance(score, (int, float)) else None

        events: list[FtlRetentionEvent] = []
        for spec in self.specs:
            value = spec.extract(peaks)
            if value <= spec.min_value:
                continue

            incumbent = self.champions.get(spec.key)
            if incumbent is not None and value <= incumbent.value:
                continue

            t_at_peak = None
            if spec.time_field is not None:
                raw = getattr(peaks, spec.time_field, None)
                t_at_peak = float(raw) if raw is not None else None

            replaced_eval = incumbent.eval_id if incumbent else None
            replaced_value = incumbent.value if incumbent else None
            self.champions[spec.key] = FtlChampion(
                eval_id=eval_id,
                value=value,
                episode=episode_str,
                score=score_f,
                t_at_peak=t_at_peak,
            )
            events.append(
                FtlRetentionEvent(
                    event="crowned",
                    metric=spec.key,
                    eval_id=eval_id,
                    value=value,
                    score=score_f,
                    t_at_peak=t_at_peak,
                    replaced_eval=replaced_eval,
                    replaced_value=replaced_value,
                )
            )
        return events

    @classmethod
    def rebuild(cls, records: Sequence[Mapping[str, Any]]) -> FtlChampionBoard:
        """Replay trajectory in eval order without emitting retention events."""
        board = cls()
        ordered = sorted(
            (r for r in records if r.get("status") == "gpu_ok"),
            key=lambda r: int(r["eval"]),
        )
        for record in ordered:
            board.consider(record)
        return board

    def to_snapshot(self) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for key, champ in self.champions.items():
            out[key] = {
                "eval": champ.eval_id,
                "value": champ.value,
                "episode": champ.episode,
                "score": champ.score,
                "t_at_peak": champ.t_at_peak,
            }
        return out

    @classmethod
    def from_snapshot(cls, data: Mapping[str, Any]) -> FtlChampionBoard:
        board = cls()
        for key, raw in data.items():
            if not isinstance(raw, Mapping):
                continue
            eval_id = raw.get("eval")
            value = raw.get("value")
            if not isinstance(eval_id, int) or not isinstance(value, (int, float)):
                continue
            board.champions[str(key)] = FtlChampion(
                eval_id=eval_id,
                value=float(value),
                episode=str(raw["episode"]) if raw.get("episode") else None,
                score=float(raw["score"]) if isinstance(raw.get("score"), (int, float)) else None,
                t_at_peak=float(raw["t_at_peak"]) if isinstance(raw.get("t_at_peak"), (int, float)) else None,
            )
        return board


def compute_keep_eval_ids(
    records: Sequence[Mapping[str, Any]],
    *,
    keep_top_score: int,
    board: FtlChampionBoard | None,
    ftl_retention_enabled: bool,
    protect_eval_ids: set[int],
) -> set[int] | None:
    """Return eval IDs to retain on disk, or ``None`` to skip pruning."""
    if keep_top_score <= 0 and not ftl_retention_enabled:
        return None

    keep_ids = set(protect_eval_ids)

    if keep_top_score > 0:
        scored: list[tuple[float, int]] = []
        for rec in records:
            score = rec.get("score")
            eval_id = rec.get("eval")
            if not isinstance(score, (int, float)) or eval_id is None:
                continue
            try:
                scored.append((float(score), int(eval_id)))
            except (TypeError, ValueError):
                continue
        keep_ids.update(eval_id for _score, eval_id in sorted(scored, reverse=True)[:keep_top_score])

    if ftl_retention_enabled and board is not None:
        keep_ids.update(board.champion_eval_ids())

    return keep_ids


def load_ftl_champions(path: Path) -> FtlChampionBoard | None:
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(data, dict):
        return None
    return FtlChampionBoard.from_snapshot(data)


def save_ftl_champions(path: Path, board: FtlChampionBoard) -> None:
    path.write_text(json.dumps(board.to_snapshot(), indent=2) + "\n", encoding="utf-8")


def append_ftl_retention_events(path: Path, events: Sequence[FtlRetentionEvent]) -> None:
    if not events:
        return
    with path.open("a", encoding="utf-8") as fh:
        for event in events:
            fh.write(json.dumps(event.to_dict(), sort_keys=True) + "\n")
