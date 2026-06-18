"""Top-K near-miss parents ranked by graded pre-GPU rejection score."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np

from ..optimize import SearchDimension
from ..optimize.candidates import _overrides_to_vector
from .types import is_graded_rejection


@dataclass
class NearMissPool:
    max_size: int = 32
    _entries: list[tuple[float, dict[str, float]]] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self._entries)

    def consider(self, record: Mapping[str, Any]) -> None:
        if self.max_size <= 0 or not is_graded_rejection(record):
            return
        score = record.get("score")
        overrides = record.get("overrides")
        if not isinstance(overrides, Mapping):
            return
        try:
            score_f = float(score)
        except (TypeError, ValueError):
            return
        if not math.isfinite(score_f):
            return
        params = {
            str(k): float(v)
            for k, v in overrides.items()
            if isinstance(v, (int, float)) and math.isfinite(float(v))
        }
        if not params:
            return
        self._entries.append((score_f, params))
        self._entries.sort(key=lambda item: item[0], reverse=True)
        if len(self._entries) > self.max_size:
            self._entries = self._entries[: self.max_size]

    def sample_params(self, rng: np.random.Generator) -> dict[str, float] | None:
        if not self._entries:
            return None
        _score, params = self._entries[int(rng.integers(len(self._entries)))]
        return dict(params)

    def param_sets(self) -> list[dict[str, float]]:
        return [dict(params) for _score, params in self._entries]

    @classmethod
    def rebuild_from_trajectory(
        cls,
        records: Sequence[Mapping[str, Any]],
        *,
        max_size: int = 32,
    ) -> NearMissPool:
        pool = cls(max_size=max_size)
        for record in records:
            pool.consider(record)
        return pool


def near_miss_vectors_from_trajectory(
    records: Sequence[Mapping[str, Any]],
    dims: Sequence[SearchDimension],
    *,
    max_size: int = 32,
) -> list[list[float]]:
    pool = NearMissPool.rebuild_from_trajectory(records, max_size=max_size)
    return [_overrides_to_vector(params, dims) for params in pool.param_sets()]
