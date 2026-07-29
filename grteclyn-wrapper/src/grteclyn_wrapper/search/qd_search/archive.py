"""MAP-Elites archive data structures."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

from ..validation_tiers import Tier, TIER_NAMES


@dataclass
class Elite:
    cell: tuple[int, int]
    score: float
    descriptors: tuple[float, float]
    params: dict[str, float]
    episode: str | None
    tier: int = int(Tier.REJECTED)
    tier_name: str = "rejected"
    objectives: dict[str, float] = field(default_factory=dict)
    descriptor_details: dict[str, float] = field(default_factory=dict)
    #: The raw optimizer vector this elite was decoded from, in genome
    #: coordinates.  ``params`` holds the *decoded* physical values, which are
    #: not the same thing on the trajectory speed axes -- reading them back as a
    #: genome is ``DebugPreGPU.md`` PG-1.  ``None`` for archives written before
    #: genomes were stored; ``sampling.elite_genome`` then reconstructs it.
    genome: list[float] | None = None


@dataclass
class QDArchive:
    bins: int
    cells: dict[tuple[int, int], Elite] = field(default_factory=dict)

    def insert(self, elite: Elite) -> bool:
        existing = self.cells.get(elite.cell)
        if existing is None or elite.score > existing.score:
            self.cells[elite.cell] = elite
            return True
        return False

    @property
    def coverage(self) -> float:
        return len(self.cells) / float(self.bins * self.bins)

    @property
    def best(self) -> Elite | None:
        if not self.cells:
            return None
        return max(self.cells.values(), key=lambda e: e.score)

    def tier_histogram(self) -> dict[str, int]:
        hist: dict[str, int] = {name: 0 for name in TIER_NAMES.values()}
        for e in self.cells.values():
            hist[e.tier_name] = hist.get(e.tier_name, 0) + 1
        return hist

    def to_dict(self) -> dict[str, Any]:
        return {
            "bins": self.bins,
            "coverage": self.coverage,
            "num_elites": len(self.cells),
            "best_score": self.best.score if self.best else None,
            "tier_histogram": self.tier_histogram(),
            "cells": [
                {
                    "cell": list(e.cell),
                    "score": e.score,
                    "descriptors": list(e.descriptors),
                    "descriptor_details": e.descriptor_details,
                    "tier": e.tier,
                    "tier_name": e.tier_name,
                    "objectives": e.objectives,
                    "params": e.params,
                    "genome": e.genome,
                    "episode": e.episode,
                }
                for e in sorted(self.cells.values(), key=lambda e: (-e.tier, -e.score))
            ],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> QDArchive:
        archive = cls(bins=int(data["bins"]))
        for cell_data in data.get("cells", []):
            cell = tuple(cell_data["cell"])
            archive.cells[cell] = Elite(
                cell=cell,
                score=float(cell_data["score"]),
                descriptors=tuple(cell_data["descriptors"]),
                params=dict(cell_data["params"]),
                genome=(
                    [float(v) for v in cell_data["genome"]]
                    if isinstance(cell_data.get("genome"), list)
                    else None
                ),
                episode=cell_data.get("episode"),
                tier=int(cell_data.get("tier", int(Tier.REJECTED))),
                tier_name=str(cell_data.get("tier_name", "rejected")),
                objectives=dict(cell_data.get("objectives", {})),
                descriptor_details=dict(cell_data.get("descriptor_details", {})),
            )
        return archive
