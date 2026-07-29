"""Shadow MAP-Elites archive for graded pre-GPU rejections."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from ..grtresna_convergence_gate import GRTresnaConvergenceConfig
from ..pre_gpu.descriptors import pre_gpu_descriptor_details
from ..pre_gpu.types import is_graded_rejection
from ..validation_tiers import Tier
from .archive import Elite, QDArchive
from .descriptors import _bin_index

PRE_GPU_ARCHIVE_FILE = "pre_gpu_archive.json"


def shadow_elite_from_record(
    record: Mapping[str, Any],
    *,
    bins: int,
    convergence_config: GRTresnaConvergenceConfig | None = None,
) -> Elite | None:
    if not is_graded_rejection(record):
        return None
    score = record.get("score")
    overrides = record.get("overrides")
    if score is None or not isinstance(overrides, Mapping):
        return None
    details = pre_gpu_descriptor_details(
        record, convergence_config=convergence_config,
    )
    d1, d2 = details["x"], details["y"]
    cell = (_bin_index(d1, bins), _bin_index(d2, bins))
    params = {str(k): float(v) for k, v in overrides.items() if isinstance(v, (int, float))}
    genome = record.get("genome")
    if not isinstance(genome, list):
        genome = None  # legacy record; sampling.elite_genome inverts params instead
    return Elite(
        cell=cell,
        score=float(score),
        descriptors=(d1, d2),
        params=params,
        genome=[float(v) for v in genome] if genome is not None else None,
        episode=record.get("episode"),
        tier=int(Tier.REJECTED),
        tier_name="rejected",
        descriptor_details=details,
    )


def load_pre_gpu_archive(path: Path, *, bins: int) -> QDArchive:
    if not path.is_file():
        return QDArchive(bins=bins)
    data = json.loads(path.read_text(encoding="utf-8"))
    archive = QDArchive.from_dict(data)
    if archive.bins != bins:
        raise ValueError(
            f"pre_gpu archive bins mismatch: file has {archive.bins}, expected {bins}"
        )
    return archive


def rebuild_pre_gpu_archive_from_trajectory(
    records: Sequence[Mapping[str, Any]],
    *,
    bins: int,
    convergence_config: GRTresnaConvergenceConfig | None = None,
) -> QDArchive:
    archive = QDArchive(bins=bins)
    for record in records:
        elite = shadow_elite_from_record(
            record, bins=bins, convergence_config=convergence_config,
        )
        if elite is not None:
            archive.insert(elite)
    return archive
