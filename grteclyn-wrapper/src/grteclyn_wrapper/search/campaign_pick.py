"""Pick elite eval ids from campaign trajectory.jsonl files."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def _pickable_record(rec: dict[str, Any], *, include_dry_run: bool = False) -> bool:
    status = rec.get("status")
    if include_dry_run and (rec.get("dry_run") or status == "dry_run"):
        return True
    if status != "gpu_ok":
        return False
    if rec.get("grtresna_rejected") or rec.get("solved_ftl_rejected"):
        return False
    if rec.get("grtresna_failed") or rec.get("preflight_rejected"):
        return False
    return True


def _gpu_ok_record(rec: dict[str, Any]) -> bool:
    return _pickable_record(rec, include_dry_run=False)


def load_gpu_ok_records(
    trajectory_path: Path,
    *,
    include_dry_run: bool = False,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    text = trajectory_path.read_text(encoding="utf-8")
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        if _pickable_record(rec, include_dry_run=include_dry_run):
            records.append(rec)
    return records


def record_rank_value(rec: dict[str, Any], rank_by: str = "score") -> float:
    if rank_by == "score":
        return float(rec.get("score") or float("-inf"))
    components = rec.get("components") or {}
    details = rec.get("descriptor_details") or {}
    for source in (details, components):
        if rank_by in source and source[rank_by] is not None:
            return float(source[rank_by])
    return float("-inf")


def pick_top_eval_ids(
    trajectory_path: Path,
    *,
    top_k: int = 1,
    rank_by: str = "score",
    include_dry_run: bool = False,
) -> list[int]:
    """Return eval ids for the top-K pickable rows sorted by ``rank_by``."""
    if top_k <= 0:
        raise ValueError(f"top_k must be positive, got {top_k}")
    records = load_gpu_ok_records(trajectory_path, include_dry_run=include_dry_run)
    if not records:
        label = "pickable" if include_dry_run else "gpu_ok"
        raise ValueError(f"no {label} records in {trajectory_path}")
    records.sort(key=lambda rec: record_rank_value(rec, rank_by=rank_by), reverse=True)
    return [int(rec["eval"]) for rec in records[:top_k]]


def pick_top_eval_id(
    trajectory_path: Path,
    *,
    rank_by: str = "score",
    include_dry_run: bool = False,
) -> int:
    return pick_top_eval_ids(
        trajectory_path,
        top_k=1,
        rank_by=rank_by,
        include_dry_run=include_dry_run,
    )[0]
