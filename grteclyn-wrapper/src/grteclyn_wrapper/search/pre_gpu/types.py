"""Status helpers for graded pre-GPU rejections."""

from __future__ import annotations

from typing import Any, Mapping

from ..trajectory_log import infer_trajectory_status

GRADED_REJECTION_STATUSES = frozenset({
    "grtresna_rejected",
    "solved_ftl_rejected",
    "postload_rejected",
})


def is_graded_rejection(record: Mapping[str, Any]) -> bool:
    status = record.get("status")
    if not isinstance(status, str):
        status = infer_trajectory_status(record)
    return status in GRADED_REJECTION_STATUSES
