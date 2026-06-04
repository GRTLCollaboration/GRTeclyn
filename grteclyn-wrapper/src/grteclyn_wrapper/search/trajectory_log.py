"""Human-readable trajectory.jsonl serialization and eval logging."""

from __future__ import annotations

import json
from typing import Any, Mapping


def infer_trajectory_status(record: Mapping[str, Any]) -> str:
    """Short status label for grep and log lines."""
    if record.get("surrogate_predicted"):
        return "surrogate_skipped"
    if record.get("preflight_rejected"):
        return "preflight_rejected"
    if record.get("grtresna_failed"):
        return "grtresna_failed"
    if record.get("grtresna_rejected"):
        return "grtresna_rejected"
    if record.get("solved_ftl_rejected"):
        return "solved_ftl_rejected"
    exit_code = record.get("exit_code")
    if exit_code is not None:
        return "gpu_ok" if exit_code == 0 else "gpu_failed"
    if record.get("dry_run"):
        return "dry_run"
    return "unknown"


# Keys emitted early in each JSON line (eval first for scanning).
_TRAJECTORY_HEAD_KEYS = (
    "eval",
    "status",
    "episode",
    "score",
    "fitness",
    "exit_code",
    "solved_ftl_rejected",
    "grtresna_rejected",
    "grtresna_failed",
    "preflight_rejected",
    "surrogate_predicted",
    "reason",
    "components",
)


def format_trajectory_line(record: Mapping[str, Any]) -> str:
    """Serialize one trajectory record with eval first (never sort_keys)."""
    status = record.get("status")
    if status is None:
        status = infer_trajectory_status(record)

    ordered: dict[str, Any] = {}
    if "eval" in record:
        ordered["eval"] = record["eval"]
    ordered["status"] = status
    for key in _TRAJECTORY_HEAD_KEYS:
        if key in ("eval", "status"):
            continue
        if key in record:
            ordered[key] = record[key]
    for key, value in record.items():
        if key not in ordered:
            ordered[key] = value
    return json.dumps(ordered, ensure_ascii=False) + "\n"


def format_eval_log_prefix(record: Mapping[str, Any], *, tag: str = "optimize") -> str:
    """Console prefix with eval id first: ``[optimize] eval 19:``."""
    ev = record.get("eval")
    if ev is None:
        return f"[{tag}] eval surrogate:"
    return f"[{tag}] eval {ev}:"


def format_eval_log_line(record: Mapping[str, Any], *, tag: str = "optimize") -> str:
    """One-line summary for stdout (eval number first)."""
    status = record.get("status") or infer_trajectory_status(record)
    parts = [format_eval_log_prefix(record, tag=tag).rstrip(":"), f"status={status}"]
    score = record.get("score")
    if score is not None:
        try:
            parts.append(f"score={float(score):.4f}")
        except (TypeError, ValueError):
            parts.append(f"score={score}")
    if record.get("episode"):
        parts.append(f"episode={record['episode']}")
    reason = record.get("reason")
    if reason:
        text = str(reason)
        if len(text) > 120:
            text = text[:117] + "..."
        parts.append(f"reason={text}")
    return " ".join(parts)
