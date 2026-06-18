"""Trajectory and eval-dir helpers for MAP-Elites campaigns."""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path
from typing import Any, Mapping, Sequence


def _load_trajectory_records(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    records: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            records.append(json.loads(line))
    return records


def _prune_eval_dirs(
    qd_dir: Path,
    records: Sequence[Mapping[str, Any]],
    *,
    keep_eval_ids: set[int] | None = None,
    keep_top: int | None = None,
    protect_eval_ids: set[int] | None = None,
) -> int:
    """Delete completed eval dirs outside ``keep_eval_ids``.

    When ``keep_eval_ids`` is ``None``, falls back to legacy top-N-by-score if
    ``keep_top > 0``. When both are unset/non-positive, pruning is disabled.
    """
    if keep_eval_ids is None:
        if keep_top is None or keep_top <= 0:
            return 0
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
        keep_eval_ids = {eval_id for _score, eval_id in sorted(scored, reverse=True)[:keep_top]}

    protected = set(protect_eval_ids or set())
    scored_ids: set[int] = set()
    terminal_unscored_ids: set[int] = set()
    for rec in records:
        eval_id = rec.get("eval")
        if eval_id is None:
            continue
        try:
            eval_id_int = int(eval_id)
        except (TypeError, ValueError):
            continue
        score = rec.get("score")
        if isinstance(score, (int, float)):
            scored_ids.add(eval_id_int)
        else:
            # Completed trajectory rows with no score (e.g. pipeline_interrupted)
            # are dead evals and should be pruned when not in the keep set.
            terminal_unscored_ids.add(eval_id_int)

    deleted = 0
    for eval_dir in qd_dir.glob("eval_*"):
        if not eval_dir.is_dir():
            continue
        suffix = eval_dir.name.removeprefix("eval_")
        if not suffix.isdigit():
            continue
        eval_id = int(suffix)
        if eval_id in keep_eval_ids or eval_id in protected:
            continue
        prunable = eval_id in scored_ids or eval_id in terminal_unscored_ids
        if not prunable:
            # No trajectory row yet — still in flight; only protected ids are kept.
            continue
        shutil.rmtree(eval_dir, ignore_errors=True)
        deleted += 1

    if deleted:
        print(
            f"[qd] pruned {deleted} eval dirs; kept {len(keep_eval_ids)} "
            f"+ protected {sorted(protected)}",
            flush=True,
        )
    return deleted


def _iterations_for_target_evals(
    *,
    target_evals: int,
    batch: int,
    completed_evals: int = 0,
    include_init: bool = True,
) -> int:
    """Return MAP-Elites batch count needed to reach ``target_evals``."""
    if target_evals <= 0:
        return 0
    if completed_evals > 0:
        remaining = target_evals - completed_evals
        if remaining <= 0:
            return 0
        return int(math.ceil(remaining / batch))
    init = batch if include_init else 0
    if target_evals <= init:
        return 0
    return int(math.ceil((target_evals - init) / batch))
