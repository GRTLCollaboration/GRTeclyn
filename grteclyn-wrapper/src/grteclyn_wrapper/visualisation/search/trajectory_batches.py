"""Load QD trajectory records and aggregate per-GPU-batch statistics."""

from __future__ import annotations

import json
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class BatchStats:
    batch: int
    eval_lo: int
    eval_hi: int
    batch_best: float
    cum_best: float
    delta: float
    cells_improved: int
    gpu_ok: int


def load_trajectory(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    if not records:
        raise ValueError(f"no records in {path}")
    return records


def batch_stats_from_records(records: list[dict[str, Any]], batch_size: int) -> list[BatchStats]:
    if batch_size < 1:
        raise ValueError("batch_size must be >= 1")

    rows: list[BatchStats] = []
    cum_best: float | None = None

    for offset in range(0, len(records), batch_size):
        chunk = records[offset : offset + batch_size]
        scores = [float(r["score"]) for r in chunk if r.get("score") is not None]
        if not scores:
            continue

        batch_best = max(scores)
        if cum_best is None:
            delta = 0.0
            cum_best = batch_best
        else:
            delta = max(0.0, batch_best - cum_best)
            cum_best = max(cum_best, batch_best)

        rows.append(
            BatchStats(
                batch=offset // batch_size + 1,
                eval_lo=int(chunk[0]["eval"]),
                eval_hi=int(chunk[-1]["eval"]),
                batch_best=batch_best,
                cum_best=cum_best,
                delta=delta,
                cells_improved=sum(1 for r in chunk if r.get("improved")),
                gpu_ok=sum(1 for r in chunk if r.get("status") == "gpu_ok"),
            )
        )
    return rows


def batch_stats_from_campaign(campaign_dir: Path, batch_size: int = 8) -> list[BatchStats]:
    traj_path = campaign_dir / "trajectory.jsonl"
    if not traj_path.is_file():
        raise FileNotFoundError(f"missing {traj_path}")
    return batch_stats_from_records(load_trajectory(traj_path), batch_size)


def rolling_mean(values: list[float], window: int) -> list[float]:
    if window < 1:
        raise ValueError("window must be >= 1")
    out: list[float] = []
    for i in range(len(values)):
        start = max(0, i - window + 1)
        out.append(statistics.mean(values[start : i + 1]))
    return out


def format_batch_summary(rows: list[BatchStats]) -> str:
    deltas = [r.delta for r in rows]
    pos = [d for d in deltas if d > 0]
    third = max(1, len(rows) // 3)

    def _slice_mean(chunk: list[BatchStats]) -> tuple[float, int]:
        ds = [r.delta for r in chunk]
        improving = sum(1 for d in ds if d > 0)
        return statistics.mean(ds), improving

    early = _slice_mean(rows[:third])
    mid = _slice_mean(rows[third : 2 * third])
    late = _slice_mean(rows[2 * third :])

    last_improve = max((r.batch for r in rows if r.delta > 0), default=0)
    lines = [
        f"batches={len(rows)}  evals={rows[0].eval_lo}-{rows[-1].eval_hi}  "
        f"best={rows[-1].cum_best:.1f}",
        f"improving batches: {len(pos)}/{len(rows)} ({100 * len(pos) / len(rows):.0f}%)",
        f"last archive gain: batch {last_improve}" if last_improve else "no archive gains",
        f"early  mean Δ={early[0]:.1f}  improving={early[1]}/{third}",
        f"mid    mean Δ={mid[0]:.1f}  improving={mid[1]}/{third}",
        f"late   mean Δ={late[0]:.1f}  improving={late[1]}/{len(rows) - 2 * third}",
    ]
    return "\n".join(lines)
