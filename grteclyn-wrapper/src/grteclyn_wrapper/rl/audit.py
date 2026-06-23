from __future__ import annotations

import json
import time
from pathlib import Path

from grteclyn_wrapper.metrics.aggregation.collector import read_episode_metrics
from grteclyn_wrapper.metrics.score.scorer import score_episode


def wait_consumer_drain(episode_path: Path, *, timeout: float = 120.0, poll_seconds: float = 0.5) -> None:
    """Block until plot consumer state stops changing or timeout."""
    state_path = episode_path / "small_data" / "consume_state.json"
    deadline = time.monotonic() + timeout
    last_mtime = None
    stable_since = time.monotonic()
    while time.monotonic() < deadline:
        if not state_path.exists():
            time.sleep(poll_seconds)
            continue
        mtime = state_path.stat().st_mtime
        if last_mtime is None:
            last_mtime = mtime
            time.sleep(poll_seconds)
            continue
        if mtime == last_mtime:
            if time.monotonic() - stable_since >= 2.0:
                return
        else:
            last_mtime = mtime
            stable_since = time.monotonic()
        time.sleep(poll_seconds)


def wait_for_frame_record(
    episode_path: Path,
    sim_time: float,
    *,
    timeout: float = 30.0,
    poll_seconds: float = 0.25,
) -> dict | None:
    score_path = episode_path / "small_data" / "score_timeseries.jsonl"
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if score_path.exists():
            for line in score_path.read_text(encoding="utf-8").splitlines():
                if not line.strip():
                    continue
                record = json.loads(line)
                if abs(float(record.get("t", -1.0)) - sim_time) < 1e-9:
                    return record
        time.sleep(poll_seconds)
    return None


def compute_audit_penalty(
    episode_path: Path,
    accumulated_dense: float,
    *,
    objective_mode: str = "general_ftl",
    target_stop_time: float | None = None,
    clip_min: float = -2000.0,
) -> float:
    metrics = read_episode_metrics(episode_path)
    full_qd = score_episode(
        metrics,
        objective_mode=objective_mode,
        target_stop_time=target_stop_time,
    )
    audit = min(0.0, full_qd.total - accumulated_dense)
    return max(audit, clip_min)
