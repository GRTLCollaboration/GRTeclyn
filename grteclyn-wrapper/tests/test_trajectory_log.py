"""Tests for trajectory.jsonl field order and eval-first logging."""

from __future__ import annotations

import json

from grteclyn_wrapper.search.trajectory_log import (
    format_eval_log_line,
    format_trajectory_line,
    infer_trajectory_status,
)


def test_trajectory_line_puts_eval_first() -> None:
    record = {
        "components": {"survival": 1.0},
        "eval": 19,
        "score": 10.47,
        "episode": "/runs/eval_000019",
        "exit_code": 0,
    }
    line = format_trajectory_line(record)
    assert line.startswith('{"eval": 19, "status": "gpu_ok"')
    parsed = json.loads(line)
    keys = list(parsed.keys())
    assert keys[0] == "eval"
    assert keys[1] == "status"
    assert parsed["status"] == "gpu_ok"


def test_infer_solved_ftl_rejected() -> None:
    record = {"eval": 22, "solved_ftl_rejected": True, "score": -50.0}
    assert infer_trajectory_status(record) == "solved_ftl_rejected"
    line = format_trajectory_line(record)
    assert '"eval": 22' in line.split(",")[0]


def test_eval_log_line_order_eval_score_status() -> None:
    line = format_eval_log_line(
        {"eval": 19, "score": 10.47, "exit_code": 0, "episode": "/e"},
    )
    assert line == "[optimize] eval 19 score=10.4700 status=gpu_ok"
    assert "episode=" not in line
