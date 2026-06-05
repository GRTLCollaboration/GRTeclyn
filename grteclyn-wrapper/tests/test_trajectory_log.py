"""Tests for trajectory.jsonl field order and eval-first logging."""

from __future__ import annotations

import json

from grteclyn_wrapper.core.evaluation import Evaluation
from grteclyn_wrapper.search.trajectory_log import (
    format_eval_log_line,
    format_trajectory_line,
    infer_trajectory_status,
    trajectory_flags_from_evaluation,
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


def test_trajectory_flags_from_evaluation_gpu_ok() -> None:
    res = Evaluation(
        score=12.0,
        components={"survival": 1.0},
        notes=[],
        episode_path="/runs/eval_000006",
        exit_code=0,
        preflight_rejected=False,
        reason=None,
        metrics={},
    )
    assert trajectory_flags_from_evaluation(res) == {"exit_code": 0}


def test_trajectory_flags_from_evaluation_rejections() -> None:
    grtresna = Evaluation(
        score=-200.0,
        components={"grtresna_rejection": -200.0},
        notes=["Hamiltonian constraint not converged"],
        episode_path="/runs/eval_000001",
        exit_code=None,
        preflight_rejected=False,
        reason="Hamiltonian constraint not converged",
        metrics={},
    )
    assert trajectory_flags_from_evaluation(grtresna) == {
        "reason": "Hamiltonian constraint not converged",
        "grtresna_rejected": True,
    }

    solved = Evaluation(
        score=-84.0,
        components={"solved_ftl_rejection": -84.0},
        notes=[],
        episode_path="/runs/eval_000007",
        exit_code=None,
        preflight_rejected=False,
        reason="solved_ftl_rejected",
        metrics={},
    )
    assert trajectory_flags_from_evaluation(solved) == {
        "reason": "solved_ftl_rejected",
        "solved_ftl_rejected": True,
    }


def test_infer_status_from_components_fallback() -> None:
    record = {
        "eval": 7,
        "score": -84.0,
        "components": {"solved_ftl_rejection": -84.0},
    }
    assert infer_trajectory_status(record) == "solved_ftl_rejected"
    line = format_trajectory_line(record)
    assert json.loads(line)["status"] == "solved_ftl_rejected"
