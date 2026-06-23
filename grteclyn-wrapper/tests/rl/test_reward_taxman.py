from __future__ import annotations

from grteclyn_wrapper.rl.audit import compute_audit_penalty
from grteclyn_wrapper.rl.reward import TaxManState, compute_dense_reward, evaluate_fences


def test_dense_reward_increases_with_ftl() -> None:
    state = TaxManState()
    low = compute_dense_reward(ftl_geo=0.05, l2_ham=0.01, min_lapse=0.3, state=state)
    high = compute_dense_reward(ftl_geo=0.20, l2_ham=0.01, min_lapse=0.3, state=TaxManState())
    assert high > low


def test_wec_fence_none_guard_inactive() -> None:
    state = TaxManState()
    result = evaluate_fences(
        min_lapse=0.2,
        l2_ham=0.01,
        wec_violation_fraction=None,
        horizon_detected=False,
        state=state,
    )
    assert not result.terminated


def test_wec_fence_triggers() -> None:
    state = TaxManState()
    result = evaluate_fences(
        min_lapse=0.2,
        l2_ham=0.01,
        wec_violation_fraction=0.2,
        horizon_detected=False,
        state=state,
    )
    assert result.terminated
    assert result.penalty == -5000.0


def test_audit_penalty_clipped() -> None:
    audit = min(0.0, -500.0 - 400.0)
    clipped = max(audit, -2000.0)
    assert clipped == -900.0
