"""Tests for splash early termination predicates."""

from __future__ import annotations

from grteclyn_wrapper.metrics.splash_early_term import evaluate_splash_early_term


def test_collapse_complete_when_lapse_drops() -> None:
    decision = evaluate_splash_early_term(
        t=4.0,
        rho=0.1,
        lapse=0.004,
        activity=0.2,
        peak_rho_so_far=0.1,
    )
    assert decision.should_stop
    assert decision.reason == "collapse_complete"


def test_no_stop_for_successful_lapse_evolution() -> None:
    decision = evaluate_splash_early_term(
        t=4.0,
        rho=0.01,
        lapse=0.05,
        activity=0.2,
        peak_rho_so_far=0.01,
        previous_activity=0.25,
    )
    assert not decision.should_stop


def test_dispersion_complete_after_peak() -> None:
    decision = evaluate_splash_early_term(
        t=8.0,
        rho=0.001,
        lapse=0.8,
        activity=0.01,
        peak_rho_so_far=0.05,
        previous_activity=0.02,
    )
    assert decision.should_stop
    assert decision.reason == "dispersion_complete"


def test_numerical_failure_on_nan() -> None:
    decision = evaluate_splash_early_term(
        t=2.0,
        rho=float("nan"),
        lapse=0.9,
        activity=0.1,
        peak_rho_so_far=0.01,
    )
    assert decision.should_stop
    assert decision.reason == "numerical_failure"
