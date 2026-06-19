"""Tests for critical_collapse objective scalarization."""

from __future__ import annotations

from grteclyn_wrapper.metrics.score.objectives import _critical_collapse_total


def test_discovery_high_peak_beats_pre_collapsed_static() -> None:
    dynamic = _critical_collapse_total(
        {
            "survival": 1.0,
            "central_energy_peak": 0.8,
            "focusing_efficiency": 2.0,
            "pre_collapsed_penalty": 0.0,
        },
        [],
        splash_mode="discovery",
    )
    static = _critical_collapse_total(
        {
            "survival": 1.0,
            "central_energy_peak": 0.2,
            "focusing_efficiency": 0.5,
            "pre_collapsed_penalty": -0.75,
        },
        [],
        splash_mode="discovery",
    )
    assert dynamic > static


def test_threshold_rewards_lapse_band() -> None:
    in_band = _critical_collapse_total(
        {
            "survival": 1.0,
            "central_energy_peak": 0.5,
            "focusing_efficiency": 1.0,
            "pre_collapsed_penalty": 0.0,
            "central_lapse_collapse": 0.8,
        },
        [],
        splash_mode="threshold",
    )
    out_band = _critical_collapse_total(
        {
            "survival": 1.0,
            "central_energy_peak": 0.5,
            "focusing_efficiency": 1.0,
            "pre_collapsed_penalty": 0.0,
            "central_lapse_collapse": -1.0,
        },
        [],
        splash_mode="threshold",
    )
    assert in_band > out_band


def test_survival_zero_zeros_score() -> None:
    total = _critical_collapse_total(
        {
            "survival": 0.0,
            "central_energy_peak": 1.0,
            "focusing_efficiency": 5.0,
        },
        [],
    )
    assert total == 0.0


def test_ignores_ftl_components() -> None:
    notes: list[str] = []
    total = _critical_collapse_total(
        {
            "survival": 1.0,
            "central_energy_peak": 0.5,
            "focusing_efficiency": 1.0,
            "operational_ftl_solved": 1.0,
            "shift_drive": 1.0,
        },
        notes,
    )
    assert total > 0.0
    assert any("critical_collapse" in note for note in notes)


def test_modest_focus_without_absolute_peak_scores_low() -> None:
    """Regression: v3 blobs scored ~1000 from focus alone with tiny rho."""
    total = _critical_collapse_total(
        {
            "survival": 0.5,
            "central_energy_peak": 0.115,
            "focusing_efficiency": 3.43,
            "wave_focusing_quality": 0.25,
            "collapse_lapse_progress": 0.47,
            "dispersion_penalty": -0.5,
            "pre_collapsed_penalty": -0.25,
        },
        [],
        splash_mode="discovery",
    )
    assert total < 200.0
