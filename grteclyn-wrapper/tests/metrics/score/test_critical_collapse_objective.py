"""Tests for critical_collapse objective scalarization."""

from __future__ import annotations

import pytest

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


def test_exotic_penalty_reduces_splash_score() -> None:
    clean = _critical_collapse_total(
        {
            "survival": 1.0,
            "central_energy_peak": 0.16,
            "focusing_efficiency": 1.0,
            "exotic_penalty": 0.0,
        },
        [],
        splash_mode="discovery",
    )
    exotic = _critical_collapse_total(
        {
            "survival": 1.0,
            "central_energy_peak": 0.16,
            "focusing_efficiency": 1.0,
            "exotic_penalty": -1.6,
        },
        [],
        splash_mode="discovery",
    )
    assert exotic < clean
    # Exotic penalty is capped at -0.3 and weighted 50x in critical_collapse
    # (previously uncapped at -1.6 weighted 200x = 320 pts drag).
    assert clean - exotic == 15.0


def test_geometric_splash_dominates_matter() -> None:
    """A geometric splash (chi-well + Psi4 wave) outscores pure matter pile-up."""
    geometric = _critical_collapse_total(
        {
            "survival": 1.0,
            "central_energy_peak": 0.1,
            "geometric_curvature_well": 1.0,
            "geometric_wave_arrival": 1.0,
            "geometric_crunch": 1.0,
        },
        [],
        splash_mode="discovery",
    )
    matter_only = _critical_collapse_total(
        {
            "survival": 1.0,
            "central_energy_peak": 0.1,
            "geometric_curvature_well": 0.0,
            "geometric_wave_arrival": 0.0,
            "geometric_crunch": 0.0,
        },
        [],
        splash_mode="discovery",
    )
    assert geometric > matter_only
    # Geometric terms contribute 800+600+300 = 1700 pts at full strength
    assert geometric - matter_only == pytest.approx(1700.0)


def test_geometric_splash_note_mentions_geometry() -> None:
    notes: list[str] = []
    _critical_collapse_total(
        {
            "survival": 1.0,
            "central_energy_peak": 0.2,
            "geometric_curvature_well": 0.5,
        },
        notes,
        splash_mode="discovery",
    )
    assert any("geometric splash" in n for n in notes)


def test_exotic_penalty_capped_for_critical_collapse() -> None:
    """Even a maximally exotic config (-1.6) loses at most 15 pts, not 320."""
    moderate = _critical_collapse_total(
        {
            "survival": 1.0,
            "central_energy_peak": 0.5,
            "focusing_efficiency": 2.0,
            "exotic_penalty": -0.3,
        },
        [],
        splash_mode="discovery",
    )
    extreme = _critical_collapse_total(
        {
            "survival": 1.0,
            "central_energy_peak": 0.5,
            "focusing_efficiency": 2.0,
            "exotic_penalty": -1.6,
        },
        [],
        splash_mode="discovery",
    )
    # Both capped at -0.3, so scores are identical
    assert moderate == extreme
