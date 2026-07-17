"""Pump-work scoring tax and collapse diagnostic integration."""

from __future__ import annotations

from pathlib import Path

from grteclyn_wrapper.metrics.diagnostics.collapse import read_collapse_metrics
from grteclyn_wrapper.metrics.score.objectives import _general_ftl_total
from grteclyn_wrapper.metrics.score.penalties import compute_penalty_components
from grteclyn_wrapper.metrics.score.types import ScoringContext
from grteclyn_wrapper.metrics.types.diagnostics import CollapseMetrics
from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics


def _ctx_with_collapse(**kwargs) -> ScoringContext:
    collapse = CollapseMetrics(
        final_time=16.0,
        min_lapse=0.9,
        min_chi=0.9,
        max_abs_k=0.1,
        max_horizon_radius=0.0,
        min_theta_plus=0.1,
        scalar_phi_range=0.1,
        scalar_pi_range=0.1,
        **kwargs,
    )
    metrics = EpisodeMetrics(
        collapse=collapse,
        constraints=None,
        stability=None,
        comoving=None,
        ftl=None,
        termination_reason="completed",
    )
    return ScoringContext(
        metrics=metrics,
        target_stop_time=16.0,
        domain_half_width=32.0,
        weights={},
        components={},
        notes=[],
        stationary_artifact=False,
    )


def test_passive_pump_work_integrates_to_zero(tmp_path: Path) -> None:
    path = tmp_path / "collapse_diagnostics.dat"
    # RadialRecipe layout: t + 13 legacy + pump_work = 15 cols.
    rows = [
        "0.0 1 1 0 0 0 0 0 0.1 0 0 0 0 0 0.0",
        "1.0 1 1 0 0 0 0 0 0.1 0 0 0 0 0 0.0",
        "2.0 1 1 0 0 0 0 0 0.1 0 0 0 0 0 0.0",
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    m = read_collapse_metrics(path)
    assert m is not None
    assert m.pump_energy == 0.0
    assert m.pump_energy_norm == 0.0


def test_active_pump_work_produces_energy_and_penalty(tmp_path: Path) -> None:
    path = tmp_path / "collapse_diagnostics.dat"
    rows = [
        "0.0 1 1 0 0 0 0 0 0.1 0 0 0 0 0 1.0",
        "1.0 1 1 0 0 0 0 0 0.1 0 0 0 0 0 1.0",
        "2.0 1 1 0 0 0 0 0 0.1 0 0 0 0 0 1.0",
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    m = read_collapse_metrics(path)
    assert m is not None
    assert m.pump_energy is not None and m.pump_energy > 0.0
    assert m.pump_energy_norm is not None and m.pump_energy_norm > 0.0

    ctx = _ctx_with_collapse(
        pump_energy=m.pump_energy, pump_energy_norm=m.pump_energy_norm
    )
    compute_penalty_components(ctx)
    assert ctx.components["pump_energy_penalty"] < 0.0


def test_general_ftl_total_decreases_with_pump_penalty(monkeypatch) -> None:
    monkeypatch.setenv("SCORE_PUMP_ENERGY_WEIGHT", "40")
    base = {
        "nontriviality_gate": 1.0,
        "ftl_geo_evolving": 0.0,
        "operational_ftl_geodesic": 0.0,
        "survival": 0.5,
        "exotic_penalty": 0.0,
        "pump_energy_penalty": 0.0,
    }
    clean = _general_ftl_total(dict(base), [])
    taxed = _general_ftl_total({**base, "pump_energy_penalty": -1.0}, [])
    assert taxed < clean
    assert taxed == clean - 40.0
