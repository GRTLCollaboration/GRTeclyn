"""Transient igniter: f_geo peak selection respects RL_PUMP_STOP_TIME."""

from __future__ import annotations

from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import (
    EvolvingGeodesicFtlReport,
    _eligible_emit_reports,
)
from grteclyn_wrapper.metrics.score.ftl import _post_pump_emit_ok
from grteclyn_wrapper.grtresna.matter.wiring import (
    evolution_overrides_from_complex_scalar,
)


def _rep(f_geo: float) -> EvolvingGeodesicFtlReport:
    return EvolvingGeodesicFtlReport(
        f_geo=f_geo,
        f_geo_frozen_peak=None,
        t_emit=0.0,
        t_arrival=10.0,
        t_flat=12.0,
        n_rays=3,
        n_reached=3,
        max_h_drift=0.0,
        h_quality_ok=True,
        max_h_rel_drift=0.0,
        notes=(),
    )


def test_eligible_emit_reports_filters_pre_stop(monkeypatch) -> None:
    monkeypatch.setenv("RL_PUMP_STOP_TIME", "4.0")
    reports = [
        (0.0, _rep(0.9)),
        (2.0, _rep(0.8)),
        (4.0, _rep(0.1)),
        (6.0, _rep(0.2)),
    ]
    kept = _eligible_emit_reports(reports)
    assert [te for te, _ in kept] == [4.0, 6.0]


def test_post_pump_emit_ok(monkeypatch) -> None:
    monkeypatch.delenv("RL_PUMP_STOP_TIME", raising=False)
    assert _post_pump_emit_ok(0.0) is True
    monkeypatch.setenv("RL_PUMP_STOP_TIME", "5.0")
    assert _post_pump_emit_ok(4.9) is False
    assert _post_pump_emit_ok(5.0) is True


def test_pump_controller_overrides_carry_stop_time(monkeypatch) -> None:
    monkeypatch.setenv("RL_PUMP_STOP_TIME", "5.0")
    overrides = evolution_overrides_from_complex_scalar(
        mass=1.0, lam=640.0, bs_omega=0.8, mu=85333.0
    )
    assert overrides["rl_pump_stop_time"] == 5.0
