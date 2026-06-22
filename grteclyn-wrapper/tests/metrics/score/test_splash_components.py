"""Tests for splash score components."""

from __future__ import annotations

import pytest

from grteclyn_wrapper.metrics.score.splash import compute_splash_components
from grteclyn_wrapper.metrics.score.types import ScoringContext
from grteclyn_wrapper.metrics.types.central import CentralFieldMetrics
from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics


def _central(**kwargs) -> CentralFieldMetrics:
    defaults = dict(
        n_frames=3,
        t=(0.0, 1.0, 2.0),
        rho_req=(1.0e-4, 5.0e-3, 2.0e-3),
        lapse=(1.0, 0.5, 0.4),
        scalar_activity=(0.01, 0.02, 0.03),
        peak_rho_req_at_origin=5.0e-3,
        peak_rho_req_time=1.0,
        initial_rho_req_at_origin=1.0e-4,
        min_lapse_at_origin=0.4,
        wave_chromaticity=0.6,
    )
    defaults.update(kwargs)
    return CentralFieldMetrics(**defaults)


def _ctx(central: CentralFieldMetrics | None) -> ScoringContext:
    return ScoringContext(
        metrics=EpisodeMetrics(
            collapse=None,
            constraints=None,
            stability=None,
            comoving=None,
            ftl=None,
            termination_reason="completed_or_partial",
            central=central,
        ),
        target_stop_time=16.0,
        domain_half_width=8.0,
        weights={},
    )


def test_zero_central_metrics_yields_zero_splash_components() -> None:
    ctx = _ctx(None)
    compute_splash_components(ctx)
    assert ctx.components["central_energy_peak"] == 0.0
    assert ctx.components["focusing_efficiency"] == 0.0
    assert ctx.components["pre_collapsed_penalty"] == 0.0


def test_central_energy_peak_normalization() -> None:
    ctx = _ctx(_central(peak_rho_req_at_origin=1.0e-2))
    compute_splash_components(ctx)
    assert ctx.components["central_energy_peak"] == 1.0


def test_focusing_efficiency_capped() -> None:
    ctx = _ctx(_central(initial_rho_req_at_origin=1.0e-6, peak_rho_req_at_origin=1.0))
    compute_splash_components(ctx)
    assert ctx.components["focusing_efficiency"] == 5.0


def test_pre_collapsed_penalty_when_initial_rho_high() -> None:
    ctx = _ctx(_central(initial_rho_req_at_origin=1.0e-2, lapse=(0.15, 0.1, 0.05)))
    compute_splash_components(ctx, splash_mode="discovery")
    assert ctx.components["pre_collapsed_penalty"] == -0.75


def test_pre_collapsed_lapse_penalty_uses_initial_not_min() -> None:
    ctx = _ctx(_central(lapse=(1.0, 0.5, 0.05), min_lapse_at_origin=0.05))
    compute_splash_components(ctx, splash_mode="discovery")
    assert ctx.components["pre_collapsed_penalty"] == 0.0


def test_threshold_mode_lapse_band_reward() -> None:
    good = _ctx(_central(min_lapse_at_origin=0.03))
    compute_splash_components(good, splash_mode="threshold")
    assert good.components["central_lapse_collapse"] > 0.0

    bad = _ctx(_central(min_lapse_at_origin=0.005))
    compute_splash_components(bad, splash_mode="threshold")
    assert bad.components["central_lapse_collapse"] == -1.0


def test_dispersion_penalty_for_fading_blob() -> None:
    ctx = _ctx(
        _central(
            peak_rho_req_at_origin=1.0e-3,
            initial_rho_req_at_origin=3.0e-4,
            scalar_activity=(0.10, 0.08, 0.06, 0.04),
            wave_chromaticity=0.6,
        )
    )
    compute_splash_components(ctx)
    assert ctx.components["dispersion_penalty"] < 0.0
    assert ctx.components["wave_focusing_quality"] < 0.6


def test_discovery_lapse_progress() -> None:
    ctx = _ctx(_central(min_lapse_at_origin=0.10))
    compute_splash_components(ctx, splash_mode="discovery")
    assert ctx.components["collapse_lapse_progress"] == pytest.approx(2.0 / 3.0)


def test_geometric_components_zero_without_geometry() -> None:
    """No chi/K/weyl4 data (legacy) -> geometric components are 0."""
    ctx = _ctx(_central())  # no chi/trace_K/weyl4
    compute_splash_components(ctx)
    assert ctx.components["geometric_curvature_well"] == 0.0
    assert ctx.components["geometric_wave_arrival"] == 0.0
    assert ctx.components["geometric_crunch"] == 0.0


def test_geometric_curvature_well_from_chi_drop() -> None:
    """chi crushing toward 0 produces a curvature-well reward."""
    ctx = _ctx(
        _central(
            chi=(1.0, 0.7, 0.4),  # min chi 0.4 -> chi_drop 0.6 == CHI_DROP_TARGET
            trace_K=(0.0, 0.5, 1.0),
            t=(0.0, 1.0, 3.0),
            weyl4=(0.0, 5.0e-3, 1.0e-2),
        )
    )
    compute_splash_components(ctx)
    # chi_drop 0.6 / target 0.20 -> saturates at 1.0
    assert ctx.components["geometric_curvature_well"] == pytest.approx(1.0)
    # peak |K| 1.0 / target 0.35 -> 1.0
    assert ctx.components["geometric_crunch"] == pytest.approx(1.0)
    # peak |weyl4| 1e-2 / target 0.20 -> 0.05
    assert ctx.components["geometric_wave_arrival"] == pytest.approx(0.05)


def test_geometric_wave_arrival_partial() -> None:
    """A weak Weyl pulse gives a fractional wave-arrival reward."""
    ctx = _ctx(
        _central(
            chi=(1.0, 0.95, 0.9),
            trace_K=(0.0, 0.1, 0.2),
            t=(0.0, 1.0, 3.0),
            weyl4=(0.0, 2.0e-3, 5.0e-3),  # peak 5e-3 at t=3 / target 0.20 = 0.025
        )
    )
    compute_splash_components(ctx)
    assert ctx.components["geometric_wave_arrival"] == pytest.approx(0.025)
    # chi only dropped to 0.9 -> chi_drop 0.1 / 0.20
    assert ctx.components["geometric_curvature_well"] == pytest.approx(0.1 / 0.20)


def test_geometric_wave_arrival_gated_before_focus_time() -> None:
    """Early A_ij spikes before t=2 earn no wave-arrival credit."""
    ctx = _ctx(
        _central(
            chi=(1.0, 0.9, 0.85),
            trace_K=(0.0, 0.1, 0.2),
            t=(0.0, 0.5, 1.5),
            weyl4=(0.0, 0.15, 0.10),  # peak 0.15 at t=0.5 — before gate
        )
    )
    compute_splash_components(ctx)
    assert ctx.components["geometric_wave_arrival"] == 0.0
    assert any("gw_wave_gated" in note for note in ctx.notes)


def test_chi_rising_at_center_scores_zero_well() -> None:
    """Initial chi below 1 that rises during evolution must not earn well credit."""
    ctx = _ctx(_central(chi=(0.75, 0.85, 0.98), trace_K=(0.1, 0.15, 0.2)))
    compute_splash_components(ctx)
    assert ctx.components["geometric_curvature_well"] == 0.0


def test_chi_drop_empty_series_is_safe() -> None:
    metrics = _central(chi=(), trace_K=())
    assert metrics.chi_drop == 0.0


def test_horizon_bonus_gated_without_center_chi_well() -> None:
    from grteclyn_wrapper.metrics.types.diagnostics import CollapseMetrics

    ctx = ScoringContext(
        metrics=EpisodeMetrics(
            collapse=CollapseMetrics(
                corroborated_trapped=True,
                first_corroborated_time=1.7,
                final_time=16.0,
                max_abs_k=1.0,
                max_horizon_radius=10.0,
                min_chi=0.1,
                min_lapse=0.01,
                min_theta_plus=-0.5,
                r_at_min_theta_plus=5.0,
                scalar_phi_range=0.5,
                scalar_pi_range=0.5,
            ),
            constraints=None,
            stability=None,
            comoving=None,
            ftl=None,
            termination_reason="completed_or_partial",
            central=_central(chi=(0.75, 0.85, 0.98), trace_K=(0.1, 0.15, 0.2)),
        ),
        target_stop_time=16.0,
        domain_half_width=8.0,
        weights={},
    )
    compute_splash_components(ctx)
    assert ctx.components["horizon_formation_time"] == 0.0
    assert any("horizon_bonus_gated" in note for note in ctx.notes)
