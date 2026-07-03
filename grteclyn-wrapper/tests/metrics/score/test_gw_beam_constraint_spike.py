"""Tests for constraint-spike detection and gw_beam scoring guards."""

from __future__ import annotations

from pathlib import Path

import pytest

from grteclyn_wrapper.metrics.diagnostics.constraints import (
    read_constraint_metrics,
    spike_info_from_rows,
)
from grteclyn_wrapper.metrics.diagnostics.psi4 import read_psi4_metrics
from grteclyn_wrapper.metrics.score.gw_beam import (
    GW_FATAL_HAM_L2,
    gw_beam_total,
    gw_health_multiplier,
)
from grteclyn_wrapper.metrics.score.health import compute_health_components
from grteclyn_wrapper.metrics.score.types import ScoringContext
from grteclyn_wrapper.metrics.types.diagnostics import ConstraintMetrics
from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics
from grteclyn_wrapper.metrics.types.psi4 import Psi4Metrics


def _episode(constraints: ConstraintMetrics | None) -> EpisodeMetrics:
    return EpisodeMetrics(
        collapse=None,
        constraints=constraints,
        stability=None,
        comoving=None,
        ftl=None,
        termination_reason="test",
    )


def test_spike_info_detects_single_step_ham_cliff() -> None:
    rows = [[float(t), 0.05, 0.002] for t in [i * 0.01 for i in range(1, 433)]]
    rows.append([4.33, 2.0e4, 0.34])
    rows.extend([[4.34 + i * 0.01, 1.0e4, 0.3] for i in range(100)])

    info = spike_info_from_rows(rows)
    assert info.has_constraint_spike
    assert info.constraint_spike_time == pytest.approx(4.33)
    assert info.max_step_ham_ratio > 1.0e5


def test_spike_info_clean_run() -> None:
    rows = [[float(t), 0.05, 0.002] for t in [i * 0.01 for i in range(1, 500)]]
    info = spike_info_from_rows(rows)
    assert not info.has_constraint_spike
    assert info.constraint_spike_time is None


def test_constraint_spike_penalty_demotes_gamed_score() -> None:
    base = {
        "gw_beam_quality": 0.2,
        "gw_peak_power": 0.0003,
        "gw_health_multiplier": 0.0,
        "survival": 0.0,
        "stability": 0.01,
        "constraint_health": 0.0,
        "horizon_penalty": -1.0,
        "exotic_penalty": -1.6,
        "instability_penalty": -0.99,
        "constraint_spike_penalty": -0.75,
    }
    total = gw_beam_total(base)
    assert total < 0.0
    assert 1000 * 0.2 + 100 * 0.0003 > 100  # raw GW would dominate without veto


def test_gw_health_multiplier_vetoes_fatal_ham() -> None:
    constraints = ConstraintMetrics(
        final_time=24.0,
        max_hamiltonian_l2=GW_FATAL_HAM_L2 + 1.0,
        max_momentum_l2=0.1,
        final_hamiltonian_l2=1.0,
        final_momentum_l2=0.1,
        min_rho_required=-1.0,
        max_rho_required=1.0,
        integral_negative_rho=90.0,
        has_constraint_spike=True,
        constraint_spike_time=4.5,
        ham_spike_ratio=1e6,
        max_step_ham_ratio=1e6,
        mom_spike_ratio=10.0,
    )
    assert gw_health_multiplier(_episode(constraints)) == 0.0


def test_gw_health_multiplier_vetoes_early_spike() -> None:
    constraints = ConstraintMetrics(
        final_time=24.0,
        max_hamiltonian_l2=50.0,
        max_momentum_l2=0.01,
        final_hamiltonian_l2=50.0,
        final_momentum_l2=0.01,
        min_rho_required=0.0,
        max_rho_required=1.0,
        integral_negative_rho=0.0,
        has_constraint_spike=True,
        constraint_spike_time=4.5,
        ham_spike_ratio=10.0,
        max_step_ham_ratio=10.0,
        mom_spike_ratio=2.0,
    )
    assert gw_health_multiplier(_episode(constraints)) == 0.0


def test_compute_health_sets_spike_penalty() -> None:
    constraints = ConstraintMetrics(
        final_time=24.0,
        max_hamiltonian_l2=2.0e7,
        max_momentum_l2=0.3,
        final_hamiltonian_l2=1.0e5,
        final_momentum_l2=0.2,
        min_rho_required=-1.0,
        max_rho_required=1.0,
        integral_negative_rho=90.0,
        ham_spike_ratio=1.0e6,
        max_step_ham_ratio=1.0e6,
        mom_spike_ratio=100.0,
        has_constraint_spike=True,
        constraint_spike_time=4.33,
    )
    ctx = ScoringContext(metrics=_episode(constraints), target_stop_time=24.0, domain_half_width=32.0, weights={})
    compute_health_components(ctx)
    assert ctx.components["constraint_spike_penalty"] == pytest.approx(-0.75)


def test_psi4_peak_capped_before_spike_time(tmp_path: Path) -> None:
    path = tmp_path / "psi4_directional.dat"
    path.write_text(
        "\n".join(
            [
                "0.0 0.01 0.005 0.5",
                "4.0 0.02 0.010 0.5",
                "5.0 1.50 0.800 0.6",
                "20.0 2.00 1.000 0.7",
            ]
        )
        + "\n"
    )
    full = read_psi4_metrics(tmp_path)
    capped = read_psi4_metrics(tmp_path, max_peak_time=4.0)
    assert full is not None and capped is not None
    assert full.peak_total_power == pytest.approx(2.0)
    assert capped.peak_total_power == pytest.approx(0.02)


@pytest.mark.skipif(
    not Path(
        "runs/grtresna_qd/gw_beam_qd100_v3/eval_000051/data/constraint_norms.dat"
    ).is_file(),
    reason="eval 51 data not present",
)
def test_eval_51_detects_spike() -> None:
    root = Path("runs/grtresna_qd/gw_beam_qd100_v3/eval_000051")
    metrics = read_constraint_metrics(root / "data" / "constraint_norms.dat")
    assert metrics is not None
    assert metrics.has_constraint_spike
    assert metrics.constraint_spike_time == pytest.approx(4.49, abs=0.02)
