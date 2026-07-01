"""Dispersion gate on the coordinate operational-FTL rewards.

Regression guard for the qball_traj_spiral_v1 failure mode: leaders banked the
full ``operational_ftl`` + ``ftl_persistence`` reward at only 28-40% confinement
because those two terms -- unlike the gauge-invariant geodesic terms -- were not
scaled by ``structural_persistence``.  A coordinate channel that only opens as
the matter disperses must be down-gated so persistent lumps + FTL out-rank
dissolving clouds.
"""

from __future__ import annotations

from grteclyn_wrapper.metrics.probes.ftl.general import GeneralFtlReport
from grteclyn_wrapper.metrics.score.ftl import compute_ftl_components
from grteclyn_wrapper.metrics.score.types import ScoringContext
from grteclyn_wrapper.metrics.types.diagnostics import (
    ConfinementMetrics,
    FtlPersistenceMetrics,
)
from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics


def _dispersed_confinement(frac: float) -> ConfinementMetrics:
    return ConfinementMetrics(
        n_frames=2,
        final_time=16.0,
        initial_rms_radius=5.0,
        final_rms_radius=10.0,
        max_rms_radius=10.0,
        initial_confined_frac=0.9,
        final_confined_frac=frac,
        min_confined_frac=frac,
        spread_ratio=2.0,
        initial_total=100.0,
        final_total=100.0,
    )


def _ctx(confinement: ConfinementMetrics | None) -> ScoringContext:
    # A strong evolved coordinate shortcut (f_op=0.1 == OP_FTL_TARGET) with a
    # sustained persistence window, but NO geodesic probe (so the coordinate
    # channel is not zeroed by the geodesic-contradiction gate).
    em = EpisodeMetrics(
        collapse=None,
        constraints=None,
        stability=None,
        comoving=None,
        ftl=None,
        termination_reason="completed",
        general_ftl_evolved=GeneralFtlReport(
            f_op=0.1,
            t_min=None,
            t_flat=1.0,
            max_local_speed=1.2,
            superluminal_fraction=0.1,
            path_offaxis=False,
            reachable=True,
            notes=(),
        ),
        ftl_persistence=FtlPersistenceMetrics(
            n_samples=3,
            f_op_min=0.1,
            f_op_median=0.1,
            f_op_last=0.1,
            max_local_speed_min=1.1,
            max_shift_max=0.2,
        ),
        confinement=confinement,
    )
    from grteclyn_wrapper.metrics.score.survival import compute_survival_components

    ctx = ScoringContext(
        metrics=em, target_stop_time=16.0, domain_half_width=None, weights={}
    )
    compute_survival_components(ctx)
    return ctx


def test_confined_channel_keeps_full_operational_ftl() -> None:
    ctx = _ctx(None)  # no confinement series -> structural_persistence == 1.0
    compute_ftl_components(ctx)
    assert ctx.components["structural_persistence"] == 1.0
    assert ctx.components["operational_ftl"] == 1.0
    assert ctx.components["ftl_persistence"] == 1.0


def test_dispersed_channel_down_gates_operational_ftl() -> None:
    ctx = _ctx(_dispersed_confinement(0.3))
    compute_ftl_components(ctx)
    assert ctx.components["structural_persistence"] == 0.3
    # Full gate (default strength 1.0): reward scaled by structural_persistence.
    assert abs(ctx.components["operational_ftl"] - 0.3) < 1e-9
    assert abs(ctx.components["ftl_persistence"] - 0.3) < 1e-9
    assert any("down-gated for dispersal" in n for n in ctx.notes)


def test_dispersion_gate_disabled_leaves_reward_ungated() -> None:
    ctx = _ctx(_dispersed_confinement(0.3))
    compute_ftl_components(ctx, dispersion_gate=0.0)
    assert ctx.components["structural_persistence"] == 0.3
    assert ctx.components["operational_ftl"] == 1.0
    assert ctx.components["ftl_persistence"] == 1.0


def test_dispersion_gate_partial_strength_interpolates() -> None:
    ctx = _ctx(_dispersed_confinement(0.4))
    compute_ftl_components(ctx, dispersion_gate=0.5)
    # multiplier = (1 - 0.5) + 0.5 * 0.4 = 0.7
    assert abs(ctx.components["operational_ftl"] - 0.7) < 1e-9
    assert abs(ctx.components["ftl_persistence"] - 0.7) < 1e-9


def test_score_episode_env_toggle(monkeypatch) -> None:
    from grteclyn_wrapper.metrics.score import score_episode

    em = _ctx(_dispersed_confinement(0.3)).metrics

    monkeypatch.setenv("SCORE_FTL_DISPERSION_GATE", "0")
    ungated = score_episode(em, objective_mode="general_ftl")
    assert ungated.components["operational_ftl"] == 1.0

    monkeypatch.setenv("SCORE_FTL_DISPERSION_GATE", "1")
    gated = score_episode(em, objective_mode="general_ftl")
    assert abs(gated.components["operational_ftl"] - 0.3) < 1e-9
    # The dispersed config scores strictly lower once the gate is active.
    assert gated.total < ungated.total
