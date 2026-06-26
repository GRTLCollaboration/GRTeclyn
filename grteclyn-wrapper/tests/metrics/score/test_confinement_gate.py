"""The confinement gate is the trustworthy 'matter dispersed' detector.

Regression guard: peak/total density is spatially blind (it can RISE under pump
injection while the lump sprays apart), so survival must be gated by the
mass-weighted confined fraction from confinement.dat.
"""

from __future__ import annotations

from grteclyn_wrapper.metrics.diagnostics.confinement import read_confinement_metrics
from grteclyn_wrapper.metrics.types.diagnostics import ConfinementMetrics
from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics
from grteclyn_wrapper.metrics.score.survival import compute_survival_components
from grteclyn_wrapper.metrics.score.types import ScoringContext


def _ctx(confinement: ConfinementMetrics | None) -> ScoringContext:
    em = EpisodeMetrics(
        collapse=None,
        constraints=None,
        stability=None,
        comoving=None,
        ftl=None,
        termination_reason="x",
        confinement=confinement,
    )
    return ScoringContext(
        metrics=em, target_stop_time=16.0, domain_half_width=None, weights={}
    )


def test_confinement_absent_leaves_persistence_ungated() -> None:
    ctx = _ctx(None)
    compute_survival_components(ctx)
    assert ctx.components["confinement_retention"] == 1.0
    assert ctx.components["structural_persistence"] == 1.0


def test_dispersed_matter_gates_persistence_and_notes() -> None:
    # 40% confined, spread almost doubled: clearly dispersed.
    cm = ConfinementMetrics(
        n_frames=2,
        final_time=3.2,
        initial_rms_radius=4.78,
        final_rms_radius=7.04,
        max_rms_radius=7.04,
        initial_confined_frac=0.83,
        final_confined_frac=0.40,
        min_confined_frac=0.40,
        spread_ratio=7.04 / 4.78,
        initial_total=208.0,
        final_total=1499.0,
    )
    ctx = _ctx(cm)
    compute_survival_components(ctx)
    assert ctx.components["confinement_retention"] == 0.40
    # density_retention is 1.0 (no constraints series), so persistence == retention.
    assert ctx.components["structural_persistence"] == 0.40
    assert any("DISPERSED" in n for n in ctx.notes)


def test_reader_round_trip(tmp_path) -> None:
    p = tmp_path / "confinement.dat"
    p.write_text(
        "# time  total_activity  peak_activity  rms_radius  confined_frac  "
        "bary_x  bary_y  bary_z  r_conf\n"
        "0.0  208.0  0.16  4.78  0.83  34.0  32.0  29.0  6.0\n"
        "3.2  1499.0  0.22  7.04  0.40  33.0  32.0  33.0  6.0\n",
        encoding="utf-8",
    )
    cm = read_confinement_metrics(p)
    assert cm is not None
    assert cm.n_frames == 2
    assert abs(cm.final_confined_frac - 0.40) < 1e-9
    assert abs(cm.spread_ratio - 7.04 / 4.78) < 1e-6
