from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from grteclyn_wrapper.metrics import read_episode_metrics, score_episode


def _write_episode(root: Path, rows: list[tuple[float, float, float, float]], areal: list[tuple[float, float]]) -> None:
    data = root / "data"
    small = root / "small_data"
    data.mkdir(parents=True)
    small.mkdir(parents=True)
    with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as handle:
        for t, lapse, chi, max_k in rows:
            handle.write(
                f"{t:g} {lapse:g} {chi:g} {max_k:g} "
                "0 0 0 0 0 0 0 0 0 0\n"
            )
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as handle:
        for t, *_ in rows:
            handle.write(f"{t:g} 1e-4 1e-4 0 0 0\n")
    with (small / "areal_radius.dat").open("w", encoding="utf-8") as handle:
        for t, radius in areal:
            handle.write(f"{t:g} {radius:g} 0.0625\n")


def test_stability_score_distinguishes_static_from_collapsing() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        stable = root / "stable"
        collapsing = root / "collapsing"

        _write_episode(
            stable,
            rows=[
                (0.0, 1.0, 0.9, 0.05),
                (1.0, 0.99, 0.895, 0.055),
                (2.0, 0.98, 0.89, 0.06),
            ],
            areal=[(0.0, 1.0), (1.0, 1.01), (2.0, 1.0)],
        )
        _write_episode(
            collapsing,
            rows=[
                (0.0, 1.0, 0.9, 0.05),
                (1.0, 0.55, 0.55, 0.7),
                (2.0, 0.2, 0.3, 1.2),
            ],
            areal=[(0.0, 1.0), (1.0, 0.7), (2.0, 0.4)],
        )

        stable_score = score_episode(read_episode_metrics(stable), target_stop_time=2.0)
        collapsing_score = score_episode(read_episode_metrics(collapsing), target_stop_time=2.0)

        assert stable_score.components["stability"] > 0.8
        assert collapsing_score.components["stability"] < 0.2
        assert stable_score.total > collapsing_score.total


def _write_episode_with_rho(
    root: Path,
    rows: list[tuple[float, float, float, float]],
    peak_rho: list[float],
) -> None:
    data = root / "data"
    data.mkdir(parents=True)
    with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as handle:
        for t, lapse, chi, max_k in rows:
            handle.write(f"{t:g} {lapse:g} {chi:g} {max_k:g} 0 0 0 0 0 0 0 0 0 0\n")
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as handle:
        for (t, *_), rho in zip(rows, peak_rho):
            handle.write(f"{t:g} 1e-4 1e-4 0 {rho:g} 0\n")


def test_survival_penalises_dissipated_structure() -> None:
    """Both runs reach the stop time; only the dissipated one loses peak rho.

    The integrator-only ``survival`` used to score both 1.0.  Gating it by
    structural persistence must drop the dissipated run while leaving the run
    that holds its structure near 1.0.
    """
    rows = [(0.0, 1.0, 0.9, 0.05), (1.0, 1.0, 0.9, 0.05), (2.0, 1.0, 0.9, 0.05)]
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        persistent = root / "persistent"
        dissipated = root / "dissipated"

        _write_episode_with_rho(persistent, rows, peak_rho=[1.0e-3, 1.0e-3, 1.0e-3])
        _write_episode_with_rho(dissipated, rows, peak_rho=[1.0e-3, 5.0e-4, 1.0e-4])

        persistent_score = score_episode(read_episode_metrics(persistent), target_stop_time=2.0)
        dissipated_score = score_episode(read_episode_metrics(dissipated), target_stop_time=2.0)

        assert persistent_score.components["numerical_survival"] == 1.0
        assert dissipated_score.components["numerical_survival"] == 1.0

        assert persistent_score.components["survival"] == 1.0
        assert dissipated_score.components["survival"] == 0.1
        assert dissipated_score.components["structural_persistence"] == 0.1


def _episode_with_ftl(final_peak_rho: float):
    from grteclyn_wrapper.metrics.probes.ftl.general import GeneralFtlReport
    from grteclyn_wrapper.metrics.types.diagnostics import ConstraintMetrics
    from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics

    # Tilted cones + frame-drag but no genuine shortcut (t_min == t_flat, f_op=0),
    # so operational_ftl stays 0 and the persistence gate is what we exercise.
    report = GeneralFtlReport(
        f_op=0.0,
        t_min=1.0,
        t_flat=1.0,
        max_local_speed=1.2,
        superluminal_fraction=0.3,
        path_offaxis=False,
        reachable=True,
        notes=(),
        max_shift=0.5,
    )
    constraints = ConstraintMetrics(
        final_time=2.0,
        max_hamiltonian_l2=1.0e-4,
        max_momentum_l2=1.0e-4,
        final_hamiltonian_l2=1.0e-4,
        final_momentum_l2=1.0e-4,
        min_rho_required=0.0,
        max_rho_required=1.0,
        integral_negative_rho=0.0,
        final_peak_rho_required=final_peak_rho,
        initial_peak_rho_required=1.0,
    )
    return EpisodeMetrics(
        collapse=None,
        constraints=constraints,
        stability=None,
        comoving=None,
        ftl=None,
        termination_reason="",
        general_ftl_evolved=report,
    )


def test_ftl_shaping_rewards_scale_with_persistence() -> None:
    """A geometry that fragments/dissipates must not bank full cone-tilt credit.

    The shaping rewards (ftl_precursor / channel_progress / shift_drive) are
    gated by structural persistence, so halving the retained peak energy density
    must halve each of them relative to the fully-persistent case.
    """
    full = score_episode(
        _episode_with_ftl(1.0), target_stop_time=2.0, objective_mode="ftl_first"
    )
    half = score_episode(
        _episode_with_ftl(0.5), target_stop_time=2.0, objective_mode="ftl_first"
    )

    assert full.components["structural_persistence"] == 1.0
    assert half.components["structural_persistence"] == 0.5

    for key in ("ftl_precursor", "channel_progress", "shift_drive"):
        assert full.components[key] > 0.0, key
        assert abs(half.components[key] - 0.5 * full.components[key]) < 1e-9, key

    # The fragmented (half-persistence) candidate must rank strictly below the
    # coherent one once the shaping credit is discounted.
    assert half.total < full.total


def _episode_with_solved_ftl(max_local_speed: float, superluminal_fraction: float):
    from grteclyn_wrapper.metrics.probes.ftl.general import GeneralFtlReport
    from grteclyn_wrapper.metrics.types.diagnostics import ConstraintMetrics
    from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics

    # Strong solved-FTL signal (f_op saturates the raw reward to 1.0); the gate
    # must then key off the geometry of the superluminal region, not f_op alone.
    solved = GeneralFtlReport(
        f_op=0.05,
        t_min=0.95,
        t_flat=1.0,
        max_local_speed=max_local_speed,
        superluminal_fraction=superluminal_fraction,
        path_offaxis=False,
        reachable=True,
        notes=(),
        max_shift=0.5,
    )
    constraints = ConstraintMetrics(
        final_time=2.0,
        max_hamiltonian_l2=1.0e-4,
        max_momentum_l2=1.0e-4,
        final_hamiltonian_l2=1.0e-4,
        final_momentum_l2=1.0e-4,
        min_rho_required=0.0,
        max_rho_required=1.0,
        integral_negative_rho=0.0,
        final_peak_rho_required=1.0,
    )
    return EpisodeMetrics(
        collapse=None,
        constraints=constraints,
        stability=None,
        comoving=None,
        ftl=None,
        termination_reason="",
        general_ftl_solved=solved,
    )


def test_solved_ftl_gated_against_uniform_coordinate_offset() -> None:
    """A near-uniform, marginally-superluminal field must not bank solved FTL.

    A genuine warp channel is localized (small superluminal fraction) with a
    real peak above c; a global lapse/coordinate offset has the whole slice at
    ~1.0x and trivially saturates superluminal_fraction. The localization +
    peak-margin gate must collapse the former and keep the latter.
    """
    artifact = score_episode(
        _episode_with_solved_ftl(max_local_speed=1.017, superluminal_fraction=1.0),
        target_stop_time=2.0,
        objective_mode="ftl_first",
    )
    channel = score_episode(
        _episode_with_solved_ftl(max_local_speed=1.20, superluminal_fraction=0.1),
        target_stop_time=2.0,
        objective_mode="ftl_first",
    )

    # The uniform offset is gated essentially to zero; the localized bubble keeps
    # most of its credit.
    assert artifact.components["operational_ftl_solved"] < 1e-6
    assert channel.components["operational_ftl_solved"] > 0.5
    assert channel.total > artifact.total


def _episode_with_coherence(structure_coherence: float | None):
    from grteclyn_wrapper.metrics.probes.ftl.general import GeneralFtlReport
    from grteclyn_wrapper.metrics.types.diagnostics import ConstraintMetrics
    from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics

    report = GeneralFtlReport(
        f_op=0.0,
        t_min=1.0,
        t_flat=1.0,
        max_local_speed=1.2,
        superluminal_fraction=0.1,
        path_offaxis=False,
        reachable=True,
        notes=(),
        max_shift=0.5,
        structure_coherence=structure_coherence,
    )
    constraints = ConstraintMetrics(
        final_time=2.0,
        max_hamiltonian_l2=1.0e-4,
        max_momentum_l2=1.0e-4,
        final_hamiltonian_l2=1.0e-4,
        final_momentum_l2=1.0e-4,
        min_rho_required=0.0,
        max_rho_required=1.0,
        integral_negative_rho=0.0,
        final_peak_rho_required=1.0,  # full density retention; only coherence varies
    )
    return EpisodeMetrics(
        collapse=None,
        constraints=constraints,
        stability=None,
        comoving=None,
        ftl=None,
        termination_reason="",
        general_ftl_evolved=report,
    )


def test_fragmentation_gates_survival_despite_full_density_retention() -> None:
    """A fragmented end-state keeps its peak density but must lose survival credit.

    Density retention alone reads 1.0 for a bubble that shatters into two equally
    dense lobes (peak rho is unchanged).  Folding morphological coherence into
    ``structural_persistence`` must drop survival/shaping for the fragmented case
    while leaving the single-lump case untouched.
    """
    coherent = score_episode(
        _episode_with_coherence(1.0), target_stop_time=2.0, objective_mode="ftl_first"
    )
    fragmented = score_episode(
        _episode_with_coherence(0.5), target_stop_time=2.0, objective_mode="ftl_first"
    )

    assert coherent.components["structural_persistence"] == 1.0
    assert fragmented.components["structural_persistence"] == 0.5
    assert fragmented.components["survival"] == 0.5
    # Shaping rewards are persistence-gated too, so the fragmented candidate ranks
    # strictly lower overall.
    assert fragmented.total < coherent.total


def test_structure_coherence_counts_comparable_lobes() -> None:
    """One blob -> 1.0; two comparable lobes -> ~0.5; faint debris ignored."""
    import numpy as np

    from grteclyn_wrapper.metrics.probes.ftl.general import structure_coherence

    n = 64
    yy, xx = np.mgrid[0:n, 0:n]

    def blob(cx: float, cy: float, amp: float, w: float = 4.0):
        return amp * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2.0 * w * w))

    single = blob(32, 32, 1.0)
    assert structure_coherence(single) == 1.0

    two_lobes = blob(20, 32, 1.0) + blob(44, 32, 1.0)
    assert abs(structure_coherence(two_lobes) - 0.5) < 1e-9

    # A dominant lump with only faint debris is still a single coherent structure.
    lump_plus_debris = blob(32, 32, 1.0) + blob(8, 8, 0.05)
    assert structure_coherence(lump_plus_debris) == 1.0


def test_precursor_is_graded_not_binary() -> None:
    """ftl_precursor must rise smoothly with cone-tilt, not snap to 1.0.

    With the old tight saturation scales, max_local_speed 1.05 and 1.30 both
    scored ~1.0.  After widening the scales the stronger tilt must score
    materially higher than the weak one (a usable gradient).
    """
    from grteclyn_wrapper.metrics.probes.ftl.general import GeneralFtlReport
    from grteclyn_wrapper.metrics.types.diagnostics import ConstraintMetrics
    from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics

    def episode(speed: float, frac: float):
        report = GeneralFtlReport(
            f_op=0.0,
            t_min=1.0,
            t_flat=1.0,
            max_local_speed=speed,
            superluminal_fraction=frac,
            path_offaxis=False,
            reachable=True,
            notes=(),
            max_shift=0.5,
        )
        constraints = ConstraintMetrics(
            final_time=2.0,
            max_hamiltonian_l2=1.0e-4,
            max_momentum_l2=1.0e-4,
            final_hamiltonian_l2=1.0e-4,
            final_momentum_l2=1.0e-4,
            min_rho_required=0.0,
            max_rho_required=1.0,
            integral_negative_rho=0.0,
            final_peak_rho_required=1.0,
        )
        return EpisodeMetrics(
            collapse=None,
            constraints=constraints,
            stability=None,
            comoving=None,
            ftl=None,
            termination_reason="",
            general_ftl_evolved=report,
        )

    weak = score_episode(episode(1.05, 0.02), target_stop_time=2.0).components["ftl_precursor"]
    strong = score_episode(episode(1.30, 0.10), target_stop_time=2.0).components["ftl_precursor"]

    assert 0.0 < weak < 0.9, weak
    assert strong > weak + 0.15, (weak, strong)


def test_superluminal_fraction_ignores_sub_margin_noise() -> None:
    """A field only marginally over c=1 must not saturate the fraction.

    The whole evolved slice sits at c~1 with O(1e-3) noise; a bare ``> 1.0``
    test marked nearly every cell superluminal.  The margin must reject a
    uniform 1.01x field (fraction 0) while still counting a genuine 1.05x one.
    """
    import numpy as np

    from grteclyn_wrapper.metrics.probes.ftl.general import operational_ftl_on_grid

    n = 9

    def report_for(eps: float):
        alpha = np.ones((n, n))
        gamma = np.zeros((n, n, 2, 2))
        gamma[:, :, 0, 0] = 1.0
        gamma[:, :, 1, 1] = 1.0
        beta = np.zeros((n, n, 2))
        beta[:, :, 0] = eps  # uniform shift -> max local speed ~ 1 + eps
        mid = n // 2
        return operational_ftl_on_grid(
            alpha, beta, gamma, spacing=(1.0, 1.0), source=(1, mid), target=(n - 2, mid)
        )

    sub_margin = report_for(0.02)  # speed ~1.02, below the 0.05 margin -> noise
    supra_margin = report_for(0.10)  # speed ~1.10, genuinely superluminal

    assert sub_margin.superluminal_fraction == 0.0
    assert supra_margin.superluminal_fraction > 0.5


if __name__ == "__main__":
    test_stability_score_distinguishes_static_from_collapsing()
    test_survival_penalises_dissipated_structure()
    test_ftl_shaping_rewards_scale_with_persistence()
    test_solved_ftl_gated_against_uniform_coordinate_offset()
    test_fragmentation_gates_survival_despite_full_density_retention()
    test_structure_coherence_counts_comparable_lobes()
    test_precursor_is_graded_not_binary()
    test_superluminal_fraction_ignores_sub_margin_noise()
    print("stability score test passed")
