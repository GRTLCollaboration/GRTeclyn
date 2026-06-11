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


if __name__ == "__main__":
    test_stability_score_distinguishes_static_from_collapsing()
    test_survival_penalises_dissipated_structure()
    test_ftl_shaping_rewards_scale_with_persistence()
    print("stability score test passed")
