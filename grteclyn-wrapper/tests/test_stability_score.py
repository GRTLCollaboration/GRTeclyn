from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from grteclyn_wrapper.metrics import read_episode_metrics
from grteclyn_wrapper.score import score_episode


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


if __name__ == "__main__":
    test_stability_score_distinguishes_static_from_collapsing()
    print("stability score test passed")
