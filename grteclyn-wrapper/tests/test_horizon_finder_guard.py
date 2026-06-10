"""Regression guard for the miscentered apparent-horizon finder.

The C++ theta_plus / apparent-horizon proxy must be measured about the physics
grid center.  When it was measured about the domain corner instead, the
minimum-expansion point landed at ``r ~ |grid_center|`` (near the domain
half-width) and produced a spurious trapped-surface verdict.  These tests pin
the parsing of ``r_at_min_theta_plus`` and the score-side guard that rejects an
off-center "horizon" as an artifact instead of silently vetoing the candidate.
"""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from grteclyn_wrapper.metrics import read_episode_metrics, score_episode


def _write_collapse(root: Path, *, ah_radius: float, min_theta: float, r_at_min_theta: float) -> None:
    """Write a collapse_diagnostics.dat with a trapped signal located at
    ``r_at_min_theta`` (column 9, about the grid center)."""
    data = root / "data"
    data.mkdir(parents=True)
    with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as handle:
        # cols: t lapse chi max_k c4 c5 c6 ah_r theta r_at_min_theta phi.. (14 total)
        for t in (0.0, 1.0, 2.0):
            handle.write(
                f"{t:g} 1.0 0.98 0.02 0 0 0 "
                f"{ah_radius:g} {min_theta:g} {r_at_min_theta:g} 0 0 0 0\n"
            )
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as handle:
        for t in (0.0, 1.0, 2.0):
            handle.write(f"{t:g} 1e-4 1e-4 0 0 0\n")


def test_r_at_min_theta_plus_is_parsed() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_collapse(root, ah_radius=110.0, min_theta=-0.03, r_at_min_theta=108.0)
        metrics = read_episode_metrics(root)
        assert metrics.collapse is not None
        assert metrics.collapse.r_at_min_theta_plus == 108.0
        assert metrics.collapse.max_horizon_radius == 110.0


def test_offcenter_horizon_is_rejected_when_domain_known() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        # L_full=128 -> half-width 64; r=108 > 0.5*64 -> miscentered artifact.
        _write_collapse(root, ah_radius=110.0, min_theta=-0.03, r_at_min_theta=108.0)
        metrics = read_episode_metrics(root)
        score = score_episode(metrics, target_stop_time=2.0, domain_half_width=64.0)
        assert score.components["horizon_penalty"] == 0.0
        assert any("miscentered" in note for note in score.notes)


def test_offcenter_horizon_penalizes_without_domain_info() -> None:
    """Backward compatible: with no domain half-width the guard is inactive."""
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_collapse(root, ah_radius=110.0, min_theta=-0.03, r_at_min_theta=108.0)
        metrics = read_episode_metrics(root)
        score = score_episode(metrics, target_stop_time=2.0)
        assert score.components["horizon_penalty"] == -1.0


def test_interior_trapped_surface_is_penalized() -> None:
    """A genuine interior horizon (small r) is still penalized."""
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_collapse(root, ah_radius=4.0, min_theta=-0.03, r_at_min_theta=4.0)
        metrics = read_episode_metrics(root)
        score = score_episode(metrics, target_stop_time=2.0, domain_half_width=64.0)
        assert score.components["horizon_penalty"] == -1.0


if __name__ == "__main__":
    test_r_at_min_theta_plus_is_parsed()
    test_offcenter_horizon_is_rejected_when_domain_known()
    test_offcenter_horizon_penalizes_without_domain_info()
    test_interior_trapped_surface_is_penalized()
    print("horizon finder guard tests passed")
