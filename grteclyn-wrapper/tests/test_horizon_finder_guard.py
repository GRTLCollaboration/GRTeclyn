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


def _write_collapse(
    root: Path,
    *,
    rows: list[tuple[float, float, float, float, float, float]],
) -> None:
    """Write collapse rows: (t, lapse, ah_r, min_theta, r_at_min_theta)."""
    data = root / "data"
    data.mkdir(parents=True)
    with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as handle:
        for t, lapse, ah_radius, min_theta, r_at_min_theta in rows:
            handle.write(
                f"{t:g} {lapse:g} 0.98 0.02 0 0 0 "
                f"{ah_radius:g} {min_theta:g} {r_at_min_theta:g} 0 0 0 0\n"
            )
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as handle:
        for t, *_rest in rows:
            handle.write(f"{t:g} 1e-4 1e-4 0 0 0\n")


def test_r_at_min_theta_plus_is_parsed() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_collapse(
            root,
            rows=[
                (0.0, 1.0, 110.0, -0.03, 108.0),
                (1.0, 1.0, 110.0, -0.03, 108.0),
                (2.0, 1.0, 110.0, -0.03, 108.0),
            ],
        )
        metrics = read_episode_metrics(root)
        assert metrics.collapse is not None
        assert metrics.collapse.r_at_min_theta_plus == 108.0
        assert metrics.collapse.max_horizon_radius == 110.0


def test_offcenter_horizon_is_rejected_when_domain_known() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_collapse(
            root,
            rows=[
                (0.0, 0.1, 110.0, -0.5, 108.0),
                (1.0, 0.1, 110.0, -0.5, 108.0),
                (2.0, 0.1, 110.0, -0.5, 108.0),
            ],
        )
        metrics = read_episode_metrics(root)
        score = score_episode(metrics, target_stop_time=2.0, domain_half_width=64.0)
        assert score.components["horizon_penalty"] == 0.0
        assert any("miscentered" in note for note in score.notes)


def test_uncorroborated_theta_is_suppressed() -> None:
    """theta+<=0 with a healthy lapse is not a corroborated trapped surface."""
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_collapse(
            root,
            rows=[
                (0.0, 1.0, 7.0, -0.01, 7.0),
                (1.0, 1.0, 8.0, -0.02, 7.5),
                (2.0, 1.0, 9.0, -0.03, 8.0),
            ],
        )
        metrics = read_episode_metrics(root)
        assert metrics.collapse is not None
        assert metrics.collapse.corroborated_trapped is False
        score = score_episode(metrics, target_stop_time=2.0, domain_half_width=32.0)
        assert score.components["horizon_penalty"] == 0.0
        assert any("uncorroborated" in note.lower() or "without lapse-collapse" in note for note in score.notes)


def test_offcenter_horizon_penalizes_without_domain_info() -> None:
    """Backward compatible: with no domain half-width the guard is inactive."""
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_collapse(
            root,
            rows=[
                (0.0, 0.1, 110.0, -0.5, 108.0),
                (1.0, 0.1, 110.0, -0.5, 108.0),
                (2.0, 0.1, 110.0, -0.5, 108.0),
            ],
        )
        metrics = read_episode_metrics(root)
        score = score_episode(metrics, target_stop_time=2.0)
        assert score.components["horizon_penalty"] == -1.0


def test_interior_trapped_surface_is_penalized() -> None:
    """A genuine interior horizon (small r, lapse collapse) is still penalized."""
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_collapse(
            root,
            rows=[
                (0.0, 0.1, 4.0, -0.5, 4.0),
                (1.0, 0.08, 4.0, -0.6, 4.0),
                (2.0, 0.05, 4.0, -0.7, 4.0),
            ],
        )
        metrics = read_episode_metrics(root)
        assert metrics.collapse is not None
        assert metrics.collapse.corroborated_trapped is True
        score = score_episode(metrics, target_stop_time=2.0, domain_half_width=64.0)
        assert score.components["horizon_penalty"] == -1.0


def test_late_only_corroborated_collapse_is_suppressed() -> None:
    """Late trailing collapse must not veto an otherwise healthy FTL window."""
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        rows = []
        for step in range(16):
            t = float(step)
            if t < 12.0:
                rows.append((t, 0.95, 8.0, -0.01, 7.0))
            else:
                rows.append((t, 0.12, 10.0, -1.5, 5.0))
        _write_collapse(root, rows=rows)
        metrics = read_episode_metrics(root)
        assert metrics.collapse is not None
        assert metrics.collapse.corroborated_trapped is True
        assert metrics.collapse.first_corroborated_time == 12.0
        score = score_episode(metrics, target_stop_time=15.0, domain_half_width=32.0)
        assert score.components["horizon_penalty"] == 0.0
        assert any("late collapse penalty suppressed" in note for note in score.notes)


if __name__ == "__main__":
    test_r_at_min_theta_plus_is_parsed()
    test_offcenter_horizon_is_rejected_when_domain_known()
    test_uncorroborated_theta_is_suppressed()
    test_offcenter_horizon_penalizes_without_domain_info()
    test_interior_trapped_surface_is_penalized()
    test_late_only_corroborated_collapse_is_suppressed()
    print("horizon finder guard tests passed")
