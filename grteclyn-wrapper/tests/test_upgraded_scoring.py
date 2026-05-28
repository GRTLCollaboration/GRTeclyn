from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory

from grteclyn_wrapper.constrained_recipe import constrained_overrides
from grteclyn_wrapper.ftl_metrics import compute_ftl_metrics
from grteclyn_wrapper.metrics import read_episode_metrics
from grteclyn_wrapper.score import score_episode
from grteclyn_wrapper.seeds import get_seed


def _write_warp_episode(root: Path) -> None:
    seed = get_seed("alcubierre_warp", velocity=0.5)
    overrides = dict(seed.overrides)
    constrained_overrides(overrides, phantom=True)
    root.mkdir(parents=True)
    (root / "metadata.json").write_text(
        json.dumps({"overrides": overrides}),
        encoding="utf-8",
    )
    data = root / "data"
    small = root / "small_data"
    data.mkdir()
    small.mkdir()
    with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as handle:
        for t in (0.0, 1.0, 2.0):
            chi = 0.9 - 0.02 * t
            handle.write(f"{t:g} 0.95 {chi:g} 0.08 0 0 0 0 0 0 0 0 0 0\n")
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as handle:
        for t in (0.0, 1.0, 2.0):
            handle.write(f"{t:g} 1e-4 1e-4 0 0 0\n")
    with (small / "shell_profiles.dat").open("w", encoding="utf-8") as handle:
        handle.write("# time chi_mean_4\n")
        for t in (0.0, 1.0, 2.0):
            handle.write(f"{t:g} {0.9 - 0.01 * t:g}\n")
    with (small / "areal_radius.dat").open("w", encoding="utf-8") as handle:
        for t in (0.0, 1.0, 2.0):
            handle.write(f"{t:g} 1.0 0.0625\n")


def test_alcubierre_beats_flat_with_upgraded_scoring() -> None:
    flat = get_seed("flat_minkowski")
    warp = get_seed("alcubierre_warp", velocity=0.5)
    flat_m = compute_ftl_metrics(flat.overrides, L=8.0)
    warp_m = compute_ftl_metrics(warp.overrides, L=8.0)

    assert flat_m.f_log == 0.0
    assert flat_m.f_asymmetry == 0.0
    assert warp_m.f_log > 0.45
    assert warp_m.f_asymmetry > 0.05
    assert warp_m.f_log > flat_m.f_log

    with TemporaryDirectory() as tmp:
        episode = Path(tmp) / "alcubierre"
        _write_warp_episode(episode)
        metrics = read_episode_metrics(episode, ftl_L=8.0)
        score = score_episode(metrics, target_stop_time=2.0)
        assert score.components["ftl_shortcut"] > 0.45
        assert score.components["expansion_asymmetry"] > 0.05
        assert score.components["comoving_stability"] > 0.5
        assert score.total > 9.5


def test_ellis_has_high_asymmetry_without_null_shortcut() -> None:
    seed = get_seed("ellis_bronnikov", num_bases=16)
    metrics = compute_ftl_metrics(seed.overrides, L=8.0)
    assert metrics.f_null == 0.0
    assert metrics.f_asymmetry > 0.95
    assert metrics.f_log > 0.95


if __name__ == "__main__":
    test_alcubierre_beats_flat_with_upgraded_scoring()
    test_ellis_has_high_asymmetry_without_null_shortcut()
    print("upgraded scoring tests passed")
