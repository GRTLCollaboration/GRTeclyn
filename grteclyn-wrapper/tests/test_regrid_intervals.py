from grteclyn_wrapper.core.params import regrid_intervals_for_max_level


def test_regrid_intervals_for_max_level_matches_amr_levels() -> None:
    assert regrid_intervals_for_max_level(0) == []
    assert regrid_intervals_for_max_level(1) == [16]
    assert regrid_intervals_for_max_level(2) == [16, 16]
    assert regrid_intervals_for_max_level(3) == [16, 16, 8]


def test_promotion_overrides_replaces_stale_cmaes_intervals() -> None:
    import importlib.util
    from pathlib import Path

    script = Path(__file__).resolve().parents[1] / "scripts/search/replay_grtresna_eval.py"
    spec = importlib.util.spec_from_file_location("replay_grtresna_eval", script)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    out = mod._promotion_overrides(
        {"max_level": 2, "regrid_interval": [16, 16]},
        n_full=256,
        l_full=128.0,
        stop_time=30.0,
        plot_interval=24,
        max_level=3,
        regrid_threshold=0.02,
    )
    assert out["max_level"] == 3
    assert out["regrid_interval"] == [16, 16, 8]
