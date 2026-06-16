from __future__ import annotations

from grteclyn_wrapper.search.qd_search.descriptors import _descriptor_details


def test_ftl_lifetime_ignores_frozen_peak_when_4d_ran_zero() -> None:
    components = {
        "ftl_geo_evolving": 0.0,
        "operational_ftl_geodesic": 0.5,
    }
    metrics = {
        "ftl_timeseries": {
            "n_frames": 10,
            "f_geo_peak": 0.05,
            "ftl_lifetime_fraction": 0.9,
        },
        "evolving_geodesic": {
            "f_geo": 0.0,
            "h_quality_ok": True,
            "n_rays": 5,
            "n_reached": 5,
        },
    }
    details = _descriptor_details(components, metrics, mode="ftl_lifetime")
    assert details["x"] == 0.0
    assert details["y"] == 0.0
    assert details["f_geo_peak"] == 0.0
    assert details["ftl_geo_timeavg"] == 0.0


def test_ftl_lifetime_uses_4d_strength_when_shortcut_found() -> None:
    components = {"ftl_geo_evolving": 0.12}
    metrics = {
        "ftl_timeseries": {"n_frames": 10, "f_geo_peak": 0.05, "ftl_lifetime_fraction": 0.1},
        "evolving_geodesic": {
            "f_geo": 0.12,
            "h_quality_ok": True,
            "n_rays": 3,
            "n_reached": 3,
        },
    }
    details = _descriptor_details(components, metrics, mode="ftl_lifetime")
    assert details["x"] > 0.0
    assert details["y"] == 1.0
    assert details["f_geo_evol"] == 0.12
