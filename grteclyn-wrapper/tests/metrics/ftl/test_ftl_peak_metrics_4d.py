from __future__ import annotations

from grteclyn_wrapper.search.ftl_peak_metrics import (
    authoritative_geo_strength,
    from_timeseries_mapping,
    peak_fields_for_descriptor_details,
)


def test_from_timeseries_mapping_zeros_frozen_when_4d_ran() -> None:
    ts = {
        "n_frames": 10,
        "f_geo_peak": 0.05,
        "t_at_f_geo_peak": 6.4,
        "ftl_lifetime_fraction": 0.9,
    }
    metrics = {
        "evolving_geodesic": {
            "f_geo": 0.0,
            "h_quality_ok": True,
            "n_rays": 5,
            "n_reached": 5,
        }
    }
    peaks = from_timeseries_mapping(ts, components={"ftl_geo_evolving": 0.0}, metrics=metrics)
    assert peaks is not None
    assert peaks.f_geo_peak == 0.0
    assert peaks.f_geo_evol == 0.0
    assert peaks.ftl_geo_timeavg == 0.0
    assert peaks.ftl_lifetime_fraction == 0.0


def test_peak_fields_report_4d_not_frozen() -> None:
    fields = peak_fields_for_descriptor_details(
        {"n_frames": 10, "f_geo_peak": 0.05, "ftl_lifetime_fraction": 0.9},
        components={"ftl_geo_evolving": 0.0, "ftl_geo_timeavg": 0.04},
        metrics={
            "evolving_geodesic": {
                "f_geo": 0.0,
                "h_quality_ok": True,
                "n_rays": 5,
                "n_reached": 5,
            }
        },
    )
    assert fields["f_geo_peak"] == 0.0
    assert fields["ftl_geo_timeavg"] == 0.0


def test_authoritative_geo_strength_ignores_frozen_component() -> None:
    strength = authoritative_geo_strength(
        {
            "evolving_geodesic": {
                "f_geo": 0.0,
                "h_quality_ok": True,
                "n_rays": 5,
                "n_reached": 5,
            }
        },
        {"operational_ftl_geodesic": 0.5},
    )
    assert strength == 0.0
