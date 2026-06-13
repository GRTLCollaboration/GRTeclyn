from __future__ import annotations

from grteclyn_wrapper.search.ftl_peak_metrics import (
    from_timeseries_mapping,
    from_trajectory_record,
    peak_fields_for_descriptor_details,
)


def test_from_timeseries_mapping_includes_speed_peaks() -> None:
    ts = {
        "n_frames": 2,
        "f_geo_peak": 0.04,
        "t_at_f_geo_peak": 6.4,
        "f_op_peak": 0.05,
        "t_at_f_op_peak": 6.4,
        "max_local_speed_peak": 1.11,
        "t_at_max_speed": 16.0,
        "superluminal_fraction_peak": 0.7,
        "t_at_superluminal_peak": 12.8,
        "ftl_lifetime_fraction": 0.5,
    }
    peaks = from_timeseries_mapping(ts, components={"ftl_geo_timeavg": 0.02})
    assert peaks is not None
    assert peaks.max_local_speed_peak == 1.11
    assert peaks.t_at_max_speed == 16.0
    assert peaks.superluminal_fraction_peak == 0.7
    assert peaks.ftl_geo_timeavg == 0.02


def test_peak_fields_for_descriptor_details() -> None:
    fields = peak_fields_for_descriptor_details(
        {
            "n_frames": 1,
            "f_geo_peak": 0.03,
            "t_at_f_geo_peak": 3.2,
            "f_op_peak": 0.02,
            "t_at_f_op_peak": 3.2,
            "max_local_speed_peak": 1.05,
            "t_at_max_speed": 3.2,
            "superluminal_fraction_peak": 0.1,
            "t_at_superluminal_peak": 3.2,
            "ftl_lifetime_fraction": 1.0,
        },
        components={"ftl_geo_timeavg": 0.015},
    )
    assert fields["f_geo_peak"] == 0.03
    assert fields["max_local_speed_peak"] == 1.05
    assert fields["ftl_geo_timeavg"] == 0.015


def test_from_trajectory_record_enriched_details() -> None:
    record = {
        "eval": 9,
        "status": "gpu_ok",
        "score": 5.0,
        "descriptor_details": {
            "f_geo_peak": 0.04,
            "t_at_f_geo_peak": 6.4,
            "f_op_peak": 0.03,
            "t_at_f_op_peak": 6.4,
            "max_local_speed_peak": 1.08,
            "t_at_max_speed": 16.0,
            "superluminal_fraction_peak": 0.55,
            "t_at_superluminal_peak": 9.6,
            "ftl_lifetime": 0.43,
            "n_frames": 7.0,
            "ftl_geo_timeavg": 0.02,
        },
        "components": {"ftl_geo_timeavg": 0.02},
    }
    peaks = from_trajectory_record(record)
    assert peaks is not None
    assert peaks.f_geo_peak == 0.04
    assert peaks.max_local_speed_peak == 1.08
    assert peaks.ftl_lifetime_fraction == 0.43
