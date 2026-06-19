"""Tests for wave_focusing MAP-Elites descriptors."""

from __future__ import annotations

from grteclyn_wrapper.search.qd_search.descriptors import _bin_index, _descriptor_details


def test_wave_focusing_x_from_central_chromaticity() -> None:
    metrics = {"central": {"wave_chromaticity": 0.72}}
    details = _descriptor_details({}, metrics, mode="wave_focusing")
    assert details["x"] == 0.72


def test_wave_focusing_x_fallback_from_omega_override() -> None:
    overrides = {"grtresna_bs_omega": 0.2}
    details = _descriptor_details({}, None, mode="wave_focusing", overrides=overrides)
    assert details["x"] == 0.5


def test_wave_focusing_y_from_profile_width_override() -> None:
    overrides = {"grtresna_bs_profile_width": 10.0}
    details = _descriptor_details({}, None, mode="wave_focusing", overrides=overrides)
    assert details["y"] == 0.5


def test_descriptor_values_clipped_to_unit_interval() -> None:
    metrics = {"central": {"wave_chromaticity": 1.5}}
    overrides = {"grtresna_bs_profile_width": 100.0}
    details = _descriptor_details({}, metrics, mode="wave_focusing", overrides=overrides)
    assert 0.0 <= details["x"] <= 1.0
    assert 0.0 <= details["y"] <= 1.0


def test_bin_index_maps_extremes() -> None:
    bins = 8
    assert _bin_index(0.0, bins) == 0
    assert _bin_index(0.999, bins) == bins - 1
