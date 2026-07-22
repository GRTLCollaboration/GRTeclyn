"""Validation tests for freely falling observer timing."""

from __future__ import annotations

import json

import numpy as np

from grteclyn_wrapper.metrics.probes.ftl.metric_field import EvolvingMetricField
from grteclyn_wrapper.metrics.probes.ftl.observer_timing import (
    compute_freefall_observer_timing,
    write_freefall_observer_timing_json,
)


def _flat_field(*, shift_x: float = 0.0) -> EvolvingMetricField:
    n = 9
    metric = np.zeros((n, n, n, 4, 4))
    metric[..., 0, 0] = -1.0 + shift_x * shift_x
    metric[..., 0, 1] = shift_x
    metric[..., 1, 0] = shift_x
    for component in range(1, 4):
        metric[..., component, component] = 1.0
    times = np.linspace(0.0, 10.0, 11)
    return EvolvingMetricField(
        g_stack=np.stack([metric] * len(times)),
        times=times,
        origin=np.array([-2.0, -2.0, -2.0]),
        spacing=(1.0, 0.5, 0.5, 0.5),
    )


def test_flat_freely_falling_observers_have_zero_advance(tmp_path):
    report = compute_freefall_observer_timing(
        _flat_field(),
        emission_tau=2.0,
        observer_step=0.05,
        ray_step=0.05,
    )

    assert report.reached
    assert report.optimizer_success
    assert report.miss_distance < report.reception_tolerance
    assert report.reception_tau is not None
    assert abs(report.reception_tau - report.flat_reception_tau) < 2.0e-3
    assert report.fractional_arrival_advance < 1.0e-4
    assert np.linalg.norm(report.emitter_displacement) < 1.0e-8
    assert np.linalg.norm(report.receiver_displacement) < 1.0e-8

    output = tmp_path / "observer_timing.json"
    write_freefall_observer_timing_json(output, report)
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["reached"] is True
    assert payload["fractional_arrival_advance"] < 1.0e-4


def test_flat_result_is_unchanged_in_constant_shift_coordinates():
    report = compute_freefall_observer_timing(
        _flat_field(shift_x=0.2),
        emission_tau=2.0,
        observer_step=0.05,
        ray_step=0.05,
    )

    assert report.reached
    assert report.reception_tau is not None
    assert abs(report.reception_tau - report.flat_reception_tau) < 2.0e-3
    assert report.fractional_arrival_advance < 1.0e-4
