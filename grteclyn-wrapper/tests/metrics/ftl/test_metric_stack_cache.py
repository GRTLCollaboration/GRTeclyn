"""Metric stack cache: persist slices before plotfile deletion."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from grteclyn_wrapper.metrics.probes import warpfactory as wf

from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import (
    compute_evolving_geodesic_ftl,
    compute_evolving_geodesic_ftl_from_metric_stack_cache,
)
from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic_options import (
    HQ_OPTIONS,
    SEARCH_OPTIONS,
)
from grteclyn_wrapper.metrics.probes.ftl.metric_field import evolving_field_from_analytic_stack
from grteclyn_wrapper.metrics.probes.ftl.metric_stack_cache import (
    evolving_field_from_metric_stack_cache,
    metric_stack_dir,
    slice_count,
)


def _write_alcubierre_slices(cache_dir: Path) -> int:
    g, spacing = wf.alcubierre_metric(
        n_space=33,
        velocity=0.5,
        half_width=4.0,
        dt=0.2,
    )
    field = evolving_field_from_analytic_stack(g, spacing)
    cache_dir.mkdir(parents=True, exist_ok=True)
    for i, t in enumerate(field.times):
        out = cache_dir / f"RadialRecipePlt{i * 24:05d}.npz"
        np.savez_compressed(
            out,
            t=np.float64(t),
            g=field.g_stack[i].astype(np.float32),
            origin=field.origin.astype(np.float64),
            spacing=np.asarray(field.spatial_spacing, dtype=np.float64),
        )
    return field.g_stack.shape[0]


def test_metric_stack_cache_round_trip(tmp_path: Path) -> None:
    cache_dir = metric_stack_dir(tmp_path)
    n_time = _write_alcubierre_slices(cache_dir)
    assert slice_count(cache_dir) == n_time
    field = evolving_field_from_metric_stack_cache(cache_dir)
    assert field is not None
    assert field.g_stack.shape[0] == n_time
    assert field.g_stack.shape[1:4] == (33, 33, 33)


def test_evolving_trace_from_metric_stack_cache(tmp_path: Path) -> None:
    cache_dir = metric_stack_dir(tmp_path)
    _write_alcubierre_slices(cache_dir)
    report = compute_evolving_geodesic_ftl_from_metric_stack_cache(
        cache_dir,
        options=HQ_OPTIONS,
    )
    assert report is not None
    assert report.h_quality_ok
    assert report.f_geo > 0.1
    assert report.n_reached == report.n_rays
    assert report.f_geo_frozen_peak is not None
    assert report.f_geo_frozen_peak > 0.0


def test_cached_field_matches_direct_analytic_stack(tmp_path: Path) -> None:
    g, spacing = wf.alcubierre_metric(n_space=33, velocity=0.5, half_width=4.0, dt=0.2)
    direct = evolving_field_from_analytic_stack(g, spacing)
    cache_dir = metric_stack_dir(tmp_path)
    _write_alcubierre_slices(cache_dir)
    loaded = evolving_field_from_metric_stack_cache(cache_dir)
    assert loaded is not None
    direct_report = compute_evolving_geodesic_ftl(direct)
    cached_report = compute_evolving_geodesic_ftl(loaded)
    assert abs(direct_report.f_geo - cached_report.f_geo) < 0.05
