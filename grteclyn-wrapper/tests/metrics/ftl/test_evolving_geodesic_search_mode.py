"""Fast search-profile 4D evolving geodesic integration."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from grteclyn_wrapper.metrics.probes import warpfactory as wf
from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import (
    compute_evolving_geodesic_ftl,
    compute_evolving_geodesic_ftl_from_metric_stack_cache,
)
from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic_options import SEARCH_OPTIONS
from grteclyn_wrapper.metrics.probes.ftl.metric_field import evolving_field_from_analytic_stack
from grteclyn_wrapper.metrics.probes.ftl.metric_stack_cache import (
    metric_stack_dir,
    subsample_slice_files,
)


def _write_alcubierre_stack(cache_dir: Path, *, n_slices: int = 8) -> None:
    g, spacing = wf.alcubierre_metric(
        n_space=33,
        velocity=0.5,
        half_width=4.0,
        dt=0.25,
    )
    field = evolving_field_from_analytic_stack(g, spacing)
    cache_dir.mkdir(parents=True, exist_ok=True)
    for i in range(min(n_slices, field.g_stack.shape[0])):
        out = cache_dir / f"RadialRecipePlt{i * 24:05d}.npz"
        np.savez_compressed(
            out,
            t=np.float64(field.times[i]),
            g=field.g_stack[i].astype(np.float32),
            origin=field.origin.astype(np.float64),
            spacing=np.asarray(field.spatial_spacing, dtype=np.float64),
        )


def test_subsample_slice_files_evenly_caps_count() -> None:
    files = [Path(f"RadialRecipePlt{i:05d}.npz") for i in range(0, 240, 24)]
    picked = subsample_slice_files(files, stride=2, max_slices=4)
    assert len(picked) == 4
    assert picked[0] == files[0]
    assert picked[-1] == files[8]


def test_search_mode_skips_frozen_peak_scan(tmp_path: Path) -> None:
    cache = metric_stack_dir(tmp_path / "small_data")
    _write_alcubierre_stack(cache, n_slices=8)
    with patch(
        "grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic._frozen_peak_from_g_slices",
    ) as frozen_mock:
        report = compute_evolving_geodesic_ftl_from_metric_stack_cache(
            cache,
            options=SEARCH_OPTIONS,
        )
    frozen_mock.assert_not_called()
    assert report is not None
    assert report.f_geo_frozen_peak is None


def test_search_mode_uses_fewer_rays_than_default(tmp_path: Path) -> None:
    cache = metric_stack_dir(tmp_path / "small_data")
    _write_alcubierre_stack(cache, n_slices=8)
    with patch(
        "grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic.compute_evolving_geodesic_ftl",
    ) as compute_mock:
        compute_mock.return_value = None
        compute_evolving_geodesic_ftl_from_metric_stack_cache(
            cache,
            options=SEARCH_OPTIONS,
        )
    assert compute_mock.call_args.kwargs["n_rays"] == SEARCH_OPTIONS.n_rays
