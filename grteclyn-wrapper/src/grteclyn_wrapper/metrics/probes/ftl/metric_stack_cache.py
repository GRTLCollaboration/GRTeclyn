"""Persist per-plotfile 4-metric slices before HDF5 plotfiles are deleted."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from .geodesic import build_metric_3d_from_plotfile
from .metric_field import EvolvingMetricField

_PLOT_INDEX_RE = re.compile(r"(\d+)\s*$")
DEFAULT_N_SPACE = 65
METRIC_STACK_SUBDIR = "metric_stack"


def metric_stack_dir(small_data_dir: Path) -> Path:
    return Path(small_data_dir) / METRIC_STACK_SUBDIR


def _plot_sort_key(path: Path) -> int:
    match = _PLOT_INDEX_RE.search(path.stem)
    return int(match.group(1)) if match else 0


def append_slice_from_plotfile(
    plotfile: str | Path,
    cache_dir: Path,
    *,
    t: float,
    n_space: int = DEFAULT_N_SPACE,
    half_width: float | None = None,
) -> Path:
    """Sample ``g_{mu nu}`` from a plotfile and write one compressed slice file."""
    plotfile = Path(plotfile)
    g, origin, spacing = build_metric_3d_from_plotfile(
        plotfile, n=n_space, half_width=half_width
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = cache_dir / f"{plotfile.name}.npz"
    np.savez_compressed(
        out,
        t=np.float64(t),
        g=g.astype(np.float32),
        origin=origin.astype(np.float64),
        spacing=np.asarray(spacing, dtype=np.float64),
    )
    return out


def list_slice_files(cache_dir: Path) -> list[Path]:
    if not cache_dir.is_dir():
        return []
    return sorted(cache_dir.glob("*.npz"), key=_plot_sort_key)


def slice_count(cache_dir: Path) -> int:
    return len(list_slice_files(cache_dir))


def subsample_slice_files(
    files: Sequence[Path],
    *,
    stride: int = 1,
    max_slices: int | None = None,
) -> list[Path]:
    """Reduce temporal resolution for fast in-loop 4D scoring."""
    if not files:
        return []
    stride = max(1, int(stride))
    picked = list(files[::stride]) if stride > 1 else list(files)
    if max_slices is not None and len(picked) > max_slices:
        idx = np.linspace(0, len(picked) - 1, int(max_slices), dtype=int)
        picked = [picked[int(i)] for i in idx]
    return picked


def evolving_field_from_metric_stack_cache(
    cache_dir: Path,
    *,
    slice_stride: int = 1,
    max_slices: int | None = None,
) -> EvolvingMetricField | None:
    """Rebuild ``EvolvingMetricField`` from cached per-plotfile slices."""
    files = subsample_slice_files(
        list_slice_files(cache_dir),
        stride=slice_stride,
        max_slices=max_slices,
    )
    if len(files) < 3:
        return None

    times: list[float] = []
    slices: list[NDArray[np.float64]] = []
    origin: NDArray[np.float64] | None = None
    spacing_xyz: tuple[float, float, float] | None = None

    for path in files:
        data = np.load(path)
        times.append(float(data["t"]))
        slices.append(np.asarray(data["g"], dtype=np.float64))
        if origin is None:
            origin = np.asarray(data["origin"], dtype=np.float64)
            sp = np.asarray(data["spacing"], dtype=np.float64)
            spacing_xyz = (float(sp[0]), float(sp[1]), float(sp[2]))

    assert origin is not None and spacing_xyz is not None
    times_arr = np.asarray(times, dtype=float)
    dt = float(np.mean(np.diff(times_arr))) if len(times_arr) > 1 else 1.0
    if not np.isfinite(dt) or dt <= 0.0:
        dt = 1.0
    g_stack = np.stack(slices, axis=0)
    return EvolvingMetricField(
        g_stack=g_stack,
        times=times_arr,
        origin=origin,
        spacing=(dt, spacing_xyz[0], spacing_xyz[1], spacing_xyz[2]),
    )
