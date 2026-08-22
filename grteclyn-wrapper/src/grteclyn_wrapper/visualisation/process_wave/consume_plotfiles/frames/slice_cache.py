"""Keep the 2-D slices, not the plotfiles, so a movie can be rescaled later.

THE PROBLEM.  Frames are rendered as plotfiles arrive and the plotfiles are
deleted immediately -- they are far too large to keep.  So the renderer can
never know how large a field will grow later in the run, and both of its
choices are wrong for a paper movie: rescaling every frame makes the colourbar
and its tick labels move in every frame, and locking the scale from the first
plotfile puts late-time features off the end of it.

THE FIX.  The thing the renderer actually draws is a single 2-D slice, and it is
already in memory when the frame is drawn.  A slice is around ten thousand times
smaller than the plotfile it came from -- at N = 512 that is a ~1 MB array out of
a ~30 GB plotfile -- so it can be kept for the whole run at negligible cost.
Afterwards the envelope over every cached slice gives the one scale that clips
nothing, and every frame is redrawn against it.

Rough cost: ``N * N * 4`` bytes per field per frame before compression, and these
fields are smooth, so ``savez_compressed`` typically takes a large bite out of
that.  For N = 512, 19 fields and 500 frames the raw figure is about 5 GB;
cache only the fields you actually want movies of if that matters.

Only the native (uniform covering grid) render path caches.  AMR levels go
through yt's fixed-resolution buffer instead and are not cached yet.
"""

from __future__ import annotations

import json
import os
import re
from typing import Iterable, Sequence

import numpy as np

#: Slices live beside the frames, under a name that cannot collide with a
#: ``<field>_<axis>`` frame directory.
CACHE_DIR_NAME = "_slice_cache"

_SLICE_RE = re.compile(r"slice_(\d+)\.npz$")


def _series_dir(frames_out_dir: str, field: str, axis: str) -> str:
    return os.path.join(frames_out_dir, CACHE_DIR_NAME, f"{field}_{axis}")


def cache_slice(
    frames_out_dir: str,
    field: str,
    axis: str,
    frame_idx: int,
    arr: np.ndarray,
    extent: Sequence[float],
    *,
    time: float,
    coord_val: float,
) -> str:
    """Store one slice, with everything needed to redraw its frame."""
    out_dir = _series_dir(frames_out_dir, field, axis)
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"slice_{frame_idx:04d}.npz")
    np.savez_compressed(
        path,
        arr=np.asarray(arr, dtype=np.float32),
        extent=np.asarray(extent, dtype=np.float64),
        time=np.float64(time),
        coord_val=np.float64(coord_val),
    )
    return path


def cached_series(frames_out_dir: str, field: str, axis: str) -> list[str]:
    """Cached slices for one field/axis, in frame order."""
    out_dir = _series_dir(frames_out_dir, field, axis)
    if not os.path.isdir(out_dir):
        return []
    found = [
        os.path.join(out_dir, n) for n in os.listdir(out_dir) if _SLICE_RE.search(n)
    ]
    return sorted(found, key=lambda p: int(_SLICE_RE.search(p).group(1)))


def cached_fields(frames_out_dir: str) -> list[tuple[str, str]]:
    """Every ``(field, axis)`` pair with cached slices."""
    root = os.path.join(frames_out_dir, CACHE_DIR_NAME)
    if not os.path.isdir(root):
        return []
    pairs = []
    for name in sorted(os.listdir(root)):
        if not os.path.isdir(os.path.join(root, name)) or "_" not in name:
            continue
        field, _, axis = name.rpartition("_")
        if field and axis:
            pairs.append((field, axis))
    return pairs


def load_slice(path: str) -> tuple[np.ndarray, list[float], float, float]:
    with np.load(path) as data:
        return (
            np.asarray(data["arr"], dtype=np.float64),
            [float(v) for v in data["extent"]],
            float(data["time"]),
            float(data["coord_val"]),
        )


def series_zlim(paths: Iterable[str], field: str) -> tuple[float, float] | None:
    """The scale that clips nothing across a whole cached series.

    Each slice is asked for the limits it would have chosen on its own -- the
    same rule the live renderer uses -- and the envelope of those is taken.  So
    no frame is scaled worse than it would have been per-frame, and the scale is
    the same in all of them.
    """
    from .zlim import _auto_zlim_from_array

    lo = hi = None
    for path in paths:
        try:
            arr, _, _, _ = load_slice(path)
        except (OSError, ValueError):
            continue
        one = _auto_zlim_from_array(arr, field)
        if one is None:
            continue
        lo = one[0] if lo is None else min(lo, one[0])
        hi = one[1] if hi is None else max(hi, one[1])
    if lo is None or hi is None:
        return None
    return (float(lo), float(hi))


def rerender_series(
    frames_out_dir: str,
    field: str,
    axis: str,
    *,
    zlim: Sequence[float] | None = None,
    corner: bool = False,
    verbose: bool = False,
) -> tuple[int, tuple[float, float] | None]:
    """Redraw every frame of one series against a single fixed scale.

    Returns ``(frames_written, zlim_used)``.  With ``zlim`` given, that scale is
    used; otherwise it is measured from the cache.
    """
    from ..config import _field_frame_config
    from .slice import draw_slice_png

    paths = cached_series(frames_out_dir, field, axis)
    if not paths:
        return (0, None)

    limits = tuple(zlim) if zlim is not None else series_zlim(paths, field)
    if limits is None:
        return (0, None)

    cfg = _field_frame_config(field)
    written = 0
    for path in paths:
        idx = int(_SLICE_RE.search(path).group(1))
        try:
            arr, extent, time, coord_val = load_slice(path)
        except (OSError, ValueError) as exc:
            print(f"WARNING: rerender skipped {os.path.basename(path)}: {exc}")
            continue
        draw_slice_png(
            arr, extent,
            field=field, cfg=cfg, axis=axis,
            coord_val=coord_val, time=time, zlim=limits,
            frames_out_dir=frames_out_dir, frame_idx=idx,
            corner=corner, verbose=verbose, note=" (cached)",
        )
        written += 1
    return (written, (float(limits[0]), float(limits[1])))


def rerender_all(
    frames_out_dir: str, *, corner: bool = False, verbose: bool = False
) -> dict[str, list[float]]:
    """Redraw every cached series, each against its own fixed scale."""
    used: dict[str, list[float]] = {}
    for field, axis in cached_fields(frames_out_dir):
        written, limits = rerender_series(
            frames_out_dir, field, axis, corner=corner, verbose=verbose
        )
        if not written or limits is None:
            print(f"[rerender] {field}_{axis}: nothing cached, skipped")
            continue
        used[f"{field}_{axis}"] = [limits[0], limits[1]]
        print(
            f"[rerender] {field}_{axis}: {written} frame(s) at a fixed "
            f"{limits[0]:.6g} .. {limits[1]:.6g}"
        )
    if used:
        record = os.path.join(frames_out_dir, CACHE_DIR_NAME, "rerender_zlims.json")
        try:
            with open(record, "w", encoding="utf-8") as fh:
                json.dump(used, fh, indent=2, sort_keys=True)
        except OSError as exc:
            print(f"WARNING: could not write {record}: {exc}")
    return used
