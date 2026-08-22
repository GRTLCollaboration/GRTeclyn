"""Tests for the slice cache — the fix for moving movie colourbars.

Plotfiles are far too large to keep, so they are deleted as they are consumed
and the renderer can never know how large a field will grow later.  Caching the
2-D slice instead — already in memory, ~1e4 times smaller — lets every frame be
redrawn afterwards against one scale measured over the whole run.

The load-bearing property is in ``test_rerendered_colourbar_does_not_move``:
identical colourbar pixels in every frame, while the picture still changes.
"""

from __future__ import annotations

import glob

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from grteclyn_wrapper.visualisation.process_wave.consume_plotfiles.frames import (
    slice_cache,
)

N = 48
EXTENT = [-12.0, 12.0, -12.0, 12.0]


def _well(depth: float) -> np.ndarray:
    x = np.linspace(-12, 12, N)
    xx, yy = np.meshgrid(x, x, indexing="xy")
    return -depth * np.exp(-(xx**2 + yy**2) / 8.0)


def _cache_deepening_series(frames_dir, n=6, first=0.004, last=0.04):
    """A well that deepens over the run: what breaks both existing strategies."""
    depths = np.linspace(first, last, n)
    for i, d in enumerate(depths):
        slice_cache.cache_slice(
            str(frames_dir), "chi_minus_1", "z", i * 10, _well(d), EXTENT,
            time=float(i), coord_val=0.0,
        )
    return depths


def test_cached_slice_round_trips(tmp_path) -> None:
    arr = _well(0.02)
    slice_cache.cache_slice(str(tmp_path), "chi", "z", 7, arr, EXTENT, time=3.5, coord_val=1.25)
    (paths,) = (slice_cache.cached_series(str(tmp_path), "chi", "z"),)
    assert len(paths) == 1
    got, extent, time, coord = slice_cache.load_slice(paths[0])
    assert got == pytest.approx(arr, abs=1e-6)      # float32 storage
    assert extent == EXTENT
    assert (time, coord) == (3.5, 1.25)


def test_series_is_returned_in_frame_order(tmp_path) -> None:
    """Sorting by name breaks at 10 frames; the envelope must not depend on it."""
    for idx in (0, 900, 40, 1000, 200):
        slice_cache.cache_slice(
            str(tmp_path), "chi", "z", idx, _well(0.01), EXTENT, time=0.0, coord_val=0.0
        )
    got = slice_cache.cached_series(str(tmp_path), "chi", "z")
    assert [int(p.rsplit("_", 1)[1].split(".")[0]) for p in got] == [0, 40, 200, 900, 1000]


def test_scale_covers_the_whole_series_not_the_first_frame(tmp_path) -> None:
    """A first-frame lock would clip everything after it."""
    from grteclyn_wrapper.visualisation.process_wave.consume_plotfiles.frames.zlim import (
        _auto_zlim_from_array,
    )

    depths = _cache_deepening_series(tmp_path)
    paths = slice_cache.cached_series(str(tmp_path), "chi_minus_1", "z")
    lo, hi = slice_cache.series_zlim(paths, "chi_minus_1")

    per_frame = [_auto_zlim_from_array(_well(d), "chi_minus_1")[1] for d in depths]
    assert hi == pytest.approx(max(per_frame), rel=1e-6)   # nothing is clipped
    assert hi > 2.0 * per_frame[0]        # and it is NOT the first frame's scale
    assert lo == pytest.approx(-hi)       # symmetric, as this field's rule requires


def test_series_zlim_is_none_when_nothing_readable(tmp_path) -> None:
    assert slice_cache.series_zlim([], "chi") is None
    bad = tmp_path / "slice_0000.npz"
    bad.write_text("not an npz")
    assert slice_cache.series_zlim([str(bad)], "chi") is None


def test_cached_fields_finds_field_names_containing_underscores(tmp_path) -> None:
    """`chi_minus_1_z` must split into (`chi_minus_1`, `z`), not (`chi`, ...)."""
    slice_cache.cache_slice(str(tmp_path), "chi_minus_1", "z", 0, _well(0.01), EXTENT,
                            time=0.0, coord_val=0.0)
    slice_cache.cache_slice(str(tmp_path), "phi", "x", 0, _well(0.01), EXTENT,
                            time=0.0, coord_val=0.0)
    assert set(slice_cache.cached_fields(str(tmp_path))) == {("chi_minus_1", "z"), ("phi", "x")}


def test_rerendered_colourbar_does_not_move(tmp_path) -> None:
    """The bug, end to end: identical bar in every frame, picture still changing."""
    from PIL import Image

    _cache_deepening_series(tmp_path, n=5)
    written, zlim = slice_cache.rerender_series(str(tmp_path), "chi_minus_1", "z")
    assert written == 5 and zlim is not None

    pngs = sorted(glob.glob(str(tmp_path / "chi_minus_1_z" / "frames" / "*.png")))
    assert len(pngs) == 5
    ims = [np.asarray(Image.open(p).convert("L")) for p in pngs]
    assert len({im.shape for im in ims}) == 1, "frame geometry must be stable too"

    w = ims[0].shape[1]
    bar = slice(w - 70, w)
    for im in ims[1:]:
        assert np.array_equal(im[:, bar], ims[0][:, bar]), "colourbar moved"
    body = slice(0, w - 90)
    assert max(np.abs(im[:, body].astype(int) - ims[0][:, body].astype(int)).max()
               for im in ims) > 40, "picture stopped evolving — scale is wrong"


def test_rerender_all_reports_nothing_when_cache_is_empty(tmp_path) -> None:
    assert slice_cache.rerender_all(str(tmp_path)) == {}
