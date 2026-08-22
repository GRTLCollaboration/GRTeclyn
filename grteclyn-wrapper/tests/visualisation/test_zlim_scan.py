"""Tests for the series colour-scale scan.

The defect this guards: fields flagged ``per_frame_zlim`` (``chi``,
``chi_minus_1``) rescaled their colourbar from every frame's own percentiles, so
the bar and its tick labels changed in every frame of the movie and the same
colour meant a different number each time.  A scanned series range must win over
that, and over the hand-set presets, without disturbing any other path.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from grteclyn_wrapper.visualisation.process_wave.consume_plotfiles.frames import zlim_scan
from grteclyn_wrapper.visualisation.process_wave.consume_plotfiles.frames.zlim import (
    _SCANNED_KEY,
    _resolve_plot_zlim,
)


def _fake_locker(per_file):
    """Stand in for the yt-backed per-plotfile measurement."""
    def locker(path, args_dict, *, include_per_frame_fields=False):
        return per_file[path]
    return locker


def test_scan_takes_the_envelope_over_the_series(monkeypatch) -> None:
    """The scale must cover every frame, not the first one."""
    per_file = {
        "plt00": {"chi": [0.90, 1.01]},
        "plt01": {"chi": [0.60, 1.02]},
        "plt02": {"chi": [0.20, 1.00]},   # the deep well arrives late
    }
    monkeypatch.setattr(zlim_scan, "_lock_frame_zlims_from_plotfile", _fake_locker(per_file))
    out = zlim_scan.scan_series_zlims(list(per_file), {})
    assert out["chi"] == [0.20, 1.02]
    assert out[_SCANNED_KEY] is True


def test_scan_always_includes_the_last_plotfile(monkeypatch) -> None:
    """Extremes live at late times; dropping the last frame would clip the movie."""
    per_file = {f"plt{i:02d}": {"chi": [1.0 - 0.1 * i, 1.0]} for i in range(5)}
    seen = []

    def locker(path, args_dict, *, include_per_frame_fields=False):
        seen.append(path)
        return per_file[path]

    monkeypatch.setattr(zlim_scan, "_lock_frame_zlims_from_plotfile", locker)
    out = zlim_scan.scan_series_zlims(list(per_file), {}, stride=3)
    assert "plt04" in seen
    assert out["chi"][0] == pytest.approx(0.6)


def test_scan_survives_an_unreadable_plotfile(monkeypatch) -> None:
    """One bad file must not lose the whole scale."""
    def locker(path, args_dict, *, include_per_frame_fields=False):
        if path == "bad":
            raise RuntimeError("corrupt")
        return {"chi": [0.5, 1.0]}

    monkeypatch.setattr(zlim_scan, "_lock_frame_zlims_from_plotfile", locker)
    out = zlim_scan.scan_series_zlims(["bad", "good"], {})
    assert out["chi"] == [0.5, 1.0]


def test_scan_returns_empty_rather_than_a_bogus_scale(monkeypatch) -> None:
    """Nothing readable -> caller must fall back, not render against garbage."""
    def locker(path, args_dict, *, include_per_frame_fields=False):
        raise RuntimeError("no")

    monkeypatch.setattr(zlim_scan, "_lock_frame_zlims_from_plotfile", locker)
    assert zlim_scan.scan_series_zlims(["a", "b"], {}) == {}
    assert zlim_scan.scan_series_zlims([], {}) == {}


def test_scan_records_what_it_measured(monkeypatch, tmp_path) -> None:
    """A movie must be traceable to the scale it was rendered against."""
    monkeypatch.setattr(
        zlim_scan, "_lock_frame_zlims_from_plotfile",
        _fake_locker({"p0": {"chi": [0.3, 1.0]}, "p1": {"chi": [0.2, 1.1]}}),
    )
    zlim_scan.scan_series_zlims(["p0", "p1"], {}, record_dir=str(tmp_path))
    rec = json.loads((tmp_path / zlim_scan.SCAN_RECORD_NAME).read_text())
    assert rec["zlims"]["chi"] == [0.2, 1.1]
    assert rec["plotfiles_scanned"] == 2
    assert _SCANNED_KEY not in rec["zlims"]


def test_scanned_range_beats_per_frame_rescaling() -> None:
    """The actual bug: chi_minus_1 is flagged per_frame_zlim and moved the bar."""
    win = np.linspace(-0.004, 0.004, 400)          # this frame is quiet
    cfg = {"zlim": (-0.9, 0.9), "per_frame_zlim": True}
    scanned = {"chi_minus_1": [-0.05, 0.05], _SCANNED_KEY: True}

    got = _resolve_plot_zlim(
        "chi_minus_1", win, cfg,
        auto_zlim=None, frame_zlims=scanned, use_global_zlim=True,
    )
    assert got == (-0.05, 0.05)


def test_without_a_scan_nothing_changes() -> None:
    """No scan -> the old per-frame behaviour, untouched."""
    win = np.linspace(-0.004, 0.004, 400)
    cfg = {"zlim": (-0.9, 0.9), "per_frame_zlim": True}

    got = _resolve_plot_zlim(
        "chi_minus_1", win, cfg,
        auto_zlim=None, frame_zlims={"chi_minus_1": [-0.05, 0.05]}, use_global_zlim=True,
    )
    assert got != (-0.05, 0.05)          # the t=0 lock must NOT hijack it
    # Per-frame scaling: a symmetric range set by this frame's own 99.5th
    # percentile, so just under the data's own maximum.
    assert got[0] == pytest.approx(-got[1])
    assert 0.0039 < got[1] <= 0.004
