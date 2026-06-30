"""The matter-confinement diagnostic must be AMR-aware.

Regression guard for the resolution-blind bug: the old extractor sampled a
level-0 covering grid, so confined_frac/peak were identical at max_level 2 and 3
(the refined soliton core was never sampled).  The AMR-aware path integrates over
all cells weighted by true cell volume, so it sees the finest data.
"""

from __future__ import annotations

import numpy as np
import pytest

yt = pytest.importorskip("yt")

from grteclyn_wrapper.visualisation.process_wave.consume_plotfiles.extraction.confinement import (
    _extract_confinement_line,
)


def _parse(row: str | None) -> dict[str, float]:
    assert row is not None
    cols = [float(x) for x in row.split()]
    keys = [
        "time", "total", "peak", "rms_radius", "confined_frac",
        "bary_x", "bary_y", "bary_z", "r_conf",
    ]
    return dict(zip(keys, cols))


def _uniform_ds():
    n = 16
    phi = np.full((n, n, n), 0.05, dtype=np.float64)
    pi = np.zeros((n, n, n), dtype=np.float64)
    data = {("stream", "phi"): phi, ("stream", "Pi"): pi}
    return yt.load_uniform_grid(
        data, (n, n, n),
        bbox=np.array([[0.0, 16.0], [0.0, 16.0], [0.0, 16.0]]),
    )


def _two_level_ds():
    """Base grid with faint matter + a refined grid holding a dense core."""
    n = 16
    base_phi = np.full((n, n, n), 0.01, dtype=np.float64)
    base_pi = np.zeros((n, n, n), dtype=np.float64)
    # Refined grid: level 1, covers the central octant [6,10]^3 (code units),
    # carrying a strong field the base grid cannot represent.
    rdim = 8  # 4 code units * refine_by(2) / (1/level1 dx) -> 8 fine cells
    ref_phi = np.full((rdim, rdim, rdim), 0.5, dtype=np.float64)
    ref_pi = np.zeros((rdim, rdim, rdim), dtype=np.float64)
    grids = [
        {
            "left_edge": [0.0, 0.0, 0.0], "right_edge": [16.0, 16.0, 16.0],
            "level": 0, "dimensions": [n, n, n],
            ("stream", "phi"): base_phi, ("stream", "Pi"): base_pi,
        },
        {
            "left_edge": [6.0, 6.0, 6.0], "right_edge": [10.0, 10.0, 10.0],
            "level": 1, "dimensions": [rdim, rdim, rdim],
            ("stream", "phi"): ref_phi, ("stream", "Pi"): ref_pi,
        },
    ]
    return yt.load_amr_grids(
        grids, [n, n, n],
        bbox=np.array([[0.0, 16.0], [0.0, 16.0], [0.0, 16.0]]),
    )


def test_uniform_grid_amr_path_matches_level0():
    ds = _uniform_ds()
    # Monkeypatch load to return our dataset regardless of path.
    import grteclyn_wrapper.visualisation.process_wave.consume_plotfiles.extraction.confinement as mod
    orig = yt.load
    yt.load = lambda *_a, **_k: ds
    try:
        amr = _parse(_extract_confinement_line("x", t=0.0, well_width=2.0, level=None))
        lvl0 = _parse(_extract_confinement_line("x", t=0.0, well_width=2.0, level=0))
    finally:
        yt.load = orig
    # Single level => identical (dV constant cancels in the fraction).
    assert amr["confined_frac"] == pytest.approx(lvl0["confined_frac"], rel=1e-9)
    assert amr["peak"] == pytest.approx(lvl0["peak"], rel=1e-9)


def test_amr_path_sees_refined_core_that_level0_misses():
    ds = _two_level_ds()
    orig = yt.load
    yt.load = lambda *_a, **_k: ds
    try:
        amr = _parse(_extract_confinement_line("x", t=0.0, well_width=2.0, level=None))
        lvl0 = _parse(_extract_confinement_line("x", t=0.0, well_width=2.0, level=0))
    finally:
        yt.load = orig
    # The dense core (phi=0.5) lives only on the refined grid: the AMR-aware peak
    # must pick it up, while the level-0 covering grid sees only the base 0.01.
    assert amr["peak"] > 10.0 * lvl0["peak"]
    assert amr["peak"] == pytest.approx(0.5, rel=1e-6)
