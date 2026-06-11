"""Tests for solved-geometry operational FTL on .gridinit files."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

from grteclyn_wrapper.grtresna.io import GridinitData, read_gridinit, write_gridinit
from grteclyn_wrapper.metrics.probes.ftl.general import operational_ftl_on_grid
from grteclyn_wrapper.metrics.probes.ftl.solved import (
    build_xz_slice_from_gridinit,
    compute_solved_geometry_ftl,
)
from grteclyn_wrapper.search.solved_ftl_gate import (
    solved_ftl_has_signal,
    solved_geometry_rejection_fitness,
)


def test_gridinit_roundtrip():
    nz, ny, nx = 8, 8, 16
    n_comp = 4
    names = ["chi", "lapse", "shift1", "shift3"]
    data = np.ones((nz, ny, nx, n_comp), dtype=np.float64)
    data[:, :, :, 0] = 4.0  # chi > 1 in channel for portal shortcut
    dx = np.array([1.0, 1.0, 1.0])
    origin = np.array([0.0, 0.0, 0.0])
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "test.gridinit"
        write_gridinit(data, names, dx, origin, path)
        loaded = read_gridinit(path)
    assert loaded.data.shape == data.shape
    assert loaded.comp_names == names
    np.testing.assert_allclose(loaded.dx_xyz, dx)
    np.testing.assert_allclose(loaded.data, data)


def test_flat_gridinit_no_ftl_signal():
    nz, ny, nx = 12, 12, 24
    names = ["chi", "h11", "h33", "lapse", "shift1", "shift3"]
    n_comp = len(names)
    data = np.zeros((nz, ny, nx, n_comp), dtype=np.float64)
    for i, name in enumerate(names):
        if name == "chi":
            data[:, :, :, i] = 1.0
        elif name.startswith("h"):
            data[:, :, :, i] = 1.0
        elif name == "lapse":
            data[:, :, :, i] = 1.0
    dx = np.array([0.5, 0.5, 0.5])
    origin = np.array([-6.0, -3.0, -3.0])
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "flat.gridinit"
        write_gridinit(data, names, dx, origin, path)
        solved = compute_solved_geometry_ftl(path, n=41, L=6.0)
    assert solved is not None
    assert not solved_ftl_has_signal(solved)
    assert solved_geometry_rejection_fitness(solved) >= 75.0


def test_compressed_channel_has_precursor_or_shortcut():
    nz, ny, nx = 16, 8, 48
    names = ["chi", "h11", "h33", "lapse", "shift1", "shift3"]
    n_comp = len(names)
    data = np.zeros((nz, ny, nx, n_comp), dtype=np.float64)
    idx = {n: i for i, n in enumerate(names)}
    data[:, :, :, idx["chi"]] = 1.0
    data[:, :, :, idx["h11"]] = 1.0
    data[:, :, :, idx["h33"]] = 1.0
    data[:, :, :, idx["lapse"]] = 1.0
    mid_j = ny // 2
    data[:, mid_j, nx // 4 : 3 * nx // 4, idx["chi"]] = 4.0
    data[:, mid_j, nx // 4 : 3 * nx // 4, idx["h11"]] = 4.0
    data[:, mid_j, nx // 4 : 3 * nx // 4, idx["h33"]] = 4.0
    dx = np.array([0.25, 0.25, 0.25])
    origin = np.array([-6.0, -1.0, -2.0])
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "channel.gridinit"
        write_gridinit(data, names, dx, origin, path)
        solved = compute_solved_geometry_ftl(path, n=61, L=8.0)
    assert solved is not None
    assert (
        solved_ftl_has_signal(solved)
        or solved.operational.max_local_speed > 1.0
        or solved.mechanisms.throat_pinch > 0.2
        or solved.mechanisms.portal_compression > 0.2
    )
