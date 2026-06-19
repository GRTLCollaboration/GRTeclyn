"""Tests for parallel Chombo→gridinit conversion."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pytest

from grteclyn_wrapper.grtresna import io
from grteclyn_wrapper.grtresna.io import chombo_to_uniform, convert_chombo_to_gridinit


def _paint_box_reference(
    data: np.ndarray,
    interior: np.ndarray,
    lo: np.ndarray,
    sz: np.ndarray,
    dx_lev: float,
    dx_target: np.ndarray,
    nx: int,
    ny: int,
    nz: int,
) -> None:
    """Original per-cell scatter the vectorized ``_paint_box`` must reproduce."""
    for kk in range(int(sz[2])):
        tk = io._target_span_slice(int(lo[2] + kk), 2, dx_lev, dx_target, nz)
        for jj in range(int(sz[1])):
            tj = io._target_span_slice(int(lo[1] + jj), 1, dx_lev, dx_target, ny)
            for ii in range(int(sz[0])):
                ti = io._target_span_slice(int(lo[0] + ii), 0, dx_lev, dx_target, nx)
                data[tk, tj, ti, :] = interior[:, kk, jj, ii]


@pytest.mark.parametrize(
    "dx_lev, dx_target, grid, lo, sz",
    [
        # finer source than target -> nearest-cell fallback (downsampling)
        (1.0, (2.0, 2.0, 2.0), (32, 32, 32), (0, 0, 0), (20, 18, 22)),
        (0.5, (2.0, 2.0, 2.0), (40, 40, 40), (10, 12, 8), (16, 16, 16)),
        # coarser source than target -> piecewise-constant upsampling
        (4.0, (2.0, 2.0, 2.0), (64, 64, 64), (2, 3, 1), (8, 7, 9)),
        # anisotropic target spacing
        (3.0, (2.0, 2.5, 1.7), (50, 40, 60), (1, 0, 4), (10, 9, 11)),
    ],
)
def test_vectorized_paint_box_matches_per_cell_loop(dx_lev, dx_target, grid, lo, sz) -> None:
    """The vectorized AMR paint must be byte-identical to the original O(sz^3)
    per-cell scatter across refinement, coarsening and anisotropic regimes."""
    nx, ny, nz = grid
    dxt = np.asarray(dx_target, dtype=np.float64)
    lo_a = np.asarray(lo, dtype=np.int64)
    sz_a = np.asarray(sz, dtype=np.int64)
    ncomp = 5
    rng = np.random.default_rng(0)
    interior = rng.standard_normal((ncomp, int(sz_a[2]), int(sz_a[1]), int(sz_a[0])))

    expected = np.zeros((nz, ny, nx, ncomp))
    actual = np.zeros((nz, ny, nx, ncomp))
    _paint_box_reference(expected, interior, lo_a, sz_a, dx_lev, dxt, nx, ny, nz)
    io._paint_box(actual, interior, lo_a, sz_a, dx_lev, dxt, nx, ny, nz)

    np.testing.assert_array_equal(actual, expected)

REPO_ROOT = Path(__file__).resolve().parents[3]
HDF5_CANDIDATES = [
    REPO_ROOT / "runs/grtresna_promote/fixture_hdf5/InitialDataFinal.3d.hdf5",
    REPO_ROOT
    / "runs/grtresna_promote/l128n256_qd_eval074/grtresna/Outputs/InitialDataFinal.3d.hdf5",
    REPO_ROOT
    / "runs/grtresna_qd/qd_20260608T175934Z/eval_000074/grtresna/Outputs/InitialDataFinal.3d.hdf5",
]


def _find_hdf5() -> Path | None:
    for path in HDF5_CANDIDATES:
        if path.is_file():
            return path
    return None


@pytest.mark.parametrize(
    "peak_shift, expect_stationary",
    [(0.20, False), (0.02, True)],
)
def test_stationarity_reads_seeded_gridinit_shift_magnitude(
    tmp_path: Path, peak_shift: float, expect_stationary: float
) -> None:
    """A GRTresna candidate has no ``recipe_beta_*`` profile, so stationarity must
    be read from the *real* seeded shift magnitude in the gridinit: a strongly
    sheared geometry is non-stationary, a near-zero one is stationary."""
    from grteclyn_wrapper.grtresna.io import write_gridinit
    from grteclyn_wrapper.metrics import STATIONARY_BETA_EPS
    from grteclyn_wrapper.metrics.diagnostics.comoving import read_comoving_metrics

    n = 16
    comp_names = ["chi", "shift1", "shift2", "shift3"]
    data = np.zeros((n, n, n, len(comp_names)), dtype=np.float64)
    data[..., 0] = 1.0  # chi
    # Localized shift bump at the grid centre, peak magnitude == peak_shift.
    data[n // 2, n // 2, n // 2, 1] = peak_shift
    dx_xyz = np.array([1.0, 1.0, 1.0])
    origin = np.array([0.0, 0.0, 0.0])
    write_gridinit(data, comp_names, dx_xyz, origin, tmp_path / "initial_data.gridinit")

    # GRTresna-style overrides: no recipe_beta_* keys.
    overrides = {"grtresna_shift_seed": peak_shift}
    cm = read_comoving_metrics(tmp_path, overrides, ftl_L=8.0)

    assert cm is not None
    assert cm.beta_mean == pytest.approx(peak_shift, abs=1e-9)
    assert cm.stationary is expect_stationary
    assert (peak_shift < STATIONARY_BETA_EPS) is expect_stationary


@pytest.mark.skipif(_find_hdf5() is None, reason="no GRTresna HDF5 fixture available")
def test_parallel_matches_serial(tmp_path: Path) -> None:
    hdf5 = _find_hdf5()
    assert hdf5 is not None

    serial, names_s, dx_s, origin_s = chombo_to_uniform(
        hdf5, nx=256, ny=256, nz=256, Lx=128.0, Ly=128.0, Lz=128.0, num_workers=1,
    )
    parallel, names_p, dx_p, origin_p = chombo_to_uniform(
        hdf5, nx=256, ny=256, nz=256, Lx=128.0, Ly=128.0, Lz=128.0, num_workers=8,
    )

    assert names_s == names_p
    np.testing.assert_allclose(dx_s, dx_p)
    np.testing.assert_allclose(origin_s, origin_p)
    np.testing.assert_allclose(serial, parallel, rtol=0.0, atol=0.0)

    out = tmp_path / "test.gridinit"
    convert_chombo_to_gridinit(
        hdf5,
        out,
        nx=256,
        ny=256,
        nz=256,
        L=128.0,
        num_workers=8,
    )
    assert out.is_file()
    assert out.stat().st_size > 0


@pytest.mark.skipif(_find_hdf5() is None, reason="no GRTresna HDF5 fixture available")
def test_cropped_export_matches_evolution_dx(tmp_path: Path) -> None:
    """Export the central evolution window at dx=0.5 (L=64, N=128)."""
    from grteclyn_wrapper.grtresna.io import read_gridinit

    hdf5 = _find_hdf5()
    assert hdf5 is not None
    out = tmp_path / "cropped.gridinit"
    convert_chombo_to_gridinit(
        hdf5,
        out,
        nx=128,
        ny=128,
        nz=128,
        Lx=64.0,
        Ly=64.0,
        Lz=64.0,
        source_origin=(32.0, 32.0, 32.0),
        target_center=(32.0, 32.0, 32.0),
        num_workers=1,
    )
    gi = read_gridinit(out)
    np.testing.assert_allclose(gi.dx_xyz, [0.5, 0.5, 0.5])
    np.testing.assert_allclose(gi.origin, [0.0, 0.0, 0.0])
    assert gi.data.shape[:3] == (128, 128, 128)


@pytest.mark.skipif(_find_hdf5() is None, reason="no GRTresna HDF5 fixture available")
def test_parallel_conversion_is_faster(tmp_path: Path) -> None:
    hdf5 = _find_hdf5()
    assert hdf5 is not None

    # Use the fixture's native resolution; large HQ runs (4096 level-0 boxes)
    # benefit more, but correctness is the hard requirement here.
    t0 = time.perf_counter()
    chombo_to_uniform(
        hdf5, nx=128, ny=128, nz=128, Lx=128.0, Ly=128.0, Lz=128.0, num_workers=1,
    )
    serial_s = time.perf_counter() - t0

    t0 = time.perf_counter()
    chombo_to_uniform(
        hdf5, nx=128, ny=128, nz=128, Lx=128.0, Ly=128.0, Lz=128.0, num_workers=16,
    )
    parallel_s = time.perf_counter() - t0

    assert parallel_s <= serial_s * 1.05, (
        f"parallel conversion slower than expected: serial={serial_s:.1f}s "
        f"parallel={parallel_s:.1f}s"
    )
