"""Tests for parallel Chombo→gridinit conversion."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pytest

from grteclyn_wrapper.grtresna.io import chombo_to_uniform, convert_chombo_to_gridinit

REPO_ROOT = Path(__file__).resolve().parents[2]
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
