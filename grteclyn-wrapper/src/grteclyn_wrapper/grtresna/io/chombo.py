"""Chombo HDF5 AMR reading and uniform-grid painting."""

from __future__ import annotations

import logging
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from .gridinit import _reflect_half_z_to_full

if TYPE_CHECKING:
    import h5py

logger = logging.getLogger(__name__)

def _infer_ghost_cells(
    n_doubles: int, box_size: NDArray, n_comp: int,
) -> int:
    """Determine the ghost-cell count from the stored data size."""
    for g in range(8):
        padded = int(np.prod(box_size + 2 * g)) * n_comp
        if padded == n_doubles:
            return g
    raise ValueError(
        f"Cannot infer ghost width: n_doubles={n_doubles}, "
        f"box_size={box_size.tolist()}, n_comp={n_comp}"
    )


def _default_gridinit_workers() -> int:
    """Thread count for AMR box painting (I/O-bound slice writes release the GIL)."""
    return min(32, os.cpu_count() or 1)


def _target_span_slice(
    src_index: int,
    axis: int,
    dx_lev: float,
    dx_target: NDArray,
    n_target: int,
    target_origin: float = 0.0,
) -> slice:
    """Target cell-center indices covered by one Chombo source cell."""
    dx = float(dx_target[axis])
    left = src_index * dx_lev
    right = (src_index + 1) * dx_lev
    start = int(np.ceil((left - target_origin) / dx - 0.5 - 1.0e-12))
    stop = int(np.floor((right - target_origin) / dx - 0.5 + 1.0e-12)) + 1
    start = max(0, min(start, n_target))
    stop = max(start, min(stop, n_target))
    if start == stop:
        nearest = int(
            np.clip(((left + right) * 0.5 - target_origin) / dx - 0.5, 0, n_target - 1)
        )
        return slice(nearest, nearest + 1)
    return slice(start, stop)


def _is_aligned_1to1(dx_lev: float, dx_target: NDArray) -> bool:
    return all(
        abs(dx_lev - float(dx_target[axis])) < 1.0e-12 * max(dx_lev, 1.0e-15)
        for axis in range(3)
    )


def _axis_source_map(
    lo: int,
    sz: int,
    axis: int,
    dx_lev: float,
    dx_target: NDArray,
    n_target: int,
    target_origin: float = 0.0,
) -> NDArray:
    """Map every target cell index on one axis to the local source-cell index
    that the piecewise-constant paint assigns it (``-1`` where untouched).

    The 3-D scatter is separable: a target voxel is painted by source cell
    ``(kk, jj, ii)`` iff its three axis indices each land in the corresponding
    1-D source span, and the coarse-to-fine loop lets the highest-index source
    win each contested cell.  Building the per-axis winner map in O(sz) and
    taking the tensor product reproduces the old per-cell loop exactly, without
    its O(sz^3) Python cost.
    """
    src_of_target = np.full(int(n_target), -1, dtype=np.int64)
    for local in range(int(sz)):
        span = _target_span_slice(
            int(lo + local), axis, dx_lev, dx_target, n_target, target_origin,
        )
        src_of_target[span] = local
    return src_of_target


def _paint_box(
    data: NDArray,
    interior: NDArray,
    lo: NDArray,
    sz: NDArray,
    dx_lev: float,
    dx_target: NDArray,
    nx: int,
    ny: int,
    nz: int,
    target_origin: NDArray | None = None,
) -> None:
    """Paint one AMR box interior onto the uniform target grid."""
    origin = np.zeros(3, dtype=np.float64) if target_origin is None else np.asarray(
        target_origin, dtype=np.float64,
    )
    use_fast_path = _is_aligned_1to1(dx_lev, dx_target) and np.allclose(origin, 0.0)
    if use_fast_path:
        ti = slice(int(lo[0]), int(lo[0]) + int(sz[0]))
        tj = slice(int(lo[1]), int(lo[1]) + int(sz[1]))
        tk = slice(int(lo[2]), int(lo[2]) + int(sz[2]))
        data[tk, tj, ti, :] = np.moveaxis(interior, 0, -1)
        return

    map_i = _axis_source_map(
        int(lo[0]), int(sz[0]), 0, dx_lev, dx_target, nx, float(origin[0]),
    )
    map_j = _axis_source_map(
        int(lo[1]), int(sz[1]), 1, dx_lev, dx_target, ny, float(origin[1]),
    )
    map_k = _axis_source_map(
        int(lo[2]), int(sz[2]), 2, dx_lev, dx_target, nz, float(origin[2]),
    )

    tk = np.nonzero(map_k >= 0)[0]
    tj = np.nonzero(map_j >= 0)[0]
    ti = np.nonzero(map_i >= 0)[0]
    if tk.size == 0 or tj.size == 0 or ti.size == 0:
        return

    gathered = interior[
        :,
        map_k[tk][:, None, None],
        map_j[tj][None, :, None],
        map_i[ti][None, None, :],
    ]
    data[np.ix_(tk, tj, ti)] = np.moveaxis(gathered, 0, -1)


def _extract_level_boxes(
    raw: NDArray,
    offsets: NDArray,
    box_records: NDArray,
    n_comp: int,
) -> list[tuple[NDArray, NDArray, NDArray]]:
    """Unpack all AMR boxes at one level from the flat Chombo buffer."""
    boxes: list[tuple[NDArray, NDArray, NDArray]] = []
    ghost: int | None = None

    for bi in range(len(box_records)):
        rec = box_records[bi]
        lo = np.array([rec["lo_i"], rec["lo_j"], rec["lo_k"]], dtype=np.int64)
        hi = np.array([rec["hi_i"], rec["hi_j"], rec["hi_k"]], dtype=np.int64)
        sz = (hi - lo + 1).astype(np.int64)
        start = int(offsets[bi])
        end = int(offsets[bi + 1]) if bi + 1 < len(offsets) else len(raw)
        n_doubles = end - start

        if ghost is None:
            ghost = _infer_ghost_cells(n_doubles, sz, n_comp)

        padded = sz + 2 * ghost
        chunk = raw[start : start + int(np.prod(padded)) * n_comp]
        arr = chunk.reshape(n_comp, int(padded[2]), int(padded[1]), int(padded[0]))
        interior = arr[
            :,
            ghost : ghost + int(sz[2]),
            ghost : ghost + int(sz[1]),
            ghost : ghost + int(sz[0]),
        ]
        boxes.append((lo, sz, interior))

    return boxes


def _paint_level(
    data: NDArray,
    f: h5py.File,
    level: int,
    n_comp: int,
    dx_target: NDArray,
    nx: int,
    ny: int,
    nz: int,
    *,
    target_origin: NDArray | None = None,
    num_workers: int = 1,
) -> None:
    """Paint one AMR level onto the uniform target grid."""
    grp = f[f"level_{level}"]
    dx_lev = float(grp.attrs["dx"])
    raw = grp["data:datatype=0"][:]
    offsets = grp["data:offsets=0"][:]
    box_records = grp["boxes"][:]
    boxes = _extract_level_boxes(raw, offsets, box_records, n_comp)

    def _paint_one(box: tuple[NDArray, NDArray, NDArray]) -> None:
        lo, sz, interior = box
        _paint_box(
            data, interior, lo, sz, dx_lev, dx_target, nx, ny, nz, target_origin,
        )

    # Parallel box painting is only safe when each source cell maps to a
    # disjoint target region (1:1 aligned spacing). Finer AMR levels use
    # nearest-cell fallback where adjacent boxes can contend for the same
    # target index; those levels must paint sequentially to preserve the
    # coarse-to-fine overwrite order.
    can_parallel = _is_aligned_1to1(dx_lev, dx_target)
    if num_workers <= 1 or len(boxes) <= 1 or not can_parallel:
        for box in boxes:
            _paint_one(box)
        return

    with ThreadPoolExecutor(max_workers=num_workers) as pool:
        list(pool.map(_paint_one, boxes))


def read_chombo_domain(
    chombo_path: str | Path,
) -> tuple[tuple[int, int, int], float, tuple[float, float, float]]:
    """Read the coarse-grid dimensions and domain size from a Chombo HDF5.

    Returns (N_cells_xyz, dx_coarse, L_xyz) where N_cells is per-axis
    cell count at level 0 and L_xyz is per-axis physical extent.
    """
    import h5py

    with h5py.File(str(chombo_path), "r") as f:
        l0 = f["level_0"]
        prob = l0.attrs["prob_domain"]
        dx = float(l0.attrs["dx"])
        ncells = (int(prob[3]) + 1, int(prob[4]) + 1, int(prob[5]) + 1)
        lengths = (ncells[0] * dx, ncells[1] * dx, ncells[2] * dx)
    return ncells, dx, lengths


def chombo_to_uniform(
    chombo_path: str | Path,
    nx: int,
    ny: int,
    nz: int,
    Lx: float | None = None,
    Ly: float | None = None,
    Lz: float | None = None,
    *,
    source_origin: Sequence[float] | tuple[float, float, float] | None = None,
    num_workers: int = 0,
) -> tuple[NDArray, list[str], NDArray, NDArray]:
    """Read a Chombo HDF5 file and flatten to a uniform grid.

    Parameters
    ----------
    chombo_path : path to InitialDataFinal.3d.hdf5
    nx, ny, nz : target uniform grid resolution per axis
    Lx, Ly, Lz : domain extents per axis. If None, auto-detected from
                  the Chombo header. For half-z domains with reflection,
                  set Lz to the full (doubled) extent.
    source_origin : minimum corner of the export window in GRTresna native
                    coordinates (default ``(0, 0, 0)``).

    Returns
    -------
    data : (nz, ny, nx, n_comp) float64 array
    comp_names : list of component names
    dx_xyz : (3,) array of per-axis cell spacings
    origin : (3,) array
    """
    import h5py

    chombo_path = Path(chombo_path)
    with h5py.File(chombo_path, "r") as f:
        n_comp = int(f.attrs["num_components"])
        num_levels = int(f.attrs["num_levels"])

        comp_names: list[str] = [""] * n_comp
        for key, val in f.attrs.items():
            if key.startswith("component_"):
                idx = int(key.split("_")[1])
                comp_names[idx] = val if isinstance(val, str) else val.decode()

        l0 = f["level_0"]
        prob_domain = l0.attrs["prob_domain"]
        dx_coarse = float(l0.attrs["dx"])

        chombo_ncells = np.array([
            int(prob_domain[3]) + 1,
            int(prob_domain[4]) + 1,
            int(prob_domain[5]) + 1,
        ])
        chombo_L = chombo_ncells * dx_coarse

        if Lx is None:
            Lx = float(chombo_L[0])
        if Ly is None:
            Ly = float(chombo_L[1])
        if Lz is None:
            Lz = float(chombo_L[2])

        dx_xyz = np.array([Lx / nx, Ly / ny, Lz / nz])
        target_origin = np.zeros(3, dtype=np.float64)
        if source_origin is not None:
            target_origin = np.asarray(source_origin, dtype=np.float64).reshape(3)

        data = np.zeros((nz, ny, nx, n_comp), dtype=np.float64)
        workers = num_workers if num_workers > 0 else _default_gridinit_workers()

        for lev in range(num_levels):
            _paint_level(
                data,
                f,
                lev,
                n_comp,
                dx_xyz,
                nx,
                ny,
                nz,
                target_origin=target_origin,
                num_workers=workers,
            )

        # z-reflection: if Chombo used half-z (reflective BC at z=0).
        #
        # _paint_level places Chombo's stored half-domain (z >= 0) in the first
        # half of `data` because Chombo physical coordinates start at zero.  The
        # .gridinit file, however, is written with origin_z = target_center_z -
        # Lz/2 so that GRTeclyn's z=0 reflective plane samples the middle of the
        # file.  Therefore the positive-z Chombo data must occupy the upper half
        # of the full reflected array, with the lower half filled by parity
        # reflection.  Putting the stored half-domain in data[:mid_k] directly
        # makes GRTeclyn's z=0 slice sample the far z-boundary tail.
        Lz_chombo = float(chombo_L[2])
        if Lz_chombo < Lz * 0.75:
            _reflect_half_z_to_full(data, comp_names)

    origin = np.array([0.0, 0.0, 0.0])
    return data, comp_names, dx_xyz, origin
