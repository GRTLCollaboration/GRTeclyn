"""Read GRTresna Chombo-HDF5 output and write a flat binary grid file
that GRTeclyn's ExternalGridInitialData can read.

The Chombo HDF5 layout stores AMR checkpoint data as per-level box
arrays with offsets.  Each box includes *ghost cells* (typically 3) so
the stored array per box is (sz+2*ghost)^3, not sz^3.  We strip
ghosts, then flatten to a single uniform Cartesian grid by painting
coarsest-to-finest so finer AMR data overwrites coarser.

Output format (.gridinit):
    Header (ASCII, newline-terminated lines):
        GRTECLYN_GRID_INIT_V1
        num_components <int>
        component_names <space-separated names>
        nx_ny_nz <int> <int> <int>
        dx <float>
        origin <float> <float> <float>
        END_HEADER
    Body (binary, C-order float64):
        data[nz][ny][nx][num_components]
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import h5py
import numpy as np
from numpy.typing import NDArray

GRTECLYN_STATE_VARS = [
    "chi",
    "h11", "h12", "h13", "h22", "h23", "h33",
    "K",
    "A11", "A12", "A13", "A22", "A23", "A33",
    "Theta",
    "Gamma1", "Gamma2", "Gamma3",
    "lapse",
    "shift1", "shift2", "shift3",
    "B1", "B2", "B3",
    "phi", "Pi",
]


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


def _paint_level(
    data: NDArray,
    f: h5py.File,
    level: int,
    n_comp: int,
    dx_target: float,
    nx: int,
    ny: int,
    nz: int,
) -> None:
    """Paint one AMR level onto the uniform target grid (vectorized)."""
    grp = f[f"level_{level}"]
    dx_lev = float(grp.attrs["dx"])
    raw = grp["data:datatype=0"][:]
    offsets = grp["data:offsets=0"][:]
    box_records = grp["boxes"][:]

    ghost = None  # infer once from first box

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
        # Chombo: [comp][z_padded][y_padded][x_padded]
        arr = chunk.reshape(n_comp, int(padded[2]), int(padded[1]), int(padded[0]))

        # Strip ghost cells to get interior data only
        g = ghost
        interior = arr[:, g:g+int(sz[2]), g:g+int(sz[1]), g:g+int(sz[0])]
        # interior shape: (n_comp, sz_z, sz_y, sz_x)

        # Physical cell-center coordinates
        ix = np.arange(sz[0])
        iy = np.arange(sz[1])
        iz = np.arange(sz[2])

        phys_x = (lo[0] + ix + 0.5) * dx_lev
        phys_y = (lo[1] + iy + 0.5) * dx_lev
        phys_z = (lo[2] + iz + 0.5) * dx_lev

        ti = np.clip((phys_x / dx_target).astype(np.int64), 0, nx - 1)
        tj = np.clip((phys_y / dx_target).astype(np.int64), 0, ny - 1)
        tk = np.clip((phys_z / dx_target).astype(np.int64), 0, nz - 1)

        TK, TJ, TI = np.meshgrid(tk, tj, ti, indexing="ij")
        for c in range(n_comp):
            data[TK, TJ, TI, c] = interior[c]


def chombo_to_uniform(
    chombo_path: str | Path,
    nx: int,
    ny: int,
    nz: int,
    L: float,
) -> tuple[NDArray, list[str], float, NDArray]:
    """Read a Chombo HDF5 file and flatten to a uniform grid.

    Parameters
    ----------
    chombo_path : path to InitialDataFinal.3d.hdf5
    nx, ny, nz : target uniform grid resolution
    L : domain side length (must match the Chombo run)

    Returns
    -------
    data : (nz, ny, nx, n_comp) float64 array  [C-order: z outermost]
    comp_names : list of component names
    dx : cell spacing
    origin : (3,) array
    """
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

        chombo_nz_cells = int(prob_domain[5]) + 1
        Lz_chombo = chombo_nz_cells * dx_coarse

        dx_target = L / nx

        data = np.zeros((nz, ny, nx, n_comp), dtype=np.float64)

        # Paint coarsest first, then finer levels overwrite
        for lev in range(num_levels):
            _paint_level(data, f, lev, n_comp, dx_target, nx, ny, nz)

        # z-reflection: if Chombo used half-z (reflective BC at z=0)
        if Lz_chombo < L * 0.75:
            mid_k = nz // 2
            reflected = data[:mid_k, :, :, :].copy()[::-1]

            z_odd_names = {
                "h13", "h23", "A13", "A23",
                "Gamma3", "shift3", "B3",
            }
            for i, name in enumerate(comp_names):
                if name in z_odd_names:
                    reflected[:, :, :, i] *= -1.0

            data[mid_k:, :, :, :] = reflected

    origin = np.array([0.0, 0.0, 0.0])
    return data, comp_names, dx_target, origin


def write_gridinit(
    data: NDArray,
    comp_names: list[str],
    dx: float,
    origin: NDArray,
    output_path: str | Path,
) -> Path:
    """Write the uniform grid to a .gridinit binary file."""
    output_path = Path(output_path)
    nz, ny, nx, n_comp = data.shape

    with open(output_path, "wb") as fout:
        header_lines = [
            "GRTECLYN_GRID_INIT_V1",
            f"num_components {n_comp}",
            f"component_names {' '.join(comp_names)}",
            f"nx_ny_nz {nx} {ny} {nz}",
            f"dx {dx:.15e}",
            f"origin {origin[0]:.15e} {origin[1]:.15e} {origin[2]:.15e}",
            "END_HEADER",
        ]
        header_bytes = ("\n".join(header_lines) + "\n").encode("ascii")
        fout.write(header_bytes)
        fout.write(data.tobytes())

    return output_path


def convert_chombo_to_gridinit(
    chombo_path: str | Path,
    output_path: str | Path,
    N: int = 64,
    L: float = 128.0,
) -> Path:
    """One-shot: read Chombo HDF5, flatten, write .gridinit.

    N is the number of cells per side (uniform cube N x N x N).
    L is the domain side length.
    """
    data, comp_names, dx, origin = chombo_to_uniform(
        chombo_path, nx=N, ny=N, nz=N, L=L,
    )
    return write_gridinit(data, comp_names, dx, origin, output_path)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <chombo.hdf5> <output.gridinit> [N] [L]")
        sys.exit(1)
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 64
    L = float(sys.argv[4]) if len(sys.argv) > 4 else 128.0
    out = convert_chombo_to_gridinit(sys.argv[1], sys.argv[2], N=N, L=L)
    print(f"Wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")
