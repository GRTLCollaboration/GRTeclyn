"""Read GRTresna Chombo-HDF5 output and write a flat binary grid file
that GRTeclyn's ExternalGridInitialData can read.

The Chombo HDF5 layout stores AMR checkpoint data as per-level box
arrays with offsets.  Each box includes *ghost cells* (typically 3) so
the stored array per box is (sz+2*ghost)^3, not sz^3.  We strip
ghosts, then flatten to a single uniform Cartesian grid by painting
coarsest-to-finest so finer AMR data overwrites coarser.

Supports non-cubic target grids: nx, ny, nz can differ, and per-axis
domain lengths Lx, Ly, Lz can be specified independently.  The output
.gridinit stores per-axis cell spacings (dx_x, dx_y, dx_z).

Output format (.gridinit v2):
    Header (ASCII, newline-terminated lines):
        GRTECLYN_GRID_INIT_V2
        num_components <int>
        component_names <space-separated names>
        nx_ny_nz <int> <int> <int>
        dx <float> <float> <float>
        origin <float> <float> <float>
        END_HEADER
    Body (binary, C-order float64):
        data[nz][ny][nx][num_components]
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

# h5py is imported lazily inside the functions that need it so that importing
# this module (e.g. to construct a GRTresnaConfig or for a dry-run) does not
# require h5py to be installed.

logger = logging.getLogger(__name__)

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
    "teo_rho", "teo_j1", "teo_j2", "teo_j3",
    "teo_S11", "teo_S12", "teo_S13", "teo_S22", "teo_S23", "teo_S33",
]

_Z_ODD_NAMES = frozenset({
    "h13", "h23", "A13", "A23",
    "Gamma3", "shift3", "B3",
    "teo_j3", "teo_S13", "teo_S23",
})


def _reflect_half_z_to_full(data: NDArray, comp_names: list[str]) -> None:
    """Move stored z>=0 data to the upper half and parity-reflect z<0 in place."""
    mid_k = data.shape[0] // 2
    positive_z = data[:mid_k, :, :, :].copy()
    reflected = positive_z[::-1].copy()
    for i, name in enumerate(comp_names):
        if name in _Z_ODD_NAMES:
            reflected[:, :, :, i] *= -1.0
    data[:mid_k, :, :, :] = reflected
    data[mid_k:mid_k + positive_z.shape[0], :, :, :] = positive_z


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
    dx_target: NDArray,
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

    ghost = None

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

        g = ghost
        interior = arr[:, g:g+int(sz[2]), g:g+int(sz[1]), g:g+int(sz[0])]

        ix = np.arange(sz[0])
        iy = np.arange(sz[1])
        iz = np.arange(sz[2])

        phys_x = (lo[0] + ix + 0.5) * dx_lev
        phys_y = (lo[1] + iy + 0.5) * dx_lev
        phys_z = (lo[2] + iz + 0.5) * dx_lev

        ti = np.clip((phys_x / dx_target[0]).astype(np.int64), 0, nx - 1)
        tj = np.clip((phys_y / dx_target[1]).astype(np.int64), 0, ny - 1)
        tk = np.clip((phys_z / dx_target[2]).astype(np.int64), 0, nz - 1)

        TK, TJ, TI = np.meshgrid(tk, tj, ti, indexing="ij")
        for c in range(n_comp):
            data[TK, TJ, TI, c] = interior[c]


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
) -> tuple[NDArray, list[str], NDArray, NDArray]:
    """Read a Chombo HDF5 file and flatten to a uniform grid.

    Parameters
    ----------
    chombo_path : path to InitialDataFinal.3d.hdf5
    nx, ny, nz : target uniform grid resolution per axis
    Lx, Ly, Lz : domain extents per axis. If None, auto-detected from
                  the Chombo header. For half-z domains with reflection,
                  set Lz to the full (doubled) extent.

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

        data = np.zeros((nz, ny, nx, n_comp), dtype=np.float64)

        for lev in range(num_levels):
            _paint_level(data, f, lev, n_comp, dx_xyz, nx, ny, nz)

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


def write_gridinit(
    data: NDArray,
    comp_names: list[str],
    dx_xyz: NDArray,
    origin: NDArray,
    output_path: str | Path,
) -> Path:
    """Write the uniform grid to a .gridinit binary file.

    dx_xyz may be a scalar (isotropic) or a 3-element array (per-axis).
    """
    output_path = Path(output_path)
    nz, ny, nx, n_comp = data.shape

    dx_arr = np.atleast_1d(np.asarray(dx_xyz, dtype=np.float64))
    if dx_arr.size == 1:
        dx_arr = np.full(3, dx_arr[0])

    with open(output_path, "wb") as fout:
        header_lines = [
            "GRTECLYN_GRID_INIT_V2",
            f"num_components {n_comp}",
            f"component_names {' '.join(comp_names)}",
            f"nx_ny_nz {nx} {ny} {nz}",
            f"dx {dx_arr[0]:.15e} {dx_arr[1]:.15e} {dx_arr[2]:.15e}",
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
    *,
    nx: int | None = None,
    ny: int | None = None,
    nz: int | None = None,
    N: int | None = None,
    Lx: float | None = None,
    Ly: float | None = None,
    Lz: float | None = None,
    L: float | None = None,
    target_center: tuple[float, float, float] | None = None,
    delete_source: bool = False,
) -> Path:
    """One-shot: read Chombo HDF5, flatten, write .gridinit.

    Grid resolution: specify ``nx, ny, nz`` individually, or ``N`` for a
    uniform cube.  If all are None, defaults to ``N=64``.

    Domain size: specify ``Lx, Ly, Lz`` individually, or ``L`` for a
    uniform cube.  If all are None, auto-detected from the Chombo header
    (with Lz doubled if the source used half-z reflection).

    If *delete_source* is True, the Chombo HDF5 file is deleted after a
    successful conversion to reclaim disk space.
    """
    if N is not None:
        nx = nx or N
        ny = ny or N
        nz = nz or N
    nx = nx or 64
    ny = ny or 64
    nz = nz or 64

    if L is not None:
        Lx = Lx or L
        Ly = Ly or L
        Lz = Lz or L

    # Auto-detect domain and double Lz if Chombo used half-z
    if Lx is None or Ly is None or Lz is None:
        _, _, chombo_lengths = read_chombo_domain(chombo_path)
        if Lx is None:
            Lx = chombo_lengths[0]
        if Ly is None:
            Ly = chombo_lengths[1]
        if Lz is None:
            # If Chombo Lz is roughly half of Lx, assume half-z reflection
            if chombo_lengths[2] < chombo_lengths[0] * 0.75:
                Lz = chombo_lengths[0]  # full domain = Lx
                logger.info(
                    "Auto-detected half-z domain (Chombo Lz=%.1f < Lx=%.1f); "
                    "setting Lz=%.1f with z-reflection",
                    chombo_lengths[2], chombo_lengths[0], Lz,
                )
            else:
                Lz = chombo_lengths[2]

    data, comp_names, dx_xyz, origin = chombo_to_uniform(
        chombo_path, nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz,
    )

    # Coordinate alignment. The GRTresna matter sits at the centre of the solve
    # domain (physical (Lx/2, Ly/2, Lz/2), stored at the centre of the gridinit
    # array). GRTeclyn's loader samples absolute coordinates px=(i+0.5)*dx and
    # indexes the gridinit via (px - origin)/file_dx. Setting
    #   origin = target_center - L/2
    # makes GRTeclyn's box centre (target_center) map onto the gridinit centre,
    # so GRTeclyn evolves the central window of the GRTresna domain rather than
    # an off-centre corner. For the half-z RadialRecipe box target_center_z = 0
    # correctly maps the z=0 reflective plane onto the gridinit centre.
    if target_center is not None:
        origin = np.array([
            float(target_center[0]) - 0.5 * float(Lx),
            float(target_center[1]) - 0.5 * float(Ly),
            float(target_center[2]) - 0.5 * float(Lz),
        ])

    result = write_gridinit(data, comp_names, dx_xyz, origin, output_path)

    if delete_source:
        src = Path(chombo_path)
        if src.exists():
            src.unlink()
            logger.info("Deleted source HDF5: %s", src)

    return result


@dataclass(frozen=True)
class GridinitData:
    """Uniform grid from a .gridinit v2 file."""

    data: NDArray  # shape (nz, ny, nx, n_comp)
    comp_names: list[str]
    dx_xyz: NDArray  # (3,)
    origin: NDArray  # (3,)


def read_gridinit(path: str | Path) -> GridinitData:
    """Read a GRTECLYN_GRID_INIT_V2 binary written by :func:`write_gridinit`."""
    path = Path(path)
    with open(path, "rb") as fin:
        lines: list[str] = []
        while True:
            raw = fin.readline()
            if not raw:
                raise ValueError(f"Unexpected EOF before END_HEADER in {path}")
            line = raw.decode("ascii").strip()
            lines.append(line)
            if line == "END_HEADER":
                break

        if not lines or lines[0] != "GRTECLYN_GRID_INIT_V2":
            raise ValueError(f"Not a GRTECLYN_GRID_INIT_V2 file: {path}")

        num_components = 0
        comp_names: list[str] = []
        nx = ny = nz = 0
        dx_xyz = np.zeros(3, dtype=np.float64)
        origin = np.zeros(3, dtype=np.float64)

        for line in lines[1:-1]:
            parts = line.split()
            if not parts:
                continue
            key = parts[0]
            if key == "num_components":
                num_components = int(parts[1])
            elif key == "component_names":
                comp_names = parts[1:]
            elif key == "nx_ny_nz":
                nx, ny, nz = int(parts[1]), int(parts[2]), int(parts[3])
            elif key == "dx":
                dx_xyz = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
            elif key == "origin":
                origin = np.array([float(parts[1]), float(parts[2]), float(parts[3])])

        if num_components <= 0 or not comp_names or nx <= 0 or ny <= 0 or nz <= 0:
            raise ValueError(f"Incomplete gridinit header in {path}")

        body = fin.read()
    expected = nz * ny * nx * num_components * 8
    if len(body) != expected:
        raise ValueError(
            f"gridinit body size mismatch in {path}: got {len(body)} bytes, "
            f"expected {expected}"
        )
    data = np.frombuffer(body, dtype=np.float64).reshape(nz, ny, nx, num_components)
    return GridinitData(
        data=data.copy(),
        comp_names=comp_names,
        dx_xyz=dx_xyz,
        origin=origin,
    )


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO)
    if len(sys.argv) < 3:
        print(
            f"Usage: {sys.argv[0]} <chombo.hdf5> <output.gridinit> "
            "[nx] [ny] [nz] [Lx] [Ly] [Lz]"
        )
        sys.exit(1)

    args = sys.argv[3:]
    kw: dict = {}
    if len(args) >= 3:
        kw["nx"], kw["ny"], kw["nz"] = int(args[0]), int(args[1]), int(args[2])
    elif len(args) == 1:
        kw["N"] = int(args[0])
    if len(args) >= 6:
        kw["Lx"], kw["Ly"], kw["Lz"] = float(args[3]), float(args[4]), float(args[5])

    out = convert_chombo_to_gridinit(
        sys.argv[1], sys.argv[2], delete_source=False, **kw,
    )
    print(f"Wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")
