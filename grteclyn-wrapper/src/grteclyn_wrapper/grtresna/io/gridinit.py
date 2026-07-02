"""Gridinit v2 binary I/O."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

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
