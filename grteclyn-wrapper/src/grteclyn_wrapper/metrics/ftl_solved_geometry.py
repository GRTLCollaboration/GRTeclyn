"""Operational FTL and mechanism descriptors on GRTresna-solved initial data.

Reads a ``.gridinit`` produced by the elliptic solve and evaluates the same
mechanism-agnostic Dijkstra probe used for evolved plotfiles, plus cheap
classifiers for how the shortcut is produced (shift-warp vs throat vs portal).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from ..grtresna.io import GridinitData, read_gridinit
from .ftl_general import GeneralFtlReport, operational_ftl_on_grid

MECHANISM_MIN_SCORE: float = 0.05


@dataclass(frozen=True)
class MechanismScores:
    """Cheap t=0 mechanism proxies in [0, 1] (not mutually exclusive)."""

    shift_warp: float
    throat_pinch: float
    portal_compression: float

    @property
    def dominant_index(self) -> int | None:
        """0=shift, 1=throat, 2=portal; None when no proxy clears the floor."""
        scores = (self.shift_warp, self.throat_pinch, self.portal_compression)
        if max(scores) < MECHANISM_MIN_SCORE:
            return None
        return int(np.argmax(scores))

    @property
    def mechanism_descriptor(self) -> float:
        """Single [0,1] behavior axis for MAP-Elites: 0 shift .. 1 portal.

        Continuous blend of the proxies that clear the floor, so distinct
        mechanism mixes land in distinct cells instead of collapsing onto the
        argmax (which clustered every candidate at one value).  Returns 0.0
        when no proxy is active (caller treats it as "no mechanism")."""
        weights = np.array(
            [self.shift_warp, self.throat_pinch, self.portal_compression],
            dtype=float,
        )
        weights = np.where(weights >= MECHANISM_MIN_SCORE, weights, 0.0)
        total = float(weights.sum())
        if total <= 0.0:
            return 0.0
        axis = np.array([0.0, 0.5, 1.0])
        return float(np.dot(weights, axis) / total)


@dataclass(frozen=True)
class SolvedGeometryFtl:
    operational: GeneralFtlReport
    mechanisms: MechanismScores
    integration_L: float


def _comp_index(names: Sequence[str], key: str) -> int | None:
    try:
        return list(names).index(key)
    except ValueError:
        return None


def build_xz_slice_from_gridinit(
    grid: GridinitData,
    *,
    n: int = 129,
    L: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], tuple[float, float]]:
    """Extract x-z midplane (y = domain centre) for operational FTL."""
    data = grid.data
    nz, ny, nx, _ = data.shape
    names = grid.comp_names
    dx, dy, dz = grid.dx_xyz
    ox, oy, oz = grid.origin

    cx = ox + 0.5 * nx * dx
    cy = oy + 0.5 * ny * dy
    cz = oz + 0.5 * nz * dz

    half_x = 0.49 * nx * dx
    half_z = 0.49 * nz * dz
    if L is not None:
        half_x = min(half_x, float(L))
        half_z = min(half_z, float(L))

    x_axis = np.linspace(cx - half_x, cx + half_x, n)
    z_axis = np.linspace(cz - half_z, cz + half_z, n)

    # Nearest-neighbour source indices, precomputed once and shared across all
    # components (vectorised: the per-cell Python loop is the search-loop
    # bottleneck, so this stays a few array ops).
    ii = np.clip(np.round((x_axis - ox) / dx - 0.5).astype(np.int64), 0, nx - 1)
    kk = np.clip(np.round((z_axis - oz) / dz - 0.5).astype(np.int64), 0, nz - 1)
    jj = int(np.clip(round((cy - oy) / dy - 0.5), 0, ny - 1))
    # Index grids with image layout [i (x), k (z)].
    ki = kk[None, :].repeat(n, axis=0)
    xi = ii[:, None].repeat(n, axis=1)

    def sample(comp: str, default: float = 0.0) -> NDArray[np.float64]:
        c = _comp_index(names, comp)
        if c is None:
            return np.full((n, n), default, dtype=np.float64)
        return data[ki, jj, xi, c]

    if _comp_index(names, "chi") is None or _comp_index(names, "lapse") is None:
        raise KeyError("gridinit must include chi and lapse")
    chi = np.clip(sample("chi"), 1.0e-10, None)
    alpha = np.clip(sample("lapse"), 1.0e-10, None)
    beta_x = sample("shift1", 0.0)
    beta_z = sample("shift3", 0.0)

    inv_chi = 1.0 / chi
    h11 = sample("h11", 1.0) * inv_chi
    h13 = sample("h13", 0.0) * inv_chi
    h33 = sample("h33", 1.0) * inv_chi

    beta = np.zeros((n, n, 2), dtype=np.float64)
    beta[:, :, 0] = beta_x
    beta[:, :, 1] = beta_z

    gamma = np.zeros((n, n, 2, 2), dtype=np.float64)
    gamma[:, :, 0, 0] = h11
    gamma[:, :, 0, 1] = gamma[:, :, 1, 0] = h13
    gamma[:, :, 1, 1] = h33

    spacing = (x_axis[1] - x_axis[0], z_axis[1] - z_axis[0])
    return alpha, beta, gamma, spacing


def _mechanism_scores(
    alpha: NDArray[np.float64],
    beta: NDArray[np.float64],
    gamma: NDArray[np.float64],
    spacing: tuple[float, float],
    *,
    integration_L: float,
) -> MechanismScores:
    """Classify shortcut style from the solved slice fields."""
    chi = 1.0 / np.clip(gamma[:, :, 0, 0], 1.0e-10, None)
    beta_mag = np.sqrt(beta[:, :, 0] ** 2 + beta[:, :, 1] ** 2)
    beta_max = float(np.max(beta_mag))
    shift_warp = float(np.clip(beta_max / 0.5, 0.0, 1.0))

    mid = chi.shape[1] // 2
    x = (np.arange(chi.shape[0]) - chi.shape[0] / 2) * spacing[0]
    chi_ax = np.clip(chi[:, mid], 1.0e-10, None)
    r_pos = np.abs(x)
    r_areal = r_pos / np.sqrt(chi_ax)
    r_asymp = integration_L / max(float(np.sqrt(chi_ax[-1])), 1.0e-8)
    f_throat = 0.0
    if r_asymp > 1.0e-8 and len(r_areal) >= 5:
        for idx in range(1, len(r_areal) - 1):
            if r_pos[idx] < 1.0e-6:
                continue
            if (
                r_areal[idx] <= r_areal[idx - 1]
                and r_areal[idx] <= r_areal[idx + 1]
                and r_areal[idx] < 0.95 * r_asymp
            ):
                pinch = max(0.0, 1.0 - float(r_areal[idx]) / r_asymp)
                f_throat = max(f_throat, pinch)
    throat_pinch = float(np.clip(f_throat, 0.0, 1.0))

    d_proper = float(np.trapezoid(1.0 / np.sqrt(chi_ax), x))
    d_coord = 2.0 * integration_L
    portal_compression = float(np.clip((d_coord - d_proper) / d_coord, 0.0, 1.0))

    return MechanismScores(
        shift_warp=shift_warp,
        throat_pinch=throat_pinch,
        portal_compression=portal_compression,
    )


def compute_solved_geometry_ftl(
    gridinit_path: str | Path,
    *,
    n: int = 129,
    L: float | None = None,
) -> SolvedGeometryFtl | None:
    """Operational FTL + mechanism scores on constraint-satisfying initial data."""
    path = Path(gridinit_path)
    if not path.is_file():
        return None
    try:
        grid = read_gridinit(path)
        alpha, beta, gamma, spacing = build_xz_slice_from_gridinit(grid, n=n, L=L)
    except (OSError, ValueError, KeyError):
        return None

    integration_L = L if L is not None else float(
        0.49 * min(grid.data.shape[2] * grid.dx_xyz[0], grid.data.shape[0] * grid.dx_xyz[2])
    )
    mid = n // 2
    operational = operational_ftl_on_grid(
        alpha,
        beta,
        gamma,
        spacing=spacing,
        source=(1, mid),
        target=(n - 2, mid),
    )
    mechanisms = _mechanism_scores(
        alpha, beta, gamma, spacing, integration_L=integration_L
    )
    return SolvedGeometryFtl(
        operational=operational,
        mechanisms=mechanisms,
        integration_L=integration_L,
    )

