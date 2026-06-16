"""Paint boson-star complex scalar fields onto a gridinit array."""

from __future__ import annotations

from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from .boson_star_profile import BosonStarProfile
from .lump_fields import _coordinate_grids

# GRTresna HDF5 names -> GRTeclyn RadialRecipe names
_COMPLEX_SCALAR_RENAMES = {
    "phi_re": "phi",
    "Pi_re": "Pi",
    "phi_im": "phi2",
    "Pi_im": "Pi2",
}


def rename_complex_scalar_components(comp_names: list[str]) -> list[str]:
    return [_COMPLEX_SCALAR_RENAMES.get(name, name) for name in comp_names]


def apply_post_solve_lapse_correction(
    data: NDArray, comp_names: Sequence[str]
) -> NDArray:
    """Divide the painted Pi2 by the solved lapse.

    GRTresna paints the imaginary momentum as ``Pi_im = -omega*phi0`` assuming
    ``alpha = 1``.  A stationary boson star needs ``Pi_im = -omega*phi0/alpha``,
    so once the constraint solve has fixed the lapse the evolved initial data
    must rescale ``Pi2 /= alpha``.  Component names must already be renamed to
    the GRTeclyn convention (``Pi2``, ``lapse``).
    """
    idx = {name: i for i, name in enumerate(comp_names)}
    if "Pi2" not in idx or "lapse" not in idx:
        return data
    out = np.array(data, dtype=np.float64, copy=True)
    lapse = out[:, :, :, idx["lapse"]]
    with np.errstate(divide="ignore", invalid="ignore"):
        out[:, :, :, idx["Pi2"]] = out[:, :, :, idx["Pi2"]] / np.maximum(
            lapse, 1.0e-12
        )
    return out


def paint_boson_star_fields_on_grid(
    data: NDArray,
    comp_names: list[str],
    dx: NDArray,
    origin: NDArray,
    profile: BosonStarProfile,
    *,
    center: tuple[float, float, float] = (0.0, 0.0, 0.0),
    apply_post_solve_pi2: bool = True,
) -> tuple[NDArray, list[str]]:
    """Paint phi, phi2, Pi, Pi2 from a radial BS profile onto the grid."""
    names = rename_complex_scalar_components(list(comp_names))
    name_to_idx = {name: i for i, name in enumerate(names)}

    required = ("phi", "Pi", "phi2", "Pi2")
    for key in required:
        if key not in name_to_idx:
            raise ValueError(f"Missing component {key!r} in gridinit names: {names}")

    nz, ny, nx, _ = data.shape
    px_grid, py_grid, pz_grid = _coordinate_grids(nx, ny, nz, dx, origin)
    cx, cy, cz = center
    rr = np.sqrt(
        (px_grid - cx) ** 2 + (py_grid - cy) ** 2 + (pz_grid - cz) ** 2
    )

    phi1 = np.asarray(profile.eval_phi0(rr), dtype=np.float64)
    phi2 = np.zeros_like(phi1)
    pi1 = np.zeros_like(phi1)

    if apply_post_solve_pi2 and "lapse" in name_to_idx:
        lapse = data[:, :, :, name_to_idx["lapse"]]
        with np.errstate(divide="ignore", invalid="ignore"):
            pi2 = -profile.omega * phi1 / np.maximum(lapse, 1.0e-12)
    else:
        pi2 = -profile.omega * phi1

    painted = np.array(data, dtype=np.float64, copy=True)
    painted[:, :, :, name_to_idx["phi"]] = phi1
    painted[:, :, :, name_to_idx["phi2"]] = phi2
    painted[:, :, :, name_to_idx["Pi"]] = pi1
    painted[:, :, :, name_to_idx["Pi2"]] = pi2

    return painted, names
