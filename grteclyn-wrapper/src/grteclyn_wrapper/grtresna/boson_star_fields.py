"""Paint boson-star complex scalar fields onto a gridinit array."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

from .boson_star_profile import BosonStarProfile
from .lump_fields import _coordinate_grids, angular_factor

# GRTeclyn RadialRecipe state layout (must match StateVariables.hpp).
_NUM_CCZ4_VARS = 25
# Canonical complex field Phi+ (slots 25-28) and phantom field Phi- (29-32).
_C_PHI_P = _NUM_CCZ4_VARS + 0
_C_PI_P = _NUM_CCZ4_VARS + 1
_C_PHI2_P = _NUM_CCZ4_VARS + 2
_C_PI2_P = _NUM_CCZ4_VARS + 3
_C_PHI_M = _NUM_CCZ4_VARS + 4
_C_PI_M = _NUM_CCZ4_VARS + 5
_C_PHI2_M = _NUM_CCZ4_VARS + 6
_C_PI2_M = _NUM_CCZ4_VARS + 7
_BICOMPLEX_NCOMP = _NUM_CCZ4_VARS + 8

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


def _raw_lump_phi_grid(
    lump: Mapping[str, Any], px: NDArray, py: NDArray, pz: NDArray,
) -> NDArray:
    """phi of one lump using the RAW amplitude (no exotic 0.25 scaling).

    Matches GRTresna's ``BosonStarParams::lump_phi1`` so the painted Phi+/Phi-
    fields agree with the per-lump-signed constraint solve.  The exotic sign is
    carried by the phantom FIELD (Phi-), not by an amplitude rescale.
    """
    amp = float(lump.get("amp", 0.0))
    if amp == 0.0:
        return np.zeros_like(px, dtype=np.float64)
    width = float(lump.get("width", 5.0))
    center = np.asarray(lump.get("center", (0.0, 0.0, 0.0)), dtype=np.float64)
    dx = px - center[0]
    dy = py - center[1]
    dz = pz - center[2]
    r2 = dx * dx + dy * dy + dz * dz
    if int(lump.get("profile", 0)) == 1:
        r = np.sqrt(r2)
        soft = 0.25 * width
        env = 0.5 * (1.0 - np.tanh((r - width) / soft))
    else:
        env = np.exp(-r2 / (2.0 * width * width))
    return amp * angular_factor(int(lump.get("mode", 0)), dx, dy, width) * env


def paint_bicomplex_fields_on_grid(
    data: NDArray,
    comp_names: list[str],
    dx: NDArray,
    origin: NDArray,
    lumps: Sequence[Mapping[str, Any]],
    *,
    omega: float,
) -> tuple[NDArray, list[str]]:
    """Paint TWO complex fields for the ``grtresna_bicomplex_scalar`` model.

    Canonical lumps (``exotic`` == 0) are superposed into the CANONICAL field
    Phi+ (state slots 25-28); exotic lumps into the PHANTOM field Phi- (slots
    29-32).  Each field is a boson star: ``phi_im = 0`` and
    ``Pi_im = -omega * phi_re / lapse`` at t=0 (the U(1) phase velocity), with
    ``Pi_re`` from any bulk kinematics (zero for at-rest lumps).  The metric
    (slots 0-24, incl. the solved lapse) is taken from GRTresna unchanged.
    """
    nz, ny, nx, ncomp = data.shape
    # Ensure room for the full canonical+phantom matter block.
    if ncomp < _BICOMPLEX_NCOMP:
        pad = np.zeros((nz, ny, nx, _BICOMPLEX_NCOMP - ncomp), dtype=np.float64)
        data = np.concatenate([data, pad], axis=3)
        ncomp = _BICOMPLEX_NCOMP

    out = np.array(data, dtype=np.float64, copy=True)

    name_to_idx = {name: i for i, name in enumerate(comp_names)}
    lapse_idx = name_to_idx.get("lapse", None)
    lapse = (
        out[:, :, :, lapse_idx]
        if lapse_idx is not None
        else np.ones((nz, ny, nx), dtype=np.float64)
    )
    inv_lapse = 1.0 / np.maximum(lapse, 1.0e-12)

    px_grid, py_grid, pz_grid = _coordinate_grids(nx, ny, nz, dx, origin)

    def _group_phi(exotic_wanted: bool) -> NDArray:
        acc = np.zeros((nz, ny, nx), dtype=np.float64)
        for lump in lumps:
            is_exotic = bool(int(lump.get("exotic", 0)))
            if is_exotic != exotic_wanted:
                continue
            acc += _raw_lump_phi_grid(lump, px_grid, py_grid, pz_grid)
        return acc

    phi_re_p = _group_phi(False)
    phi_re_m = _group_phi(True)

    # Zero out the whole matter block first, then write the two fields.
    out[:, :, :, _NUM_CCZ4_VARS:_BICOMPLEX_NCOMP] = 0.0
    # Canonical Phi+ (Pi_re from kinematics is ~0 for at-rest lumps).
    out[:, :, :, _C_PHI_P] = phi_re_p
    out[:, :, :, _C_PHI2_P] = 0.0
    out[:, :, :, _C_PI_P] = 0.0
    out[:, :, :, _C_PI2_P] = -omega * phi_re_p * inv_lapse
    # Phantom Phi-.
    out[:, :, :, _C_PHI_M] = phi_re_m
    out[:, :, :, _C_PHI2_M] = 0.0
    out[:, :, :, _C_PI_M] = 0.0
    out[:, :, :, _C_PI2_M] = -omega * phi_re_m * inv_lapse

    names = list(comp_names)
    matter_names = [
        "phi", "Pi", "phi2", "Pi2",
        "phi_lump1", "Pi_lump1", "phi_lump2", "Pi_lump2",
    ]
    for slot, nm in zip(range(_NUM_CCZ4_VARS, _BICOMPLEX_NCOMP), matter_names):
        if slot < len(names):
            names[slot] = nm
        else:
            names.append(nm)
    return out, names
