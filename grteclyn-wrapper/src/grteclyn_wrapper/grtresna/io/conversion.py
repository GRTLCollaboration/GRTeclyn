"""Chombo HDF5 to .gridinit conversion pipeline."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .chombo import chombo_to_uniform, read_chombo_domain
from .gridinit import write_gridinit

logger = logging.getLogger(__name__)

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
    source_origin: Sequence[float] | tuple[float, float, float] | None = None,
    delete_source: bool = False,
    lumps: Sequence[Mapping[str, Any]] | None = None,
    shift_seed: float = 0.0,
    num_workers: int = 0,
    matter_model: str = "",
    bs_omega: float = 0.0,
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
        chombo_path,
        nx=nx,
        ny=ny,
        nz=nz,
        Lx=Lx,
        Ly=Ly,
        Lz=Lz,
        source_origin=source_origin,
        num_workers=num_workers,
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

    is_bicomplex = matter_model == "grtresna_bicomplex_scalar"

    if is_bicomplex:
        # Two-complex-field (canonical + phantom) model.  Repaint the matter
        # block analytically from the lumps grouped by sign; the metric (incl.
        # the solved lapse) comes from GRTresna unchanged.
        from ..fields.boson_star import (
            paint_bicomplex_fields_on_grid,
            rename_complex_scalar_components,
        )
        from ..fields.lump import (
            paint_shift_seed_on_grid,
            shift_lump_centers_for_gridinit,
        )

        comp_names = rename_complex_scalar_components(comp_names)
        shifted_lumps = (
            shift_lump_centers_for_gridinit(
                lumps,
                grid_center=(
                    float(origin[0]) + 0.5 * float(Lx),
                    float(origin[1]) + 0.5 * float(Ly),
                    float(origin[2]) + 0.5 * float(Lz),
                ),
            )
            if lumps
            else []
        )
        if shift_seed and shifted_lumps:
            data = paint_shift_seed_on_grid(
                data, comp_names, dx_xyz, origin, shifted_lumps,
                shift_seed=float(shift_seed),
            )
        data, comp_names = paint_bicomplex_fields_on_grid(
            data, comp_names, dx_xyz, origin, shifted_lumps, omega=bs_omega,
        )
        result = write_gridinit(data, comp_names, dx_xyz, origin, output_path)
        if delete_source:
            src = Path(chombo_path)
            if src.exists():
                src.unlink()
                logger.info("Deleted source HDF5: %s", src)
        return result

    if lumps:
        from ..fields.lump import (
            paint_lump_fields_on_grid,
            paint_shift_seed_on_grid,
            shift_lump_centers_for_gridinit,
        )

        shifted_lumps = shift_lump_centers_for_gridinit(
            lumps,
            grid_center=(
                float(origin[0]) + 0.5 * float(Lx),
                float(origin[1]) + 0.5 * float(Ly),
                float(origin[2]) + 0.5 * float(Lz),
            ),
        )
        # Seed the initial shift before appending the diagnostic phi_lump
        # channels so the shift columns line up with GRTeclyn's state vars.
        if shift_seed:
            data = paint_shift_seed_on_grid(
                data, comp_names, dx_xyz, origin, shifted_lumps,
                shift_seed=float(shift_seed),
            )
        data, comp_names = paint_lump_fields_on_grid(
            data, comp_names, dx_xyz, origin, shifted_lumps,
        )

    complex_scalar_names = ("phi_re", "phi_im", "Pi_re", "Pi_im", "phi2", "Pi2")
    if any(name in comp_names for name in complex_scalar_names):
        from ..fields.boson_star import (
            apply_post_solve_lapse_correction,
            rename_complex_scalar_components,
        )

        # GRTresna emits phi/Pi/phi2/Pi2 already; the rename is a no-op for that
        # layout but also accepts the phi_re/phi_im legacy naming.
        comp_names = rename_complex_scalar_components(comp_names)
        data = apply_post_solve_lapse_correction(data, comp_names)

    result = write_gridinit(data, comp_names, dx_xyz, origin, output_path)

    if delete_source:
        src = Path(chombo_path)
        if src.exists():
            src.unlink()
            logger.info("Deleted source HDF5: %s", src)

    return result
