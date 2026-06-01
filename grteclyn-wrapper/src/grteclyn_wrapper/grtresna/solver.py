"""Orchestrate GRTresna to produce constraint-satisfying initial data.

Given a parameter vector theta (BH masses, spins, scalar field config),
this module:
1. Writes a GRTresna params.txt
2. Runs GRTresna via MPI
3. Converts the Chombo HDF5 output to a .gridinit file for GRTeclyn
4. Cleans up heavy intermediate files to keep disk alive in a conveyor

The resulting .gridinit can be loaded by GRTeclyn's
ExternalGridInitialData to start an evolution from a fully
constraint-satisfied initial state.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

from ..core.config import REPO_ROOT
from .io import convert_chombo_to_gridinit

logger = logging.getLogger(__name__)

_DEFAULT_GRTRESNA_ROOT = Path(
    os.environ.get("GRTRESNA_ROOT", str(REPO_ROOT.parent / "GRTresna"))
)
_DEFAULT_EXAMPLE = "ScalarFieldBH"


@dataclass
class GRTresnaConfig:
    """All knobs for a GRTresna solve."""

    grtresna_root: Path = _DEFAULT_GRTRESNA_ROOT
    example: str = _DEFAULT_EXAMPLE
    mpi_ranks: int = 8
    max_NL_iterations: int = 50
    write_diagnostics: bool = False
    timeout: int = 3600

    # Grid — GRTresna uses half-z by default (reflective lo_boundary z=1)
    N: tuple[int, int, int] = (64, 64, 32)
    L: float = 128.0
    max_level: int = 3
    block_factor: int = 16
    max_grid_size: int = 16

    # BH params
    bh1_bare_mass: float = 1.0
    bh1_spin: tuple[float, float, float] = (0.0, 0.0, 0.5)
    bh1_momentum: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bh1_offset: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bh2_bare_mass: float = 0.0
    bh2_spin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bh2_momentum: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bh2_offset: tuple[float, float, float] = (0.0, 0.0, 0.0)

    # Scalar field
    G_Newton: float = 1.0
    phi_0: float = 0.0
    dphi: float = 0.1
    dphi_length: float = 5.0
    pi_0: float = 0.0
    dpi: float = 0.1
    dpi_length: float = 5.0
    scalar_mass: float = 0.1

    # Conformal / K
    regularised_part_psi: float = 1.0
    sign_of_K: int = 1

    # Boundary conditions
    use_compact_Vi_ansatz: int = 1
    hi_boundary: tuple[int, int, int] = (0, 0, 0)
    lo_boundary: tuple[int, int, int] = (0, 0, 1)

    # Output target resolution for .gridinit (per-axis)
    gridinit_nx: int = 64
    gridinit_ny: int = 64
    gridinit_nz: int = 64

    # Cleanup: remove HDF5 / pout / intermediate files after conversion
    cleanup: bool = True
    # Remove the entire work_dir after producing the gridinit
    cleanup_workdir: bool = False


def _fmt(v: Any) -> str:
    """Format a value for GRTresna params.txt."""
    if isinstance(v, (tuple, list)):
        return " ".join(str(x) for x in v)
    return str(v)


def write_grtresna_params(cfg: GRTresnaConfig, path: Path) -> None:
    """Write a GRTresna params.txt."""
    lines = [
        f"output_path = Outputs/",
        f"output_filename = InitialDataFinal.3d.hdf5",
        f"error_filename = Ham_and_Mom_errors",
        f"write_diagnostics = {1 if cfg.write_diagnostics else 0}",
        f"",
        f"N = {_fmt(cfg.N)}",
        f"L = {cfg.L}",
        f"max_level = {cfg.max_level}",
        f"block_factor = {cfg.block_factor}",
        f"max_grid_size = {cfg.max_grid_size}",
        f"regrid_radius = 10",
        f"",
        f"is_periodic = 0 0 0",
        f"use_compact_Vi_ansatz = {cfg.use_compact_Vi_ansatz}",
        f"hi_boundary = {_fmt(cfg.hi_boundary)}",
        f"lo_boundary = {_fmt(cfg.lo_boundary)}",
        f"",
        f"G_Newton = {cfg.G_Newton}",
        f"phi_0 = {cfg.phi_0}",
        f"dphi = {cfg.dphi}",
        f"dphi_length = {cfg.dphi_length}",
        f"pi_0 = {cfg.pi_0}",
        f"dpi = {cfg.dpi}",
        f"dpi_length = {cfg.dpi_length}",
        f"scalar_mass = {cfg.scalar_mass}",
        f"regularised_part_psi = {cfg.regularised_part_psi}",
        f"sign_of_K = {cfg.sign_of_K}",
        f"",
        f"bh1_bare_mass = {cfg.bh1_bare_mass}",
        f"bh1_spin = {_fmt(cfg.bh1_spin)}",
        f"bh1_momentum = {_fmt(cfg.bh1_momentum)}",
        f"bh1_offset = {_fmt(cfg.bh1_offset)}",
        f"bh2_bare_mass = {cfg.bh2_bare_mass}",
        f"bh2_spin = {_fmt(cfg.bh2_spin)}",
        f"bh2_momentum = {_fmt(cfg.bh2_momentum)}",
        f"bh2_offset = {_fmt(cfg.bh2_offset)}",
        f"",
        f"max_NL_iterations = {cfg.max_NL_iterations}",
        f"deactivate_zero_mode = 0",
    ]
    path.write_text("\n".join(lines) + "\n")


def _find_executable(cfg: GRTresnaConfig) -> Path:
    """Locate the compiled GRTresna executable."""
    example_dir = cfg.grtresna_root / "Examples" / cfg.example
    candidates = list(example_dir.glob("Main_*3d.*.ex"))
    if not candidates:
        raise FileNotFoundError(
            f"No GRTresna executable found in {example_dir}. "
            "Build it first with `make all`."
        )
    return candidates[0]


def _cleanup_workdir(work_dir: Path, keep_gridinit: Path | None) -> None:
    """Remove heavy files from the work directory."""
    outputs_dir = work_dir / "Outputs"
    removed_bytes = 0

    # Remove all HDF5 files (the big ones: 300-400 MB each)
    for hdf5 in outputs_dir.glob("*.hdf5"):
        sz = hdf5.stat().st_size
        hdf5.unlink()
        removed_bytes += sz

    # Remove NL iteration diagnostics
    for diag in outputs_dir.glob("NL_iteration_*.hdf5"):
        if diag.exists():
            removed_bytes += diag.stat().st_size
            diag.unlink()

    # Remove pout (MPI rank logs, can be large with many ranks)
    pout_dir = work_dir / "pout"
    if pout_dir.is_dir():
        for pf in pout_dir.iterdir():
            if pf.is_file():
                removed_bytes += pf.stat().st_size
                pf.unlink()

    if removed_bytes > 0:
        logger.info(
            "Cleanup: removed %.1f MB from %s",
            removed_bytes / 1e6,
            work_dir,
        )


def parse_convergence(work_dir: Path) -> dict[str, float] | None:
    """Read the final Ham/Mom errors from the GRTresna error file."""
    err_file = work_dir / "Ham_and_Mom_errors.txt"
    if not err_file.exists():
        return None
    lines = err_file.read_text().strip().splitlines()
    if len(lines) < 2:
        return None
    last = lines[-1].split()
    if len(last) < 3:
        return None
    try:
        return {
            "iteration": int(last[0]),
            "ham_pct": float(last[1]),
            "mom_pct": float(last[2]),
        }
    except (ValueError, IndexError):
        return None


def solve(
    cfg: GRTresnaConfig,
    work_dir: Path | None = None,
    gridinit_path: Path | None = None,
) -> Path:
    """Run GRTresna and return the path to the .gridinit file.

    Parameters
    ----------
    cfg : solver configuration
    work_dir : directory to run in (created if None)
    gridinit_path : where to write the .gridinit (defaults to work_dir)

    Returns
    -------
    Path to the .gridinit file
    """
    if work_dir is None:
        work_dir = Path(tempfile.mkdtemp(prefix="grtresna_"))
    work_dir.mkdir(parents=True, exist_ok=True)

    exe = _find_executable(cfg)

    outputs_dir = work_dir / "Outputs"
    outputs_dir.mkdir(exist_ok=True)
    (work_dir / "pout").mkdir(exist_ok=True)

    params_path = work_dir / "params.txt"
    write_grtresna_params(cfg, params_path)

    cmd_parts = ["mpirun"]
    if cfg.mpi_ranks > 1:
        cmd_parts.append("--oversubscribe")
    cmd_parts.extend(["-np", str(cfg.mpi_ranks), str(exe), str(params_path)])

    logger.info(
        "Running GRTresna: %d ranks, max %d iterations, grid %s, L=%.0f",
        cfg.mpi_ranks, cfg.max_NL_iterations, cfg.N, cfg.L,
    )

    result = subprocess.run(
        cmd_parts,
        cwd=str(work_dir),
        capture_output=True,
        text=True,
        timeout=cfg.timeout,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"GRTresna failed (rc={result.returncode}):\n"
            f"stdout: {result.stdout[-2000:]}\n"
            f"stderr: {result.stderr[-2000:]}"
        )

    convergence = parse_convergence(work_dir)
    if convergence:
        logger.info(
            "GRTresna converged: iter=%d, Ham=%.4f%%, Mom=%.4f%%",
            convergence["iteration"], convergence["ham_pct"], convergence["mom_pct"],
        )

    chombo_hdf5 = outputs_dir / "InitialDataFinal.3d.hdf5"
    if not chombo_hdf5.exists():
        raise FileNotFoundError(
            f"GRTresna did not produce {chombo_hdf5}. "
            f"Check {work_dir / 'Ham_and_Mom_errors.txt'}."
        )

    if gridinit_path is None:
        gridinit_path = work_dir / "initial_data.gridinit"

    convert_chombo_to_gridinit(
        chombo_hdf5,
        gridinit_path,
        nx=cfg.gridinit_nx,
        ny=cfg.gridinit_ny,
        nz=cfg.gridinit_nz,
        L=cfg.L,
        delete_source=cfg.cleanup,
    )

    if cfg.cleanup:
        _cleanup_workdir(work_dir, keep_gridinit=gridinit_path)

    if cfg.cleanup_workdir:
        # Move gridinit out before nuking the dir
        if gridinit_path.parent == work_dir:
            import tempfile as _tf
            fd, tmp = _tf.mkstemp(suffix=".gridinit")
            os.close(fd)
            shutil.move(str(gridinit_path), tmp)
            shutil.rmtree(work_dir, ignore_errors=True)
            shutil.move(tmp, str(gridinit_path))
        else:
            shutil.rmtree(work_dir, ignore_errors=True)

    return gridinit_path
