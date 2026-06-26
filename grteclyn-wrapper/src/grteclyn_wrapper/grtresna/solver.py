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
import math
import os
import signal
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

from ..core.config import REPO_ROOT
from .io import convert_chombo_to_gridinit
from .gridinit_export import GridinitExportSpec
from .matter_models import finalize_solver_config
from .matter_wiring import write_matter_metadata

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
    # Adaptive early-exit for the outer nonlinear solve. nl_stall_tolerance stops
    # the solve once the per-iteration relative improvement of both Ham and Mom
    # errors drops below it (residual at its floor) without losing solution
    # quality; nl_exit_tolerance stops once both relative errors (in %) are below
    # it. 0 disables a check.
    nl_exit_tolerance: float = 1.0
    nl_stall_tolerance: float = 0.02
    write_diagnostics: bool = False
    timeout: int = 3600

    # Path to the ``mpirun`` that matches the GRTresna build environment. When
    # None it is auto-resolved (see ``_resolve_mpirun``) so ``--grtresna`` works
    # without the caller having to prefix PATH manually.
    mpirun: str | None = None

    # Grid — GRTresna uses half-z by default (reflective lo_boundary z=1)
    N: tuple[int, int, int] = (64, 64, 32)
    L: float = 128.0
    max_level: int = 3
    block_factor: int = 16
    max_grid_size: int = 16
    refine_threshold: float = 0.5
    regrid_radius: float = 10.0
    coefficient_average_type: str = "harmonic"

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
    scalar_lambda: float = 0.0
    scalar_mu: float = 0.0

    # Matter model selection (empty => legacy independent real scalars via ScalarFieldBH)
    matter_model: str = ""

    # Boson-star parameters (BosonStarBH example)
    bs_phi_c: float = 0.08
    bs_profile_width: float = 8.0
    bs_omega: float = 0.0
    # Phantom sign for the complex scalar: +1 canonical (WEC-satisfying),
    # -1 phantom (flips entire T_ab, giving negative energy density while
    # preserving U(1) charge conservation).
    scalar_sign: int = 1

    # Momentum-carrying scalar "cloud" (a single localised lump). All default
    # to off so a config that does not set them reproduces the legacy spherical
    # data. The lump's conjugate momentum is built so the configuration carries
    # net linear momentum (lump_velocity) and/or angular momentum (lump_omega
    # with lump_mode >= 1); GRTresna then solves the momentum constraint for it.
    # NOTE: this single-lump API is kept for convenience; for a richer matter
    # distribution set ``lumps`` (a list of dicts) instead -- when non-empty it
    # takes precedence and a num_lumps/indexed block is written.
    lump_amp: float = 0.0
    lump_width: float = 5.0
    lump_center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    lump_velocity: tuple[float, float, float] = (0.0, 0.0, 0.0)
    lump_omega: float = 0.0
    lump_mode: int = 0
    # Phantom/ghost flag for the legacy single lump (0 => canonical positive
    # energy, 1 => exotic negative-energy source). The multi-lump ``lumps``
    # basis carries its own per-entry ``exotic`` key.
    lump_exotic: int = 0

    # Multi-lump scalar basis. Each entry is a dict with keys:
    #   amp, width, center (3-tuple), velocity (3-tuple), omega, mode.
    # A superposition of lumps paints an arbitrary energy/momentum distribution
    # which GRTresna turns into constraint-satisfying initial data. When this is
    # non-empty it overrides the single-lump fields above.
    lumps: list[dict] = field(default_factory=list)

    # Seed an initial coordinate shift in the exported gridinit, aligned with the
    # matter momentum density (beta^i = shift_seed * S_i / max|S|, peak |beta| =
    # |shift_seed|). GRTresna emits zero initial shift, so this is what makes the
    # warp/channel mechanism (shift_drive, channel_progress) reachable from t=0.
    # 0 disables the seed (GRTresna's native zero-shift gauge).
    shift_seed: float = 0.0

    # Conformal / K
    regularised_part_psi: float = 1.0
    sign_of_K: int = 1

    # Slicing + nonlinear-solve robustness. The default CTTKHybrid ansatz
    # K = sign*sqrt(24 pi G rho + ...) is imaginary wherever rho < 0, so it
    # cannot solve EXOTIC (negative-energy) initial data. ``maximal_slicing``
    # switches GRTresna to a K = 0 York/Lichnerowicz solve that sources the
    # matter energy elliptically and handles either sign of rho. For the
    # indefinite operators that strong exotic matter produces, ``psi_relaxation``
    # under-relaxes the Newton step and ``psi_floor`` clamps psi > 0 so the
    # solve does not diverge to NaN. Defaults below reproduce the original
    # (canonical-only) behaviour; the search enables the exotic-safe path.
    maximal_slicing: bool = False
    psi_relaxation: float = 1.0
    psi_floor: float = -1.0
    maximal_jacobian_cap: float = -1.0

    # Boundary conditions
    use_compact_Vi_ansatz: int = 1
    hi_boundary: tuple[int, int, int] = (0, 0, 0)
    lo_boundary: tuple[int, int, int] = (0, 0, 1)

    # Output target resolution for .gridinit (per-axis)
    gridinit_nx: int = 64
    gridinit_ny: int = 64
    gridinit_nz: int = 64
    # Parallel AMR box painting during Chombo→gridinit conversion (0 = auto).
    gridinit_workers: int = 0

    # Centre of the GRTeclyn evolution box that will LOAD this gridinit. The
    # GRTresna matter sits at the centre of the (large) solve domain; the
    # gridinit ``origin`` is written so that this matter centre maps onto
    # ``target_center`` in GRTeclyn coordinates (GRTeclyn then evolves the
    # central window of the GRTresna domain at its own resolution). The default
    # matches the RadialRecipe template (L_full=64, half-z => centre (32,32,0)).
    target_center: tuple[float, float, float] = (32.0, 32.0, 0.0)

    # When set (typically by ``GRTresnaDomainConfig.apply_to_solver``), controls
    # the physical export window and spacing so gridinit dx matches evolution dx.
    gridinit_export: GridinitExportSpec | None = None

    # Cleanup: remove HDF5 / pout / intermediate files after conversion
    cleanup: bool = True
    # Remove the entire work_dir after producing the gridinit
    cleanup_workdir: bool = False


def apply_exotic_safe_solver(cfg: GRTresnaConfig) -> GRTresnaConfig:
    """Switch ``cfg`` to the K=0 maximal-slicing York/Lichnerowicz path in place.

    The default CTTKHybrid ansatz K = sign*sqrt(24 pi G rho + ...) is imaginary
    wherever rho < 0, so any configuration containing an EXOTIC (negative-energy)
    lump produces NaN residuals unless it is solved on the maximal-slicing path
    (matter sourced elliptically, with Newton under-relaxation and a psi-floor
    for the indefinite operator). Purely canonical configs must NOT use this --
    the standard ansatz is markedly more robust for positive multi-lump matter.

    Only fills values the caller left at their canonical defaults, so explicit
    overrides win. This is the single source of truth shared by the search
    driver (``build_grtresna_config``) and standalone solves/smokes.
    """
    if not cfg.maximal_slicing:
        cfg.maximal_slicing = True
    if cfg.psi_relaxation == 1.0:
        cfg.psi_relaxation = 0.6
    if cfg.psi_floor <= 0.0:
        cfg.psi_floor = 0.1
    if cfg.maximal_jacobian_cap <= 0.0:
        cfg.maximal_jacobian_cap = 25.0
    if cfg.coefficient_average_type == "harmonic":
        cfg.coefficient_average_type = "arithmetic"
    return cfg


def config_has_exotic_lump(cfg: GRTresnaConfig) -> bool:
    """True when any configured matter is flagged exotic (negative energy).

    Covers independent-scalar exotic lumps and phantom complex scalars.
    """
    if cfg.matter_model == "grtresna_complex_scalar" and cfg.scalar_sign < 0:
        return True
    if cfg.lumps:
        return any(int(lump.get("exotic", 0)) for lump in cfg.lumps)
    return bool(cfg.lump_exotic)


def _fmt(v: Any) -> str:
    """Format a value for GRTresna params.txt."""
    if isinstance(v, (tuple, list)):
        return " ".join(str(x) for x in v)
    return str(v)


def _lump_lines(cfg: GRTresnaConfig) -> list[str]:
    """Render the momentum-carrying matter block.

    If ``cfg.lumps`` is non-empty, write a ``num_lumps`` + indexed
    ``lump{k}_*`` basis; otherwise write the legacy single-lump keys.
    """
    if cfg.lumps:
        lines = [f"num_lumps = {len(cfg.lumps)}"]
        for k, lump in enumerate(cfg.lumps):
            lines.extend([
                f"lump{k}_amp = {lump.get('amp', 0.0)}",
                f"lump{k}_width = {lump.get('width', 5.0)}",
                f"lump{k}_center = {_fmt(tuple(lump.get('center', (0.0, 0.0, 0.0))))}",
                f"lump{k}_velocity = {_fmt(tuple(lump.get('velocity', (0.0, 0.0, 0.0))))}",
                f"lump{k}_omega = {lump.get('omega', 0.0)}",
                f"lump{k}_mode = {int(lump.get('mode', 0))}",
                f"lump{k}_exotic = {int(lump.get('exotic', 0))}",
                f"lump{k}_profile = {int(lump.get('profile', 0))}",
            ])
        return lines
    return [
        f"lump_amp = {cfg.lump_amp}",
        f"lump_width = {cfg.lump_width}",
        f"lump_center = {_fmt(cfg.lump_center)}",
        f"lump_velocity = {_fmt(cfg.lump_velocity)}",
        f"lump_omega = {cfg.lump_omega}",
        f"lump_mode = {cfg.lump_mode}",
        f"lump_exotic = {int(cfg.lump_exotic)}",
    ]


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
        f"refine_threshold = {cfg.refine_threshold}",
        f"block_factor = {cfg.block_factor}",
        f"max_grid_size = {cfg.max_grid_size}",
        f"regrid_radius = {cfg.regrid_radius}",
        f"coefficient_average_type = {cfg.coefficient_average_type}",
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
        f"scalar_lambda = {cfg.scalar_lambda}",
    ]
    if cfg.matter_model in (
        "grtresna_complex_scalar",
        "grtresna_bicomplex_scalar",
    ):
        lines.extend([
            f"scalar_sign = {cfg.scalar_sign}",
            f"bs_phi_c = {cfg.bs_phi_c}",
            f"bs_profile_width = {cfg.bs_profile_width}",
            f"bs_omega = {cfg.bs_omega}",
        ])
    lines.extend(_lump_lines(cfg))
    lines.extend([
        f"regularised_part_psi = {cfg.regularised_part_psi}",
        f"sign_of_K = {cfg.sign_of_K}",
        f"maximal_slicing = {1 if cfg.maximal_slicing else 0}",
        f"psi_relaxation = {cfg.psi_relaxation}",
        f"psi_floor = {cfg.psi_floor}",
        f"maximal_jacobian_cap = {cfg.maximal_jacobian_cap}",
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
        f"NL_exit_tolerance = {cfg.nl_exit_tolerance}",
        f"NL_stall_tolerance = {cfg.nl_stall_tolerance}",
        f"deactivate_zero_mode = 0",
    ])
    path.write_text("\n".join(lines) + "\n")


def _candidate_env_prefixes() -> list[Path]:
    """Ordered list of conda/micromamba env prefixes that may hold mpirun."""
    name = os.environ.get("GRTRESNA_ENV_NAME", "grtresna")
    prefixes: list[Path] = []
    explicit = os.environ.get("GRTRESNA_ENV")
    if explicit:
        prefixes.append(Path(explicit))
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        prefixes.append(Path(conda_prefix))
    home = Path.home()
    for root in (
        home / ".mlspace" / "envs",
        home / "micromamba" / "envs",
        home / "miniconda3" / "envs",
        home / "miniforge3" / "envs",
        home / "anaconda3" / "envs",
        Path("/opt/conda/envs"),
    ):
        prefixes.append(root / name)
    return prefixes


def _resolve_mpirun(cfg: GRTresnaConfig) -> tuple[str, dict[str, str]]:
    """Resolve an ``mpirun`` and the environment to run GRTresna under.

    Resolution order:
      1. ``cfg.mpirun`` (explicit override)
      2. ``$GRTRESNA_MPIRUN`` (explicit override)
      3. ``bin/mpirun`` inside a candidate conda/micromamba env prefix
         (``$GRTRESNA_ENV``, ``$CONDA_PREFIX``, common env roots for
         ``$GRTRESNA_ENV_NAME`` which defaults to "grtresna")
      4. ``mpirun`` already on PATH

    Returns the mpirun path and a subprocess env with that mpirun's directory
    prepended to PATH (so its helper binaries, e.g. prted/orted, resolve).
    """
    candidates: list[str] = []
    if cfg.mpirun:
        candidates.append(cfg.mpirun)
    env_override = os.environ.get("GRTRESNA_MPIRUN")
    if env_override:
        candidates.append(env_override)
    for prefix in _candidate_env_prefixes():
        candidates.append(str(prefix / "bin" / "mpirun"))

    mpirun_path: str | None = None
    mpirun_prefix: Path | None = None
    for cand in candidates:
        if cand and Path(cand).is_file():
            mpirun_path = cand
            cand_path = Path(cand).resolve()
            if cand_path.parent.name == "bin":
                mpirun_prefix = cand_path.parent.parent
            break
    if mpirun_path is None:
        which = shutil.which("mpirun")
        if which:
            mpirun_path = which
    if mpirun_path is None:
        raise FileNotFoundError(
            "Could not locate 'mpirun' for GRTresna. Set $GRTRESNA_MPIRUN to its "
            "absolute path, or $GRTRESNA_ENV / $GRTRESNA_ENV_NAME to the conda "
            "environment that built GRTresna, or add mpirun to PATH."
        )

    env = dict(os.environ)
    mpirun_dir = str(Path(mpirun_path).resolve().parent)
    env["PATH"] = mpirun_dir + os.pathsep + env.get("PATH", "")
    if mpirun_prefix is not None:
        env["CONDA_PREFIX"] = str(mpirun_prefix)
        lib_dir = mpirun_prefix / "lib"
        if lib_dir.exists():
            env["LD_LIBRARY_PATH"] = (
                str(lib_dir) + os.pathsep + env.get("LD_LIBRARY_PATH", "")
            )
    return mpirun_path, env


def _find_executable(cfg: GRTresnaConfig) -> Path:
    """Locate the compiled GRTresna executable."""
    example_dir = cfg.grtresna_root / "Examples" / cfg.example
    candidates = list(example_dir.glob("Main_*3d.*.ex"))
    if not candidates:
        raise FileNotFoundError(
            f"No GRTresna executable found in {example_dir}. "
            "Build it first with `make all`."
        )
    serial = sorted(path for path in candidates if ".MPI." not in path.name)
    mpi = sorted(path for path in candidates if ".MPI." in path.name)
    if cfg.mpi_ranks <= 1 and serial:
        return serial[0]
    if cfg.mpi_ranks > 1 and mpi:
        return mpi[0]
    return sorted(candidates)[0]


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
        ham_pct = float(last[1])
        mom_pct = float(last[2])
    except (ValueError, IndexError):
        return None
    # Static (momentum-free) matter has an identically zero momentum source, so
    # GRTresna's *relative* momentum error is 0/0 = NaN even when the solve is
    # perfectly converged.  A finite Hamiltonian residual proves the fields are
    # finite, so a NaN momentum error in that case means "no momentum to
    # violate" -> the momentum constraint is trivially satisfied (0%).  Without
    # this, every static-matter candidate is falsely rejected as "nonfinite
    # convergence" (a genuine divergence shows up as NaN in *both* Ham and Mom).
    if math.isfinite(ham_pct) and not math.isfinite(mom_pct):
        mom_pct = 0.0
    return {
        "iteration": int(last[0]),
        "ham_pct": ham_pct,
        "mom_pct": mom_pct,
    }


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
    caller_cwd = Path.cwd()
    if work_dir is None:
        work_dir = Path(tempfile.mkdtemp(prefix="grtresna_"))
    else:
        work_dir = Path(work_dir).expanduser()
        if not work_dir.is_absolute():
            work_dir = caller_cwd / work_dir
        work_dir = work_dir.resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    finalize_solver_config(cfg)
    exe = _find_executable(cfg)

    outputs_dir = work_dir / "Outputs"
    outputs_dir.mkdir(exist_ok=True)
    (work_dir / "pout").mkdir(exist_ok=True)

    params_path = work_dir / "params.txt"
    write_grtresna_params(cfg, params_path)

    run_env = dict(os.environ)
    if cfg.mpi_ranks <= 1 and ".MPI." not in exe.name:
        cmd_parts = [str(exe), str(params_path)]
        launch_desc = "serial"
    else:
        mpirun_path, run_env = _resolve_mpirun(cfg)
        cmd_parts = [mpirun_path]
        launch_desc = mpirun_path
        cmd_parts.append("--oversubscribe")
        cmd_parts.extend(["-np", str(cfg.mpi_ranks), str(exe), str(params_path)])

    logger.info(
        "Running GRTresna: %d ranks, max %d iterations, grid %s, L=%.0f (launcher=%s)",
        cfg.mpi_ranks, cfg.max_NL_iterations, cfg.N, cfg.L, launch_desc,
    )

    proc = subprocess.Popen(
        cmd_parts,
        cwd=str(work_dir),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=run_env,
        start_new_session=True,
    )
    try:
        stdout, stderr = proc.communicate(timeout=cfg.timeout)
    except subprocess.TimeoutExpired as exc:
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            stdout, stderr = proc.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            stdout, stderr = proc.communicate()
        raise RuntimeError(
            f"GRTresna timed out after {cfg.timeout}s:\n"
            f"stdout: {stdout[-2000:]}\n"
            f"stderr: {stderr[-2000:]}"
        ) from exc

    result = subprocess.CompletedProcess(cmd_parts, proc.returncode, stdout, stderr)
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
    else:
        gridinit_path = Path(gridinit_path).expanduser()
        if not gridinit_path.is_absolute():
            gridinit_path = caller_cwd / gridinit_path
        gridinit_path = gridinit_path.resolve()

    export = cfg.gridinit_export
    convert_kwargs: dict[str, object] = {
        "nx": cfg.gridinit_nx,
        "ny": cfg.gridinit_ny,
        "nz": cfg.gridinit_nz,
        "target_center": cfg.target_center,
        "delete_source": cfg.cleanup,
        "lumps": cfg.lumps or None,
        "shift_seed": cfg.shift_seed,
        "num_workers": cfg.gridinit_workers,
        "matter_model": cfg.matter_model,
        "bs_omega": cfg.bs_omega,
    }
    if export is not None:
        convert_kwargs.update(
            {
                "Lx": export.lx,
                "Ly": export.ly,
                "Lz": export.lz,
                "source_origin": export.source_origin,
            }
        )
    else:
        convert_kwargs["L"] = cfg.L

    convert_chombo_to_gridinit(
        chombo_hdf5,
        gridinit_path,
        **convert_kwargs,
    )

    if cfg.lumps or cfg.matter_model == "grtresna_complex_scalar":
        write_matter_metadata(gridinit_path.with_suffix(".matter.json"), cfg)

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
