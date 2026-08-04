"""Run GRTresna via MPI and produce a .gridinit file."""

from __future__ import annotations

import logging
import os
import shutil
import signal
import subprocess
import tempfile
from pathlib import Path

from ..io.conversion import convert_chombo_to_gridinit
from ..matter.models import finalize_solver_config
from ..matter.sign_consistency import enforce_id_evolution_sign_consistency
from ..matter.wiring import write_matter_metadata
from .config import GRTresnaConfig
from .convergence import parse_convergence
from .params import write_grtresna_params

logger = logging.getLogger(__name__)


def _candidate_env_prefixes() -> list[Path]:
    """Prefixes that may contain ``bin/mpirun`` for GRTresna.

    Host-specific env roots are not guessed — set ``GRTRESNA_ENV`` in ``.env``.
    """

    prefixes: list[Path] = []
    explicit = os.environ.get("GRTRESNA_ENV")
    if explicit:
        prefixes.append(Path(explicit).expanduser())
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        prefixes.append(Path(conda_prefix))
    return prefixes


def _resolve_mpirun(cfg: GRTresnaConfig) -> tuple[str, dict[str, str]]:
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
    del keep_gridinit
    outputs_dir = work_dir / "Outputs"
    removed_bytes = 0

    for hdf5 in outputs_dir.glob("*.hdf5"):
        sz = hdf5.stat().st_size
        hdf5.unlink()
        removed_bytes += sz

    for diag in outputs_dir.glob("NL_iteration_*.hdf5"):
        if diag.exists():
            removed_bytes += diag.stat().st_size
            diag.unlink()

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


def solve(
    cfg: GRTresnaConfig,
    work_dir: Path | None = None,
    gridinit_path: Path | None = None,
) -> Path:
    """Run GRTresna and return the path to the .gridinit file."""
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
    # Refuse to solve when ID exotic signs disagree with the evolution model
    # (single-complex + per-lump phantom is the eval-118 failure mode).
    enforce_id_evolution_sign_consistency(cfg)
    exe = _find_executable(cfg)

    outputs_dir = work_dir / "Outputs"
    outputs_dir.mkdir(exist_ok=True)
    (work_dir / "pout").mkdir(exist_ok=True)

    params_path = work_dir / "params.txt"
    write_grtresna_params(cfg, params_path)

    run_env = dict(os.environ)
    if cfg.mpi_ranks <= 1:
        # One rank never needs a launcher.  An MPI-built executable started
        # directly enters MPI singleton mode, which is why the old ".MPI." name
        # test is gone: on a node whose PRRTE is unusable (mpirun segfaults at
        # DVM start-up) the only GRTresna binary present is the MPI one, and
        # routing it through mpirun turned every candidate into a bogus
        # grtresna_failed rejection.  Verified against an 8-rank solve of the
        # same params: Ham 0.63% / Mom 0.65% singleton vs 0.93% / 0.74% MPI.
        cmd_parts = [str(exe), str(params_path)]
        launch_desc = "singleton" if ".MPI." in exe.name else "serial"
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

    if cfg.lumps or cfg.matter_model in (
        "grtresna_complex_scalar",
        "grtresna_bicomplex_scalar",
    ):
        write_matter_metadata(gridinit_path.with_suffix(".matter.json"), cfg)

    if cfg.cleanup:
        _cleanup_workdir(work_dir, keep_gridinit=gridinit_path)

    if cfg.cleanup_workdir:
        if gridinit_path.parent == work_dir:
            fd, tmp = tempfile.mkstemp(suffix=".gridinit")
            os.close(fd)
            shutil.move(str(gridinit_path), tmp)
            shutil.rmtree(work_dir, ignore_errors=True)
            shutil.move(tmp, str(gridinit_path))
        else:
            shutil.rmtree(work_dir, ignore_errors=True)

    return gridinit_path
