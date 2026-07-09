"""CMA-ES iteration loop for geometry-first matter adjustment.

Given a GeometryMotif target, this module repeatedly adjusts scalar-lump
parameters (amp, width, center, velocity, omega) and invokes GRTresna to find
the matter configuration that best reproduces the target geometry while
satisfying the Hamiltonian and momentum constraints.

The loop is cheap (GRTresna-only, no GPU evolution) and fits between the
geometry-first scout and the existing QD/CMA-ES/HQ campaign stages.
"""

from __future__ import annotations

import dataclasses
import json
import logging
import shutil
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from ..grtresna.fit.motif import (
    FittedMatter,
    MAX_LUMP_AMP,
    MAX_LUMP_WIDTH,
    MAX_VELOCITY,
    MIN_LUMP_WIDTH,
    build_grtresna_config_from_fitted,
    fitted_matter_to_dict,
    write_fitted_matter_json,
)
from ..grtresna.solver import GRTresnaConfig, solve
from ..grtresna.solver.convergence import parse_convergence
from ..initial_data.motif import GeometryMotif, MomentumTarget
from ..projection.motif_preservation import compare_motif_preservation
from .mismatch import GATE_FITNESS, MismatchReport, compute_mismatch

logger = logging.getLogger(__name__)

# Per-lump parameter order (fixed lump count / exotic / mode flags).
_LUMP_PARAMS = ("amp", "width", "cx", "cy", "cz", "vx", "vy", "vz", "omega")
_PARAMS_PER_LUMP = len(_LUMP_PARAMS)

# Bounds for each per-lump parameter.
_LOWER = np.array([0.01, MIN_LUMP_WIDTH, -64.0, -64.0, -64.0, -MAX_VELOCITY, -MAX_VELOCITY, -MAX_VELOCITY, -0.5])
_UPPER = np.array([MAX_LUMP_AMP, MAX_LUMP_WIDTH, 64.0, 64.0, 64.0, MAX_VELOCITY, MAX_VELOCITY, MAX_VELOCITY, 0.5])


@dataclass(frozen=True)
class IterationConfig:
    """Configuration for the iterative adjustment loop."""

    max_evals: int = 80
    popsize: int = 8
    sigma0: float = 0.2
    max_concurrent_grtresna: int = 6
    mpi_ranks: int = 8
    plateau_generations: int = 5
    plateau_tolerance: float = 1e-4
    retention_threshold: float = 0.5
    L: float | None = None
    gridinit_n: int = 64
    grtresna_L: float = 128.0
    seed: int | None = None


@dataclass
class IterationResult:
    """Output summary of the iteration loop."""

    best_fitness: float
    best_vector: list[float]
    best_fitted_matter: FittedMatter
    best_gridinit_path: Path | None
    best_preservation_score: float
    total_evals: int
    generations: int
    converged: bool
    fitness_history: list[float]
    notes: list[str]


def _vector_from_lumps(lumps: Sequence[Mapping[str, Any]]) -> np.ndarray:
    """Encode lump parameters into a flat vector (warm-start)."""
    parts = []
    for lump in lumps:
        parts.extend([
            float(lump.get("amp", 0.1)),
            float(lump.get("width", 3.0)),
            float(lump.get("center", (0.0, 0.0, 0.0))[0]),
            float(lump.get("center", (0.0, 0.0, 0.0))[1]),
            float(lump.get("center", (0.0, 0.0, 0.0))[2]),
            float(lump.get("velocity", (0.0, 0.0, 0.0))[0]),
            float(lump.get("velocity", (0.0, 0.0, 0.0))[1]),
            float(lump.get("velocity", (0.0, 0.0, 0.0))[2]),
            float(lump.get("omega", 0.0)),
        ])
    return np.array(parts, dtype=float)


def _lumps_from_vector(
    vector: np.ndarray,
    *,
    template_lumps: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Decode a flat vector back into lump dicts (preserving exotic/mode flags)."""
    n_lumps = len(template_lumps)
    lumps = []
    for i in range(n_lumps):
        offset = i * _PARAMS_PER_LUMP
        v = vector[offset: offset + _PARAMS_PER_LUMP]
        template = template_lumps[i]
        lumps.append({
            "amp": float(v[0]),
            "width": float(v[1]),
            "center": (float(v[2]), float(v[3]), float(v[4])),
            "velocity": (float(v[5]), float(v[6]), float(v[7])),
            "omega": float(v[8]),
            "mode": int(template.get("mode", 0)),
            "exotic": int(template.get("exotic", 0)),
        })
    return lumps


def _clip_vector(vector: np.ndarray, n_lumps: int) -> np.ndarray:
    """Clip vector to per-lump bounds."""
    clipped = vector.copy()
    for i in range(n_lumps):
        offset = i * _PARAMS_PER_LUMP
        clipped[offset: offset + _PARAMS_PER_LUMP] = np.clip(
            clipped[offset: offset + _PARAMS_PER_LUMP], _LOWER, _UPPER
        )
    return clipped


def _evaluate_candidate(
    *,
    vector: np.ndarray,
    motif: GeometryMotif,
    template_lumps: Sequence[Mapping[str, Any]],
    base_fitted: FittedMatter,
    base_cfg: GRTresnaConfig,
    work_dir: Path,
    L: float | None,
) -> tuple[float, MismatchReport, Path | None]:
    """Single candidate evaluation: vector → lumps → solve → mismatch."""
    lumps = _lumps_from_vector(vector, template_lumps=template_lumps)

    # Build updated FittedMatter (preserving non-lump fields)
    fitted = FittedMatter(
        lumps=tuple(lumps),
        scalar_mass=base_fitted.scalar_mass,
        maximal_slicing=base_fitted.maximal_slicing,
        static_lens_only=base_fitted.static_lens_only,
        momentum_target=base_fitted.momentum_target,
        notes=base_fitted.notes,
    )

    cfg = build_grtresna_config_from_fitted(fitted, base=base_cfg)
    gridinit_path = work_dir / "initial_data.gridinit"
    grtresna_dir = work_dir / "grtresna"

    try:
        solve(cfg, work_dir=grtresna_dir, gridinit_path=gridinit_path)
    except (RuntimeError, FileNotFoundError, OSError) as exc:
        logger.debug("GRTresna solve failed for %s: %s", work_dir.name, exc)
        report = MismatchReport(
            fitness=GATE_FITNESS,
            chi_l2=0.0,
            beta_l2=0.0,
            exotic_penalty=0.0,
            convergence_penalty=0.0,
            solve_failed=True,
            notes=(f"solve failed: {exc}",),
        )
        return GATE_FITNESS, report, None

    report = compute_mismatch(
        motif,
        gridinit_path,
        lumps=lumps,
        grtresna_work_dir=grtresna_dir,
        L=L,
    )
    return report.fitness, report, gridinit_path


def run_iterate(
    motif: GeometryMotif,
    initial_fitted: FittedMatter,
    *,
    out_dir: Path,
    config: IterationConfig | None = None,
    base_grtresna_cfg: GRTresnaConfig | None = None,
) -> IterationResult:
    """Run the CMA-ES iteration loop to adjust matter for a geometry target.

    Parameters
    ----------
    motif : GeometryMotif
        Target geometry (carries overrides for profile generation).
    initial_fitted : FittedMatter
        Starting point from one-shot fit (warm-start for CMA-ES).
    out_dir : Path
        Output directory for logs, best gridinit, and fitted_matter JSONs.
    config : IterationConfig, optional
        Loop parameters (defaults are sensible for 1-3 lumps at n=64).
    base_grtresna_cfg : GRTresnaConfig, optional
        Base solver config (overridden by fitted-matter fields each eval).
    """
    try:
        import cma
    except ImportError as exc:
        raise ImportError(
            "Iteration loop requires the 'cma' package. Install: pip install cma"
        ) from exc

    cfg = config or IterationConfig()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    n_lumps = len(initial_fitted.lumps)
    if n_lumps == 0:
        raise ValueError("initial_fitted has no lumps to optimize")

    template_lumps = list(initial_fitted.lumps)
    x0 = _vector_from_lumps(template_lumps)
    x0 = _clip_vector(x0, n_lumps)

    # Base GRTresna config
    if base_grtresna_cfg is None:
        n_grid = cfg.gridinit_n
        base_grtresna_cfg = GRTresnaConfig(
            mpi_ranks=cfg.mpi_ranks,
            N=(n_grid, n_grid, n_grid // 2),
            L=cfg.grtresna_L,
            gridinit_nx=n_grid,
            gridinit_ny=n_grid,
            gridinit_nz=n_grid,
            scalar_mass=initial_fitted.scalar_mass,
            dphi=0.0,
            dpi=0.0,
            bh1_bare_mass=0.0,
            bh1_spin=(0.0, 0.0, 0.0),
        )

    # CMA-ES bounds
    lower = np.tile(_LOWER, n_lumps)
    upper = np.tile(_UPPER, n_lumps)
    bounds = [lower.tolist(), upper.tolist()]

    opts: dict[str, Any] = {
        "popsize": cfg.popsize,
        "bounds": bounds,
        "maxfevals": cfg.max_evals,
        "verbose": -9,  # suppress CMA-ES stdout
        "CMA_stds": (cfg.sigma0 * (upper - lower)).tolist(),
    }
    if cfg.seed is not None:
        opts["seed"] = cfg.seed

    es = cma.CMAEvolutionStrategy(x0.tolist(), cfg.sigma0, opts)

    # Logging
    log_path = out_dir / "iteration_log.jsonl"
    log_handle = open(log_path, "w", encoding="utf-8")

    best_fitness = GATE_FITNESS
    best_vector = x0.copy()
    best_gridinit: Path | None = None
    best_dir = out_dir / "best"
    best_dir.mkdir(exist_ok=True)
    fitness_history: list[float] = []
    notes: list[str] = []
    total_evals = 0
    generation = 0
    plateau_count = 0
    prev_best = GATE_FITNESS

    try:
        while not es.stop() and total_evals < cfg.max_evals:
            generation += 1
            solutions = es.ask()
            n_pop = len(solutions)
            clipped = [_clip_vector(np.array(s), n_lumps) for s in solutions]

            # Parallel GRTresna evaluations
            fitnesses = [GATE_FITNESS] * n_pop
            with ThreadPoolExecutor(max_workers=cfg.max_concurrent_grtresna) as pool:
                futures = {}
                for idx, vec in enumerate(clipped):
                    eval_id = total_evals + idx
                    eval_dir = out_dir / f"eval_{eval_id:06d}"
                    eval_dir.mkdir(exist_ok=True)
                    fut = pool.submit(
                        _evaluate_candidate,
                        vector=vec,
                        motif=motif,
                        template_lumps=template_lumps,
                        base_fitted=initial_fitted,
                        base_cfg=dataclasses.replace(base_grtresna_cfg),
                        work_dir=eval_dir,
                        L=cfg.L,
                    )
                    futures[fut] = idx

                for fut in as_completed(futures):
                    idx = futures[fut]
                    eval_id = total_evals + idx
                    try:
                        fitness, report, gridinit_path = fut.result()
                    except Exception as exc:
                        logger.warning("eval %d raised: %s", eval_id, exc)
                        fitness = GATE_FITNESS
                        report = MismatchReport(
                            fitness=GATE_FITNESS,
                            chi_l2=0.0,
                            beta_l2=0.0,
                            exotic_penalty=0.0,
                            convergence_penalty=0.0,
                            solve_failed=True,
                            notes=(f"exception: {exc}",),
                        )
                        gridinit_path = None

                    fitnesses[idx] = fitness

                    # Log this eval
                    log_entry = {
                        "eval_id": eval_id,
                        "generation": generation,
                        "fitness": fitness,
                        **report.to_dict(),
                    }
                    log_handle.write(json.dumps(log_entry) + "\n")
                    log_handle.flush()

                    # Track best
                    if fitness < best_fitness:
                        best_fitness = fitness
                        best_vector = clipped[idx].copy()
                        # Copy best gridinit to stable location
                        if gridinit_path is not None and gridinit_path.exists():
                            dst = best_dir / "initial_data.gridinit"
                            shutil.copy2(gridinit_path, dst)
                            best_gridinit = dst

            total_evals += n_pop
            es.tell(solutions, fitnesses)

            fitness_history.append(best_fitness)
            logger.info(
                "gen %d: best_fitness=%.6f, evals=%d/%d",
                generation, best_fitness, total_evals, cfg.max_evals,
            )

            # Plateau detection
            if abs(prev_best - best_fitness) < cfg.plateau_tolerance:
                plateau_count += 1
            else:
                plateau_count = 0
            prev_best = best_fitness

            if plateau_count >= cfg.plateau_generations:
                notes.append(
                    f"plateau stop after {cfg.plateau_generations} generations "
                    f"with no improvement > {cfg.plateau_tolerance}"
                )
                break

            # Prune non-best eval dirs to save disk
            for idx in range(n_pop):
                eval_id = total_evals - n_pop + idx
                eval_dir = out_dir / f"eval_{eval_id:06d}"
                if eval_dir.exists() and eval_dir != best_dir:
                    gridinit_in_dir = eval_dir / "initial_data.gridinit"
                    if not gridinit_in_dir.exists() or fitnesses[idx] > best_fitness:
                        shutil.rmtree(eval_dir, ignore_errors=True)

    finally:
        log_handle.close()

    # Build best FittedMatter for output
    best_lumps = _lumps_from_vector(best_vector, template_lumps=template_lumps)
    best_fitted = FittedMatter(
        lumps=tuple(best_lumps),
        scalar_mass=initial_fitted.scalar_mass,
        maximal_slicing=initial_fitted.maximal_slicing,
        static_lens_only=initial_fitted.static_lens_only,
        momentum_target=initial_fitted.momentum_target,
        notes=initial_fitted.notes + (f"iterate: {total_evals} evals, best_fitness={best_fitness:.6f}",),
    )
    write_fitted_matter_json(best_fitted, out_dir / "best_fitted_matter.json")

    # Preservation check on best gridinit
    preservation_score = 0.0
    if best_gridinit is not None and best_gridinit.exists():
        try:
            preservation = compare_motif_preservation(
                motif, best_gridinit, ftl_L=cfg.L,
            )
            preservation_score = preservation.retention_score
        except Exception as exc:
            notes.append(f"preservation check failed: {exc}")

    converged = preservation_score >= cfg.retention_threshold

    # Write summary
    summary = {
        "best_fitness": best_fitness,
        "best_preservation_score": preservation_score,
        "total_evals": total_evals,
        "generations": generation,
        "converged": converged,
        "fitness_history": fitness_history,
        "notes": notes,
    }
    (out_dir / "iterate_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    return IterationResult(
        best_fitness=best_fitness,
        best_vector=best_vector.tolist(),
        best_fitted_matter=best_fitted,
        best_gridinit_path=best_gridinit,
        best_preservation_score=preservation_score,
        total_evals=total_evals,
        generations=generation,
        converged=converged,
        fitness_history=fitness_history,
        notes=notes,
    )
