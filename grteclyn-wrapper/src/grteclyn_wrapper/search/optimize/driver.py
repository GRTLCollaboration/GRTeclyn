"""CMA-ES optimization loop driver."""

from __future__ import annotations

import math
import random
import threading
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from ...core.config import ExampleConfig, resolve_example
from ...core.episode import write_json
from ...core.scratch import purge_orphan_scratch
from ..grtresna_convergence_gate import (
    DEFAULT_GRTRESNA_CONVERGENCE_CONFIG,
    GRTresnaConvergenceConfig,
)
from ..ftl_retention import (
    FTL_CHAMPIONS_FILE,
    FTL_RETENTION_LOG,
    FtlChampionBoard,
    append_ftl_retention_events,
    compute_keep_eval_ids,
    save_ftl_champions,
)
from ..surrogate import RBFSurrogate, screen_candidates
from ..trajectory_log import format_trajectory_line
from .candidates import (
    _clip_vector,
    _force_exotic_template,
    _jitter_vector,
    _load_warm_start_vectors,
    _random_vector,
)
from .dimension import SearchDimension
from .eval import (
    ObjectiveMode,
    _collect_training,
    _evaluate_generation_parallel,
    _evaluate_generation_pipelined,
    _objective,
    _track_trajectory,
)
from .spaces import DEFAULT_SEARCH_SPACE, GRTRESNA_SEARCH_SPACE

try:
    import cma
except ImportError:
    cma = None  # type: ignore[assignment]

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]

@dataclass(frozen=True)
class OptimizeResult:
    best_params: dict[str, float]
    best_score: float
    best_episode: str
    generations: int
    evaluations: int
    trajectory: list[dict[str, Any]]


def _apply_optimize_retention(
    *,
    opt_dir: Path,
    trajectory: Sequence[Mapping[str, Any]],
    new_records: Sequence[Mapping[str, Any]],
    board: "FtlChampionBoard | None",
    keep_top_eval_dirs: int,
    ftl_retention_enabled: bool,
    champions_path: Path,
    retention_path: Path,
) -> None:
    """Keep top-N-by-score eval dirs (union FTL champions); prune the rest.

    Mirrors the QD ``_ingest`` retention so a CMA-ES run also saves only the
    best eval dirs on disk plus one champion dir per FTL peak metric
    (``ftl_champions.json`` / ``ftl_retention.jsonl``), deleting everything
    else as the conveyor advances.
    """
    # Lazy import: qd_search.io <-> optimize would otherwise form an import cycle
    # (qd_search/__init__ imports optimize's search-space defs).
    from ..qd_search.io import _prune_eval_dirs

    protect_eval_ids: set[int] = set()
    events = []
    for rec in new_records:
        eval_id = rec.get("eval")
        if isinstance(eval_id, int):
            protect_eval_ids.add(eval_id)
        if ftl_retention_enabled and board is not None and rec.get("status") == "gpu_ok":
            events.extend(board.consider(rec))
    if ftl_retention_enabled and board is not None:
        append_ftl_retention_events(retention_path, events)
        save_ftl_champions(champions_path, board)
    keep_ids = compute_keep_eval_ids(
        trajectory,
        keep_top_score=int(keep_top_eval_dirs),
        board=board if ftl_retention_enabled else None,
        ftl_retention_enabled=ftl_retention_enabled,
        protect_eval_ids=protect_eval_ids,
    )
    _prune_eval_dirs(
        opt_dir,
        trajectory,
        keep_eval_ids=keep_ids,
        protect_eval_ids=protect_eval_ids,
    )


def run_optimize(
    *,
    runs_dir: Path,
    executable: ExecutableConfig | None = None,
    max_generations: int = 50,
    population_size: int | None = None,
    sigma0: float = 0.3,
    seed: int | None = None,
    base_overrides: Mapping[str, Any] | None = None,
    search_space: Sequence[SearchDimension] | None = None,
    template: Path | None = None,
    example: ExampleConfig | str = "RadialRecipe",
    name: str | None = None,
    dry_run: bool = False,
    constrained: bool = True,
    phantom: bool = True,
    use_preflight: bool = True,
    cuda_devices: str | None = None,
    gpu_ids: Sequence[int] | None = None,
    check_params: bool = True,
    score_weights: Mapping[str, float] | None = None,
    objective_mode: ObjectiveMode = "weighted",
    ftl_L: float | None = None,
    x0: Sequence[float] | None = None,
    consume_plotfiles: bool = True,
    consumer_radii: Sequence[float] = (4.0, 8.0),
    consumer_keep_last: int = 1,
    surrogate: bool = False,
    surrogate_keep_fraction: float = 0.5,
    surrogate_warmup: int | None = None,
    grtresna: bool = False,
    grtresna_config: "GRTresnaConfig | None" = None,
    grtresna_solved_ftl_gate: bool | None = None,
    solved_ftl_gate_config: Any | None = None,
    grtresna_convergence_config: GRTresnaConvergenceConfig | None = None,
    grtresna_postload_gate: bool = False,
    postload_gate_config: Any | None = None,
    warm_start_trajectories: Sequence[Path] = (),
    warm_start_top_k: int = 8,
    warm_start_jitter: float = 0.08,
    random_injection_fraction: float = 0.0,
    exotic_injection_fraction: float = 0.0,
    keep_top_eval_dirs: int = 0,
    ftl_retention_enabled: bool = False,
    warm_start_include_near_miss: bool | None = None,
    warm_start_near_miss_k: int = 4,
    target_evals: int | None = None,
    use_pipeline: bool = True,
    slots_per_gpu: int = 1,
    cluster_cpu_fraction: float | None = None,
    pipeline_share: float | None = None,
    max_concurrent_grtresna_override: int = 0,
    reserve_cores: int | None = None,
    grtresna_mpi_ranks: int = 8,
) -> OptimizeResult:
    """Run multi-GPU CMA-ES optimization loop.

    Parameters
    ----------
    gpu_ids : list of int, optional
        Available GPU indices for parallel evaluation. If None, uses
        cuda_devices for sequential mode. If provided, each member of
        the CMA-ES population is assigned a GPU round-robin.
    """
    if cma is None:
        raise ImportError(
            "CMA-ES optimization requires the 'cma' package. "
            "Install it with: pip install cma"
        )
    if np is None:
        raise ImportError("numpy is required for optimization.")

    example_cfg = example if isinstance(example, ExampleConfig) else resolve_example(example)
    dims = list(search_space or (GRTRESNA_SEARCH_SPACE if grtresna else DEFAULT_SEARCH_SPACE))
    tpl = template or example_cfg.template
    base = dict(base_overrides or {})
    base.setdefault("N_full", 64)
    base.setdefault("max_level", 2)
    base.setdefault("stop_time", 2.0)
    base.setdefault("plot_interval", 10)
    base.setdefault("checkpoint_interval", -1)
    base.setdefault("dt_multiplier", 0.02)
    base.setdefault("regrid_threshold", 0.01)
    max_lvl = int(base["max_level"])
    if "regrid_interval" not in base and max_lvl > 0:
        intervals = [16] * min(max_lvl, 2) + [8] * max(0, max_lvl - 2)
        base["regrid_interval"] = intervals

    target_stop_time = float(base["stop_time"])
    solved_ftl_gate = grtresna if grtresna_solved_ftl_gate is None else grtresna_solved_ftl_gate

    if name is None:
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        name = f"optimize_{timestamp}"
    opt_dir = (runs_dir / name).expanduser().resolve()
    opt_dir.mkdir(parents=True, exist_ok=False)

    # Plotfiles live on the node's own disk now, which is the small one.  A
    # campaign killed mid-flight leaves directories nothing will return for;
    # reclaim them before adding a few hundred more.
    purge_orphan_scratch()

    # FTL champion retention + top-N disk pruning (same logic as the QD loop):
    # keep only the best eval dirs on disk plus one champion dir per FTL peak.
    ftl_board = FtlChampionBoard() if ftl_retention_enabled else None
    ftl_champions_path = opt_dir / FTL_CHAMPIONS_FILE
    ftl_retention_path = opt_dir / FTL_RETENTION_LOG
    retention_active = bool(ftl_retention_enabled or keep_top_eval_dirs > 0)

    warm_vectors = _load_warm_start_vectors(
        warm_start_trajectories,
        dims,
        max(0, warm_start_top_k),
        include_near_miss=(
            grtresna if warm_start_include_near_miss is None else warm_start_include_near_miss
        ),
        near_miss_k=warm_start_near_miss_k,
    )

    if x0 is not None:
        initial = list(x0)
    elif warm_vectors:
        initial = list(warm_vectors[0])
    else:
        initial = [d.center for d in dims]

    bounds = [[d.lower for d in dims], [d.upper for d in dims]]

    effective_popsize = population_size
    if effective_popsize is None and gpu_ids is not None:
        effective_popsize = len(gpu_ids)

    opts: dict[str, Any] = {
        "maxiter": max_generations,
        "bounds": bounds,
        "CMA_stds": [sigma0 * d.range for d in dims],
        "verbose": -1,
    }
    if effective_popsize is not None:
        opts["popsize"] = effective_popsize
    if seed is not None:
        opts["seed"] = seed

    eval_counter = [0]
    trajectory: list[dict[str, Any]] = []
    trajectory_path = opt_dir / "trajectory.jsonl"
    trajectory_lock = threading.Lock()

    gpu_pool = None
    cpu_admission = None
    pipeline_sizing: dict[str, int | float] = {}
    effective_gpu_ids = list(gpu_ids) if gpu_ids else []
    if use_pipeline and effective_gpu_ids and not dry_run:
        from ..gpu_pool import CpuAdmissionController, build_pipeline_resources

        gpu_pool, cpu_admission, pipeline_sizing = build_pipeline_resources(
            effective_gpu_ids,
            slots_per_gpu=slots_per_gpu,
            cluster_cpu_fraction=cluster_cpu_fraction,
            pipeline_share=pipeline_share,
            mpi_ranks=grtresna_mpi_ranks,
            reserve_cores=reserve_cores,
        )
        if max_concurrent_grtresna_override > 0:
            cpu_admission = CpuAdmissionController(max_concurrent_grtresna_override)
            pipeline_sizing["max_grtresna"] = max_concurrent_grtresna_override

    write_json(opt_dir / "metadata.json", {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "example": example_cfg.name,
        "max_generations": max_generations,
        "population_size": effective_popsize,
        "sigma0": sigma0,
        "seed": seed,
        "constrained": constrained,
        "phantom": phantom,
        "use_preflight": use_preflight,
        "dry_run": dry_run,
        "gpu_ids": list(gpu_ids) if gpu_ids else None,
        "consume_plotfiles": consume_plotfiles,
        "grtresna": grtresna,
        "grtresna_solved_ftl_gate": solved_ftl_gate,
        "grtresna_max_ham_pct": (
            grtresna_convergence_config.max_ham_pct
            if grtresna_convergence_config is not None
            else DEFAULT_GRTRESNA_CONVERGENCE_CONFIG.max_ham_pct
        ),
        "grtresna_max_mom_pct": (
            grtresna_convergence_config.max_mom_pct
            if grtresna_convergence_config is not None
            else DEFAULT_GRTRESNA_CONVERGENCE_CONFIG.max_mom_pct
        ),
        "objective_mode": objective_mode,
        "warm_start_trajectories": [str(p) for p in warm_start_trajectories],
        "warm_start_top_k": warm_start_top_k,
        "warm_start_jitter": warm_start_jitter,
        "warm_start_include_near_miss": (
            grtresna if warm_start_include_near_miss is None else warm_start_include_near_miss
        ),
        "warm_start_near_miss_k": warm_start_near_miss_k,
        "random_injection_fraction": random_injection_fraction,
        "exotic_injection_fraction": exotic_injection_fraction,
        "keep_top_eval_dirs": keep_top_eval_dirs,
        "ftl_retention_enabled": ftl_retention_enabled,
        "target_evals": target_evals,
        "use_pipeline": use_pipeline,
        "slots_per_gpu": slots_per_gpu,
        "max_concurrent_grtresna": pipeline_sizing.get("max_grtresna"),
        "gpu_slots": pipeline_sizing.get("gpu_slots"),
        "cluster_cpu_fraction": cluster_cpu_fraction,
        "pipeline_share": pipeline_share,
        "base_overrides": base,
        "search_space": [
            {"key": d.param_key, "lower": d.lower, "upper": d.upper, "initial": d.center}
            for d in dims
        ],
    })

    es = cma.CMAEvolutionStrategy(initial, sigma0, opts)
    rng = random.Random(seed)

    best_score = -math.inf
    best_params: dict[str, float] = {}
    best_episode = ""
    gen = 0

    print(f"[optimize] Starting CMA-ES: {len(dims)}D, popsize={es.popsize}, "
          f"max_gen={max_generations}, GPUs={gpu_ids or cuda_devices}, "
          f"pipeline={use_pipeline and gpu_pool is not None}")
    if target_evals is not None:
        print(f"[optimize] target_evals={target_evals}")
    if warm_vectors:
        print(f"[optimize] loaded {len(warm_vectors)} warm-start vectors")

    metadata_path = opt_dir / "metadata.json"
    common_eval_kwargs = dict(
        dims=dims,
        base_overrides=base,
        opt_dir=opt_dir,
        example=example_cfg,
        template=tpl,
        executable=executable,
        eval_counter=eval_counter,
        constrained=constrained,
        phantom=phantom,
        use_preflight=use_preflight,
        check_params=check_params,
        trajectory=trajectory,
        target_stop_time=target_stop_time,
        score_weights=score_weights,
        objective_mode=objective_mode,
        ftl_L=ftl_L,
        consume_plotfiles=consume_plotfiles,
        consumer_radii=consumer_radii,
        consumer_keep_last=consumer_keep_last,
        grtresna=grtresna,
        grtresna_base=grtresna_config,
        grtresna_solved_ftl_gate=solved_ftl_gate,
        solved_ftl_gate_config=solved_ftl_gate_config,
        grtresna_convergence_config=grtresna_convergence_config,
        grtresna_postload_gate=grtresna_postload_gate,
        postload_gate_config=postload_gate_config,
        live_trajectory_path=trajectory_path if not dry_run else None,
    )

    def _evaluate_subset(subset: list) -> list[float]:
        if use_pipeline and gpu_pool is not None and cpu_admission is not None and not dry_run:
            return _evaluate_generation_pipelined(
                subset,
                gpu_pool=gpu_pool,
                cpu_admission=cpu_admission,
                gpu_ids=effective_gpu_ids,
                metadata_path=metadata_path,
                **common_eval_kwargs,
            )
        if gpu_ids is not None and len(gpu_ids) > 1 and not dry_run:
            return _evaluate_generation_parallel(
                subset,
                gpu_ids=gpu_ids,
                **common_eval_kwargs,
            )
        out: list[float] = []
        for sol in subset:
            out.append(_objective(
                sol,
                cuda_devices=cuda_devices,
                dry_run=dry_run,
                **common_eval_kwargs,
            ))
        return out

    warmup = surrogate_warmup if surrogate_warmup is not None else 2 * es.popsize

    while not es.stop():
        if target_evals is not None and eval_counter[0] >= target_evals:
            break
        gen += 1
        solutions = [list(sol) for sol in es.ask()]
        n_solutions = len(solutions)
        if target_evals is not None:
            remaining = target_evals - eval_counter[0]
            if remaining <= 0:
                break
            if remaining < n_solutions:
                solutions = solutions[:remaining]
                n_solutions = len(solutions)

        if gen == 1 and warm_vectors:
            n_warm = min(len(warm_vectors), n_solutions)
            for i in range(n_warm):
                source = warm_vectors[i]
                if i == 0:
                    solutions[i] = _clip_vector(source, dims)
                else:
                    solutions[i] = _jitter_vector(
                        source, dims, rng, warm_start_jitter
                    )

        n_random = max(0, min(
            n_solutions,
            int(round(n_solutions * max(0.0, random_injection_fraction))),
        ))
        n_exotic = max(0, min(
            n_solutions - n_random,
            int(round(n_solutions * max(0.0, exotic_injection_fraction))),
        ))
        inject_pos = n_solutions - n_random - n_exotic
        for j in range(n_exotic):
            base_vec = (
                _jitter_vector(warm_vectors[j % len(warm_vectors)], dims, rng, warm_start_jitter)
                if warm_vectors else _random_vector(dims, rng)
            )
            solutions[inject_pos + j] = _force_exotic_template(
                base_vec, dims, rng, template_index=gen + j
            )
        for j in range(n_random):
            solutions[n_solutions - n_random + j] = _random_vector(dims, rng)

        traj_before = len(trajectory)

        x_train, y_train = _collect_training(trajectory, dims)
        use_surrogate_now = (
            surrogate
            and x_train is not None
            and x_train.shape[0] >= warmup
        )

        if use_surrogate_now:
            model = RBFSurrogate(
                lower=np.asarray([d.lower for d in dims], dtype=float),
                upper=np.asarray([d.upper for d in dims], dtype=float),
            ).fit(x_train, y_train)
            sols_arr = np.asarray(solutions, dtype=float)
            eval_idx, predicted = screen_candidates(
                model, sols_arr,
                keep_fraction=surrogate_keep_fraction,
                min_eval=1,
            )
            subset = [solutions[i] for i in eval_idx]
            sub_fit = _evaluate_subset(subset)
            fitnesses = [-float(predicted[i]) for i in range(len(solutions))]
            for k, i in enumerate(eval_idx):
                fitnesses[i] = sub_fit[k]
            n_skipped = len(solutions) - len(eval_idx)
            for i in range(len(solutions)):
                if i not in set(eval_idx):
                    _track_trajectory(trajectory, {
                        "eval": None,
                        "surrogate_predicted": True,
                        "score": float(predicted[i]),
                        "overrides": {d.param_key: float(solutions[i][j]) for j, d in enumerate(dims)},
                    }, trajectory_lock=trajectory_lock, live_trajectory_path=trajectory_path if not dry_run else None)
            print(f"[optimize] gen {gen}: surrogate screened "
                  f"{len(eval_idx)}/{len(solutions)} evaluated on GPU "
                  f"({n_skipped} predicted)")
        else:
            fitnesses = _evaluate_subset(list(solutions))

        es.tell(solutions, fitnesses)

        evaluated_records = [
            rec for rec in trajectory[traj_before:]
            if not rec.get("surrogate_predicted")
        ]
        gen_scores = [rec.get("score", -math.inf) for rec in evaluated_records]
        gen_best = max(gen_scores) if gen_scores else -math.inf
        gen_mean = sum(gen_scores) / len(gen_scores) if gen_scores else 0.0

        for rec in evaluated_records:
            sc = rec.get("score", -math.inf)
            if sc > best_score:
                best_score = sc
                best_params = dict(rec.get("overrides", {}))
                best_episode = rec.get("episode", "")

        print(f"[optimize] gen {gen}/{max_generations}: "
              f"best={gen_best:.4f} mean={gen_mean:.4f} "
              f"all-time-best={best_score:.4f} evals={eval_counter[0]}")

        with (opt_dir / "trajectory.jsonl").open("w", encoding="utf-8") as fh:
            for rec in trajectory:
                fh.write(format_trajectory_line(rec))

        if retention_active and not dry_run:
            _apply_optimize_retention(
                opt_dir=opt_dir,
                trajectory=trajectory,
                new_records=evaluated_records,
                board=ftl_board,
                keep_top_eval_dirs=keep_top_eval_dirs,
                ftl_retention_enabled=ftl_retention_enabled,
                champions_path=ftl_champions_path,
                retention_path=ftl_retention_path,
            )

    result = OptimizeResult(
        best_params=best_params,
        best_score=best_score,
        best_episode=best_episode,
        generations=gen,
        evaluations=eval_counter[0],
        trajectory=trajectory,
    )

    write_json(opt_dir / "result.json", {
        "best_params": result.best_params,
        "best_score": result.best_score,
        "best_episode": result.best_episode,
        "generations": result.generations,
        "evaluations": result.evaluations,
    })

    print(f"\n[optimize] Done. {gen} generations, {eval_counter[0]} evaluations.")
    print(f"[optimize] Best score: {best_score:.4f}")
    print(f"[optimize] Best episode: {best_episode}")
    print(f"[optimize] Results: {opt_dir}")

    return result
