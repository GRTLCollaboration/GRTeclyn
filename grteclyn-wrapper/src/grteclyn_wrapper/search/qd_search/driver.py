"""MAP-Elites search driver."""

from __future__ import annotations

import json
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from ...core.config import ExampleConfig, ExecutableConfig, resolve_example
from ...core.episode import write_json
from ...core.params import regrid_intervals_for_max_level
from ...core.evaluation import (
    Evaluation,
    _run_cpu_grtresna_gates,
    _run_gpu_session,
)
from ..eval_pipeline import (
    EvalJob,
    EvalPipeline,
    resume_eval_counter,
    scan_partial_eval_dirs,
)
from ..gpu_pool import (
    CpuAdmissionController,
    build_pipeline_resources,
)
from ..scoring_pool import ScoringPool, scoring_workers_from_env
from ..optimize import DEFAULT_SEARCH_SPACE, SearchDimension
from ..optimize.candidates import _vector_to_overrides
from ..trajectory_log import (
    format_trajectory_line,
    infer_trajectory_status,
    trajectory_flags_from_evaluation,
)
from ..validation_tiers import (
    DEFAULT_TIER_CONFIG,
    Tier,
    TierConfig,
    build_survivors,
    convergence_signals,
    evaluate_tiers,
    survivor_front,
)
from .. import qd_search
from .archive import Elite, QDArchive
from .descriptors import _bin_index, _descriptor_details
from ..ftl_retention import (
    FTL_CHAMPIONS_FILE,
    FTL_RETENTION_LOG,
    FtlChampionBoard,
    append_ftl_retention_events,
    compute_keep_eval_ids,
    save_ftl_champions,
)
from .io import _iterations_for_target_evals, _load_trajectory_records, _prune_eval_dirs
from ..pre_gpu import NearMissPool, attach_pre_gpu_metrics_to_record, is_graded_rejection
from .pre_gpu_archive import (
    PRE_GPU_ARCHIVE_FILE,
    load_pre_gpu_archive,
    rebuild_pre_gpu_archive_from_trajectory,
    shadow_elite_from_record,
)
from .sampling import QdSamplingConfig, sample_next_candidate, _sample_random


def run_qd_search(
    *,
    runs_dir: Path,
    executable: ExecutableConfig | None = None,
    iterations: int = 10,
    batch_size: int | None = None,
    bins: int = 8,
    init_random: int | None = None,
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
    ftl_L: float | None = None,
    consume_plotfiles: bool = True,
    consumer_radii: Sequence[float] = (4.0, 8.0),
    consumer_keep_last: int = 1,
    tier_config: TierConfig = DEFAULT_TIER_CONFIG,
    survivor_min_tier: int = int(Tier.OPERATIONAL),
    descriptor_mode: str = "legacy",
    objective_mode: str = "weighted",
    grtresna: bool = False,
    grtresna_config: Any | None = None,
    grtresna_solved_ftl_gate: bool = False,
    solved_ftl_gate_config: Any | None = None,
    grtresna_convergence_config: Any | None = None,
    grtresna_postload_gate: bool = False,
    postload_gate_config: Any | None = None,
    resume: bool = False,
    target_evals: int | None = None,
    seed_overrides: Sequence[Mapping[str, Any]] | None = None,
    keep_top_eval_dirs: int = 0,
    ftl_retention_enabled: bool = True,
    slots_per_gpu: int = 1,
    cluster_cpu_fraction: float = 0.30,
    pipeline_share: float = 1.0,
    max_concurrent_grtresna_override: int = 0,
    reserve_cores: int = 4,
    grtresna_mpi_ranks: int = 8,
    remove_partial: bool = False,
    use_pipeline: bool = True,
    prune_interval: int = 10,
    pre_gpu_learning: bool | None = None,
    near_miss_pool_size: int = 32,
) -> QDArchive:
    example_cfg = example if isinstance(example, ExampleConfig) else resolve_example(example)
    dims = list(search_space or DEFAULT_SEARCH_SPACE)
    tpl = template or example_cfg.template
    base = dict(base_overrides or {})
    base.setdefault("N_full", 64)
    base.setdefault("max_level", 1)
    base.setdefault("stop_time", 2.0)
    base.setdefault("plot_interval", 10)
    base.setdefault("checkpoint_interval", -1)
    base.setdefault("dt_multiplier", 0.02)
    max_lvl = int(base["max_level"])
    if "regrid_interval" not in base and max_lvl > 0:
        base["regrid_interval"] = regrid_intervals_for_max_level(max_lvl)
    target_stop_time = float(base["stop_time"])
    pre_gpu_enabled = grtresna if pre_gpu_learning is None else bool(pre_gpu_learning and grtresna)

    rng = np.random.default_rng(seed)
    effective_gpu_ids = list(gpu_ids) if gpu_ids else (
        [int(cuda_devices)] if cuda_devices is not None else [0]
    )
    batch = batch_size or len(effective_gpu_ids)
    n_init = init_random if init_random is not None else batch

    if name is None:
        if resume:
            raise ValueError("resume requires --name pointing at an existing qd_<timestamp> directory")
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        name = f"qd_{timestamp}"
    qd_dir = (runs_dir / name).expanduser().resolve()

    trajectory_path = qd_dir / "trajectory.jsonl"
    metadata_path = qd_dir / "metadata.json"
    conv_history: list[dict[str, Any]] = []
    completed_evals = 0

    if resume:
        if not qd_dir.is_dir():
            raise FileNotFoundError(f"resume directory not found: {qd_dir}")
        archive_path = qd_dir / "archive.json"
        if not archive_path.exists():
            raise FileNotFoundError(f"resume requires archive.json in {qd_dir}")
        archive = QDArchive.from_dict(json.loads(archive_path.read_text(encoding="utf-8")))
        if archive.bins != bins:
            raise ValueError(
                f"resume bins mismatch: archive has {archive.bins}, requested {bins}"
            )
        trajectory = _load_trajectory_records(trajectory_path)
        traj_ids = {int(rec["eval"]) for rec in trajectory if "eval" in rec}
        partial_records = scan_partial_eval_dirs(
            qd_dir,
            trajectory_eval_ids=traj_ids,
            remove_partial=remove_partial,
        )
        for record in partial_records:
            trajectory.append(record)
            with trajectory_path.open("a", encoding="utf-8") as fh:
                fh.write(format_trajectory_line(record))
        completed_evals = resume_eval_counter(
            metadata_path=metadata_path,
            trajectory_eval_ids=traj_ids,
        )
        validation_path = qd_dir / "validation.json"
        if validation_path.exists():
            validation = json.loads(validation_path.read_text(encoding="utf-8"))
            conv_history = list(validation.get("history", []))
        if metadata_path.exists() and seed is None:
            meta = json.loads(metadata_path.read_text(encoding="utf-8"))
            seed = meta.get("seed")
        if target_evals is not None:
            iterations = _iterations_for_target_evals(
                target_evals=target_evals,
                batch=batch,
                completed_evals=completed_evals,
            )
        print(
            f"[qd] resume {qd_dir.name}: completed={completed_evals} "
            f"elites={len(archive.cells)} planned_batches={iterations}"
            + (f" target_evals={target_evals}" if target_evals is not None else "")
        )
        if target_evals is not None and iterations == 0:
            print(f"[qd] target_evals={target_evals} already reached; nothing to do.")
            return archive
    else:
        qd_dir.mkdir(parents=True, exist_ok=False)
        archive = QDArchive(bins=bins)
        trajectory = []
        if target_evals is not None:
            iterations = _iterations_for_target_evals(
                target_evals=target_evals,
                batch=batch,
            )

    eval_counter = [completed_evals]
    counter_lock = threading.Lock()
    trajectory_lock = threading.Lock()
    archive_lock = threading.Lock()
    ftl_retention_path = qd_dir / FTL_RETENTION_LOG
    ftl_champions_path = qd_dir / FTL_CHAMPIONS_FILE
    ftl_board = (
        FtlChampionBoard.rebuild(trajectory) if ftl_retention_enabled else FtlChampionBoard()
    )

    pre_gpu_archive_path = qd_dir / PRE_GPU_ARCHIVE_FILE
    if pre_gpu_enabled:
        try:
            pre_gpu_archive = load_pre_gpu_archive(pre_gpu_archive_path, bins=bins)
        except (json.JSONDecodeError, KeyError, ValueError):
            pre_gpu_archive = rebuild_pre_gpu_archive_from_trajectory(
                trajectory,
                bins=bins,
                convergence_config=grtresna_convergence_config,
            )
        near_miss_pool = NearMissPool.rebuild_from_trajectory(
            trajectory, max_size=near_miss_pool_size,
        )
    else:
        pre_gpu_archive = QDArchive(bins=bins)
        near_miss_pool = NearMissPool(max_size=0)
    sampling_config = QdSamplingConfig()

    gpu_pool, cpu_admission, sizing = build_pipeline_resources(
        effective_gpu_ids,
        slots_per_gpu=slots_per_gpu,
        cluster_cpu_fraction=cluster_cpu_fraction,
        pipeline_share=pipeline_share,
        mpi_ranks=grtresna_mpi_ranks,
        reserve_cores=reserve_cores,
    )
    if max_concurrent_grtresna_override > 0:
        cpu_admission = CpuAdmissionController(max_concurrent_grtresna_override)
        sizing["max_grtresna"] = max_concurrent_grtresna_override

    # Scoring runs off the GPU-lease path in its own bounded process pool so the
    # CPU-heavy yt analysis overlaps with GPU evolutions instead of holding a
    # GPU slot (default: one scorer per GPU slot; override via SCORING_WORKERS).
    scoring_workers = scoring_workers_from_env(gpu_pool.total_slots)
    scoring_pool = ScoringPool(scoring_workers)
    sizing["scoring_workers"] = scoring_workers

    metadata = {
        "mode": "qd_map_elites",
        "example": example_cfg.name,
        "iterations": iterations,
        "batch_size": batch,
        "bins": bins,
        "seed": seed,
        "gpu_ids": list(effective_gpu_ids),
        "descriptor_mode": descriptor_mode,
        "base_overrides": base,
        "search_space": [
            {"key": d.param_key, "lower": d.lower, "upper": d.upper} for d in dims
        ],
        "grtresna": grtresna,
        "grtresna_solved_ftl_gate": grtresna_solved_ftl_gate,
        "keep_top_eval_dirs": keep_top_eval_dirs,
        "ftl_retention_enabled": ftl_retention_enabled,
        "grtresna_max_ham_pct": getattr(
            grtresna_convergence_config, "max_ham_pct", 5.0,
        ),
        "grtresna_max_mom_pct": getattr(
            grtresna_convergence_config, "max_mom_pct", 5.0,
        ),
        "slots_per_gpu": slots_per_gpu,
        "max_concurrent_grtresna": sizing["max_grtresna"],
        "scoring_workers": sizing.get("scoring_workers"),
        "gpu_slots": sizing["gpu_slots"],
        "cluster_cpu_fraction": cluster_cpu_fraction,
        "pipeline_share": pipeline_share,
        "last_eval_counter": completed_evals,
        "use_pipeline": use_pipeline,
        "pre_gpu_learning": pre_gpu_enabled,
        "near_miss_pool_size": near_miss_pool_size if pre_gpu_enabled else 0,
    }
    if target_evals is not None:
        metadata["target_evals"] = target_evals
    if resume:
        metadata["resumed_at"] = datetime.now(timezone.utc).isoformat()
        metadata["resume_from_evals"] = completed_evals
        metadata["resume_planned_batches"] = iterations
        if metadata_path.exists():
            existing = json.loads(metadata_path.read_text(encoding="utf-8"))
            metadata["created_at"] = existing.get(
                "created_at", datetime.now(timezone.utc).isoformat()
            )
        else:
            metadata["created_at"] = datetime.now(timezone.utc).isoformat()
    else:
        metadata["created_at"] = datetime.now(timezone.utc).isoformat()
    write_json(metadata_path, metadata)

    completions_since_prune = [0]
    completions_since_report = [0]
    stop_event = threading.Event()

    def _append_live_trajectory(record: Mapping[str, Any]) -> None:
        with trajectory_lock:
            with trajectory_path.open("a", encoding="utf-8") as fh:
                fh.write(format_trajectory_line(record))

    def _build_record_and_elite(
        *,
        idx: int,
        x: list[float],
        res: Evaluation | None,
    ) -> tuple[dict[str, Any], Elite | None]:
        overrides = _vector_to_overrides(x, dims, base)
        if res is None or res.preflight_rejected:
            record: dict[str, Any] = {
                "eval": idx,
                "preflight_rejected": True,
                "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
            }
            if res is not None:
                record.update(trajectory_flags_from_evaluation(res))
            record["status"] = infer_trajectory_status(record)
            return record, None

        descriptor_details = _descriptor_details(
            res.components,
            res.metrics,
            mode=descriptor_mode,
            overrides={d.param_key: overrides.get(d.param_key) for d in dims},
        )
        d1, d2 = descriptor_details["x"], descriptor_details["y"]
        cell = (_bin_index(d1, bins), _bin_index(d2, bins))
        params = {d.param_key: float(overrides[d.param_key]) for d in dims}
        assessment = evaluate_tiers(res.components, metrics=res.metrics, config=tier_config)
        flags = trajectory_flags_from_evaluation(res)
        status = infer_trajectory_status({**flags, "components": res.components})

        elite = None
        if status == "gpu_ok":
            elite = Elite(
                cell=cell,
                score=res.score,
                descriptors=(d1, d2),
                params=params,
                episode=res.episode_path,
                tier=assessment.tier,
                tier_name=assessment.tier_name,
                objectives=assessment.objectives,
                descriptor_details=descriptor_details,
            )

        record = {
            "eval": idx,
            "score": res.score,
            "cell": list(cell),
            "descriptors": [d1, d2],
            "descriptor_details": descriptor_details,
            "tier": assessment.tier,
            "tier_name": assessment.tier_name,
            "improved": None,
            "episode": res.episode_path,
            "components": res.components,
            "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
        }
        record.update(flags)
        record["status"] = status
        attach_pre_gpu_metrics_to_record(record, res.metrics)
        return record, elite

    def _sample_next_candidate(worker_rng: np.random.Generator) -> list[float]:
        if not pre_gpu_enabled:
            elites = list(archive.cells.values())
            if elites and worker_rng.random() < 0.85:
                from .sampling import _mutate_elite
                parent = elites[int(worker_rng.integers(len(elites)))]
                return _mutate_elite(parent, dims, worker_rng)
            if elites:
                from .sampling import _feasible_bounds, _sample_feasible_box
                return _sample_feasible_box(
                    dims, _feasible_bounds(dims, elites), worker_rng,
                )
            return _sample_random(dims, worker_rng)
        return sample_next_candidate(
            dims=dims,
            rng=worker_rng,
            archive=archive,
            pre_gpu_archive=pre_gpu_archive,
            near_miss_pool=near_miss_pool,
            config=sampling_config,
        )

    def _batched_prune() -> None:
        in_flight = pipeline.in_flight_eval_ids if use_pipeline else set()
        keep_ids = compute_keep_eval_ids(
            trajectory,
            keep_top_score=int(keep_top_eval_dirs),
            board=ftl_board if ftl_retention_enabled else None,
            ftl_retention_enabled=ftl_retention_enabled,
            protect_eval_ids=in_flight,
        )
        _prune_eval_dirs(
            qd_dir,
            trajectory,
            keep_eval_ids=keep_ids,
            protect_eval_ids=in_flight,
        )
        metadata["last_eval_counter"] = eval_counter[0]
        write_json(metadata_path, metadata)

    def _write_validation() -> dict[str, Any]:
        items = [
            {
                "label": f"cell_{e.cell[0]}_{e.cell[1]}",
                "tier": e.tier,
                "score": e.score,
                "objectives": e.objectives,
                "episode": e.episode,
            }
            for e in archive.cells.values()
        ]
        survivors = build_survivors(items, min_tier=survivor_min_tier)
        front = survivor_front(survivors)
        best = archive.best
        snapshot = {
            "iteration": len(conv_history),
            "coverage": archive.coverage,
            "best_score": best.score if best else 0.0,
            "front_labels": [s.label for s in front],
        }
        conv_history.append(snapshot)
        signals = convergence_signals(conv_history)
        write_json(qd_dir / "validation.json", {
            "tier_histogram": archive.tier_histogram(),
            "survivor_min_tier": survivor_min_tier,
            "num_survivors": len(survivors),
            "survivor_front": [
                {
                    "label": s.label,
                    "tier": s.tier,
                    "tier_name": next(
                        (e.tier_name for e in archive.cells.values()
                         if f"cell_{e.cell[0]}_{e.cell[1]}" == s.label),
                        "",
                    ),
                    "score": s.score,
                    "objectives": s.objectives,
                    "episode": s.episode,
                }
                for s in front
            ],
            "convergence": signals,
            "history": conv_history,
        })
        return signals

    def _evaluate_job(job: EvalJob[list[float]]) -> Evaluation:
        overrides = _vector_to_overrides(job.payload, dims, base)
        name_eval = f"eval_{job.eval_id:06d}"
        if not grtresna:
            return qd_search.evaluate_overrides(
                overrides,
                out_dir=qd_dir,
                name=name_eval,
                example=example_cfg,
                template=tpl,
                executable=executable,
                constrained=constrained,
                phantom=phantom,
                use_preflight=use_preflight,
                check_params=check_params,
                dry_run=dry_run,
                target_stop_time=target_stop_time,
                score_weights=score_weights,
                objective_mode=objective_mode,
                ftl_L=ftl_L,
                consume_plotfiles=consume_plotfiles,
                consumer_radii=consumer_radii,
                consumer_keep_last=consumer_keep_last,
                gpu_pool=gpu_pool,
                scoring_pool=scoring_pool,
            )

        with cpu_admission.admit():
            cpu_phase = _run_cpu_grtresna_gates(
                overrides=overrides,
                out_dir=qd_dir,
                name=name_eval,
                example=example_cfg,
                grtresna_base=grtresna_config,
                grtresna_convergence_config=grtresna_convergence_config,
                grtresna_solved_ftl_gate=grtresna_solved_ftl_gate,
                solved_ftl_gate_config=solved_ftl_gate_config,
                ftl_L=ftl_L,
                dry_run=dry_run,
            )
        if isinstance(cpu_phase, Evaluation):
            return cpu_phase
        cpu_result = cpu_phase

        # The GPU slot is leased inside _run_gpu_session for the evolution binary
        # only; scoring runs in scoring_pool after the slot is released.
        return _run_gpu_session(
            cpu_result=cpu_result,
            example=example_cfg,
            template=tpl,
            executable=executable,
            check_params=check_params,
            dry_run=dry_run,
            target_stop_time=target_stop_time,
            score_weights=score_weights,
            objective_mode=objective_mode,
            ftl_L=ftl_L,
            consume_plotfiles=consume_plotfiles,
            consumer_radii=consumer_radii,
            consumer_keep_last=consumer_keep_last,
            consumer_ftl_timeseries=True,
            consumer_evolving_geodesic=None,
            grtresna_postload_gate=grtresna_postload_gate,
            postload_gate_config=postload_gate_config,
            gpu_pool=gpu_pool,
            scoring_pool=scoring_pool,
        )

    _worker_rng_local = threading.local()
    _worker_rng_seeds = list(
        np.random.SeedSequence(seed if seed is not None else 0).spawn(
            gpu_pool.total_slots + cpu_admission.max_concurrent + 4
        )
    )
    _worker_rng_assign_lock = threading.Lock()
    _worker_rng_next_idx = [0]

    def _get_worker_rng() -> np.random.Generator:
        rng = getattr(_worker_rng_local, "rng", None)
        if rng is None:
            with _worker_rng_assign_lock:
                idx = _worker_rng_next_idx[0]
                _worker_rng_next_idx[0] += 1
                child_seed = _worker_rng_seeds[idx % len(_worker_rng_seeds)]
            _worker_rng_local.rng = np.random.default_rng(child_seed)
            rng = _worker_rng_local.rng
        return rng

    def _on_eval_complete(job: EvalJob[list[float]], res: Evaluation) -> None:
        record, elite = _build_record_and_elite(idx=job.eval_id, x=job.payload, res=res)
        retention_events: list[dict[str, Any]] = []
        with archive_lock:
            if elite is not None:
                record["improved"] = archive.insert(elite)
            if pre_gpu_enabled and is_graded_rejection(record):
                near_miss_pool.consider(record)
                shadow = shadow_elite_from_record(
                    record,
                    bins=bins,
                    convergence_config=grtresna_convergence_config,
                )
                if shadow is not None:
                    record["shadow_improved"] = pre_gpu_archive.insert(shadow)
            trajectory.append(record)
            if ftl_retention_enabled and record.get("status") == "gpu_ok":
                retention_events = ftl_board.consider(record)
        _append_live_trajectory(record)
        if ftl_retention_enabled:
            if retention_events:
                append_ftl_retention_events(ftl_retention_path, retention_events)
            if record.get("status") == "gpu_ok":
                save_ftl_champions(ftl_champions_path, ftl_board)

        completions_since_prune[0] += 1
        completions_since_report[0] += 1
        if completions_since_prune[0] >= prune_interval:
            completions_since_prune[0] = 0
            _batched_prune()
        if completions_since_report[0] >= batch:
            completions_since_report[0] = 0
            write_json(qd_dir / "archive.json", archive.to_dict())
            if pre_gpu_enabled:
                write_json(pre_gpu_archive_path, pre_gpu_archive.to_dict())
            signals = _write_validation()
            best = archive.best
            print(
                f"[qd] report: elites={len(archive.cells)} coverage={archive.coverage:.2f} "
                f"best={best.score if best else float('nan'):.4f} "
                f"evals={eval_counter[0]}/{total_target}"
                + (
                    f" near_miss={len(near_miss_pool)} "
                    f"shadow_cov={pre_gpu_archive.coverage:.2f}"
                    if pre_gpu_enabled
                    else ""
                )
            )
            if signals.get("converged") and (
                target_evals is None or eval_counter[0] >= total_target
            ):
                stop_event.set()

        if eval_counter[0] < total_target and not stop_event.is_set():
            with archive_lock:
                new_x = _sample_next_candidate(_get_worker_rng())
            pipeline.submit(new_x)

    # Keep a job resident at every stage (GRTresna solve, GPU evolve, score) so
    # no stage starves while another is busy. Scoring now runs off the GPU-lease
    # path, so it needs its own share of in-flight slots.
    max_in_flight = (
        gpu_pool.total_slots + cpu_admission.max_concurrent + scoring_pool.max_workers
    )
    pipeline = EvalPipeline(
        gpu_pool=gpu_pool,
        cpu_admission=cpu_admission,
        evaluate_fn=_evaluate_job,
        on_complete=_on_eval_complete,
        counter_lock=counter_lock,
        eval_counter=eval_counter,
        metadata_path=metadata_path,
        max_in_flight=max_in_flight,
    )

    start_iter = len(conv_history)
    n_seed = 0 if resume else len(seed_overrides or [])
    total_target = (
        target_evals
        if target_evals is not None
        else completed_evals + (0 if resume else n_init) + n_seed + iterations * batch
    )
    print(
        f"[qd] MAP-Elites: {len(dims)}D, bins={bins}x{bins}, batch={batch}, "
        f"iters={iterations}, GPUs={effective_gpu_ids}, pipeline={use_pipeline}, "
        f"evals={completed_evals}->{total_target}"
    )

    if not use_pipeline:
        # Legacy batch mode scores inline; the pipeline scoring pool is unused.
        scoring_pool.shutdown(wait=False)
        return _run_legacy_batch_mode(
            qd_dir=qd_dir,
            archive=archive,
            trajectory=trajectory,
            trajectory_path=trajectory_path,
            eval_counter=eval_counter,
            counter_lock=counter_lock,
            dims=dims,
            base=base,
            bins=bins,
            example_cfg=example_cfg,
            tpl=tpl,
            executable=executable,
            constrained=constrained,
            phantom=phantom,
            use_preflight=use_preflight,
            cuda_devices=cuda_devices,
            effective_gpu_ids=effective_gpu_ids,
            check_params=check_params,
            dry_run=dry_run,
            target_stop_time=target_stop_time,
            score_weights=score_weights,
            objective_mode=objective_mode,
            ftl_L=ftl_L,
            consume_plotfiles=consume_plotfiles,
            consumer_radii=consumer_radii,
            consumer_keep_last=consumer_keep_last,
            grtresna=grtresna,
            grtresna_config=grtresna_config,
            grtresna_solved_ftl_gate=grtresna_solved_ftl_gate,
            solved_ftl_gate_config=solved_ftl_gate_config,
            grtresna_convergence_config=grtresna_convergence_config,
            grtresna_postload_gate=grtresna_postload_gate,
            postload_gate_config=postload_gate_config,
            resume=resume,
            seed_overrides=seed_overrides,
            n_init=n_init,
            iterations=iterations,
            batch=batch,
            descriptor_mode=descriptor_mode,
            tier_config=tier_config,
            survivor_min_tier=survivor_min_tier,
            rng=rng,
            ftl_retention_enabled=ftl_retention_enabled,
            ftl_board=ftl_board,
            ftl_retention_path=ftl_retention_path,
            ftl_champions_path=ftl_champions_path,
            keep_top_eval_dirs=keep_top_eval_dirs,
            conv_history=conv_history,
            start_iter=start_iter,
            total_target=total_target,
            pre_gpu_enabled=pre_gpu_enabled,
            near_miss_pool=near_miss_pool,
            pre_gpu_archive=pre_gpu_archive,
            sampling_config=sampling_config,
            pre_gpu_archive_path=pre_gpu_archive_path,
        )

    if not resume:
        seed_vectors = [
            [float(s.get(d.param_key, d.center)) for d in dims]
            for s in (seed_overrides or [])
        ]
        # Fill every pipeline stage on startup (previously only n_init≈batch,
        # which left GPUs idle while the batch drained through CPU stages).
        remaining = max(0, total_target - eval_counter[0])
        n_fill = min(remaining, max_in_flight)
        n_random = max(n_init, n_fill - len(seed_vectors))
        n_random = min(n_random, max(0, remaining - len(seed_vectors)))
        init_vectors = seed_vectors + [
            _sample_random(dims, rng) for _ in range(n_random)
        ]
        pipeline.submit_many(init_vectors)
    else:
        remaining = max(0, total_target - eval_counter[0])
        refill = min(remaining, max_in_flight)
        if refill > 0:
            with archive_lock:
                vectors = [_sample_next_candidate(_get_worker_rng()) for _ in range(refill)]
            pipeline.submit_many(vectors)

    pipeline.wait_until_idle()
    pipeline.shutdown()
    scoring_pool.shutdown()
    _batched_prune()
    write_json(qd_dir / "archive.json", archive.to_dict())
    if pre_gpu_enabled:
        write_json(pre_gpu_archive_path, pre_gpu_archive.to_dict())
    _write_validation()
    final_hist = archive.tier_histogram()
    print(f"[qd] Done. elites={len(archive.cells)} coverage={archive.coverage:.2f} "
          f"tiers={final_hist} dir={qd_dir}")
    return archive


def _run_legacy_batch_mode(
    *,
    qd_dir: Path,
    archive: QDArchive,
    trajectory: list[dict[str, Any]],
    trajectory_path: Path,
    eval_counter: list[int],
    counter_lock: threading.Lock,
    dims: Sequence[SearchDimension],
    base: dict[str, Any],
    bins: int,
    example_cfg: ExampleConfig,
    tpl: Path,
    executable: ExecutableConfig | None,
    constrained: bool,
    phantom: bool,
    use_preflight: bool,
    cuda_devices: str | None,
    effective_gpu_ids: Sequence[int],
    check_params: bool,
    dry_run: bool,
    target_stop_time: float,
    score_weights: Mapping[str, float] | None,
    objective_mode: str,
    ftl_L: float | None,
    consume_plotfiles: bool,
    consumer_radii: Sequence[float],
    consumer_keep_last: int,
    grtresna: bool,
    grtresna_config: Any | None,
    grtresna_solved_ftl_gate: bool,
    solved_ftl_gate_config: Any | None,
    grtresna_convergence_config: Any | None,
    grtresna_postload_gate: bool,
    postload_gate_config: Any | None,
    resume: bool,
    seed_overrides: Sequence[Mapping[str, Any]] | None,
    n_init: int,
    iterations: int,
    batch: int,
    descriptor_mode: str,
    tier_config: TierConfig,
    survivor_min_tier: int,
    rng: np.random.Generator,
    ftl_retention_enabled: bool,
    ftl_board: FtlChampionBoard,
    ftl_retention_path: Path,
    ftl_champions_path: Path,
    keep_top_eval_dirs: int,
    conv_history: list[dict[str, Any]],
    start_iter: int,
    total_target: int,
    pre_gpu_enabled: bool = False,
    near_miss_pool: NearMissPool | None = None,
    pre_gpu_archive: QDArchive | None = None,
    sampling_config: QdSamplingConfig | None = None,
    pre_gpu_archive_path: Path | None = None,
) -> QDArchive:
    trajectory_lock = threading.Lock()

    def _append_live_trajectory(record: Mapping[str, Any]) -> None:
        with trajectory_lock:
            with trajectory_path.open("a", encoding="utf-8") as fh:
                fh.write(format_trajectory_line(record))

    def _write_trajectory() -> None:
        with trajectory_lock:
            with trajectory_path.open("w", encoding="utf-8") as fh:
                for rec in trajectory:
                    fh.write(format_trajectory_line(rec))

    def _record_result(
        *,
        idx: int,
        x: list[float],
        res: Evaluation | None,
        insert_archive: bool,
    ) -> dict[str, Any]:
        overrides = _vector_to_overrides(x, dims, base)
        if res is None or res.preflight_rejected:
            record: dict[str, Any] = {
                "eval": idx,
                "preflight_rejected": True,
                "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
            }
            if res is not None:
                record.update(trajectory_flags_from_evaluation(res))
            record["status"] = infer_trajectory_status(record)
            return record

        descriptor_details = _descriptor_details(
            res.components,
            res.metrics,
            mode=descriptor_mode,
            overrides={d.param_key: overrides.get(d.param_key) for d in dims},
        )
        d1, d2 = descriptor_details["x"], descriptor_details["y"]
        cell = (_bin_index(d1, bins), _bin_index(d2, bins))
        params = {d.param_key: float(overrides[d.param_key]) for d in dims}
        assessment = evaluate_tiers(res.components, metrics=res.metrics, config=tier_config)
        flags = trajectory_flags_from_evaluation(res)
        status = infer_trajectory_status({**flags, "components": res.components})
        improved = None
        if insert_archive and status == "gpu_ok":
            elite = Elite(
                cell=cell,
                score=res.score,
                descriptors=(d1, d2),
                params=params,
                episode=res.episode_path,
                tier=assessment.tier,
                tier_name=assessment.tier_name,
                objectives=assessment.objectives,
                descriptor_details=descriptor_details,
            )
            improved = archive.insert(elite)
        record = {
            "eval": idx,
            "score": res.score,
            "cell": list(cell),
            "descriptors": [d1, d2],
            "descriptor_details": descriptor_details,
            "tier": assessment.tier,
            "tier_name": assessment.tier_name,
            "improved": improved,
            "episode": res.episode_path,
            "components": res.components,
            "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
        }
        record.update(flags)
        record["status"] = status
        if res is not None:
            attach_pre_gpu_metrics_to_record(record, res.metrics)
        return record

    def _evaluate_batch(vectors: list[list[float]]) -> list[tuple[int, list[float], Evaluation | None] | None]:
        results: list[tuple[int, list[float], Evaluation | None] | None] = [None] * len(vectors)

        def _eval_one(i: int, x: list[float], gpu: str | None) -> None:
            with counter_lock:
                eval_counter[0] += 1
                idx = eval_counter[0]
            overrides = _vector_to_overrides(x, dims, base)
            res = qd_search.evaluate_overrides(
                overrides,
                out_dir=qd_dir,
                name=f"eval_{idx:06d}",
                example=example_cfg,
                template=tpl,
                executable=executable,
                constrained=constrained,
                phantom=phantom,
                use_preflight=use_preflight,
                cuda_devices=gpu if gpu is not None else cuda_devices,
                check_params=check_params,
                dry_run=dry_run,
                target_stop_time=target_stop_time,
                score_weights=score_weights,
                objective_mode=objective_mode,
                ftl_L=ftl_L,
                consume_plotfiles=consume_plotfiles,
                consumer_radii=consumer_radii,
                consumer_keep_last=consumer_keep_last,
                grtresna=grtresna,
                grtresna_base=grtresna_config,
                grtresna_solved_ftl_gate=grtresna_solved_ftl_gate,
                solved_ftl_gate_config=solved_ftl_gate_config,
                grtresna_convergence_config=grtresna_convergence_config,
                grtresna_postload_gate=grtresna_postload_gate,
                postload_gate_config=postload_gate_config,
            )
            results[i] = (idx, x, res)
            _append_live_trajectory(_record_result(idx=idx, x=x, res=res, insert_archive=False))

        if len(effective_gpu_ids) > 1 and not dry_run:
            threads = []
            for i, x in enumerate(vectors):
                gpu = str(effective_gpu_ids[i % len(effective_gpu_ids)])
                t = threading.Thread(target=_eval_one, args=(i, x, gpu))
                t.start()
                threads.append(t)
            for t in threads:
                t.join()
        else:
            for i, x in enumerate(vectors):
                _eval_one(i, x, cuda_devices)
        return results

    def _ingest(results: list[tuple[int, list[float], Evaluation | None] | None]) -> None:
        protected_eval_ids: set[int] = set()
        retention_events = []
        for item in results:
            if item is None:
                continue
            idx, x, res = item
            protected_eval_ids.add(idx)
            record = _record_result(idx=idx, x=x, res=res, insert_archive=True)
            if pre_gpu_enabled and near_miss_pool is not None and pre_gpu_archive is not None:
                if is_graded_rejection(record):
                    near_miss_pool.consider(record)
                    shadow = shadow_elite_from_record(
                        record,
                        bins=bins,
                        convergence_config=grtresna_convergence_config,
                    )
                    if shadow is not None:
                        pre_gpu_archive.insert(shadow)
            trajectory.append(record)
            if ftl_retention_enabled and record.get("status") == "gpu_ok":
                retention_events.extend(ftl_board.consider(record))
        _write_trajectory()
        if ftl_retention_enabled:
            append_ftl_retention_events(ftl_retention_path, retention_events)
            save_ftl_champions(ftl_champions_path, ftl_board)
        keep_ids = compute_keep_eval_ids(
            trajectory,
            keep_top_score=int(keep_top_eval_dirs),
            board=ftl_board if ftl_retention_enabled else None,
            ftl_retention_enabled=ftl_retention_enabled,
            protect_eval_ids=protected_eval_ids,
        )
        _prune_eval_dirs(
            qd_dir,
            trajectory,
            keep_eval_ids=keep_ids,
            protect_eval_ids=protected_eval_ids,
    )

    if not resume:
        seed_vectors = [
            [float(s.get(d.param_key, d.center)) for d in dims]
            for s in (seed_overrides or [])
        ]
        init_vectors = seed_vectors + [_sample_random(dims, rng) for _ in range(n_init)]
        _ingest(_evaluate_batch(init_vectors))
        write_json(qd_dir / "archive.json", archive.to_dict())

    for it in range(1, iterations + 1):
        vectors: list[list[float]] = []
        for _ in range(batch):
            if pre_gpu_enabled and pre_gpu_archive is not None and near_miss_pool is not None:
                vectors.append(
                    sample_next_candidate(
                        dims=dims,
                        rng=rng,
                        archive=archive,
                        pre_gpu_archive=pre_gpu_archive,
                        near_miss_pool=near_miss_pool,
                        config=sampling_config or QdSamplingConfig(),
                    )
                )
            else:
                from .sampling import _feasible_bounds, _mutate_elite, _sample_feasible_box
                elites = list(archive.cells.values())
                feasible_bounds = _feasible_bounds(dims, elites) if elites else None
                if elites and rng.random() < 0.85:
                    parent = elites[int(rng.integers(len(elites)))]
                    vectors.append(_mutate_elite(parent, dims, rng))
                elif feasible_bounds is not None:
                    vectors.append(_sample_feasible_box(dims, feasible_bounds, rng))
                else:
                    vectors.append(_sample_random(dims, rng))
        _ingest(_evaluate_batch(vectors))
        write_json(qd_dir / "archive.json", archive.to_dict())
        if pre_gpu_enabled and pre_gpu_archive_path is not None and pre_gpu_archive is not None:
            write_json(pre_gpu_archive_path, pre_gpu_archive.to_dict())
        print(f"[qd] legacy iter {it}/{iterations} evals={eval_counter[0]}/{total_target}")

    write_json(qd_dir / "archive.json", archive.to_dict())
    return archive
