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
from ...core.evaluation import Evaluation
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
from .io import _iterations_for_target_evals, _load_trajectory_records, _prune_eval_dirs
from .sampling import (
    _ELITE_FRACTION,
    _feasible_bounds,
    _mutate_elite,
    _sample_feasible_box,
    _sample_random,
)


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
) -> QDArchive:
    example_cfg = example if isinstance(example, ExampleConfig) else resolve_example(example)
    dims = list(search_space or DEFAULT_SEARCH_SPACE)
    tpl = template or example_cfg.template
    base = dict(base_overrides or {})
    base.setdefault("N_full", 64)
    base.setdefault("max_level", 2)
    base.setdefault("stop_time", 2.0)
    base.setdefault("plot_interval", 10)
    base.setdefault("checkpoint_interval", -1)
    base.setdefault("dt_multiplier", 0.02)
    target_stop_time = float(base["stop_time"])

    rng = np.random.default_rng(seed)
    batch = batch_size or (len(gpu_ids) if gpu_ids else 4)
    n_init = init_random if init_random is not None else batch

    if name is None:
        if resume:
            raise ValueError("resume requires --name pointing at an existing qd_<timestamp> directory")
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        name = f"qd_{timestamp}"
    qd_dir = (runs_dir / name).expanduser().resolve()

    trajectory_path = qd_dir / "trajectory.jsonl"
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
        if trajectory:
            completed_evals = max(int(rec["eval"]) for rec in trajectory if "eval" in rec)
        validation_path = qd_dir / "validation.json"
        if validation_path.exists():
            validation = json.loads(validation_path.read_text(encoding="utf-8"))
            conv_history = list(validation.get("history", []))
        metadata_path = qd_dir / "metadata.json"
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

    metadata = {
        "mode": "qd_map_elites",
        "example": example_cfg.name,
        "iterations": iterations,
        "batch_size": batch,
        "bins": bins,
        "seed": seed,
        "gpu_ids": list(gpu_ids) if gpu_ids else None,
        "descriptor_mode": descriptor_mode,
        "base_overrides": base,
        "search_space": [
            {"key": d.param_key, "lower": d.lower, "upper": d.upper} for d in dims
        ],
        "grtresna": grtresna,
        "grtresna_solved_ftl_gate": grtresna_solved_ftl_gate,
        "keep_top_eval_dirs": keep_top_eval_dirs,
        "grtresna_max_ham_pct": getattr(
            grtresna_convergence_config, "max_ham_pct", 5.0,
        ),
        "grtresna_max_mom_pct": getattr(
            grtresna_convergence_config, "max_mom_pct", 5.0,
        ),
    }
    if target_evals is not None:
        metadata["target_evals"] = target_evals
    if resume:
        metadata["resumed_at"] = datetime.now(timezone.utc).isoformat()
        metadata["resume_from_evals"] = completed_evals
        metadata["resume_planned_batches"] = iterations
        existing_meta_path = qd_dir / "metadata.json"
        if existing_meta_path.exists():
            existing = json.loads(existing_meta_path.read_text(encoding="utf-8"))
            metadata["created_at"] = existing.get(
                "created_at", datetime.now(timezone.utc).isoformat()
            )
        else:
            metadata["created_at"] = datetime.now(timezone.utc).isoformat()
    else:
        metadata["created_at"] = datetime.now(timezone.utc).isoformat()
    write_json(qd_dir / "metadata.json", metadata)

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
            res.components, res.metrics, mode=descriptor_mode,
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
            _append_live_trajectory(
                _record_result(idx=idx, x=x, res=res, insert_archive=False)
            )

        if gpu_ids and len(gpu_ids) > 1 and not dry_run:
            threads = []
            for i, x in enumerate(vectors):
                gpu = str(gpu_ids[i % len(gpu_ids)])
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
        for item in results:
            if item is None:
                continue
            idx, x, res = item
            protected_eval_ids.add(idx)
            trajectory.append(
                _record_result(idx=idx, x=x, res=res, insert_archive=True)
            )
        _write_trajectory()
        _prune_eval_dirs(
            qd_dir,
            trajectory,
            keep_top=int(keep_top_eval_dirs),
            protect_eval_ids=protected_eval_ids,
        )

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

    start_iter = len(conv_history)
    n_seed = 0 if resume else len(seed_overrides or [])
    total_target = (
        target_evals
        if target_evals is not None
        else completed_evals + (0 if resume else n_init) + n_seed + iterations * batch
    )
    print(
        f"[qd] MAP-Elites: {len(dims)}D, bins={bins}x{bins}, batch={batch}, "
        f"iters={iterations}, GPUs={gpu_ids or cuda_devices}, "
        f"evals={completed_evals}->{total_target}"
    )

    if not resume:
        seed_vectors = [
            [float(s.get(d.param_key, d.center)) for d in dims]
            for s in (seed_overrides or [])
        ]
        init_vectors = seed_vectors + [
            _sample_random(dims, rng) for _ in range(n_init)
        ]
        _ingest(_evaluate_batch(init_vectors))
        write_json(qd_dir / "archive.json", archive.to_dict())
        _write_validation()

    for it in range(1, iterations + 1):
        vectors: list[list[float]] = []
        elites = list(archive.cells.values())
        feasible_bounds = _feasible_bounds(dims, elites) if elites else None
        for _ in range(batch):
            if elites and rng.random() < _ELITE_FRACTION:
                parent = elites[int(rng.integers(len(elites)))]
                vectors.append(_mutate_elite(parent, dims, rng))
            elif feasible_bounds is not None:
                vectors.append(_sample_feasible_box(dims, feasible_bounds, rng))
            else:
                vectors.append(_sample_random(dims, rng))
        _ingest(_evaluate_batch(vectors))

        write_json(qd_dir / "archive.json", archive.to_dict())
        signals = _write_validation()
        best = archive.best
        n_front = len(conv_history[-1]["front_labels"])
        global_iter = start_iter + it
        print(
            f"[qd] iter {it}/{iterations} (total {global_iter}): "
            f"elites={len(archive.cells)} coverage={archive.coverage:.2f} "
            f"best={best.score if best else float('nan'):.4f} "
            f"front(T>={survivor_min_tier})={n_front} "
            f"converged={signals.get('converged')} evals={eval_counter[0]}/{total_target}"
        )
        if signals.get("converged"):
            print(f"[qd] archive converged at iter {it}: "
                  f"coverage_delta={signals.get('coverage_delta'):.3f} "
                  f"best_delta={signals.get('best_score_delta'):.4f} "
                  f"front_stable={signals.get('front_stable')}")

    write_json(qd_dir / "archive.json", archive.to_dict())
    final_hist = archive.tier_histogram()
    print(f"[qd] Done. elites={len(archive.cells)} coverage={archive.coverage:.2f} "
          f"tiers={final_hist} dir={qd_dir}")
    return archive
