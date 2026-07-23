"""MAP-Elites driver for the pure-geometry stationary atlas."""

from __future__ import annotations

import json
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from ..qd_search.archive import Elite, QDArchive
from ..validation_tiers import Tier, TIER_NAMES
from .config import GeometryAtlasConfig
from .genome import GeometryGenome, GeometryGenomeConfig, mutate_genome, sample_genome, zero_genome
from .refine import seed_alcubierre_genome
from .render import RenderConfig, render_and_write
from .score import GeometryAtlasEvaluation, evaluate_genome, evaluation_summary

logger = logging.getLogger(__name__)

# Stationary scout thresholds (screening only; not evolved F_op).
_NONTRIVIAL_F_GEO = 1.0e-3
_NONTRIVIAL_F_FF = 1.0e-3


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _make_run_dir(cfg: GeometryAtlasConfig) -> Path:
    name = cfg.name or f"geometry_atlas_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
    root = Path(cfg.runs_dir) / name
    root.mkdir(parents=True, exist_ok=True)
    (root / "evals").mkdir(exist_ok=True)
    (root / "elites").mkdir(exist_ok=True)
    return root


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _append_jsonl(path: Path, record: dict) -> None:
    with path.open("a", encoding="utf-8") as fout:
        fout.write(json.dumps(record) + "\n")


def _campaign_genome_config(cfg: GeometryAtlasConfig) -> GeometryGenomeConfig:
    """Align genome box metadata with the campaign render domain."""
    return GeometryGenomeConfig(
        **{**asdict(cfg.genome), "box_half_width": 0.5 * cfg.render.L}
    )


def _with_campaign_genome(
    genome: GeometryGenome, cfg: GeometryAtlasConfig
) -> GeometryGenome:
    return GeometryGenome(
        coeffs=genome.coeffs,
        centers=genome.centers,
        config=_campaign_genome_config(cfg),
    )


def _tier_from_evaluation(ev: GeometryAtlasEvaluation) -> tuple[int, str]:
    """Map Stage-1 screening outcomes onto the shared tier ladder.

    Geometry atlas never claims evolved persistence: at most NONTRIVIAL for a
    measurable stationary shortcut signal, else CONSTRUCTED for a valid metric.
    """
    if ev.rejected:
        return int(Tier.REJECTED), TIER_NAMES[int(Tier.REJECTED)]
    nontrivial = (ev.f_geo >= _NONTRIVIAL_F_GEO) or (
        ev.freefall_reached and ev.f_ff >= _NONTRIVIAL_F_FF
    )
    tier = Tier.NONTRIVIAL if nontrivial else Tier.CONSTRUCTED
    return int(tier), TIER_NAMES[int(tier)]


def _elite_from_evaluation(
    ev: GeometryAtlasEvaluation, *, episode: str | None = None
) -> Elite:
    tier, tier_name = _tier_from_evaluation(ev)
    return Elite(
        cell=ev.cell,
        score=ev.score,
        descriptors=ev.descriptors,
        params={"eval_id": float(ev.eval_id)},
        episode=episode if episode is not None else ev.gridinit_path,
        tier=tier,
        tier_name=tier_name,
        objectives={
            "f_geo": ev.f_geo,
            "f_ff": ev.f_ff,
            "integral_negative_rho": ev.integral_negative_rho,
            "min_rho": ev.min_rho,
        },
        descriptor_details={
            "shift_fraction": ev.descriptors[0],
            "log_exotic_energy": ev.descriptors[1],
        },
    )


def _genome_from_elite_payload(payload: dict) -> GeometryGenome:
    return GeometryGenome.from_dict(payload["genome"])


def _load_warm_start_genomes(cfg: GeometryAtlasConfig) -> list[GeometryGenome]:
    """Load JSON genomes (e.g. a CMA elite) to inject into the archive at start."""
    loaded: list[GeometryGenome] = []
    for gpath in cfg.warm_start_genomes:
        try:
            payload = _load_json(Path(gpath))
            gdata = payload.get("genome", payload) if isinstance(payload, dict) else payload
            loaded.append(_with_campaign_genome(GeometryGenome.from_dict(gdata), cfg))
        except Exception as exc:  # noqa: BLE001 — bad warm-start must not abort the run
            logger.warning("Skipping warm-start genome %s: %s", gpath, exc)
    return loaded


def run_geometry_atlas(cfg: GeometryAtlasConfig) -> tuple[Path, QDArchive, dict[str, Any]]:
    """Run (or resume) a pure-geometry MAP-Elites campaign.

    Returns ``(run_dir, archive, summary)``.
    """
    root = _make_run_dir(cfg)
    state_path = root / "state.json"
    archive_path = root / "archive.json"
    traj_path = root / "trajectory.jsonl"
    meta_path = root / "metadata.json"

    rng = np.random.default_rng(cfg.seed)
    archive = QDArchive(bins=cfg.bins)
    next_eval_id = 0
    elite_genomes: dict[tuple[int, int], GeometryGenome] = {}

    if cfg.resume and state_path.exists():
        state = _load_json(state_path)
        next_eval_id = int(state.get("next_eval_id", 0))
        if archive_path.exists():
            archive = QDArchive.from_dict(_load_json(archive_path))
        # Reload elite genomes from per-eval JSON.
        for cell_key, elite in list(archive.cells.items()):
            eval_id = int(elite.params.get("eval_id", -1))
            eval_json = root / "evals" / f"eval_{eval_id:06d}.json"
            if eval_json.exists():
                payload = _load_json(eval_json)
                if payload.get("genome"):
                    elite_genomes[cell_key] = _genome_from_elite_payload(payload)
        logger.info(
            "Resumed geometry atlas at eval_id=%d with %d elites",
            next_eval_id,
            len(archive.cells),
        )
    else:
        _write_json(meta_path, {"created_at": _utc_now(), "config": cfg.to_dict()})

    # Seed with Minkowski (flat baseline) and a moderate Alcubierre warm-start.
    if next_eval_id == 0:
        seeds: list[GeometryGenome] = [_with_campaign_genome(zero_genome(cfg.genome), cfg)]
        if cfg.genome.enable_alcubierre:
            seeds.append(
                _with_campaign_genome(seed_alcubierre_genome(_campaign_genome_config(cfg)), cfg)
            )
        seeds.extend(_load_warm_start_genomes(cfg))
        for seed in seeds:
            if next_eval_id >= cfg.target_evals:
                break
            ev = _evaluate_and_store(
                seed,
                eval_id=next_eval_id,
                cfg=cfg,
                root=root,
                keep_gridinit=True,
            )
            next_eval_id += 1
            if not ev.rejected:
                elite_path = _elite_gridinit_path(root, ev.cell)
                if archive.insert(_elite_from_evaluation(ev, episode=str(elite_path))):
                    elite_genomes[ev.cell] = seed
                    _persist_elite_gridinit(ev, seed, cfg.render, root)
            _append_jsonl(traj_path, evaluation_summary(ev))
        _checkpoint(root, archive, next_eval_id, cfg)

    while next_eval_id < cfg.target_evals:
        batch_n = min(cfg.batch_size, cfg.target_evals - next_eval_id)
        genomes: list[tuple[int, GeometryGenome]] = []
        for _ in range(batch_n):
            eval_id = next_eval_id
            next_eval_id += 1
            if archive.cells and rng.random() > cfg.random_fraction:
                # Mutate a random elite.
                cell = list(archive.cells.keys())[int(rng.integers(0, len(archive.cells)))]
                parent = elite_genomes.get(cell)
                if parent is None:
                    parent = sample_genome(rng, cfg.genome)
                child = mutate_genome(parent, rng)
            else:
                child = sample_genome(rng, cfg.genome)
            genomes.append((eval_id, _with_campaign_genome(child, cfg)))

        results = _evaluate_batch(genomes, cfg, root)
        for ev, genome in results:
            _append_jsonl(traj_path, evaluation_summary(ev))
            if ev.rejected:
                continue
            elite_path = _elite_gridinit_path(root, ev.cell)
            elite = _elite_from_evaluation(ev, episode=str(elite_path))
            if archive.insert(elite):
                elite_genomes[ev.cell] = genome
                _persist_elite_gridinit(ev, genome, cfg.render, root)
                logger.info(
                    "New elite cell=%s score=%.4f f_geo=%.4f f_ff=%.4f E-=%.4g",
                    ev.cell,
                    ev.score,
                    ev.f_geo,
                    ev.f_ff,
                    ev.integral_negative_rho,
                )
        _checkpoint(root, archive, next_eval_id, cfg)

    qd_score_f_geo = float(
        sum(float(e.objectives.get("f_geo", 0.0)) for e in archive.cells.values())
    )
    best_f_geo = max(
        (float(e.objectives.get("f_geo", 0.0)) for e in archive.cells.values()),
        default=0.0,
    )
    summary = {
        "finished_at": _utc_now(),
        "evals": next_eval_id,
        "coverage": archive.coverage,
        "num_elites": len(archive.cells),
        "qd_score_f_geo": qd_score_f_geo,
        "best_f_geo": best_f_geo,
        "best_score": archive.best.score if archive.best else None,
        "best_cell": list(archive.best.cell) if archive.best else None,
        "best_objectives": archive.best.objectives if archive.best else None,
    }
    _write_json(root / "summary.json", summary)
    logger.info(
        "Geometry atlas done: evals=%d elites=%d coverage=%.2f",
        next_eval_id,
        len(archive.cells),
        archive.coverage,
    )
    return root, archive, summary


def _evaluate_and_store(
    genome: GeometryGenome,
    *,
    eval_id: int,
    cfg: GeometryAtlasConfig,
    root: Path,
    keep_gridinit: bool,
) -> GeometryAtlasEvaluation:
    ev = evaluate_genome(
        genome,
        eval_id=eval_id,
        render_cfg=cfg.render,
        out_dir=root / "evals",
        bins=cfg.bins,
        n_rays=cfg.n_rays,
        compute_ff=cfg.compute_ff,
        keep_gridinit=keep_gridinit,
        localise_probe=cfg.localise_probe,
    )
    _write_json(root / "evals" / f"eval_{eval_id:06d}.json", ev.to_dict())
    return ev


def _evaluate_batch(
    genomes: list[tuple[int, GeometryGenome]],
    cfg: GeometryAtlasConfig,
    root: Path,
) -> list[tuple[GeometryAtlasEvaluation, GeometryGenome]]:
    """Evaluate a batch, optionally in parallel threads (CPU-bound numpy)."""
    results: list[tuple[GeometryAtlasEvaluation, GeometryGenome]] = []
    # Einstein-source evaluation is memory-heavy; keep concurrency modest.
    workers = max(1, min(cfg.batch_size, 2))

    def _one(item: tuple[int, GeometryGenome]) -> tuple[GeometryAtlasEvaluation, GeometryGenome]:
        eval_id, genome = item
        ev = _evaluate_and_store(
            genome,
            eval_id=eval_id,
            cfg=cfg,
            root=root,
            keep_gridinit=False,
        )
        return ev, genome

    if workers == 1 or len(genomes) == 1:
        for item in genomes:
            results.append(_one(item))
        return results

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(_one, item) for item in genomes]
        for fut in as_completed(futs):
            results.append(fut.result())
    results.sort(key=lambda pair: pair[0].eval_id)
    return results


def _elite_gridinit_path(root: Path, cell: tuple[int, int]) -> Path:
    return root / "elites" / f"cell_{cell[0]}_{cell[1]}.gridinit"


def _persist_elite_gridinit(
    ev: GeometryAtlasEvaluation,
    genome: GeometryGenome,
    render_cfg: RenderConfig,
    root: Path,
) -> Path:
    """Write / refresh the elite cell's gridinit for Stage-2 handoff."""
    path = _elite_gridinit_path(root, ev.cell)
    cell_tag = path.stem
    render_and_write(genome, path, render_cfg)
    meta = {
        "cell": list(ev.cell),
        "eval_id": ev.eval_id,
        "score": ev.score,
        "f_geo": ev.f_geo,
        "f_ff": ev.f_ff,
        "integral_negative_rho": ev.integral_negative_rho,
        "gridinit_path": str(path),
        "genome": genome.to_dict(),
        "diagnostics": ev.diagnostics,
    }
    _write_json(root / "elites" / f"{cell_tag}.json", meta)
    # Also update the evaluation record with the elite gridinit path.
    eval_json = root / "evals" / f"eval_{ev.eval_id:06d}.json"
    if eval_json.exists():
        payload = _load_json(eval_json)
        payload["gridinit_path"] = str(path)
        _write_json(eval_json, payload)
    return path


def _checkpoint(
    root: Path,
    archive: QDArchive,
    next_eval_id: int,
    cfg: GeometryAtlasConfig,
) -> None:
    _write_json(root / "archive.json", archive.to_dict())
    _write_json(
        root / "state.json",
        {
            "next_eval_id": next_eval_id,
            "updated_at": _utc_now(),
            "bins": cfg.bins,
            "target_evals": cfg.target_evals,
            "seed": cfg.seed,
        },
    )


__all__ = ["run_geometry_atlas"]
