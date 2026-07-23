"""Focused CMA-ES refinement of a geometry-atlas elite.

MAP-Elites explores diversity; this module hill-climbs frozen ``f_geo`` (with
the atlas rank score) from a seed genome before any larger archive campaign.
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

import numpy as np

from .config import GeometryAtlasConfig
from .ansatz import ANALYTIC_PARAMS
from .genome import (
    PARAMS_PER_CENTER,
    GeometryGenome,
    GeometryGenomeConfig,
    genome_bounds,
    zero_genome,
)
from .render import RenderConfig, render_and_write
from .score import evaluate_genome, evaluation_summary

logger = logging.getLogger(__name__)

ObjectiveName = Literal["score", "f_geo"]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _append_jsonl(path: Path, record: dict) -> None:
    with path.open("a", encoding="utf-8") as fout:
        fout.write(json.dumps(record) + "\n")


def _campaign_genome_config(cfg: GeometryAtlasConfig) -> GeometryGenomeConfig:
    return GeometryGenomeConfig(
        **{**asdict(cfg.genome), "box_half_width": 0.5 * cfg.render.L}
    )


def seed_alcubierre_genome(cfg: GeometryGenomeConfig) -> GeometryGenome:
    """Warm-start near a strong Alcubierre bubble (RBF coeffs near zero)."""
    g = zero_genome(cfg)
    coeffs = g.coeffs.copy()
    a0 = cfg.n_centers * PARAMS_PER_CENTER
    # Start in the trusted h-quality window (v≲1.2 on N~48); CMA can climb.
    coeffs[a0] = float(min(1.1, 0.85 * cfg.alc_velocity_max))
    coeffs[a0 + 1] = float(
        np.clip(4.0, cfg.alc_radius_min, cfg.alc_radius_max)
    )
    coeffs[a0 + 2] = float(
        np.clip(1.5, cfg.alc_sigma_min, cfg.alc_sigma_max)
    )
    coeffs[a0 + 3] = 0.0  # +x
    return GeometryGenome(coeffs=coeffs, centers=g.centers.copy(), config=cfg)


def _f_geo_objective(ev) -> float:
    """Maximise raw frozen shortcut; tiny penalty for null-H drift."""
    raw = float(ev.diagnostics.get("f_geo_raw", ev.f_geo))
    drift = float(ev.diagnostics.get("geo_max_h_rel_drift", 0.0))
    # Soft cost keeps CMA from ignoring integration quality entirely.
    return raw - 0.05 * max(0.0, drift - 1.0e-2)


def _fitness(ev, *, objective: ObjectiveName) -> float:
    """CMA-ES minimises this; rejected candidates get a large penalty."""
    if ev.rejected:
        return 1.0e6
    if objective == "f_geo":
        return -_f_geo_objective(ev)
    return -float(ev.score)


def _is_better(ev, best_value: float, *, objective: ObjectiveName) -> bool:
    if ev.rejected:
        return False
    if objective == "f_geo":
        return _f_geo_objective(ev) > best_value
    return float(ev.score) > best_value


def _best_value(ev, *, objective: ObjectiveName) -> float:
    return _f_geo_objective(ev) if objective == "f_geo" else float(ev.score)


def run_geometry_cmaes(
    cfg: GeometryAtlasConfig,
    *,
    seed_genome: GeometryGenome | None = None,
    max_evals: int | None = None,
    sigma0: float = 0.25,
    population_size: int | None = None,
    objective: ObjectiveName = "f_geo",
    alc_only: bool = False,
) -> tuple[Path, GeometryGenome, dict[str, Any]]:
    """Run CMA-ES; return ``(run_dir, best_genome, summary)``.

    ``objective='f_geo'`` maximises frozen shortcut strength directly.
    ``alc_only=True`` freezes the RBF block and optimises only the trailing
    analytic topologies (shift tube, tunnel, lens, throat).
    """
    try:
        import cma
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pycma is required for geometry CMA-ES refinement") from exc

    genome_cfg = _campaign_genome_config(cfg)
    seed = seed_genome or seed_alcubierre_genome(genome_cfg)
    if seed.config != genome_cfg:
        seed = GeometryGenome(
            coeffs=seed.coeffs, centers=seed.centers, config=genome_cfg
        )
    if alc_only:
        # Freeze RBF deformations; keep only the analytic topology tail free.
        frozen = zero_genome(genome_cfg).coeffs.copy()
        frozen[-ANALYTIC_PARAMS:] = seed.analytic_coeffs
        seed = GeometryGenome(
            coeffs=frozen, centers=seed.centers.copy(), config=genome_cfg
        )

    name = cfg.name or f"geometry_cmaes_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
    root = Path(cfg.runs_dir) / name
    root.mkdir(parents=True, exist_ok=True)
    (root / "evals").mkdir(exist_ok=True)
    traj_path = root / "trajectory.jsonl"
    _write_json(
        root / "metadata.json",
        {
            "created_at": _utc_now(),
            "config": cfg.to_dict(),
            "mode": "cmaes",
            "objective": objective,
            "alc_only": alc_only,
        },
    )

    lo_full, hi_full = genome_bounds(genome_cfg)
    n_evals = int(max_evals if max_evals is not None else cfg.target_evals)

    if alc_only:
        a0 = genome_cfg.n_centers * PARAMS_PER_CENTER
        x0 = np.clip(seed.coeffs[a0:].copy(), lo_full[a0:], hi_full[a0:])
        lo = lo_full[a0:]
        hi = hi_full[a0:]
        popsize = int(population_size) if population_size is not None else 8
    else:
        x0 = np.clip(seed.coeffs.copy(), lo_full, hi_full)
        lo, hi = lo_full, hi_full
        popsize = (
            int(population_size)
            if population_size is not None
            else max(8, min(16, x0.size // 4))
        )

    opts = cma.CMAOptions()
    opts.set("bounds", [lo.tolist(), hi.tolist()])
    opts.set("seed", int(cfg.seed))
    opts.set("maxfevals", n_evals)
    opts.set("verbose", -9)
    opts.set("verb_log", 0)
    opts.set("popsize", popsize)

    es = cma.CMAEvolutionStrategy(x0.tolist(), float(sigma0), opts)
    best_genome = seed
    best_metric = -1.0e300
    best_ev: dict[str, Any] | None = None
    eval_id = 0

    def _embed(x: np.ndarray) -> GeometryGenome:
        if alc_only:
            coeffs = seed.coeffs.copy()
            coeffs[-ANALYTIC_PARAMS:] = np.asarray(x, dtype=np.float64)
        else:
            coeffs = np.asarray(x, dtype=np.float64)
        return GeometryGenome(
            coeffs=coeffs,
            centers=seed.centers.copy(),
            config=genome_cfg,
        ).clipped()

    while not es.stop() and eval_id < n_evals:
        xs = es.ask()
        fitness: list[float] = []
        for x in xs:
            if eval_id >= n_evals:
                fitness.append(1.0e9)
                continue
            genome = _embed(x)
            ev = evaluate_genome(
                genome,
                eval_id=eval_id,
                render_cfg=cfg.render,
                out_dir=root / "evals",
                bins=cfg.bins,
                n_rays=cfg.n_rays,
                compute_ff=cfg.compute_ff,
                keep_gridinit=False,
                localise_probe=cfg.localise_probe,
            )
            _write_json(root / "evals" / f"eval_{eval_id:06d}.json", ev.to_dict())
            _append_jsonl(traj_path, evaluation_summary(ev))
            fitness.append(_fitness(ev, objective=objective))
            if _is_better(ev, best_metric, objective=objective):
                best_metric = _best_value(ev, objective=objective)
                best_genome = genome
                best_ev = evaluation_summary(ev)
                logger.info(
                    "CMA-ES new best eval=%d objective=%s value=%.6f f_geo=%.4f",
                    eval_id,
                    objective,
                    best_metric,
                    ev.f_geo,
                )
            eval_id += 1
        es.tell(xs, fitness)

    best_path = root / "best_genome.json"
    _write_json(best_path, best_genome.to_dict())
    gridinit_path = root / "best.gridinit"
    render_and_write(best_genome, gridinit_path, cfg.render)

    summary = {
        "finished_at": _utc_now(),
        "evals": eval_id,
        "objective": objective,
        "alc_only": alc_only,
        "best_metric": best_metric if best_metric > -1.0e299 else None,
        "best_score": None if best_ev is None else best_ev.get("score"),
        "best_f_geo": None if best_ev is None else best_ev.get("f_geo"),
        "best": best_ev,
        "best_genome": str(best_path),
        "best_gridinit": str(gridinit_path),
        "cma_stop": dict(es.stop()),
    }
    _write_json(root / "summary.json", summary)
    logger.info(
        "Geometry CMA-ES done: evals=%d best_f_geo=%s",
        eval_id,
        summary["best_f_geo"],
    )
    return root, best_genome, summary


__all__ = ["run_geometry_cmaes", "seed_alcubierre_genome"]
