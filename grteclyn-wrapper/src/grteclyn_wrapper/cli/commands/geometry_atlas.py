"""CLI command: pure-geometry MAP-Elites atlas / calibration / CMA-ES refine."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ...search.geometry_atlas import (
    GeometryAtlasConfig,
    GeometryGenome,
    GeometryGenomeConfig,
    RenderConfig,
    calibrate_atlas_probe,
    run_geometry_atlas,
    run_geometry_cmaes,
    seed_alcubierre_genome,
)
from ...search.geometry_atlas.calibrate import DEFAULT_CAND146_GRIDINIT


def _genome_config(args: argparse.Namespace) -> GeometryGenomeConfig:
    return GeometryGenomeConfig(
        n_centers=int(args.n_centers),
        support_radius=float(args.support_radius),
        rbf_width=float(args.rbf_width),
        alpha_amp=float(args.alpha_amp),
        beta_amp=float(args.beta_amp),
        log_metric_amp=float(args.log_metric_amp),
        kij_amp=float(args.kij_amp),
        mutation_sigma=float(args.mutation_sigma),
        box_half_width=0.5 * float(args.L),
        alc_velocity_max=float(args.alc_velocity_max),
        enable_alcubierre=not bool(args.no_alcubierre),
    )


def _atlas_config(args: argparse.Namespace) -> GeometryAtlasConfig:
    return GeometryAtlasConfig(
        runs_dir=Path(args.runs_dir),
        name=args.name,
        target_evals=int(args.target_evals),
        bins=int(args.bins),
        seed=int(args.seed),
        batch_size=int(args.batch_size),
        n_rays=int(args.n_rays),
        compute_ff=not bool(args.no_ff),
        resume=bool(args.resume),
        random_fraction=float(args.random_fraction),
        localise_probe=not bool(args.fullbox_probe),
        warm_start_genomes=tuple(
            p.strip() for p in str(args.seed_genome or "").split(",") if p.strip()
        ),
        genome=_genome_config(args),
        render=RenderConfig(n=int(args.n), L=float(args.L)),
    )


def run_geometry_atlas_command(args: argparse.Namespace) -> int:
    mode = str(getattr(args, "mode", "map_elites") or "map_elites")
    cfg = _atlas_config(args)

    if mode == "calibrate":
        out = Path(args.runs_dir) / (args.name or "geometry_atlas_calibration")
        if args.cand146 == "":
            cand146 = DEFAULT_CAND146_GRIDINIT
        elif args.cand146.lower() in {"none", "skip"}:
            cand146 = None
        else:
            cand146 = Path(args.cand146)
        report = calibrate_atlas_probe(
            out,
            n_rays=int(args.n_rays),
            alc_n=int(args.calibrate_n),
            alc_L=float(args.calibrate_L),
            cand146_path=cand146,
            localise=not bool(args.fullbox_probe),
        )
        print(json.dumps({"calibration": str(out), "report": report}, indent=2))
        return 0 if report["verdict"]["alcubierre_localised_ok"] else 2

    if mode == "cmaes":
        seed = None
        if args.seed_genome:
            seed = GeometryGenome.from_dict(
                json.loads(Path(args.seed_genome).read_text(encoding="utf-8"))
            )
        elif not args.no_alcubierre:
            seed = seed_alcubierre_genome(_genome_config(args))
        root, _best, summary = run_geometry_cmaes(
            cfg,
            seed_genome=seed,
            max_evals=int(args.target_evals),
            sigma0=float(args.cma_sigma0),
            population_size=int(args.cma_popsize) if args.cma_popsize else None,
            objective=str(args.cma_objective),
            alc_only=bool(args.alc_only),
        )
        print(json.dumps({"geometry_cmaes": str(root), "summary": summary}, indent=2))
        return 0

    root, archive, summary = run_geometry_atlas(cfg)
    print(
        json.dumps(
            {
                "geometry_atlas": str(root),
                "summary": summary,
                "archive_elites": len(archive.cells),
            },
            indent=2,
        )
    )
    return 0
