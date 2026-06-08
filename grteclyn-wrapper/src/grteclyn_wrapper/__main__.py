from __future__ import annotations

import argparse
import json
import os
import random
from pathlib import Path
from typing import Any, Mapping

from .core.config import default_runs_dir, resolve_example, resolve_executable
from .core.episode import create_episode, update_metadata, write_json
from .core.params import write_params
from .core.runner import run_episode
from .initial_data.candidates import resolve_initial_data_overrides
from .initial_data.constrained_recipe import constrained_overrides
from .initial_data.seeds import get_seed, list_seeds
from .initial_data.validate_guesser import run_validation
from .metrics.episode_metrics import dataclass_to_dict, read_episode_metrics
from .metrics.score import score_episode
from .search.atlas import run_atlas
from .search.optimize import run_optimize
from .projection.postload_gate import PostLoadGateConfig
from .search.grtresna_convergence_gate import GRTresnaConvergenceConfig
from .search.solved_ftl_gate import SolvedFtlGateConfig


def _grtresna_postload_gate_enabled(args: argparse.Namespace, *, use_grtresna: bool) -> bool:
    if not use_grtresna:
        return False
    if getattr(args, "grtresna_postload_gate", False):
        return True
    return os.environ.get("POSTLOAD_GATE", "0") == "1"


def _postload_gate_config_from_args(args: argparse.Namespace) -> PostLoadGateConfig:
    return PostLoadGateConfig(
        max_hamiltonian_l2=getattr(args, "postload_max_ham_l2", 1.0e-2),
        max_momentum_l2=getattr(args, "postload_max_mom_l2", 1.0e-2),
    )


SWEEP_RANGES = {
    "wormhole_phi_perturbation_amplitude": (-0.04, 0.04),
    "wormhole_support_strength": (0.2, 1.0),
    "wormhole_phi_perturbation_width": (0.25, 1.0),
}


def _parse_override(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError(f"Override must be key=value, got {value!r}")
    key, raw = value.split("=", 1)
    key = key.strip()
    raw = raw.strip()
    if not key:
        raise argparse.ArgumentTypeError("Override key cannot be empty")
    return key, raw


def _coerce_value(raw: str) -> Any:
    for caster in (int, float):
        try:
            return caster(raw)
        except ValueError:
            pass
    if raw.lower() in {"true", "false"}:
        return raw.lower() == "true"
    return raw


def _collect_overrides(pairs: list[tuple[str, str]]) -> dict[str, Any]:
    return {key: _coerce_value(value) for key, value in pairs}


def _parse_score_weights(pairs: list[tuple[str, str]]) -> dict[str, float]:
    weights: dict[str, float] = {}
    for key, raw in pairs:
        weights[key] = float(raw)
    return weights


def _finalize_score(
    episode_dir: Path,
    target_stop_time: float | None,
    *,
    score_weights: Mapping[str, float] | None = None,
    ftl_L: float | None = None,
) -> int:
    metrics = read_episode_metrics(episode_dir, ftl_L=ftl_L)
    score = score_episode(metrics, target_stop_time=target_stop_time, weights=score_weights)
    write_json(
        episode_dir / "score.json",
        {
            "score": dataclass_to_dict(score),
            "metrics": dataclass_to_dict(metrics),
        },
    )
    print(json.dumps({"episode": str(episode_dir), "score": score.total}, indent=2))
    return 0


def _run_single(args: argparse.Namespace, overrides: dict[str, Any]) -> int:
    example = resolve_example(args.example)
    phantom = getattr(args, "phantom", False) or getattr(args, "phantom_default", False)
    if getattr(args, "constrained", False) and example.name == "RadialRecipe":
        constrained_overrides(overrides, phantom=phantom)
    # Evolve exotic (phantom) matter when the recipe was solved for it, so the
    # geometry is faithfully sourced instead of mismatched canonical matter.
    if phantom and example.name == "RadialRecipe":
        overrides.setdefault("recipe_exotic_matter", 1)
    runs_dir = Path(args.runs_dir).expanduser().resolve()
    metadata: dict[str, Any] = {
        "mode": args.command,
        "example": example.name,
        "overrides": overrides,
    }
    if getattr(args, "initial_data_source", None):
        metadata["initial_data_source"] = args.initial_data_source
    episode = create_episode(
        runs_dir,
        name=args.name,
        metadata=metadata,
    )
    template = Path(args.template).expanduser().resolve() if args.template else example.template
    write_params(
        template,
        episode.params_path,
        episode_dir=episode.path,
        example=example,
        overrides=overrides,
    )

    if args.dry_run:
        update_metadata(episode, {"dry_run": True})
        print(f"Wrote dry-run episode: {episode.path}")
        return 0

    executable = resolve_executable(
        args.executable,
        example=example,
        mpi_ranks=args.mpi_ranks,
        comp=args.comp,
        cuda=not args.no_cuda,
        debug=args.debug,
    )
    result = run_episode(
        episode,
        executable,
        check_params=not args.skip_check_params,
        cuda_devices=args.cuda_devices,
        consume_plotfiles=args.consume_plotfiles,
        consumer_radii=args.consumer_radii,
        consumer_delete=args.consumer_delete,
        consumer_keep_last=getattr(args, "consumer_keep_last", 1),
    )
    _finalize_score(
        episode.path,
        target_stop_time=overrides.get("stop_time"),
        score_weights=getattr(args, "score_weights", None),
        ftl_L=getattr(args, "ftl_L", None),
    )
    return result.returncode


def _sample_overrides(base: dict[str, Any], rng: random.Random) -> dict[str, Any]:
    overrides = dict(base)
    for key, (lo, hi) in SWEEP_RANGES.items():
        overrides.setdefault(key, rng.uniform(lo, hi))
    return overrides


def _run_sweep(args: argparse.Namespace, base_overrides: dict[str, Any]) -> int:
    rng = random.Random(args.seed)
    status = 0
    for index in range(args.count):
        overrides = _sample_overrides(base_overrides, rng)
        name = args.name or f"sweep_{index + 1:06d}"
        per_args = argparse.Namespace(**vars(args))
        per_args.name = name
        per_args.command = "sweep"
        rc = _run_single(per_args, overrides)
        status = status or rc
        if rc != 0 and args.stop_on_failure:
            break
    return status


def _run_optimize_command(args: argparse.Namespace, base_overrides: dict[str, Any]) -> int:
    # The GRTresna gridinit loader (ExternalGridInitialData) only exists in the
    # RadialRecipe example. If --grtresna is set but the example was left at the
    # default, switch to RadialRecipe so the loaded initial data is actually used.
    if getattr(args, "grtresna", False) and args.example == "SupportedWormholeCollapse":
        print("[optimize] --grtresna requires the RadialRecipe example; "
              "switching --example to RadialRecipe.")
        args.example = "RadialRecipe"
    example = resolve_example(args.example)
    executable = None
    if not args.dry_run:
        executable = resolve_executable(
            args.executable,
            example=example,
            mpi_ranks=args.mpi_ranks,
            comp=args.comp,
            cuda=not args.no_cuda,
            debug=args.debug,
        )

    template = Path(args.template).expanduser().resolve() if args.template else None

    # Non-spherical search: extend the radial space with gauge (lapse/shift)
    # angular modes and activate them via the constant base overrides.
    from .search.optimize import build_search_space, ANGULAR_BASE_OVERRIDES
    nonspherical = getattr(args, "nonspherical", False)
    use_grtresna = getattr(args, "grtresna", False)
    grtresna_lumps = getattr(args, "grtresna_lumps", 5)
    grtresna_ansatz = getattr(args, "grtresna_ansatz", "free")
    grtresna_shell_profile = getattr(args, "grtresna_shell_profile", "compact")
    search_space = build_search_space(
        nonspherical=nonspherical, grtresna=use_grtresna,
        grtresna_lumps=grtresna_lumps,
        grtresna_ansatz=grtresna_ansatz,
        grtresna_shell_profile=grtresna_shell_profile,
    )
    if use_grtresna and grtresna_ansatz == "ring":
        base_overrides = {
            **base_overrides,
            "grtresna_ring_lumps": grtresna_lumps,
        }
    if use_grtresna and grtresna_ansatz == "shell":
        base_overrides = {
            **base_overrides,
            "grtresna_shell_lumps": grtresna_lumps,
        }
    grtresna_full_z = bool(getattr(args, "grtresna_full_z", False))
    grtresna_domain = None
    if use_grtresna:
        from .grtresna.domain import GRTresnaDomainConfig

        grtresna_domain = GRTresnaDomainConfig(
            full_z=grtresna_full_z,
            l_full=getattr(args, "grtresna_evolution_l_full", 64.0),
            n_full=getattr(args, "grtresna_evolution_n_full", 64),
            grtresna_l=getattr(args, "grtresna_domain_l", 128.0),
            grtresna_nx=getattr(args, "grtresna_domain_nx", 64),
            grtresna_ny=getattr(args, "grtresna_domain_ny", 64),
            grtresna_nz=getattr(args, "grtresna_domain_nz", None),
            gridinit_nx=getattr(args, "grtresna_gridinit_nx", 64),
            gridinit_ny=getattr(args, "grtresna_gridinit_ny", 64),
            gridinit_nz=getattr(args, "grtresna_gridinit_nz", 64),
        )
        base_overrides = {**base_overrides, **grtresna_domain.evolution_overrides()}
    if nonspherical and not use_grtresna:
        base_overrides = {**ANGULAR_BASE_OVERRIDES, **base_overrides}

    # GRTresna-in-the-loop: build the fixed solver config (grid / MPI / cleanup
    # / BH content); the searched grtresna_lump_* dims are layered per-candidate.
    grtresna_config = None
    if use_grtresna:
        from .grtresna.solver import GRTresnaConfig

        grtresna_config = GRTresnaConfig(
            mpi_ranks=getattr(args, "grtresna_ranks", 8),
            max_NL_iterations=getattr(args, "grtresna_iterations", 50),
            timeout=getattr(args, "grtresna_timeout", 3600),
            max_level=getattr(args, "grtresna_max_level", 3),
            refine_threshold=getattr(args, "grtresna_refine_threshold", 0.5),
            regrid_radius=getattr(args, "grtresna_regrid_radius", 0.0),
            coefficient_average_type=getattr(
                args, "grtresna_coefficient_average_type", "harmonic"
            ),
            psi_relaxation=getattr(args, "grtresna_psi_relaxation", 1.0),
            psi_floor=getattr(args, "grtresna_psi_floor", -1.0),
            maximal_jacobian_cap=getattr(args, "grtresna_jacobian_cap", -1.0),
            # The searched lump basis is the matter source; leaving the legacy
            # radial scalar profile on mixes in an unrelated canonical cloud.
            dphi=0.0,
            dpi=0.0,
            # focus on moving/rotating MATTER: no black holes by default
            bh1_bare_mass=0.0,
            bh1_spin=(0.0, 0.0, 0.0),
            cleanup=True,
        )
        if grtresna_domain is not None:
            grtresna_config = grtresna_domain.apply_to_solver(grtresna_config)

    x0 = None
    if getattr(args, "seed_name", None):
        seed_obj = get_seed(args.seed_name)
        x0 = []
        for dim in search_space:
            seed_value = float(seed_obj.overrides.get(dim.param_key, dim.center))
            x0.append(max(dim.lower, min(dim.upper, seed_value)))

    gpu_ids = None
    if getattr(args, "gpu_ids", None):
        gpu_ids = args.gpu_ids

    use_constrained = getattr(args, "constrained", False)
    use_phantom = getattr(args, "phantom", False)
    use_preflight = getattr(args, "preflight", False)
    if args.command == "optimize":
        use_constrained = True if not hasattr(args, "_no_constrained") else use_constrained
        # The search defaults to phantom-supported constrained data, but
        # --no-phantom switches to a *normal-matter* (rho >= 0) constrained
        # solve so the search targets FTL achievable WITHOUT exotic matter:
        # geometries that demand negative energy then fail the Hamiltonian
        # constraint and are pruned / penalized instead of being force-fed
        # phantom support.
        use_phantom = not getattr(args, "no_phantom", False)
        use_preflight = True if not hasattr(args, "_no_preflight") else use_preflight

    solved_ftl_gate_config = SolvedFtlGateConfig(
        f_op_floor=getattr(args, "solved_ftl_f_op_floor", 1.0e-4),
        near_luminal_speed_floor=getattr(
            args, "solved_ftl_near_luminal_speed_floor", 0.99
        ),
        superluminal_speed_floor=getattr(
            args, "solved_ftl_superluminal_speed_floor", 1.01
        ),
        superluminal_fraction_floor=getattr(
            args, "solved_ftl_superluminal_fraction_floor", 0.02
        ),
        max_physical_coord_speed=getattr(
            args, "solved_ftl_max_physical_coord_speed", 8.0
        ),
        max_physical_f_op=getattr(args, "solved_ftl_max_physical_f_op", 0.85),
        rejection_speed_target=getattr(
            args, "solved_ftl_rejection_speed_target", 1.01
        ),
    )
    grtresna_convergence_config = (
        GRTresnaConvergenceConfig(
            max_ham_pct=getattr(args, "grtresna_max_ham_pct", 5.0),
            max_mom_pct=getattr(args, "grtresna_max_mom_pct", 5.0),
        )
        if use_grtresna
        else None
    )

    result = run_optimize(
        search_space=search_space,
        runs_dir=Path(args.runs_dir).expanduser().resolve(),
        executable=executable,
        max_generations=args.max_generations,
        population_size=args.population_size,
        sigma0=args.sigma0,
        seed=args.seed,
        base_overrides=base_overrides,
        template=template,
        example=example,
        name=args.name,
        dry_run=args.dry_run,
        constrained=use_constrained,
        phantom=use_phantom,
        use_preflight=use_preflight,
        cuda_devices=args.cuda_devices,
        gpu_ids=gpu_ids,
        check_params=not args.skip_check_params,
        x0=x0,
        consume_plotfiles=not getattr(args, "no_consume_plotfiles", False),
        consumer_radii=getattr(args, "consumer_radii", [4.0, 8.0]),
        consumer_keep_last=getattr(args, "consumer_keep_last", 1),
        score_weights=getattr(args, "score_weights", None),
        objective_mode=getattr(args, "objective_mode", "weighted"),
        ftl_L=getattr(args, "ftl_L", None),
        surrogate=getattr(args, "surrogate", False),
        surrogate_keep_fraction=getattr(args, "surrogate_keep_fraction", 0.5),
        grtresna=use_grtresna,
        grtresna_config=grtresna_config,
        solved_ftl_gate_config=solved_ftl_gate_config,
        grtresna_convergence_config=grtresna_convergence_config,
        grtresna_postload_gate=_grtresna_postload_gate_enabled(args, use_grtresna=use_grtresna),
        postload_gate_config=(
            _postload_gate_config_from_args(args) if use_grtresna else None
        ),
        warm_start_trajectories=[
            Path(p).expanduser().resolve()
            for p in getattr(args, "warm_start_trajectory", [])
        ],
        warm_start_top_k=getattr(args, "warm_start_top_k", 8),
        warm_start_jitter=getattr(args, "warm_start_jitter", 0.08),
        random_injection_fraction=getattr(args, "random_injection_fraction", 0.0),
        exotic_injection_fraction=getattr(args, "exotic_injection_fraction", 0.0),
    )
    print(json.dumps({
        "best_score": result.best_score,
        "best_episode": result.best_episode,
        "best_params": result.best_params,
        "generations": result.generations,
        "evaluations": result.evaluations,
    }, indent=2))
    return 0


def _run_qd_command(args: argparse.Namespace, base_overrides: dict[str, Any]) -> int:
    from .search.qd_search import run_qd_search

    if getattr(args, "grtresna", False) and args.example == "SupportedWormholeCollapse":
        print("[qd] --grtresna requires the RadialRecipe example; switching --example to RadialRecipe.")
        args.example = "RadialRecipe"
    example = resolve_example(args.example)
    executable = None
    if not args.dry_run:
        executable = resolve_executable(
            args.executable,
            example=example,
            mpi_ranks=args.mpi_ranks,
            comp=args.comp,
            cuda=not args.no_cuda,
            debug=args.debug,
        )
    template = Path(args.template).expanduser().resolve() if args.template else None

    from .search.optimize import ANGULAR_BASE_OVERRIDES, build_search_space
    nonspherical = getattr(args, "nonspherical", False)
    use_grtresna = getattr(args, "grtresna", False)
    grtresna_lumps = getattr(args, "grtresna_lumps", 5)
    grtresna_ansatz = getattr(args, "grtresna_ansatz", "free")
    grtresna_shell_profile = getattr(args, "grtresna_shell_profile", "compact")
    search_space = build_search_space(
        nonspherical=nonspherical,
        grtresna=use_grtresna,
        grtresna_lumps=grtresna_lumps,
        grtresna_ansatz=grtresna_ansatz,
        grtresna_shell_profile=grtresna_shell_profile,
    )
    if use_grtresna and grtresna_ansatz == "ring":
        base_overrides = {**base_overrides, "grtresna_ring_lumps": grtresna_lumps}
    if use_grtresna and grtresna_ansatz == "shell":
        base_overrides = {**base_overrides, "grtresna_shell_lumps": grtresna_lumps}
    if nonspherical and not use_grtresna:
        base_overrides = {**ANGULAR_BASE_OVERRIDES, **base_overrides}

    grtresna_config = None
    solved_ftl_gate_config = None
    if use_grtresna:
        from .grtresna.domain import GRTresnaDomainConfig
        from .grtresna.solver import GRTresnaConfig

        grtresna_domain = GRTresnaDomainConfig(
            full_z=bool(getattr(args, "grtresna_full_z", False)),
            l_full=getattr(args, "grtresna_evolution_l_full", 64.0),
            n_full=getattr(args, "grtresna_evolution_n_full", 64),
            grtresna_l=getattr(args, "grtresna_domain_l", 128.0),
            grtresna_nx=getattr(args, "grtresna_domain_nx", 64),
            grtresna_ny=getattr(args, "grtresna_domain_ny", 64),
            grtresna_nz=getattr(args, "grtresna_domain_nz", None),
            gridinit_nx=getattr(args, "grtresna_gridinit_nx", 64),
            gridinit_ny=getattr(args, "grtresna_gridinit_ny", 64),
            gridinit_nz=getattr(args, "grtresna_gridinit_nz", 64),
        )
        base_overrides = {**base_overrides, **grtresna_domain.evolution_overrides()}
        grtresna_config = GRTresnaConfig(
            mpi_ranks=getattr(args, "grtresna_ranks", 8),
            max_NL_iterations=getattr(args, "grtresna_iterations", 50),
            timeout=getattr(args, "grtresna_timeout", 3600),
            max_level=getattr(args, "grtresna_max_level", 3),
            refine_threshold=getattr(args, "grtresna_refine_threshold", 0.5),
            regrid_radius=getattr(args, "grtresna_regrid_radius", 0.0),
            coefficient_average_type=getattr(args, "grtresna_coefficient_average_type", "harmonic"),
            psi_relaxation=getattr(args, "grtresna_psi_relaxation", 1.0),
            psi_floor=getattr(args, "grtresna_psi_floor", -1.0),
            maximal_jacobian_cap=getattr(args, "grtresna_jacobian_cap", -1.0),
            bh1_bare_mass=0.0,
            bh1_spin=(0.0, 0.0, 0.0),
            bh2_bare_mass=0.0,
            dphi=0.0,
            dpi=0.0,
            cleanup=True,
        )
        grtresna_config = grtresna_domain.apply_to_solver(grtresna_config)
        solved_ftl_gate_config = SolvedFtlGateConfig(
            f_op_floor=getattr(args, "solved_ftl_f_op_floor", 1.0e-4),
            near_luminal_speed_floor=getattr(args, "solved_ftl_near_luminal_speed_floor", 0.99),
            superluminal_speed_floor=getattr(args, "solved_ftl_superluminal_speed_floor", 1.01),
            superluminal_fraction_floor=getattr(args, "solved_ftl_superluminal_fraction_floor", 0.02),
            max_physical_coord_speed=getattr(args, "solved_ftl_max_physical_coord_speed", 8.0),
            max_physical_f_op=getattr(args, "solved_ftl_max_physical_f_op", 0.85),
            rejection_speed_target=getattr(args, "solved_ftl_rejection_speed_target", 1.01),
        )
        grtresna_convergence_config = GRTresnaConvergenceConfig(
            max_ham_pct=getattr(args, "grtresna_max_ham_pct", 5.0),
            max_mom_pct=getattr(args, "grtresna_max_mom_pct", 5.0),
        )
    else:
        grtresna_convergence_config = None

    archive = run_qd_search(
        runs_dir=Path(args.runs_dir).expanduser().resolve(),
        executable=executable,
        iterations=args.iterations,
        batch_size=args.batch_size,
        bins=args.bins,
        init_random=args.init_random,
        seed=args.seed,
        base_overrides=base_overrides,
        search_space=search_space,
        template=template,
        example=example,
        name=args.name,
        dry_run=args.dry_run,
        constrained=not use_grtresna,
        phantom=getattr(args, "phantom", False) or getattr(args, "phantom_default", False),
        use_preflight=(False if use_grtresna else getattr(args, "preflight", False)),
        cuda_devices=args.cuda_devices,
        gpu_ids=getattr(args, "gpu_ids", None),
        check_params=not args.skip_check_params,
        score_weights=getattr(args, "score_weights", None),
        objective_mode=getattr(args, "objective_mode", "weighted"),
        ftl_L=getattr(args, "ftl_L", None),
        consume_plotfiles=getattr(args, "consume_plotfiles", True),
        consumer_radii=getattr(args, "consumer_radii", [4.0, 8.0]),
        consumer_keep_last=getattr(args, "consumer_keep_last", 1),
        descriptor_mode=getattr(args, "descriptor_mode", "legacy"),
        grtresna=use_grtresna,
        grtresna_config=grtresna_config,
        grtresna_solved_ftl_gate=use_grtresna,
        solved_ftl_gate_config=solved_ftl_gate_config,
        grtresna_convergence_config=grtresna_convergence_config,
        grtresna_postload_gate=_grtresna_postload_gate_enabled(args, use_grtresna=use_grtresna),
        postload_gate_config=(
            _postload_gate_config_from_args(args) if use_grtresna else None
        ),
    )
    best = archive.best
    print(json.dumps({
        "num_elites": len(archive.cells),
        "coverage": archive.coverage,
        "best_score": best.score if best else None,
        "best_episode": best.episode if best else None,
    }, indent=2))
    return 0


def _run_pareto_command(args: argparse.Namespace) -> int:
    from .search.pareto import front_to_dict, load_trajectory_points, pareto_front

    points = load_trajectory_points(Path(args.trajectory).expanduser().resolve())
    front = pareto_front(points)
    payload = front_to_dict(front)
    if args.output:
        write_json(Path(args.output).expanduser().resolve(), payload)
    print(json.dumps(payload, indent=2))
    return 0


def _run_warpfactory_command(args: argparse.Namespace) -> int:
    from .metrics import warpfactory as wf

    if getattr(args, "convergence", False):
        result = wf.convergence_order(
            velocity=args.velocity, half_width=args.half_width, dt=args.dt
        )
        print(json.dumps(result, indent=2))
        return 0

    if args.metric == "minkowski":
        g, spacing = wf.minkowski_metric(
            half_width=args.half_width, n_space=args.n_space, dt=args.dt
        )
    elif args.metric == "alcubierre":
        g, spacing = wf.alcubierre_metric(
            velocity=args.velocity,
            bubble_radius=args.bubble_radius,
            sigma=args.sigma,
            half_width=args.half_width,
            n_space=args.n_space,
            dt=args.dt,
        )
    else:  # pragma: no cover - argparse choices guard this
        raise ValueError(f"unknown metric {args.metric!r}")

    report = wf.evaluate_four_metric(
        g,
        spacing,
        n_directions=args.n_directions,
        n_speeds=args.n_speeds,
        max_speed=args.max_speed,
    )
    payload = dataclass_to_dict(report)
    if args.output:
        write_json(Path(args.output).expanduser().resolve(), payload)
    print(json.dumps(payload, indent=2))
    return 0


def _run_atlas_command(args: argparse.Namespace, base_overrides: dict[str, Any]) -> int:
    example = resolve_example(args.example)
    executable = None
    if not args.dry_run:
        executable = resolve_executable(
            args.executable,
            example=example,
            mpi_ranks=args.mpi_ranks,
            comp=args.comp,
            cuda=not args.no_cuda,
            debug=args.debug,
        )

    template = Path(args.template).expanduser().resolve() if args.template else example.template
    paths, records, summary = run_atlas(
        runs_dir=Path(args.runs_dir).expanduser().resolve(),
        executable=executable,
        count=args.count,
        seed=args.seed,
        base_overrides=base_overrides,
        template=template,
        example=example,
        name=args.name,
        dry_run=args.dry_run,
        stop_on_failure=args.stop_on_failure,
        cuda_devices=args.cuda_devices,
        check_params=not args.skip_check_params,
        constrained=getattr(args, "constrained", False),
        phantom=getattr(args, "phantom", False),
        preflight=getattr(args, "preflight", False),
    )
    print(json.dumps({"atlas": str(paths.root), "records": len(records), "summary": summary}, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run isolated GRTeclyn example episodes.")
    parser.add_argument("--runs-dir", default=str(default_runs_dir()), help="Directory for episode folders.")
    parser.add_argument(
        "--example",
        default="SupportedWormholeCollapse",
        choices=["SupportedWormholeCollapse", "RadialRecipe", "RotatingWormholeCollapse"],
        help="GRTeclyn example to run.",
    )
    parser.add_argument(
        "--template",
        default=None,
        help="Source params template. Defaults to the selected example template.",
    )
    parser.add_argument("--executable", default=None, help="Executable path. Defaults to the selected example binary name.")
    parser.add_argument("--mpi-ranks", type=int, default=1, help="MPI ranks; >1 selects the MPI executable name.")
    parser.add_argument("--comp", default="gnu", help="Compiler tag in the executable name.")
    parser.add_argument("--debug", action="store_true", help="Select DEBUG executable naming.")
    parser.add_argument("--no-cuda", action="store_true", help="Select non-CUDA executable naming.")
    parser.add_argument("--cuda-devices", default=None, help="CUDA_VISIBLE_DEVICES value for the run.")
    parser.add_argument("--skip-check-params", action="store_true", help="Skip the check_params=1 preflight.")
    parser.add_argument("--consume-plotfiles", action="store_true", help="Run consume_plotfiles as a side process.")
    parser.add_argument("--consumer-radii", nargs="+", type=float, default=[8.0, 16.0], help="Radii for consume_plotfiles.")
    parser.add_argument("--consumer-delete", action="store_true", help="Let consume_plotfiles delete old plotfiles.")
    parser.add_argument(
        "--consumer-keep-last",
        type=int,
        default=1,
        help="Plotfiles to retain when consuming (>=3 lets evolved-FTL/effective-EC score the evolved geometry).",
    )
    parser.add_argument("--dry-run", action="store_true", help="Create episode files without launching GRTeclyn.")
    parser.add_argument("--constrained", action="store_true", help="Derive phi from chi to satisfy Hamiltonian constraint (RadialRecipe only).")
    parser.add_argument("--phantom", action="store_true", help="Use phantom scalar field coupling (negative kinetic term) in constrained mode.")
    parser.add_argument("--no-phantom", dest="no_phantom", action="store_true", help="Force the constrained solve to use NORMAL matter (rho>=0) during optimize: search for FTL without exotic matter.")
    parser.add_argument("--preflight", action="store_true", help="Pre-flight constraint check; reject bad candidates before GPU launch (RadialRecipe only).")
    parser.add_argument("--seed-name", default=None, choices=list_seeds(), help="Start from a known-solution seed (RadialRecipe only).")
    parser.add_argument(
        "--candidate-id",
        default=None,
        help="Validation candidate ID from validate_guesser (e.g. random_000, bubble_wall_016).",
    )
    parser.add_argument(
        "--nonspherical-id",
        default=None,
        help="Non-spherical ray-validation candidate (e.g. quadrupole_bubble_001).",
    )
    parser.add_argument(
        "--validation-seed",
        type=int,
        default=42,
        help="RNG seed for deterministic candidate lookup.",
    )
    parser.add_argument("--set", action="append", type=_parse_override, default=[], metavar="KEY=VALUE", help="Params override.")
    parser.add_argument(
        "--score-weight",
        action="append",
        type=_parse_override,
        default=[],
        metavar="KEY=VALUE",
        help="Override score component weight (e.g. ftl_shortcut=5.0).",
    )
    parser.add_argument(
        "--ftl-L",
        type=float,
        default=None,
        dest="ftl_L",
        help="Half-axis length for t=0 FTL profile integration.",
    )
    parser.add_argument("--name", default=None, help="Episode folder name.")

    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("reproduce", help="Run one episode from the template plus overrides.")

    sweep = subparsers.add_parser("sweep", help="Run a random sweep over current wormhole parameters.")
    sweep.add_argument("--count", type=int, default=1, help="Number of random episodes.")
    sweep.add_argument("--seed", type=int, default=None, help="Random seed.")
    sweep.add_argument("--stop-on-failure", action="store_true", help="Stop sweep at first non-zero run.")

    atlas = subparsers.add_parser("atlas", help="Run a low-resolution failure-atlas batch.")
    atlas.add_argument("--count", type=int, default=10, help="Number of sampled episodes.")
    atlas.add_argument("--seed", type=int, default=None, help="Random seed.")
    atlas.add_argument("--stop-on-failure", action="store_true", help="Stop atlas at first solver failure.")

    opt = subparsers.add_parser("optimize", help="CMA-ES optimization over RadialRecipe coefficients.")
    opt.add_argument("--max-generations", type=int, default=50, help="Maximum CMA-ES generations.")
    opt.add_argument("--population-size", type=int, default=None, help="CMA-ES population size (default: auto, or num GPUs).")
    opt.add_argument("--sigma0", type=float, default=0.3, help="Initial CMA-ES step size.")
    opt.add_argument("--seed", type=int, default=None, help="Random seed for CMA-ES.")
    opt.add_argument("--gpu-ids", nargs="+", type=int, default=None, help="GPU indices for parallel eval (e.g. 0 1 2 3 4 5 6 7).")
    opt.add_argument(
        "--objective-mode",
        choices=["weighted", "ftl_first"],
        default="weighted",
        help="Scoring scalarization: weighted legacy score or FTL-first ordering.",
    )
    opt.add_argument(
        "--warm-start-trajectory",
        action="append",
        default=[],
        help="Previous trajectory.jsonl to seed the first generation; repeatable.",
    )
    opt.add_argument(
        "--warm-start-top-k",
        type=int,
        default=8,
        help="Top prior candidates loaded from warm-start trajectories.",
    )
    opt.add_argument(
        "--warm-start-jitter",
        type=float,
        default=0.08,
        help="Fraction of each parameter range used when jittering warm-starts.",
    )
    opt.add_argument(
        "--random-injection-fraction",
        type=float,
        default=0.0,
        help="Fraction of every generation replaced by random bounded candidates.",
    )
    opt.add_argument(
        "--exotic-injection-fraction",
        type=float,
        default=0.0,
        help="Fraction of every generation replaced by forced exotic templates.",
    )
    opt.add_argument("--surrogate", action="store_true", help="Enable RBF surrogate pre-screening to skip low-value GPU evaluations.")
    opt.add_argument("--surrogate-keep-fraction", type=float, default=0.5, help="Fraction of each generation evaluated on GPU when surrogate is active.")
    opt.add_argument(
        "--no-consume-plotfiles",
        action="store_true",
        help="Disable the plotfile consumer. By default optimize streams each "
             "plotfile into radial/frame products and DELETES the raw plotfile, "
             "so disk stays bounded and per-eval frames are produced. Pass this "
             "to keep raw plotfiles instead (uses a lot of disk).",
    )
    opt.add_argument(
        "--nonspherical",
        action="store_true",
        help="Open up non-spherical geometries: add axisymmetric Legendre angular "
             "modes on the lapse and shift (gauge sector, constraint-preserving, "
             "no extra exotic matter) so the search can sculpt directional FTL channels.",
    )
    opt.add_argument(
        "--grtresna",
        action="store_true",
        help="GRTresna-in-the-loop: produce constraint-satisfying, momentum-carrying "
             "scalar-field initial data via the 3D elliptic solver each evaluation "
             "(searches grtresna_lump_* dims), then evolve it on GPU. Replaces the 1D "
             "radial recipe search space.",
    )
    opt.add_argument(
        "--grtresna-ranks", type=int, default=8,
        help="MPI ranks for each GRTresna solve (default: 8).",
    )
    opt.add_argument(
        "--grtresna-lumps", type=int, default=5,
        help="Number of scalar lumps in the GRTresna matter basis (default: 5). "
             "Each lump adds 11 searched dimensions (amp/width/center/velocity/"
             "omega/mode/exotic).",
    )
    opt.add_argument(
        "--grtresna-ansatz",
        choices=["free", "ring", "shell"],
        default="free",
        help="GRTresna matter parameterization. 'free' searches every lump "
             "independently (11*K dimensions). 'ring' searches a reduced rotating "
             "counterflow/exotic ring template (14D, planar/equatorial) and "
             "expands it into K lumps. 'shell' is the full-sphere discovery "
             "ansatz (16D): lumps cover the whole 2-sphere with an arbitrary "
             "orientation axis and poloidal+toroidal currents, reaching 3D "
             "configurations the planar ring cannot.",
    )
    opt.add_argument(
        "--grtresna-shell-profile",
        choices=["middle", "compact", "outer_precursor", "inner_shift"],
        default="compact",
        help="Shell ansatz bounds preset. 'compact' caps lump width (default); "
             "'middle' restores the pre-170329Z bounds; 'outer_precursor' and "
             "'inner_shift' bias toward the eval-128 / eval-57 leader regimes.",
    )
    opt.add_argument(
        "--grtresna-full-z",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Use a full-z GRTresna/GRTeclyn domain instead of the legacy reflective z=0 half-domain.",
    )
    opt.add_argument(
        "--grtresna-evolution-l-full",
        type=float,
        default=64.0,
        help="GRTeclyn evolution box length used for grtresna domain mapping.",
    )
    opt.add_argument(
        "--grtresna-evolution-n-full",
        type=int,
        default=64,
        help="GRTeclyn evolution grid size override for grtresna search episodes.",
    )
    opt.add_argument(
        "--grtresna-domain-l",
        type=float,
        default=128.0,
        help="Physical side length of the GRTresna solve domain.",
    )
    opt.add_argument("--grtresna-domain-nx", type=int, default=64)
    opt.add_argument("--grtresna-domain-ny", type=int, default=64)
    opt.add_argument(
        "--grtresna-domain-nz",
        type=int,
        default=None,
        help="GRTresna solve z cells; default is nx for full-z and nx/2 for half-z.",
    )
    opt.add_argument("--grtresna-gridinit-nx", type=int, default=64)
    opt.add_argument("--grtresna-gridinit-ny", type=int, default=64)
    opt.add_argument("--grtresna-gridinit-nz", type=int, default=64)
    opt.add_argument(
        "--grtresna-iterations", type=int, default=50,
        help="Max non-linear iterations per GRTresna solve (default: 50).",
    )
    opt.add_argument(
        "--grtresna-timeout", type=int, default=3600,
        help="Wall-clock timeout in seconds for each GRTresna solve.",
    )
    opt.add_argument(
        "--grtresna-max-level", type=int, default=3,
        help="Max AMR level for each GRTresna solve (default: 3).",
    )
    opt.add_argument(
        "--grtresna-refine-threshold", type=float, default=0.5,
        help="Refinement threshold for GRTresna AMR tagging (default: 0.5).",
    )
    opt.add_argument(
        "--grtresna-regrid-radius", type=float, default=0.0,
        help="Forced GRTresna regrid radius; 0 uses source-based tagging (default: 0).",
    )
    opt.add_argument(
        "--grtresna-coefficient-average-type",
        choices=["harmonic", "arithmetic"],
        default="harmonic",
        help="Coefficient averaging for GRTresna AMR multigrid; exotic candidates override harmonic to arithmetic.",
    )
    opt.add_argument(
        "--grtresna-psi-relaxation", type=float, default=1.0,
        help="Under-relaxation for GRTresna nonlinear psi updates; exotic candidates override 1.0 to a safer value.",
    )
    opt.add_argument(
        "--grtresna-psi-floor", type=float, default=-1.0,
        help="Positive psi floor for GRTresna nonlinear updates; exotic candidates override non-positive values.",
    )
    opt.add_argument(
        "--grtresna-jacobian-cap", type=float, default=-1.0,
        help="Optional absolute cap for the maximal-slicing psi Jacobian; exotic candidates set a safe default.",
    )
    opt.add_argument(
        "--grtresna-max-ham-pct", type=float, default=5.0,
        help="Reject GRTresna solves above this Hamiltonian residual (%%).",
    )
    opt.add_argument(
        "--grtresna-max-mom-pct", type=float, default=5.0,
        help="Reject GRTresna solves above this momentum residual (%%).",
    )
    opt.add_argument(
        "--solved-ftl-f-op-floor", type=float, default=1.0e-4,
        help="Solved-geometry F_op threshold that passes a candidate to GPU.",
    )
    opt.add_argument(
        "--solved-ftl-near-luminal-speed-floor", type=float, default=0.99,
        help="Exploratory solved-geometry max_c floor for subluminal near misses (requires max_c < 1).",
    )
    opt.add_argument(
        "--solved-ftl-superluminal-speed-floor", type=float, default=1.01,
        help="Solved-geometry max_c floor for genuine superluminal local-speed precursors.",
    )
    opt.add_argument(
        "--solved-ftl-superluminal-fraction-floor", type=float, default=0.02,
        help="Minimum solved-geometry superluminal area fraction for speed-precursor passes.",
    )
    opt.add_argument(
        "--solved-ftl-max-physical-coord-speed", type=float, default=8.0,
        help="Reject solved slices above this max_c as near-degenerate numerical artifacts.",
    )
    opt.add_argument(
        "--solved-ftl-max-physical-f-op", type=float, default=0.85,
        help="Reject solved slices above this F_op as near-degenerate numerical artifacts.",
    )
    opt.add_argument(
        "--solved-ftl-rejection-speed-target", type=float, default=1.01,
        help="max_c target used to grade solved-FTL rejection fitness.",
    )
    opt.add_argument(
        "--grtresna-postload-gate",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Run a short GRTeclyn post-load Ham/Mom check before full evolution.",
    )
    opt.add_argument("--postload-max-ham-l2", type=float, default=1.0e-2)
    opt.add_argument("--postload-max-mom-l2", type=float, default=1.0e-2)

    qd = subparsers.add_parser("qd", help="MAP-Elites quality-diversity search (Spacetime Failure Atlas).")
    qd.add_argument("--iterations", type=int, default=10, help="Number of MAP-Elites batches after the initial fill.")
    qd.add_argument("--batch-size", type=int, default=None, help="Candidates per batch (default: number of GPUs or 4).")
    qd.add_argument("--bins", type=int, default=8, help="Behavior-descriptor grid resolution per axis.")
    qd.add_argument("--init-random", type=int, default=None, help="Random candidates in the initial fill (default: batch size).")
    qd.add_argument("--seed", type=int, default=None, help="Random seed.")
    qd.add_argument("--gpu-ids", nargs="+", type=int, default=None, help="GPU indices for parallel eval.")
    qd.add_argument(
        "--nonspherical",
        action="store_true",
        help="Use the geometry-first nonspherical lapse/shift angular search space.",
    )
    qd.add_argument(
        "--descriptor-mode",
        choices=["legacy", "channel"],
        default="legacy",
        help="MAP-Elites descriptors: legacy FTL/mechanism grid, or channel path-closeness/mechanism-balance grid.",
    )
    qd.add_argument(
        "--objective-mode",
        choices=["weighted", "ftl_first"],
        default="weighted",
        help="Scoring scalarization used as elite quality.",
    )
    qd.add_argument(
        "--grtresna",
        action="store_true",
        help="Evaluate MAP-Elites candidates through the GRTresna constraint solve before GPU evolution.",
    )
    qd.add_argument("--grtresna-ranks", type=int, default=8)
    qd.add_argument("--grtresna-lumps", type=int, default=5)
    qd.add_argument(
        "--grtresna-ansatz",
        choices=["free", "ring", "shell"],
        default="free",
        help="GRTresna matter parameterization for QD candidates.",
    )
    qd.add_argument(
        "--grtresna-shell-profile",
        choices=["middle", "compact", "outer_precursor", "inner_shift"],
        default="compact",
        help="Shell ansatz bounds preset when --grtresna-ansatz=shell.",
    )
    qd.add_argument("--grtresna-full-z", action=argparse.BooleanOptionalAction, default=False)
    qd.add_argument("--grtresna-evolution-l-full", type=float, default=64.0)
    qd.add_argument("--grtresna-evolution-n-full", type=int, default=64)
    qd.add_argument("--grtresna-domain-l", type=float, default=128.0)
    qd.add_argument("--grtresna-domain-nx", type=int, default=64)
    qd.add_argument("--grtresna-domain-ny", type=int, default=64)
    qd.add_argument("--grtresna-domain-nz", type=int, default=None)
    qd.add_argument("--grtresna-gridinit-nx", type=int, default=64)
    qd.add_argument("--grtresna-gridinit-ny", type=int, default=64)
    qd.add_argument("--grtresna-gridinit-nz", type=int, default=64)
    qd.add_argument("--grtresna-iterations", type=int, default=50)
    qd.add_argument("--grtresna-timeout", type=int, default=3600, help="Wall-clock timeout in seconds for each GRTresna solve.")
    qd.add_argument("--grtresna-max-level", type=int, default=3)
    qd.add_argument("--grtresna-refine-threshold", type=float, default=0.5)
    qd.add_argument("--grtresna-regrid-radius", type=float, default=0.0)
    qd.add_argument(
        "--grtresna-coefficient-average-type",
        choices=["harmonic", "arithmetic"],
        default="harmonic",
    )
    qd.add_argument("--grtresna-psi-relaxation", type=float, default=1.0)
    qd.add_argument("--grtresna-psi-floor", type=float, default=-1.0)
    qd.add_argument("--grtresna-jacobian-cap", type=float, default=-1.0)
    qd.add_argument(
        "--grtresna-max-ham-pct", type=float, default=5.0,
        help="Reject GRTresna solves above this Hamiltonian residual (%%).",
    )
    qd.add_argument(
        "--grtresna-max-mom-pct", type=float, default=5.0,
        help="Reject GRTresna solves above this momentum residual (%%).",
    )
    qd.add_argument("--solved-ftl-f-op-floor", type=float, default=1.0e-4)
    qd.add_argument("--solved-ftl-near-luminal-speed-floor", type=float, default=0.99)
    qd.add_argument("--solved-ftl-superluminal-speed-floor", type=float, default=1.01)
    qd.add_argument("--solved-ftl-superluminal-fraction-floor", type=float, default=0.02)
    qd.add_argument("--solved-ftl-max-physical-coord-speed", type=float, default=8.0)
    qd.add_argument("--solved-ftl-max-physical-f-op", type=float, default=0.85)
    qd.add_argument("--solved-ftl-rejection-speed-target", type=float, default=1.01)
    qd.add_argument(
        "--grtresna-postload-gate",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Run a short GRTeclyn post-load Ham/Mom check before full evolution.",
    )
    qd.add_argument("--postload-max-ham-l2", type=float, default=1.0e-2)
    qd.add_argument("--postload-max-mom-l2", type=float, default=1.0e-2)

    pareto = subparsers.add_parser("pareto", help="Extract the multi-objective Pareto front from a trajectory.jsonl.")
    pareto.add_argument("--trajectory", required=True, help="Path to an optimizer trajectory.jsonl.")
    pareto.add_argument("--output", default=None, help="Optional path to write the front JSON.")

    wfp = subparsers.add_parser(
        "warpfactory",
        help="Warp Factory-style multi-observer energy-condition evaluation of an analytic 4-metric.",
    )
    wfp.add_argument("--metric", choices=["minkowski", "alcubierre"], default="alcubierre", help="Analytic metric to evaluate.")
    wfp.add_argument("--velocity", type=float, default=0.5, help="Bubble velocity (Alcubierre).")
    wfp.add_argument("--bubble-radius", type=float, default=2.0, help="Bubble radius (Alcubierre).")
    wfp.add_argument("--sigma", type=float, default=2.0, help="Wall steepness (Alcubierre).")
    wfp.add_argument("--half-width", type=float, default=4.0, help="Spatial half-extent of the grid.")
    wfp.add_argument("--n-space", type=int, default=22, help="Grid points per spatial axis.")
    wfp.add_argument("--dt", type=float, default=0.2, help="Time spacing for d_t finite differences.")
    wfp.add_argument("--n-directions", type=int, default=60, help="Sampled observer directions on the sphere.")
    wfp.add_argument("--n-speeds", type=int, default=4, help="Sampled timelike observer speeds.")
    wfp.add_argument("--max-speed", type=float, default=0.9, help="Maximum timelike observer speed (<1).")
    wfp.add_argument("--convergence", action="store_true", help="Run the finite-difference convergence-order self-test instead of an EC report.")
    wfp.add_argument("--output", default=None, help="Optional path to write the report JSON.")

    validate = subparsers.add_parser(
        "validate",
        help="Batch-validate the metric guesser on synthetic candidates (no optimizer).",
    )
    validate.add_argument("--seed", type=int, default=42, help="Random seed for candidate generation.")
    validate.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_out"),
        help="Directory for guesser_validation.csv and summary JSON.",
    )
    validate.add_argument("--no-write", action="store_true", help="Print summary only; skip file output.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.phantom_default = False
    args.initial_data_source = None
    args.score_weights = _parse_score_weights(args.score_weight) or None
    overrides = _collect_overrides(args.set)
    source = None
    if (
        getattr(args, "seed_name", None)
        or getattr(args, "candidate_id", None)
        or getattr(args, "nonspherical_id", None)
    ):
        base_overrides, phantom_default, source = resolve_initial_data_overrides(
            seed_name=getattr(args, "seed_name", None),
            candidate_id=getattr(args, "candidate_id", None),
            nonspherical_id=getattr(args, "nonspherical_id", None),
            validation_seed=getattr(args, "validation_seed", 42),
        )
        overrides = {**base_overrides, **overrides}
        args.phantom_default = phantom_default
        args.initial_data_source = source
    if args.command == "optimize":
        return _run_optimize_command(args, overrides)
    if args.command == "qd":
        return _run_qd_command(args, overrides)
    if args.command == "pareto":
        return _run_pareto_command(args)
    if args.command == "warpfactory":
        return _run_warpfactory_command(args)
    if args.command == "validate":
        output_dir = None if args.no_write else args.output_dir
        run_validation(seed=args.seed, output_dir=output_dir)
        return 0
    if args.command == "atlas":
        return _run_atlas_command(args, overrides)
    if args.command == "sweep":
        return _run_sweep(args, overrides)
    return _run_single(args, overrides)


if __name__ == "__main__":
    raise SystemExit(main())
