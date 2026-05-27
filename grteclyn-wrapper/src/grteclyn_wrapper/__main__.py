from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Any

from .atlas import run_atlas
from .config import default_runs_dir, resolve_example, resolve_executable
from .constrained_recipe import constrained_overrides
from .episode import create_episode, update_metadata, write_json
from .metrics import dataclass_to_dict, read_episode_metrics
from .optimize import run_optimize
from .params import write_params
from .seeds import get_seed, list_seeds
from .runner import run_episode
from .score import score_episode
from .validate_guesser import run_validation


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


def _finalize_score(episode_dir: Path, target_stop_time: float | None) -> int:
    metrics = read_episode_metrics(episode_dir)
    score = score_episode(metrics, target_stop_time=target_stop_time)
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
    if getattr(args, "constrained", False) and example.name == "RadialRecipe":
        constrained_overrides(overrides, phantom=getattr(args, "phantom", False))
    runs_dir = Path(args.runs_dir).expanduser().resolve()
    episode = create_episode(
        runs_dir,
        name=args.name,
        metadata={"mode": args.command, "example": example.name, "overrides": overrides},
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
    )
    _finalize_score(episode.path, target_stop_time=overrides.get("stop_time"))
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

    x0 = None
    if getattr(args, "seed_name", None):
        seed_obj = get_seed(args.seed_name)
        from .optimize import DEFAULT_SEARCH_SPACE
        x0 = []
        for dim in DEFAULT_SEARCH_SPACE:
            x0.append(float(seed_obj.overrides.get(dim.param_key, dim.center)))

    result = run_optimize(
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
        constrained=getattr(args, "constrained", False),
        phantom=getattr(args, "phantom", False),
        use_preflight=getattr(args, "preflight", False),
        cuda_devices=args.cuda_devices,
        check_params=not args.skip_check_params,
        x0=x0,
    )
    print(json.dumps({
        "best_score": result.best_score,
        "best_episode": result.best_episode,
        "best_params": result.best_params,
        "generations": result.generations,
        "evaluations": result.evaluations,
    }, indent=2))
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
        choices=["SupportedWormholeCollapse", "RadialRecipe"],
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
    parser.add_argument("--dry-run", action="store_true", help="Create episode files without launching GRTeclyn.")
    parser.add_argument("--constrained", action="store_true", help="Derive phi from chi to satisfy Hamiltonian constraint (RadialRecipe only).")
    parser.add_argument("--phantom", action="store_true", help="Use phantom scalar field coupling (negative kinetic term) in constrained mode.")
    parser.add_argument("--preflight", action="store_true", help="Pre-flight constraint check; reject bad candidates before GPU launch (RadialRecipe only).")
    parser.add_argument("--seed-name", default=None, choices=list_seeds(), help="Start from a known-solution seed (RadialRecipe only).")
    parser.add_argument("--set", action="append", type=_parse_override, default=[], metavar="KEY=VALUE", help="Params override.")
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
    opt.add_argument("--population-size", type=int, default=None, help="CMA-ES population size (default: auto).")
    opt.add_argument("--sigma0", type=float, default=0.3, help="Initial CMA-ES step size.")
    opt.add_argument("--seed", type=int, default=None, help="Random seed for CMA-ES.")

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
    overrides = _collect_overrides(args.set)
    if getattr(args, "seed_name", None):
        seed = get_seed(args.seed_name)
        seed_overrides = dict(seed.overrides)
        seed_overrides.update(overrides)
        overrides = seed_overrides
    if args.command == "optimize":
        return _run_optimize_command(args, overrides)
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
