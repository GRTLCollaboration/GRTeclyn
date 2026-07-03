"""CLI argument parser for grteclyn-wrapper."""

from __future__ import annotations

import argparse
from pathlib import Path

from ..core.config import default_runs_dir
from ..initial_data.seeds import list_seeds
from .args import parse_override
from .grtresna_args import (
    add_grtresna_matter_selection_args,
    add_grtresna_solved_ftl_gate_arg,
    add_grtresna_speed_args,
)
from .pipeline_args import add_pipeline_args


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
    parser.add_argument("--set", action="append", type=parse_override, default=[], metavar="KEY=VALUE", help="Params override.")
    parser.add_argument(
        "--score-weight",
        action="append",
        type=parse_override,
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
    opt.add_argument(
        "--target-evals",
        type=int,
        default=None,
        help="Stop after this many GPU/CPU evaluations (overrides max_generations when reached first).",
    )
    opt.add_argument("--population-size", type=int, default=None, help="CMA-ES population size (default: auto, or num GPUs).")
    opt.add_argument("--sigma0", type=float, default=0.3, help="Initial CMA-ES step size.")
    opt.add_argument("--seed", type=int, default=None, help="Random seed for CMA-ES.")
    opt.add_argument("--gpu-ids", nargs="+", type=int, default=None, help="GPU indices for parallel eval (e.g. 0 1 2 3 4 5 6 7).")
    opt.add_argument(
        "--objective-mode",
        choices=["weighted", "ftl_first", "robust_ftl", "general_ftl", "critical_collapse"],
        default="weighted",
        help="Scoring scalarization: weighted legacy score, FTL-first ordering, "
        "robust_ftl (FTL-first tilted toward persistent/healthy/low-exotic), "
        "or general_ftl (gauge-invariant shortcut only; no warp-motor shaping).",
    )
    opt.add_argument(
        "--pin-dimension",
        action="append",
        default=[],
        dest="pin_dimension",
        metavar="KEY=VALUE",
        help="Pin a search dimension to a constant (repeatable). The key is "
        "removed from the optimizer search space and forced via base "
        "overrides, e.g. --pin-dimension grtresna_matter_layout=2 "
        "--pin-dimension grtresna_shell_toroidal_velocity=0.",
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
        "--warm-start-include-near-miss",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Also warm-start from top graded pre-GPU rejections in trajectory files (default: on with --grtresna).",
    )
    opt.add_argument(
        "--warm-start-near-miss-k",
        type=int,
        default=4,
        help="Max near-miss vectors to add when warm-starting CMA-ES.",
    )
    opt.add_argument(
        "--keep-top-eval-dirs",
        type=int,
        default=0,
        help="After each generation, delete completed eval_* dirs outside the top N "
        "scored records. trajectory.jsonl stays intact; 0 disables.",
    )
    opt.add_argument(
        "--ftl-retention",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Also retain eval dirs that hold the run-best FTL peak per metric "
        "(f_geo, speed, lifetime, ...); writes ftl_retention.jsonl + ftl_champions.json.",
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
        choices=["free", "ring", "shell", "sh", "trajectory", "boson_star", "splash"],
        default="free",
        help="GRTresna matter parameterization. 'free' searches every lump "
             "independently (11*K dimensions). 'ring' searches a reduced rotating "
             "counterflow/exotic ring template (14D, planar/equatorial) and "
             "expands it into K lumps. 'shell' is the full-sphere discovery "
             "ansatz (16D): lumps cover the whole 2-sphere with an arbitrary "
             "orientation axis and poloidal+toroidal currents, reaching 3D "
             "configurations the planar ring cannot. "
             "'trajectory' enables trajectory-guided FTL geometry survey. "
             "Legacy alias: 'boson_star' sets --grtresna-matter-sector boson_star.",
    )
    add_grtresna_matter_selection_args(opt)
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
    add_grtresna_speed_args(opt)
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
        "--grtresna-keep-source", action="store_true", default=False,
        help="Keep the GRTresna Chombo HDF5 + workdir per eval (disables cleanup). "
             "Use for conversion validation/debugging; consumes much more disk.",
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
    add_grtresna_solved_ftl_gate_arg(opt)
    opt.add_argument(
        "--grtresna-postload-gate",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Run a short GRTeclyn post-load Ham/Mom check before full evolution.",
    )
    opt.add_argument("--postload-max-ham-l2", type=float, default=1.0e-2)
    opt.add_argument("--postload-max-mom-l2", type=float, default=1.0e-2)
    add_pipeline_args(opt)

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
        "--seed-eval-dirs",
        nargs="+",
        default=None,
        metavar="EVAL_DIR",
        help="Warm-start the initial MAP-Elites population from these prior eval "
        "directories (reads each metadata.json 'overrides'). Use to seed a fresh "
        "campaign from promoted survivors.",
    )
    qd.add_argument(
        "--keep-top-eval-dirs",
        type=int,
        default=0,
        help="After each batch, delete completed eval_* directories outside the "
        "top N scored records. trajectory.jsonl/archive.json stay intact; 0 disables.",
    )
    qd.add_argument(
        "--ftl-retention",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Also retain eval dirs that hold the campaign-best FTL peak per metric "
        "(f_geo, speed, lifetime, etc.); writes ftl_retention.jsonl + ftl_champions.json.",
    )
    qd.add_argument(
        "--descriptor-mode",
        choices=["legacy", "channel", "speed_horizon", "speed_super", "ftl_lifetime", "wave_focusing", "gw_beam"],
        default="legacy",
        help="MAP-Elites descriptors: legacy FTL/mechanism grid, channel "
        "path-closeness/mechanism-balance grid (needs shift>0), speed_horizon "
        "cone-tilt(max_local_speed) vs horizon-free(min_theta_plus) grid that "
        "illuminates the fast-but-not-trapped niche without needing shift, "
        "speed_super recalibrated cone-tilt vs superluminal_fraction grid "
        "(localized vs widespread superluminal region) that stays discriminating "
        "in the nontrivial-but-not-operational regime, ftl_lifetime "
        "peak gauge-invariant strength vs FTL-lifetime fraction (time-resolved) "
        "grid that separates transient shortcuts from sustained ones, or "
        "gw_beam (log total GW power vs Z-axis beam ratio).",
    )
    qd.add_argument(
        "--objective-mode",
        choices=["weighted", "ftl_first", "robust_ftl", "general_ftl", "critical_collapse", "gw_beam"],
        default="weighted",
        help="Scoring scalarization used as elite quality. robust_ftl tilts "
        "ftl_first toward persistent/healthy/low-exotic geometries; "
        "general_ftl rewards gauge-invariant shortcuts only; gw_beam "
        "rewards strong directional gravitational-wave emission.",
    )
    qd.add_argument(
        "--pin-dimension",
        action="append",
        default=[],
        dest="pin_dimension",
        metavar="KEY=VALUE",
        help="Pin a search dimension to a constant (repeatable). The key is "
        "removed from the optimizer search space and forced via base "
        "overrides, e.g. --pin-dimension grtresna_matter_layout=2 "
        "--pin-dimension grtresna_shell_toroidal_velocity=0.",
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
        choices=["free", "ring", "shell", "sh", "trajectory", "boson_star", "splash"],
        default="free",
        help="GRTresna geometry ansatz for scalar lumps (shell/ring/free/trajectory). "
             "'trajectory' enables the trajectory-guided FTL geometry survey. "
             "Legacy: 'boson_star' alias sets --grtresna-matter-sector boson_star.",
    )
    add_grtresna_matter_selection_args(qd)
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
    add_grtresna_speed_args(qd)
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
        "--grtresna-keep-source", action="store_true", default=False,
        help="Keep the GRTresna Chombo HDF5 + workdir per eval (disables cleanup). "
             "Use for conversion validation/debugging; consumes much more disk.",
    )
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
    add_grtresna_solved_ftl_gate_arg(qd)
    qd.add_argument(
        "--grtresna-postload-gate",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Run a short GRTeclyn post-load Ham/Mom check before full evolution.",
    )
    qd.add_argument("--postload-max-ham-l2", type=float, default=1.0e-2)
    qd.add_argument("--postload-max-mom-l2", type=float, default=1.0e-2)
    qd.add_argument(
        "--resume",
        action="store_true",
        help="Continue an existing MAP-Elites campaign in --name (loads archive + trajectory).",
    )
    qd.add_argument(
        "--target-evals",
        type=int,
        default=None,
        help="Total evaluations to reach (sets batch count from remaining evals on resume).",
    )
    qd.add_argument(
        "--pre-gpu-learning",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Learn from graded GRTresna rejections via near-miss pool + shadow archive (default: on with --grtresna).",
    )
    qd.add_argument(
        "--near-miss-pool-size",
        type=int,
        default=32,
        help="Top-K near-miss parents kept for MAP-Elites mutation.",
    )
    add_pipeline_args(qd)

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
