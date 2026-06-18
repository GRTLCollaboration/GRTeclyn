"""GRTresna-related CLI argument helpers."""

from __future__ import annotations

import argparse
import os
from typing import Any

from ..projection.postload_gate import PostLoadGateConfig
from ..grtresna.matter_models import (
    MATTER_COUPLING_CANONICAL,
    MATTER_COUPLING_EXOTIC,
    MATTER_SECTOR_BOSON_STAR,
    MATTER_SECTOR_SCALAR,
    MatterSelection,
    resolve_matter_selection,
)
from ..search.solved_ftl_gate import SolvedFtlGateConfig


def grtresna_postload_gate_enabled(args: argparse.Namespace, *, use_grtresna: bool) -> bool:
    if not use_grtresna:
        return False
    if getattr(args, "grtresna_postload_gate", False):
        return True
    return os.environ.get("POSTLOAD_GATE", "0") == "1"


def add_grtresna_solved_ftl_gate_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--grtresna-solved-ftl-gate",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Reject GRTresna candidates with no solved-geometry coordinate FTL "
            "precursor before GPU evolution. Default: on with --grtresna for "
            "warp-discovery objectives; off for --objective-mode general_ftl "
            "(gauge-invariant null-geodesic scoring)."
        ),
    )


def grtresna_solved_ftl_gate_enabled(
    args: argparse.Namespace,
    *,
    use_grtresna: bool,
    objective_mode: str = "weighted",
) -> bool:
    """Whether to apply the pre-GPU solved operational-FTL filter."""
    if not use_grtresna:
        return False
    explicit = getattr(args, "grtresna_solved_ftl_gate", None)
    if explicit is not None:
        return bool(explicit)
    # general_ftl rewards 4D null-geodesic shortcuts, not coordinate precursors.
    if objective_mode == "general_ftl":
        return False
    return True


def postload_gate_config_from_args(args: argparse.Namespace) -> PostLoadGateConfig:
    return PostLoadGateConfig(
        max_hamiltonian_l2=getattr(args, "postload_max_ham_l2", 1.0e-2),
        max_momentum_l2=getattr(args, "postload_max_mom_l2", 1.0e-2),
    )


def add_grtresna_matter_selection_args(parser: argparse.ArgumentParser) -> None:
    """Matter sector (scalar vs boson star) and energy coupling (canonical vs exotic)."""
    parser.add_argument(
        "--grtresna-matter-sector",
        choices=[MATTER_SECTOR_SCALAR, MATTER_SECTOR_BOSON_STAR],
        default=MATTER_SECTOR_SCALAR,
        help=(
            "Matter sector for GRTresna initial data. "
            "'scalar' = independent real scalar lumps (shell/ring/free ansatz). "
            "'boson_star' = single complex U(1) mini boson star (7-D search)."
        ),
    )
    parser.add_argument(
        "--grtresna-matter-coupling",
        choices=[MATTER_COUPLING_CANONICAL, MATTER_COUPLING_EXOTIC],
        default=MATTER_COUPLING_CANONICAL,
        help=(
            "Energy coupling. 'canonical' = positive-energy matter; "
            "'exotic' = phantom (−rho): all scalar lumps exotic, or boson scalar_sign=-1."
        ),
    )


def resolve_grtresna_matter_from_args(args: argparse.Namespace) -> MatterSelection:
    """Resolve matter selection from CLI, honoring legacy ``--grtresna-ansatz boson_star``."""
    sector = getattr(args, "grtresna_matter_sector", MATTER_SECTOR_SCALAR)
    if getattr(args, "grtresna_ansatz", "free") == MATTER_SECTOR_BOSON_STAR:
        sector = MATTER_SECTOR_BOSON_STAR
    coupling = getattr(args, "grtresna_matter_coupling", MATTER_COUPLING_CANONICAL)
    return resolve_matter_selection(sector, coupling)


def add_grtresna_speed_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--grtresna-nl-exit-tolerance",
        type=float,
        default=1.0,
        help=(
            "Stop GRTresna when both Ham and Mom residuals fall below this %% "
            "(0 disables absolute early exit)."
        ),
    )
    parser.add_argument(
        "--grtresna-nl-stall-tolerance",
        type=float,
        default=0.02,
        help=(
            "Stop when per-iteration Ham/Mom relative improvement drops below "
            "this fraction (0 disables stall early exit)."
        ),
    )
    parser.add_argument(
        "--grtresna-gridinit-workers",
        type=int,
        default=0,
        help=(
            "Threads for parallel Chombo→gridinit box painting "
            "(0 = auto, min(32, cpu_count))."
        ),
    )


def grtresna_speed_kwargs(args: argparse.Namespace) -> dict[str, float | int]:
    return {
        "nl_exit_tolerance": getattr(args, "grtresna_nl_exit_tolerance", 1.0),
        "nl_stall_tolerance": getattr(args, "grtresna_nl_stall_tolerance", 0.02),
        "gridinit_workers": getattr(args, "grtresna_gridinit_workers", 0),
    }


def solved_ftl_gate_config_from_args(args: argparse.Namespace) -> SolvedFtlGateConfig:
    return SolvedFtlGateConfig(
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


def ensure_radial_recipe_for_grtresna(args: argparse.Namespace, *, label: str) -> None:
    if getattr(args, "grtresna", False) and args.example == "SupportedWormholeCollapse":
        print(f"[{label}] --grtresna requires the RadialRecipe example; "
              "switching --example to RadialRecipe.")
        args.example = "RadialRecipe"
