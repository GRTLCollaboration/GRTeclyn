"""GRTresna-related CLI argument helpers."""

from __future__ import annotations

import argparse
import os
from typing import Any

from ..projection.postload_gate import PostLoadGateConfig
from ..search.solved_ftl_gate import SolvedFtlGateConfig


def grtresna_postload_gate_enabled(args: argparse.Namespace, *, use_grtresna: bool) -> bool:
    if not use_grtresna:
        return False
    if getattr(args, "grtresna_postload_gate", False):
        return True
    return os.environ.get("POSTLOAD_GATE", "0") == "1"


def postload_gate_config_from_args(args: argparse.Namespace) -> PostLoadGateConfig:
    return PostLoadGateConfig(
        max_hamiltonian_l2=getattr(args, "postload_max_ham_l2", 1.0e-2),
        max_momentum_l2=getattr(args, "postload_max_mom_l2", 1.0e-2),
    )


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
