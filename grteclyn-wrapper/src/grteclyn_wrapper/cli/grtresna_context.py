"""Build GRTresna search space and solver config from CLI args."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Any

from ..grtresna.domain import GRTresnaDomainConfig
from ..grtresna.matter.models import matter_selection_base_overrides
from ..grtresna.solver import GRTresnaConfig
from ..search.grtresna_convergence_gate import GRTresnaConvergenceConfig
from ..search.optimize import ANGULAR_BASE_OVERRIDES, SearchDimension, build_search_space
from ..search.optimize.spaces import restrict_retrograde_omega
from ..search.solved_ftl_gate import SolvedFtlGateConfig
from .grtresna_args import (
    grtresna_speed_kwargs,
    resolve_grtresna_matter_from_args,
    solved_ftl_gate_config_from_args,
)


@dataclass
class GRTresnaSearchContext:
    search_space: list[SearchDimension]
    base_overrides: dict[str, Any]
    use_grtresna: bool
    grtresna_config: GRTresnaConfig | None
    solved_ftl_gate_config: SolvedFtlGateConfig | None
    grtresna_convergence_config: GRTresnaConvergenceConfig | None


def _grtresna_domain_from_args(args: argparse.Namespace) -> GRTresnaDomainConfig:
    return GRTresnaDomainConfig(
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


def _retrograde_only(overrides: dict[str, Any]) -> bool:
    """Mirrors the decoder's reading of ``trajectory_retrograde_only``."""
    try:
        return bool(int(float(overrides.get("trajectory_retrograde_only", 0))))
    except (TypeError, ValueError):
        return False


def build_grtresna_search_context(
    args: argparse.Namespace,
    base_overrides: dict[str, Any],
) -> GRTresnaSearchContext:
    nonspherical = getattr(args, "nonspherical", False)
    use_grtresna = getattr(args, "grtresna", False)
    grtresna_lumps = getattr(args, "grtresna_lumps", 5)
    grtresna_ansatz = getattr(args, "grtresna_ansatz", "free")
    grtresna_shell_profile = getattr(args, "grtresna_shell_profile", "compact")
    matter = resolve_grtresna_matter_from_args(args)

    # Detect whether shell_static is pinned (from --pin-dimension or base).
    # If pinned to 1, velocity dims are excluded from the search space.
    shell_static = True
    boson_allow_exotic = False
    pin_specs = getattr(args, "pin_dimension", []) or []
    for spec in pin_specs:
        key, sep, val = spec.partition("=")
        if key.strip() == "grtresna_shell_static" and sep:
            shell_static = int(round(float(val))) >= 1

    import os

    if os.environ.get("GRTRESNA_BOSON_ALLOW_EXOTIC", "0").lower() in ("1", "true", "yes"):
        boson_allow_exotic = True

    boson_matched_bounds = os.environ.get(
        "GRTRESNA_BOSON_MATCHED_BOUNDS", "0"
    ).lower() in ("1", "true", "yes")

    search_space = build_search_space(
        nonspherical=nonspherical,
        grtresna=use_grtresna,
        grtresna_lumps=grtresna_lumps,
        grtresna_ansatz=grtresna_ansatz,
        grtresna_shell_profile=grtresna_shell_profile,
        grtresna_matter_sector=matter.sector,
        grtresna_shell_static=shell_static,
        grtresna_boson_allow_exotic=boson_allow_exotic,
        grtresna_boson_matched_bounds=boson_matched_bounds,
    )
    # A retrograde-only campaign folds positive omega back at decode, so the
    # upper half of every omega axis is a mirror of the lower half and the
    # optimizer is fitting a shape that does not exist.  Search the half that
    # survives decoding.  DebugPreGPU.md PG-7.
    search_space = restrict_retrograde_omega(
        search_space,
        enabled=_retrograde_only(base_overrides),
    )
    overrides = dict(base_overrides)
    if use_grtresna and matter.is_scalar and grtresna_ansatz == "ring":
        overrides = {**overrides, "grtresna_ring_lumps": grtresna_lumps}
    if use_grtresna and grtresna_ansatz == "shell":
        overrides = {**overrides, "grtresna_shell_lumps": grtresna_lumps}
    if use_grtresna and grtresna_ansatz == "trajectory":
        overrides = {
            **overrides,
            "trajectory_mode": 1,
            "trajectory_num_lumps": grtresna_lumps,
        }
    if use_grtresna and grtresna_ansatz == "sh":
        overrides = {**overrides, "grtresna_sh_lumps": grtresna_lumps}
    if use_grtresna:
        # Matter-selection defaults first; campaign EXTRA_SETS / --set win so an
        # explicit grtresna_matter_model=grtresna_bicomplex_scalar is preserved.
        overrides = {**matter_selection_base_overrides(matter), **overrides}
    if boson_allow_exotic:
        overrides["grtresna_boson_allow_exotic"] = 1.0
    if nonspherical and not use_grtresna:
        overrides = {**ANGULAR_BASE_OVERRIDES, **overrides}

    pins: dict[str, float] = {}
    for spec in getattr(args, "pin_dimension", []) or []:
        key, sep, val = spec.partition("=")
        key = key.strip()
        if not key or not sep:
            raise ValueError(f"--pin-dimension expects KEY=VALUE, got {spec!r}")
        try:
            pins[key] = float(val)
        except ValueError as exc:
            raise ValueError(
                f"--pin-dimension value for {key!r} must be numeric, got {val!r}"
            ) from exc

    if pins:
        unknown = pins.keys() - {d.param_key for d in search_space}
        if unknown:
            raise ValueError(
                f"--pin-dimension keys not in search space: {sorted(unknown)}"
            )
        search_space = [d for d in search_space if d.param_key not in pins]
        overrides = {**overrides, **pins}

    grtresna_config = None
    solved_ftl_gate_config = None
    grtresna_convergence_config = None

    if use_grtresna:
        grtresna_domain = _grtresna_domain_from_args(args)
        overrides = {**overrides, **grtresna_domain.evolution_overrides()}
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
            bh1_bare_mass=0.0,
            bh1_spin=(0.0, 0.0, 0.0),
            bh2_bare_mass=0.0,
            dphi=0.0,
            dpi=0.0,
            cleanup=not getattr(args, "grtresna_keep_source", False),
            **grtresna_speed_kwargs(args),
        )
        grtresna_config = grtresna_domain.apply_to_solver(grtresna_config)
        solved_ftl_gate_config = solved_ftl_gate_config_from_args(args)
        grtresna_convergence_config = GRTresnaConvergenceConfig(
            max_ham_pct=getattr(args, "grtresna_max_ham_pct", 5.0),
            max_mom_pct=getattr(args, "grtresna_max_mom_pct", 5.0),
        )

    return GRTresnaSearchContext(
        search_space=search_space,
        base_overrides=overrides,
        use_grtresna=use_grtresna,
        grtresna_config=grtresna_config,
        solved_ftl_gate_config=solved_ftl_gate_config,
        grtresna_convergence_config=grtresna_convergence_config,
    )
