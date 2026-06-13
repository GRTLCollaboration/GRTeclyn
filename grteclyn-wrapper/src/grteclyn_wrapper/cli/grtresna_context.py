"""Build GRTresna search space and solver config from CLI args."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Any

from ..grtresna.domain import GRTresnaDomainConfig
from ..grtresna.solver import GRTresnaConfig
from ..search.grtresna_convergence_gate import GRTresnaConvergenceConfig
from ..search.optimize import ANGULAR_BASE_OVERRIDES, SearchDimension, build_search_space
from ..search.solved_ftl_gate import SolvedFtlGateConfig
from .grtresna_args import grtresna_speed_kwargs, solved_ftl_gate_config_from_args


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


def build_grtresna_search_context(
    args: argparse.Namespace,
    base_overrides: dict[str, Any],
) -> GRTresnaSearchContext:
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
    overrides = dict(base_overrides)
    if use_grtresna and grtresna_ansatz == "ring":
        overrides = {**overrides, "grtresna_ring_lumps": grtresna_lumps}
    if use_grtresna and grtresna_ansatz == "shell":
        overrides = {**overrides, "grtresna_shell_lumps": grtresna_lumps}
    if nonspherical and not use_grtresna:
        overrides = {**ANGULAR_BASE_OVERRIDES, **overrides}

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
