"""Aggregate all episode metrics from diagnostics and probes."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np

from ..diagnostics import (
    read_collapse_metrics,
    read_comoving_metrics,
    read_constraint_metrics,
    read_curvature_invariant_metrics,
    read_energy_condition_metrics,
    read_ftl_timeseries_metrics,
    read_growth_metrics,
    read_qei_metrics,
    read_stability_metrics,
    read_transport_metrics,
)
from ..io.plotfile_lock import PLOTFILE_READ_LOCK
from ..probes.boundary import read_boundary_flux_metrics
from ..probes.ftl.analytic import compute_ftl_metrics, load_overrides_from_episode
from ..probes.ftl.general import (
    compute_general_ftl,
    compute_general_ftl_from_plotfile,
    find_latest_plotfile,
    find_recent_plotfiles,
    matter_coherence_from_plotfile,
    wait_for_plotfile_complete,
)
from ..probes.ftl.geodesic import (
    compute_geodesic_ftl_from_plotfile,
    compute_trajectory_qei_from_plotfile,
)
from ..probes.physical import compute_physical_metrics
from ..probes.warpfactory import effective_energy_conditions_from_plotfiles
from ..types.diagnostics import EffectiveEnergyConditionMetrics, FtlPersistenceMetrics
from ..types.episode import EpisodeMetrics
from .context import EpisodeContext, build_episode_context


def read_episode_metrics(
    episode_dir: Path,
    *,
    ftl_L: float | None = None,
) -> EpisodeMetrics:
    ctx = build_episode_context(episode_dir, ftl_L=ftl_L)

    collapse = read_collapse_metrics(ctx.collapse_path)
    constraints = read_constraint_metrics(ctx.constraint_path)
    stability = read_stability_metrics(ctx.collapse_path, ctx.areal_path)
    growth = read_growth_metrics(ctx.collapse_path, ctx.constraint_path)
    energy_conditions = read_energy_condition_metrics(ctx.energy_conditions_path)
    curvature = read_curvature_invariant_metrics(ctx.curvature_path)

    general_ftl_evolved = None
    try:
        plotfile = find_latest_plotfile(ctx.episode_dir, complete_only=False)
        if plotfile is not None:
            wait_for_plotfile_complete(plotfile, timeout=30.0)
        plotfile = find_latest_plotfile(ctx.episode_dir)
        if plotfile is not None:
            with PLOTFILE_READ_LOCK:
                general_ftl_evolved = compute_general_ftl_from_plotfile(
                    plotfile, n=129, L=ctx.ftl_L
                )
                coherence = matter_coherence_from_plotfile(plotfile, n=64)
            general_ftl_evolved = replace(
                general_ftl_evolved, structure_coherence=coherence
            )
    except Exception:
        general_ftl_evolved = None

    ftl_persistence = None
    try:
        recent = find_recent_plotfiles(ctx.episode_dir, count=5)
        if len(recent) >= 2:
            f_ops: list[float] = []
            speeds: list[float] = []
            shifts: list[float] = []
            for plotfile in recent:
                try:
                    with PLOTFILE_READ_LOCK:
                        rep = compute_general_ftl_from_plotfile(
                            plotfile, n=97, L=ctx.ftl_L
                        )
                except Exception:
                    continue
                f_ops.append(float(rep.f_op))
                speeds.append(float(rep.max_local_speed))
                shifts.append(float(rep.max_shift))
            if len(f_ops) >= 2:
                ftl_persistence = FtlPersistenceMetrics(
                    n_samples=len(f_ops),
                    f_op_min=float(min(f_ops)),
                    f_op_median=float(np.median(f_ops)),
                    f_op_last=float(f_ops[-1]),
                    max_local_speed_min=float(min(speeds)),
                    max_shift_max=float(max(shifts)),
                )
    except Exception:
        ftl_persistence = None

    geodesic_ftl = None
    try:
        plotfile = find_latest_plotfile(ctx.episode_dir)
        if plotfile is not None and general_ftl_evolved is not None:
            if (
                general_ftl_evolved.f_op > 1.0e-3
                or general_ftl_evolved.max_local_speed > 1.0
            ):
                with PLOTFILE_READ_LOCK:
                    geodesic_ftl = compute_geodesic_ftl_from_plotfile(
                        plotfile, n=65, half_width=ctx.ftl_L
                    )
    except Exception:
        geodesic_ftl = None

    boundary_flux = read_boundary_flux_metrics(ctx.boundary_flux_path)
    if boundary_flux is None:
        boundary_flux = read_boundary_flux_metrics(ctx.boundary_flux_fallback_path)

    # Time-resolved FTL: the gauge-invariant shortcut peaks mid-run and diffuses,
    # so the final frame is half-blind.  This stream (written in-flight by the
    # consumer as it processes+deletes each plotfile) lets the scorer average a
    # composite FTL x stability score over the whole run.
    ftl_timeseries = None
    try:
        ftl_timeseries = read_ftl_timeseries_metrics(ctx.ftl_timeseries_path)
    except Exception:
        ftl_timeseries = None

    effective_ec = None
    try:
        recent = find_recent_plotfiles(ctx.episode_dir, count=5)
        if len(recent) >= 3:
            with PLOTFILE_READ_LOCK:
                rep = effective_energy_conditions_from_plotfiles(
                    [str(p) for p in recent], n_space=32, half_width=ctx.ftl_L
                )
            effective_ec = EffectiveEnergyConditionMetrics(
                nec_min=rep.nec_min,
                wec_min=rep.wec_min,
                nec_slack_min=rep.nec_slack_min,
                rho_eulerian_min=rep.rho_eulerian_min,
                wec_violation_fraction=rep.wec_violation_fraction,
                s_energy_conditions=rep.s_energy_conditions,
                n_points=rep.n_points,
            )
    except Exception:
        effective_ec = None

    ftl = None
    comoving = None
    physical = None
    general_ftl = None
    overrides = load_overrides_from_episode(ctx.episode_dir)
    if overrides:
        try:
            ftl = compute_ftl_metrics(overrides, L=ctx.ftl_L)
        except Exception:
            ftl = None
        try:
            general_ftl = compute_general_ftl(overrides, L=ctx.ftl_L, n=97)
        except Exception:
            general_ftl = None
        try:
            comoving = read_comoving_metrics(ctx.episode_dir, overrides, ftl_L=ctx.ftl_L)
        except Exception:
            comoving = None
        try:
            physical = compute_physical_metrics(overrides, L=ctx.ftl_L)
        except Exception:
            physical = None

    general_ftl_solved = None
    mechanism_descriptor = None
    if ctx.gridinit_path.is_file():
        try:
            from ..probes.ftl.solved import compute_solved_geometry_ftl
            from ...search.solved_ftl_gate import DEFAULT_SOLVED_FTL_GATE_POLICY

            solved = compute_solved_geometry_ftl(ctx.gridinit_path, L=ctx.ftl_L)
            if solved is not None and not DEFAULT_SOLVED_FTL_GATE_POLICY.is_degenerate(solved):
                general_ftl_solved = solved.operational
                mechanism_descriptor = solved.mechanisms.mechanism_descriptor
        except Exception:
            pass

    trajectory_qei = None
    try:
        plotfile = find_latest_plotfile(ctx.episode_dir)
        if plotfile is not None:
            with PLOTFILE_READ_LOCK:
                trajectory_qei = compute_trajectory_qei_from_plotfile(
                    plotfile, n=49, half_width=ctx.ftl_L
                )
    except Exception:
        trajectory_qei = None

    qei = read_qei_metrics(physical, trajectory_violation=trajectory_qei)
    transport = read_transport_metrics(ctx.collapse_path)

    if collapse is None and constraints is None:
        reason = "missing_diagnostics"
    else:
        reason = "completed_or_partial"

    return EpisodeMetrics(
        collapse=collapse,
        constraints=constraints,
        stability=stability,
        comoving=comoving,
        ftl=ftl,
        termination_reason=reason,
        growth=growth,
        physical=physical,
        energy_conditions=energy_conditions,
        curvature=curvature,
        general_ftl=general_ftl,
        general_ftl_evolved=general_ftl_evolved,
        general_ftl_solved=general_ftl_solved,
        mechanism_descriptor=mechanism_descriptor,
        effective_ec=effective_ec,
        geodesic_ftl=geodesic_ftl,
        boundary_flux=boundary_flux,
        qei=qei,
        transport=transport,
        ftl_persistence=ftl_persistence,
        ftl_timeseries=ftl_timeseries,
    )
