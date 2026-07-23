"""Stationary geometry scoring: frozen f_geo, stationary f_ff, energy cost."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from ...metrics.probes.ftl.geodesic import (
    GeodesicFtlReport,
    build_metric_3d_from_gridinit,
    geodesic_report_best_direction,
    geodesic_report_from_metric_g,
)
from ...metrics.probes.ftl.metric_field import EvolvingMetricField
from ...metrics.probes.ftl.observer_timing import (
    FreeFallTimingReport,
    compute_freefall_observer_timing,
)
from .genome import GeometryGenome
from .render import RenderConfig, RenderedGeometry, render_and_write, render_genome


@dataclass
class GeometryAtlasEvaluation:
    """Full Stage-1 evaluation of one stationary geometry genome."""

    eval_id: int
    score: float
    descriptors: tuple[float, float]
    cell: tuple[int, int]
    rejected: bool
    reject_reason: str | None
    f_geo: float
    f_ff: float
    integral_negative_rho: float
    min_rho: float
    h_quality_ok: bool
    freefall_reached: bool
    diagnostics: dict[str, float] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)
    genome: dict | None = None
    gridinit_path: str | None = None

    def to_dict(self) -> dict:
        return {
            "eval_id": self.eval_id,
            "score": self.score,
            "descriptors": list(self.descriptors),
            "cell": list(self.cell),
            "rejected": self.rejected,
            "reject_reason": self.reject_reason,
            "f_geo": self.f_geo,
            "f_ff": self.f_ff,
            "integral_negative_rho": self.integral_negative_rho,
            "min_rho": self.min_rho,
            "h_quality_ok": self.h_quality_ok,
            "freefall_reached": self.freefall_reached,
            "diagnostics": self.diagnostics,
            "notes": list(self.notes),
            "genome": self.genome,
            "gridinit_path": self.gridinit_path,
        }


def stationary_field_from_g(
    g: NDArray[np.float64],
    origin: NDArray[np.float64],
    spacing: Sequence[float],
    *,
    n_time: int = 5,
    dt: float = 1.0,
) -> EvolvingMetricField:
    """Wrap a frozen 4-metric as a time-independent evolving field stack."""
    g_stack = np.stack([g] * n_time, axis=0)
    times = np.arange(n_time, dtype=float) * dt
    return EvolvingMetricField(
        g_stack=g_stack,
        times=times,
        origin=np.asarray(origin, dtype=float),
        spacing=(float(dt), float(spacing[0]), float(spacing[1]), float(spacing[2])),
    )


def probe_half_length(
    *,
    support_radius: float,
    box_length: float,
    scale: float = 1.5,
) -> float:
    """Localise emitter/detector about the support rather than the full box.

    Caps at 40% of the box so rays remain inside the domain margin band.
    """
    desired = max(float(support_radius) * float(scale), 4.0)
    return float(min(desired, 0.4 * float(box_length)))


def score_metric_g(
    g: NDArray[np.float64],
    origin: NDArray[np.float64],
    spacing: Sequence[float],
    *,
    n_rays: int = 5,
    directions: Sequence[str] = ("x", "y", "z"),
    emission_tau: float = 2.0,
    compute_ff: bool = True,
    half_length: float | None = None,
) -> tuple[GeodesicFtlReport, FreeFallTimingReport | None]:
    """Compute frozen null and (optional) stationary free-fall scores."""
    if len(directions) > 1:
        geo = geodesic_report_best_direction(
            g,
            origin,
            spacing,
            n_rays=n_rays,
            h_tol=1.0e-6,
            directions=directions,
            half_length=half_length,
        )
    else:
        axis = {"x": 0, "y": 1, "z": 2}.get(directions[0], 0)
        geo = geodesic_report_from_metric_g(
            g,
            origin,
            spacing,
            n_rays=n_rays,
            h_tol=1.0e-6,
            axis=axis,
            half_length=half_length,
        )

    ff: FreeFallTimingReport | None = None
    if compute_ff:
        # Need enough coordinate time for emission + transit on a stationary stack.
        t_flat = float(geo.t_flat) if geo.t_flat else float(spacing[0] * g.shape[0])
        dt = max(0.5, t_flat / 8.0)
        n_time = max(5, int(np.ceil((emission_tau + 1.5 * t_flat) / dt)) + 2)
        field = stationary_field_from_g(g, origin, spacing, n_time=n_time, dt=dt)
        try:
            ff = compute_freefall_observer_timing(field, emission_tau=emission_tau)
        except Exception as exc:  # noqa: BLE001 — probe failures become soft zeros
            ff = FreeFallTimingReport(
                emission_tau=emission_tau,
                emission_t=0.0,
                reception_tau=None,
                reception_t=None,
                initial_proper_separation=0.0,
                flat_reception_tau=0.0,
                receiver_clock_advance=None,
                fractional_arrival_advance=0.0,
                miss_distance=1.0e30,
                reception_tolerance=0.0,
                reached=False,
                max_null_relative_drift=0.0,
                emitter_mass_shell_drift=0.0,
                receiver_mass_shell_drift=0.0,
                emitter_displacement=(0.0, 0.0, 0.0),
                receiver_displacement=(0.0, 0.0, 0.0),
                optimizer_success=False,
                optimizer_evaluations=0,
                notes=(f"freefall failed: {exc}",),
            )
    return geo, ff


def descriptor_axes(
    morphology: float,
    integral_negative_rho: float,
    *,
    energy_floor: float = 1.0e-6,
    energy_ceiling: float = 1.0e2,
) -> tuple[float, float]:
    """MAP-Elites axes: geometry morphology (shift fraction) vs log exotic cost.

    ``morphology`` is the shift-channel fraction in ``[0, 1]`` so that distinct
    geometry families (shift tubes vs conformal lenses) occupy separate cells.
    Shortcut strength (``f_geo``) is the objective, not a diversity axis.
    """
    x = float(np.clip(morphology, 0.0, 1.0))
    e = max(float(integral_negative_rho), energy_floor)
    log_e = np.log10(e)
    log_lo = np.log10(energy_floor)
    log_hi = np.log10(energy_ceiling)
    y = float(np.clip((log_e - log_lo) / (log_hi - log_lo), 0.0, 1.0))
    return x, y


def cell_from_descriptors(
    descriptors: tuple[float, float], bins: int
) -> tuple[int, int]:
    ix = min(bins - 1, max(0, int(descriptors[0] * bins)))
    iy = min(bins - 1, max(0, int(descriptors[1] * bins)))
    return ix, iy


def _rank_score(
    *,
    f_ff: float,
    f_geo: float,
    freefall_reached: bool,
    h_quality_ok: bool,
    diagnostics: dict[str, float],
) -> tuple[float, bool, str | None, list[str]]:
    """Within-cell ranking score with hard rejection gates."""
    notes: list[str] = []
    # Hard rejects.
    if diagnostics.get("min_lapse", 1.0) <= 1.0e-3:
        return -1.0e9, True, "nonpositive_lapse", notes
    if diagnostics.get("min_gamma_eig", 1.0) <= 0.0:
        return -1.0e9, True, "metric_not_spd", notes
    if diagnostics.get("boundary_max_abs_alpha_m1", 0.0) > 0.05:
        return -1.0e9, True, "boundary_not_flat_alpha", notes
    if diagnostics.get("boundary_max_abs_beta", 0.0) > 0.05:
        return -1.0e9, True, "boundary_not_flat_beta", notes
    if diagnostics.get("ham_l2", 0.0) > 5.0:
        return -1.0e9, True, "constraint_inconsistent", notes

    if not h_quality_ok:
        notes.append("null_h_drift_high")
    if not freefall_reached:
        notes.append("freefall_missed")

    # Primary: valid free-fall fractional advance; secondary: frozen f_geo.
    ff_term = f_ff if freefall_reached else 0.0
    geo_term = f_geo if h_quality_ok else 0.5 * f_geo
    # Soft penalty for large constraint residuals / roughness.
    ham_pen = 0.01 * float(diagnostics.get("ham_l2", 0.0))
    mom_pen = 0.01 * float(diagnostics.get("mom_l2", 0.0))
    score = 1000.0 * ff_term + 100.0 * geo_term - ham_pen - mom_pen
    return score, False, None, notes


def evaluate_genome(
    genome: GeometryGenome,
    *,
    eval_id: int,
    render_cfg: RenderConfig,
    out_dir: Path | None = None,
    bins: int = 8,
    n_rays: int = 3,
    compute_ff: bool = True,
    keep_gridinit: bool = False,
    localise_probe: bool = True,
) -> GeometryAtlasEvaluation:
    """Render, score, and descriptor-bin one genome."""
    rendered: RenderedGeometry | None = None
    gridinit_path: str | None = None
    notes: list[str] = []

    try:
        if out_dir is not None and keep_gridinit:
            out_dir.mkdir(parents=True, exist_ok=True)
            path = out_dir / f"eval_{eval_id:06d}.gridinit"
            rendered, written = render_and_write(genome, path, render_cfg)
            gridinit_path = str(written)
            g, origin, spacing = build_metric_3d_from_gridinit(written)
        else:
            rendered = render_genome(genome, render_cfg)
            # Build g directly from ADM without a temp file when possible.
            g, origin, spacing = _metric_from_rendered(rendered, render_cfg)
            if out_dir is not None:
                # Always write a lightweight JSON later; optional gridinit for elites.
                pass
    except Exception as exc:  # noqa: BLE001
        descriptors = (0.0, 0.0)
        return GeometryAtlasEvaluation(
            eval_id=eval_id,
            score=-1.0e9,
            descriptors=descriptors,
            cell=cell_from_descriptors(descriptors, bins),
            rejected=True,
            reject_reason=f"render_failed: {exc}",
            f_geo=0.0,
            f_ff=0.0,
            integral_negative_rho=0.0,
            min_rho=0.0,
            h_quality_ok=False,
            freefall_reached=False,
            diagnostics={},
            notes=[str(exc)],
            genome=genome.to_dict(),
            gridinit_path=None,
        )

    assert rendered is not None
    diag = dict(rendered.diagnostics)
    half_length = None
    if localise_probe:
        half_length = probe_half_length(
            support_radius=float(genome.config.support_radius),
            box_length=float(render_cfg.L),
        )
        diag["probe_half_length"] = float(half_length)
    geo, ff = score_metric_g(
        g,
        origin,
        spacing,
        n_rays=n_rays,
        compute_ff=compute_ff,
        half_length=half_length,
    )
    # Always keep the raw timing advantage.  Hard-zeroing on h-drift previously
    # hid strong Alcubierre channels (v≳1) behind flat transverse axes.
    f_geo_raw = float(geo.f_geo)
    h_ok = bool(geo.h_quality_ok)
    f_geo = f_geo_raw if h_ok else 0.5 * f_geo_raw
    if not h_ok and f_geo_raw > 0.0:
        diag["f_geo_raw_untrusted"] = f_geo_raw
        notes.extend(geo.notes)
        notes.append("f_geo_softened_for_h_drift")

    f_ff = float(ff.fractional_arrival_advance) if ff is not None else 0.0
    freefall_reached = bool(ff.reached) if ff is not None else False
    if ff is not None:
        diag["freefall_miss_distance"] = float(ff.miss_distance)
        diag["freefall_max_null_rel_drift"] = float(ff.max_null_relative_drift)
        notes.extend(ff.notes)

    diag["f_geo"] = f_geo
    diag["f_geo_raw"] = f_geo_raw
    diag["f_ff"] = f_ff
    diag["geo_n_reached"] = float(geo.n_reached)
    diag["geo_max_h_rel_drift"] = float(geo.max_h_rel_drift)

    descriptors = descriptor_axes(
        float(diag.get("shift_fraction", 0.0)),
        float(diag.get("integral_negative_rho", 0.0)),
    )
    cell = cell_from_descriptors(descriptors, bins)
    score, rejected, reason, rank_notes = _rank_score(
        f_ff=f_ff,
        f_geo=f_geo,
        freefall_reached=freefall_reached,
        h_quality_ok=h_ok,
        diagnostics={**diag, "min_lapse": float(rendered.grid.metrics["min_lapse"])},
    )
    notes.extend(rank_notes)

    return GeometryAtlasEvaluation(
        eval_id=eval_id,
        score=score,
        descriptors=descriptors,
        cell=cell,
        rejected=rejected,
        reject_reason=reason,
        f_geo=f_geo,
        f_ff=f_ff,
        integral_negative_rho=float(diag.get("integral_negative_rho", 0.0)),
        min_rho=float(diag.get("min_rho", 0.0)),
        h_quality_ok=bool(geo.h_quality_ok),
        freefall_reached=freefall_reached,
        diagnostics=diag,
        notes=notes,
        genome=genome.to_dict(),
        gridinit_path=gridinit_path,
    )


def _metric_from_rendered(
    rendered: RenderedGeometry,
    cfg: RenderConfig,
) -> tuple[NDArray[np.float64], NDArray[np.float64], tuple[float, float, float]]:
    """Build geodesic metric array directly from rendered ADM fields."""
    alpha = rendered.alpha  # (nz, ny, nx)
    beta = rendered.beta
    gamma = rendered.gamma
    # Transpose to (nx, ny, nz, ...)
    alpha_t = np.transpose(alpha, (2, 1, 0))
    beta_t = np.transpose(beta, (2, 1, 0, 3))
    gamma_t = np.transpose(gamma, (2, 1, 0, 3, 4))

    beta_low = np.einsum("...ij,...j->...i", gamma_t, beta_t)
    beta_sq = np.einsum("...i,...i->...", beta_low, beta_t)
    g = np.zeros(alpha_t.shape + (4, 4), dtype=float)
    g[..., 0, 0] = beta_sq - alpha_t * alpha_t
    g[..., 0, 1:] = beta_low
    g[..., 1:, 0] = beta_low
    g[..., 1:, 1:] = gamma_t

    origin = np.asarray(cfg.origin, dtype=float) + 0.5 * cfg.dx
    spacing = (cfg.dx, cfg.dx, cfg.dx)
    return g, origin, spacing


def evaluation_summary(ev: GeometryAtlasEvaluation) -> dict:
    """Compact summary for archive / logging."""
    return {
        "eval_id": ev.eval_id,
        "score": ev.score,
        "cell": list(ev.cell),
        "descriptors": list(ev.descriptors),
        "f_geo": ev.f_geo,
        "f_ff": ev.f_ff,
        "integral_negative_rho": ev.integral_negative_rho,
        "min_rho": ev.min_rho,
        "rejected": ev.rejected,
        "reject_reason": ev.reject_reason,
        "h_quality_ok": ev.h_quality_ok,
        "freefall_reached": ev.freefall_reached,
        "gridinit_path": ev.gridinit_path,
    }


__all__ = [
    "GeometryAtlasEvaluation",
    "cell_from_descriptors",
    "descriptor_axes",
    "evaluate_genome",
    "evaluation_summary",
    "probe_half_length",
    "score_metric_g",
    "stationary_field_from_g",
]
