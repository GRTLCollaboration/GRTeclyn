from __future__ import annotations

from .splash import SURVIVAL_FORGIVENESS_WELL
from .types import ScoringContext
from .weights import HEALTH_COMPONENTS


def compute_total(
    ctx: ScoringContext,
    objective_mode: str,
    nontriviality: float,
    *,
    splash_mode: str = "discovery",
    exotic_penalty_weight: float = 1.0,
) -> float:
    components = ctx.components
    notes = ctx.notes
    w = ctx.weights

    if objective_mode == "ftl_first":
        return _ftl_first_total(components, notes, exotic_penalty_weight=exotic_penalty_weight)
    if objective_mode == "robust_ftl":
        return _robust_ftl_total(components, notes, exotic_penalty_weight=exotic_penalty_weight)
    if objective_mode == "general_ftl":
        return _general_ftl_total(components, notes, exotic_penalty_weight=exotic_penalty_weight)
    if objective_mode == "critical_collapse":
        return _critical_collapse_total(components, notes, splash_mode=splash_mode)
    return _weighted_total(components, w, nontriviality)


def _ftl_first_total(
    components: dict[str, float],
    notes: list[str],
    *,
    exotic_penalty_weight: float = 1.0,
) -> float:
    # Lexicographic-style scalarization: path-level FTL dominates; precursor
    # and shift are shaping only.  Trapped-surface proxies (horizon_penalty
    # < 0) veto high local-speed artifacts that are not traversable channels.
    health_gate = components.get("nontriviality_gate", 0.0)
    horizon = components.get("horizon_penalty", 0.0)
    trapped_surface = horizon < -0.05
    # Two tiers of FTL signal:
    #   * VALIDATED FTL (geodesic shortcut, evolved operational advantage,
    #     sustained persistence, constraint-solved t=0 shortcut) -- the
    #     actual goal, kept dominant so a real result always wins.
    #   * SHAPING GRADIENTS (channel_progress, ftl_precursor, shift_drive) --
    #     coordinate cone-tilt heuristics that merely point toward FTL.  At
    #     the old weights (channel 350, precursor 60) they dwarfed survival
    #     (20), so a marginally-stable geometry that fragments but tilts its
    #     light cones out-ranked a coherent, persistent survivor -- the
    #     opposite of the stated priority.  They are cut to ~40% so they
    #     still guide the search out of flat space without out-voting the
    #     health terms.
    # HEALTH terms (survival/stability) are correspondingly boosted so a
    # coherent, persistent, stable structure is rewarded on par with a
    # strong shaping signal.  Survival already folds in density retention
    # *and* morphological coherence (see structural_persistence), so a
    # fragmenting end-state cannot bank the larger survival weight.
    # SHAPING-vs-VALIDATED balance.  The honest geodesic recalibration
    # (GEO_FTL_TARGET=2e-1) deliberately makes a realistic few-percent
    # gauge-invariant shortcut score modestly (~160 pts for 3%).  The
    # coordinate-cone shaping rewards (channel_progress + operational_ftl_solved
    # + precursor + shift) must therefore stay *below* that, or a t=0
    # coordinate precursor that yields no gauge-invariant shortcut outranks a
    # validated one (observed in v9: coord-only candidates took #2/#3 over a
    # genuine 3.3% shortcut).  operational_ftl_solved is a *t=0* constraint-
    # solved coordinate signal -- pure shaping -- so it is the most strongly
    # capped; channel_progress is trimmed to keep total shaping subordinate.
    total = (
        1000.0 * components.get("ftl_geo_evolving", 0.0)
        + 1000.0 * components.get("operational_ftl_geodesic", 0.0)
        + 400.0 * components.get("operational_ftl", 0.0)
        + 300.0 * components.get("ftl_persistence", 0.0)
        + 100.0 * components.get("channel_progress", 0.0)
        + 50.0 * components.get("operational_ftl_solved", 0.0)
        + 30.0 * components.get("ftl_precursor", 0.0)
        + 20.0 * components.get("shift_drive", 0.0)
        + 50.0 * components.get("ftl_shortcut", 0.0)
        + 5.0 * components.get("nontrivial_geometry", 0.0)
        + health_gate * (
            70.0 * components.get("survival", 0.0)
            + 6.0 * components.get("constraint_health", 0.0)
            + 10.0 * components.get("stability", 0.0)
            + 15.0 * components.get("instability_penalty", 0.0)
            + 8.0 * components.get("comoving_stability", 0.0)
            + 40.0 * components.get("energy_condition", 0.0)
        )
        # Physicality pressure.  Every v10 top-5 elite was a maximally-exotic
        # warp bubble (energy_condition ~0.05, near-floor exotic_penalty):
        # at the old weights (exotic 1, energy 2) the exotic cost was ~-0.5
        # against a ~150-point shortcut, so the search had *no* incentive to
        # leave the easy exotic corner.  energy_condition is boosted (2 -> 40)
        # to materially reward a genuinely NEC-respecting geometry, and the
        # graded exotic penalty is boosted (1 -> 40): a fully-exotic geometry
        # now costs ~-64 pts (typical ~-20), a real ~20-30% dent against a
        # realistic ~150-250-pt shortcut -- enough that a cleaner candidate of
        # comparable magnitude wins, but NOT enough to swamp the FTL reward.
        # (A first cut at weight 100 floored the whole population negative
        # before any shortcut appeared -- a -160 full-exotic penalty made
        # "minimise exotic" dominate "find FTL", inverting the FTL-first
        # priority.)  The penalty stays graded (0..-1.6), so the QD gradient
        # is preserved.
        + 40.0 * exotic_penalty_weight * components.get("exotic_penalty", 0.0)
        # Moderate stationary penalty: the FTL shaping rewards are already
        # zeroed for a stationary artifact (see the stationary-artifact gate)
        # and the geodesic shortcut is reliability-gated, so a static lens
        # can no longer masquerade as FTL.  The penalty therefore only needs
        # to keep a shift-free geometry ranked below a genuine shift-driven
        # one -- a catastrophic weight (it used to be 250) instead floored the
        # entire population negative and erased the QD gradient.  At this
        # weight a healthy static structure still scores mildly positive
        # while any real shift (graded by |beta|) lifts it further.
        + 8.0 * components.get("stationary_artifact_penalty", 0.0)
        + 500.0 * horizon
    )
    if trapped_surface:
        notes.append(
            f"trapped-surface proxy active (horizon_penalty={horizon:.3f}); "
            "local precursor/shift alone cannot rank this candidate highly"
        )
    notes.append(
        "objective_mode=ftl_first: channel/operational FTL dominate precursor/shift"
    )
    return total


def _robust_ftl_total(
    components: dict[str, float],
    notes: list[str],
    *,
    exotic_penalty_weight: float = 1.0,
) -> float:
    # Robustness-tilted FTL scalarization (option B).  Same FTL-first priority
    # as ``ftl_first`` -- the gauge-invariant geodesic shortcut
    # (operational_ftl_geodesic, 1000x) stays dominant so a genuine result
    # still wins -- but the rest of the budget is rebalanced toward
    # *persistent, healthy, low-exotic* geometries (the OBSERVER_EC survivors)
    # rather than the highest-magnitude-but-exotic/transient peak:
    #   * coordinate-only signals trimmed (operational_ftl 400->200, shaping cut),
    #   * lasting FTL rewarded (ftl_persistence 300->500),
    #   * health boosted (survival 70->150, comoving 8->20, stability 10->25,
    #     energy_condition 40->60),
    #   * exotic penalty hardened 40->70 (a fully-exotic geometry now costs
    #     ~-112 pts, a real ~20-30% dent, kept below the ~100 weight that the
    #     ftl_first comments note floored the whole population negative).
    health_gate = components.get("nontriviality_gate", 0.0)
    horizon = components.get("horizon_penalty", 0.0)
    trapped_surface = horizon < -0.05
    total = (
        1000.0 * components.get("ftl_geo_evolving", 0.0)
        + 1000.0 * components.get("operational_ftl_geodesic", 0.0)
        + 200.0 * components.get("operational_ftl", 0.0)
        + 500.0 * components.get("ftl_persistence", 0.0)
        + 60.0 * components.get("channel_progress", 0.0)
        + 30.0 * components.get("operational_ftl_solved", 0.0)
        + 20.0 * components.get("ftl_precursor", 0.0)
        + 15.0 * components.get("shift_drive", 0.0)
        + 30.0 * components.get("ftl_shortcut", 0.0)
        + 5.0 * components.get("nontrivial_geometry", 0.0)
        + health_gate * (
            150.0 * components.get("survival", 0.0)
            + 10.0 * components.get("constraint_health", 0.0)
            + 25.0 * components.get("stability", 0.0)
            + 15.0 * components.get("instability_penalty", 0.0)
            + 20.0 * components.get("comoving_stability", 0.0)
            + 60.0 * components.get("energy_condition", 0.0)
        )
        + 70.0 * exotic_penalty_weight * components.get("exotic_penalty", 0.0)
        + 8.0 * components.get("stationary_artifact_penalty", 0.0)
        + 500.0 * horizon
    )
    if trapped_surface:
        notes.append(
            f"trapped-surface proxy active (horizon_penalty={horizon:.3f}); "
            "local precursor/shift alone cannot rank this candidate highly"
        )
    notes.append(
        "objective_mode=robust_ftl: persistent/healthy/low-exotic geometry "
        "prioritized; gauge-invariant shortcut stays dominant"
    )
    return total


def _general_ftl_total(
    components: dict[str, float],
    notes: list[str],
    *,
    exotic_penalty_weight: float = 1.0,
) -> float:
    # General-FTL profile: the gauge-invariant 4D null shortcut is the ONLY
    # FTL reward.  All coordinate/warp-motor shaping (shift_drive,
    # channel_progress, operational_ftl_solved, ftl_precursor) is removed so a
    # static wormhole / ring / portal is not penalised for beta^i = 0.  Health
    # and stationarity are boosted because a *persistent, stable* throat is the
    # goal, not a translating bubble.
    health_gate = components.get("nontriviality_gate", 0.0)
    horizon = components.get("horizon_penalty", 0.0)
    total = (
        1000.0 * components.get("ftl_geo_evolving", 0.0)
        + 1000.0 * components.get("operational_ftl_geodesic", 0.0)
        + 600.0 * components.get("ftl_persistence", 0.0)
        + 200.0 * components.get("operational_ftl", 0.0)
        + 5.0 * components.get("nontrivial_geometry", 0.0)
        + health_gate * (
            150.0 * components.get("survival", 0.0)
            + 10.0 * components.get("constraint_health", 0.0)
            + 40.0 * components.get("stability", 0.0)
            + 30.0 * components.get("comoving_stability", 0.0)
            + 15.0 * components.get("instability_penalty", 0.0)
            + 60.0 * components.get("energy_condition", 0.0)
        )
        + 70.0 * exotic_penalty_weight * components.get("exotic_penalty", 0.0)
        + 500.0 * horizon
    )
    notes.append(
        "objective_mode=general_ftl: gauge-invariant shortcut only; "
        "warp-motor shaping and stationary penalty disabled"
    )
    return total


def _critical_collapse_total(
    components: dict[str, float],
    notes: list[str],
    *,
    splash_mode: str = "discovery",
) -> float:
    survival = float(components.get("survival", 0.0))
    well = float(components.get("geometric_curvature_well", 0.0))
    horizon_bonus = float(components.get("horizon_formation_time", 0.0))

    # Successful collapse terminates early.  Forgive survival only on geometric
    # terms — matter corroboration still scales with actual persistence.
    if well > SURVIVAL_FORGIVENESS_WELL or horizon_bonus > 0.0:
        geometric_survival = 1.0
        notes.append("survival_forgiven: critical_collapse_achieved (geometric terms only)")
    else:
        geometric_survival = survival
    matter_survival = survival

    if geometric_survival <= 0.0 and matter_survival <= 0.0:
        notes.append(
            f"objective_mode=critical_collapse splash_mode={splash_mode}: "
            "zero survival zeros score"
        )
        return 0.0
    peak = float(components.get("central_energy_peak", 0.0))
    focus = min(float(components.get("focusing_efficiency", 0.0)), 5.0)
    # Relative growth only counts when absolute peak density is meaningful.
    focus_effective = focus * peak
    lapse_progress = float(components.get("collapse_lapse_progress", 0.0))
    lapse_term = float(components.get("central_lapse_collapse", 0.0))
    dispersion = float(components.get("dispersion_penalty", 0.0))
    pre_penalty = float(components.get("pre_collapsed_penalty", 0.0))
    constraint_quality = float(components.get("constraint_quality", 0.0))

    exotic = float(components.get("exotic_penalty", 0.0))
    # For critical collapse the evolved geometry naturally develops rho < 0
    # regions near the focusing core (gauge artifacts + strong curvature).
    # The full exotic penalty (-1.6 saturated, costing -320 pts) is
    # counterproductive here: it punishes the very dynamics we seek.  Cap
    # the penalty contribution so it provides a mild preference for cleaner
    # data without dominating the score budget.
    exotic_capped = max(exotic, -0.3)

    # --- Spacetime-splash (geometric) terms ---------------------------------
    # The splash we want is a *geometric* event -- a converging gravitational
    # wave that crushes the conformal factor (chi -> 0), spikes |K|, and emits
    # a Weyl/Psi4 pulse at the center -- NOT a slow matter pile-up.  These
    # geometric signatures are therefore the primary reward; the matter
    # (rho) terms below act as secondary corroboration only.
    curvature_well = float(components.get("geometric_curvature_well", 0.0))
    wave_arrival = float(components.get("geometric_wave_arrival", 0.0))
    crunch = float(components.get("geometric_crunch", 0.0))

    total = (
        # Primary: spacetime curvature concentration + GW focusing.
        800.0 * curvature_well * geometric_survival
        + 600.0 * wave_arrival * geometric_survival
        + 300.0 * crunch * geometric_survival
        # Secondary: matter density corroboration of the focusing event.
        + 400.0 * peak * matter_survival
        + 200.0 * focus_effective * matter_survival
        + 100.0 * pre_penalty
        + 100.0 * dispersion
        - 50.0 * (1.0 - constraint_quality) * geometric_survival
        + 50.0 * exotic_capped * geometric_survival
    )
    if splash_mode == "threshold":
        total += 500.0 * lapse_term
    else:
        total += 200.0 * lapse_progress * geometric_survival
        if horizon_bonus > 0.0:
            total += 500.0 * horizon_bonus * geometric_survival

    notes.append(
        f"objective_mode=critical_collapse splash_mode={splash_mode}: "
        "geometric splash (chi-well + Psi4 wave + K-crunch) primary; "
        "rho/focus/lapse secondary; FTL terms ignored"
    )
    return total


def _weighted_total(
    components: dict[str, float],
    weights: dict[str, float],
    nontriviality: float,
) -> float:
    total = 0.0
    for key, value in components.items():
        if key == "nontriviality_gate":
            continue
        gate = nontriviality if key in HEALTH_COMPONENTS else 1.0
        total += weights.get(key, 0.0) * value * gate
    return total
