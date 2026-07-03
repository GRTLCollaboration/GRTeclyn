"""Assess candidates against the falsification tier ladder."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from .types import (
    DEFAULT_TIER_CONFIG,
    FAIL,
    PASS,
    UNAVAILABLE,
    GateResult,
    Tier,
    TierAssessment,
    TierConfig,
    TIER_NAMES,
)


def _f(components: Mapping[str, float], key: str, default: float = 0.0) -> float:
    try:
        return float(components.get(key, default))
    except (TypeError, ValueError):
        return default


def survivor_objectives(components: Mapping[str, float]) -> dict[str, float]:
    """Trade-off axes used for the survivor Pareto front."""
    return {
        "operational_ftl": max(
            _f(components, "operational_ftl"),
            _f(components, "ftl_persistence"),
            _f(components, "operational_ftl_solved"),
        ),
        "anec_condition": _f(components, "anec_condition"),
        "constraint_growth": _f(components, "constraint_growth"),
        "tidal_comfort": _f(components, "tidal_comfort"),
    }


def evaluate_tiers(
    components: Mapping[str, float],
    *,
    metrics: Mapping[str, Any] | None = None,
    extra: Mapping[str, Any] | None = None,
    config: TierConfig = DEFAULT_TIER_CONFIG,
) -> TierAssessment:
    """Assess one candidate against the falsification ladder."""
    del metrics  # reserved for future metric-aware gates
    extra = extra or {}
    gates: list[GateResult] = []

    def add(tier: Tier, status: str, detail: str = "") -> None:
        gates.append(GateResult(int(tier), TIER_NAMES[tier], status, detail))

    if not components:
        add(Tier.CONSTRUCTED, FAIL, "no score components (rejected / preflight)")
    else:
        q = max(_f(components, "initial_constraint_quality"), _f(components, "constraint_health"))
        if q >= config.constructed_floor:
            add(Tier.CONSTRUCTED, PASS, f"constraint_quality={q:.3g}")
        else:
            add(Tier.CONSTRUCTED, FAIL, f"constraint_quality={q:.3g} < {config.constructed_floor}")

    signal = max(
        _f(components, "operational_ftl_solved"),
        _f(components, "operational_ftl"),
        _f(components, "ftl_precursor"),
        _f(components, "ftl_shortcut"),
        _f(components, "nonflat_geometry"),
        _f(components, "ftl_persistence"),
        _f(components, "expansion_asymmetry"),
    )
    # GW-beam campaigns emit no FTL/collapse signal by construction; their
    # "nontrivial" signature is genuine gravitational-wave emission.  Let a
    # meaningful gw_beam_quality (log-scaled power) clear the NONTRIVIAL bar
    # so GW emitters are not uniformly tier-rejected.
    gw_signal = _f(components, "gw_beam_quality")
    if gw_signal >= config.gw_emitter_floor:
        signal = max(signal, gw_signal)
    if signal >= config.nontrivial_floor:
        add(Tier.NONTRIVIAL, PASS, f"max_signal={signal:.3g}")
    else:
        add(Tier.NONTRIVIAL, FAIL, f"max_signal={signal:.3g} < {config.nontrivial_floor} (flat/trivial)")

    has_evolved = ("operational_ftl" in components) or ("ftl_persistence" in components)
    if not has_evolved:
        add(Tier.OPERATIONAL, UNAVAILABLE, "no evolved FTL diagnostic (not evolved yet)")
    else:
        op = max(_f(components, "operational_ftl"), _f(components, "ftl_persistence"))
        horizon = _f(components, "horizon_penalty")
        constr = _f(components, "constraint_health")
        if op < config.operational_floor:
            add(Tier.OPERATIONAL, FAIL, f"evolved F_op={op:.3g} < {config.operational_floor}")
        elif horizon < config.max_trapped:
            add(Tier.OPERATIONAL, FAIL, f"trapped-surface crossing (horizon={horizon:.3g})")
        elif constr < config.constraint_floor:
            add(Tier.OPERATIONAL, FAIL, f"constraint_health={constr:.3g} < {config.constraint_floor}")
        else:
            if "ftl_geo_evolving" in components:
                geo = _f(components, "ftl_geo_evolving")
                geo_label = "f_geo_evolving"
            else:
                geo = _f(components, "operational_ftl_geodesic")
                geo_label = "f_geo"
            if (
                "ftl_geo_evolving" not in components
                and "operational_ftl_geodesic" not in components
            ):
                add(
                    Tier.OPERATIONAL,
                    PASS,
                    f"evolved F_op={op:.3g}, geodesic unchecked, horizon={horizon:.3g}",
                )
            elif op < config.operational_floor:
                add(Tier.OPERATIONAL, FAIL, f"evolved F_op={op:.3g} < {config.operational_floor}")
            elif geo < config.geodesic_floor and op >= config.operational_floor:
                add(
                    Tier.OPERATIONAL,
                    FAIL,
                    f"coordinate FTL without gauge-invariant shortcut "
                    f"(F_op={op:.3g}, {geo_label}={geo:.3g})",
                )
            elif geo < config.geodesic_floor:
                add(
                    Tier.OPERATIONAL,
                    FAIL,
                    f"{geo_label}={geo:.3g} < {config.geodesic_floor}",
                )
            else:
                add(
                    Tier.OPERATIONAL,
                    PASS,
                    f"evolved F_op={op:.3g}, {geo_label}={geo:.3g}, horizon={horizon:.3g}",
                )

    surv = _f(components, "survival")
    stab = max(_f(components, "stability"), _f(components, "comoving_stability"))
    grow = _f(components, "constraint_growth")
    persist = _f(components, "ftl_persistence", _f(components, "operational_ftl"))
    if surv < config.survival_floor:
        add(Tier.PERSISTENT, FAIL, f"survival={surv:.3g} < {config.survival_floor}")
    elif stab < config.stability_floor:
        add(Tier.PERSISTENT, FAIL, f"stability={stab:.3g} < {config.stability_floor}")
    elif grow < config.growth_floor:
        add(Tier.PERSISTENT, FAIL, f"constraint_growth={grow:.3g} < {config.growth_floor}")
    elif persist < config.persistence_floor:
        add(Tier.PERSISTENT, FAIL, f"channel did not persist (F_op_ev={persist:.3g})")
    else:
        add(Tier.PERSISTENT, PASS, f"survival={surv:.3g}, stability={stab:.3g}, growth={grow:.3g}")

    if "energy_condition" not in components:
        add(Tier.OBSERVER_EC, UNAVAILABLE, "no observer energy-condition diagnostic")
    else:
        ec = _f(components, "energy_condition")
        exotic = _f(components, "exotic_penalty")
        qei = _f(components, "qei_penalty")
        if ec < config.ec_floor:
            add(Tier.OBSERVER_EC, FAIL, f"energy_condition={ec:.3g} < {config.ec_floor}")
        elif exotic < -config.exotic_cap:
            add(Tier.OBSERVER_EC, FAIL, f"exotic cost too large (penalty={exotic:.3g})")
        elif qei < -config.qei_cap:
            add(Tier.OBSERVER_EC, FAIL, f"QEI violation too large (penalty={qei:.3g})")
        else:
            add(
                Tier.OBSERVER_EC,
                PASS,
                f"energy_condition={ec:.3g}, exotic_penalty={exotic:.3g}, qei_penalty={qei:.3g}",
            )

    rc = extra.get("resolution_converged")
    if rc is None:
        add(Tier.CONVERGED, UNAVAILABLE, "no resolution-ladder replay")
    else:
        add(Tier.CONVERGED, PASS if rc else FAIL, f"resolution_converged={bool(rc)}")

    af = extra.get("analytic_form")
    if af is None:
        add(Tier.ANALYTIC, UNAVAILABLE, "no closed-form back-derivation")
    else:
        add(Tier.ANALYTIC, PASS if af else FAIL, f"analytic_form={bool(af)}")

    reached = Tier.REJECTED
    for g in gates:
        if g.status == PASS:
            reached = Tier(g.tier)
        else:
            break

    return TierAssessment(
        tier=int(reached),
        tier_name=TIER_NAMES[reached],
        gates=gates,
        objectives=survivor_objectives(components),
    )


def tier_histogram(assessments: Sequence[TierAssessment]) -> dict[str, int]:
    """Count how many candidates reached each tier."""
    hist: dict[str, int] = {name: 0 for name in TIER_NAMES.values()}
    for a in assessments:
        hist[a.tier_name] += 1
    return hist
