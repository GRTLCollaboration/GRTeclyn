"""Falsification-tier layer for the MAP-Elites Spacetime Failure Atlas.

A quality-diversity archive answers *"what is the best candidate of each
kind?"*, but the score is only a proxy: a high-scoring elite may be a gauge
artifact, a short-stop-time survivor, or a degenerate-metric numerical blow-up.
The search therefore cannot certify that any elite is *the* FTL solution.  What
it can do is record, for every elite, **how far up an increasingly demanding
falsification ladder it survived**, and track when the archive (and its front of
survivors) has stopped improving.

The ladder is intentionally honest: each tier needs the previous one *and* its
own gate.  Gates whose diagnostics are not present (e.g. no resolution-ladder
replay yet) are reported ``UNAVAILABLE`` rather than silently passed, so a
single short episode climbs only as far as its evidence allows.

Tiers
-----
``T0 CONSTRUCTED``        constraint-satisfying data exists (not rejected).
``T1 NONTRIVIAL``         carries some FTL / non-flat signal (not Minkowski).
``T2 OPERATIONAL``        evolved operational shortcut beats the flat control,
                          bounded constraints, no trapped-surface crossing.
``T3 PERSISTENT``         survives long evolution: stable, non-growing, channel
                          persists (kills short-stop-time exploiters).
``T4 OBSERVER_EC``        observer-robust energy-condition accounting available
                          and exotic cost bounded.
``T5 CONVERGED``          resolution-convergence evidence (external replay flag).
``T6 ANALYTIC``           closed-form back-derivation reproduced the geometry.

A candidate that clears T0..T4 from one campaign and then T5/T6 from promotion
runs is no longer "a high-scoring run" -- it is a result.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from enum import IntEnum
from pathlib import Path
from typing import Any, Mapping, Sequence

from .pareto import ParetoPoint, pareto_front


class Tier(IntEnum):
    """Falsification ladder.  ``REJECTED`` means even T0 failed."""

    REJECTED = -1
    CONSTRUCTED = 0
    NONTRIVIAL = 1
    OPERATIONAL = 2
    PERSISTENT = 3
    OBSERVER_EC = 4
    CONVERGED = 5
    ANALYTIC = 6


TIER_NAMES: dict[int, str] = {
    Tier.REJECTED: "rejected",
    Tier.CONSTRUCTED: "constructed",
    Tier.NONTRIVIAL: "nontrivial",
    Tier.OPERATIONAL: "operational",
    Tier.PERSISTENT: "persistent",
    Tier.OBSERVER_EC: "observer_ec",
    Tier.CONVERGED: "converged",
    Tier.ANALYTIC: "analytic",
}

PASS = "pass"
FAIL = "fail"
UNAVAILABLE = "unavailable"


# Survivor front objectives, all framed as *maximize* (higher reward = better),
# matching the bounded-reward convention in ``metrics/score.py``.
SURVIVOR_OBJECTIVES: tuple[str, ...] = (
    "operational_ftl",   # evolved operational shortcut benefit
    "anec_condition",    # less exotic energy along the axis (reward)
    "constraint_growth", # slower constraint growth (reward)
    "tidal_comfort",     # lower tidal curvature (reward)
)


@dataclass(frozen=True)
class TierConfig:
    """Thresholds for the falsification gates.

    All FTL/constraint floors are deliberately permissive for *discovery*; the
    point of the ladder is to sort candidates, not to pre-judge them.  Tighten
    these for a promotion/validation pass.
    """

    constructed_floor: float = 0.05   # min constraint-quality reward to count as built
    nontrivial_floor: float = 1.0e-3  # min FTL/non-flat signal
    operational_floor: float = 1.0e-3 # min evolved operational FTL benefit
    constraint_floor: float = 0.20    # min constraint-health reward at T2
    max_trapped: float = -0.5         # horizon_penalty must be >= this (no trapped crossing)
    survival_floor: float = 0.99      # full-run survival
    stability_floor: float = 0.25     # min stability reward
    growth_floor: float = 0.50        # min constraint-growth reward (non-growing)
    persistence_floor: float = 1.0e-4 # evolved channel must persist
    ec_floor: float = 1.0e-3          # min observer energy-condition margin
    exotic_cap: float = 6.0           # |exotic_penalty| ceiling at T4


DEFAULT_TIER_CONFIG = TierConfig()


@dataclass(frozen=True)
class GateResult:
    tier: int
    name: str
    status: str
    detail: str = ""


@dataclass(frozen=True)
class TierAssessment:
    tier: int
    tier_name: str
    gates: list[GateResult]
    objectives: dict[str, float]

    def to_dict(self) -> dict[str, Any]:
        return {
            "tier": int(self.tier),
            "tier_name": self.tier_name,
            "gates": [
                {"tier": g.tier, "name": g.name, "status": g.status, "detail": g.detail}
                for g in self.gates
            ],
            "objectives": dict(self.objectives),
        }


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
    """Assess one candidate against the falsification ladder.

    ``components`` is the post-scoring component map (``Evaluation.components``).
    ``extra`` may carry promotion-run evidence that a single episode lacks:
    ``resolution_converged`` (bool) and ``analytic_form`` (truthy) drive T5/T6.
    """
    extra = extra or {}
    gates: list[GateResult] = []

    def add(tier: Tier, status: str, detail: str = "") -> None:
        gates.append(GateResult(int(tier), TIER_NAMES[tier], status, detail))

    # --- T0 CONSTRUCTED -------------------------------------------------
    if not components:
        add(Tier.CONSTRUCTED, FAIL, "no score components (rejected / preflight)")
    else:
        q = max(_f(components, "initial_constraint_quality"), _f(components, "constraint_health"))
        if q >= config.constructed_floor:
            add(Tier.CONSTRUCTED, PASS, f"constraint_quality={q:.3g}")
        else:
            add(Tier.CONSTRUCTED, FAIL, f"constraint_quality={q:.3g} < {config.constructed_floor}")

    # --- T1 NONTRIVIAL --------------------------------------------------
    signal = max(
        _f(components, "operational_ftl_solved"),
        _f(components, "operational_ftl"),
        _f(components, "ftl_precursor"),
        _f(components, "ftl_shortcut"),
        _f(components, "nonflat_geometry"),
        _f(components, "ftl_persistence"),
        _f(components, "expansion_asymmetry"),
    )
    if signal >= config.nontrivial_floor:
        add(Tier.NONTRIVIAL, PASS, f"max_ftl_signal={signal:.3g}")
    else:
        add(Tier.NONTRIVIAL, FAIL, f"max_ftl_signal={signal:.3g} < {config.nontrivial_floor} (flat/trivial)")

    # --- T2 OPERATIONAL -------------------------------------------------
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
            add(Tier.OPERATIONAL, PASS, f"evolved F_op={op:.3g}, horizon={horizon:.3g}")

    # --- T3 PERSISTENT --------------------------------------------------
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

    # --- T4 OBSERVER_EC -------------------------------------------------
    if "energy_condition" not in components:
        add(Tier.OBSERVER_EC, UNAVAILABLE, "no observer energy-condition diagnostic")
    else:
        ec = _f(components, "energy_condition")
        exotic = _f(components, "exotic_penalty")  # <= 0
        if ec < config.ec_floor:
            add(Tier.OBSERVER_EC, FAIL, f"energy_condition={ec:.3g} < {config.ec_floor}")
        elif exotic < -config.exotic_cap:
            add(Tier.OBSERVER_EC, FAIL, f"exotic cost too large (penalty={exotic:.3g})")
        else:
            add(Tier.OBSERVER_EC, PASS, f"energy_condition={ec:.3g}, exotic_penalty={exotic:.3g}")

    # --- T5 CONVERGED (external) ---------------------------------------
    rc = extra.get("resolution_converged")
    if rc is None:
        add(Tier.CONVERGED, UNAVAILABLE, "no resolution-ladder replay")
    else:
        add(Tier.CONVERGED, PASS if rc else FAIL, f"resolution_converged={bool(rc)}")

    # --- T6 ANALYTIC (external) ----------------------------------------
    af = extra.get("analytic_form")
    if af is None:
        add(Tier.ANALYTIC, UNAVAILABLE, "no closed-form back-derivation")
    else:
        add(Tier.ANALYTIC, PASS if af else FAIL, f"analytic_form={bool(af)}")

    # Highest contiguously-passed tier; stop at the first non-PASS gate.
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


@dataclass(frozen=True)
class SurvivorPoint:
    label: str
    tier: int
    score: float
    objectives: dict[str, float]
    episode: str | None


def survivor_front(
    survivors: Sequence[SurvivorPoint],
    *,
    objectives: Sequence[str] = SURVIVOR_OBJECTIVES,
) -> list[SurvivorPoint]:
    """Non-dominated set among gauntlet survivors (the candidate solutions)."""
    points = [
        ParetoPoint(label=s.label, objectives=s.objectives, total=s.score, episode=s.episode)
        for s in survivors
    ]
    front_labels = {p.label for p in pareto_front(points, objectives=objectives)}
    front = [s for s in survivors if s.label in front_labels]
    front.sort(key=lambda s: (-s.tier, -s.score))
    return front


def build_survivors(
    items: Sequence[Mapping[str, Any]],
    *,
    min_tier: int = int(Tier.OPERATIONAL),
) -> list[SurvivorPoint]:
    """Collect survivors at or above ``min_tier`` from assessed records.

    Each item must carry ``label``, ``tier``, ``score``, ``objectives`` and may
    carry ``episode``.
    """
    out: list[SurvivorPoint] = []
    for it in items:
        if int(it.get("tier", Tier.REJECTED)) < min_tier:
            continue
        out.append(SurvivorPoint(
            label=str(it.get("label")),
            tier=int(it.get("tier", Tier.REJECTED)),
            score=float(it.get("score", 0.0)),
            objectives=dict(it.get("objectives", {})),
            episode=it.get("episode"),
        ))
    return out


def assess_campaign(
    campaign_dir: Path,
    *,
    min_tier: int = int(Tier.OPERATIONAL),
    config: TierConfig = DEFAULT_TIER_CONFIG,
    write: bool = True,
) -> dict[str, Any]:
    """Assess an existing campaign offline from its ``eval_*/score.json`` files.

    Works for both the CMA-ES optimizer and the MAP-Elites driver, since both
    write per-evaluation ``score.json`` with a ``score.components`` block.  The
    result mirrors the live ``validation.json`` (tier histogram + survivor
    front) and is written into ``campaign_dir`` unless ``write=False``.
    """
    campaign_dir = Path(campaign_dir).expanduser().resolve()
    assessments: list[TierAssessment] = []
    items: list[dict[str, Any]] = []
    for score_path in sorted(campaign_dir.glob("eval_*/score.json")):
        try:
            payload = json.loads(score_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        components = (payload.get("score") or {}).get("components") or {}
        if not components:
            continue
        metrics = payload.get("metrics") or {}
        a = evaluate_tiers(components, metrics=metrics, config=config)
        assessments.append(a)
        items.append({
            "label": score_path.parent.name,
            "tier": a.tier,
            "score": float((payload.get("score") or {}).get("total", 0.0)),
            "objectives": a.objectives,
            "episode": str(score_path.parent),
        })

    survivors = build_survivors(items, min_tier=min_tier)
    front = survivor_front(survivors)
    front_labels = {s.label for s in front}
    result = {
        "campaign": str(campaign_dir),
        "num_evaluated": len(items),
        "tier_histogram": tier_histogram(assessments),
        "survivor_min_tier": min_tier,
        "num_survivors": len(survivors),
        "survivor_front": [
            {
                "label": s.label,
                "tier": s.tier,
                "tier_name": TIER_NAMES.get(Tier(s.tier), str(s.tier)),
                "score": s.score,
                "objectives": s.objectives,
                "episode": s.episode,
            }
            for s in front
        ],
        "candidates": sorted(
            (
                {**it, "tier_name": TIER_NAMES.get(Tier(it["tier"]), str(it["tier"])),
                 "on_front": it["label"] in front_labels}
                for it in items
            ),
            key=lambda d: (-int(d["tier"]), -float(d["score"])),
        ),
    }
    if write:
        (campaign_dir / "validation.json").write_text(
            json.dumps(result, indent=2, sort_keys=True), encoding="utf-8"
        )
    return result


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    union = a | b
    if not union:
        return 1.0
    return len(a & b) / len(union)


def convergence_signals(
    history: Sequence[Mapping[str, Any]],
    *,
    window: int = 3,
    score_eps: float = 1.0e-3,
) -> dict[str, Any]:
    """Decide whether the archive has stopped improving.

    ``history`` is an ordered list of per-iteration snapshots, each with keys
    ``coverage`` (float), ``best_score`` (float) and ``front_labels`` (list of
    survivor-front labels).  The archive is declared *converged* when, over the
    last ``window`` iterations: coverage did not grow, best score improved by
    less than ``score_eps``, and the survivor front membership was stable.
    """
    n = len(history)
    if n == 0:
        return {"converged": False, "reason": "no history"}
    last = history[-1]
    if n < window + 1:
        return {
            "converged": False,
            "reason": f"need {window + 1} iters, have {n}",
            "coverage": last.get("coverage"),
            "best_score": last.get("best_score"),
        }

    recent = history[-(window + 1):]
    coverages = [float(h.get("coverage", 0.0)) for h in recent]
    scores = [float(h.get("best_score", 0.0)) for h in recent]
    coverage_delta = coverages[-1] - coverages[0]
    best_delta = scores[-1] - scores[0]

    front_sets = [set(h.get("front_labels", []) or []) for h in recent]
    pair_stability = [
        _jaccard(front_sets[i], front_sets[i + 1]) for i in range(len(front_sets) - 1)
    ]
    front_stable = bool(pair_stability) and min(pair_stability) >= 0.999

    coverage_stalled = coverage_delta <= 0.0
    score_stalled = best_delta < score_eps
    converged = coverage_stalled and score_stalled and front_stable

    return {
        "converged": bool(converged),
        "coverage": coverages[-1],
        "coverage_delta": coverage_delta,
        "best_score": scores[-1],
        "best_score_delta": best_delta,
        "front_stability": min(pair_stability) if pair_stability else 1.0,
        "coverage_stalled": coverage_stalled,
        "score_stalled": score_stalled,
        "front_stable": front_stable,
        "window": window,
    }
