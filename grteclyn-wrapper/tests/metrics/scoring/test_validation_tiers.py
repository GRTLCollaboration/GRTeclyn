import json

from grteclyn_wrapper.search.validation_tiers import (
    Tier,
    assess_campaign,
    build_survivors,
    convergence_signals,
    evaluate_tiers,
    survivor_front,
    tier_histogram,
)


# --- Tier ladder -------------------------------------------------------


def test_empty_components_is_rejected() -> None:
    a = evaluate_tiers({})
    assert a.tier == int(Tier.REJECTED)
    assert a.tier_name == "rejected"


def test_flat_minkowski_constructed_but_trivial() -> None:
    # Constraint-clean flat space: built, but no FTL signal -> stops at T0.
    flat = {
        "constraint_health": 1.0,
        "initial_constraint_quality": 1.0,
        "operational_ftl": 0.0,
        "ftl_persistence": 0.0,
        "ftl_shortcut": 0.0,
        "nonflat_geometry": 0.0,
        "survival": 1.0,
    }
    a = evaluate_tiers(flat)
    assert a.tier == int(Tier.CONSTRUCTED)


def test_evolved_operational_candidate_reaches_operational() -> None:
    # A candidate with a mild evolved shortcut, clean constraints, no trapped
    # surface should climb to T2 but not T3 (it is not stable/persistent here).
    cand = {
        "constraint_health": 0.8,
        "initial_constraint_quality": 0.6,
        "operational_ftl": 0.046,
        "ftl_persistence": 0.046,
        "ftl_shortcut": 0.0,
        "nonflat_geometry": 0.3,
        "horizon_penalty": 0.0,
        "survival": 1.0,
        "stability": 0.05,   # unstable -> blocks T3
        "constraint_growth": 0.1,
    }
    a = evaluate_tiers(cand)
    assert a.tier == int(Tier.OPERATIONAL)
    # The blocking gate at T3 should be a FAIL, not UNAVAILABLE.
    t3 = next(g for g in a.gates if g.tier == int(Tier.PERSISTENT))
    assert t3.status == "fail"


def test_trapped_surface_blocks_operational() -> None:
    cand = {
        "constraint_health": 0.8,
        "initial_constraint_quality": 0.6,
        "operational_ftl": 0.1,
        "horizon_penalty": -1.0,  # trapped surface crossed
        "survival": 1.0,
    }
    a = evaluate_tiers(cand)
    assert a.tier == int(Tier.NONTRIVIAL)


def test_persistent_candidate_reaches_t3() -> None:
    cand = {
        "constraint_health": 0.9,
        "initial_constraint_quality": 0.8,
        "operational_ftl": 0.02,
        "ftl_persistence": 0.02,
        "horizon_penalty": 0.0,
        "survival": 1.0,
        "stability": 0.85,
        "comoving_stability": 0.85,
        "constraint_growth": 0.9,
    }
    a = evaluate_tiers(cand)
    assert a.tier == int(Tier.PERSISTENT)
    # T4 needs EC diagnostics which are absent -> UNAVAILABLE (not FAIL).
    t4 = next(g for g in a.gates if g.tier == int(Tier.OBSERVER_EC))
    assert t4.status == "unavailable"


def test_external_evidence_drives_t5_t6() -> None:
    cand = {
        "constraint_health": 0.9,
        "initial_constraint_quality": 0.8,
        "operational_ftl": 0.02,
        "ftl_persistence": 0.02,
        "horizon_penalty": 0.0,
        "survival": 1.0,
        "stability": 0.85,
        "constraint_growth": 0.9,
        "energy_condition": 0.5,
        "exotic_penalty": -0.1,
    }
    a = evaluate_tiers(cand, extra={"resolution_converged": True, "analytic_form": "chi(r)=..."})
    assert a.tier == int(Tier.ANALYTIC)


def test_tier_histogram_counts() -> None:
    a1 = evaluate_tiers({})
    a2 = evaluate_tiers({"constraint_health": 1.0, "survival": 1.0})
    hist = tier_histogram([a1, a2])
    assert hist["rejected"] == 1
    assert hist["constructed"] == 1


# --- Survivor front ----------------------------------------------------


def test_survivor_front_filters_by_tier_and_dominance() -> None:
    items = [
        # below survivor tier -> excluded
        {"label": "a", "tier": int(Tier.NONTRIVIAL), "score": 9.0,
         "objectives": {"operational_ftl": 0.5, "anec_condition": 0.5,
                        "constraint_growth": 0.5, "tidal_comfort": 0.5}},
        # dominated by c on every axis -> excluded from front
        {"label": "b", "tier": int(Tier.OPERATIONAL), "score": 5.0,
         "objectives": {"operational_ftl": 0.1, "anec_condition": 0.1,
                        "constraint_growth": 0.1, "tidal_comfort": 0.1}},
        {"label": "c", "tier": int(Tier.OPERATIONAL), "score": 8.0,
         "objectives": {"operational_ftl": 0.9, "anec_condition": 0.9,
                        "constraint_growth": 0.9, "tidal_comfort": 0.9}},
    ]
    survivors = build_survivors(items, min_tier=int(Tier.OPERATIONAL))
    assert {s.label for s in survivors} == {"b", "c"}
    front = survivor_front(survivors)
    assert {s.label for s in front} == {"c"}


# --- Archive convergence ----------------------------------------------


def test_convergence_needs_enough_history() -> None:
    sig = convergence_signals([{"coverage": 0.1, "best_score": 1.0, "front_labels": []}])
    assert sig["converged"] is False


def test_convergence_detected_when_stalled() -> None:
    # Coverage flat, best score flat, front stable across the window -> converged.
    history = [
        {"coverage": 0.5, "best_score": 12.0, "front_labels": ["c"]},
        {"coverage": 0.5, "best_score": 12.0, "front_labels": ["c"]},
        {"coverage": 0.5, "best_score": 12.0, "front_labels": ["c"]},
        {"coverage": 0.5, "best_score": 12.0, "front_labels": ["c"]},
    ]
    sig = convergence_signals(history, window=3)
    assert sig["converged"] is True


def test_convergence_false_when_still_improving() -> None:
    history = [
        {"coverage": 0.2, "best_score": 8.0, "front_labels": ["a"]},
        {"coverage": 0.4, "best_score": 10.0, "front_labels": ["a", "b"]},
        {"coverage": 0.6, "best_score": 12.0, "front_labels": ["c"]},
        {"coverage": 0.8, "best_score": 14.0, "front_labels": ["d"]},
    ]
    sig = convergence_signals(history, window=3)
    assert sig["converged"] is False
    assert sig["coverage_delta"] > 0


# --- Offline campaign assessment --------------------------------------


def test_assess_campaign_reads_score_json(tmp_path) -> None:
    eval_dir = tmp_path / "eval_000001"
    eval_dir.mkdir()
    (eval_dir / "score.json").write_text(json.dumps({
        "score": {
            "total": 16.75,
            "components": {
                "constraint_health": 0.9,
                "initial_constraint_quality": 0.8,
                "operational_ftl": 0.02,
                "ftl_persistence": 0.02,
                "horizon_penalty": 0.0,
                "survival": 1.0,
                "stability": 0.85,
                "constraint_growth": 0.9,
            },
        },
        "metrics": {},
    }))

    result = assess_campaign(tmp_path, min_tier=int(Tier.OPERATIONAL))

    assert result["num_evaluated"] == 1
    assert result["tier_histogram"]["persistent"] == 1
    assert result["num_survivors"] == 1
    assert result["survivor_front"][0]["label"] == "eval_000001"
    assert (tmp_path / "validation.json").exists()
