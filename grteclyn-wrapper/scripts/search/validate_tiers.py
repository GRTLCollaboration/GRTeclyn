#!/usr/bin/env python3
"""Assess a search campaign against the falsification-tier ladder.

Reads every ``eval_*/score.json`` in a campaign directory, climbs each
candidate up the tier ladder (CONSTRUCTED -> NONTRIVIAL -> OPERATIONAL ->
PERSISTENT -> OBSERVER_EC -> CONVERGED -> ANALYTIC), and prints the tier
histogram plus the Pareto front of survivors -- the candidate solutions. Writes
``validation.json`` into the campaign directory.

Examples
--------
    uv run python scripts/search/validate_tiers.py runs/grtresna_search/optimize_20260602T181607Z
    uv run python scripts/search/validate_tiers.py runs/qd/qd_2026... --min-tier 3
"""

from __future__ import annotations

import argparse
from pathlib import Path

from grteclyn_wrapper.search.validation_tiers import (
    DEFAULT_TIER_CONFIG,
    Tier,
    TierConfig,
    assess_campaign,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("campaign_dir", type=Path, help="campaign directory containing eval_*/score.json")
    parser.add_argument(
        "--min-tier", type=int, default=int(Tier.OPERATIONAL),
        help="lowest tier counted as a survivor for the Pareto front (default: 2 = OPERATIONAL)",
    )
    parser.add_argument("--operational-floor", type=float, default=DEFAULT_TIER_CONFIG.operational_floor)
    parser.add_argument("--stability-floor", type=float, default=DEFAULT_TIER_CONFIG.stability_floor)
    parser.add_argument("--no-write", action="store_true", help="do not write validation.json")
    args = parser.parse_args()

    config = TierConfig(
        operational_floor=args.operational_floor,
        stability_floor=args.stability_floor,
    )
    result = assess_campaign(
        args.campaign_dir,
        min_tier=args.min_tier,
        config=config,
        write=not args.no_write,
    )

    print(f"campaign       : {result['campaign']}")
    print(f"evaluated      : {result['num_evaluated']}")
    print("tier histogram :")
    for name, count in result["tier_histogram"].items():
        if count:
            print(f"  {name:<12} {count}")
    print(f"survivors (T>={args.min_tier}): {result['num_survivors']}")
    print(f"survivor front : {len(result['survivor_front'])} non-dominated candidate(s)")
    for s in result["survivor_front"]:
        obj = s["objectives"]
        print(
            f"  {s['label']:<16} tier={s['tier_name']:<11} score={s['score']:.3f} "
            f"F_op={obj.get('operational_ftl', 0.0):.4g} "
            f"anec={obj.get('anec_condition', 0.0):.3g} "
            f"stab={obj.get('constraint_growth', 0.0):.3g}"
        )
    if not result["survivor_front"]:
        print("  (none yet reached the survivor tier -- promote top candidates and re-run)")


if __name__ == "__main__":
    main()
