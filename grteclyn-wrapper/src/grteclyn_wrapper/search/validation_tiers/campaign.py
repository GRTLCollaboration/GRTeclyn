"""Offline campaign assessment from per-eval score files."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .evaluate import evaluate_tiers, tier_histogram
from .survivors import build_survivors, survivor_front
from .types import DEFAULT_TIER_CONFIG, Tier, TierAssessment, TierConfig, TIER_NAMES


def assess_campaign(
    campaign_dir: Path,
    *,
    min_tier: int = int(Tier.OPERATIONAL),
    config: TierConfig = DEFAULT_TIER_CONFIG,
    write: bool = True,
) -> dict[str, Any]:
    """Assess an existing campaign offline from its ``eval_*/score.json`` files."""
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
