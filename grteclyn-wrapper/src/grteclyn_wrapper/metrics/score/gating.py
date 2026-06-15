from __future__ import annotations

from .types import ScoringContext


def compute_nontriviality_gate(ctx: ScoringContext) -> float:
    components = ctx.components
    notes = ctx.notes

    # Non-triviality gate: a trivial flat spacetime aces every "health" reward
    # (survival, stability, clean constraints, no exotic energy) while scoring
    # zero on all exoticity terms, which previously let it out-rank genuine warp
    # / wormhole geometries (the flat-space attractor).  We gate the health
    # rewards by how non-trivial the geometry actually is, so they only count
    # once there is real structure to be healthy about.  The exoticity/FTL terms
    # are never gated -- they supply the gradient out of flatness.
    # Note: exotic energy is deliberately NOT a non-triviality signal anymore --
    # we do not want negative-energy content to "earn" the health rewards.  FTL
    # and curvature structure supply the gradient out of flatness instead.
    nontriviality = max(
        components.get("nonflat_geometry", 0.0),
        components.get("expansion_asymmetry", 0.0),
        components.get("nontrivial_geometry", 0.0),
        components.get("curvature_activity", 0.0),
        components.get("shift_drive", 0.0),
        components.get("channel_progress", 0.0),
        components.get("operational_ftl", 0.0),
        components.get("operational_ftl_geodesic", 0.0),
        components.get("ftl_persistence", 0.0),
        components.get("ftl_precursor", 0.0),
        components.get("operational_ftl_solved", 0.0),
        components.get("ftl_shortcut", 0.0),
    )
    nontriviality = float(min(max(nontriviality, 0.0), 1.0))
    components["nontriviality_gate"] = nontriviality
    if nontriviality < 0.05:
        notes.append(
            "near-trivial geometry: health rewards gated out (flat-space attractor guard)"
        )
    return nontriviality
