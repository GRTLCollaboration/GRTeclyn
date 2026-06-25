from __future__ import annotations

import math

from .types import ScoringContext

# Below this peak matter energy density there is effectively no structure to
# persist, so the structural-persistence ratio is left undefined (survival is
# not persistence-gated) rather than dividing by numerical noise.
_RHO_PERSISTENCE_FLOOR: float = 1.0e-8


def compute_survival_components(ctx: ScoringContext) -> None:
    metrics = ctx.metrics
    components = ctx.components
    notes = ctx.notes
    final_time = ctx.final_time

    # Numerical survival: did the integrator reach the configured stop time
    # without crashing?  This is a *necessary* condition but NOT sufficient --
    # empty/dissipated space is the easiest thing to march to the end, so on its
    # own this perversely rewards junk (a candidate whose matter completely
    # disperses still scores 1.0).  We therefore gate it by structural
    # persistence below.
    if ctx.target_stop_time and ctx.target_stop_time > 0 and final_time is not None:
        numerical_survival = min(final_time / ctx.target_stop_time, 1.0)
    elif final_time is not None:
        numerical_survival = 1.0
    else:
        numerical_survival = 0.0
        notes.append("no time-series diagnostics were found")

    # Structural persistence has two independent failure modes, both of which
    # must cost the candidate its "survived" credit:
    #
    #   1. DENSITY RETENTION -- fraction of the peak matter energy density still
    #      present at the final step.  ``max_rho_required`` is the peak over the
    #      whole run, ``final_peak_rho_required`` is the value at the stop time.
    #      A configuration that dissipates/disperses sees its peak rho collapse
    #      toward 0; one that holds (or concentrates) keeps this near 1.0.
    #
    #   2. MORPHOLOGICAL COHERENCE -- whether the matter that *is* still dense
    #      remains a single coherent structure or has fragmented into several
    #      lobes.  Density retention alone is blind to this: a coherent bubble
    #      that shatters into two equally-dense counter-rotating lobes keeps its
    #      peak density (retention ~1.0) yet is no longer the localized warp
    #      structure we are searching for.  ``structure_coherence`` (computed on
    #      the evolved matter slice) is ~1/k for a structure split into k
    #      comparable pieces, so it gates the survived-credit accordingly.
    #
    # Persistence is the product of the two and defaults to 1.0 when the
    # underlying series/slice is unavailable, leaving survival un-gated.
    density_retention = 1.0
    if metrics.constraints is not None:
        peak_rho = metrics.constraints.max_rho_required
        final_rho = metrics.constraints.final_peak_rho_required
        if peak_rho is not None and final_rho is not None and peak_rho > _RHO_PERSISTENCE_FLOOR:
            ratio = final_rho / peak_rho
            if not math.isfinite(ratio):
                ratio = 0.0
            density_retention = float(min(max(ratio, 0.0), 1.0))
        else:
            notes.append("matter-density time series unavailable; survival not density-gated")

    coherence = 1.0
    if (
        metrics.general_ftl_evolved is not None
        and metrics.general_ftl_evolved.structure_coherence is not None
    ):
        raw_coherence = metrics.general_ftl_evolved.structure_coherence
        if not math.isfinite(raw_coherence):
            raw_coherence = 0.0
            notes.append("structure_coherence was NaN/inf; treating as fully fragmented")
        coherence = float(min(max(raw_coherence, 0.0), 1.0))
        if coherence < 0.95:
            notes.append(
                f"matter fragmented into ~{round(1.0 / max(coherence, 1e-6))} lobes "
                f"(coherence={coherence:.2f}); survival/shaping rewards gated down"
            )

    structural_persistence = density_retention * coherence

    components["numerical_survival"] = numerical_survival
    components["structural_persistence"] = structural_persistence
    components["survival"] = numerical_survival * structural_persistence
    if structural_persistence < 0.5:
        notes.append(
            f"matter structure dissipated (retained {structural_persistence:.0%} of peak energy density)"
        )
