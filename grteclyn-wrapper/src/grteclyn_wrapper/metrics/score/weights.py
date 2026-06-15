from __future__ import annotations

DEFAULT_WEIGHTS: dict[str, float] = {
    # t=0 FTL signals are now only faint "looked promising at t=0" hints -- they
    # are downweighted hard because a t=0 shortcut that does not survive the
    # evolution is a gauge/initial-data artifact (see operational_ftl, which
    # rewards the *evolved* advantage instead).
    "ftl_shortcut": 1.0,
    "expansion_asymmetry": 0.75,
    "nonflat_geometry": 1.0,
    "comoving_stability": 2.5,
    "survival": 1.5,
    "constraint_health": 2.0,
    "lapse_health": 1.0,
    "horizon_penalty": 1.5,
    "energy_condition": 2.0,
    "stability": 0.5,
    "instability_penalty": 8.0,
    "nontrivial_geometry": 0.25,
    "initial_constraint_quality": 0.5,
    # Proposed extensions (Sec. "Proposed Extensions").
    "constraint_growth": 2.0,
    "anec_condition": 1.5,
    "tidal_comfort": 1.0,
    # Mechanism-agnostic operational FTL, measured on the EVOLVED spacetime
    # (Dijkstra arrival-time advantage that survives the dynamics) -- this is the
    # primary goal and dominates the reward budget.  A t=0-only shortcut scores
    # nothing here.  Curvature activity rewards genuinely non-trivial geometry
    # (and keeps flat space out of the running).
    "operational_ftl": 3.0,
    "operational_ftl_geodesic": 12.0,
    # Continuous FTL *precursor* (shaping gradient).  operational_ftl is a hard
    # gate -- it stays exactly 0 until a connected superluminal channel beats
    # the flat baseline end-to-end -- so an "almost there" geometry (light cones
    # tilted, a few locally superluminal cells, but no through-channel yet) gets
    # no signal and the optimizer cannot climb toward FTL.  The precursor rewards
    # the cone-tilting itself (max local coordinate light speed past 1 + the
    # superluminal-cell fraction), which rises smoothly *before* the gate fires.
    # Flat space has max_local_speed == 1 exactly, so the precursor is 0 in
    # flatness: it is a clean slope out of the flat-space basin toward FTL.
    "ftl_precursor": 3.0,
    # Direct frame-dragging / light-cone-tilt drive.  This gives the optimizer a
    # signal before the end-to-end FTL path exists, but unlike stability terms it
    # points at the mechanism we actually need.
    "shift_drive": 2.0,
    # Coupled path-progress shaping: rewards evolved geometries whose fastest
    # causal path is close to beating the flat baseline *and* carry both cone
    # opening (precursor) and frame-drag (shift) together -- not either alone.
    "channel_progress": 2.5,
    # Constraint-satisfying GRTresna geometry at t=0 (mechanism-agnostic Dijkstra).
    # Shaping gradient before GPU evolution; evolved operational_ftl remains primary.
    "operational_ftl_solved": 4.0,
    "curvature_activity": 0.5,
    # Exotic-matter PENALTY.  The objective is FTL *without* exotic matter, so
    # any negative-energy requirement is punished -- both the matter sourced at
    # t=0 by the constrained solve and the geometric exotic energy of the
    # evolved spacetime (effective NEC violation, T^eff = G/8pi).  The component
    # is negative; this weight sets how hard exoticity is penalized.
    "exotic_penalty": 8.0,
    "qei_penalty": 4.0,
    "boundary_penalty": 6.0,
    "transport_objective": 2.0,
}


# "Health/niceness" rewards that a trivial flat spacetime maximises by default.
# They are gated by a non-triviality factor (see ``score_episode``) so flat
# Minkowski cannot bank them and out-rank a genuinely exotic geometry.
HEALTH_COMPONENTS: frozenset[str] = frozenset({
    "survival",
    "constraint_health",
    "energy_condition",
    "initial_constraint_quality",
    "lapse_health",
    "stability",
    "comoving_stability",
    "constraint_growth",
    "anec_condition",
    "tidal_comfort",
})
