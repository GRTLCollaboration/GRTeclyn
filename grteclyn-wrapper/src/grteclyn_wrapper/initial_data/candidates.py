"""Resolve RadialRecipe initial-data sources (seeds or validation candidates)."""

from __future__ import annotations

from typing import Any

from .nonspherical_guesser import NonSphericalCandidate, nonspherical_to_overrides
from .seeds import get_seed
from .validate_guesser import CandidateSpec, generate_candidates, spec_to_overrides

# GPU long-run promotion ladder (from smoke + t=5 stability checks).
PROMOTED_SURVIVORS: tuple[str, ...] = ("ellis_bronnikov", "bubble_wall_016")

# Do not schedule extended resolution / long-time GPU runs for these IDs.
REJECTED_LONG_RUN: tuple[str, ...] = ("random_000",)


def lookup_candidate(candidate_id: str, *, seed: int = 42) -> CandidateSpec:
    for spec in generate_candidates(seed):
        if spec.candidate_id == candidate_id:
            return spec
    raise KeyError(
        f"Unknown candidate_id={candidate_id!r} (validation seed={seed}). "
        "Run validate_guesser.generate_candidates() to list IDs."
    )


def lookup_nonspherical_candidate(
    candidate_id: str, *, seed: int = 42
) -> NonSphericalCandidate:
    from .nonspherical_guesser import generate_nonspherical_candidates

    for spec in generate_nonspherical_candidates(seed):
        if spec.candidate_id == candidate_id:
            return spec
    raise KeyError(
        f"Unknown non-spherical candidate_id={candidate_id!r} (seed={seed})."
    )


def candidate_requires_phantom(spec: CandidateSpec) -> bool:
    return spec.field_type == "phantom"


def resolve_initial_data_overrides(
    *,
    seed_name: str | None = None,
    candidate_id: str | None = None,
    nonspherical_id: str | None = None,
    validation_seed: int = 42,
) -> tuple[dict[str, Any], bool, str]:
    """Return (overrides, phantom, source_label) for RadialRecipe episodes."""
    if sum(x is not None for x in (seed_name, candidate_id, nonspherical_id)) > 1:
        raise ValueError(
            "Specify only one of seed_name, candidate_id, or nonspherical_id."
        )

    if nonspherical_id is not None:
        spec = lookup_nonspherical_candidate(nonspherical_id, seed=validation_seed)
        return nonspherical_to_overrides(spec), True, nonspherical_id

    if candidate_id is not None:
        spec = lookup_candidate(candidate_id, seed=validation_seed)
        return spec_to_overrides(spec), candidate_requires_phantom(spec), candidate_id

    if seed_name is not None:
        seed = get_seed(seed_name)
        phantom = seed_name == "ellis_bronnikov" or bool(
            getattr(seed, "phantom", False)
        )
        return dict(seed.overrides), phantom, seed_name

    raise ValueError("One of seed_name, candidate_id, or nonspherical_id is required.")
