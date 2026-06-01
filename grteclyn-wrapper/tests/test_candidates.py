"""Tests for candidate / seed resolution."""

from __future__ import annotations

import pytest

from grteclyn_wrapper.initial_data.candidates import (
    PROMOTED_SURVIVORS,
    REJECTED_LONG_RUN,
    resolve_initial_data_overrides,
)
from grteclyn_wrapper.initial_data.nonspherical_guesser import (
    generate_nonspherical_candidates,
    nonspherical_to_overrides,
)
from grteclyn_wrapper.initial_data.validate_guesser import lookup_candidate, spec_to_overrides


def test_lookup_bubble_wall_016():
    spec = lookup_candidate("bubble_wall_016", seed=42)
    assert spec.field_type == "phantom"
    overrides = spec_to_overrides(spec)
    assert overrides["recipe_num_bases"] == 12


def test_resolve_seed():
    overrides, phantom, source = resolve_initial_data_overrides(
        seed_name="ellis_bronnikov"
    )
    assert source == "ellis_bronnikov"
    assert phantom is True
    assert "recipe_chi_coeff_0" in overrides


def test_resolve_candidate():
    overrides, phantom, source = resolve_initial_data_overrides(
        candidate_id="random_000"
    )
    assert source == "random_000"
    assert phantom is True


def test_nonspherical_overrides():
    spec = generate_nonspherical_candidates(42)[1]
    overrides = nonspherical_to_overrides(spec)
    assert overrides["recipe_num_chi_angular_modes"] == len(spec.modes)
    assert overrides["recipe_chi_mode_ell_0"] == spec.modes[0].ell


def test_promotion_lists():
    assert "ellis_bronnikov" in PROMOTED_SURVIVORS
    assert "random_000" in REJECTED_LONG_RUN
