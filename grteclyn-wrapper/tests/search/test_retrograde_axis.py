"""A retrograde-only campaign should not search the half of omega it folds away.

``trajectory_retrograde_only`` negates any positive ``omega_rot`` during decode,
so genome ``+x`` and ``-x`` describe the same spacetime.  Left unrestricted, the
upper half of every omega axis is a mirror image the optimizer cannot tell from
signal: CMA-ES fits covariance across the fold, and QD mutations that cross zero
land back where they started.
"""

from __future__ import annotations

import math

from grteclyn_wrapper.search.optimize.candidates import _vector_to_overrides
from grteclyn_wrapper.search.optimize.dimension import SearchDimension
from grteclyn_wrapper.search.optimize.spaces import (
    grtresna_trajectory_search_space,
    restrict_retrograde_omega,
)

_OMEGA = "trajectory_lump0_omega_rot"


def _by_key(dims) -> dict[str, SearchDimension]:
    return {d.param_key: d for d in dims}


def test_the_omega_axis_is_two_sided_by_default() -> None:
    dims = _by_key(grtresna_trajectory_search_space(num_lumps=2))
    assert dims[_OMEGA].lower < 0.0 < dims[_OMEGA].upper


def test_a_retrograde_campaign_searches_only_the_retrograde_half() -> None:
    dims = grtresna_trajectory_search_space(num_lumps=2)
    restricted = _by_key(restrict_retrograde_omega(dims, enabled=True))
    for key in (_OMEGA, "trajectory_lump1_omega_rot"):
        assert restricted[key].upper == 0.0
        assert restricted[key].lower < 0.0
        assert restricted[key].center <= 0.0


def test_nothing_else_moves() -> None:
    dims = grtresna_trajectory_search_space(num_lumps=2)
    restricted = restrict_retrograde_omega(dims, enabled=True)
    for original, updated in zip(dims, restricted):
        if original.param_key.endswith("_omega_rot"):
            continue
        assert original == updated


def test_disabled_is_a_no_op() -> None:
    dims = grtresna_trajectory_search_space(num_lumps=2)
    assert restrict_retrograde_omega(dims, enabled=False) == list(dims)


def test_the_restricted_axis_never_triggers_the_fold() -> None:
    """Every genome in the narrowed box decodes without being reflected.

    That is the property the fold was destroying: two distinct genomes used to
    decode to one physical omega, and the stored value could never be positive.
    """
    dims = restrict_retrograde_omega(
        grtresna_trajectory_search_space(num_lumps=2), enabled=True
    )
    base = {"trajectory_retrograde_only": 1}
    omega_keys = [d.param_key for d in dims if d.param_key.endswith("_omega_rot")]

    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        x = [d.lower + fraction * d.range for d in dims]
        folded = _vector_to_overrides(x, dims, base)
        unfolded = _vector_to_overrides(x, dims, {})
        for key in omega_keys:
            # Identical with and without the fold => the fold did nothing.
            assert math.isclose(
                float(folded[key]), float(unfolded[key]), rel_tol=0.0, abs_tol=1e-15
            )
            assert float(folded[key]) <= 0.0
