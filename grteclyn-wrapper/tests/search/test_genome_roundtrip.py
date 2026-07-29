"""Genome encode/decode round-trip invariance.

WHY THIS FILE EXISTS
--------------------
``test_trajectory_speed_cap.py`` already tests the *forward* transform
(genome -> physical) in eighteen different ways, and every one of them passes.
None of them tests whether the transform can be **undone**.

It cannot.  ``_clamp_trajectory_speed`` writes the physical value back into the
same dict key the genome coordinate arrived in, and the search loop then reads
that dict back as if it were a genome (``sampling._mutate_params``,
``candidates._overrides_to_vector`` for CMA-ES warm start, the near-miss pool,
``--seed-eval-dirs``, and the surrogate's training set).  Each pass multiplies
rotation by ``v_max / R0`` again -- a 9-19x contraction per round trip against
mutation noise that stays at full scale, so rotation and radial drift are not
inherited at all.  Measured on
``runs/grtresna_qd/qball_traj_bicomplex_8lump_v1``: the implied genome
coordinate drops from 0.490 (uniform random init) to 0.241 (pure noise) the
moment mutation takes over, and stays there for 175 evaluations.

See ``research/neuralspacetime/DebugPreGPU.md`` PG-1.

THE INVARIANT
-------------
Whatever a stage stores as a candidate's genome must decode to the same physics
when fed back through the decoder.  Formally, with ``D`` = decode and ``E`` =
whatever the search loop uses to recover a vector::

    D(E(D(x))) == D(x)          for every x in the search space

That is the honest statement of the requirement: the *physics* must be a fixed
point.

``E(D(x)) == x`` is NOT available for every ``x``, and it would be wrong to
demand it.  ``D`` legitimately projects: the speed disk clamp, the spiral-radius
floor and the retrograde fold each discard information on purpose, because the
configurations they discard are unphysical.  The correct strong form is stated
against that projection ``P = E . D``::

    E(D(x)) == P(x)             exactly, for every x
    P(P(x)) == P(x)             P is idempotent
    P(x)    == x                whenever no cap binds

Written that way the tests are *stricter* than a blanket identity, not looser:
they pin which coordinates the projection may move and to exactly what value, so
a silent contraction cannot hide behind "something got clamped."  ``P`` is
:func:`grteclyn_wrapper.search.optimize.candidates.project_vector`.

These tests were written to FAIL on the pre-fix code.  A fix that makes them
pass without changing the forward transform is a correct fix.
"""

from __future__ import annotations

import math

import pytest

from grteclyn_wrapper.search.optimize.candidates import (
    TRAJECTORY_V_MAX_DEFAULT,
    _overrides_to_vector,
    _vector_to_overrides,
    project_vector,
)
from grteclyn_wrapper.search.optimize.dimension import SearchDimension
from grteclyn_wrapper.search.qd_search.archive import Elite
from grteclyn_wrapper.search.qd_search.sampling import _mutate_elite, _mutate_params

# A representative slice of the qball_trajectory space: two lumps, each with the
# radius plus the two speed axes that carry the campaign's stated mechanism
# (counter-rotating lumps -> frame-dragging shear).
TRAJECTORY_DIMS = [
    SearchDimension("trajectory_lump0_R0", 1.5, 8.0, 5.0),
    SearchDimension("trajectory_lump0_omega_rot", -1.0, 1.0, 0.0),
    SearchDimension("trajectory_lump0_v_rad", -1.0, 1.0, 0.0),
    SearchDimension("trajectory_lump1_R0", 1.5, 8.0, 3.0),
    SearchDimension("trajectory_lump1_omega_rot", -1.0, 1.0, 0.0),
    SearchDimension("trajectory_lump1_v_rad", -1.0, 1.0, 0.0),
]

# Genome vectors chosen to exercise real rotation, not the degenerate x=0 case
# that hides the contraction.  (A zero genome is a fixed point of the bug.)
#
# UNCAPPED lie strictly inside every physical limit, so the decoder is a pure
# bijection on them and the identity E(D(x)) == x must hold exactly.
UNCAPPED_GENOMES = [
    pytest.param([4.0, -0.60, -0.20, 2.5, -0.45, 0.10], id="moderate-retrograde"),
    pytest.param([3.0, 0.25, 0.00, 6.0, -0.05, 0.02], id="small-rotation"),
]

# This one deliberately trips a cap: lump1 sits at R0 = 1.5 and asks for an
# inward drift of -0.30, which would carry it inside TRAJECTORY_R_MIN before
# stop_time.  The decoder projects it to the boundary.  That is the clamp doing
# its job, not the PG-1 contraction -- see test_the_projection_moves_only_capped_axes,
# which pins the resulting value exactly.
CAPPED_GENOME = [8.0, -0.95, 0.05, 1.5, 0.80, -0.30]

GENOMES = UNCAPPED_GENOMES + [pytest.param(CAPPED_GENOME, id="near-cap-mixed-sign")]

TRAJECTORY_SPEED_KEYS = [
    "trajectory_lump0_omega_rot",
    "trajectory_lump0_v_rad",
    "trajectory_lump1_omega_rot",
    "trajectory_lump1_v_rad",
]


def _physics(overrides):
    """The quantities that actually reach GRTresna, in physical units."""
    out = {}
    for k in (0, 1):
        r0 = float(overrides[f"trajectory_lump{k}_R0"])
        omega = float(overrides[f"trajectory_lump{k}_omega_rot"])
        v_rad = float(overrides[f"trajectory_lump{k}_v_rad"])
        out[f"lump{k}_v_tan"] = abs(omega) * r0
        out[f"lump{k}_v_rad"] = v_rad
    return out


# ---------------------------------------------------------------------------
# The core invariant
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("genome", UNCAPPED_GENOMES)
def test_encode_decode_recovers_the_genome(genome) -> None:
    """E(D(x)) == x, for genomes inside every physical limit.

    Failed before the fix: the decoded overrides carry physical omega, and
    ``_overrides_to_vector`` read it back with no inverse transform, so the
    recovered vector was contracted by v_max/R0 (~9-19x).
    """
    overrides = _vector_to_overrides(genome, TRAJECTORY_DIMS, {})
    recovered = _overrides_to_vector(overrides, TRAJECTORY_DIMS)

    for dim, want, got in zip(TRAJECTORY_DIMS, genome, recovered):
        assert got == pytest.approx(want, abs=1e-12), (
            f"{dim.param_key}: genome {want} was stored/recovered as {got} "
            f"(ratio {got / want if want else float('nan'):.4g}). "
            "The decoder wrote a physical value into the genome's key."
        )


@pytest.mark.parametrize("genome", GENOMES)
def test_encode_decode_lands_on_the_projection(genome) -> None:
    """E(D(x)) == P(x) exactly, and P is idempotent -- including under a cap."""
    projected = project_vector(genome, TRAJECTORY_DIMS, {})
    recovered = _overrides_to_vector(
        _vector_to_overrides(genome, TRAJECTORY_DIMS, {}), TRAJECTORY_DIMS
    )
    assert recovered == pytest.approx(projected, abs=1e-12)

    twice = project_vector(projected, TRAJECTORY_DIMS, {})
    assert twice == pytest.approx(projected, abs=1e-12), (
        "the projection is not idempotent: a stored candidate moves again every "
        "time it is reloaded, which is the PG-1 failure mode however small the step."
    )


def test_the_projection_moves_only_capped_axes() -> None:
    """Pin exactly which coordinate a cap may move, and to exactly what value.

    Without this, ``E(D(x)) == P(x)`` would be vacuous -- any contraction at all
    could hide inside P.  The one legitimate move here is the spiral-radius floor
    on lump1: v_rad_phys is held at (r_min - R0)/stop_time, and nothing else in
    the vector is allowed to budge.
    """
    projected = project_vector(CAPPED_GENOME, TRAJECTORY_DIMS, {})

    r0, r_min, stop_time = 1.5, 0.1, 16.0
    expected = ((r_min - r0) / stop_time) / TRAJECTORY_V_MAX_DEFAULT
    assert projected[5] == pytest.approx(expected, abs=1e-12)

    moved = [
        d.param_key
        for i, d in enumerate(TRAJECTORY_DIMS)
        if not math.isclose(projected[i], CAPPED_GENOME[i], abs_tol=1e-12)
    ]
    assert moved == ["trajectory_lump1_v_rad"], (
        f"the decoder moved {moved}; only the documented spiral-radius floor on "
        "lump1 may move, and only that axis."
    )


@pytest.mark.parametrize("genome", GENOMES)
def test_decode_is_idempotent_in_physics(genome) -> None:
    """D(E(D(x))) == D(x) -- the weak form.

    Even if a fix keeps genome and physical values sharing a key name, running
    a stored candidate back through the decoder must not move the physics.
    Today a second pass shrinks tangential speed by another factor v_max/R0.
    """
    once = _vector_to_overrides(genome, TRAJECTORY_DIMS, {})
    twice = _vector_to_overrides(
        _overrides_to_vector(once, TRAJECTORY_DIMS), TRAJECTORY_DIMS, {}
    )

    p1, p2 = _physics(once), _physics(twice)
    for key in p1:
        assert p2[key] == pytest.approx(p1[key], rel=1e-9, abs=1e-12), (
            f"{key}: {p1[key]:.6g} -> {p2[key]:.6g} on a second decode pass "
            "(a stored candidate is not a fixed point of the decoder)."
        )


def test_contraction_factor_is_not_silently_applied() -> None:
    """Pin the specific failure: one round trip must not scale omega by v_max/R0.

    This is the number that reached the paper -- research.tex:1101 describes
    candidate 146 as 'roughly an order of magnitude quieter than champion 87',
    and mean R0/v_max for that candidate is 13.5x.
    """
    r0 = 4.0
    genome = [r0, -0.8, 0.0, 3.0, 0.0, 0.0]
    overrides = _vector_to_overrides(genome, TRAJECTORY_DIMS, {})
    recovered = _overrides_to_vector(overrides, TRAJECTORY_DIMS)

    predicted_contraction = TRAJECTORY_V_MAX_DEFAULT / r0  # ~0.075
    omega_in, omega_out = genome[1], recovered[1]
    assert not math.isclose(
        omega_out, omega_in * predicted_contraction, rel_tol=1e-6
    ), (
        "omega was contracted by exactly v_max/R0 on one round trip -- this is "
        "the PG-1 signature."
    )


# ---------------------------------------------------------------------------
# Heritability: the search-loop consequence
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("genome", GENOMES)
def test_mutation_inherits_the_parent(genome) -> None:
    """A zero-noise mutation of a stored elite must reproduce that elite.

    ``_mutate_params`` is how ~60% of QD samples are drawn.  With sigma=0 it is
    the identity on the genome -- unless what was stored is not a genome.
    """
    import numpy as np

    stored = _vector_to_overrides(genome, TRAJECTORY_DIMS, {})
    params = {d.param_key: float(stored[d.param_key]) for d in TRAJECTORY_DIMS}

    child = _mutate_params(
        params, TRAJECTORY_DIMS, np.random.default_rng(0), sigma=0.0
    )
    want_vec = project_vector(genome, TRAJECTORY_DIMS, {})

    for dim, want, got in zip(TRAJECTORY_DIMS, want_vec, child):
        assert got == pytest.approx(want, abs=1e-12), (
            f"{dim.param_key}: noiseless mutation of an elite moved it "
            f"{want} -> {got}. Rotation is not heritable."
        )


@pytest.mark.parametrize("genome", GENOMES)
def test_a_stored_genome_is_reproduced_exactly(genome) -> None:
    """The robust half of the fix: an elite carrying its genome needs no inverse.

    ``_mutate_params`` has to reconstruct the genome from decoded physics, which
    works but depends on the transform staying invertible.  An ``Elite`` that
    recorded the optimizer vector reproduces it under a noiseless mutation with
    no reconstruction at all -- including for the capped genome, where the raw
    vector is deliberately *not* the projection.
    """
    import numpy as np

    stored = _vector_to_overrides(genome, TRAJECTORY_DIMS, {})
    elite = Elite(
        cell=(0, 0),
        score=1.0,
        descriptors=(0.0, 0.0),
        params={d.param_key: float(stored[d.param_key]) for d in TRAJECTORY_DIMS},
        genome=list(genome),
        episode=None,
    )

    child = _mutate_elite(
        elite, TRAJECTORY_DIMS, np.random.default_rng(0), sigma=0.0
    )
    assert child == pytest.approx(genome, abs=1e-12)


def test_rotation_signal_survives_mutation_noise() -> None:
    """The statistical form of the field failure, as a unit test.

    Take a strongly-rotating elite, mutate it many times at the production
    sigma, and require the children's mean rotation to retain the parent's
    signal.  Under PG-1 the parent contributes ~2% and the children are pure
    noise about zero -- exactly what trajectory.jsonl shows (0.490 -> 0.241).
    """
    import numpy as np

    rng = np.random.default_rng(1234)
    parent_genome = [4.0, -0.85, 0.0, 4.0, -0.85, 0.0]
    stored = _vector_to_overrides(parent_genome, TRAJECTORY_DIMS, {})
    params = {d.param_key: float(stored[d.param_key]) for d in TRAJECTORY_DIMS}

    omega_idx = [
        i
        for i, d in enumerate(TRAJECTORY_DIMS)
        if d.param_key.endswith("_omega_rot")
    ]
    children = [
        _mutate_params(params, TRAJECTORY_DIMS, rng, sigma=0.15) for _ in range(400)
    ]
    mean_abs = np.mean([abs(c[i]) for c in children for i in omega_idx])

    # Pure reflected noise at sigma = 0.15 * range = 0.3 gives E|x| ~ 0.239.
    # A heritable axis centred on |0.85| must sit well above that.
    assert mean_abs > 0.5, (
        f"mean |omega| across mutated children = {mean_abs:.3f}; pure noise "
        "would give ~0.239 and the parent was at 0.85. The parent's rotation "
        "is not being inherited."
    )


# ---------------------------------------------------------------------------
# Generalisation: no dimension may be lossy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("genome", UNCAPPED_GENOMES)
def test_no_search_dimension_is_lossy(genome) -> None:
    """Guard the whole space, not just the two axes we know about.

    Any future transform added to ``_vector_to_overrides`` that writes back into
    a searched key will trip this.  Restricted to genomes inside every physical
    limit so that a genuine clamp is not mistaken for a lossy transform; the
    capped case is pinned exactly by
    ``test_the_projection_moves_only_capped_axes``.
    """
    overrides = _vector_to_overrides(genome, TRAJECTORY_DIMS, {})
    recovered = _overrides_to_vector(overrides, TRAJECTORY_DIMS)

    lossy = [
        d.param_key
        for d, want, got in zip(TRAJECTORY_DIMS, genome, recovered)
        if not math.isclose(want, got, abs_tol=1e-12)
    ]
    assert not lossy, (
        f"{len(lossy)} of {len(TRAJECTORY_DIMS)} searched dimensions do not "
        f"survive a store/reload round trip: {lossy}"
    )
