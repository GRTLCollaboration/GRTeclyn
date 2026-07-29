"""Three small defects that each cost the search something quietly.

PG-5: worker RNG streams were spawned from a count that disagreed with the
number of workers, so the surplus workers wrapped round and re-drew streams
already in use -- two threads mutating with identical noise.

PG-6: the exotic-injection knob matched only ``grtresna_lump{k}_*``, so for
trajectory campaigns it selected nothing while the logs reported injections.

PG-8: a NaN descriptor raised out of the archive binning, and the exception
escaped the completion callback before it could refill the pipeline -- one
in-flight slot lost per occurrence, permanently.
"""

from __future__ import annotations

import math
import random
from pathlib import Path

import numpy as np

import grteclyn_wrapper.search.qd_search.driver
from grteclyn_wrapper.search.optimize.candidates import _force_exotic_template
from grteclyn_wrapper.search.optimize.dimension import SearchDimension
from grteclyn_wrapper.search.qd_search.descriptors import _bin_index


# --- PG-5 -----------------------------------------------------------------


def test_every_worker_gets_its_own_random_stream() -> None:
    """More workers than the old fixed pool: none may share a stream."""
    sequence = np.random.SeedSequence(7)
    draws = []
    for _ in range(40):
        (child,) = sequence.spawn(1)
        draws.append(np.random.default_rng(child).normal(size=8).tobytes())
    assert len(set(draws)) == len(draws)


def test_on_demand_spawning_matches_the_old_bulk_spawn() -> None:
    """The streams existing workers get must not shift under them."""
    bulk = [
        np.random.default_rng(s).normal(size=4).tolist()
        for s in np.random.SeedSequence(7).spawn(5)
    ]
    sequence = np.random.SeedSequence(7)
    on_demand = []
    for _ in range(5):
        (child,) = sequence.spawn(1)
        on_demand.append(np.random.default_rng(child).normal(size=4).tolist())
    assert on_demand == bulk


# --- PG-6 -----------------------------------------------------------------


def _trajectory_dims(n: int) -> list[SearchDimension]:
    dims: list[SearchDimension] = []
    for k in range(n):
        dims.append(SearchDimension(f"trajectory_lump{k}_R0", 1.5, 8.0, 3.0))
        dims.append(SearchDimension(f"trajectory_lump{k}_omega_rot", -1.0, 1.0, 0.0))
        dims.append(SearchDimension(f"trajectory_lump{k}_exotic", 0.0, 1.0, 0.0))
    return dims


def test_the_exotic_template_reaches_trajectory_lumps() -> None:
    dims = _trajectory_dims(4)
    x = [d.center for d in dims]
    rng = random.Random(0)

    exotic_at = {
        d.param_key: i for i, d in enumerate(dims) if d.param_key.endswith("_exotic")
    }

    flagged_per_template = []
    for template_index in range(4):
        vec = _force_exotic_template(x, dims, rng, template_index)
        flagged = {
            k
            for k in range(4)
            if vec[exotic_at[f"trajectory_lump{k}_exotic"]] >= 0.5
        }
        flagged_per_template.append(flagged)

    # Every template must flag at least one lump -- the whole failure was that
    # none of them flagged any.
    assert all(flagged for flagged in flagged_per_template)
    # And the four templates must not all be the same pattern.
    assert len({frozenset(f) for f in flagged_per_template}) > 1


def test_the_exotic_template_leaves_the_orbit_alone() -> None:
    """The orbit is searched; a template that rewrote it would erase the search."""
    dims = _trajectory_dims(3)
    x = [d.center for d in dims]
    vec = _force_exotic_template(x, dims, random.Random(0), 1)
    for i, dim in enumerate(dims):
        if not dim.param_key.endswith("_exotic"):
            assert vec[i] == x[i]


# --- PG-8 -----------------------------------------------------------------


def test_a_nan_descriptor_does_not_raise() -> None:
    assert _bin_index(float("nan"), 8) == 0


def test_infinite_descriptors_land_at_the_ends() -> None:
    assert _bin_index(float("inf"), 8) == 7
    assert _bin_index(float("-inf"), 8) == 0


def test_ordinary_descriptors_bin_exactly_as_before() -> None:
    bins = 8
    for value in (0.0, 0.01, 0.125, 0.5, 0.99, 1.0, -0.3, 1.7):
        expected = int(min(bins - 1, max(0, math.floor(value * bins))))
        assert _bin_index(value, bins) == expected


def test_the_pipeline_refill_survives_a_failed_recording() -> None:
    """The refill must sit in a ``finally``, not after the recording.

    The callback runs inside a live search with real archives and GPU leases,
    so this pins the structure rather than the behaviour: whatever the recording
    does, the submit that keeps the pipeline full still happens.
    """
    source = (
        Path(grteclyn_wrapper.search.qd_search.driver.__file__)
    ).read_text(encoding="utf-8")
    body = source.split("def _on_eval_complete", 1)[1].split("\n    def ", 1)[0]
    assert "try:" in body
    assert "finally:" in body
    finally_block = body.split("finally:", 1)[1]
    assert "pipeline.submit(new_x)" in finally_block
    # And nowhere else -- a second, unprotected refill would defeat the point.
    assert source.count("pipeline.submit(new_x)") == 1
