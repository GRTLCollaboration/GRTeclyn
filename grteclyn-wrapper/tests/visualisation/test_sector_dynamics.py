"""Tests for the Bondi scrutiny stream (core position, momentum, gauge).

These cover the pure numerics; the yt/plotfile plumbing is exercised by the
consumer's own integration path.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from grteclyn_wrapper.visualisation.process_wave.consume_plotfiles.extraction.sector_dynamics import (
    SECTOR_DYNAMICS_COLUMNS,
    _NAN_CORE,
    _core_moments,
    _momentum_x,
    _proper_separation,
)

N = 48
DX = 0.5
AXES = [(np.arange(N) + 0.5) * DX for _ in range(3)]


def _grid():
    return np.meshgrid(*AXES, indexing="ij")


def _blob(x, y, z, *, center, width=1.0, amp=1.0):
    r2 = (x - center[0]) ** 2 + (y - center[1]) ** 2 + (z - center[2]) ** 2
    return amp * np.exp(-r2 / (2.0 * width * width))


def test_core_tracker_ignores_the_halo() -> None:
    """The point of the stream: a big shed halo must not drag the core.

    The reference pair run loses ~35% of its tracked activity to a halo that
    doubles the rms radius by t=60, which is exactly what makes the plain
    barycentre unreliable late.  A core-weighted centroid must stay on the star.
    """
    x, y, z = _grid()
    # On a cell centre, so the sampled peak is the continuum peak.
    centre = (8.25, 12.25, 12.25)
    core = _blob(x, y, z, center=centre, width=1.0, amp=1.0)
    # A broad, off-centre halo carrying MORE total weight than the core.
    halo = _blob(x, y, z, center=(18.25, 12.25, 12.25), width=6.0, amp=0.12)
    assert halo.sum() > core.sum()

    plain_barycentre_x = float((x * (core + halo)).sum() / (core + halo).sum())
    tracked = _core_moments(core + halo, x, y, z)

    # The plain first moment is dragged well off the star; the core tracker is not.
    assert plain_barycentre_x > centre[0] + 2.0
    assert tracked["x"] == pytest.approx(centre[0], abs=0.2)
    assert tracked["y"] == pytest.approx(centre[1], abs=0.2)
    assert tracked["peak"] == pytest.approx(1.0, rel=0.05)


def test_absent_sector_reads_as_absent() -> None:
    """An empty sector must be nan, never a core sitting at the origin."""
    x, y, z = _grid()
    out = _core_moments(np.zeros_like(x), x, y, z)
    assert math.isnan(out["x"]) and math.isnan(out["peak"])
    assert out["ix"] is None


def test_momentum_sign_follows_the_motion() -> None:
    """A rightward-translating field carries positive x momentum.

    With the evolution's convention (d_t phi = alpha * Pi, so Pi = d_t phi for
    alpha = 1) a profile f(x - v t) has Pi = -v f'(x), and the momentum density
    S_x = -(Pi d_x phi) = +v (f')^2 > 0.
    """
    x, y, z = _grid()
    phi1 = _blob(x, y, z, center=(12.0, 12.0, 12.0), width=1.5)
    dphi_dx = np.gradient(phi1, DX, axis=0)
    zero = np.zeros_like(phi1)
    sqrt_gamma = np.ones_like(phi1)

    v = 0.1
    px_right = _momentum_x(phi1, -v * dphi_dx, zero, zero, sqrt_gamma, DX)
    px_left = _momentum_x(phi1, +v * dphi_dx, zero, zero, sqrt_gamma, DX)

    assert px_right > 0.0
    assert px_left == pytest.approx(-px_right, rel=1e-9)
    # Momentum scales linearly with the velocity.
    px_double = _momentum_x(phi1, -2.0 * v * dphi_dx, zero, zero, sqrt_gamma, DX)
    assert px_double == pytest.approx(2.0 * px_right, rel=1e-9)


def test_momentum_at_rest_is_zero() -> None:
    """A stationary real profile carries no momentum (the t=0 pair state)."""
    x, y, z = _grid()
    phi1 = _blob(x, y, z, center=(12.0, 12.0, 12.0), width=1.5)
    zero = np.zeros_like(phi1)
    px = _momentum_x(phi1, zero, zero, zero, np.ones_like(phi1), DX)
    assert px == pytest.approx(0.0, abs=1e-12)


def test_momentum_uses_the_proper_volume_element() -> None:
    """sqrt(gamma) weighting must actually be applied."""
    x, y, z = _grid()
    phi1 = _blob(x, y, z, center=(12.0, 12.0, 12.0), width=1.5)
    dphi_dx = np.gradient(phi1, DX, axis=0)
    pi1 = -0.1 * dphi_dx
    zero = np.zeros_like(phi1)
    flat = _momentum_x(phi1, pi1, zero, zero, np.ones_like(phi1), DX)
    curved = _momentum_x(phi1, pi1, zero, zero, 2.0 * np.ones_like(phi1), DX)
    assert curved == pytest.approx(2.0 * flat, rel=1e-9)


def test_proper_separation_reads_the_metric() -> None:
    """gamma_xx = h11/chi: proper length must scale as sqrt(gamma_xx)."""
    shape = (N, N, N)
    chi = np.ones(shape)
    core_a = {"x": 15.0, "y": 6.25, "z": 6.25}
    core_b = {"x": 5.0, "y": 6.25, "z": 6.25}

    flat = _proper_separation(np.ones(shape), chi, core_a, core_b, AXES)
    assert flat == pytest.approx(10.0, rel=1e-12)

    # gamma_xx = 4 everywhere -> proper distance doubles.
    stretched = _proper_separation(4.0 * np.ones(shape), chi, core_a, core_b, AXES)
    assert stretched == pytest.approx(2.0 * flat, rel=1e-12)

    # chi scales the other way: gamma_xx = h11/chi.
    squashed = _proper_separation(np.ones(shape), 4.0 * chi, core_a, core_b, AXES)
    assert squashed == pytest.approx(0.5 * flat, rel=1e-12)

    # A missing core cannot produce a separation.
    assert math.isnan(
        _proper_separation(np.ones(shape), chi, core_a, dict(_NAN_CORE), AXES)
    )


def test_proper_separation_resolves_below_one_cell() -> None:
    """The whole point of the 2026-08-21 repair.

    The old integrator ran between the two cores' integer PEAK CELL indices, so
    every answer was a multiple of dx and a pair holding its separation to a
    tenth of a cell read as perfectly constant -- or jumped a whole cell.  Both
    endpoints are sub-cell now, and flat space must return the coordinate
    separation exactly, whatever fraction of a cell it lands on.
    """
    shape = (N, N, N)
    chi = np.ones(shape)
    h11 = np.ones(shape)
    for xa, xb in ((15.0043, 5.0044), (12.26, 11.74), (15.0, 5.0)):
        core_a = {"x": xa, "y": 6.25, "z": 6.25}
        core_b = {"x": xb, "y": 6.25, "z": 6.25}
        got = _proper_separation(h11, chi, core_a, core_b, AXES)
        assert got == pytest.approx(abs(xa - xb), abs=1e-12), f"{xa} - {xb}"


def test_proper_separation_interpolates_across_the_row() -> None:
    """Transverse position is sub-cell too, so the row must be interpolated."""
    shape = (N, N, N)
    _, y, _ = _grid()
    h11 = 1.0 + 0.01 * (y - 6.0)          # varies only with y
    core_a = {"x": 15.0, "y": 6.25, "z": 6.25}   # 6.25 is a cell EDGE: worst case
    core_b = {"x": 5.0, "y": 6.25, "z": 6.25}
    got = _proper_separation(h11, np.ones(shape), core_a, core_b, AXES)
    assert got == pytest.approx(10.0 * math.sqrt(1.0 + 0.01 * 0.25), rel=1e-12)


def test_proper_separation_refuses_to_run_off_the_grid() -> None:
    """Clamping at the boundary would silently return a truncated length."""
    shape = (N, N, N)
    core_a = {"x": 1.0e3, "y": 6.25, "z": 6.25}
    core_b = {"x": 5.0, "y": 6.25, "z": 6.25}
    assert math.isnan(
        _proper_separation(np.ones(shape), np.ones(shape), core_a, core_b, AXES)
    )


def test_column_contract_is_stable() -> None:
    """The stream is positional; guard the order the analysis scripts read."""
    assert SECTOR_DYNAMICS_COLUMNS[0] == "time"
    assert SECTOR_DYNAMICS_COLUMNS.index("core_x_canon") == 1
    assert SECTOR_DYNAMICS_COLUMNS.index("core_x_phantom") == 5
    assert SECTOR_DYNAMICS_COLUMNS.index("coord_sep") == 9
    assert SECTOR_DYNAMICS_COLUMNS.index("proper_sep") == 10
    assert SECTOR_DYNAMICS_COLUMNS.index("px_total") == 13
    assert len(SECTOR_DYNAMICS_COLUMNS) == 17
