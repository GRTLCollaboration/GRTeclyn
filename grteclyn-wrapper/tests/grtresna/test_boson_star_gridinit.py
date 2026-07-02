"""Tests for boson-star field painting on gridinit arrays."""

from __future__ import annotations

import numpy as np

from grteclyn_wrapper.grtresna.fields.boson_star import (
    apply_post_solve_lapse_correction,
    paint_boson_star_fields_on_grid,
)
from grteclyn_wrapper.grtresna.profiles.boson_star import solve_mini_boson_star


def test_paint_boson_star_fields_sets_phase_velocity() -> None:
    profile = solve_mini_boson_star(mass=0.1, phi_c=0.05, r_max=40.0, n_points=100)
    nz, ny, nx = 8, 8, 8
    names = ["chi", "lapse", "phi", "Pi", "phi2", "Pi2"]
    data = np.zeros((nz, ny, nx, len(names)), dtype=np.float64)
    data[:, :, :, names.index("lapse")] = 0.85
    dx = np.array([2.0, 2.0, 2.0])
    origin = np.array([-8.0, -8.0, -8.0])

    painted, out_names = paint_boson_star_fields_on_grid(
        data, names, dx, origin, profile, center=(0.0, 0.0, 0.0),
    )

    phi_idx = out_names.index("phi")
    pi2_idx = out_names.index("Pi2")
    center_val = painted[nz // 2, ny // 2, nx // 2, phi_idx]
    pi2_val = painted[nz // 2, ny // 2, nx // 2, pi2_idx]
    assert center_val > 0.0
    assert pi2_val < 0.0
    assert abs(pi2_val + profile.omega * center_val / 0.85) < 1.0e-10


def test_post_solve_lapse_correction_divides_pi2() -> None:
    """GRTresna paints Pi2 at alpha=1; the conversion must divide by lapse."""
    names = ["chi", "lapse", "phi", "Pi", "phi2", "Pi2"]
    nz, ny, nx = 4, 4, 4
    data = np.zeros((nz, ny, nx, len(names)), dtype=np.float64)
    data[:, :, :, names.index("lapse")] = 0.5
    data[:, :, :, names.index("Pi2")] = -0.2  # = -omega*phi0 at alpha=1

    corrected = apply_post_solve_lapse_correction(data, names)

    assert np.allclose(corrected[:, :, :, names.index("Pi2")], -0.4)
    # Original array is not mutated.
    assert np.allclose(data[:, :, :, names.index("Pi2")], -0.2)


def test_post_solve_lapse_correction_noop_without_columns() -> None:
    names = ["chi", "phi", "Pi"]
    data = np.ones((2, 2, 2, len(names)), dtype=np.float64)
    out = apply_post_solve_lapse_correction(data, names)
    assert np.array_equal(out, data)
