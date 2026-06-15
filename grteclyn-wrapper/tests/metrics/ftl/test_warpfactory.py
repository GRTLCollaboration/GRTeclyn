"""Tests for the Warp Factory-style energy-condition evaluator.

These reproduce the canonical results of Helmerich et al. (arXiv:2404.03095):
flat Minkowski satisfies every energy condition exactly, while the Alcubierre
warp bubble violates the NEC/WEC with violation growing with velocity.
"""

from __future__ import annotations

import numpy as np

from grteclyn_wrapper import warpfactory as wf


def test_fibonacci_sphere_unit_vectors():
    pts = wf.fibonacci_sphere(64)
    assert pts.shape == (64, 3)
    norms = np.linalg.norm(pts, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-12)


def test_flat_minkowski_satisfies_all_conditions():
    g, spacing = wf.minkowski_metric(half_width=4.0, n_space=12, dt=0.25)
    report = wf.evaluate_four_metric(g, spacing, n_directions=30, n_speeds=3)
    # Constant metric -> exactly vanishing stress-energy.
    assert abs(report.nec_min) < 1e-9
    assert abs(report.wec_min) < 1e-9
    assert report.nec_violation_fraction == 0.0
    assert report.wec_violation_fraction == 0.0
    assert report.s_energy_conditions > 0.99


def test_alcubierre_violates_energy_conditions():
    g, spacing = wf.alcubierre_metric(
        velocity=0.5, bubble_radius=2.0, sigma=2.0, half_width=4.0, n_space=22, dt=0.2
    )
    report = wf.evaluate_four_metric(g, spacing, n_directions=60, n_speeds=4, max_speed=0.9)
    # Canonical warp-drive result: NEC and WEC are violated by some observer.
    assert report.nec_min < -1e-4
    assert report.wec_min < -1e-4
    assert report.wec_violation_fraction > 0.0
    # The Eulerian (at-rest) observer already sees negative energy density.
    assert report.rho_eulerian_min < 0.0
    # Exotic matter -> the bounded energy-condition reward is small.
    assert report.s_energy_conditions < 0.1
    assert any("NEC violated" in note for note in report.notes)


def test_violation_grows_with_velocity():
    def worst_nec(v: float) -> float:
        g, spacing = wf.alcubierre_metric(
            velocity=v, bubble_radius=2.0, sigma=2.0, half_width=4.0, n_space=20, dt=0.2
        )
        return wf.evaluate_four_metric(
            g, spacing, n_directions=40, n_speeds=3, max_speed=0.9
        ).nec_min

    slow = worst_nec(0.3)
    fast = worst_nec(0.9)
    # A faster bubble requires a deeper energy-condition violation.
    assert fast < slow < 0.0


def test_einstein_tensor_vanishes_for_flat_metric():
    g, spacing = wf.minkowski_metric(half_width=3.0, n_space=10, dt=0.25)
    G = wf.einstein_tensor(g, spacing)
    assert np.allclose(G, 0.0, atol=1e-10)


def test_flat_metric_is_type_I_vacuum():
    g, spacing = wf.minkowski_metric(half_width=4.0, n_space=12, dt=0.25)
    report = wf.evaluate_four_metric(g, spacing, n_directions=20, n_speeds=2)
    # Vacuum is Type I everywhere with all slacks exactly zero (on the boundary).
    assert report.type_I_fraction == 1.0
    assert report.nec_slack_min == 0.0
    assert report.s_energy_conditions > 0.99


def test_hawking_ellis_eigenvalue_test_flags_alcubierre_violation():
    # At sufficient velocity the bubble wall is cleanly Type I and the exact
    # eigenvalue NEC/WEC slacks are negative -- an observer-independent verdict.
    g, spacing = wf.alcubierre_metric(velocity=0.9, half_width=4.0, n_space=28, dt=0.2)
    report = wf.evaluate_four_metric(g, spacing, n_directions=40, n_speeds=3)
    assert report.nec_slack_min is not None and report.nec_slack_min < 0.0
    assert report.wec_slack_min is not None and report.wec_slack_min < 0.0
    # The invariant (frame-independent) energy density goes negative.
    assert report.rho_invariant_min is not None and report.rho_invariant_min < 0.0
    assert any("exactly" in note for note in report.notes)


def test_hybrid_margin_never_certifies_clean_when_observers_see_violation():
    # Even when the eigenvalue test only catches clean bulk points, the observer
    # bound must keep the score from falsely reading "all conditions satisfied".
    g, spacing = wf.alcubierre_metric(velocity=0.3, half_width=4.0, n_space=28, dt=0.2)
    report = wf.evaluate_four_metric(g, spacing, n_directions=40, n_speeds=3)
    assert report.nec_min < 0.0
    assert report.s_energy_conditions < 0.5


def test_convergence_order_is_better_than_second_order():
    result = wf.convergence_order(velocity=0.5, resolutions=(20, 28, 40, 56))
    # Fourth-order stencils should push the observed order well above 2.
    assert result["order_estimate"] > 2.0
    assert len(result["functional"]) == 4
