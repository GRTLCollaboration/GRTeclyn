from __future__ import annotations

import numpy as np

from grteclyn_wrapper.constrained_recipe import RecipeBasis, fit_gaussian_basis
from grteclyn_wrapper.ftl_metrics import compute_ftl_metrics
from grteclyn_wrapper.seeds import get_seed


def _alcubierre_shape(r: np.ndarray, *, radius: float, sigma: float) -> np.ndarray:
    normalizer = 2.0 * np.tanh(sigma * radius)
    return (
        np.tanh(sigma * (r + radius))
        - np.tanh(sigma * (r - radius))
    ) / normalizer


def test_flat_minkowski_has_zero_ftl() -> None:
    seed = get_seed("flat_minkowski")
    metrics = compute_ftl_metrics(seed.overrides, L=8.0)
    assert metrics.f_null == 0.0
    assert metrics.f_portal == 0.0
    assert metrics.f_throat == 0.0
    assert metrics.f_shortcut == 0.0
    assert metrics.s_nonflat < 0.05


def test_chi_compression_bump_has_portal_shortcut() -> None:
    overrides = get_seed("flat_minkowski").overrides.copy()
    overrides["recipe_chi_coeff_0"] = 0.4
    overrides["recipe_chi_coeff_1"] = 0.2
    overrides["recipe_basis_width"] = 1.0
    metrics = compute_ftl_metrics(overrides, L=8.0)
    assert metrics.s_nonflat > 0.1
    assert metrics.f_portal > 0.0
    assert metrics.f_shortcut > 0.0


def test_ellis_bronnikov_has_throat_pinch() -> None:
    seed = get_seed("ellis_bronnikov", b0=0.5, num_bases=16)
    metrics = compute_ftl_metrics(seed.overrides, L=8.0)
    assert metrics.s_nonflat > 0.1
    assert metrics.f_throat > 0.0
    assert metrics.f_shortcut > 0.0


def test_alcubierre_warp_has_null_shortcut() -> None:
    seed = get_seed("alcubierre_warp", velocity=0.5, bubble_radius=2.0)
    metrics = compute_ftl_metrics(seed.overrides, L=8.0)
    assert metrics.s_nonflat > 0.05
    assert metrics.f_null > 0.0
    assert metrics.f_shortcut > 0.0
    assert metrics.path_valid


def test_singular_lapse_invalidates_null_ray() -> None:
    overrides = get_seed("flat_minkowski").overrides.copy()
    overrides["recipe_alpha_asymptotic"] = 0.05
    overrides["recipe_alpha_coeff_0"] = -0.08
    overrides["recipe_alpha_coeff_1"] = -0.08
    metrics = compute_ftl_metrics(overrides, L=8.0)
    assert metrics.f_null == 0.0
    assert not metrics.path_valid


if __name__ == "__main__":
    test_flat_minkowski_has_zero_ftl()
    test_chi_compression_bump_has_portal_shortcut()
    test_ellis_bronnikov_has_throat_pinch()
    test_alcubierre_warp_has_null_shortcut()
    test_singular_lapse_invalidates_null_ray()
    print("ftl metrics tests passed")
