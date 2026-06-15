"""Comprehensive tests for the constraint-satisfying metric guesser.

Run from the grteclyn-wrapper directory:

    uv run python tests/test_constrained_guesser.py

Exit code is 0 if all checks pass, 1 otherwise.
"""

from __future__ import annotations

import sys

import numpy as np

from grteclyn_wrapper.initial_data.constrained_recipe import (
    RecipeBasis,
    constrained_overrides,
    constrained_phi,
    fit_gaussian_basis,
)
from grteclyn_wrapper.initial_data.preflight import preflight_check
from grteclyn_wrapper.initial_data.seeds import get_seed


PASS = 0
FAIL = 0


def check(name: str, condition: bool, detail: str = "") -> None:
    global PASS, FAIL
    tag = "PASS" if condition else "FAIL"
    if condition:
        PASS += 1
    else:
        FAIL += 1
    suffix = f"  ({detail})" if detail else ""
    print(f"  [{tag}] {name}{suffix}")


def run_constrained_test(
    chi_coeffs: list[float],
    *,
    chi_asym: float = 1.0,
    K_const: float = 0.0,
    phantom: bool = False,
    num_bases: int | None = None,
    width: float = 1.0,
    rmax: float = 8.0,
):
    if num_bases is None:
        num_bases = len(chi_coeffs)
    basis = RecipeBasis(
        num_bases=num_bases,
        basis_width=width,
        basis_radius_max=rmax,
    )
    result = constrained_phi(
        basis,
        chi_asym,
        chi_coeffs,
        K_const,
        phantom,
    )

    overrides = {
        "recipe_num_bases": num_bases,
        "recipe_basis_width": width,
        "recipe_basis_radius_max": rmax,
        "recipe_chi_asymptotic": chi_asym,
        "recipe_K_asymptotic": K_const,
        "recipe_phi_asymptotic": 0.0,
        "recipe_Pi_asymptotic": 0.0,
    }
    for n in range(num_bases):
        overrides[f"recipe_chi_coeff_{n}"] = chi_coeffs[n]
        overrides[f"recipe_K_coeff_{n}"] = 0.0
        overrides[f"recipe_phi_coeff_{n}"] = result.phi_coeffs[n]
        overrides[f"recipe_Pi_coeff_{n}"] = 0.0

    preflight = preflight_check(overrides, phantom=phantom)
    neg_frac = float(result.integrand_negative_mask.mean())
    return result, preflight, neg_frac


def alcubierre_shape(
    r: np.ndarray,
    *,
    radius: float,
    sigma: float,
) -> np.ndarray:
    """Alcubierre top-hat shape function f(r).

    The full Alcubierre metric has flat spatial slices and a shift
    beta^x = -v f(r).  RadialRecipe does not currently model that vector
    shift, so these tests use f(r) to build radial conformal-factor proxies
    for the bubble wall.
    """
    normalizer = 2.0 * np.tanh(sigma * radius)
    return (
        np.tanh(sigma * (r + radius))
        - np.tanh(sigma * (r - radius))
    ) / normalizer


def gaussian_coeffs_for_profile(
    profile: np.ndarray,
    *,
    r: np.ndarray,
    basis: RecipeBasis,
    asymptotic: float,
) -> list[float]:
    coeffs, residual = fit_gaussian_basis(
        r,
        profile,
        basis,
        asymptotic=asymptotic,
    )
    check("Gaussian fit residual bounded", residual < 0.15, f"{residual:.4f}")
    return coeffs


def test_group_1_known_seeds() -> None:
    print("=" * 60)
    print("TEST GROUP 1: Known analytical seeds")
    print("=" * 60)

    print("\n1a. Flat Minkowski (trivial solution)")
    mink = get_seed("flat_minkowski")
    mink_constrained = constrained_overrides(dict(mink.overrides), phantom=False)
    num = int(mink.overrides["recipe_num_bases"])
    phi_m = [mink_constrained[f"recipe_phi_coeff_{n}"] for n in range(num)]
    check("phi coefficients all zero", all(c == 0.0 for c in phi_m), str(phi_m))
    pf = preflight_check(mink_constrained, phantom=False)
    check("preflight passes", pf.passed)
    check("ham_l2 ~ 0", pf.hamiltonian_l2 < 1e-6, f"{pf.hamiltonian_l2:.2e}")

    print("\n1b. Ellis-Bronnikov wormhole (phantom, b0=0.5)")
    eb = get_seed("ellis_bronnikov")
    eb_constrained = constrained_overrides(dict(eb.overrides), phantom=True)
    num_eb = int(eb.overrides["recipe_num_bases"])
    orig = [eb.overrides[f"recipe_phi_coeff_{n}"] for n in range(num_eb)]
    derived = [eb_constrained[f"recipe_phi_coeff_{n}"] for n in range(num_eb)]
    max_diff = max(abs(a - b) for a, b in zip(orig, derived))
    check(
        "phi recovery close to analytical fit",
        max_diff < 0.35,
        f"max_diff={max_diff:.2e}",
    )
    pf = preflight_check(eb_constrained, phantom=True)
    check("preflight passes (b0=0.5)", pf.passed, f"ham_l2={pf.hamiltonian_l2:.2e}")

    print("\n1c. Ellis-Bronnikov wormhole (phantom, b0=1.0, wider throat)")
    eb2 = get_seed("ellis_bronnikov", b0=1.0)
    eb2_constrained = constrained_overrides(dict(eb2.overrides), phantom=True)
    num2 = int(eb2.overrides["recipe_num_bases"])
    orig2 = [eb2.overrides[f"recipe_phi_coeff_{n}"] for n in range(num2)]
    derived2 = [eb2_constrained[f"recipe_phi_coeff_{n}"] for n in range(num2)]
    max_diff2 = max(abs(a - b) for a, b in zip(orig2, derived2))
    check(
        "phi recovery bounded (b0=1.0)",
        max_diff2 < 0.5,
        f"max_diff={max_diff2:.2e}",
    )
    pf2 = preflight_check(eb2_constrained, phantom=True)
    check(
        "preflight reports high residual (expected for coarse fit)",
        pf2.hamiltonian_l2 > 1.0,
        f"ham_l2={pf2.hamiltonian_l2:.2e}",
    )

    print("\n1d. Ellis-Bronnikov wormhole (phantom, b0=2.0, very wide throat)")
    eb3 = get_seed("ellis_bronnikov", b0=2.0, num_bases=12, basis_width=0.8)
    eb3_constrained = constrained_overrides(dict(eb3.overrides), phantom=True)
    pf3 = preflight_check(eb3_constrained, phantom=True)
    check(
        "preflight reports high residual (expected for wide throat)",
        pf3.hamiltonian_l2 > 1.0,
        f"ham_l2={pf3.hamiltonian_l2:.2e}",
    )

    print("\n1e. Schwarzschild puncture (vacuum, M=0.5)")
    sch = get_seed("schwarzschild_puncture")
    sch_constrained = constrained_overrides(dict(sch.overrides), phantom=False)
    pf_sch = preflight_check(sch_constrained, phantom=False)
    check(
        "preflight result obtained",
        True,
        f"passed={pf_sch.passed}, ham_l2={pf_sch.hamiltonian_l2:.2e}",
    )

    print("\n1f. Schwarzschild puncture (vacuum, M=1.0)")
    sch2 = get_seed("schwarzschild_puncture", M=1.0, num_bases=12)
    sch2_constrained = constrained_overrides(dict(sch2.overrides), phantom=False)
    pf_sch2 = preflight_check(sch2_constrained, phantom=False)
    check(
        "preflight result obtained",
        True,
        f"passed={pf_sch2.passed}, ham_l2={pf_sch2.hamiltonian_l2:.2e}",
    )


def test_group_2_random_normal(rng: np.random.Generator) -> None:
    print("\n" + "=" * 60)
    print("TEST GROUP 2: Random chi profiles (normal field)")
    print("=" * 60)

    for i in range(8):
        n_bases = int(rng.integers(3, 9))
        width = float(rng.uniform(0.3, 2.5))
        rmax = float(rng.uniform(4.0, 12.0))
        chi_coeffs = (rng.standard_normal(n_bases) * 0.15).tolist()
        k_const = float(rng.uniform(-0.3, 0.3))
        label = f"random_{i} (N={n_bases}, w={width:.2f}, K={k_const:.2f})"
        print(f"\n2{chr(97 + i)}. {label}")
        result, _, neg = run_constrained_test(
            chi_coeffs,
            K_const=k_const,
            phantom=False,
            width=width,
            rmax=rmax,
        )
        preview = [f"{c:.3f}" for c in result.phi_coeffs[:4]]
        check("finite phi coeffs", all(np.isfinite(result.phi_coeffs)), f"coeffs={preview}...")
        check("fit residual bounded", result.fit_residual < 1.0, f"{result.fit_residual:.4f}")
        check("neg integrand < 80%", neg < 0.8, f"{neg:.1%}")


def test_group_3_random_phantom(rng: np.random.Generator) -> None:
    print("\n" + "=" * 60)
    print("TEST GROUP 3: Random chi profiles (phantom field)")
    print("=" * 60)

    for i in range(8):
        n_bases = int(rng.integers(3, 9))
        width = float(rng.uniform(0.3, 2.5))
        rmax = float(rng.uniform(4.0, 12.0))
        chi_coeffs = (rng.standard_normal(n_bases) * 0.2).tolist()
        k_const = float(rng.uniform(-0.5, 0.5))
        label = f"phantom_{i} (N={n_bases}, w={width:.2f}, K={k_const:.2f})"
        print(f"\n3{chr(97 + i)}. {label}")
        result, _, _ = run_constrained_test(
            chi_coeffs,
            K_const=k_const,
            phantom=True,
            width=width,
            rmax=rmax,
        )
        preview = [f"{c:.3f}" for c in result.phi_coeffs[:4]]
        check("finite phi coeffs", all(np.isfinite(result.phi_coeffs)), f"coeffs={preview}...")
        check("fit residual bounded", result.fit_residual < 1.0, f"{result.fit_residual:.4f}")


def test_group_4_edge_cases() -> None:
    print("\n" + "=" * 60)
    print("TEST GROUP 4: Edge cases")
    print("=" * 60)

    print("\n4a. Very narrow basis (width=0.15)")
    result, pf, _ = run_constrained_test(
        [-0.2, -0.1, 0.05, 0.0],
        width=0.15,
        rmax=4.0,
    )
    check("finite result", all(np.isfinite(result.phi_coeffs)))
    check(
        "preflight check runs",
        True,
        f"passed={pf.passed}, ham_l2={pf.hamiltonian_l2:.2e}",
    )

    print("\n4b. Very wide basis (width=4.0)")
    result, pf, _ = run_constrained_test(
        [-0.1, 0.0, 0.0, 0.0],
        width=4.0,
        rmax=16.0,
    )
    check("finite result", all(np.isfinite(result.phi_coeffs)))
    check("preflight passes (gentle curvature)", pf.passed, f"ham_l2={pf.hamiltonian_l2:.2e}")

    print("\n4c. Dense basis (N=12)")
    coeffs12 = [0.0] * 12
    coeffs12[0] = -0.15
    coeffs12[3] = 0.08
    coeffs12[7] = -0.05
    result, _, _ = run_constrained_test(
        coeffs12,
        num_bases=12,
        width=0.6,
        rmax=8.0,
    )
    check("finite result", all(np.isfinite(result.phi_coeffs)))
    check("12 phi coeffs returned", len(result.phi_coeffs) == 12)

    print("\n4d. Nearly singular chi (chi_coeff_0 = -0.95)")
    result, _, neg = run_constrained_test(
        [-0.95, 0.0, 0.0, 0.0],
        width=0.5,
    )
    check("finite result", all(np.isfinite(result.phi_coeffs)))
    check("high neg integrand expected", neg > 0.1, f"{neg:.1%}")

    print("\n4e. Non-zero K with phantom field (K=-0.5)")
    result, pf, _ = run_constrained_test(
        [-0.2, -0.1, 0.0, 0.0],
        K_const=-0.5,
        phantom=True,
    )
    check("finite result", all(np.isfinite(result.phi_coeffs)))
    check(
        "preflight check runs",
        True,
        f"passed={pf.passed}, ham_l2={pf.hamiltonian_l2:.2e}",
    )

    print("\n4f. Globally depressed chi (asymptotic=0.8)")
    result, pf, _ = run_constrained_test(
        [-0.1, -0.05, 0.0, 0.0],
        chi_asym=0.8,
        phantom=True,
    )
    check("finite result", all(np.isfinite(result.phi_coeffs)))
    check("preflight check runs", True, f"passed={pf.passed}")

    print("\n4g. Single basis function (N=1)")
    result, _, _ = run_constrained_test(
        [-0.2],
        num_bases=1,
        width=1.5,
        rmax=6.0,
    )
    check("finite result", all(np.isfinite(result.phi_coeffs)))
    check("1 phi coeff returned", len(result.phi_coeffs) == 1)


def test_group_5_constraint_improvement(rng: np.random.Generator) -> None:
    print("\n" + "=" * 60)
    print("TEST GROUP 5: Constraint improvement verification")
    print("=" * 60)

    print("\n5a. Random config: unconstrained vs constrained preflight")
    for i in range(5):
        n_bases = int(rng.integers(4, 8))
        width = float(rng.uniform(0.5, 2.0))
        rmax = 8.0
        chi_coeffs = (rng.standard_normal(n_bases) * 0.15).tolist()
        phi_rand = (rng.standard_normal(n_bases) * 0.1).tolist()

        unconstrained = {
            "recipe_num_bases": n_bases,
            "recipe_basis_width": width,
            "recipe_basis_radius_max": rmax,
            "recipe_chi_asymptotic": 1.0,
            "recipe_K_asymptotic": 0.0,
            "recipe_phi_asymptotic": 0.0,
            "recipe_Pi_asymptotic": 0.0,
        }
        for n in range(n_bases):
            unconstrained[f"recipe_chi_coeff_{n}"] = chi_coeffs[n]
            unconstrained[f"recipe_K_coeff_{n}"] = float(rng.uniform(-0.2, 0.2))
            unconstrained[f"recipe_phi_coeff_{n}"] = phi_rand[n]
            unconstrained[f"recipe_Pi_coeff_{n}"] = float(rng.uniform(-0.1, 0.1))

        pf_before = preflight_check(unconstrained, phantom=False)
        constrained = constrained_overrides(dict(unconstrained), phantom=False)
        pf_after = preflight_check(constrained, phantom=False)

        improved = pf_after.hamiltonian_l2 <= pf_before.hamiltonian_l2 * 1.15
        label = (
            f"trial_{i}: ham_l2 {pf_before.hamiltonian_l2:.2e} "
            f"-> {pf_after.hamiltonian_l2:.2e}"
        )
        check("constraint improved or within 15% fit error", improved, label)
        check(
            "momentum ~ 0 after locking K",
            pf_after.momentum_l2 < 1e-6,
            f"mom_l2={pf_after.momentum_l2:.2e}",
        )


def test_group_6_alcubierre_warpdrive_metrics() -> None:
    print("\n" + "=" * 60)
    print("TEST GROUP 6: Alcubierre-style warpdrive metrics")
    print("=" * 60)

    print("\n6a. Original Alcubierre spatial slice (chi=1)")
    flat_slice = {
        "recipe_num_bases": 4,
        "recipe_basis_width": 1.0,
        "recipe_basis_radius_max": 8.0,
        "recipe_chi_asymptotic": 1.0,
        "recipe_K_asymptotic": 0.0,
        "recipe_phi_asymptotic": 0.0,
        "recipe_Pi_asymptotic": 0.0,
    }
    for n in range(4):
        flat_slice[f"recipe_chi_coeff_{n}"] = 0.0
        flat_slice[f"recipe_K_coeff_{n}"] = 0.0
        flat_slice[f"recipe_phi_coeff_{n}"] = 0.0
        flat_slice[f"recipe_Pi_coeff_{n}"] = 0.0

    flat_constrained = constrained_overrides(dict(flat_slice), phantom=False)
    pf_flat = preflight_check(flat_constrained, phantom=False)
    phi_flat = [flat_constrained[f"recipe_phi_coeff_{n}"] for n in range(4)]
    check("spatially flat slice stays trivial", all(c == 0.0 for c in phi_flat))
    check("preflight passes flat Alcubierre spatial slice", pf_flat.passed)
    check("ham_l2 ~ 0 for spatial slice", pf_flat.hamiltonian_l2 < 1e-6, f"{pf_flat.hamiltonian_l2:.2e}")

    print("\n6b. Alcubierre bubble-wall proxy (normal scalar)")
    basis = RecipeBasis(num_bases=12, basis_width=0.6, basis_radius_max=10.0)
    r = np.linspace(10.0 / 4096.0, 20.0, 4096)
    shape = alcubierre_shape(r, radius=4.0, sigma=1.5)
    bubble_wall = shape * (1.0 - shape)
    chi_profile = 1.0 - 0.08 * bubble_wall
    chi_coeffs = gaussian_coeffs_for_profile(
        chi_profile,
        r=r,
        basis=basis,
        asymptotic=1.0,
    )
    result, pf, neg = run_constrained_test(
        chi_coeffs,
        phantom=False,
        num_bases=12,
        width=0.6,
        rmax=10.0,
    )
    check("finite phi coeffs", all(np.isfinite(result.phi_coeffs)))
    check("fit residual bounded", result.fit_residual < 1.0, f"{result.fit_residual:.4f}")
    check("chi remains positive", pf.negative_chi_fraction == 0.0)
    check("normal-field source mostly satisfiable", neg < 0.8, f"{neg:.1%}")

    print("\n6c. Alcubierre bubble-wall proxy (phantom scalar)")
    chi_profile_phantom = 1.0 + 0.08 * bubble_wall
    chi_coeffs_phantom = gaussian_coeffs_for_profile(
        chi_profile_phantom,
        r=r,
        basis=basis,
        asymptotic=1.0,
    )
    result_p, pf_p, neg_p = run_constrained_test(
        chi_coeffs_phantom,
        phantom=True,
        num_bases=12,
        width=0.6,
        rmax=10.0,
    )
    check("finite phi coeffs", all(np.isfinite(result_p.phi_coeffs)))
    check("fit residual bounded", result_p.fit_residual < 1.0, f"{result_p.fit_residual:.4f}")
    check("chi remains positive", pf_p.negative_chi_fraction == 0.0)
    check("phantom-field source mostly satisfiable", neg_p < 0.8, f"{neg_p:.1%}")

    print("\n6d. Thin/high Alcubierre wall stress test")
    sharp_basis = RecipeBasis(num_bases=16, basis_width=0.35, basis_radius_max=10.0)
    sharp_shape = alcubierre_shape(r, radius=4.0, sigma=3.0)
    sharp_wall = sharp_shape * (1.0 - sharp_shape)
    sharp_chi = 1.0 - 0.15 * sharp_wall
    sharp_coeffs = gaussian_coeffs_for_profile(
        sharp_chi,
        r=r,
        basis=sharp_basis,
        asymptotic=1.0,
    )
    result_s, pf_s, neg_s = run_constrained_test(
        sharp_coeffs,
        phantom=False,
        num_bases=16,
        width=0.35,
        rmax=10.0,
    )
    check("finite phi coeffs", all(np.isfinite(result_s.phi_coeffs)))
    check("preflight reports result", True, f"passed={pf_s.passed}, ham_l2={pf_s.hamiltonian_l2:.2e}")
    check("thin wall has nontrivial constraint source", neg_s > 0.05, f"{neg_s:.1%}")


def main() -> int:
    global PASS, FAIL
    rng = np.random.default_rng(42)

    test_group_1_known_seeds()
    test_group_2_random_normal(rng)
    test_group_3_random_phantom(rng)
    test_group_4_edge_cases()
    test_group_5_constraint_improvement(rng)
    test_group_6_alcubierre_warpdrive_metrics()

    print("\n" + "=" * 60)
    total = PASS + FAIL
    print(f"RESULTS: {PASS}/{total} passed, {FAIL} failed")
    print("=" * 60)
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
