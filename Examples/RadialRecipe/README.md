# RadialRecipe

A GRTeclyn example that builds initial data from a parameterised "recipe"
(radial basis functions + optional angular/`Y_lm` modes) for `chi`, lapse,
shift, `K`, and the scalar field, then evolves it with CCZ4 + matter. It is the
simulation backend for the closed-loop spacetime search in `grteclyn-wrapper`.

## Initial-data recipe

`RadialRecipeInitialData` superposes Gaussian-shell radial bases on each field
and adds optional non-spherical structure:

- `recipe_<field>_coeff_n`, `recipe_basis_width`, `recipe_basis_radius_max` —
  radial profile of `chi`, `alpha` (lapse), `beta` (shift), `K`, `phi`, `Pi`.
- `recipe_num_chi_Ylm_modes` + `recipe_chi_Ylm_*` — spherical-harmonic deformation
  of `chi`.
- `recipe_num_{chi,lapse,beta,K}_angular_modes` + `recipe_<field>_mode_*` —
  Legendre `P_l(cos theta)` angular modes on each field, so the search can build
  general non-spherical lapse/shift/K geometries (not just axial warp bubbles).

## In-situ diagnostics

Enabled in `params.txt` and written under `<output_path>/<data_subpath>/`:

| Flag | Output | Contents |
|------|--------|----------|
| `calculate_constraint_norms` | `constraint_norms.dat` | `L2_Ham`, `L2_Mom`, `min/max_rho_req`, `integral_neg_rho` |
| `calculate_energy_conditions` | `energy_conditions.dat` | observer-sampled `matter_min_{NEC,WEC,SEC,DEC}`, `matter_integral_NEC_violation` |
| `calculate_curvature_invariants` | `curvature_invariants.dat` | `max_abs_ricci_scalar`, `max_ricci_tensor_sq` (`R_ij R^ij`), `max_Kij_sq` (`K_ij K^ij`), `L2_ricci_scalar` |

Collapse/horizon diagnostics (`collapse_diagnostics.dat`) are always written on
level 0. The energy conditions use a Warp-Factory-style observer sampling of the
3+1 stress-energy `(rho, j_i, S_ij)` reconstructed inline via
`compute_emtensor`; the rigorous Hawking–Ellis eigenvalue refinement and the
geometry-sourced effective stress energy (`T^eff = G / 8pi`) are evaluated
post-hoc from plotfiles by the Python `warpfactory.py` module.

## Matter sector: canonical vs. exotic (phantom)

The reconstruction logic is "propose geometry → reconstruct the matter that
sources it." A wormhole/warp geometry generally requires **exotic** (phantom,
`rho <= 0`) matter to satisfy the Hamiltonian constraint, which is why
`min_rho_req` goes negative for such seeds.

Two matter models are available (both in `Source/Matter/`):

- `ScalarField` — canonical, `rho = 1/2 Pi^2 + 1/2 |grad phi|^2 + V >= 0`.
- `ExoticScalarField` — phantom, every stress-energy component scaled by
  `-recipe_support_strength`, giving `rho <= 0` (NEC/WEC violating).

Select at runtime:

```
recipe_exotic_matter    = 1     # evolve ExoticScalarField (phantom)
recipe_support_strength = 1.0   # magnitude of the phantom support
```

When set, the exotic matter is used consistently in the **RHS evolution**, the
**Hamiltonian/momentum constraints**, and the **energy-condition diagnostic**, so
the matter that is actually evolved matches the matter the geometry was
reconstructed for. With `recipe_exotic_matter = 0` the canonical `ScalarField` is
used. The `grteclyn-wrapper` sets `recipe_exotic_matter = 1` automatically
whenever `--phantom` is passed.

## Two findings worth your attention

1. **Matter-sector EC is a null result by construction with a canonical field.**
   A canonical `ScalarField` has `rho >= 0`, so its NEC/WEC are `~0`; `--phantom`
   alone only shapes the initial data. With `recipe_exotic_matter = 1` the evolved
   matter is genuinely exotic and the `matter_*` columns go negative exactly where
   the geometry demands it (verified on Ellis–Bronnikov: `matter_min_NEC ~ -0.07`,
   integrated NEC violation `~2.1`). The curvature invariants and the general FTL
   measure read the geometry directly and are meaningful either way; the effective
   stress energy `T^eff = G / 8pi` is what the `matter_*` columns cannot see and is
   computed post-hoc by `warpfactory.py`.
2. **Evolved-data FTL / effective-EC needs more plot vars.** Plotfiles currently
   store only `chi`, `K`, `lapse` — not the full `h_ij`/shift — so the general FTL
   measure (`ftl_general.py`) runs on the `t=0` reconstructed slice. To run it (and
   the effective EC) on evolved spacetimes, add the metric components to
   `amr.plot_vars` and feed the plotfile grid into `operational_ftl_on_grid` /
   `warpfactory.py`.

## Build & run

```bash
# Build (single GPU, no MPI)
make COMP=gnu USE_CUDA=TRUE USE_MPI=FALSE CUDA_ARCH=90 -j

# Run directly
./main3d.gnu.CUDA.ex params.txt

# Or via the closed-loop wrapper (auto-sets exotic matter for phantom seeds)
BUILD=0 SEED_NAME=ellis_bronnikov \
  bash ../../grteclyn-wrapper/scripts/run_radialrecipe_gpu_smoke.sh
```
