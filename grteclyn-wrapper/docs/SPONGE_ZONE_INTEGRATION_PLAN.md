# Numerical Sponge Zone — Integration Plan

> Status: planning document. No code changes have been made yet.

## Feasibility Summary

Yes, a numerical sponge zone can be added to this GRTeclyn codebase.

- **No existing sponge zone**: A search for "sponge" only returns a note in `research/RL/research.md` describing the concept; there is no implementation.
- **Boundary conditions**: the production run (`Examples/SupportedWormholeCollapse/params_2gpu.txt`) already uses non-periodic boundaries (`isPeriodic = 0 0 0`), Sommerfeld on the high faces (`hi_boundary = 1 1 1`), and reflective on the low faces (`lo_boundary = 2 2 2`). The `RadialRecipe` example (`Examples/RadialRecipe/params.txt`) uses `isPeriodic = 0 0 0`, `hi_boundary = 1 1 1`, `lo_boundary = 1 1 2` (Sommerfeld on all faces except reflective on the low-z face).
- **Dissipation infrastructure**: Kreiss-Oliger dissipation is already applied uniformly via the `sigma` parameter in `CCZ4RHS` / `CCZ4RHSWithMatter` (`Source/Grids/FourthOrderDerivatives.hpp::add_dissipation`).
- **Insertion point**: the RHS is assembled in `specificEvalRHS` (per example) and boundary conditions are enforced in `GRAMRLevel::advance` after the RHS call. Either location can host the sponge layer.
- **Target campaign**: the immediate driver is the `RadialRecipe` MAP-ElitesDynamics campaign (`research/neuralspacetime/MapElitesDynamics.md`), which runs `RadialRecipe` at `N=128 / L=64 / max_level=1` (QD) and `N=256 / L=128 / max_level=3` (HQ) with Sommerfeld boundaries and physics centered at `recipe_params.grid_center`.

The sponge zone can therefore be implemented as a spatially-ramped extra dissipation (or a damping term) in a finite outer layer, without rewriting the CCZ4 or boundary-condition machinery.

---

## Where to Add the Sponge Zone

### Option A — Per-example, in `SupportedWormholeLevel::specificEvalRHS`

Add a second `ParallelFor` after the existing `ccz4rhs` call that computes a sponge damping increment on `a_rhs` inside the outer layer.

**Files to modify**
- `Examples/SupportedWormholeCollapse/SupportedWormholeLevel.cpp`
- `Examples/SupportedWormholeCollapse/SimulationParameters.hpp` (read new knobs)
- optionally `Examples/SupportedWormholeCollapse/params_2gpu.txt` (turn on / tune)

**Pros**
- Least invasive; no shared headers change.
- Easy to test against the existing paper run.
- Can be tuned without affecting other examples.

**Cons**
- Other examples (BinaryBH, ScalarField, etc.) do not get the sponge.
- Slightly more boilerplate if later copied to other examples.

**When to choose this**
If the immediate goal is to improve the `SupportedWormholeCollapse` production runs and validate the idea before generalizing.

---

### Option B — Shared, inside `CCZ4RHSWithMatter`

Make `CCZ4RHSWithMatter` accept an optional `SpongeZone` object and use a position-dependent dissipation coefficient inside `add_dissipation`.

**Files to modify**
- New: `Source/Grids/SpongeZone.hpp` (and `Make.package`)
- `Source/Matter/CCZ4RHSWithMatter.hpp` / `.impl.hpp` (add sponge member / constructor overload)
- `Source/CCZ4/CCZ4RHS.hpp` / `.impl.hpp` (possibly expose a spatial sigma hook)
- `Source/GRTeclynCore/SimulationParametersBase.hpp` (read shared knobs)
- Every example that constructs `CCZ4RHSWithMatter` (e.g. `SupportedWormholeLevel.cpp`, `RotatingWormholeCollapse`, etc.)

**Pros**
- Reusable across all matter examples.
- The sponge is part of the physical RHS computation, which is the most consistent place.

**Cons**
- Touches shared template code; changes ripple to all matter examples.
- Requires careful handling of the non-matter `CCZ4RHS` base class.

**When to choose this**
If the sponge is intended to become a standard feature for all CCZ4+Matter evolutions.

---

### Option C — Generic post-RHS step in `GRAMRLevel::advance`

Apply a sponge operator after `specificEvalRHS` returns and after `m_boundaries.apply_sommerfeld_boundaries` is called.

**Files to modify**
- New: `Source/Grids/SpongeZone.hpp` (and `Make.package`)
- `Source/GRTeclynCore/GRAMRLevel.cpp` (add the post-RHS call)
- `Source/GRTeclynCore/SimulationParametersBase.hpp` or `AMReXParameters.hpp` (read knobs)

**Pros**
- Automatically covers every example, including matter-free ones.
- Sponge is treated as a numerical boundary-layer fix, not a physics change.

**Cons**
- Applied after Sommerfeld RHS enforcement; may partially undo or conflict with boundary treatment.
- Less obvious to example authors where the extra damping comes from.

**When to choose this**
If the goal is a global numerical fix that should be transparent to all examples.

---

## Recommended Insertion for `RadialRecipe` — `RadialRecipeMatter::eval_rhs`

For the MAP-ElitesDynamics campaign the cleanest place is the shared matter dispatch in `Examples/RadialRecipe/RadialRecipeMatterDispatch.hpp`. This function already centralizes all five matter variants (`GRTresnaIndependentScalars`, `ComplexScalarField`, `BiComplexScalarField`, `ExoticScalarField`, and `ScalarField`), so a single post-`ccz4rhs` sponge call covers every campaign configuration without duplicating code.

**Files to modify**
- New: `Source/Grids/SpongeZone.hpp` (and `Source/Grids/Make.package`)
- `Examples/RadialRecipe/SimulationParameters.hpp` (read `sponge_zone.*` knobs)
- `Examples/RadialRecipe/RadialRecipeMatterDispatch.hpp` (apply sponge after each `ccz4rhs` ParallelFor)
- `Examples/RadialRecipe/params.txt` (add example tuning block)
- Optional: `grteclyn-wrapper/scripts/campaigns/lib/*.sh` or input templates used by the MAP-ElitesDynamics campaign so the sponge is on by default for new runs.

**How to wire it**

1. Read the sponge parameters in `SimulationParameters::read_shared_params` or `read_recipe_params`:
   ```cpp
   pp.load("sponge_zone.enabled", sponge_zone.enabled, false);
   pp.load("sponge_zone.inner_radius", sponge_zone.inner_radius, 0.0);
   pp.load("sponge_zone.outer_radius", sponge_zone.outer_radius, 0.0);
   pp.load("sponge_zone.sigma_factor", sponge_zone.sigma_factor, 1.0);
   pp.load("sponge_zone.damping_sigma", sponge_zone.damping_sigma, 0.0);
   pp.load("sponge_zone.ramp_power", sponge_zone.ramp_power, 2);
   pp.load("sponge_zone.center", sponge_zone.center, recipe_params.grid_center);
   ```

2. Pass the configured `SpongeZone` object to `RadialRecipeMatter::eval_rhs`, e.g.:
   ```cpp
   inline void eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                        const SimulationParameters &params, double dx,
                        const std::array<double, AMREX_SPACEDIM> &center,
                        double time)
   {
       SpongeZone sponge(params.sponge_zone, dx);
       ...
   }
   ```

3. After each `ccz4rhs` ParallelFor, add a second `ParallelFor` over `a_rhs` that calls the sponge operator:
   ```cpp
   amrex::ParallelFor(
       a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
       { sponge.apply(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
   ```
   This adds the spatially-ramped extra dissipation (and optional damping) to the RHS already computed by `CCZ4RHSWithMatter`.

**Why this is the best fit for `RadialRecipe`**
- `RadialRecipeMatter::eval_rhs` is the single dispatch point for every matter model used in the campaign, so one change covers all search variants.
- The sponge is applied after the physical RHS but before the Sommerfeld boundary RHS enforcement in `GRAMRLevel::advance`, so it supplements rather than fights the boundary conditions.
- No changes to `CCZ4RHSWithMatter` or `GRAMRLevel` are required; the risk to other examples is zero.

**Pros**
- Minimal, campaign-focused change.
- Works for all matter sectors in `RadialRecipe`.
- Easy to A/B test by toggling `sponge_zone.enabled`.

**Cons**
- Not automatically available to `SupportedWormholeCollapse` or `BinaryBH` unless the same code is duplicated or later moved to `CCZ4RHSWithMatter` (Option B).

---

## Sponge Mechanism Options

### Mechanism 1 — Spatially Varying Kreiss-Oliger Dissipation

In the outer layer, multiply the existing dissipation coefficient `sigma` by a profile `f(r)` that ramps from 0 in the interior to 1 at the boundary.

```
sigma_eff(r) = sigma * (1 + (sponge_sigma_factor - 1) * f(r))
```

or, equivalently,

```
sigma_eff(r) = sigma_inner * (1 - f(r)) + sigma_outer * f(r)
```

A smooth ramp is typically:

```
f(r) = 0                                          for r <= r_sponge_in
f(r) = ((r - r_sponge_in) / (r_sponge_out - r_sponge_in))^p   for r_sponge_in < r < r_sponge_out
f(r) = 1                                          for r >= r_sponge_out
```

**Pros**
- Directly reuses the existing `FourthOrderDerivatives::add_dissipation` / `SixthOrderDerivatives::add_dissipation` machinery.
- Conceptually matches the note in `research/RL/research.md`.

**Cons**
- High dissipation can slow waves and alter phase; needs tuning.
- Must be applied to variables that already receive KO dissipation.

---

### Mechanism 2 — Damping Toward Asymptotic Values

Add a forcing term to the RHS in the outer layer:

```
∂t u -> ∂t u - sigma_sponge(r) * (u - u_asymptotic)
```

where `u_asymptotic` is the flat-space / background value (e.g. `chi=1`, `h_ij=δ_ij`, `lapse=1`, others=0).

**Pros**
- Very effective at killing incoming reflections.
- Explicitly drives variables toward the same values used by the Sommerfeld boundary condition.

**Cons**
- Requires a per-variable asymptotic value map.
- Can violate constraints if applied too aggressively or too close to the source.

---

### Mechanism 3 — Combined Dissipation + Damping

Use the ramped Kreiss-Oliger dissipation to smooth high-frequency reflections and add a weak damping term toward asymptotic values to absorb low-frequency junk.

**Pros**
- Best suppression across a broad frequency range.
- Robust for long runs.

**Cons**
- More parameters to tune.
- Slightly more code.

---

## Which Variables to Act On

### Scope A — All State Variables

Apply the sponge to every component of `a_rhs` (or every component that receives KO dissipation).

**Pros**
- Simplest implementation.
- Reflections can appear in any variable, so blanket treatment is safest.

**Cons**
- May over-damp the scalar field if it is supposed to carry physical waves out of the domain.

---

### Scope B — Geometry/Gauge Only, Exclude Scalar Matter

Sponge the CCZ4 metric/gauge variables (`chi`, `h_ij`, `A_ij`, `K`, `Theta`, `Gamma^i`, `lapse`, `shift_i`) but leave the scalar field (`phi`, `Pi`) untouched.

**Pros**
- Preserves the physical scalar evolution in the bulk and near the boundary.

**Cons**
- Scalar radiation can still reflect unless the metric is sufficiently damped.
- Requires a variable list; more code.

---

### Scope C — Configurable Per-Variable List

Read a list of variable indices from the parameter file (e.g. `sponge_vars = chi h11 h22 h33 lapse`) and only sponge those components.

**Pros**
- Maximum flexibility.

**Cons**
- More parameters and more error handling.

---

## Suggested Shared Implementation: `Source/Grids/SpongeZone.hpp`

A reusable header-only class is recommended regardless of which insertion point is chosen.

```cpp
struct SpongeZoneParams
{
    bool enabled = false;
    double inner_radius = 0.0;     // radius where ramp starts
    double outer_radius = 0.0;     // radius where ramp ends (usually domain boundary)
    double sigma_factor = 1.0;     // multiplier at outer_radius for dissipation-only
    double damping_sigma = 0.0;    // damping coefficient at outer_radius
    int ramp_power = 2;            // polynomial order of the ramp
    std::array<double, AMREX_SPACEDIM> center{};

    void read(GRParmParse& pp);
    AMREX_GPU_DEVICE double profile(double r) const;
    AMREX_GPU_DEVICE double effective_sigma(double base_sigma, double r) const;
    AMREX_GPU_DEVICE double damping_coefficient(double r) const;
};
```

The profile is computed once per cell from the cell-centered coordinate:

```cpp
r = sqrt((x - cx)^2 + (y - cy)^2 + (z - cz)^2)
f = clamp(((r - r_in) / (r_out - r_in))^p, 0, 1)
```

Helper functions:
- `effective_sigma(base_sigma, r)` returns `base_sigma * (1 + (sigma_factor - 1) * f(r))`.
- `damping_coefficient(r)` returns `damping_sigma * f(r)`.

This class can be used in any of the three insertion points above.

---

## Parameter File Additions

```
# Sponge zone
sponge_zone.enabled = 1
sponge_zone.inner_radius = 25.0       # ~80% of L/2 for L=64
sponge_zone.outer_radius = 32.0       # domain boundary
sponge_zone.sigma_factor = 10.0       # 10x KO dissipation at boundary
sponge_zone.damping_sigma = 0.1       # optional damping term
sponge_zone.ramp_power = 2
sponge_zone.center = 0.0 0.0 0.0
```

The existing `sigma` parameter remains the base dissipation coefficient; the sponge zone only modifies it in the outer layer.

### Suggested defaults for `RadialRecipe` / MAP-ElitesDynamics

For `L_full = 64.0` with physics centered at `recipe_params.grid_center` (typically the domain center for the campaign), the sponge can start at ~80% of the half-domain and end at the boundary:

```
sponge_zone.enabled = 1
sponge_zone.inner_radius = 25.0       # r ~ 0.8 * (L/2)
sponge_zone.outer_radius = 32.0       # domain boundary / sqrt(3) for box corner
sponge_zone.sigma_factor = 10.0         # 10x base KO dissipation at boundary
sponge_zone.damping_sigma = 0.0         # start with dissipation only; add damping if needed
sponge_zone.ramp_power = 2
sponge_zone.center = 32.0 32.0 32.0   # match recipe_params.grid_center
```

For the HQ replay (`L_full = 128.0`), scale the radii proportionally:

```
sponge_zone.inner_radius = 50.0
sponge_zone.outer_radius = 64.0
```

If the sponge is meant to be a box-shaped layer based on the maximum coordinate distance to any face rather than a spherical shell, use a box profile instead:

```
r = max(|x - cx|, |y - cy|, |z - cz|)
```

This avoids leaving unsponged corners and is often preferred for Cartesian boxes.

---

## Recommended Minimal Path

For the fastest, lowest-risk validation on the target `RadialRecipe` / MAP-ElitesDynamics campaign:

1. **Create** `Source/Grids/SpongeZone.hpp` with a GPU-ramp profile (box or spherical).
2. **Register** it in `Source/Grids/Make.package`.
3. **Add** `SpongeZoneParams` to `Examples/RadialRecipe/SimulationParameters.hpp`.
4. **Modify** `Examples/RadialRecipe/RadialRecipeMatterDispatch.hpp::eval_rhs` to apply the sponge as a second `ParallelFor` after each `ccz4rhs` call.
5. **Add** a `sponge_zone.*` block to `Examples/RadialRecipe/params.txt` and to the campaign input templates.
6. **Test** with a short QD or HQ run and compare `Psi_4` / curvature-invariant / energy-condition diagnostics for late-time boundary reflections.
7. **Generalize** to `SupportedWormholeCollapse` or Option B (shared inside `CCZ4RHSWithMatter`) only after the RadialRecipe campaign tuning is stable.

For a `SupportedWormholeCollapse` validation the same steps apply, but the insertion point is `SupportedWormholeLevel::specificEvalRHS` and the sponge center should match the wormhole center (typically the origin for the octant run).

---

## Risks and Considerations

- **Stability**: high `sigma_factor` or `damping_sigma` can violate the Courant condition; start with small factors and increase gradually.
- **AMR**: the sponge must be applied on every level, using the correct `dx` and physical coordinates. A level-local operation naturally handles this.
- **Interaction with Sommerfeld**: the sponge should be located inside the domain, well before the mathematical boundary, so it does not replace the Sommerfeld condition but supplements it.
- **Constraint violation**: strong damping can drive the metric away from the constraint-satisfying evolution; monitor Hamiltonian and momentum constraints.
- **Performance**: extra spatial profile adds a few flops per cell; negligible compared to the CCZ4 RHS.
- **Low-face reflective BC**: `SupportedWormholeCollapse` uses reflective boundaries on all low faces; `RadialRecipe` uses reflective only on the low-z face (`lo_boundary = 1 1 2`). A sponge layer is usually paired with Sommerfeld boundaries on all faces; consider switching `lo_boundary` to `1 1 1` if reflections from the low faces are a problem.
- **Matter pump interaction**: `RadialRecipe` can inject a time-dependent matter pump (`RLMatterPumpParams`) that may excite the scalar field near the boundaries. If the sponge is too strong or too close to the source, it can fight the pump and suppress the physical signal. Keep the sponge well outside the lump region (e.g. `r > 2 * recipe_basis_radius_max`).
- **Campaign templates**: the MAP-ElitesDynamics campaign inputs are generated by wrapper scripts. Ensure the sponge parameters are added to the templates (e.g. `grteclyn-wrapper/scripts/campaigns/lib/*.sh` or GRTresna-generated input files) so they are consistent across QD, CMA-ES, and HQ stages.

---

## Verification Steps

### `RadialRecipe` / MAP-ElitesDynamics

1. Build the `RadialRecipe` executable.
2. Run a short QD-style test (`L_full = 64`, `N_full = 128`, `max_level = 1`, `stop_time = 5.0`) with the sponge off and on.
3. Compare `Psi_4` extraction (if `GRTECLYN_PSI4=1`) at `r = 8`, `12`, and `24` from `recipe_params.grid_center` for late-time reflections.
4. Plot `L2` Hamiltonian / momentum constraints from `data/constraint_norms.dat`; ensure they do not grow.
5. Check `data/energy_conditions.dat` and `data/curvature_invariants.dat`; the sponge should not introduce new regions of strong EC violation or curvature spikes.
6. Visualize `chi`, `lapse`, or `phi` near the box boundaries to confirm the sponge layer is active.
7. Sweep `sponge_zone.sigma_factor` (e.g. `2.0 → 20.0`) and `sponge_zone.inner_radius` to find the minimal effective setting.

### `SupportedWormholeCollapse`

1. Build the `SupportedWormholeCollapse` executable.
2. Run a short test (e.g. `stop_time = 5.0`) with the sponge off and on.
3. Compare `Psi_4` extraction at `r = 8` and `r = 16` for late-time reflections.
4. Plot `L2` Hamiltonian / momentum constraints; ensure they do not grow.
5. Visualize `chi` or `lapse` near the boundaries to confirm the sponge layer is active.
6. Sweep `sponge_zone.sigma_factor` and `sponge_zone.inner_radius` to find the minimal effective setting.

---

## Open Decisions

- Insertion point (A, B, C, or the RadialRecipe-specific dispatch)
- Sponge profile: spherical shell (`r = sqrt(...)`) or box-aligned layer (`r = max(|x-cx|, ...)`)
- Mechanism (1, 2, or 3)
- Variable scope (A, B, or C)
- Whether to also switch `lo_boundary` from reflective to Sommerfeld in `RadialRecipe` and `SupportedWormholeCollapse`
- Default values for the production `params_2gpu.txt` and for the MAP-ElitesDynamics campaign templates
