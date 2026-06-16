# RadialRecipe — Debug notes

## Update (8 June 2026): matter–geometry decoupling in the matter-first (GRTresna) path

**Verdict: the matter-first / MAP-Elites "shell" results (paper Sec. 4b and the
MAP-Elites shell HQ-promotion table) cannot be trusted as constraint-satisfying
physical evolutions for any candidate that contains exotic and/or overlapping
lumps — which is all the headline ones (`eval_000063` = 3 exotic lumps of 5; the
promoted shell elites `eval_000106/117/094/011/...` all use an exotic wedge plus
overlapping Fibonacci lumps).** The geometry GRTeclyn evolves is sourced by a
*different* stress-energy than the one GRTresna used to satisfy the constraints,
so the evolution starts off-constraint and the "evolved operational FTL"
(`F_op^ev`) is a constraint-relaxation / gauge transient, not a sustained
physical channel. This is the same root cause (inconsistent source, `G_ab ≠
8πT_ab`) documented as **P2** for the RotatingWormhole route.

Important: this is **not** a `.gridinit` column-misalignment bug. For
RadialRecipe the loader is index-based but harmless — GRTresna writes exactly the
27 `chi…Pi` components (`../GRTresna/Source/Variables/GRChomboVariables.hpp`) in
the same order as `Examples/RadialRecipe/StateVariables.hpp` (`c_phi=25`,
`c_Pi=26`), so geometry **and** matter `phi/Pi` load into the right slots
(`Examples/RadialRecipe/ExternalGridInitialData.hpp:69-92`). The decoupling is in
the **matter model**, not the loader.

### The inconsistency (3 independent mechanisms)

GRTresna solves `G_ab = 8π T_ab^solve`, writes the geometry **and** the combined
field `phi = phi_bg + Σ_k lump_phi_k`, `Pi = Pi_bg + Σ_k lump_pi_k`
(`../GRTresna/Examples/ScalarFieldBH/MyMatterFunctions.cpp:19-41`). GRTeclyn then
re-derives `T_ab^evolve` from that `(phi, Pi)` with a **canonical, massless**
scalar field. `T_ab^solve ≠ T_ab^evolve` because:

1. **Exotic sign flip (dominant).** In `T_ab^solve` each exotic lump's
   kinetic + momentum density enters with `sign = -1`
   (`../GRTresna/Source/Matter/ScalarField.cpp:107-114`). The written `phi/Pi`
   carry **no** sign flip. GRTeclyn evolves canonical `ScalarField`
   (`rho = ½Π² + ½χ|∇φ|² + V ≥ 0`,
   `Examples/RadialRecipe/RadialRecipeLevel.cpp:415-443`), so in every exotic
   region the Hamiltonian source has the **wrong sign** (off by ~`2·|rho_exotic|`).
2. **Dropped cross terms.** `T_ab^solve` is an *independent-field* sum:
   `Σ_k |∇φ_k|²` with inter-lump cross terms dropped
   (`ScalarField.cpp:66-72,100-115`). GRTeclyn uses the true combined field, so
   `|∇φ|² = |Σ_k ∇φ_k|²` **includes** all cross terms. Mismatch is finite wherever
   lumps overlap — and the production "compact" shell deliberately overlaps lumps.
3. **Potential mismatch.** `T_ab^solve` adds `V = ½(m·φ)²` with
   `scalar_mass` (default `0.1`, `grteclyn-wrapper/.../grtresna/solver.py:82`,
   `MyMatterFunctions.cpp:14-17`). GRTeclyn's `DefaultPotential` is `V≡0`
   (`Source/Matter/DefaultPotential.hpp:26-29`). Mismatch `∝ m²φ²`.

A single real scalar with one global flag **cannot** represent mechanism (1):
mixed per-lump signs are not a real field configuration. So the inconsistency is
structural, not just a missing flag.

### Why the search never caught it

`recipe_exotic_matter` defaults to `false`
(`Examples/RadialRecipe/SimulationParameters.hpp:36`) and is **never set on the
grtresna path** (`core/evaluation.py:72-73` guards it with `and not grtresna`;
the QD launcher `scripts/campaigns/qd/run.sh` does not set it). The
search gates only on:
- GRTresna's **own** internal Ham/Mom residual (`parse_convergence`), and
- a **geometry-only** FTL probe on the `.gridinit`
  (`metrics/ftl_solved_geometry.py`, `core/evaluation.py:152-197`).

Neither sees GRTeclyn's post-load `T_ab`. The paper's quoted "GRTresna Ham/Mom
residuals 0.65–1.53%" are GRTresna-internal, **not** what GRTeclyn evolves.

### Cluster verification (the runs are not in this repo)

GRTeclyn already computes its own `t=0` constraint with the *evolution* matter
model and writes it (`RadialRecipeLevel.cpp:495-645` →
`<output>/data/constraint_norms.dat`, columns `L2_Ham L2_Mom …`). To confirm:

1. For an exotic shell elite (e.g. `l128n256_qd_eval106`), read the **first row**
   (`t=0`) of `constraint_norms.dat`. Expectation: `L2_Ham` is **orders of
   magnitude larger** than GRTresna's reported ~1% (it should reflect the wrong
   matter source), and it does **not** decrease to the GRTresna value.
2. Compare against a **purely canonical, single, well-separated, massless** lump
   candidate: its `t=0 L2_Ham` should be small and stable — isolating the effect.
3. Optional offline check (no GPU): `read_gridinit()` the `.gridinit`, recompute
   `rho = ½Π² + ½χ|∇φ|² + V` from `phi/Pi` and compare to the vacuum
   `rho_req = H_vac/16π` from the loaded geometry. They should disagree in the
   exotic/overlap regions.

### Fix options (decide before re-running)

- **(A) Minimal "make results honest" gate.** Add a GRTeclyn-side `t=0`
  constraint gate to the search: after load, reject candidates whose `L2_Ham/Mom`
  (with the *actual* evolution matter model) exceed a threshold. Catches the
  decoupling instead of fixing it; cheapest, no physics change.
- **(B) Restrict to a representable matter basis.** Forbid mixed signs and large
  lump overlap; set `scalar_mass = 0` in the GRTresna solve **or** add the mass
  potential to GRTeclyn; set `recipe_exotic_matter` consistently
  (all-canonical or all-phantom with matching `recipe_support_strength`). Then a
  single field reproduces the source up to FD error. Shrinks the search space.
- **(C) Evolve the consistent stress-energy.** Make GRTresna emit the effective
  `(rho, S_i, S_ij)` it actually solved and have GRTeclyn evolve a *genuinely
  co-evolving* matter model that reproduces it (NOT a frozen prescribed source —
  that is P2 and blows up, see RotatingWormhole Debug.md). Largest change; only
  this lets exotic+mixed configs be physical.

Recommended order: (A) immediately to quantify/quarantine existing results, then
(B) for a trustworthy reduced search, with (C) as the real long-term fix.

### Key code references

- GRTresna source EM tensor (sign flip, cross terms, +V):
  `../GRTresna/Source/Matter/ScalarField.cpp:51-123`
- GRTresna painted field/momentum + potential:
  `../GRTresna/Examples/ScalarFieldBH/MyMatterFunctions.cpp`,
  `../GRTresna/Examples/ScalarFieldBH/MatterParams.hpp`
- GRTresna output var order (27 `chi…Pi`):
  `../GRTresna/Source/Variables/GRChomboVariables.hpp`, `Source/Tools/WriteOutput.H`
- GRTeclyn loader / state vars / matter RHS / potential:
  `Examples/RadialRecipe/ExternalGridInitialData.hpp`,
  `Examples/RadialRecipe/StateVariables.hpp`,
  `Examples/RadialRecipe/RadialRecipeLevel.cpp`, `Source/Matter/DefaultPotential.hpp`
- Search wiring (no GRTeclyn constraint gate, exotic flag unset):
  `grteclyn-wrapper/src/grteclyn_wrapper/core/evaluation.py`,
  `grteclyn-wrapper/src/grteclyn_wrapper/search/optimize.py` (`build_grtresna_config`),
  `grteclyn-wrapper/scripts/campaigns/qd/run.sh`
