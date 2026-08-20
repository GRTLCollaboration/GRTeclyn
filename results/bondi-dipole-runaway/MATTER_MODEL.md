# Matter model, initial data, and the code modifications

Paths are repo-relative. `GRTeclyn/` is the evolution code, `GRTresna/` the
constraint solver, `grteclyn-wrapper/` the Python driver that configures both.
Diffs of everything described here are in [`patches/`](patches/).

---

## 1. The bicomplex scalar model

Two independent complex scalar fields share one spacetime:

- **Φ₊** — the *canonical* sector, gravitational sign +1
- **Φ₋** — the *phantom* sector, gravitational sign −1

Both obey the **same** potential,

```
V(|Φ|²) = ½ m² |Φ|² − ¼ λ |Φ|⁴ + ⅙ μ |Φ|⁶ ,        m = 1
```

The sextic term (μ > 0 with λ > 0) is required: a purely attractive quartic is
the critical case in 3D and either collapses or disperses. The *same* potential
is compiled into both codes — `GRTeclyn/Source/Matter/ComplexScalarPotential.hpp`
and `GRTresna/Source/Matter/ComplexScalarField.cpp::potential_value` — because
if the solved metric backed a different soliton than the one evolved, the lump
would relax and radiate at t = 0.

### Where the sign flip lives — and where it does not

This is the physical heart of the model, and it is one sign in one function.

**Stress-energy (sign flips):** `GRTeclyn/Source/Matter/BiComplexScalarField.impl.hpp`
accumulates each field's stress-energy through a helper that takes the field's
gravitational sign `fs`:

```cpp
out.rho    += fs * (Pi_k * Pi_k + 0.5 * Vt_k);
out.j[i]   += fs * (-Pi_k * dphi_i);
out.S[i][j]+= fs * (-0.5 * vars.h(i,j) * Vt_k / chi + dphi_i * dphi_j);
out.rho    += fs * V_of_phi;
```

called once with `+1.0` for Φ₊ and once with `−1.0` for Φ₋.

**Field equations (sign does NOT flip):** in `add_matter_rhs` both fields get an
identical Klein–Gordon right-hand side —

```cpp
rhs[c_phi]   = lapse * Pi1p + advec.phi1p();      // canonical
rhs[c_Pi]    = lapse * (K * Pi1p − dVp1) + …
rhs[c_phi_m] = lapse * Pi1m + advec.phi1m();      // phantom: same form
rhs[c_Pi_m]  = lapse * (K * Pi1m − dVm1) + …
```

So both sectors have **positive inertial and passive gravitational mass** and
**signed active gravitational mass** — exactly Bondi's premise, and the reason
the pair self-accelerates instead of simply repelling. The phantom sector is a
genuine NEC-violating source, not a sign trick in the initial data.

State layout (`GRTeclyn/Examples/RadialRecipe/StateVariables.hpp`): Φ₊ occupies
`c_phi, c_Pi, c_phi2, c_Pi2`; Φ₋ reuses the `phi_lump1/2` slots as
`c_phi_m, c_Pi_m, c_phi2_m, c_Pi2_m`. This is why plot variables named
`phi_lump1` carry the phantom field, and why `movie_phi_lump0_z.mp4` shows only
the canonical blob.

On the solver side, `GRTresna/Source/Matter/ComplexScalarField.cpp`
`compute_emtensor` applies the same idea **per lump**, so one constraint solve
can carry mixed-sign matter:

```cpp
const Real sign = (L.exotic != 0) ? -1.0 : 1.0;
rho_total += sign * (0.5*pi1_k*pi1_k + 0.5*pi2_k*pi2_k + 0.5*chi*(grad1_sq+grad2_sq) + V_k);
Si[i]     += sign * (-(pi1_k*dphi1_k[i] + pi2_k*dphi2_k[i]));
```

Cross terms between lumps are dropped (each lump is an independent field),
mirroring the real-scalar path.

## 2. Coupling rung

All rungs in the campaign keep λ²/μ = 4.8, fixing ω_min = √(1 − 3λ²/16μ) = 0.316,
so ω = 0.55 always means 45 % binding. What changes is the amplitude scale
√(3λ/4μ):

| rung | λ | μ | core amplitude | dressed ω = 0.55 star exists? |
|---|---|---|---|---|
| strong | 640 | 85 333 | 0.075 | — (collapses) |
| weak | 2 560 | 1 365 333 | 0.0375 | **no** |
| mid | 1 280 | 341 333 | 0.053 | **no** |
| **ultraweak** | **10 240** | **21 845 333** | **0.0197** | **yes** |

The ultraweak rung is the first at which a gravitationally dressed ω = 0.55 star
exists at all — the discovery that unblocked the campaign
([`DEBUGGING.md`](../../research/bondi_dipole/docs/DEBUGGING.md) §2).

## 3. Dressed-star initial data

### 3.1 The star ODE

`grteclyn-wrapper/src/grteclyn_wrapper/grtresna/profiles/boson_star_ode.py`
integrates a static, spherically symmetric, maximally sliced star in isotropic
coordinates with state `[φ₀, φ₀′, ψ, ψ′, α, α′]`:

```
ψ″ + (2/r) ψ′ = −2π ψ⁵ ρ                                  (Hamiltonian constraint)
α″ + (2/r + 2ψ′/ψ) α′ = 4π ψ⁴ α (ρ + S)                   (maximal slicing)
φ₀″ + (2/r + α′/α + 2ψ′/ψ) φ₀′ = ψ⁴ (dU − ω²/α² φ₀)       (Klein–Gordon)
ρ     = ½ ω²/α² φ₀² + ½ ψ⁻⁴ φ₀′² + U
ρ + S = 2 ω²/α² φ₀² − 2U
```

**The phantom star is obtained by `gravity_sign = −1`, which multiplies `ρ` and
`ρ + S` — and nothing else:**

```python
# Only the METRIC sources flip for the phantom sector; the Klein-Gordon
# equation below keeps its canonical form.
rho = gravity_sign * rho
rho_plus_s = gravity_sign * rho_plus_s
```

The Klein–Gordon line is untouched, exactly mirroring the evolution code. The
result is a soliton bound by its own self-interaction whose *gravity is
repulsive*: ψ < 1, α > 1 in the core (against ψ > 1, α < 1 for the canonical
star), and negative ADM mass.

### 3.2 Solving at fixed frequency

Naïve ω-shooting fails on this branch. Heavy-branch sextic stars sit with φ_c a
fraction of a percent below the effective-potential top, so at fixed φ_c the
eigenvalue in ω is an exponentially thin needle — a shooter steps over it and
lands on the light branch (observed: asking for the weak-rung ω = 0.55 star
returned a 5× lighter ω = 0.75 star). `solve_selfgrav_at_omega` therefore
mirrors the flat-space parameterisation: **bisect φ_c at fixed integration-frame
ω, then outer-iterate ω_int so the rescaled ω_phys = ω_int / α_∞ matches the
request.** φ_c becomes a solved *output*, and the returned eigenvalue matches
the request to ~1 part in 10⁵ (see `stars/star_family.csv`).

Results at ω = 0.550: canonical φ_c = 0.019695, ADM = +0.063951, α(0) = 0.9773;
phantom φ_c = 0.019589, ADM = −0.076962, α(0) = 1.0249.

### 3.3 Painting the star, and the two subtleties that mattered

1. **The amplitude must be the star's own φ_c.** Both painters rescale the
   tabulated profile by `amp/φ_c`, so any other amplitude de-tunes the
   eigenstate. The writer back-writes `lump["amp"] = profile.phi_c`.
2. **A boson star is stationary, not static.** The U(1) momentum must be
   painted with the star's own lapse,

   ```
   Π_im = −(ω / α(r)) φ₀(r)
   ```

   which is why the emitted profile table carries a third column α(r) and why
   the post-solve repaint divides Π₂ by the solved lapse. Painting with α = 1
   gives the wrong kinetic energy and the seed radiates.

Mixed-sector pairs get **one table per sector** — `qball_profile.dat` and
`qball_profile_exotic.dat` — referenced per lump via `lump{k}_profile_path`,
because sharing the canonical table would seed the phantom off-equilibrium.

### 3.4 Constraint solve and evolution hand-off

GRTresna solves CTTKHybrid with these per-lump signed sources, then the wrapper
converts the output to a GRTeclyn gridinit, **repainting the matter block
analytically** from the same star tables (canonical lumps → Φ₊, phantom lumps →
Φ₋) while keeping GRTresna's solved metric. Solve tolerances were tightened for
this campaign to `NL_exit_tolerance = 0.1 %`, `NL_stall_tolerance = 0.002`; at
the old 1 % gate the residual radiated away as a visible χ ring that sloshed the
star's envelope.

## 4. Per-lump frequency (the equal-mass capability)

The equal-|ADM| cell needs the phantom star at ω = 0.56598 while the canonical
star stays at ω = 0.550 — but the painter had exactly one global U(1) phase
velocity. Seeding the 0.56598 profile while painting momentum at 0.550 would
reproduce the off-eigenstate bug the campaign had just eliminated, so the
frequency was threaded through per lump:

| layer | file | change |
|---|---|---|
| campaign knob | `search/optimize/config.py` | `trajectory_lump{k}_bs_omega` → that lump's `qball_omega` (0 ⇒ global) |
| solve + emit | `grtresna/solver/params.py` | group lumps by (sector, ω), solve one star per group, back-write per-lump `amp` / `bs_omega` / `profile_path`, emit `lump{k}_bs_omega` |
| Python repaint | `grtresna/fields/boson_star.py` | `_boosted_lump_fields` uses the lump's own `bs_omega` for Π_im (at-rest and boosted paths) |
| C++ struct + parser | `GRTresna/Source/Matter/BosonStarParams.hpp` | new `boson_lump_t::bs_omega`, parsed as `lump{k}_bs_omega`; `total_pi2` uses it |
| C++ painter + source | `GRTresna/Source/Matter/ComplexScalarField.cpp` | non-winding paint and `compute_emtensor` use `omega_k = bs_omega > 0 ? bs_omega : global` |

`bs_omega` is deliberately a **new field**, not a reuse of the existing
`boson_lump_t::omega`: that one is the rigid *rotation* rate used by
`lump_pi1`, and overloading it would have silently spun the stars.

Regression test: `grteclyn-wrapper/tests/grtresna/test_selfgrav_profile_wiring.py::test_selfgrav_pair_mixed_frequency_equal_mass`
asserts the per-lump solve, the `params.txt` contract, and that the Python
repaint's Π_im uses the phantom's own eigenvalue. Suite: 14 tests green.

Validation in production: the equal-mass cell was born at phantom
total = 17.083 / rms = 5.162 against the 17.05 / 5.16 prediction from the
0.56598 star — i.e. the new star was painted, at its own frequency, correctly.

## 5. Diagnostics added for this campaign

| diagnostic | file | what it gives |
|---|---|---|
| per-sector barycentres | `small_data/sector_barycenters.dat` (env `GRTECLYN_SECTOR_BARYCENTERS=1`) | separate integrals, centroids and rms radii for Φ₊ and Φ₋ — **the trajectory record** |
| confinement / grip | `small_data/confinement.dat` | peak amplitude, confined fraction, min χ |
| radiation | `small_data/psi4_*.dat` (env `GRTECLYN_PSI4=1`) | Ψ₄ modes and directionality at R = 8, 16 |

Plotfiles are consumed and purged during the run, so the barycentre stream is
the only recoverable trajectory record — worth knowing before designing a rerun.

## 6. Patch inventory

| patch | contents |
|---|---|
| `patches/grtresna-matter.diff` | the C++ matter modifications: per-lump exotic sign in `compute_emtensor`, tabulated star loader with lapse column, per-lump `profile_path`, per-lump `bs_omega` |
| `patches/wrapper-per-lump-omega.diff` | the per-lump frequency wiring + its regression test |
| `patches/wrapper-phantom-star-solver.diff` | phantom-star solving (`gravity_sign`), per-sector table emission, exotic-veto lift |
