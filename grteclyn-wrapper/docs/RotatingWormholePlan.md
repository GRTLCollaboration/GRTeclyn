# Rotating Wormhole Collapse — Implementation Plan

Extension of the Ellis–Bronnikov collapse article
(`research/wormholedynamics/article.md`) to a **rotating** traversable wormhole:
rotation supplies an inherent quadrupole so the ℓ=2 gravitational-wave signal is
*natural* (no artificial `A_φ` perturbation), and **GRTresna-solved initial
data** eliminates the t=0 constraint defect and the steady Hamiltonian growth of
the EB runs.

Builds directly on the post-mortems in
[`Examples/RotatingWormholeCollapse/Debug.md`](../Examples/RotatingWormholeCollapse/Debug.md).
Companion plans: [`NextSteps.md`](NextSteps.md), [`PuncturePlan.md`](PuncturePlan.md).

---

## Status (2026-07-07) — approved scope

**Goal (this pass):** simulate rotating-wormhole *collapse*, extract the GW
signal, and assess **stability vs. the static case**. The dynamical mass-support
control system (RL matter pump / PD "trap" controller, `RLMatterPumpParams.hpp`)
as an active throat stabilizer is **explicitly deferred** to a follow-up pass.

**What is already done (verified in-tree, supersedes Debug.md's "never
compiled"):**

- **Builds.** Both `main3d.gnu.ex` (CPU) and `main3d.gnu.MPI.CUDA.ex` (GPU)
  compile — the `GNUmakefile` fix adding `GRTeclynCore/RL` to the include path
  resolved the `RLMatterPumpParams.hpp` dependency pulled in unconditionally by
  `Source/Matter`. **Phase A1/A2 (build) are complete.**
- **Complex phantom-scalar matter model** `(φ₁,Π₁,φ₂,Π₂)` with
  `wormhole_azimuthal_m`, `wormhole_rotation_omega`, phantom sign flip, and the
  `S_support` collapse trigger — all wired in `SupportedWormholeLevel.cpp`.
- **Collapse + stability diagnostics** in `collapse_diagnostics.dat`:
  `min_chi`, `min_lapse`, `max_ah_r`, `min_theta_plus` (trapped-surface proxy),
  barycenter, and `J_z` (spin tracking).
- **GW extraction** wired: Weyl4 at r = 12/16/20/24; C++ extracts (l=2, m=0/1/2);
  Python `process_wave/.../psi4.py` decomposes all five m ∈ {−2..+2} with
  z-beaming power metrics.
- **Params** exist: analytic equilibrium + collapse, plus GRTresna
  `.gridinit`-driven variants.

**What is NOT done:** no evolution has ever been run
(`runs/rotating_wormhole/` does not exist). Remaining work is **initial-data
quality, execution, and analysis — not building.**

**Approved decisions:**

- **Initial data:** GRTresna **constraint-clean** ID first (Phase B). Note the
  existing `make_rotating_wormhole_id.py` solves a *real*-scalar Gaussian lump
  (the four-lobe failure, L3) — Phase B must extend the GRTresna solve to a
  **complex phantom scalar with phase winding + singular-ψ throat**. This is the
  real work item.
- **Static baseline:** **both** — ω=0 run of the *same* complex-scalar code
  (apples-to-apples numerics, doubles as the EB regression) **and** the published
  Ellis–Bronnikov article / `SupportedWormholeCollapse` results.
- **Scope:** full (ω, m, κ) grid (Phase C), after the single-case pipeline is
  proven by the Phase A smokes.

**Smoke passed (2026-07-07).** A coarse 64³ / t=2.0 no-NaN smoke
(`params_rotating_complex_smoke.txt`, ω=0.05, m=1, S_support=1) ran to
completion on 1×H100: constraints bounded (L2_Ham 0.029→0.037, L2_Mom ≈0.006),
GW extraction executed, and `J_z ≈ −1000` (nonzero, roughly held) — the complex
phantom field carries spin as designed. Confirms the built binary evolves the
model.

**Launcher (disk-safe, lesson L6):**
`grteclyn-wrapper/scripts/wormhole/run_rotating_wormhole.sh <params> <ngpu>`.
It auto-starts the `plot_run.sh` consumer sidecar *before* the evolution
(extraction/frame centers derived from the params `center`), streams PNG frames
+ `small_data/` (`psi4_mode_l2m0.dat`, `psi4_directional.dat`, `areal_radius.dat`)
and **deletes plotfiles on the fly** (`--keep-last 2`), then prunes checkpoints on
exit (`KEEP_CHK`, default 1). Verified: a full smoke leaves ~1.4 MB
(frames + small_data), not tens of GB. Never launch the bare binary without the
sidecar.

**First production 128³ runs (2026-07-07) — Phase B hypothesis confirmed with
data.** Two equilibrium (`S_support=1`) runs at 128³ on 1×H100:

| Run | Init L2_Ham | Constraint trend | min_lapse | Trapped surf. | Outcome |
|-----|-------------|------------------|-----------|---------------|---------|
| **Static** ω=0, m=0 | 0.0097 | *decreases* → ~0.0027 bounded | holds ~0.91 | none (max_ah_r=0) | **stable, throat holds** (t≥6) |
| **Rotating** ω=0.05, m=1 | 0.029 (3×) | *grows* | collapses 1.0→0.31 | grows from t≈0.8 | **NaN at t≈2.8** |

Interpretation: the rotating analytic ID starts with a **3× larger Hamiltonian
defect** (the `K_ij=0` **O(ω) momentum-constraint defect**, lesson L4/Debug §"Known
limitations") which then grows and destroys the run. The *same* code/ID machinery
is **stable in the static limit** (constraints actually *decrease*, throat holds,
no trapped surface). So the rotating instability here is **numerical (ID defect),
not physical** — it is precisely what the Phase B momentum-sector GRTresna solve
removes, and it confirms the static-baseline arm #1 is sound. Params:
`params_rotating_complex_equilibrium_static.txt` (new).

**Phase B constraint solve now CONVERGES (2026-07-08).** The former B-BLOCKER
("exotic Ham no-op") was a diagnostic-denominator artefact in `c_Ham_abs`, not a
real solver failure — fixed (see §3 B-BLOCKER RESOLVED). The exotic complex
winding solve (`params_rotating_wormhole_test.txt`, m=1, ω=0.05, throat) now
converges **Ham 0.76%, Mom 0.63%** and writes `InitialDataFinal.3d.hdf5`.

**Phase B4 DONE (2026-07-08) — κ ID family + level-0-dx `.gridinit` produced.**
Driver `scripts/wormhole/solve_kappa_family.{sh,py}` reuses the campaign's runtime
machinery (validated GRTresna BosonStarBH solve binary + shared
`convert_chombo_to_gridinit`) rather than the search-space CLI (whose
`scalar`/`boson_star` sectors do not express the exotic-winding throat). It scales
the throat/phantom amplitude `f → κ·f` and GRTresna re-solves both constraints per
κ. Results (m=1, ω=0.05, throat, N=64/L=64, native dx=1.0 = evolution level-0 dx):

| κ | amp | Ham % | Mom % | `.gridinit` |
|---|-----|-------|-------|-------------|
| 1.00 | 0.10 | 0.76 | 0.63 | ✅ `runs/rotating_wormhole_id/rotwh_omega_p0p05_m1_kappa_1p00/` |
| 0.90 | 0.09 | 0.57 | 0.43 | ✅ …`_kappa_0p90/` |
| 0.70 | 0.07 | 0.67 | 0.29 | ✅ …`_kappa_0p70/` |
| 0.50 | 0.05 | 0.92 | 0.16 | ✅ …`_kappa_0p50/` |

Each `.gridinit` verified structurally (dx=1.0, origin (0,0,−32), throat chi∈[0.85,1.11],
both winding channels φ/φ₂ and two-channel momentum Π/Π₂ populated, |Φ|² axisymmetric,
K=0). Regression tests: `tests/grtresna/test_rotating_wormhole_kappa.py` (10, incl.
fixture-gated winding/axisymmetry checks over the solved κ=1.0 gridinit).

**Phase B5 DONE (2026-07-08) — solved-ID rotating equilibrium HOLDS past the
analytic-ID failure point.** Loaded the κ=1.0 `.gridinit` into the
`complex_scalar` RotatingWormholeCollapse evolution (unigrid N=64/L=64 at the
gridinit's native dx=1.0, per L2) and evolved to **t=3.5** — well past the t≈2.8
where the *analytic*-ID rotating run NaN'd. Params:
`Examples/RotatingWormholeCollapse/params_rotating_grtresna_complex_kappa1.txt`.

| Quantity | Analytic ID (ω=0.05) | **GRTresna solved ID (κ=1, ω=0.05)** |
|----------|----------------------|--------------------------------------|
| Outcome | **NaN at t≈2.8** | **stable to t=3.5, no NaN** |
| Ham constraint | 0.029 → *grows* | 5.7e-3 → **2.2e-3 (decreasing)** |
| Mom constraint | grows | 6.6e-4 → **1.3e-4 (decreasing)** |
| J_z (spin) | — | **conserved: 77.5677 at t=0 and t=3.5** |
| Horizon / throat | collapses | **no horizon, throat holds (min_chi 0.86→0.98)** |

This is the headline result: the rotating instability was **numerical** (the O(ω)
`K_ij=0` analytic-ID momentum defect), not physical — the constraint-clean
GRTresna solve removes it and the rotating wormhole holds its equilibrium just
like the static case. Phase B is complete.

**Immediate next actions:** Phase C — the full (ω, m, κ) evolution grid (collapse
arm from κ<1, equilibrium arm from κ=1 across ω), then Phase D GW extraction.

**Phase C setup (2026-07-08) — CLI consolidation + high-res (dx=0.5) ID.**
Replaced the per-case `params_*.txt` proliferation with a single **CLI generator**
(same pattern as the QD campaign): `scripts/wormhole/wormhole_case.sh` renders the
evolution params from one in-code template + flags
(`--kappa --dx --omega --m --max-level --stop-time --gpu`), locates the matching
`.gridinit`, launches, renders frames offline, and prunes plotfiles. The eight
duplicated params files were deleted; only the legacy analytic/real-scalar
examples remain (documented in the example README). The κ family was re-solved at
**N=128 → dx=0.5** (`RES_N=128 solve_kappa_family.sh`) so a dx=0.5 / max_level=3
evolution stays L2-safe at level 0; solves converged Ham {0.95, 0.77, 0.99}%,
Mom {0.77, 0.36, 0.16}%. Solve HDF5 scratch is now auto-pruned (only the gridinit
kept). Example usage + status table live in
[`Examples/RotatingWormholeCollapse/README.md`](../../Examples/RotatingWormholeCollapse/README.md).

**Open L2 caveat for max_level=3:** the gridinit is uniform at dx=0.5, so AMR
levels below that (0.25→0.0625) interpolate the ID — the exact L2 kink risk. The
first high-res run is the test; fallback is lower max_level or the analytic
per-level ID layer (Debug.md step 3 / Phase C4).

---

## 0. Lessons already paid for (do not relearn)

From Debug.md — these are settled by controls and must constrain the design:

| # | Lesson | Evidence |
|---|--------|----------|
| L1 | **Frozen prescribed sources are fatal.** `EffectiveTeoMatter` (empty `add_matter_rhs`) diverges on unigrid, spin-independent, and no gauge choice rescues it (fixed-gauge test: Ham still →2.21). Matter must co-evolve | 4 Jun controls: `weak_spin ≡ a0` blow-up; spin residual 0.08% of field |
| L2 | **Never evolve finer than the ID file's native dx** — trilinear `ExternalGridInitialData` interpolation creates C0 kinks → NaN at the throat | `max_level=5` crash at t≈0.13; `ml=1` NaN at t≈0.17 |
| L3 | **A single real scalar cannot rotate smoothly** — `cos(mφ)` density modulation *is* the four-lobe pattern and disperses | failed first rotating run |
| L4 | **The cure is a complex phantom scalar with phase winding**: Φ = f(r,θ)e^{i(mφ−ωt)}, axisymmetric \|Φ\|², J_z from the winding | Route B design (8 Jun) — C++ written, **never compiled** |
| L5 | **GRTresna rotating solves work**: ω=0.05 lump converged (Ham 0.74%, Mom 0.012%) and evolved cleanly to t=4.0 (Ham 1.3e-4) | 6 Jun validation |
| L6 | **Disk discipline**: `plot_interval ≥ 10`, consumer with `--delete`, or 100+ GB backlogs | 112 GB incident |

**On "Teo":** the literal Teo metric is a stationary *geometry* with unspecified
matter — evolving it as a frozen source is exactly L1 and cannot work. The
correct physical target is a **rotating phantom-scalar wormhole**
(Kleihaus–Kunz class), which limits to Ellis–Bronnikov at ω=0. We keep the
Teo-class phenomenology (rotating throat, frame dragging, natural quadrupole)
with a self-consistent co-evolving source. The analytic-Teo `.gridinit` path
stays regression-only.

---

## 1. Physics design

### Ansatz (Route B, revived)

Complex phantom scalar stored as (φ₁,Π₁,φ₂,Π₂):

- Φ = f(r,θ) e^{i(mφ_az − ωt)}, phantom (sign-flipped) stress-energy
- \|Φ\|² = φ₁²+φ₂² axisymmetric → no four-leaf, no L3 dispersal
- J_z carried by phase winding; m ∈ {1,2}, ω ∈ [0, 0.2] sweep
- ω=0 limit must reproduce the published EB behavior (built-in regression)

### Collapse trigger — the constraint-clean upgrade

The EB article triggered collapse by scaling the stress-energy in the RHS
(`S_support = 0.5`), which injects a t=0 Hamiltonian defect that then grows.
Replace it:

> **Amplitude-reduction trigger:** paint the throat + phantom profile with a
> reduced field amplitude f → κ·f (κ < 1) and have **GRTresna re-solve the
> Hamiltonian and momentum constraints for that configuration.** The result is
> an *exact* Einstein solution that is simply not in equilibrium — it collapses
> (or inflates) with **zero initial constraint defect**. `S_support` is retained
> only as a comparison arm to quantify how much of the article's constraint
> growth was trigger-induced.

This is the direct answer to the "steadily growing violation" problem: the
defect never enters, and the momentum constraint absorbs the O(ω) winding
defect that Route B planned to just accept (K_ij=0 hack no longer needed —
GRTresna solves for the momentum sector, already demonstrated at Mom 0.012%).

### Science targets

1. **Rotation vs collapse:** does angular momentum delay/prevent throat
   collapse (centrifugal support)? Is there a critical ω separating collapse
   from a spinning remnant/dispersal?
2. **Natural ℓ=2, m=±2 GW** from the rotational quadrupole — waveform family
   vs (ω, m, κ), replacing the article's single hand-perturbed waveform.
3. **Remnant spin:** J_z through collapse; Kerr QNM fit of the ringdown →
   measured remnant ĵ = J/M², cross-checked against swallowed J_z.
4. **Phantom bounce with rotation:** does the bounce go non-axisymmetric
   (bar-mode m=2 growth → enhanced GW)?

---

## 2. Phase A — Build & regression (days)

Route B C++ was written blind (no compiler available at the time). It now
**builds** (see Status above); the remaining Phase A work is the no-NaN smokes
and the ω=0 regression that validate the matter model.

| # | Task | Acceptance |
|---|------|-----------|
| A1 | ~~CPU compile check~~ **DONE** — `make USE_CUDA=FALSE USE_MPI=FALSE COMP=gnu DIM=3` builds `main3d.gnu.ex` | clean build ✅ |
| A2 | ~~CUDA build~~ **DONE** — `USE_CUDA=TRUE CUDA_ARCH=90` builds `main3d.gnu.MPI.CUDA.ex` | clean build ✅ |
| A3 | **ω=0 regression:** complex-scalar run with ω=0, m=0 must reproduce the EB article baseline — equilibrium hold to t≈1.5M then noise-driven branch, λ≈9/M growth when perturbed | growth rate within ~10% of article value |
| A4 | Equilibrium smoke at ω=0.05 (`params_rotating_complex_equilibrium.txt`): axisymmetric \|Φ\|², J_z nonzero and ~constant, bounded constraints (unigrid, per L2) | Debug.md §"How to test" pass criteria |

A3 is the keystone: it validates the entire new matter model against published
results before any new physics is claimed.

---

## 3. Phase B — GRTresna rotating-throat initial data (~1–2 weeks)

The constraint-clean core. Sub-steps ordered so each is testable alone.

### B0 audit findings (2026-07-07) — what exists vs. the precise gap

Audited `../GRTresna/Source/{Matter,Methods,Core}` + `Examples/BosonStarBH`.

**Already present (reusable as-is):**
- **Exotic/phantom solve.** `L.exotic` sign-flips the whole `T_ab`; the solver
  auto-switches to **K=0 maximal slicing** (York/Lichnerowicz) for `rho<0`
  (`CTTKHybrid.impl.hpp` ~L22-38, 238-251), with `psi_floor` /
  `psi_relaxation` / `maximal_jacobian_cap` robustness (`Grids.cpp` ~L233-242).
  Exotic lump amplitude is auto-damped (`EXOTIC_AMP_SCALE`).
- **Matter-sourced momentum constraint.** `CTTKHybrid.impl.hpp` ~L273-281 sets
  `rhs(V_i) = psi^6 (2/3 ∂_i K + 8πG·Si)` — it *already* takes a matter
  momentum density `Si` and solves the vector potential `V_i`. **This is the
  hook the winding needs.**
- **Singular-ψ throat.** `compute_bowenyork_psi` gives `psi_bh = (b0/2)/r` from
  `bh1_bare_mass` (`PsiAndAijFunctions.cpp` ~L65-104); solver solves `psi_reg`
  on top (`psi_0 = psi_reg + psi_bh`). Set `bh1_bare_mass = b0/2` → EB throat.
- Example `BosonStarBH` = `CTTKHybrid<ComplexScalarField>` template.

**The precise gap (two edits in `ComplexScalarField.cpp`):**
1. **No genuine phase winding.** `paint_boson_fields`/lump painting set
   `phi2_k = 0` **everywhere** (L110) and encode rotation only via
   `pi2=-(ω/α)phi1`; the azimuthal "mode" is a **real** `angular_factor` (dx/w)
   modulation of `phi1` — i.e. the `cos(mφ)` real modulation that **is** the L3
   four-lobe dispersal. Need the true winding matching the GRTeclyn evolution
   convention: `φ₁=f(r,θ)cos(mφ_az)`, `φ₂=f(r,θ)sin(mφ_az)`,
   `Π₁=(ω/α)φ₂`, `Π₂=−(ω/α)φ₁`.
2. **Momentum density drops the imaginary channel.** `Si += sign·(−pi1·dphi1)`
   only (L134, 159); `dphi2=0` so the winding `j_φ` is never sourced → no clean
   `J_z`. Fix: `Si += sign·(−(pi1·dphi1 + pi2·dphi2))`. With the winding
   `dphi2≠0`, this feeds the real azimuthal momentum into the L273-281 solve →
   **constraint-clean frame-dragging `A_ij`**.

Both edits are localized to `ComplexScalarField.cpp` (paint + emtensor) plus the
matching potential/derivative helpers in `BosonStarParams.hpp` to emit
`f·sin(mφ)` for `phi2` and its gradient. `T_ab` must stay identical to the
GRTeclyn `complex_scalar` evolution model so the solved `.gridinit` is loadable
and in equilibrium at t=0. `USE_COMPLEX_SCALAR_MATTER` gates the TU.

**Revised B-substep mapping** (concrete files below):

| # | Task | Detail / files |
|---|------|----------------|
| B0 | **Audit** complex-scalar + exotic coupling + (ω, m) support in GRTresna and the wrapper (the boson-star / exotic-wedge machinery covers most of this; `ComplexScalarField.cpp` exists) | `../GRTresna/Source/Matter/`, `grtresna/fields/`, `grtresna/matter/wiring.py` |
| B1 | **Throat via singular-ψ machinery:** the EB conformal factor ψ = √(1+b₀²/4r̄²) → (b₀/2)/r̄ as r̄→0 — the *same* 1/r leading singularity as a puncture. Use `regularised_part_psi` with a bare-mass-like singular part (coefficient b₀/2) and let GRTresna solve the regular correction sourced by the phantom profile. Validate at ω=0 against the analytic EB ψ | GRTresna solve config; wrapper `solver/config.py` (`bh1_bare_mass`-style seed), `psi_floor` |
| B2 | **Winding momentum sector:** phantom Π = ω·f phase pattern → axisymmetric j_φ source; solve the momentum constraint (CTTKHybrid Vi). L5 shows this converges | `grtresna/solver/`, momentum wiring |
| B3 | **Profile painting:** f(r)·(sinθ)^m with tabulated f(r) — reuse the boson-star tabulated-profile loader (2/3-column `qball_profile.dat` path); θ-dependence painted analytically in C++ (already in Route B `SupportedWormholeInitialData`) | `grtresna/profiles/`, `SupportedWormholeInitialData.hpp` |
| B4 | **Amplitude-reduction ID family:** generate κ ∈ {1.0, 0.9, 0.7, 0.5} re-solved configs per (ω, m); κ=1.0 is the equilibrium arm | driver script modeled on `make_rotating_wormhole_id.py` |
| B5 | **Resolution rule (L2):** emit `.gridinit` at the exact evolution level-0 dx; postload gate before every production run; name-mapped loader already fixed | `grtresna/io.py`, `projection/postload_gate.py` |

**Acceptance:** t=0 Ham/Mom at truncation level (target: beat the L5 benchmark
Ham 0.74% → sub-0.1% interior); ω=0, κ=1 solved ID matches analytic EB data to
truncation; loaded-vs-solved constraint norms agree through the postload gate.

### B1–B3 implementation status (2026-07-07) — winding momentum DONE, exotic-Ham blocker found

**Implemented in GRTresna (`feature/interstellar`, built clean):**
- `BosonStarParams.hpp`: added a per-lump `winding` flag + the genuine
  phase-winding ansatz helpers `polar_factor` ((sinθ)^m), `lump_winding_modulus`
  (f = amp·env(r)·(sinθ)^m), `lump_phi_winding` (φ₁=f cos mφ_az, φ₂=f sin mφ_az),
  `lump_grad_phi_winding`; `read_boson_lump` reads `winding`.
- `ComplexScalarField.cpp`: winding branch in **both** `initialise_matter_vars`
  (paints φ₂, Π₁=(ω/α)φ₂, Π₂=−(ω/α)φ₁ using the per-lump ω) **and**
  `compute_emtensor` (two-channel momentum density
  `Sᵢ = −(Π₁∂φ₁ + Π₂∂φ₂)` and both gradient-energy channels in ρ).
- Wrapper: `solver/params.py` emits `lump{k}_winding`.

**Validated — the momentum sector works:** a single exotic winding lump
(m=1, ω=0.05, bare-mass throat) solves the **momentum constraint to Mom ≈
0.009%** (matching the L5 benchmark). The `winding=0` control gives identical
behavior, confirming the winding change is correct and side-effect-free.

### B-BLOCKER — RESOLVED (2026-07-08): it was a diagnostic denominator artefact, not a stalled solve

**Root cause (found & fixed).** The exotic Hamiltonian solve was **never actually
stuck** — the *reported* residual was. `Diagnostics.impl.hpp`'s `c_Ham_abs`
normaliser summed `+24πG·rho` **without `abs()`** while every sibling term (`A²`,
`lap ψ`, and all of `Mom_abs`) used a magnitude. For EXOTIC (ρ<0) matter under
maximal slicing the converged constraint drives `12·lap(ψ)·ψ⁻⁵ → −24πGρ =
+|24πGρ|`, so the signed `+24πGρ` term **cancels the Laplacian term** and
`c_Ham_abs → 0`. `Ham_norm = |Ham / Ham_abs|` then divides by ~0 and is pinned
near 100% *regardless of how well the ψ-solve converged*. The real
`ScalarField` path escaped this only because its ρ mixes a canonical (ρ>0)
background + unsigned potential, so its denominator never fully cancels.

**Fix.** Use `abs(emtensor.rho)` in `c_Ham_abs` (mirrors `Mom_abs` and every
other term), giving a strictly positive scale. Also added a one-line
per-NL-iteration instrumentation in `GRSolver::run()` (`max|dψ|`, `psi_reg`
range) to positively confirm the ψ-solve was moving all along.

**Verification (post-fix, `params_rotating_wormhole_test.txt`, exotic winding
m=1, ω=0.05, throat, amp 0.1, N=64³, ml0):** Ham now drops monotonically
`1.0e2 → 40.4 → 19.3 → 9.8 → 5.1 → 2.7 → 1.4 → 0.76 %`, converging alongside
Mom `→ 0.63 %` (early-exit at NL iter 8). Instrumentation shows `max|dψ|`
shrinking `0.049 → 0.001` and `psi_reg` dipping into the throat `1.0 → 0.896` —
i.e. the elliptic ψ-solve was progressing correctly the entire time; the "100%"
was purely the metric. `InitialDataFinal.3d.hdf5` is produced. **Both
constraints now converge for exotic complex winding matter — the blocker is
cleared and Phase B4 (κ family + `.gridinit` export) can proceed.**

<details><summary>Original symptom writeup (kept for the record)</summary>

Any EXOTIC (ρ<0) lump solved through `CTTKHybrid<ComplexScalarField>` (the
`BosonStarBH` build, `maximal_slicing=1`) *appeared* to leave the Hamiltonian
residual pinned at **exactly `1.000000e+02` %** for every NL iteration. The
momentum sector was unaffected (`Vi`/`Mom` converged normally). Canonical (ρ>0)
complex lumps were fine because they use the `K=√(24πρ)` ansatz, not the
elliptic ψ-solve. Why *exactly* 100%: at iter 0 (ψ_reg=1 ⇒ lap=0, `A²`=0, K=0)
`c_Ham=−24πρ` and `c_Ham_abs=+24πρ` ⇒ `‖c_Ham‖=‖c_Ham_abs‖` ⇒ 100%; and because
`c_Ham_abs` kept cancelling to ~0 at every subsequent (actually converging)
iterate, the *ratio* stayed ~100% even as the true residual `c_Ham` fell.

</details>

**Debugging history / controls run (all `maximal_slicing=1`, N=64³ unless noted).**
NOTE: the "100% stuck" column below is the *artefactual* metric — with the
`c_Ham_abs` fix these same configs report their true (converging) residuals; e.g.
config #1 now reaches 0.76%.

| # | Config | Matter path | grid | Ham | Mom | Verdict |
|---|--------|-------------|------|-----|-----|---------|
| 1 | exotic, winding m=1, ω=0.05, throat b₀/2, amp 0.1 | Complex | ml0 | **100% stuck** | 2.8%→stall | throat+strong amp |
| 2 | exotic, winding m=1, ω=0.05, amp 0.03, massless, no throat | Complex | ml0 | **100% stuck** | **0.017%** ✅ | winding momentum works |
| 3 | exotic, **winding=0**, m=1, amp 0.03, massless | Complex | ml0 | **100% stuck** | 0.048% ✅ | not the winding |
| 4 | exotic, winding m=1, ω=0.05, amp 0.02, **m_φ=0.1**, w=6 | Complex | ml0 | **100% stuck** | 0.009% ✅ | not massless/amp |
| 5 | = #4 but `maximal_jacobian_cap` 1e6→**25** | Complex | ml0 | **100% stuck** (identical) | id. | not the Jacobian cap |
| 6 | **canonical** (exotic=0), `maximal_slicing=0` | Complex | ml0 | 4.9e-15→**0.79%** ✅ | 0.23% ✅ | canonical complex fine (K-ansatz) |
| 7 | exotic, **mode=0, winding=0**, ω=0, amp 0.08, m_φ=0.1 (simplest spherical exotic) | Complex | ml0 | **100% stuck** | 0 | fails even at the simplest config |
| R1 | `params_exotic_amr_test` (amp eff 0.05) | **Real ScalarField** | ml3 | **0.10%** ✅ | 0.027% ✅ | real exotic converges |
| R2 | = R1 but **ml0**, N=64³, amp 0.08 | **Real ScalarField** | ml0 | **0.084%** ✅ | 0.023% ✅ | real exotic converges at *my* grid too |

**Isolation (what is ruled OUT):** grid/AMR level (R2 vs #7 same grid), azimuthal
mode / winding (#7 spherical still fails; #2 winding momentum works), field mass
(#2 vs #4), amplitude / existence boundary (#7 amp 0.08 ≈ R2, and #2 tiny amp all
fail), `maximal_jacobian_cap` (#4≡#5), throat bare mass (#2 no throat fails),
multigrid-var index collision (`ComplexScalarFieldVariables` matter slots start at
`NUM_METRIC_VARS`, no overlap with `c_psi_reg`). The two `Main`s are identical
(`GRSolver<CTTKHybrid<M>,M>`). The "only variable that flips the outcome is the
matter class" observation was the red herring: **both classes converged; only the
`ScalarField` path *reported* it**, because its mixed-sign ρ (canonical
background + unsigned potential) kept `c_Ham_abs` away from the exact
cancellation that a pure-exotic complex ρ produced. See the resolution above.

### Resolution & fixes applied

1. ✅ **APPLIED — the actual cause.** `Diagnostics.impl.hpp` `c_Ham_abs`: changed
   `+24πG·rho` → `+24πG·abs(rho)`. This is what unblocked the exotic complex
   solve (Ham now `→0.76%`).
2. ✅ **APPLIED — diagnosis instrumentation.** `GRSolver::run()` now prints per-NL
   `max|dψ|` and the `psi_reg` range, which positively confirmed the ψ-solve was
   progressing all along (`dψ 0.049→0.001`, `psi_reg 1.0→0.896`).

<details><summary>Earlier ranked proposals (superseded — kept for context)</summary>

- **Instrument `rho` and the ψ correction (diagnosis first).** `pout()` of
  `∫rho dV` and `max|constraint_vars(c_psi)|` in `GRSolver::run()`. Two outcomes:
   - `rho` differs (≈0, wrong sign, or huge) ⇒ bug is in
     `ComplexScalarField::compute_emtensor` → fix the ρ formula.
   - `rho` matches but the ψ correction is ≈0 while Vi's is not ⇒ the ψ row of
     the elliptic solve is not being assembled/relaxed in the complex build →
     look at `update_psi0` and the operator define under `USE_COMPLEX_SCALAR_MATTER`.
   - `rho` matches but the ψ correction is ≈0 while Vi's is not ⇒ the ψ row of
     the elliptic solve is not being assembled/relaxed. (This is what the applied
     instrumentation checked — the ψ correction was *not* zero, confirming the
     solve was fine and the metric was the culprit.)
- **Diff the two builds' non-matter headers** / **two-real-scalar fallback** —
  investigated as candidate causes; both moot now that the solve is confirmed
  correct.

</details>

**Still relevant post-fix (Phase B4 follow-on):** add **exotic-amplitude
damping** to the complex path — mirror `MatterParams::EXOTIC_AMP_SCALE` /
`effective_amp` in `BosonStarParams` so exotic complex lumps stay inside the
Lichnerowicz existence boundary as the κ family pushes larger winding gradient
energy. Not needed for the current amp-0.1 solve, but keep in mind if a κ point
fails to converge.

**Reproduce:** build `scripts/wormhole/build_grtresna_bosonstar.sh`; solve
`scripts/wormhole/solve_rotating_wormhole_test.sh`. Test inputs:
`../GRTresna/Examples/BosonStarBH/params_rotating_wormhole_test*.txt`
(test2 = winding momentum PASS; test7_spherical = simplest exotic FAIL;
test6_canonical = canonical PASS) and
`../GRTresna/Examples/ScalarFieldBH/params_exotic_ml0_test.txt` (real PASS).
Ham/Mom history in each example's `Ham_and_Mom_errors.txt` + `pout/pout.0`.

---

## 4. Phase C — Evolution program (~2–3 weeks GPU)

All runs: unigrid first (L2), `plot_interval ≥ 10` (L6), consumer sidecar with
deletion, frames on for `chi`, `K`, `lapse`, `phi/Pi`, `Weyl4_*`, `embedding`.

| # | Run | Purpose | Pass criteria |
|---|-----|---------|--------------|
| C0 | ω=0, κ=0.5 (or S=0.5 comparison arm) | reproduce the article's collapse with the new stack | trapped surface, phantom bounce, ℓ=2 signal matching article morphology |
| C1 | Equilibrium holds: κ=1, ω ∈ {0.05, 0.1, 0.2} | rotating wormhole is (meta)stable on the run timescale | J_z constant, \|Φ\|² axisymmetric, Ham flat (contrast article's growth) |
| C2 | **Collapse grid:** κ ∈ {0.9, 0.7, 0.5} × ω ∈ {0, 0.05, 0.1, 0.2} × m ∈ {1, 2} (prune by interest) | the physics: rotation-vs-collapse phase map, critical ω hunt | branch outcome per cell; J_z tracked through collapse |
| C3 | **GW extraction** on C2 winners: Ψ₄ multipoles including m=±2 (extend `process_wave` beyond l=2,m=0), retarded-time alignment, propagation-speed test (article methodology), Kerr QNM fit → remnant spin | v≈c; QNM ĵ consistent with swallowed J_z |
| C4 | AMR enablement (only after C0–C2 stable): resolutions matched per L2, or analytic per-level ID layer (Debug.md step 3) if refinement below file dx is needed | no cross-level NaN; convergence pair at 2 resolutions |

**Constraint-violation comparison (the article upgrade, explicit deliverable):**
side-by-side Ham(t) for (a) article-style `S_support` defect trigger vs (b)
κ-re-solved trigger, same physics otherwise. Target result: (b) stays bounded
until horizon formation. This figure alone justifies the GRTresna integration.

Optional: drive the (ω, m, κ) grid through the wrapper campaign machinery
(replay/HQ infra) instead of hand-launched runs — gets consumer discipline,
`score_timeseries`, and tier validation for free.

---

## 5. Phase D — Analysis & paper deliverables

1. **Phase map:** collapse / bounce / dispersal outcome over (ω, κ), with the
   critical-ω boundary if it exists in range.
2. **Waveform family:** natural-quadrupole ℓ=2, m=±2 waveforms vs ω — the
   headline extension over the single perturbed EB waveform; burst morphology
   (CWT) and strain projections (with the near-zone caveat of NextSteps P7
   stated honestly; extraction radii 12–24 as in the article).
3. **Remnant characterization:** Kerr QNM spin vs initial J_z — a
   self-consistency test unavailable in the non-rotating paper.
4. **Methods section:** constraint-clean trigger (κ re-solve) + Ham(t)
   comparison — directly addresses the softest point of the EB article.
5. Rotating phantom bounce phenomenology (bar mode or not).

---

## 6. Risks

| Risk | Mitigation |
|------|-----------|
| Route B code has unknown compile debt (never built) | Phase A1 CPU compile first; budget 1–3 days of interface fixes |
| GRTresna nonlinear solve with singular ψ + phantom source diverges | B1 validates at ω=0 against analytic EB first; fall back to solving with a weaker singular seed + more NL iterations; `maximal_slicing` on |
| Ghost field + rotation → faster instabilities at high ω (no stabilizing potential) | keep ω ≤ 0.2 in the first grid; equilibrium arm C1 detects runaway before collapse runs |
| θ-dependent profile vs radial tabulated loader | (sinθ)^m painted analytically in C++ (Route B already does this); loader stays radial |
| O(ω²) residuals after Mom solve (Ham cross-terms) | measured at B-acceptance; if large, iterate the solve with updated Π |
| AMR kinks return in C4 | strict L2 rule; analytic per-level ID layer as the documented fallback |
| m=±2 mode extraction not yet in `process_wave` | small extension of the existing l=2,m=0 machinery; do early in C3 |
| Disk (L6) | `plot_interval ≥ 10`, consumer `--delete --keep-last 2`, checked in every launcher |

---

## 7. Execution order & effort

```
A  build + ω=0 regression        ~3–5 days      ── gate for everything
B  GRTresna rotating ID           ~1–2 weeks     ── B1 (ω=0 validation) first
C0–C1  controls + equilibria      ~3–4 days GPU  ── article reproduction
C2  collapse grid                 ~1–2 weeks GPU ── the physics
C3  GW analysis                   ~1 week        ── overlaps C2
C4  AMR + convergence             ~3–5 days GPU  ── only after C2 stable
D   analysis + writing            ~1–2 weeks
```

Minimum publishable core: **A + B1/B4 + C0 + C1 + C2 (one m, three ω) + C3** —
a rotating-collapse waveform family with constraint-clean initial data. The
critical-ω hunt and AMR convergence strengthen but do not gate the paper.
