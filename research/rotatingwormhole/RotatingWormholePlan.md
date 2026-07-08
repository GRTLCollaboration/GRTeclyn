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

**Immediate next actions (updated 2026-07-08 after the boundary study, see
§"Boundary study DONE" below):** the t≈5 dispersal is now understood — it is a
massless-ghost instability whose *removal timing* is boundary-set. The gating
priority is therefore the **massive-field re-solve** (`phantom_mass > 0`, which
confines the phantom cloud into a bound soliton) so the equilibrium arm actually
holds, *then* the full (ω, m, κ) grid and Phase D GW extraction.

**Phase C setup (2026-07-08) — CLI consolidation + high-res (dx=0.5) ID.**
Replaced the per-case `params_*.txt` proliferation with a single **CLI generator**
(same pattern as the QD campaign): `scripts/wormhole/wormhole_case.sh` renders the
evolution params from one in-code template + flags
(`--kappa --dx --omega --m --max-level --stop-time --gpu`), locates the matching
`.gridinit`, launches, streams frames + Ψ₄ and deletes plotfiles *during* the run
(a `consume_plotfiles --watch --delete` sidecar, L6), then drains the tail. The eight
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

**Phase C high-res DONE (2026-07-08) — dx=0.5, max_level=3 AMR is stable; L2
caveat resolved.** κ ∈ {1.0, 0.7, 0.5} (m=1, ω=0.05) each evolved to **t=8 on a
4-level AMR grid (dx=0.5→0.0625) with no NaN and no cross-level kinks** — the L2
interpolation risk did **not** materialise once `regrid_interval` was given one
entry per level (AMReX `ParmParse` aborts otherwise; the low-res runs were unigrid
so it never triggered). Constraints stay **bounded and decreasing** for every κ:

| κ (dx=0.5, ml=3, t=0→8) | L2_Ham | L2_Mom | min_chi (throat) | horizon `max_ah_r` / `min_θ₊` | matter `rho_sum` | J_z |
|---|---|---|---|---|---|---|
| 1.0 | 8.4e-3 → **2.0e-3** ↓ | 1.5e-4 → 1.1e-4 | 0.57 → 0.91 | **0 / >0 (none)** | 774 → 1.5 | −0.031 → ~0 |
| 0.7 | 6.5e-3 → **2.2e-3** ↓ | 1.1e-4 → 1.2e-4 | 0.45 → 0.79 | **0 / >0 (none)** | 375 → 2.1 | −0.018 → ~0 |
| 0.5 | 5.8e-3 → **2.3e-3** ↓ | 1.0e-4 → 1.3e-4 | 0.41 → 0.73 | **0 / >0 (none)** | 191 → 1.3 | −0.009 → ~0 |

**Headline (and a caveat that changes the low-res narrative).** At high resolution
with constraint-clean GRTresna ID, **none of the κ arms collapse to a black hole**
(`max_ah_r=0`, `min_θ₊>0` throughout). Instead all three do the *same* thing: the
throat opens toward flat (min_chi → ~1), and the supporting **phantom cloud
disperses at t≈4.5–5 — at essentially the same time for every κ** (`rho_sum` and
`J_z` both decay to ~0). The Ψ₄ ℓ=2,m=0 signal grows monotonically
(Re @ R=12: −4e-4 → −0.029), i.e. the rotational quadrupole *is* radiating. This
**supersedes the low-res "monotonic collapse trend with decreasing κ"**: at dx=1.0
the arms looked like they were marching toward collapse, but at dx=0.5/ml=3 they
disperse without forming a horizon, and the outcome is κ-independent.

**Why the caution:** the κ-independence *and* the simultaneity of the t≈5 dispersal
are the signature of a **numerical/boundary effect, not κ-driven collapse physics**.
The measured outer-boundary flux (`boundary_flux.dat`) stays small (~1–6e-5) with no
spike at t≈5, so the matter is not obviously leaving through the Sommerfeld faces —
it looks like the phantom cloud *spreads* (amplitude drops as it delocalises), and
`rho_sum` (a refined-region volume sum) is also sensitive to AMR de-refinement as
the cloud disperses. This must be pinned down before any physical claim.

Diagnostics + movies (README scripts): per-run `plots/constraints_plot.*`,
`plots/psi4_*.*` (`plot_diagnostic.sh` + `plot_extracted_psi4.py`) and
`movies/movie_{chi,K,lapse,phi,Pi,Weyl4_*}_z.mp4` (`make_movies.sh`); κ-family
comparison in `runs/rotating_wormhole/kappa_family_{diagnostics,constraints}.png`.
Note the shared collapse-diagnostics plotter caps at ≤14 columns while the
RotatingWormhole diag writes 19 (min_lapse…rho_sum, J_z) — hence the custom
comparison figure.

**Pipeline hardening (this pass):** `wormhole_case.py` now (i) emits
`regrid_interval` with one value per level, (ii) runs the plotfile consumer as a
**live sidecar** (frames+Ψ₄ rendered and plotfiles deleted *during* the run, only
the newest few kept), (iii) drains the tail with `--keep-existing-frames` (the
consumer wipes frames at startup by default — a post-run drain must not), and
(iv) defaults `--frame-jobs 4` (`-j`≫4 on NFS thrashes and hangs). Result: each
finished run is ~15 MB (frames + `.dat` only); **zero** plotfile backlog.

**Next steps (Phase C follow-up + D):**
1. **Resolution/boundary study on the t≈5 dispersal (gating).** Re-run κ=1.0 at
   dx=1.0 (N=64) and at dx=0.5 with a **larger box** (e.g. L=96/128, boundary
   farther out); if the dispersal time moves with resolution/boundary it is
   numerical, if invariant it is physical. Also add a `--central-timeseries` /
   `rho` volume integral that is AMR-robust to confirm whether matter spreads vs.
   leaves.
2. **ω=0 static baseline with the *same* CLI** (arm #1 apples-to-apples): does the
   non-rotating complex-scalar throat also disperse at t≈5, or hold? This
   discriminates a rotation-specific instability from a generic ID/boundary issue.
3. **Phase D GW analysis** once (1)–(2) settle the interpretation: Ψ₄ ℓ=2 m=0/±2
   multipoles, retarded-time alignment, strain, spectra; QNM/remnant-spin only if
   a horizon ever forms (none did here).

---

## Boundary study DONE (2026-07-08) — dispersal is BOTH physical (unstable ghost) AND boundary-timed

Ran the box-size convergence test from step 1 above: **κ=1.0, ω=0.05, m=1,
dx=0.5, ml=3, L=96** (boundary at r=48 from the throat vs r=32 for L=64),
evolved to **t=46** with the new `--box-size` flag. This overturns the
"purely numerical" caveat above with a sharper, two-part answer.

**Result — the dispersal time moves with the box.**

| Quantity | L=64 (r_bdy=32) | **L=96 (r_bdy=48)** |
|----------|-----------------|---------------------|
| `rho_sum` holds until | t≈4.0 | t≈4.5 |
| `rho_sum` half-gone | t≈4.7 | **t≈5.1** |
| `rho_sum` ≈ 0 (fully dispersed) | t≈5.0 | **t≈6.0** |
| Constraints (Ham/Mom) | bounded, ↓ | bounded, ↓ (1.0e-3 / 7.7e-5 at t=46) |
| Horizon | none | none (`max_ah_r=0`, `min_θ₊>0`) |
| min_chi (throat) | → 0.91 | → 0.91 (opens toward flat) |

The dispersal **time shifts later with the larger box** (fully gone at t≈6 vs
t≈5), and `boundary_flux.dat` now shows a **clear net-outward-flux peak at
t≈5.6** (3.9e-4) coincident with the `rho_sum` collapse — the earlier "no flux
spike" reading was the coarser L=64 diagnostic cadence. So the matter **does**
leave through the Sommerfeld faces, and *when* it fully leaves is boundary-set.

**Why the matter disperses (physics, not a bug).** Two compounding causes, both
expected:

1. **Phantom-scalar wormhole equilibria are dynamically unstable.** The
   Ellis–Bronnikov / massless-ghost throat has an unstable fundamental radial
   mode (Gonzalez–Guzmán–Sarbach 2009; Shinkai–Hayward 2002). GRTresna gives a
   genuinely constraint-satisfying slice, but it sits on an unstable saddle: any
   perturbation grows. Here the seed is a **gauge transient** — the ID is
   K=0/lapse=1 (maximal-slicing solve) loaded into a 1+log / Γ-driver evolution,
   so at t≈0.5 `max_K` spikes to ~2.2 and the lapse dips to ~0.73 as the gauge
   relaxes. That kick is enough to push the throat off the unstable equilibrium.
2. **No confining potential (`phantom_mass = 0`).** The ghost scalar is
   *massless*, so even absent the instability it has nothing to bind it — a
   massless field spreads as dispersive radiation. There is no `V=μ²|Φ|²`
   potential well (as a boson star / Kleihaus–Kunz *massive*-field wormhole would
   have) to hold the cloud together. Once perturbed it delocalises, its amplitude
   drops, the front propagates out at ~c, and the finite box then absorbs it.

The throat *geometry* is fine throughout (min_chi rises smoothly to 0.91,
constraints decrease); what fails is **matter confinement**, not the metric
solve. The rotating equilibrium doesn't collapse — it **evaporates its own
support** and relaxes toward flat space.

**What we can do about it (ranked).**

1. **Give the scalar a mass (`phantom_mass > 0`) — the primary fix.** A mass term
   creates a bound-state potential well; the phantom lump becomes a (meta)stable
   soliton (the massive-complex-field / Kleihaus–Kunz rotating wormhole class,
   which is the correct physical target anyway). Requires re-solving the κ family
   with the *same* nonzero mass so the ID stays in equilibrium, then evolving with
   `phantom_mass` matched. **Do this next** — it is the single largest lever and
   converts "dispersing lump" into "confined rotating throat."
2. **Self-interacting potential (Q-ball type),** `V = μ²|Φ|² − g|Φ|⁴ + …`, for a
   deeper, more robust well if a bare mass is marginal. Reuses the existing
   `qball_profile` machinery.
3. **Active support — the deferred RL matter pump / PD "trap" controller**
   (`RLMatterPumpParams.hpp`). This is exactly the mechanism to counter dispersal:
   feed matter back to hold the throat. Explicitly deferred in the approved scope;
   the boundary study is the quantitative motivation for turning it on.
4. **Reduce the seed kick:** load the ID lapse/shift consistently and/or ramp the
   gauge gently so the t≈0.5 `max_K` transient doesn't excite the unstable mode.
   Delays but does not cure the intrinsic instability (only #1–#3 confine matter).
5. **Sponge / higher-order absorbing outer BC:** removes boundary *reflection*
   only — it does **not** stop the physical spreading, so it will not prevent
   `rho_sum → 0`. Useful for clean GW extraction, not for confinement.
6. **Embrace it as the physics result:** a massless rotating phantom wormhole is
   unstable and radiates its support away; the transient ℓ=2 Ψ₄ burst as the
   throat relaxes *is* a publishable signature. Valid fallback if the massive-field
   branch proves hard to stabilise.

**Reproduce:**
```bash
EVO_L=96 RES_N=192 bash grteclyn-wrapper/scripts/wormhole/id/solve_kappa_family.sh 1.0 4
bash grteclyn-wrapper/scripts/wormhole/run/wormhole_case.sh --kappa 1.0 --dx 0.5 --box-size 96 --gpu 0
```
Run dir `runs/rotating_wormhole/evo_omega_p0p05_m1_kappa_1p00_dx0p5_ml3_L96/`
(84 MB final: frames + `.dat` only, sidecar-cleaned). Plots in `output/plots/`
(`constraints_plot`, `collapse_diagnostics_plot`, `psi4_analysis_M30_D10`).

**Revised next step:** the gating question is settled (dispersal = unstable
massless ghost + boundary-timed removal). Priority is now **the massive-field
re-solve (#1)** — solve the κ family and evolve with `phantom_mass > 0` — before
the (ω, m, κ) collapse grid, so the equilibrium arm actually holds.

---

## Massive-field test DONE (2026-07-08) — a bare mass does NOT confine; the profile must be a bound eigenstate

Implemented the confining-potential lever (#1 above): threaded a field mass μ
through both codes — `MASS` env → GRTresna `scalar_mass`, `--mass` flag →
GRTeclyn `phantom_mass` — with the shared convention `V = ½μ²|Φ|²` (verified:
`PhantomDecayPotential` == `ComplexScalarField::potential_value`, λ=μ_sextic=0).
Tag suffix `_mass<μ>` keeps massive/massless runs separate (backward compatible).
Solved κ=1.0, μ=0.1 (`ω=0.05 < μ`, bound regime; Ham 0.98% / Mom 0.81%) and
evolved at L=64, dx=0.5, ml=3.

**Result — the mass only *delays* dispersal, it does not confine.**

| t | massless `rho_sum` | **μ=0.1 `rho_sum`** |
|---|--------------------|---------------------|
| 0 | 633 | 633 |
| 4 | ~dispersing | 745 |
| 5 | ~1.5 | 457 |
| 6 | — | **4.9** |
| 8 | 1.5 | 1.7 |

The cloud still evaporates (`rho_sum → ~1` by t≈6–8), just ~1 M later than the
massless case. Constraints stay clean throughout (Ham 1.8e-3, Mom 6.5e-5); the
throat geometry is fine (min_chi → 0.93). So the mass is not the missing
ingredient by itself.

**Why a bare mass fails — the profile is not a bound state.** The painted lump
uses `lump0_profile = 0` (an *analytic Gaussian* f(r,θ)), which is **not** a
stationary eigenstate of the massive Klein–Gordon operator. A Gaussian is a
superposition of many radial modes; the non-eigenstate part simply radiates away
(disperses) whether or not there is a mass term. Adding μ shifts the dispersion
relation (hence the ~1 M delay) but cannot turn a non-eigenstate lump into a
bound one. **Confinement requires the *profile* to be a genuine bound state**, of
which the code already has two (`grtresna/solver/params.py`,
`profiles/boson_star.py`):

- `PROFILE_ODE_BOUND` (3): flat-space **Q-ball**, bound by self-interaction —
  needs `scalar_lambda`/`scalar_mu > 0` and a solved radial ODE profile f(r).
- `PROFILE_SELFGRAV_BOUND` (4): **self-gravitating boson star**, bound by gravity;
  its frequency ω is the gravitational eigenvalue (solved, written back to
  `bs_omega`), and the painter sets `Π₂ = −ω·f/α` consistently.

Both emit a tabulated f(r) through the C++ `profile==3` loader — i.e. the
machinery to paint a *real* bound state already exists (it is what the
self-grav boson-star seed campaign uses). The rotating-wormhole lump just wasn't
using it; it was still on the analytic Gaussian.

**Deeper caveat (phantom bound states).** Even with a bound-state profile there
is a genuine physics question: the standard Q-ball / boson-star eigenstates are
solutions for a *canonical* (positive-energy) field. Here the field sources
gravity with the phantom sign (`−support_strength`), and Ellis–Bronnikov-class
throats have an unstable radial mode independent of the matter profile. A bound
profile removes the *dispersive* channel (radiation of non-eigenstate modes) but
may not remove the *instability* channel. Expect it to help substantially, but
budget for the possibility that only **active support** (option #3, the deferred
RL pump / PD trap) or **embracing the instability as the GW signal** (option #6)
fully closes the gap.

**Reproduce:**
```bash
MASS=0.1 RES_N=128 bash grteclyn-wrapper/scripts/wormhole/id/solve_kappa_family.sh 1.0 4
bash grteclyn-wrapper/scripts/wormhole/run/wormhole_case.sh --kappa 1.0 --dx 0.5 --mass 0.1 --gpu 1
```
Run dir `runs/rotating_wormhole/evo_omega_p0p05_m1_kappa_1p00_dx0p5_ml3_mass0p1/`.

**Revised next step (updated):** the confinement lever is now understood — it is
the **bound-state profile**, not merely a mass. Switch the rotating-wormhole lump
to `PROFILE_ODE_BOUND` (Q-ball; set `scalar_lambda`/`scalar_mu > 0` and use the
ODE profile) or `PROFILE_SELFGRAV_BOUND` (boson star), matching the potential in
the evolution. If a bound profile still disperses/destabilises, that is the
physical instability of the phantom throat → pivot to active support (RL pump) or
document the transient GW burst as the result.

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
