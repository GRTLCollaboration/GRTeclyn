# MAP-Elites FTL Discovery — Matter-First Metric Discovery

> Quality-Diversity (MAP-Elites) search over **matter configurations** that, once
> turned into a self-consistent spacetime and evolved, exhibit faster-than-light
> (FTL) precursors. We do **not** hand-design a warp metric and ask "what matter
> supports it?" — we propose matter, solve Einstein's constraints for the
> geometry it actually produces, evolve that geometry, and let the search
> *discover* which matter distributions yield FTL signatures.

## Table of contents

- [The idea: matter-first, not metric-first](#the-idea-matter-first-not-metric-first)
- [The pipeline](#the-pipeline)
  - [Diagram — end-to-end overview](#diagram--end-to-end-overview)
  - [Diagram — matter-first vs metric-first](#diagram--matter-first-vs-metric-first)
  - [Stage 0 — Quality-Diversity proposer (MAP-Elites)](#stage-0--quality-diversity-proposer-map-elites)
  - [Stage 1 — Initial data (GRTresna, CPU/MPI)](#stage-1--initial-data-grtresna-cpumpi)
  - [Stage 2 — Evolution (GRTeclyn, GPU)](#stage-2--evolution-grteclyn-gpu)
  - [Stage 3 — Metrics & probes (scoring)](#stage-3--metrics--probes-scoring)
  - [Stage 4 — Archive update & feedback](#stage-4--archive-update--feedback)
- [The hard consistency rule](#the-hard-consistency-rule)
- [Behavior descriptors (the "diversity" axes)](#behavior-descriptors-the-diversity-axes)
- [Scoring model (the "quality" axis)](#scoring-model-the-quality-axis)
  - [Plain-English glossary: every metric & penalty](#plain-english-glossary-every-metric--penalty)
- [Code map (where everything lives)](#code-map-where-everything-lives)
- [Building the binaries (GRTresna + GRTeclyn)](#building-the-binaries-grtresna--grteclyn)
- [How to run a campaign](#how-to-run-a-campaign)
- [Campaign log / runs analysis](#campaign-log--runs-analysis)

## The idea: matter-first, not metric-first

Classic warp-drive analysis is **metric-first**: write down a target metric
(Alcubierre, Natário, …), then read off the stress-energy `T_ab = G_ab/8π` it
requires — which is always exotic, often pathological, and never solved as a
*consistent initial-value problem*. The "FTL" in those constructions is baked
into the chosen coordinates and need not survive a real evolution.

This project inverts that. The pipeline is **matter-first**:

1. The search proposes a **matter configuration** (a set of massive scalar-field
   "lumps" with positions, boosts, rotations, mass, and exotic flags).
2. GRTresna **solves the Einstein constraint equations** for the conformal factor
   and shift that this matter actually sources — a genuine, on-constraint initial
   data set, not a coordinate guess.
3. GRTeclyn **evolves** that spacetime forward in time on GPU.
4. Probes measure whether an FTL signature (coordinate-speed channel, sustained
   superluminal region, and ultimately a **gauge-invariant null-geodesic
   shortcut**) emerges and *persists*.
5. MAP-Elites stores the best candidate per behavior cell and mutates elites to
   explore the space.

The payoff: any FTL signal that survives this loop is a property of a
**self-consistent, evolved spacetime**, not of a hand-picked metric. The
open scientific question the campaign attacks is *which matter distributions, if
any, produce a real (gauge-invariant) shortcut* — and the metrics themselves are
iteratively hardened so the leaderboard cannot be gamed by coordinate artifacts
(see the campaign log for the bug-fix history that got us here).

## The pipeline

### Diagram — end-to-end overview

The closed loop: search proposes matter → physics builds and evolves a real
spacetime → metrics discover FTL signatures → archive feeds the next proposal.

![MAP-Elites end-to-end overview](mapelites-end-to-end.svg)

### Diagram — matter-first vs metric-first

![Matter-first vs metric-first diagram](mapelites-matter-first.svg)

---

### Stage 0 — Quality-Diversity proposer (MAP-Elites)

MAP-Elites maintains an 8×8 behavior archive keyed by a 2-D descriptor. Each
batch it either **mutates an existing elite** (boundary-reflected Gaussian
perturbation, σ=0.15, ~85% of draws) or **samples inside the feasible box** of
known elites (falling back to uniform sampling of the full space before any
elites exist). Proposals are points in the `grtresna_shell` search space — a 19-D
parameterization of the lump field (amplitudes, mode, thickness, toroidal /
poloidal / radial velocities, shift seed, exotic phase, **scalar mass**,
**static toggle**, …).

> **Not CMA-ES.** The QD campaign is pure MAP-Elites: it illuminates the whole
> behavior grid via elite-mutation + feasible-box sampling and never runs a
> CMA-ES generation. The CMA-ES driver (`run_optimize` in `optimize.py`) is a
> *separate*, single-objective optimizer used by the non-QD `run_grtresna_search.sh`;
> the `qd` command never calls it. The two only share the `SearchDimension`
> search-space definitions, which happen to live in `optimize.py` (see Code map).

### Stage 1 — Initial data (GRTresna, CPU/MPI)

The sibling repo `../GRTresna` paints the proposed lumps into `(φ, Π)` and runs a
Lichnerowicz/York elliptic **constraint solve** for the conformal factor and
vector potential. This is CPU/MPI-bound and is the main throughput bottleneck:
poorly-conditioned proposals fail to converge below the Ham/Mom gate and are
**rejected before reaching the GPU** (logged but not inserted into the archive).

### Stage 2 — Evolution (GRTeclyn, GPU)

Converged initial data is handed to GRTeclyn (`RadialRecipe` example), which
evolves the BSSN/CCZ4 system + matter forward to `STOP_TIME` on GPU, dumping
plotfiles every `PLOT_INTERVAL` steps.

### Stage 3 — Metrics & probes (scoring)

Plotfiles are post-processed (`yt`) into diagnostics: constraint norms,
apparent-horizon `theta_plus`, comoving / shift statistics, matter energy density,
and the FTL probe family — coordinate-speed (`operational_ftl_solved`), sustained
evolved speed (`operational_ftl`, `ftl_persistence`), 3-D morphological coherence,
and the gauge-invariant **null-geodesic shortcut** (`operational_ftl_geodesic`).

### Stage 4 — Archive update & feedback

`score.py` collapses the diagnostics into a single `ftl_first` fitness. Only
`gpu_ok` candidates compete for a behavior cell; the new candidate replaces the
incumbent if it scores higher. Rejected GRTresna solves are still written to
`trajectory.jsonl` but never pollute the archive grid. The updated archive feeds
the next proposer step.

## The hard consistency rule

`T_ab` used in the GRTresna **constraint solve** must equal `T_ab` used in the
GRTeclyn **evolution**. Otherwise the run starts off-constraint and any apparent
"FTL" is a constraint-relaxation transient, not physics.

The `grtresna_independent_scalars` matter path exists precisely to keep both
sides identical; any new matter sector (mass search, `λφ⁴`, complex field) must be
added to **both** sides with matching analytic forms. (Root cause documented in
`Examples/RadialRecipe/Debug.md`; see also
[Matter model — reference & future directions](#matter-model--reference--future-directions-2026-06-10).)

## Behavior descriptors (the "diversity" axes)

The current descriptor is **`ftl_lifetime`** (`qd_search/descriptors.py`), added
in [v15](#v15-time-resolved-intermediate-ftl-scoring-2026-06-13):

- **x** — peak gauge-invariant strength: the run's maximum trustworthy `f_geo`,
  ramped `(f_geo − 1e-3)/(2e-1 − 1e-3)`.
- **y** — **FTL-lifetime fraction**: the share of plotfile frames in which a
  trustworthy shortcut is alive. This is the transient-vs-sustained axis — a
  one-frame Alcubierre-like spike and a stable warp of equal peak strength land
  in different cells.

Earlier descriptors are kept for back-compat: **`speed_super`** (x = recalibrated
cone-tilt from `max_local_speed`, y = `superluminal_fraction`) was the v14
default; `speed_horizon` was retired after the `theta_plus` centering bug made
its y-axis degenerate (see campaign log).

## Scoring model (the "quality" axis)

Fitness is the `ftl_first` objective in `metrics/score.py`. Headline structure:

- **Gauge-invariant FTL is king — now time-averaged.** `operational_ftl_geodesic`
  (null-ray shortcut) carries the largest weight, but only when its reliability
  gate passes (`h_quality_ok` via relative Hamiltonian drift, and all rays reach
  the detector). Otherwise it scores 0 — integration noise is never trusted.
  Since [v15](#v15-time-resolved-intermediate-ftl-scoring-2026-06-13) this term is
  the **mean over the whole run** of the per-frame magnitude (each frame gated on
  its own reliability), not a single final-frame reading — the shortcut peaks
  mid-run and diffuses, so the average is what separates a sustained warp from a
  transient collapse.
- **Dynamical, sustained signal next.** Evolved `operational_ftl` + `ftl_persistence`
  outweigh the one-shot coordinate-speed `operational_ftl_solved`, which is
  treated as a *precursor/shaping* term (localization-gated). `operational_ftl`
  is itself **zeroed when a trustworthy geodesic probe finds no shortcut** — the
  gauge-invariant arbiter wins, so a coordinate cone-tilt artifact cannot bank
  the evolved-FTL reward (see the [`v9` review](#ftl_discovery_v9-review--shaping-rebalance--hq-promotion-2026-06-11)).
- **Persistence-honest health block.** `survival = numerical_survival ×
  structural_persistence`, where `structural_persistence = density_retention ×
  morphological_coherence`. A configuration that dissipates or fragments can no
  longer bank FTL credit — the shaping rewards are multiplied by persistence.
- **Vetoes / penalties.** Horizon formation, energy-condition / exotic use,
  instability, and stationary "warp-lens" coordinate artifacts are penalized or
  zeroed so the leaderboard tracks genuine, evolving, bound geometries.

Weight hierarchy (approximate): `operational_ftl_geodesic` ×1000 ≫
`operational_ftl` ×400 (geodesic-gated) ≫ `ftl_persistence` ×300 ≫
`channel_progress` ×100 > `operational_ftl_solved` ×50 (coordinate shaping, kept
below a realistic geodesic shortcut so it cannot out-vote a validated one) ≫
health block ×(survival 70 + stability 10 + …).

The geodesic component is a linear ramp `(f_geo − 1e-3)/(2e-1 − 1e-3)` clamped to
`[0, 1]`: full marks require a **dramatic ~20% null arrival-time shortcut**, so the
scalar tracks shortcut *magnitude* rather than saturating the moment the floor is
crossed (see the [geodesic-reward recalibration](#geodesic-reward-recalibration--ftl_discovery_v9-2026-06-11)).

The full per-component table lives in
`grteclyn-wrapper/src/grteclyn_wrapper/metrics/README.md`.

### Plain-English glossary: every metric & penalty

What each score component actually means, in three buckets: **did it make FTL?**,
**is it a real, healthy structure?**, and **penalties for cheating or breaking
physics**. (Two scoring modes exist: `weighted` is a plain weighted sum;
`ftl_first` — what every campaign uses — makes the validated FTL evidence dominate
and keeps the shaping hints subordinate.)

**1. The goal — did it make a faster-than-light shortcut?**

- **`operational_ftl_geodesic`** — *the gold standard.* Fires a real light ray
  through the evolved spacetime and checks whether it genuinely arrives earlier
  than it could through empty flat space. Gauge-invariant (can't be faked by a
  coordinate trick), so it carries the largest weight. Only counts when its
  reliability gate passes (rays stay on the constraint surface and all reach the
  detector) **and** the matter structure producing it persists (persistence-gated).
- **`operational_ftl`** — a shortcut measured on the evolved geometry via
  coordinate light-speed. Weaker evidence; it is **zeroed when a trustworthy
  geodesic probe finds no shortcut** (then it's just a coordinate illusion).
- **`ftl_persistence`** — rewards a shortcut that *lasts* across the final frames,
  not a one-frame flicker.
- **`operational_ftl_solved`** — a shortcut seen in the initial (t=0)
  constraint-solved data, before evolution. A hint, not proof (localization-gated).
- **`ftl_precursor`, `channel_progress`, `shift_drive`** — "you're getting warmer"
  gradients: light cones starting to tilt and frame-dragging starting to appear
  *before* a full shortcut exists, so the search has a slope to climb out of flat
  space. Gated down if the matter fragments.
- **`ftl_shortcut`** — a faint t=0 hint, barely weighted.

**2. Is it a real, healthy structure? (rewards)**

- **`numerical_survival`** — did the simulation run to the end without crashing?
  (Necessary but not sufficient — empty space also "survives".)
- **`structural_persistence`** — the honest survival measure, a product of two
  things: *density retention* (did the matter keep its energy density or fizzle
  out?) × *morphological coherence* (did it stay one blob or shatter into lobes?).
- **`survival`** = `numerical_survival × structural_persistence` — "it actually
  lasted as a coherent thing".
- **`stability` / `comoving_stability`** — how little the geometry drifts/wobbles.
- **`constraint_health`, `constraint_growth`, `initial_constraint_quality`** — is
  it actually solving Einstein's equations well, or accumulating numerical error?
- **`lapse_health`** — is the time-slicing well-behaved?
- **`energy_condition` / `anec_condition`** — rewards for respecting physical
  energy rules (positive energy). Weighted heavily so physical solutions are favored.
- **`tidal_comfort`** — would a passenger survive the tidal forces?
- **`curvature_activity`, `nontrivial_geometry`, `nonflat_geometry`,
  `expansion_asymmetry`** — rewards for the geometry being genuinely *warped*, so
  flat empty space can't win by default.

**3. Penalties — for cheating or breaking physics (negative)**

- **`exotic_penalty`** — the big one. We want FTL *without* exotic (negative-energy)
  matter; the more negative energy the solution needs, the larger the penalty.
  Graded (0..−1.6) and calibrated so it's a real ~20–30% dent against a shortcut
  without swamping it (a too-hot weight floors the whole population negative).
- **`horizon_penalty`** — it collapsed into a black hole / trapped surface. Heavily
  punished.
- **`instability_penalty`** — it blew up or went numerically wild.
- **`qei_penalty`** — violates the quantum energy inequality (a deeper limit on
  negative energy).
- **`boundary_penalty`** — junk reflecting off the edges of the box contaminated the
  result.
- **`stationary_artifact_penalty`** — catches a *static* lens pretending to be FTL
  (a frozen geometry that bends light but transports nothing). Graded so
  "almost moving" is penalized less than "perfectly frozen".

**The non-triviality gate** wraps the whole thing: the health/niceness rewards
(survival, stability, …) are switched off for flat empty space, so Minkowski
vacuum — which is maximally healthy and stable — can't farm points and beat a
genuinely warped geometry.

## Code map (where everything lives)

| Concern | Path |
|---------|------|
| **MAP-Elites** QD loop, archive, elite-mutation, descriptors | `grteclyn-wrapper/src/grteclyn_wrapper/search/qd_search.py` |
| Search-space (`SearchDimension`) defs + param overrides — *also* hosts the separate, unused-by-QD CMA-ES `run_optimize` | `grteclyn-wrapper/src/grteclyn_wrapper/search/optimize.py` |
| Scoring / fitness | `grteclyn-wrapper/src/grteclyn_wrapper/metrics/score.py` |
| Metric aggregation | `grteclyn-wrapper/src/grteclyn_wrapper/metrics/aggregation/collector.py` |
| FTL probes (general / geodesic) | `grteclyn-wrapper/src/grteclyn_wrapper/metrics/probes/ftl/` |
| Diagnostic dataclasses | `grteclyn-wrapper/src/grteclyn_wrapper/metrics/types/diagnostics.py` |
| Plotfile → frames | `grteclyn-wrapper/src/grteclyn_wrapper/visualisation/process_wave/consume_plotfiles.py` |
| Matter (evolution side) | `Source/Matter/GRTresnaIndependentScalars.{hpp,impl.hpp}`, `Examples/RadialRecipe/` |
| Matter (initial-data side) | `../GRTresna/Examples/ScalarFieldBH/` |
| Campaign launcher | `grteclyn-wrapper/scripts/search/run_grtresna_qd_search.sh` |

## Building the binaries (GRTresna + GRTeclyn)

Two separate toolchains. **GRTresna** (initial data) is **Chombo + conda-OpenMPI,
CPU/MPI**; **GRTeclyn** (evolution) is **AMReX + CUDA, GPU**. The recurring pain is
the GRTresna MPI/HDF5 toolchain: its Chombo `Make.defs.local` resolves
`mpicxx`/`gfortran`/`g++` and `$CONDA_PREFIX` from the **conda env**, so the build
(and link, and run) fail unless that env is on `PATH`/`LD_LIBRARY_PATH` and
`CONDA_PREFIX` points at it.

### One env to set first (every shell)

```bash
# Adjust these two paths to your machine; everything below is derived.
export GRTRESNA_ENV=/home/jovyan/.mlspace/envs/grtresna
export SIM_ROOT=/home/jovyan/nachevsky/test/simulation

export CHOMBO_HOME="${SIM_ROOT}/Chombo/lib"
export CONDA_PREFIX="${GRTRESNA_ENV}"
export PATH="${GRTRESNA_ENV}/bin:${PATH}"
export LD_LIBRARY_PATH="${GRTRESNA_ENV}/lib:${LD_LIBRARY_PATH:-}"
```

Shortcut: `source grteclyn-wrapper/scripts/lib/env.sh` sets the same vars
(`GRTRESNA_ENV`, `PATH`, `LD_LIBRARY_PATH`, `CONDA_PREFIX`). The conda env must
already contain `openmpi openmpi-mpicxx "hdf5=*=mpi_openmpi*" gcc/gxx/gfortran_linux-64`
plus the `mpicxx`/`g++`/`gcc`/`gfortran` symlinks in `$CONDA_PREFIX/bin` — see
`grteclyn-wrapper/src/grteclyn_wrapper/README.md` ("Installing GRTresna from
scratch") for the full micromamba recipe.

### Build GRTresna (initial-data solver, MPI)

```bash
cd "${SIM_ROOT}/GRTresna/Examples/ScalarFieldBH"
PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
  make all -j4 CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE
# -> Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex
```

First time only, build Chombo once: `cd "${CHOMBO_HOME}" && make lib -j"$(nproc)"`.

**Header-only edits** (e.g. `MatterParams.hpp`, `ScalarField.hpp`) often don't
trigger a recompile — force the relink:

```bash
cd "${SIM_ROOT}/GRTresna/Examples/ScalarFieldBH"
EXE=Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex
rm -f "${EXE}" \
  o/3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI/{Main_ScalarFieldBH,Grids,ScalarField,MyMatterFunctions}.o
PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
  make all -j4 CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE
```

### Build GRTeclyn (evolution, GPU)

`AMREX_HOME` defaults to `${SIM_ROOT}/amrex` (sibling of GRTeclyn), so no extra env
is needed beyond a CUDA toolchain and `COMP=gnu`.

```bash
cd "${SIM_ROOT}/GRTeclyn/Examples/RadialRecipe"
# Single-GPU (what the smoke pipeline builds):
make COMP=gnu USE_CUDA=TRUE USE_MPI=FALSE CUDA_ARCH=90 -j"$(nproc)"
#   -> main3d.gnu.CUDA.ex
# Multi-GPU (mpirun comes from the conda/local OpenMPI on PATH):
make COMP=gnu USE_CUDA=TRUE USE_MPI=TRUE  CUDA_ARCH=90 -j"$(nproc)"
#   -> main3d.gnu.MPI.CUDA.ex
```

Set `CUDA_ARCH` to your GPU: `90` = Hopper/H100, `80` = A100, `70` = V100.

### Common failures → fixes (the MPI/conda gotcha)

| Symptom | Cause | Fix |
|---------|-------|-----|
| `mpicxx: command not found` / `gfortran: command not found` during GRTresna build | conda env not on `PATH` | `export PATH="${GRTRESNA_ENV}/bin:${PATH}"` (or `source .../scripts/lib/env.sh`) |
| `cannot find -lhdf5` / HDF5 headers missing | `CONDA_PREFIX` unset, so Chombo `HDF*FLAGS` point nowhere | `export CONDA_PREFIX="${GRTRESNA_ENV}"` before `make` |
| Builds, but `mpirun ./...ex` fails with `libmpi.so`/`libhdf5.so not found` | runtime loader can't see the conda libs | `export LD_LIBRARY_PATH="${GRTRESNA_ENV}/lib:${LD_LIBRARY_PATH}"` |
| `CHOMBO_HOME` undefined / Chombo libs not found | env not exported | `export CHOMBO_HOME="${SIM_ROOT}/Chombo/lib"` (and build Chombo once) |
| A C++ header edit "did nothing" | Chombo make skips header-only deps | force-relink (rm the `.ex` + the listed `.o` files) |
| `/bin/csh: No such file` during Chombo/GRTresna make | Chombo hardcodes `/bin/csh` | point it at conda `tcsh` (README "Installing GRTresna from scratch", step 6) |

## How to run a campaign

```bash
cd grteclyn-wrapper
QD_NAME=ftl_discovery_vN QD_ITERATIONS=10 BINS=8 STOP_TIME=8.0 \
  GPU_IDS="0 1 2 3 4 5 6 7" RANKS=8 LUMPS=5 SHELL_PROFILE=compact \
  GRTRESNA_MAX_HAM_PCT=5.0 GRTRESNA_MAX_MOM_PCT=5.0 \
  nohup bash scripts/search/run_grtresna_qd_search.sh \
  > ../runs/qd_ftl_discovery_vN.launch.log 2>&1 &
```

Results land in `runs/grtresna_qd/ftl_discovery_vN/` (`trajectory.jsonl`, per-eval
`score.json`, `frames/`). The campaign log below records what each `vN` changed
and what it found.

## Campaign log / runs analysis

Reverse-chronological isn't enforced below; entries were appended as work
happened. Quick index (most consequential first):

| Campaign / section | Date | Headline |
|--------------------|------|----------|
| [**v15: time-resolved FTL scoring**](#v15-time-resolved-intermediate-ftl-scoring-2026-06-13) | 06-13 | Final-frame scoring was half-blind: the gauge-invariant shortcut **peaks mid-run and diffuses**, so the last frame both under-credits a real transient and can't tell a sustained warp from an Alcubierre-like collapse. Adds an in-flight per-plotfile FTL stream (`ftl_timeseries.dat`, process+delete), retargets the headline `operational_ftl_geodesic` to the **time-average** over the run, and adds an `ftl_lifetime` MAP-Elites axis. Validated on eval 231: f_geo rises 2.7%→**7.43% peak at t=9.6**→5.24% (t=16); the old final frame saw only 5.24%. QD now runs at HQ refinement (ml=3, dx=0.5) |
| [**v14 results & analytics**](#v14-campaign-results--analytics-2026-06-12-completed) | 06-12 | 504 evals, 351 gpu_ok, 51.6% archive coverage. Top: Eval 231 (f_geo=5.30%, ring+top-hat, score 551). 5 operational, 3 observer_ec. Ring layout dominates top-5; exotic fraction 90–99% universal. Full Alcubierre comparison — our best is 17% of Alcubierre's shortcut but is self-consistent & evolvable |
| [`v14` launch setup → matter profile + cloud layout](#v14-launch-setup-matter-profile-and-cloud-layout-2026-06-12) | 06-12 | Adds per-lump matter profile (Gaussian / smoothed top-hat "ball") + quasi-random cloud layout (`matter_layout=4`); search space 21→23 dims. GRTresna rebuilt; 182 tests pass |
| [`v12` review → λφ⁴ + FTL geometry layouts → `v13`](#ftl_discovery_v12-review--lambda-phi4--ftl-geometry-layouts--ftl_discovery_v13-2026-06-12) | 06-12 | 278 evals, zero geodesic FTL; top scores were coordinate-shaping artifacts (eval 197 scored 130 with f_geo=0). Adds searchable `grtresna_scalar_lambda` + `grtresna_matter_layout` (sphere/channel/bipolar/ring), zeros shaping when geodesic contradicts, pytest gate before QD launch |
| [Alcubierre positive control → metric-first vs matter-first verdict](#alcubierre-positive-control--metric-first-vs-matter-first-verdict-2026-06-12) | 06-12 | Prescribed Alcubierre metric → our probes detect a 32% shortcut + flag exotic matter (probes validated). Metric-first is fine for *analysis*, impossible for *dynamics*; matter-first is correct. QD-res H-gate rejects even textbook Alcubierre → added 129³ mid-res re-probe |
| [`v10` review → persistence-gate + physicality pressure → `v11`](#ftl_discovery_v10-review--persistence-gate--physicality-pressure--ftl_discovery_v11-2026-06-12) | 06-12 | 400 evals, top-5 all dynamic exotic bubbles with transient 2–3% shortcuts (same HQ-death band); #1 fragmented (persistence 0.46) yet ranked top. Persistence-gate the geodesic reward + raise exotic/energy weights |
| [HQ verdict: shortcuts did not survive refinement → `v10`](#hq-verdict-shortcuts-did-not-survive-refinement--ftl_discovery_v10-2026-06-11) | 06-11 | All 3 promoted shortcuts collapsed at HQ (`f_geo` 2–3% → 0 / 0.29%); pipeline honestly rejected its own elites. Extend QD to t=16 + add static-matter toggle |
| [`v9` review + shaping rebalance → HQ promotion](#ftl_discovery_v9-review--shaping-rebalance--hq-promotion-2026-06-11) | 06-11 | Geodesic gate fires for all evals (5 real shortcuts found); coordinate precursor out-voted validated ones → rebalanced; top 3 promoted HQ |
| [Geodesic-reward recalibration → `v9`](#geodesic-reward-recalibration--ftl_discovery_v9-2026-06-11) | 06-11 | `v8` eval 11 scored 1066 off a real but *modest* 3.3% shortcut; rescaled so the scalar reflects magnitude |
| [Null-geodesic reliability fix → `v8`](#null-geodesic-reliability-fix-2026-06-11-post-v7) | 06-11 | Forward-launch rays + relative-drift gate; gauge-invariant term finally fires |
| [`ftl_discovery_v7` review](#campaign-ftl_discovery_v7--finished-run--critical-leaderboard-review-2026-06-11) | 06-11 | Persistence/coherence honest; geodesic still blind (the bottleneck) |
| [`ftl_discovery_v4`](#campaign-ftl_discovery_v4--persistence-honest-scoring--bound-matter-2026-06-11) | 06-11 | Persistence-gated `survival`; searched `scalar_mass`; capped boosts |
| [Matter model reference & future directions](#matter-model--reference--future-directions-2026-06-10) | 06-10 | What the lumps are; how they interact; roadmap |
| [`ftl_discovery_v2`](#campaign-ftl_discovery_v2--first-healthy-run--a-scoring-concern-2026-06-10) | 06-10 | First sane scoring run; gauge-invariance under-weighted |
| [Scoring fix: stationary warp-lens artifacts](#scoring-fix-stationary-warp-lens-artifacts-2026-06-10-after-90-evals) | 06-10 | Reliability-gate + stationary-artifact gate |
| [Navigation overhaul](#navigation-overhaul-2026-06-10) | 06-10 | `speed_super` descriptor; feasible-box sampling |
| [Status / reset (theta_plus bug)](#map-elites-ftl-discovery-status) | 06-10 | `theta_plus` re-centered on `grid_center` |

---

## v15: time-resolved (intermediate) FTL scoring (2026-06-13)

**The problem.** Every campaign through v14 scored a candidate from a **single
plotfile — the last one**. The geodesic probe samples one static slice (time
derivatives zeroed) and reports `f_geo` there. But the HQ-promotion forensics on
the v14 top-5 showed the gauge-invariant shortcut is a **transient**: it forms,
**peaks mid-run, and diffuses**. Reading only the final frame is therefore
half-blind in two ways:

1. It **under-credits** a genuine channel whose peak already passed (eval 231:
   `f_geo` peaks at **7.43% near t=9.6** but is **5.24%** by t=16 — and was ~0 by
   the t=30 HQ frame).
2. It **cannot distinguish a sustained warp from a collapsing one.** A single
   snapshot of Alcubierre looks stable and interstellar; a single snapshot of an
   unstable bubble at its peak looks identical. Only the *time history* tells
   them apart. (This is exactly why the v10–v14 "winners" died at HQ.)

**The fix — measure the whole curve, in flight.** As the plotfile consumer
streams and **deletes** each plotfile (no hoarding — disk stays bounded), it now
also runs the FTL probes on that frame and appends one row to
`small_data/ftl_timeseries.dat`, alongside the other per-plotfile diagnostics:

```
# time  f_op  f_geo  geo_trustworthy  max_local_speed  superluminal_fraction  max_shift  structure_coherence  reachable  n_rays  n_reached  max_h_rel_drift
```

The cheap operational-FTL probe runs on every frame; the expensive
gauge-invariant geodesic probe is **gated** (only fires when the cheap probe sees
a coordinate channel), so per-frame cost stays ~3–7 s. The collector parses this
into `FtlTimeSeriesMetrics` (per-frame arrays + peak, time-of-peak, FTL-lifetime
fraction).

**The new headline score (avg / sum / divide).** `operational_ftl_geodesic` — the
×1000 dominant term — is now the **mean over all frames** of the per-frame
trustworthy gauge-invariant magnitude `(f_geo − 1e-3)/(2e-1 − 1e-3)`, still scaled
by end-state `structural_persistence`. Each frame is gated on its *own*
reliability (`h_quality_ok`, full ray bundle). This is the persistence-honest
average: a channel that is FTL for most of the run scores high; a one-frame
Alcubierre-like spike averages toward zero; a diffused final frame no longer
zeroes a genuinely sustained shortcut. Falls back to the single final-frame value
when no per-frame stream exists (backward-compatible). Worked example (identical
peak, different lifetime): persistent 10/10 frames → **296**, transient 5/10 →
148, one-frame spike → **30**.

**New MAP-Elites axis — `ftl_lifetime`** (now the campaign default): **x** = peak
gauge-invariant strength, **y** = FTL-lifetime fraction (share of frames the
shortcut is alive). The archive now explicitly separates *transient* shortcuts
from *sustained* ones at the same strength — the transient/stable distinction is
a first-class diversity dimension instead of being invisible.

**QD now runs at HQ refinement.** Resolution gap was the last suspect behind the
QD→HQ collapse, so the QD loop itself now evolves at **`max_level=3`, `dx=0.5`**
(`N_full=128`, `L_full=64`), `t=16` — the same grid an elite is promoted at. The
control experiments (dx=0.5 alone, ml=3 alone) had already shown resolution does
*not* kill the shortcut at t=16; the real killer is **time** (diffusion by t≈30),
which the time-average now captures directly.

### Validation (eval 231 replay, gridinit reused, GRTresna skipped)

Re-evolving eval 231 at the new QD grid (dx=0.5, ml=3, t=16) with the per-frame
stream on, processing + deleting plotfiles:

| t | `f_geo` | trustworthy |
|------|---------|-------------|
| 0.0  | 2.70%   | yes (5/5) |
| 3.2  | 4.86%   | yes |
| 6.4  | 7.18%   | yes |
| 9.6  | **7.43% (peak)** | yes |
| 12.8 | 6.34%   | yes |
| 16.0 | 5.24%   | yes |

Time-averaged magnitude = **0.275** vs the single final frame's 0.258 — the
mid-run insight *raises* the score, as it should. `ftl_lifetime = 100%`
(FTL in every frame → sustained, top lifetime cell). Confirms the mechanism
end-to-end: in-flight extraction → bounded disk → time-resolved aggregation →
time-averaged scoring → lifetime descriptor.

### Where it lives

- **In-flight extraction:** `visualisation/process_wave/consume_plotfiles/extraction/ftl.py`
  (`_extract_ftl_timeseries_line`); wired via `--ftl-timeseries`/`--ftl-l` through
  `core/plot_consumer.py` → `core/runner.py` → `core/evaluation.py` (auto-on
  whenever plotfiles are consumed).
- **Aggregation:** `metrics/diagnostics/ftl_timeseries.py`
  (`read_ftl_timeseries_metrics` → `FtlTimeSeriesMetrics`); collected in
  `metrics/aggregation/collector.py`.
- **Scoring:** `metrics/score.py` (time-averaged `operational_ftl_geodesic`).
- **Descriptor:** `search/qd_search/descriptors.py` (`ftl_lifetime` mode);
  CLI choice in `cli/parser.py`.

---

## v14 launch setup: matter profile and cloud layout (2026-06-12)

`v14` keeps the entire v13 stack intact and adds two matter-distribution knobs.
Nothing here is tuned to the FTL score or shaped toward a warp metric — the goal
is to let MAP-Elites explore a broader, more physical matter palette.

### New in v14 (matter distribution)

1. **Per-lump matter profile.** Each lump carries a `profile` selector: `0` =
   Gaussian envelope `exp(−r²/2w²)` (default → v13 unchanged), `1` = smoothed
   top-hat "ball" `½(1 − tanh((r − w)/0.25w))` (near-uniform, more volumetric,
   soft edge). In the shell ansatz it is driven by
   `grtresna_shell_profile_fraction ∈ [0,1]` (+ `grtresna_shell_profile_phase`),
   assigning the top-hat to a searchable subset of lumps (mirrors the
   `exotic_fraction` idiom); fraction `0` keeps every lump Gaussian. The envelope
   is **byte-identical** in GRTresna `MatterParams.hpp::lump_phi` (constraint
   solve + painter) and the Python gridinit painter `lump_fields.py`, so the
   solved metric and the evolved φ/Π stay consistent. Gradients/Π adapt
   automatically (finite differences). The free/indexed ansatz honours a per-lump
   `grtresna_lump{k}_profile`.

2. **Quasi-random cloud layout** (`grtresna_matter_layout = 4`). Deterministic
   low-discrepancy (R3) scatter of the ≤5 lumps inside a bounded ball, oriented by
   the searched axis frame — a reproducible, fully asymmetric distribution the
   symmetric sphere/channel/bipolar/ring layouts cannot represent.

Search space: **23 dimensions** (was 21).

### Carried over from v13 (independently reviewed, sound)

3. **λφ⁴ self-interaction** (`grtresna_scalar_lambda ∈ [0, 0.1]`, default 0):
   `V = ½(m·Σφ)² − (λ/4)(Σφ)⁴`, `dV/dφ = m²Σφ − λ(Σφ)³`. Threaded identically
   through GRTresna (`scalar_lambda`) and GRTeclyn (`recipe_scalar_lambda`,
   `GRTresnaScalarPotential.hpp`); λ=0 recovers v12.
4. **Matter layout topologies** `grtresna_matter_layout`: `0`=sphere (v12
   Fibonacci shell, bit-for-bit preserved), `1`=channel, `2`=bipolar, `3`=ring —
   now extended with `4`=cloud (above).
5. **Geodesic contradiction gate**: when a trustworthy geodesic probe reports
   `f_geo ≤ floor`, the FTL shaping terms (`operational_ftl_solved`,
   `ftl_precursor`, `channel_progress`, `shift_drive`) are zeroed, stopping
   coordinate-artifact leaderboard inflation.

### Build / rebuild status

- **GRTresna `ScalarFieldBH` recompiled** with the new `profile` field (verified
  clean build + relink; `MatterParams.hpp`, `ScalarField.cpp`,
  `MyMatterFunctions.cpp` all rebuilt). A stale binary silently ignores `profile`.
- **No GRTeclyn `Source/` change** for v14 — the evolution loads φ/Π from the
  gridinit, so the existing `main3d.gnu.CUDA.ex` is current for this change.
- Preflight: `uv run pytest` → 182 passed; C++/Python envelope agreement verified
  to 1e-9.

### Launch command

```bash
cd grteclyn-wrapper
QD_NAME=ftl_discovery_v14 QD_ITERATIONS=10 BINS=8 STOP_TIME=16.0 \
  GPU_IDS="0 1 2 3 4 5 6 7" RANKS=8 LUMPS=5 SHELL_PROFILE=compact \
  GRTRESNA_MAX_HAM_PCT=5.0 GRTRESNA_MAX_MOM_PCT=5.0 \
  nohup bash scripts/search/run_grtresna_qd_search.sh \
  > ../runs/qd_ftl_discovery_v14.launch.log 2>&1 &
```

The launcher runs the pytest preflight gate (`SKIP_QD_PREFLIGHT_TESTS=1` bypasses)
before allocating GPUs.

**Open gap:** resolved — the live 500-eval campaign below confirmed the full pipeline
end-to-end (λ + all layouts + top-hat → solve → evolve → geodesic gate firing).

### v14 campaign results & analytics (2026-06-12, completed)

504 evaluations (target 500), 351 completed on GPU. The pre-GPU rejection rate fell
to **30.4%** (v13 was ~82%), confirming the navigation-overhaul sampling fixes.
Archive: 33 elites across the 8×8 grid (**51.6% coverage**).

| Stat | Value |
|------|-------|
| Total evals | 504 |
| gpu_ok | 351 (69.6%) |
| solved_ftl_rejected | 92 (18.3%) |
| grtresna_rejected | 43 (8.5%) |
| grtresna_failed | 18 (3.6%) |
| f_geo > 0 | 55 (15.7% of gpu_ok) |
| f_geo ≥ 5% | 1 |
| operational tier | 5 |
| observer_ec tier | 3 |
| nontrivial tier | 130 |
| constructed tier | 213 |
| Median score | −18.8 |
| Archive coverage | 33/64 cells (51.6%) |

#### Top 5 candidates (ranked by `ftl_first` score)

Raw `f_geo` is the actual gauge-invariant null-geodesic arrival-time shortcut
(vs flat-space expectation); `f_op_ev` is the sustained evolved coordinate-speed
superluminal fraction. The `operational_ftl_geodesic` score component is the
ramped `[0,1]` scalar `(f_geo − 1e-3)/(0.2 − 1e-3)`.

| # | Eval | Score | Tier | f_geo (raw) | f_op_ev | Layout | Profile | λ | Exotic frac |
|---|------|-------|------|-------------|---------|--------|---------|---|-------------|
| 1 | 231 | 551.0 | observer_ec | **5.30%** | 5.30% | ring | top-hat (0.66) | 0.066 | 90.6% |
| 2 | 369 | 513.4 | observer_ec | 2.40% | **5.76%** | ring | Gaussian (0.24) | 0.054 | 91.9% |
| 3 | 483 | 441.2 | operational | 3.42% | 5.18% | bipolar | mixed (0.45) | 0.088 | 91.6% |
| 4 | 489 | 340.6 | operational | 3.57% | 4.26% | ring | Gaussian (0.15) | 0.056 | 98.4% |
| 5 | 91 | 276.6 | operational | 4.22% | 3.11% | sphere | top-hat (0.79) | 0.081 | 98.9% |

**Eval 231** (top): ring layout, moderate top-hat fraction, mid-range λ. Generates
the campaign's largest geodesic shortcut (5.30%) alongside a 5.30% sustained
coordinate-speed channel. Scored `observer_ec` — the energy-condition violation is
detected (NEC_min = −1.5e-4, WEC_min = −8.5e-4) but the FTL signal dominates.
Described in cell [2,0] (low cone-tilt 0.36, zero superluminal fraction), meaning
its FTL is localized rather than domain-filling.

**Eval 369** (2nd): ring layout, Gaussian-dominant profile. Highest evolved FTL
fraction in the campaign (f_op_ev = 5.76%) with strong persistence (39.8% scaled).
Cell [7,7] — saturates both descriptor axes, indicating domain-wide superluminal
structure. The geodesic shortcut is more modest (2.40%) but the coordinate-speed
channel is the strongest found.

**Eval 483** (3rd): only bipolar-layout candidate in the top 5. Mixed profile
(neither strongly Gaussian nor top-hat). Highest λ (0.088) among top candidates,
suggesting scalar self-interaction helps sustain the channel. Cell [4,3] —
moderate cone-tilt, moderate superluminal fraction.

**Eval 489** (4th): ring layout, strongly Gaussian. Highest exotic fraction
(98.4%). Cell [2,0] — like eval 231, FTL is localized. Survival is slightly
degraded (0.95).

**Eval 91** (5th): the only sphere-layout candidate in the top 5. Strongly
top-hat (0.79) with high λ (0.081) and near-total exotic fraction (98.9%).
Highest raw f_geo among operational-tier candidates (4.22%), but its coordinate
FTL fraction is the lowest (3.11%) and persistence is weak (1.9%), producing a
lower total score despite the stronger geodesic signal.

#### Comparison to Alcubierre (positive control)

| Metric | Alcubierre (129³) | v14 best (Eval 231) | Ratio |
|--------|-------------------|---------------------|-------|
| f_geo (raw shortcut) | **31.5%** | 5.30% | 16.8% |
| f_op_ev | N/A (prescribed metric) | 5.30% | — |
| Exotic matter | Yes (NEC < 0) | Yes (90.6% exotic frac) | — |
| Evolvable | No (frozen background) | **Yes** (scalar fields) | — |
| Self-consistent ID | No (metric-first) | **Yes** (constraint-solved) | — |

The v14 best achieves only ~17% of Alcubierre's shortcut magnitude, but it does
so as a **self-consistent, evolvable, constraint-solved scalar-field spacetime**
— not a prescribed kinematic background. The three reasons the search cannot
reach full Alcubierre-level shortcuts are the same features documented in the
[Alcubierre verdict](#alcubierre-positive-control--metric-first-vs-matter-first-verdict-2026-06-12):
(1) physical scalar matter cannot source the extreme exotic `T`, (2) the
moving-puncture gauge damps seeded warp shift, and (3) exotic matter is penalized.

#### Patterns & takeaways

- **Ring layout dominates.** 3 of the top 5 use the ring layout (v13's layout=3),
  and layout=3 also dominates the archive (109 of 351 completed evals). The ring
  appears to be the most productive topology for seeding FTL channels.
- **Exotic matter is universal.** Every top candidate requires 90–99% exotic
  fraction. The exotic penalty floors many candidates negative (median score −18.8)
  but the top performers overcome it with genuine geodesic shortcuts.
- **Top-hat vs Gaussian is mixed.** Both profiles appear in the top 5; the
  search did not converge on one as clearly superior.
- **λφ⁴ is actively used.** All top candidates have λ > 0 (0.054–0.088), and
  eval 483's high λ (0.088) correlates with strong channel persistence. The search
  is genuinely exploring the self-interaction dimension.
- **observer_ec vs operational.** The top 2 candidates are `observer_ec` (energy
  condition violation detected but FTL dominates); the next 3 are `operational`
  (genuine dynamical FTL). The tier distinction is primarily about whether the
  energy-condition penalty exceeds the FTL reward.
- **No cloud-layout (layout=4) breakthrough.** Despite being a new v14 feature,
  the quasi-random cloud layout did not produce any top-10 candidate. The
  structured layouts (ring, bipolar, sphere) were more productive.
- **Archive diversity is healthy.** 33 of 64 cells filled, spanning all tier
  levels from constructed through observer_ec. The 51.6% coverage with only
  351 gpu_ok evals confirms the descriptor axes are well-calibrated.

## MAP-Elites FTL Discovery Status

Status: **reset**. The previous QD/HQ campaign results were discarded.

### Why the reset

The `theta_plus` diagnostic measured radius from the coordinate origin instead of
`grid_center`, producing false trapped-surface detections and a spurious −1.0
horizon penalty. Fixed by re-centering on `grid_center` in `RadialRecipeLevel.cpp`
and related files; binary rebuilt. Fresh campaign `ftl_discovery_postfix` (8×8
`speed_horizon`, target 400) confirmed the fix — `theta_plus` strictly positive,
`horizon_penalty = 0` on all scored candidates. But ~93% pre-GPU rejection exposed
the navigation defects fixed next.

## Navigation Overhaul (2026-06-10)

Two navigation defects (distinct from the `theta_plus` physics bug) were fixed:

1. **Behavior-space collapse** — `speed_horizon` y-axis degenerated to a single row
   after the centering fix (coverage ~0.078). Replaced with **`speed_super`**
   descriptor: x = recalibrated cone-tilt (`max_local_speed`, floor 0.95 / target
   1.20), y = `superluminal_fraction` (share of domain with local speed > c).
   `speed_horizon` kept for back-compat.
2. **~82% pre-GPU rejection** — GRTresna solve stalled on pathological corners of
   the shell bounds. Fixes: tightened shell bounds to feasible basin, boundary
   reflection instead of hard clipping in mutation, elite-mutation fraction 0.70→0.85,
   feasible-box sampling (`_sample_feasible_box`) for random draws, harder
   GRTresna solve (50 iterations, stall tolerance 0.005).

Launched as `ftl_discovery_nav` (8×8 `speed_super`, target 400 evals). First 16
evals confirmed GPU-reach ~40% (vs ~18%) and x-axis spread. Two follow-up fixes:
y-axis reads *solved* report (not evolved, which decays to ~0) and rescaled by
observed ceiling 0.30; `chi`/`chi-1` frames fixed with per-frame percentile
auto-scaling.

### Scoring fix: stationary warp-lens artifacts (2026-06-10, after ~90 evals)

The run converged into a degenerate basin — all 15 retained elites were stationary
zero-shift geometries (static "warp-lens" artifacts). Top candidate eval_000083
scored 1164 from `operational_ftl_geodesic` computed on a geodesic report flagged
`h_quality_ok=False` (integration noise at full weight). Two fixes in
`metrics/score.py`:

1. **Reliability-gate `operational_ftl_geodesic`** — only certified when
   `h_quality_ok` AND `n_reached == n_rays`; otherwise reward 0.
2. **Stationary-artifact gate** — when stationary AND no dynamical FTL, all
   shaping rewards (`operational_ftl_solved`, `ftl_precursor`, `channel_progress`,
   `shift_drive`) are zeroed.

Effect: eval_000083 1164→−247, eval_000065 270→−246, eval_000094 194→−247.
Regression tests: `test_unreliable_geodesic_shortcut_is_not_rewarded`,
`test_stationary_warp_lens_artifact_ranks_below_genuine_candidate`.

## Matter model — reference & future directions (2026-06-10)

This is a map of *what the lumps actually are* and how to extend them, so we
don't have to re-trace the code each time. The `shell`/MAP-Elites campaign
evolves **N independent massive real scalar fields** ("lumps"), not a single
field and not massless matter.

### What the campaign actually runs

Each eval's `params.txt` sets:

```
recipe_matter_model     = grtresna_independent_scalars
recipe_num_scalar_fields = 5
recipe_scalar_mass       = 0.1
recipe_exotic_matter     = 0     # per-lump exotic flags are used instead
```

So GRTeclyn dispatches to `GRTresnaIndependentScalars` (5 fields, mass m=0.1),
**not** the `ScalarField<DefaultPotential>` / `ExoticScalarField` paths.

### Where the matter lives (file map)

**GRTeclyn (evolution side):**
- `Examples/RadialRecipe/RadialRecipeMatterDispatch.hpp` — runtime model
  selection. `uses_independent_scalars()` picks `GRTresnaIndependentScalars`
  when `recipe_matter_model == "grtresna_independent_scalars"`; otherwise
  `ExoticScalarField<DefaultPotential>` (if `recipe_exotic_matter`) or
  `ScalarField<DefaultPotential>`.
- `Source/Matter/GRTresnaIndependentScalars.{hpp,impl.hpp}` — the per-lump
  fields. `compute_emtensor` sums `sign[k]·T_ab(φ_k)` over lumps;
  `add_matter_rhs` is the curved-space Klein–Gordon RHS per lump.
- `Source/Matter/GRTresnaScalarPotential.hpp` — the **only** non-trivial
  potential currently wired in: `V = ½ m² (Σφ_k)²` (massive, shared sum).
- `Source/Matter/DefaultPotential.hpp` — `V=0` (massless); used by the legacy
  `ScalarField`/`ExoticScalarField` paths. `recipe_scalar_mass` is **ignored**
  on those paths — only `GRTresnaIndependentScalars` honors the mass.
- `Source/Matter/{ScalarField,ExoticScalarField}.{hpp,impl.hpp}` — canonical
  (ρ≥0) and phantom (ρ≤0, T scaled by `-recipe_support_strength`) single-field
  models. Templated on the potential class.
- `Examples/RadialRecipe/StateVariables.hpp` — state layout: `c_phi, c_Pi`
  plus `phi_lump0..4 / Pi_lump0..4` (`2 * GRTRESNA_MAX_INDEPENDENT_SCALARS`).
- `Examples/ScalarField/Potential.hpp` — reference `½ m² φ²` potential class
  (template signature to copy when adding new potentials).

**GRTresna (initial-data side, sibling repo `../GRTresna`):**
- `Examples/ScalarFieldBH/MatterParams.hpp` — the **lump definition**. `lump_t`
  = {amp, width, center, velocity (boost→linear momentum), omega
  (rotation→L_z), mode (0=axisym, 1=dipole, 2=quadrupole), exotic flag}.
  `lump_phi` = `amp·angular_factor(mode)·exp(-r²/2w²)`; `lump_pi` =
  `-(v·∇φ) - omega·∂_azimuth φ` (so each cloud carries net P_i ~ v and L_z ~ omega).
  Exotic lumps are damped by `EXOTIC_AMP_SCALE = 0.25` to stay in the
  Lichnerowicz/York convergent regime.
- `Examples/ScalarFieldBH/MyMatterFunctions.cpp` — paints initial `(phi, Pi)` =
  background + `Σ_k lump_phi/lump_pi`. `my_potential_function` uses
  `½ (scalar_mass·φ)²` (must match GRTeclyn's potential).

### How lumps interact / change shape (answers we keep re-deriving)

- **Shape:** live wave fields under Klein–Gordon — they oscillate, disperse,
  deform. Not rigid. Initial shape is the analytic Gaussian cloud above.
- **Interaction:** only two indirect channels — **(a) gravity** (all lumps
  source one shared metric) and **(b) a shared mass potential** `½m²(Σφ_k)²`
  whose gradient `m²·Σφ_j` is applied to *every* lump's `Pi` RHS
  (`GRTresnaIndependentScalars.impl.hpp`), so overlapping lumps cross-couple.
  There is **no** direct gradient-gradient or self-interaction term.
- **Why they "fly away":** initial boosts are O(1) (`lumpK_velocity ~ |1|`),
  amplitudes are tiny (`amp ~ 0.09` → negligible self-gravity), and the mass is
  light (m=0.1 → Compton wavelength 1/m=10 ≈ box). Nothing binds them, so they
  free-stream out. This is physically expected, not a bug — there is no soliton
  / bound-state mechanism in the current matter sector.

### Future directions (ranked: leverage vs. effort)

1. **Search the mass + cap boosts (config-only, do first).** *Done
   (2026-06-11).* `m=0.1` and O(1) velocities were why lumps dispersed/flew
   away. `grtresna_scalar_mass` is now a searched QD dimension (range
   `0.3–1.5`, default `0.6`) wired into `grtresna_shell_search_space` and
   applied to `cfg.scalar_mass` in `_apply_grtresna_overrides`
   (`search/optimize.py`); it already propagated to both the GRTresna solve
   (`solver.py`) and the GRTeclyn evolution (`matter_wiring.py`), so the
   consistency rule holds with no new physics code. The toroidal current (warp
   motor) keeps its full range, while the net-outflow velocities are capped
   (`poloidal ±1.5→±0.8`, `radial ±0.8→±0.3`) so matter stays bound. Heavier
   mass binds within the shell width (Compton wavelength `1/m`), letting the
   search settle on persistent, near-static matter.
2. **Add `λφ⁴` (or φ⁶) self-interaction ⇒ oscillons / Q-balls** (genuinely bound,
   long-lived). Extend `GRTresnaScalarPotential` (GRTeclyn) **and**
   `my_potential_function` (GRTresna `MyMatterFunctions.cpp`) identically, plus
   one `lambda` param threaded through both param files.
3. **Complex scalar field ⇒ boson stars / Q-balls** with conserved U(1) charge —
   the textbook *stationary, non-dispersing* lump. Larger lift: two real
   components per field, and constraint-solver support on the GRTresna side.
4. **Per-lump independent mass** (currently one shared `scalar_mass`) to mix a
   heavy "anchor" lump with light mobile ones.

### Hard consistency rule (do not violate)

`T_ab` used in the GRTresna **constraint solve** must equal `T_ab` used in the
GRTeclyn **evolution**, or the run starts off-constraint and any "FTL" is a
constraint-relaxation transient (root cause documented in
`Examples/RadialRecipe/Debug.md`). The `grtresna_independent_scalars` path
exists precisely to keep them identical. Any new matter (mass search, `λφ⁴`,
complex field) must be added to **both** sides with matching analytic forms.

## Campaign `ftl_discovery_v2` — first healthy run + a scoring concern (2026-06-10)

Fresh 8-GPU campaign (`runs/grtresna_qd/ftl_discovery_v2/`, `speed_super` 8×8,
`ftl_first`) launched on the vectorized-conversion + corrected-stationarity +
moderated-penalty fixes. This is the first run that scores sanely: gpu_ok
candidates get gradient-rich positives (14→405), `stationary_artifact_penalty`
is 0 on survivors (graded, e.g. −0.955 on the near-stationary eval_011), the
horizon veto fires correctly (eval_026/038 → −1.0), and the archive spreads
across both axes (y-bins 0→7, cells up to [7,7]) — the earlier y-collapse is
fixed.

### Validation of the top elite (eval_000036, score 405, cell [7,7])

Checked in detail because it carries the run's first nonzero
`operational_ftl_solved` (0.717). **Verdict: best *physical precursor* so far,
but NOT confirmed gauge-invariant FTL.** Evidence from `score.json`:

- **Genuinely physical (unlike the old eval_083 lens):** `comoving.stationary =
  false`, `beta_mean = 0.371` (real net shift), final Ham/Mom L2 = 8.2e-4 /
  1.8e-4 (stays on-constraint — not a relaxation transient), `min_theta_plus =
  0.118` (no horizon), survives to t=2.
- **But the 405 is driven by a *coordinate*-speed metric.** `operational_ftl_solved`
  derives from `general_ftl_solved.max_local_speed = 1.287`, explicitly noted as
  a "superluminal **coordinate** channel" — gauge-dependent, a precursor not a proof.
- **The one gauge-invariant test failed its reliability gate:** `geodesic_ftl`
  has `h_quality_ok = False` (null-constraint drift max H=5.4e-4, only 4/5 rays
  reached), so `operational_ftl_geodesic` was correctly zeroed.
- **No sustained signal:** evolved `operational_ftl = 0`, `ftl_persistence = 0`
  (max local speed decays 1.287 → 1.092 over the evolution).
- **Frames** (`shift1_z`, `local_speed_z`, t=0 vs t=2): a real ±0.05 shift dipole
  and a ~1.05–1.10 (on z-slice) superluminal region, both **localized near
  center and weakening in place** — they do not propagate as a warp channel.

So the scoring behaved honestly here: it labeled the result as precursor +
channel progress, *not* certified FTL.

### Concern + plan: gauge-invariant signal is under-weighted

`operational_ftl_solved` (a **coordinate** speed, gauge-dependent) currently
dominates the `ftl_first` objective, while the only gauge-invariant measure
(`geodesic_ftl`, null rays) is frequently rejected for unreliability and
contributes 0. This risks ranking gauge/shift-channel artifacts above genuine
shortcuts. Plan (future work, in priority order):

1. **Fix null-geodesic reliability** so `h_quality_ok` can pass on good
   candidates: tighten ray integration (smaller step / higher-order),
   investigate the constraint-drift source, and ensure all rays reach the
   detector. Until this passes, we have *no* trustworthy gauge-invariant verdict.
   Code: `metrics/null_geodesic.py` (`GeodesicFtlReport`, `h_quality_ok`,
   `n_reached/n_rays`, `max_h_drift`).
2. **Re-weight once geodesic is reliable:** gauge-invariant `operational_ftl_geodesic`
   should outweigh the coordinate-based `operational_ftl_solved`, so the
   leaderboard tracks real shortcuts rather than coordinate channels. Code:
   `metrics/score.py` (`ftl_first` component weights).
3. **Keep `operational_ftl_solved` as a precursor/shaping term only** (lower
   weight), gated as now by non-stationarity, so it guides the search toward
   shift channels without certifying them as FTL.

## Campaign `ftl_discovery_v4` — persistence-honest scoring + bound matter (2026-06-11)

`v2`/`v3` exposed that the leaderboard rewarded structures that *dissipate or
fly apart*: `survival` was `1.0` for HQ elites whose frames showed the lump
gone by the stop time, and the FTL shaping rewards banked cone-tilt off a
structure that no longer existed. Root cause on the matter side: with
`scalar_mass` pinned at `0.1` (Compton wavelength `1/m=10 ≈ box`) nothing binds
the lump, so it free-streams out even at zero boost. The following fixes ship in
`v4` (`runs/grtresna_qd/ftl_discovery_v4/`, `speed_super` 8×8, `ftl_first`,
8 GPUs, `STOP_TIME=8`, `PLOT_INTERVAL=40`).

### Scoring fixes

- **Done — `survival` is now persistence-gated.** Split the old single metric
  into `numerical_survival` (did the integrator reach the stop time — a
  completion gate, weighted lightly) and `structural_persistence` (fraction of
  peak matter energy density retained at the final step,
  `final_peak_rho_required / max_rho_required`). `survival = numerical_survival ×
  structural_persistence`, so a configuration that dissipates can no longer
  score as "stable". Code: `metrics/score.py`, `metrics/types/diagnostics.py`
  (`final_peak_rho_required`), `metrics/diagnostics/constraints.py` (reads it
  from `constraint_norms.dat`).
- **Done — FTL shaping rewards gated by persistence** (commit `be93469`).
  `ftl_precursor` / `channel_progress` / `shift_drive` are multiplied by
  `structural_persistence`, so a fragmented end-state can no longer out-rank a
  coherent survivor by banking cone-tilt/frame-drag credit off a structure that
  has shattered. Persistence defaults to `1.0` when the matter-density series is
  unavailable, leaving the rewards untouched. Regression test:
  `test_ftl_shaping_rewards_scale_with_persistence`.

### Search-parameter fixes

- **Done — `STOP_TIME 2 → 8`, `PLOT_INTERVAL 10 → 40`** (commit `be44e02`). At
  `t=2` a dissipating configuration is indistinguishable from a persistent one
  (both retain ~90% of peak ρ), so the persistence gate could not bite; by `t=8`
  a genuine survivor plateaus (~0.7) while dissipators keep falling (<0.55).
  Plot interval scaled with stop time so the FTL-probe/plotting cost per eval
  stays ~constant.
- **Done — `scalar_mass` is a searched dimension; fly-away velocities capped**
  (commit `e9c193d`, also recorded in *Future directions* #1 above).
  `grtresna_scalar_mass ∈ [0.3, 1.5]` (default `0.6`) lets the search pick heavy,
  bound matter (binds within the shell width). Toroidal current (warp motor)
  keeps its full `±1.2` range; the net-outflow components are tightened
  (`poloidal ±1.5→±0.8`, `radial ±0.8→±0.3`) so heavy matter stays put. Propagated
  consistently to both the GRTresna solve and the GRTeclyn evolution. Tests:
  `test_scalar_mass_flows_into_grtresna_config`,
  `test_shell_fly_away_velocities_are_capped` (search space now 18D).

### Watch-list (open, not yet verified on `v4`)

- **Heavy-mass convergence.** Larger `m` raises the shared potential
  `½m²(Σφ)²` in the GRTresna elliptic solve, which may shift the convergent
  basin. Watch the pre-GPU rejection rate over the first ~20 evals; if it spikes
  at the high-mass end, narrow the top of the mass range or raise solver
  iterations.
- **Gauge-invariant signal still under-weighted** (carried over from `v2`): the
  null-geodesic reliability gate (`h_quality_ok`) still rejects most candidates,
  so `operational_ftl_geodesic` contributes 0 and the coordinate-speed
  `operational_ftl_solved` dominates `ftl_first`. Unchanged in `v4`.

## Campaign `ftl_discovery_v7` — finished run + critical leaderboard review (2026-06-11)

Run `runs/grtresna_qd/ftl_discovery_v7/` (`speed_super` 8×8, `ftl_first`, 8 GPUs,
`STOP_TIME=8`, `PLOT_INTERVAL=80`, `scalar_mass` searched, de-saturated precursor
+ 3D coherence + rebalanced weights). Finished: 88 candidates, 53 `gpu_ok`,
coverage 0.3125 (20/64 cells), best score 606. **All 20 retained elites are tier
`nontrivial` — zero reached the certified `operational` tier** (no candidate
passed the gauge-invariant geodesic gate).

### Do the metrics work this time? — yes on persistence/coherence, still blind on gauge-invariance

The two bugs that made `v2`/`v3` leaderboards lie are fixed and verified on real
episodes:

- **`survival` is now honest and differentiating.** Top-3 read `survival = 1.0 /
  0.72 / 1.0` (eval 71/79/64) instead of a flat `1.0`. The persistence gate bites
  (eval 79's structure partially decays → 0.72). Frames confirm it: eval 71's
  `Pi_lump_sum` is a strong coherent dipole at `t=8.02` (amplitude held at ~0.01,
  same as `t=0`), **not** the dissipated nothing we saw in `v3`. The heavier
  searched mass (binds within the shell) is doing its job.
- **`structure_coherence = 1.0`** on all top-3, matching the frames (single
  connected blob with an internal dipole, no fragmentation into separated lobes).
  The 3D covering-grid count is no longer fooled by slice orientation.
- **`superluminal_fraction` de-saturated**: top-3 read 0.11 / 0.17 / 0.20 (was
  pinned at 1.0 in `v6`). The descriptor carries real signal again.
- **Score composition is now persistence-honest.** eval 71's 606 breaks down as
  `operational_ftl` (evolved, sustained) **205 pts (34%)** + `ftl_persistence`
  **155 pts (26%)** + `operational_ftl_solved` 106 (17%, itself down-gated to
  0.59) + `channel_progress` 65 + health/survival block ~49. So **~60% comes from
  the *dynamical* evolved+sustained channel**, not from a one-shot coordinate
  speed and not from integration noise. Contrast `v2`, where the top score was
  83% from an `h_quality_ok=False` geodesic artifact.

### What was found

A reproducible **family** of non-stationary, exotic-supported, bound lump
configurations that sustain a *localized* superluminal **coordinate** channel
coherently out to `t=8`. Top-3 (all cell `[7,7]`):

| eval | score | survival | coherence | β_mean | max c (evolved) | superlum. frac | op_ftl | persistence |
|------|-------|----------|-----------|--------|-----------------|----------------|--------|-------------|
| 71   | 606   | 1.00     | 1.0       | 0.444  | 1.108           | 0.111          | 0.513  | 0.517       |
| 79   | 375   | 0.72     | 1.0       | 0.466  | 1.118           | 0.166          | 0.213  | 0.216       |
| 64   | 348   | 1.00     | 1.0       | 0.505  | 1.103           | 0.198          | 0.188  | 0.160       |

All three are genuinely non-stationary (`beta_mean ≈ 0.44–0.51`, real net shift),
stay on-constraint (final Ham/Mom L2 `~1e-3 / 1e-4`), are horizon-free
(`min_theta_plus > 0`), and **require exotic matter** (`wec_violation_fraction
0.76–0.85`, `exotic_penalty −0.53 to −1.0`).

### Critical caveat — these are precursors, NOT certified FTL

`operational_ftl_geodesic = 0` for **every** top candidate: the only
gauge-invariant test (null-ray shortcut) is still rejected as unreliable
(`h_quality_ok=False`, e.g. eval 71 `max H=3.6e-4`, 4/5 rays reached). So the
leaderboard is ranking a sustained, coherent, persistent *coordinate-speed /
shift* channel (`max c ≈ 1.10–1.12` localized near center) — a strong physical
precursor — but nothing here is a proven gauge-invariant shortcut. This is the
same open watch-list item from `v2`/`v4`, now the clear bottleneck: **fixing
null-geodesic reliability is the next high-leverage step**, because until it
passes we cannot promote anything past the `nontrivial` tier.

Secondary note: top scorers cluster in cell `[7,7]` (both descriptor axes near
max), so coverage (0.31) is moderate but the *quality* front is narrow — the
search is exploiting one basin. Worth widening exploration or adding a descriptor
that separates these once the geodesic gate is trustworthy.

## Null-geodesic reliability fix (2026-06-11, post-v7)

The `v7` bottleneck above — `operational_ftl_geodesic = 0` for every elite
because `h_quality_ok=False` — was traced to **two real bugs in the ray tracer**
(`metrics/probes/ftl/geodesic.py`), not to the candidates being non-FTL.

### Bug 1 — rays launched backward (the dominant failure)

The initial ray momentum was built as `k^μ = (1,1,0,0)` (contravariant guess),
converted to covariant `k_μ = g_{μν}k^ν`, then passed through a generic
`project_null`. That projection rescales the spatial momentum to restore the null
condition but **does not pin the propagation direction** — it could (and did)
select the *backward* null root. Instrumenting the central ray of `eval_000071`
showed `dx/dλ` pointing in **−x** from step 0, so the ray left the −x boundary
after 57 steps. Almost no rays reached the detector (`reached 0/5`, `0/5`, `2/5`
on the top-3), so the gate could never pass.

*Fix:* `future_null_cov(g, n_hat)` builds the momentum on the contravariant side
— `k^μ = (k^t, n_hat)` with the null condition solved for the future-directed
`k^t` (the unique positive root, since `g_tt<0` and `n·γ·n>0` give opposite-sign
roots). This guarantees `dx^i/dλ ∝ n_hat = +x` and `dt/dλ > 0`. `project_null`
also gained an optional `dx_ref` so in-flight re-projection keeps the ray's
direction. Result: **all retained elites now reach 5/5**.

### Bug 2 — reliability gated on an unreachable absolute drift

With rays reaching, `h_quality_ok` was still `False` because it tested
`max|H| ≤ 1e-5` (`H = ½ g^{μν}k_μk_ν`, exactly 0 for null). But the metric is
only **C0 (trilinear) interpolated**, so the per-step pre-projection drift has an
interpolation-set floor: a convergence study on `eval_000071` (step
`0.05→0.025→0.0125`, grid `65→97→129`) showed `max|H|` falling only as `~O(ds)`
(`1.2e-3 → 3.2e-4`) and never approaching `1e-5`, while `f_geo` stayed **stable
at 0.0077–0.0081** (the shortcut is a real, converged signal). The absolute
threshold was therefore impossible to satisfy by construction.

*Fix:* gate on the **relative** drift `|H| / ½(|g_tt(k^t)²| + |k_iγ^{ij}k_j|)`
(scale-free "distance off the null cone"), tolerance `H_REL_TOL = 1e-2`. The
v7 elites sit at `5e-4–2e-3` relative drift — an order of magnitude inside the
bar — while the absolute drift and `f_geo` magnitude are still reported.

### Effect (re-scored on the retained v7 plotfiles)

`operational_ftl_geodesic` now fires correctly and **discriminates**:

| eval | f_geo | rel-H | h_ok | geo reward (×1500) |
|------|-------|-------|------|--------------------|
| 88   | 0.0083| 1.4e-3| ✅   | +223 |
| 15   | 0.0078| 1.2e-3| ✅   | +208 |
| 71   | 0.0077| 1.3e-3| ✅   | +206 |
| 79   | 0.0058| 2.0e-3| ✅   | +147 |
| 64   | 0.0047| 9.4e-4| ✅   | +113 |
| 80/81/83/85 | 0.0 | <1.5e-3 | ✅ | +0 (no shortcut) |

The gauge-invariant term now contributes, **re-ranking the board** (eval 88 and
15 leap up next to eval 71) and rewarding genuine null shortcuts while still
zeroing no-shortcut geometries (via `GEO_FTL_FLOOR=1e-3`). Regression tests:
`test_future_null_cov_propagates_forward_under_shift`,
`test_integrate_null_ray_reaches_detector_under_shift`.

### Campaign `ftl_discovery_v8`

Relaunched on the geodesic fix (same search space / descriptor / 8 GPUs,
`STOP_TIME=8`, `PLOT_INTERVAL=80`). Success criterion: at least some elites reach
the **`operational` tier** (certified gauge-invariant shortcut), which `v7` could
never do. Watch whether the search now climbs `f_geo` toward the `5e-2` target
instead of plateauing on the coordinate-speed channel.

## Geodesic-reward recalibration → `ftl_discovery_v9` (2026-06-11)

`v8` worked as intended — the reliability gate fired and **eval 11 reached the
`operational` tier**, the first stable, certified, multi-probe-corroborated
candidate. But it scored **1066**, ~5–18× the rest of the board, which read as a
"breakthrough" when the underlying physics was modest.

### Validation of `v8` eval 11 (score 1066, cell [2,0])

The candidate is genuinely good *and* genuinely small:

- **Real, reliability-gated shortcut.** `f_geo = 0.033`, `h_quality_ok=True`,
  `max_h_rel_drift = 3.2e-4` (≪ `1e-2` tol), **5/5 rays reached** the detector.
- **Corroborated by three independent probes.** geodesic `f_geo = 0.033`, evolved
  coordinate `f_op = 0.032` (`max c = 1.044`), `t=0` solved `f_op = 0.013`,
  persistence `f_op_median = 0.032` sustained across plotfiles — not a single-probe
  artifact.
- **Stable & coherent.** `stability = 0.944`, `structure_coherence = 1.0`; `chi`
  stays a single connected blob from t=0→8 (does not fly away or fragment). The
  local-speed map shows a coherent, localized bump — the persistence/coherence
  fixes are working.
- **But near-flat and modest.** `chi` within 2–3% of 1, `curvature_activity = 0.032`,
  `nontrivial_geometry = 0.040`; the shortcut is a physically small ~3.3%.

### Diagnosis — the reward saturated, so magnitude was lost

The geodesic component is `(f_geo − GEO_FTL_FLOOR)/(GEO_FTL_TARGET − GEO_FTL_FLOOR)`
clamped to 1. With `GEO_FTL_TARGET = 5e-2`, a marginal 3.3% shortcut already mapped
to **0.654**, and at weight ×1500 that is **981 pts (92% of the total)**. A genuinely
dramatic 30% shortcut would have scored *the same*: the scale was saturated and blind
to magnitude.

### Fix

In `metrics/score.py`:

- `GEO_FTL_TARGET`: **`5e-2 → 2e-1`** — full marks now require a dramatic ~20%
  null arrival-time shortcut; the component ramps linearly with magnitude. (This
  corrects the component value itself, so *every* objective mode benefits.)
- `ftl_first` geodesic weight: **×1500 → ×1000** — keeps the trusted gauge-invariant
  term the dominant FTL signal (2.5× the coordinate `operational_ftl` ×400) without
  inflating the scalar.

Effect on eval 11: component `0.654 → 0.161`, geodesic reward `981 → ~161 pts`,
total **`1066 → ~250`** — still clearly #1 over the ~60-pt field, but now reads as
"modest, stable, real" rather than "breakthrough." A future genuine 20% shortcut
lands at ~1000+, reserving big numbers for big physics.

Regression tests updated (`test_geodesic_confirmation_dominates_coordinate_shortcut`,
`test_stationary_warp_lens_artifact_ranks_below_genuine_candidate`) to use a *strong*
`f_geo = 0.2` so "a dramatic shortcut dominates" still holds while a marginal one
intentionally no longer does. 40/40 scoring + geodesic tests pass.

### Campaign `ftl_discovery_v9`

Relaunched on the recalibrated scoring (same search space / descriptor / 8 GPUs,
`STOP_TIME=8`, `PLOT_INTERVAL=80`). Success criterion: the leaderboard now ranks by
shortcut *magnitude* — watch whether the search pushes `f_geo` upward (toward the new
`2e-1` target) rather than parking on the first modest near-flat shortcut.

## `ftl_discovery_v9` review + shaping rebalance → HQ promotion (2026-06-11)

`v9` finished (54/88 evals `gpu_ok`). The headline is **the measurement metrics now
work**, but a *scoring-balance* regression introduced by the geodesic recalibration
let coordinate artifacts out-rank validated shortcuts. Fixed, then the genuine top 3
were promoted to HQ.

### What works (validated in production)

- **The null-geodesic reliability gate fires for every eval.** `h_quality_ok=True`,
  **5/5 rays reach** the detector on all 54 successes — the [`v7` blindness](#null-geodesic-reliability-fix-2026-06-11-post-v7)
  is gone, confirmed at scale (not just on replayed elites).
- **Honest magnitude scaling holds.** `f_geo = 0.033 → 0.16`, `0.026 → 0.12`,
  `0.012 → 0.055` — linear, faithful, no saturation.
- **Five genuine, reliability-certified gauge-invariant shortcuts** were discovered:
  eval 11 (`f_geo=3.3%`), 40 (2.6%), 80 (2.3%), 25 (1.7%), 57 (1.2%). The
  persistence/coherence and gate-notes are all accurate.

### The bug — coordinate precursor out-voted validated FTL

The [`v9` recalibration](#geodesic-reward-recalibration--ftl_discovery_v9-2026-06-11)
deliberately shrank the genuine-FTL reward (a real 3% shortcut → ~160 pts) **but the
coordinate-shaping weights were never trimmed to match**. Net result: the *t=0
coordinate precursor* `operational_ftl_solved` (×180, a "shaping gradient" by design)
out-scored validated shortcuts. The strongest real shortcut ranked **#8**:

| eval | old score | gauge-invariant geodesic | driver |
|------|-----------|--------------------------|--------|
| 40 | 326 | **2.6% (real)** | geodesic + healthy structure |
| 72 | 315 | none (`f_geo=0`) | t=0 coordinate precursor |
| 74 | 312 | none (`f_geo=0`) | t=0 coordinate precursor |
| 13 | 296 | none (flagged *gauge artifact*) | evolved coordinate channel |
| … | | | |
| 80 | 222 | **2.3% (real)** | geodesic |
| 11 | 203 | **3.3% (strongest!)** | geodesic |

Separately, eval 13 was *flagged* `coordinate FTL channel is a gauge artifact
(f_geo=0)` yet still banked `operational_ftl` points — the diagnostic was right but
the score ignored it.

### Fix (`metrics/score.py`, `ftl_first`)

- `channel_progress` **×150 → ×100**, `operational_ftl_solved` **×180 → ×50** —
  coordinate-cone shaping is now clearly subordinate to a realistic gauge-invariant
  shortcut (~160 pts), so a t=0 precursor can no longer out-vote a validated one.
- `operational_ftl` is **zeroed when a *trustworthy* geodesic probe finds no shortcut**
  (`geo_trustworthy and f_geo ≤ GEO_FTL_FLOOR`). The gauge-invariant probe is the
  arbiter; a contradicted coordinate channel is the artifact the note already names.

### Effect (re-scored offline; simulator matches the live scorer to 0.000)

The leaderboard inverts to the physically correct order — **the top 3 are all genuine
gauge-invariant shortcuts**:

| rank | eval | new score | geodesic |
|------|------|-----------|----------|
| 1 | 40 | 266 | **2.6%** |
| 2 | 11 | 200 | **3.3%** |
| 3 | 80 | 185 | **2.3%** |
| 4–7 | 74 / 72 / 15 / 13 | 174–178 | coordinate-only |

Tests: **64/64 pass**. `test_strong_evolved_shortcut_gets_operational_ftl_reward`
threshold `300 → 250` (coordinate-only is *intentionally* worth less now);
`test_geodesic_zero_flags_gauge_artifact` strengthened to assert `operational_ftl == 0`
(locks in the eval-13 demotion).

> The descriptors (BD axes: `speed_tilt`, `superluminal_fraction`) are unchanged, so
> the fix only affects which candidate wins per archive cell and the search gradient —
> meaning a future re-run steers toward genuine shortcuts instead of coordinate tilt.

### HQ promotion of the genuine top 3

With the corrected ranking, the three real shortcuts (eval 40, 11, 80) were promoted
via fresh GRTresna solve + framed GPU evolution at HQ resolution:

```bash
QD_RUN=runs/grtresna_qd/ftl_discovery_v9 NAME_PREFIX=ftl_discovery_v9 \
  CANDIDATES="40 0 11 1 80 2" \
  bash scripts/search/run_promote_qd_batch.sh
# L=128 N=256 max_level=3 t=30, Ham/Mom ≤ 10%, frames on the fly
```

Outputs in `runs/grtresna_promote/l128n256t30_ftl_discovery_v9_qd_eval0000{40,11,80}/`.
Success criterion: at HQ resolution and longer evolution (t=30), the gauge-invariant
`f_geo` survives and ideally *grows* — confirming the 2–3% shortcuts are physical, not
discretization-limited.

## HQ verdict: shortcuts did not survive refinement → `ftl_discovery_v10` (2026-06-11)

All three HQ promotions completed (L=128, N=256, max_level=3, t=30). The success
criterion above was **not met**: the gauge-invariant shortcut collapsed under
continuum refinement in every case.

| eval | QD `f_geo` | **HQ `f_geo`** | geodesic reliable? (HQ) | QD survival | HQ survival | HQ exotic penalty |
|------|-----------|----------------|--------------------------|-------------|-------------|-------------------|
| 11   | 3.30%     | **0.0%**       | ✅ 5/5, h_rel≈3e-4        | 1.00        | 0.64        | −0.83             |
| 40   | 2.58%     | **0.0%**       | ✅ 5/5, h_rel≈7e-4        | 0.99        | 0.61        | −0.92             |
| 80   | 2.33%     | **0.29%**      | ✅ 5/5, h_rel≈1e-3        | 0.56        | **0.18**    | −0.39             |

**Did we discover FTL? No.** The 2–3% shortcuts that topped the `v9` leaderboard
were resolution/gauge artifacts. At the correct resolution they vanish (evals 11,
40) or shrink to a near-noise 0.29% on a *dissipating* structure (eval 80, survival
0.56 → 0.18). Every candidate also still **requires exotic matter** (negative energy
density). The coordinate-channel signals (`operational_ftl_solved`, precursor, shift)
remained non-zero at HQ and were **correctly gated out** by the `v9` rebalance — e.g.
eval 11's `operational_ftl` was zeroed with the note *"trustworthy geodesic probe
found no gauge-invariant shortcut … coordinate channel is a gauge artifact."*

**Why this is the pipeline working, not failing.** The point of the HQ promotion +
reliable geodesic probe is to be a *rejection filter*. It did its job: it killed our
own best candidates under refinement rather than rubber-stamping them. A trustworthy
"no" on a coarse near-miss is the prerequisite for ever trusting a future "yes".

### Two diagnoses → two `v10` changes

1. **The QD window was too short.** Several QD "survivors" (evals 11, 40 at survival
   ≈ 1.0 by t=8) went on to lose geometry and drift at HQ (survival 0.6, large
   areal-radius drift over t=30). The QD loop never saw the second half of the
   evolution, so it promoted late-unstable candidates. **Fix:** extend QD
   `STOP_TIME` from 8 → **16** so the persistence/stability metrics observe the
   later instability and reject it *inside the QD loop* (`PLOT_INTERVAL` scaled
   80 → 160 to hold the ~6-frame post-processing budget).

2. **Only momentum-carrying matter was ever tested.** Every shell candidate carried
   toroidal/poloidal/radial currents + spin; we never sampled purely *static* matter,
   where the only FTL channel is the gauge/shift seed and the geometry itself rather
   than frame-dragging. **Fix:** add a searched `grtresna_shell_static` toggle
   (rounded to int, starts 0 = moving matter). When set, all lump velocities and
   `omega` are forced to zero, so the search explores static lumps alongside the
   moving-matter family.

Code: `grteclyn-wrapper/src/grteclyn_wrapper/search/optimize.py`
(`grtresna_shell_search_space` + `build_grtresna_config` shell branch);
`scripts/search/run_grtresna_qd_search.sh` (`STOP_TIME`, `PLOT_INTERVAL`);
tests in `tests/test_grtresna_shell_ansatz.py` (19-dim space + static-toggle test).

`ftl_discovery_v10` launched with these two changes.

## `ftl_discovery_v10` review → persistence-gate + physicality pressure → `ftl_discovery_v11` (2026-06-12)

`v10` ran the full **400 evals** (240 `gpu_ok` = 60%; the `Mom=nan` static-matter
convergence fix held — static/dynamic split ~40/60, both families evaluating
cleanly). The run converged/stalled: behavior-grid coverage plateaued at 37.5%
and best score froze at 332 for the last 3 windows.

**Top 5 (all dynamic):**

| # | eval | score | `f_geo` | structural_persistence | exotic_fraction | energy_condition |
|---|------|-------|---------|------------------------|-----------------|------------------|
| 1 | 258 | 332 | ~3.3% | **0.46** ⚠ | 0.99 | 0.05 |
| 2 | 283 | 304 | ~2.2% | 0.77 | 0.94 | 0.03 |
| 3 | 335 | 301 | ~1.7% | 0.70 | 0.97 | 0.06 |
| 4 | 356 | 291 | ~2.0% | 0.58 | 1.00 | 0.03 |
| 5 | 301 | 254 | ~3.1% | 0.84 | 0.91 | 0.02 |

**Metrics health:** no bugs. The geodesic recalibration is honest (3% → ~160 pts,
not the old inflated ~1000); gating fires correctly (geodesic-driven #1 has
`operational_ftl_solved=0`; #5's non-persistent shortcut yields `ftl_persistence=0`);
the t=16 window catches late collapse (eval 88 scored **−432** on a horizon +
instability penalty). The static toggle is explored and correctly scores ~0 on
every FTL channel (no currents → no frame-drag → no shift).

**Three problems the elites expose:**

1. **Shortcuts sit in the HQ-death band.** Every top-5 `f_geo` is 1.7–3.3% — the
   exact magnitude that collapsed to ~0 under continuum refinement in the `v9` HQ
   promotions. At QD resolution (N=128) this is where gauge/discretization
   artifacts live.
2. **Transient shortcut on a fragmenting structure ranks #1.** eval 258 banked the
   strongest shortcut while its lump shattered into turbulent lobes
   (`structural_persistence=0.46`, visibly broken by t=16). The geodesic reward was
   *not* persistence-gated, so a disintegrating end-state out-ranked the more
   coherent candidates.
3. **All elites are maximally-exotic warp bubbles.** `exotic_fraction` 0.91–1.00,
   `energy_condition` ~0.02–0.06 across the whole top 5. At the old weights
   (`exotic` 1, `energy` 2) the exotic cost was ~−0.5 against a ~150-pt shortcut, so
   the search had no incentive to leave the easy exotic corner.

### Two metric changes → `v11`

1. **Persistence-gate the geodesic reward.** `operational_ftl_geodesic` is now
   multiplied by `structural_persistence` (the same gate already applied to the
   shaping rewards), so a gauge-invariant shortcut only counts if the structure
   producing it actually holds together. Demotes eval-258-style transients below
   coherent survivors of comparable magnitude.
   (`score.py`, geodesic block.)
2. **Raise physicality pressure.** `energy_condition` weight 2 → **40** and the
   graded `exotic_penalty` weight 1 → **100** in the `ftl_first` objective, so a
   fully-exotic shortcut (~−30..−60 pts) is selected against relative to a cleaner
   one — forcing a real trade between shortcut strength and exotic content. The
   penalty stays graded (0..−1.6), preserving the QD gradient.
   (`score.py`, `ftl_first` weight block.)

Tests: `tests/test_grtresna_integration.py::test_geodesic_reward_gated_by_structural_persistence`
(geodesic reward scales linearly with persistence; fragmenting end-state ranks lower).

`ftl_discovery_v11` launched with these two changes (400 evals, t=16, static toggle).

## `ftl_discovery_v12` review → λφ⁴ + FTL geometry layouts → `ftl_discovery_v13` (2026-06-12)

**v12 verdict (278 evals):** metrics honest, zero `operational_ftl_geodesic` across all
survivors. Top score 130 (eval 197) came entirely from coordinate shaping
(`channel_progress`, `ftl_precursor`, gated `operational_ftl_solved`) with
`f_geo=0` — the geodesic probe correctly rejected every coordinate superluminal
bubble. Archive stuck in descriptor column 7 (max cone tilt). Not promotable.

**Root cause:** matter disperses (no self-binding) and the Fibonacci-sphere layout
only explores one topology. The search was gaming coordinate precursors, not
discovering gauge-invariant shortcuts.

### Three changes → `v13`

1. **Searchable λφ⁴ self-interaction** (`grtresna_scalar_lambda ∈ [0, 0.1]`,
   default 0). Potential on the shared field sum:
   `V = ½(m·Σφ)² − (λ/4)(Σφ)⁴`. Wired identically in GRTresna
   (`scalar_lambda`) and GRTeclyn (`recipe_scalar_lambda`). λ=0 recovers v12.

2. **Searchable matter layout** (`grtresna_matter_layout` 0–3):
   - 0 = sphere (Fibonacci shell, v12 baseline)
   - 1 = channel (directed corridor along the polar axis)
   - 2 = bipolar (two clusters on opposite sides)
   - 3 = ring (planar loop orthogonal to the axis)

   Not Alcubierre/warp-bubble mimicry — diverse FTL-oriented initial geometries.

3. **Geodesic contradiction gate on shaping.** When a trustworthy geodesic probe
   finds `f_geo=0`, all FTL shaping terms (`channel_progress`, `ftl_precursor`,
   `shift_drive`, `operational_ftl_solved`) are zeroed — stops v12-style
   leaderboard inflation.

**Cluster safety:** `run_grtresna_qd_search.sh` runs a pytest preflight gate
(`test_scalar_lambda_potential`, `test_grtresna_shell_ansatz`,
`test_matter_geometry_consistency`, `test_grtresna_integration`) before launch.
Set `SKIP_QD_PREFLIGHT_TESTS=1` to bypass (not recommended).

**Rebuild required:** C++ changes to `GRTresnaScalarPotential.hpp` and GRTresna
`MyMatterFunctions.cpp` need GRTeclyn + GRTresna recompiled before v13 evals.

Search space: **21 dimensions** (was 19).

### Matter-distribution enrichment → per-lump profile + cloud layout

Two follow-on changes broaden the *matter* the search can place. These are
richer initial data only — not warp-bubble mimicry and not tuned to the FTL
score:

4. **Per-lump matter profile.** Each lump carries a `profile` selector:
   `0` = Gaussian envelope `exp(−r²/2w²)` (default, recovers v13), `1` = smoothed
   top-hat "ball" `½(1 − tanh((r − w)/0.25w))` — a near-uniform, more volumetric
   blob with a soft edge. In the shell ansatz it is driven by
   `grtresna_shell_profile_fraction ∈ [0,1]` (+ `grtresna_shell_profile_phase`),
   which assigns the top-hat to a searchable subset of lumps (mirrors the
   `exotic_fraction` idiom); fraction `0` keeps every lump Gaussian. The envelope
   is byte-identical in GRTresna `MatterParams.hpp::lump_phi` (constraint solve +
   painter) and the Python gridinit painter `lump_fields.py`, so the solved
   metric and the evolved φ/Π stay mutually consistent. Gradients and Π adapt
   automatically (finite differences of `lump_phi`).

5. **Quasi-random cloud layout** (`grtresna_matter_layout = 4`). A deterministic
   low-discrepancy (R3) scatter of the ≤5 lumps inside a bounded ball, oriented
   by the searched axis frame — a reproducible, fully asymmetric distribution the
   symmetric sphere/channel/bipolar/ring layouts cannot represent.

**Rebuild required:** GRTresna `ScalarFieldBH` must be recompiled (new `profile`
field in `lump_t`); a stale binary silently ignores it. No GRTeclyn `Source/`
change is needed — the evolution loads φ/Π from the gridinit.

Search space: **23 dimensions** (was 21).

---

## Alcubierre positive control → metric-first vs matter-first verdict (2026-06-12)

**Motivating question:** *why can't the search even find the Alcubierre metric — a
known superluminal solution? Are we doing something wrong?* To answer it without
guessing, we built a **metric-first positive control**: prescribe the textbook
Alcubierre metric analytically and run the *exact same probes* the QD campaign
uses on it.

- Builder + report: `grteclyn-wrapper/scripts/validation/alcubierre_metric_validation.py`
- Permanent regression tests: `grteclyn-wrapper/tests/test_alcubierre_validation.py`

Alcubierre (bubble at origin, moving +x): `α=1`, `γ_ij=δ_ij`,
`β^x = -v_s f(r_s)` with the top-hat shape `f`. A +x null ray then has coordinate
speed `1 + v_s f` (superluminal in the wall), so it should register a shortcut.

### Result: the probes are validated

| config | `f_geo` (shortcut) | rel-H drift | reliability gate | min NEC |
|--------|--------------------|-------------|------------------|---------|
| control `v_s=0` | 0.000 | 0.0 | — | +0.000 |
| `v_s=2`, 65³ (≈ QD res) | **0.315** | 2.16e‑2 | ✗ FAIL | −0.40 |
| `v_s=2`, 97³ | 0.317 | 1.07e‑2 | ✗ FAIL | −0.46 |
| `v_s=2`, 129³ (≈ HQ res) | 0.315 | 5.04e‑3 | ✓ PASS | −0.47 |
| `v_s=2`, 97³ soft wall | 0.335 | 2.50e‑3 | ✓ PASS | −0.11 |
| `v_s=2`, 129³ soft wall | 0.337 | 1.12e‑3 | ✓ PASS | −0.11 |

- **Geodesic probe works.** It measures a stable ~32% gauge-invariant shortcut on
  Alcubierre (`f_geo` converged across 65³→129³), and the flat control stays dark.
- **Energy-condition probe works.** `stress_energy` (`T=G/8π`, already implemented
  in `warpfactory.py`) extracts the required matter and reports **min NEC < 0** —
  the exotic matter Alcubierre demands, correctly flagged.
- **So the campaign's null results are trustworthy:** on a genuine superluminal
  metric our metrics light up; their silence on physical-scalar candidates is real.

### Bonus finding: the QD-resolution reliability gate rejects even Alcubierre

The H-drift gate (`H_REL_TOL=1e-2`) **fails at QD resolution** (65³) on a textbook
Alcubierre bubble and only passes near HQ resolution (129³). The drift is a pure
discretization artifact (halves per refinement) while `f_geo` is unchanged — i.e.
a *genuine* sharp-walled shortcut found by the QD search would have its geodesic
reward **zeroed at QD resolution** and only certified after HQ refinement. This is
exactly the `h_quality_ok=False` we kept hitting in v7–v9.

**Fixed (2026-06-12):** `compute_geodesic_ftl_from_plotfile` now does a
**reliability re-probe**. It traces the ray fan at the cheap base resolution
(65³); if that finds a coordinate shortcut (`f_geo > GEO_REFINE_FLOOR = 1e-3`) but
the H gate fails, it re-traces *once* at `GEO_REFINE_N = 129³` (> 96, the first
grid that reliably certifies even a sharp `σ=2` Alcubierre wall) and returns that
report. The re-probe never fires for no-shortcut or already-reliable candidates,
so the extra cost is paid only on the rare warp-candidate that would otherwise be
silently discarded. Regression tests live in `tests/test_null_geodesic.py`
(`test_reliability_reprobe_*`).

### Verdict: is metric-first (geometry-first) a bad idea, or are we doing matter-first wrong?

**Neither philosophy is "wrong" — they answer different questions, and we are
already using both correctly. The distinction is *analysis* vs *dynamics*:**

- **Metric-first as an *analysis / validation* tool — good, and we now rely on it.**
  Prescribing a metric and reading off `T=G/8π` is fully implemented
  (`stress_energy`) and the geodesic probe is matter-agnostic. This is how we just
  *validated* the FTL metrics against Alcubierre. Keep it.
- **Metric-first as a *discovery engine* — dead end.** You cannot evolve a
  prescribed metric self-consistently with physical matter. GRTeclyn evolves
  *matter fields* (scalar `φ,Π` via Klein–Gordon) and *computes* `T` from them; it
  cannot be handed an arbitrary prescribed `T_μν(t)`. Holding the metric at the
  Alcubierre form forever means evolving a **fixed background** — kinematics, not
  physics — and its "FTL" is baked into the coordinates. That is the very thing
  this project was built to avoid.
- **The Alcubierre matter is not the stress-energy of any evolvable field.** Its
  `T` is **NEC-violating**; a canonical scalar field's `T` *always* satisfies the
  NEC. So no scalar configuration reproduces it even at `t=0`, and even a forced
  match would be dragged away instantly by the field's own dynamics. This is why
  Alcubierre needs the "by-fiat" frozen geometry.
- **We are *not* doing matter-first wrong.** Matter-first is the only route to a
  self-consistent, evolvable spacetime. The campaign's failure to find a surviving
  shortcut from physical matter is an **honest physics result** — consistent with
  the NEC / Ford–Roman quantum-inequality no-go theorems — not a pipeline bug.
- **The one genuine gap was instrumentation, not philosophy:** the QD-resolution
  geodesic reliability gate was too strict for sharp walls (shown above). Fixed
  with a 129³ mid-res re-probe (see above); the matter-first design is unchanged.

**Bottom line:** the search "can't find Alcubierre" for three correct reasons —
(1) physical scalar matter cannot source its exotic `T`, (2) the moving-puncture
gauge damps any seeded warp shift, and (3) we deliberately penalize exotic matter.
All three are features. The metric-first control proves our detectors are sound;
the matter-first search is asking the harder, meaningful question (*does a
**physical, evolvable** FTL exist?*), and a clean null is a real answer.
