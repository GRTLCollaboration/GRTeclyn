# MAP-Elites + CMA-ES FTL Discovery — Matter-First Metric Discovery

> Three-stage pipeline: **MAP-Elites** (wide survey) finds where good warps live in
> a matter-sector-specific search space (**23-D shell** for dispersing real scalars,
> **7-D boson star** for complex U(1) matter); **CMA-ES** (local refinement)
> hill-climbs around the best *healthy* survivors; **HQ promotion** re-runs the top
> QD + CMA-ES elites at full resolution and extended time with incremental scoring.
> All stages share the same matter-first loop — propose matter → GRTresna constraint
> solve → GRTeclyn GPU evolution → time-resolved FTL probes → score — but differ in
> proposer, matter sector, resolution, and stop time.

## Known limitations of the lump-based matter ansatz

The shell/ring search spaces parameterize matter as **N identical Gaussian lumps** placed
at deterministic lattice positions. This ansatz was designed for dimensionality reduction
(22-D shell vs ~50-D per-lump), but it imposes structural constraints that **limit the
geometric freedom available to the optimizer** and may prevent discovery of complex warp
geometries. See also the [static shell post-mortem](#static-shell-post-mortem-scalar_shell_ftl_rl_v1-2026-06-24)
for how these limitations interact with the `shell_static=1` pin.

### What is constrained

All N=5 lumps share one global value for each property — the optimizer cannot express
per-lump variation beyond low-order Legendre modulation (ℓ ≤ 2):

| Per-lump property | Independent per lump? | What actually happens |
|---|---|---|
| **Position** | No — deterministic Fibonacci lattice | Only `radius`, `thickness`, `axis_theta/phi` move all lumps together |
| **Amplitude** | No — one `shell_amp` for all | Dipole (ℓ=1) + quadrupole (ℓ=2) Legendre modulation only |
| **Width** | No — one `shell_width` for all | Same ℓ=2 modulation, clamped \[1.0, 6.0\] |
| **Velocity** | No — one (v_tor, v_pol, v_rad) for all | Every lump gets the same velocity decomposed into its local frame |
| **Omega (spin)** | No — one `shell_omega` for all | Same rigid rotation for every lump |
| **Exotic flag** | Coarse wedge — `exotic_fraction/phase` | Binary 0/1, contiguous azimuthal wedge, not independently searchable |
| **Profile** | Same coarse wedge mechanism | `profile_fraction/phase` selects Gaussian vs top-hat |
| **Azimuthal mode** | No — one `shell_mode ∈ {0,1,2}` | All lumps get the same angular harmonic |

### What cannot be expressed

Configurations the current ansatz **structurally cannot represent**:

1. **Asymmetric matter distributions** — e.g. two heavy lumps near the axis + three
   light ones at the equator. Only ℓ ≤ 2 modulation exists.
2. **Per-lump differential motion** — e.g. one lump orbiting fast, another stationary,
   a third falling inward. Critical for counter-rotating pairs and differential-ω
   rotating warp mechanisms.
3. **Nested / multi-shell structures** — inner exotic shell + outer canonical shell at
   a different radius.
4. **Non-Gaussian density profiles** — torus cross-sections, saddle-shaped distributions,
   ring-within-ring, thick vs thin shell substructures.
5. **Arbitrary lump positions** — cannot cluster matter where the geometry needs it
   (e.g. at a warp channel throat). Positions are deterministic lattice functions of
   layout + radius + thickness.

### Proposed alternatives for richer geometric search

The existing infrastructure already supports per-lump keys
([`config.py`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/config.py):
`grtresna_lump{k}_{param}`) so per-lump freedom is structurally available —
the shell ansatz deliberately collapses it for dimensionality reduction.

#### Option A — Neural field decoder (learned geometry)

Replace the lump-painting step with a small coordinate-conditioned neural network
(e.g. MLP or low-rank grid decoder) that outputs `φ(x,y,z)` and `Π(x,y,z)` directly
on the GRTresna grid. The network's latent vector (16–32 D) is the genome that
MAP-Elites/CMA-ES searches; the decoder provides continuous, unconstrained geometric
freedom.

| Property | Detail |
|----------|--------|
| Search dims | **16–32 D** latent (vs 22-D shell, but richer output) |
| Geometric freedom | Arbitrary: toruses, nested shells, saddles, asymmetric clusters |
| Integration point | Replaces lump painting → feeds `gridinit` to GRTresna constraint solve |
| Constraint safety | GRTresna still validates — non-physical proposals rejected as today |
| Training | Pre-train decoder on existing elite `gridinit` fields; fine-tune during search |
| Risk | Decoder must produce constraint-solvable fields (not guaranteed without training) |

#### Option B — Per-lump sequential RL placement

An RL agent places lumps one at a time with per-lump geometric freedom:

| MDP element | Content |
|-------------|---------|
| **State** | Current matter distribution (partial field + GRTresna constraint residual) |
| **Action** | Next lump's (position_xyz, amp, width, velocity_xyz, omega, exotic) |
| **Dense reward** | Cheap GRTresna pre-solve Ham/Mom residual after each lump placement |
| **Terminal reward** | Full `general_ftl` score after GRTeclyn evolution |

This gives per-lump geometric freedom and uses constraint-solve feedback to learn
where matter should go. Avoids the ~50-D curse because the agent learns a placement
*policy*, not a flat vector — generalization across placement steps reduces the
effective dimensionality.

| Property | Detail |
|----------|--------|
| Search dims | ~10 D per action step × N steps (but policy generalizes) |
| Geometric freedom | Full per-lump control: position, size, motion, exotic |
| Constraint feedback | Dense signal from GRTresna after each placement |
| Risk | Sequential placement imposes an arbitrary ordering; may need permutation-invariant architecture |

#### Option C — Spherical harmonic field expansion

Parameterize the matter density as a truncated spherical harmonic expansion with
radial basis functions:

```
φ(r,θ,ϕ) = Σ_ℓm a_ℓm f_ℓ(r) Y_ℓm(θ,ϕ)
```

Search over the coefficients `a_ℓm`. Naturally smooth, rotationally natural, and
the number of DOF scales as `(ℓ_max + 1)²` — 25 D for ℓ_max = 4, 49 D for
ℓ_max = 6.

| Property | Detail |
|----------|--------|
| Search dims | **(ℓ_max+1)²** — 25 D (ℓ=4) to 49 D (ℓ=6) |
| Geometric freedom | All angular structure up to ℓ_max; higher ℓ = finer spatial detail |
| Integration point | Coefficients → evaluate on grid → paint `gridinit` → GRTresna solve |
| No training needed | Purely physics-based; no neural network required |
| Risk | Radial basis choice matters; may need ~3–5 radial DOF per ℓ for multi-shell support |

### Assessment

The current lump ansatz is likely **too rigid** for general FTL discovery — it can only
express highly symmetric, globally uniform matter distributions. The best FTL found
([eval 166](#paired-shell-ftl-comparison-boson-vs-scalar-2026-06-23), score 869) is a
static exotic lens with weak 4D signal (0.10), not a genuine dynamical warp. Enabling
dynamics (`shell_static=0`) is the cheapest immediate fix, but even with dynamics on,
the uniform-lump constraint still prevents asymmetric, multi-scale, differentially-moving
structures that real warp geometries may require.

**Recommended path:** start with Option C (spherical harmonics) — no ML training needed,
physically principled, and 25–49 D is within reach of MAP-Elites with pre-GPU filtering.
Option A (neural decoder) is the highest-ceiling approach but requires training data and
integration work. Option B (sequential RL placement) gives the most natural per-lump
freedom and a defensible role for RL in this pipeline.

### Option C implementation — `sh` ansatz (2025-06-24)

The spherical-harmonic ansatz is now implemented as `--grtresna-ansatz sh`.

**Architecture.** 24 real-SH modulation coefficients (ℓ_max=4, monopole absorbed
into base amplitude) multiplicatively modulate the per-lump amplitude on a Fibonacci
sphere:

```
A_k = base_amp × max(0, 1 + Σ_{idx>0} c_idx Y_ℓm(θ_k, φ_k))
```

Lump positions, velocity decomposition (toroidal/poloidal/radial), and exotic-wedge
selection reuse the shell ansatz infrastructure. Dynamics are ON by default
(`sh_static` starts at 0), fixing the warp-motor issue identified in the static
post-mortem.

**Search space (38 D):**

| Block | Dims | Keys |
|-------|------|------|
| Base amplitude | 1 | `grtresna_sh_amp` [0.06, 0.22] |
| SH modulation | 24 | `grtresna_sh_c1`..`c24` [−0.8, 0.8] |
| Geometry | 3 | `sh_radius`, `sh_width`, `sh_thickness` |
| Physics | 2 | `scalar_mass`, `scalar_lambda` |
| Kinematics | 4 | `sh_toroidal_velocity`, `sh_poloidal_velocity`, `sh_radial_velocity`, `sh_omega` |
| Exotic | 2 | `sh_exotic_fraction`, `sh_exotic_phase` |
| Gauge/toggle | 2 | `shift_seed`, `sh_static` |

**Files added / modified:**

| File | Role |
|------|------|
| `grtresna/sh_fields.py` (new) | Real-SH evaluation, modulation helpers |
| `search/optimize/spaces.py` | `grtresna_sh_search_space()` + `build_search_space` dispatch |
| `search/optimize/config.py` | `_expand_sh_lumps_from_overrides()` — SH → lump expansion |
| `cli/parser.py` | `"sh"` added to ansatz choices |
| `cli/grtresna_context.py` | `"sh"` ansatz context wiring |
| `search/optimize/__init__.py` | Exports |
| `tests/grtresna/test_grtresna_sh_ansatz.py` (new) | 19 tests (SH math, search space, config expansion, CLI) |

**Smoke test (scalar_sh_ftl_smoke_v1, 8 evals):** 1 `gpu_ok` (score 12.0, full pipeline),
5 `grtresna_rejected`, 1 `postload_rejected`, 1 `grtresna_failed`. The successful eval
used all 24 SH coefficients, producing lumps with amplitudes 0.018–0.088 (3.5× range),
confirming per-lump variation from the SH expansion.

**Launch:**
```bash
cd grteclyn-wrapper
GRTRESNA_ANSATZ=sh PIN_DIMS="" \
QD_NAME=scalar_sh_ftl_v1 QD_TARGET_EVALS=200 \
GPU_IDS="0 1 2 3 4 5 6 7" \
  bash scripts/campaigns/qd/run.sh
```

---

## Quick start — running campaigns

All commands run from **`grteclyn-wrapper/`**. Binaries must be built first — see
[Building the binaries](#building-the-binaries-grtresna--grteclyn).

| Stage | What it does | Launcher | Output directory |
|-------|--------------|----------|------------------|
| **0 — QD** | MAP-Elites survey (8×8 archive) | `scripts/campaigns/qd/run.sh` | `runs/grtresna_qd/<name>/` |
| **1 — CMA-ES** | Local hill-climb from QD elites | `scripts/campaigns/cmaes/run.sh` | `runs/grtresna_cmaes/<name>/` |
| **2 — HQ** | Full-res replay + frames | `scripts/campaigns/hq/run_batch.sh` | `runs/grtresna_promote/<prefix>_hq_eval*/` |

Shared search defaults (grid, gates, 4D probe, objective) live in
`scripts/campaigns/lib/search_common.sh`. HQ defaults in `scripts/campaigns/lib/promote_common.sh`.

### Stage 0 — MAP-Elites (QD)

**Generic launch** (background, 8 GPUs, pipelined GRTresna + GPU):

```bash
cd grteclyn-wrapper
QD_NAME=my_campaign_v1 \
QD_TARGET_EVALS=200 \
OBJECTIVE_MODE=ftl_first \
GPU_IDS="0 1 2 3 4 5 6 7" \
GPU_SLOTS_PER_DEVICE=1 \
MAX_CONCURRENT_GRTRESNA=3 \
  nohup bash scripts/campaigns/qd/run.sh \
  > ../runs/my_campaign_v1.launch.log 2>&1 &
```

**`general_ftl` wormhole** (current production path — pins 15-D subspace, `OBJECTIVE_MODE=general_ftl`,
dirs `x y z`, pre-GPU learning):

```bash
cd grteclyn-wrapper
BRANCH=wormhole \
QD_NAME=general_ftl_wormhole_v22 \
QD_TARGET_EVALS=200 \
GPU_IDS="0 1 2 3 4 5 6 7" \
GPU_SLOTS_PER_DEVICE=1 \
MAX_CONCURRENT_GRTRESNA=3 \
  nohup bash scripts/campaigns/general_ftl/run_all.sh \
  > ../runs/general_ftl_wormhole_v22.launch.log 2>&1 &
```

**Resume** an existing QD run (same `QD_NAME`, same dir):

```bash
cd grteclyn-wrapper
QD_NAME=general_ftl_wormhole_v21 QD_RESUME=1 \
  bash scripts/campaigns/qd/run.sh
```

**Boson star** (complex scalar / U(1) matter, **7-D** search — Phase 1 single centered
Gaussian; geometry ansatz unchanged from default shell but matter sector is orthogonal):

```bash
cd grteclyn-wrapper
QD_NAME=boson_star_v1 \
QD_TARGET_EVALS=80 \
GRTRESNA_MATTER_SECTOR=boson_star \
GRTRESNA_MATTER_COUPLING=canonical \
GPU_IDS="0 1 2 3" GPU_SLOTS_PER_DEVICE=1 \
  bash scripts/campaigns/boson_star/run.sh
```

Shortcut launcher sets `GRTRESNA_MATTER_SECTOR=boson_star`, `OBJECTIVE_MODE=ftl_first`,
and `GRTRESNA_FULL_Z=1`. CMA-ES / HQ: `scripts/campaigns/boson_star/{cmaes_run,hq_run}.sh`.
Exotic (phantom) boson: `GRTRESNA_MATTER_COUPLING=exotic` (`scalar_sign=-1`).

**Bosonic shell + FTL (RL chassis)** — same **BosonStarBH** superposed shell as
[grlab splash](../grlab/README.md) (`grtresna_complex_scalar`), but scored with
**`general_ftl`** + 4D geodesic (not `critical_collapse`). Search space =
`grtresna_boson_shell_search_space()` (~18-D shell geometry + boson mass/λ/ω + exotic wedge).
**Exotic wedge search ON** by default — see [Exotic matter](#exotic-matter-by-sector) below.

```bash
cd grteclyn-wrapper
uv sync   # h5py>=3.10 required for GRTresna Chombo→gridinit bridge

QD_NAME=boson_shell_ftl_rl_v1 \
QD_TARGET_EVALS=200 \
QD_ITERATIONS=30 \
STOP_TIME=16.0 \
PLOT_INTERVAL=320 \
GRTECLYN_FRAMES=1 \
GPU_IDS="0 1 2 3 4 5 6 7" \
GPU_SLOTS_PER_DEVICE=1 \
MAX_CONCURRENT_GRTRESNA=5 \
BATCH_SIZE=8 \
  nohup bash scripts/campaigns/boson_star/ftl_shell_run.sh \
  > ../runs/boson_shell_ftl_rl_v1.launch.log 2>&1 &
```

Launcher: `scripts/campaigns/boson_star/ftl_shell_run.sh`. Pins static shell only
(`grtresna_shell_static=1`); **exotic wedge search ON** (`GRTRESNA_BOSON_ALLOW_EXOTIC=1`).
**`PLOT_INTERVAL=320`** → **6** plotfiles over t=16 (dt≈0.01). Output: `runs/grtresna_qd/boson_shell_ftl_rl_v1/`.
RL handoff: [research/RL/LabJournal.md](../RL/LabJournal.md).

| Knob | Boson shell FTL default | Notes |
|------|-------------------------|-------|
| Matter | `boson_star` + `shell` | `grtresna_complex_scalar` |
| Coupling | **`canonical`** | Per-lump exotic via **`shell_exotic_fraction/phase`** |
| Objective | `general_ftl` | dirs `x y z`, 4D `search` profile |
| Exotic search dims | **`shell_exotic_fraction`, `shell_exotic_phase`** | **`GRTRESNA_BOSON_ALLOW_EXOTIC=1`** (launcher default) |
| Frames | **on** (`GRTECLYN_FRAMES=1`) | projections `scalar_activity`, `phi` |
| Stop / plots | t=16, interval 320 | 6 dumps per eval |

| Knob | Default (search) | Notes |
|------|------------------|-------|
| Grid | N=128, L=64, ml=1 | GRTresna solve on 128³ domain |
| Stop time | t=16 | `STOP_TIME=16.0` |
| Archive | 8×8 bins | `BINS=8`, descriptor `ftl_lifetime` |
| Frames | **off** | `GRTECLYN_FRAMES=0` (speed) |
| 4D geodesic | `search` profile | cheap stack for leaderboard |
| Prune eval dirs | top 3 + FTL peaks | `QD_KEEP_TOP_EVAL_DIRS=3` |

**Outputs:** `trajectory.jsonl`, `archive.json`, `eval_*/`, `ftl_champions.json`.

**Monitor:**

```bash
tail -f runs/grtresna_qd/<name>/trajectory.jsonl
cat runs/grtresna_qd/<name>/ftl_champions.json
```

### Stage 1 — CMA-ES optimization

Warm-start from a QD (or prior CMA-ES) trajectory. **Must match** the source campaign's
`OBJECTIVE_MODE`, grid, `STOP_TIME`, pins, and 4D profile — do not switch objectives mid-handoff.

**Wormhole example** (seed QD score record, local σ=0.05 refinement):

```bash
cd grteclyn-wrapper
RUN_NAME=general_ftl_wormhole_cmaes_v1 \
OBJECTIVE_MODE=general_ftl \
WARM_START_TRAJECTORY="${GRTECLYN_ROOT}/runs/grtresna_qd/general_ftl_wormhole_v21/trajectory.jsonl" \
WARM_START_TOP_K=1 WARM_START_JITTER=0.05 SIGMA0=0.05 \
TARGET_EVALS=150 MAX_GENERATIONS=50 KEEP_TOP_EVAL_DIRS=3 \
GPU_IDS="0 1 2 3 4 5 6 7" GPU_SLOTS_PER_DEVICE=1 MAX_CONCURRENT_GRTRESNA=3 \
PIN_DIMS="$(bash -c 'source scripts/campaigns/lib/general_ftl_pins.sh && ftl_general_ftl_wormhole_pins')" \
  nohup bash scripts/campaigns/cmaes/run.sh \
  > ../runs/general_ftl_wormhole_cmaes_v1.launch.log 2>&1 &
```

| Knob | Typical | Notes |
|------|---------|-------|
| Population | = GPU count | `POPULATION` defaults to `#GPU_IDS` |
| σ₀ | 0.05–0.08 | local basin width |
| Warm-start | top-K elites | `WARM_START_TOP_K`, `WARM_START_JITTER` |
| Target | eval budget | `TARGET_EVALS` or `MAX_GENERATIONS × pop` |

**Outputs:** `runs/grtresna_cmaes/<RUN_NAME>/` — same layout as QD (`trajectory.jsonl`, `eval_*/`).

**Monitor:** `tail -f runs/grtresna_cmaes/<RUN_NAME>/trajectory.jsonl`

### Stage 2 — HQ promotion

Replays elite genomes at **N=256, L=128, ml=3, t=30** with fresh GRTresna solve,
**frames on**, 4D geodesic in **`hq`** verify mode, incremental `score_timeseries.jsonl`.

**Single candidate** (`CANDIDATES` = eval/gpu **pairs**):

```bash
cd grteclyn-wrapper
SOURCE_RUN="${GRTECLYN_ROOT}/runs/grtresna_cmaes/general_ftl_wormhole_cmaes_v1" \
NAME_PREFIX=general_ftl_wormhole_cmaes_v1 \
OBJECTIVE_MODE=general_ftl \
CANDIDATES="46 0" \
  bash scripts/campaigns/hq/run_batch.sh
```

**Top-K auto-pick** from trajectory (score-sorted `gpu_ok` rows):

```bash
cd grteclyn-wrapper
SOURCE_RUN="${GRTECLYN_ROOT}/runs/grtresna_cmaes/general_ftl_wormhole_cmaes_v1" \
NAME_PREFIX=general_ftl_wormhole_cmaes_v1 \
OBJECTIVE_MODE=general_ftl \
TOP_K=1 \
  bash scripts/campaigns/hq/run_batch.sh
```

| Knob | Search (QD/CMA-ES) | HQ |
|------|-------------------|-----|
| Grid | 128³, L=64, ml=1–2 | **256³, L=128, ml=3** |
| Stop time | 16 | **30** |
| 4D mode | `search` | **`hq`** (full stack) |
| Geodesic dirs | `x y z` (`general_ftl`) | **`x y z`** (same) |
| Frames | off | **on** (`GRTECLYN_FRAMES=1`) |
| Objective | campaign-specific | match source (`general_ftl` for wormhole) |

**Monitor:**

```bash
tail -f runs/grtresna_promote/<name>.log
tail -f runs/grtresna_promote/<name>/small_data/score_timeseries.jsonl
ls runs/grtresna_promote/<name>/frames/
bash scripts/plot/make_movies.sh runs/grtresna_promote/<name> --framerate 10
# → writes runs/grtresna_promote/<name>/movies/movie_<field>_<axis>.mp4
```

**Incremental scoring note:** mid-run HQ totals are not comparable to search finals until the
end-of-run 4D trace completes — only `ftl_geo_evolving` earns geodesic credit incrementally.

### Rules (do not skip)

1. **CMA-ES must mirror QD** — same `OBJECTIVE_MODE`, pins, grid, `STOP_TIME`, and geodesic config.
2. **`general_ftl` needs `GRTECLYN_GEO_DIRECTIONS=x y z`** — wormhole shortcuts live on **z**; x-only
   scoring replays elites at the wrong fitness (see [v22 CMA-ES](#v22-cma-es-wormhole-refinement-general_ftl_wormhole_cmaes_v1-2026-06-18)).
3. **HQ `CANDIDATES` is eval/gpu pairs** — e.g. `"46 0 39 1"` not a bare eval list.
4. **Search turns frames off, HQ turns them on** — by design (`search_common.sh` vs `promote_common.sh`).

Campaign results and tuning history → [campaign log](#campaign-log--runs-analysis).

### End-to-end orchestrator (one folder)

Runs **QD → CMA-ES → HQ** sequentially under `runs/campaigns/<CAMPAIGN_NAME>/`:

```bash
cd grteclyn-wrapper
CAMPAIGN_NAME=general_ftl_wormhole_v22 \
  bash scripts/campaigns/run_full_campaign.sh
```

See `scripts/campaigns/run_full_campaign.sh` and `scripts/campaigns/README.md`.

## Table of contents

**Quick reference**

- [Quick start — running campaigns](#quick-start--running-campaigns)

**Reference**

- [The idea: matter-first, not metric-first](#the-idea-matter-first-not-metric-first)
- [Matter selector: scalar vs boson star](#matter-selector-scalar-vs-boson-star)
- [The pipeline](#the-pipeline)
- [The hard consistency rule](#the-hard-consistency-rule)
- [Behavior descriptors](#behavior-descriptors-the-diversity-axes)
- [Scoring model](#scoring-model-the-quality-axis)
- [Code map](#code-map-where-everything-lives)
- [Building the binaries](#building-the-binaries-grtresna--grteclyn)
- [How to run a campaign](#how-to-run-a-campaign)

**Campaign log** (reverse-chronological)

- [Quick index](#campaign-log--runs-analysis)
- [Boson star: unpinned QD + frames (v3)](#boson-star-unpinned-qd--frames-v3-2026-06-18)
- [Boson star: integration + GPU smoke](#boson-star-matter-sector-integration--gpu-smoke-2026-06-18)
- [v22: Pre-GPU rejection learning + v21 resume](#v22-pre-gpu-rejection-learning--v21-resume-2026-06-18)
  - [Final results: top 3 + FTL champions](#v22-final-results-top-3--ftl-champions)
  - [Physical validation & next steps](#v22-physical-validation--next-steps)
- [v22 CMA-ES: wormhole refinement](#v22-cma-es-wormhole-refinement-general_ftl_wormhole_cmaes_v1-2026-06-18)
  - [HQ eval 046: final results (t=30)](#hq-eval-046-final-results-t30)
  - [HQ eval 046: matter dynamics](#hq-eval-046-matter-dynamics)
- [v21: Pipelined QD + GPU tenancy tuning](#v21-pipelined-qd--gpu-tenancy-tuning-2026-06-17)
- [v10–v20: pipeline evolution & runs](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17)
- [Foundational entries (2026-06-10)](#foundational-entries-2026-06-10)

## The idea: matter-first, not metric-first

Classic warp-drive analysis is **metric-first**: write down a target metric
(Alcubierre, Natário, …), then read off the stress-energy `T_ab = G_ab/8π` it
requires — which is always exotic, often pathological, and never solved as a
*consistent initial-value problem*. The "FTL" in those constructions is baked
into the chosen coordinates and need not survive a real evolution.

This project inverts that. The pipeline is **matter-first**:

1. The search proposes a **matter configuration** — either **multi-lump real scalars**
   (shell / ring / free ansätze) or a **single centered boson star** (complex scalar).
2. GRTresna **solves the Einstein constraint equations** for conformal factor and shift.
3. GRTeclyn **evolves** that spacetime forward on GPU.
4. Probes measure whether an FTL signature emerges and *persists*.
5. MAP-Elites stores the best candidate per behavior cell and mutates elites.

The payoff: any FTL signal that survives this loop is a property of a
**self-consistent, evolved spacetime**, not a hand-picked metric. Metrics are
iteratively hardened so the leaderboard cannot be gamed by coordinate artifacts
(see the [campaign log](#campaign-log--runs-analysis)).

## Matter selector: scalar vs boson star

The pipeline now selects matter **orthogonally** to the geometry ansatz (shell / ring / free).
Default remains the existing **real-scalar lump** path; boson star is opt-in.

| Knob | Values | Default |
|------|--------|---------|
| `GRTRESNA_MATTER_SECTOR` / `--grtresna-matter-sector` | `scalar`, `boson_star` | `scalar` |
| `GRTRESNA_MATTER_COUPLING` / `--grtresna-matter-coupling` | `canonical`, `exotic` | `canonical` |

| Mode | Matter model | GRTresna example | Search space | Notes |
|------|--------------|------------------|--------------|-------|
| **scalar + canonical** | `grtresna_independent_scalars` | `ScalarFieldBH` | 23-D shell (etc.) | Production wormhole path — unchanged |
| **scalar + exotic** | same, all lumps `exotic=1` | `ScalarFieldBH` | 23-D shell | Maximal-slicing exotic path |
| **boson_star + canonical** | `grtresna_complex_scalar` | `BosonStarBH` | **7-D** (`spaces.py`) | U(1) charge conserved; `scalar_sign=+1` |
| **boson_star + exotic** | `grtresna_complex_scalar` | `BosonStarBH` | **7-D** | Phantom boson; `scalar_sign=-1` |

**Exotic matter by sector**

| Sector / launcher | Exotic in search? | How to enable |
|-------------------|-------------------|---------------|
| **Real scalar shell** (`qd/run.sh`, `general_ftl/`) | **Yes** — `grtresna_shell_exotic_fraction`, `grtresna_lump{k}_exotic` | Default search dims; wormhole QD explores exotic mixes |
| **Boson shell FTL** (`boson_star/ftl_shell_run.sh`) | **Yes** — `grtresna_shell_exotic_fraction`, `grtresna_shell_exotic_phase` | **`GRTRESNA_BOSON_ALLOW_EXOTIC=1`** (default in launcher); per-lump exotic via shell wedge |
| **Boson splash** (`splash/run.sh`) | **No** — same boson-shell space + `PIN_DIMS` sign=1 | `critical_collapse`, not FTL |
| **Centered 7-D boson** (`boson_star/run.sh`) | **Optional global sign** — dim exists but often pinned | `GRTRESNA_MATTER_COUPLING=exotic` → phantom boson (`scalar_sign=-1`); **not** per-lump exotic fraction |

**`boson_shell_ftl_rl_v1`** uses **canonical coupling** with **searchable exotic shell wedges**
(`GRTRESNA_BOSON_ALLOW_EXOTIC=1` in `ftl_shell_run.sh`). To disable wedge search and pin
canonical-only (splash-style):

```bash
GRTRESNA_BOSON_ALLOW_EXOTIC=0 \
PIN_DIMS="grtresna_scalar_sign=1 grtresna_shell_static=1" \
  bash scripts/campaigns/boson_star/ftl_shell_run.sh
```

Global phantom boson (whole-field `scalar_sign=-1`, no exotic-fraction dims):

```bash
GRTRESNA_MATTER_COUPLING=exotic \
PIN_DIMS="grtresna_shell_static=1" \
  bash scripts/campaigns/boson_star/ftl_shell_run.sh
```

**7-D boson search dimensions:** `grtresna_scalar_mass`, `grtresna_scalar_lambda`,
`grtresna_bs_phi_c`, `grtresna_bs_profile_width`, `grtresna_bs_omega` (pinned 0),
`grtresna_scalar_sign`, `grtresna_shift_seed`.

**Phase 1 limitation:** boson star is a **single centered Gaussian** — it does **not** place
boson stars at shell lump positions. Shell geometry ansätze and boson matter are separate
search spaces, not combined. Multi-site boson painting is future work.

**Why boson star:** real scalar lumps **disperse** at \(t \gg 50\); complex fields with
conserved U(1) charge (`Pi2 = -\omega\phi/\alpha`) persist under evolution — the motivation
for [action plan step 3](#iv-strategic-action-plan) below.

**Wrapper wiring:** `grteclyn_wrapper/grtresna/matter_models.py` (selection + overrides),
`matter_wiring.py` (`.matter.json`), `search/optimize/spaces.py` (`grtresna_boson_star_search_space`),
`scripts/campaigns/lib/search_common.sh` (frame fields: `phi2`, `Pi2`, `scalar_activity`).

## The pipeline

### Diagram — end-to-end overview

The closed loop: search proposes matter → physics builds and evolves a real
spacetime → metrics discover FTL signatures → archive feeds the next proposal.

![MAP-Elites end-to-end overview](mapelites-end-to-end.svg)

### Diagram — matter-first vs metric-first

![Matter-first vs metric-first diagram](mapelites-matter-first.svg)

### Stage 0 — Quality-Diversity proposer (MAP-Elites)

MAP-Elites maintains an 8×8 behavior archive keyed by a 2-D descriptor. Each
batch it either **mutates an existing elite** (boundary-reflected Gaussian
perturbation, σ=0.15, ~85% of draws) or **samples inside the feasible box** of
known elites. Proposals are points in the active search space — **23-D shell** (default
real scalars) or **7-D boson star** when `GRTRESNA_MATTER_SECTOR=boson_star`.

> **QD vs CMA-ES.** QD (`campaigns/qd/run.sh`) illuminates the 8×8 grid.
> CMA-ES (`campaigns/cmaes/run.sh` → `optimize`) hill-climbs from QD survivors.
> Same search space (`search/optimize/spaces.py`); only the proposer differs.
> See [CMA-ES refinement](#cma-es-refinement-after-map-elites).

### Stage 1 — Initial data (GRTresna, CPU/MPI)

`../GRTresna` paints lumps into `(φ, Π)` and runs a Lichnerowicz/York **constraint
solve**. CPU/MPI-bound; poorly-conditioned proposals fail the Ham/Mom gate and are
**rejected before the GPU**.

### Stage 2 — Evolution (GRTeclyn, GPU)

Converged initial data evolves via `RadialRecipe` to `STOP_TIME`, dumping plotfiles
every `PLOT_INTERVAL` steps.

### Stage 3 — Metrics & probes (scoring)

Plotfiles → diagnostics: constraints, `theta_plus`, comoving/shift stats, matter
density, and FTL probes — coordinate (`operational_ftl_solved`), evolved
(`operational_ftl`, `ftl_persistence`), gauge-invariant geodesic shortcut
(`operational_ftl_geodesic`), and (opt-in) **4D evolving** end-to-end null trace
(`ftl_geo_evolving`; authoritative in search when enabled — see
[4D probe](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17)).

### Stage 4 — Archive update & feedback

`metrics/score/` collapses diagnostics into `ftl_first` fitness. Only `gpu_ok`
candidates compete for a behavior cell. Rejected GRTresna solves log to
`trajectory.jsonl` but never enter the archive.

## The hard consistency rule

`T_ab` used in the GRTresna **constraint solve** must equal `T_ab` used in the
GRTeclyn **evolution**. Otherwise the run starts off-constraint and any apparent
"FTL" is a constraint-relaxation transient, not physics.

The `grtresna_independent_scalars` matter path exists precisely to keep both
sides identical; any new matter sector must be added to **both** sides with
matching analytic forms. Boson star uses **`grtresna_complex_scalar`** on both
GRTresna (`BosonStarBH` + `ComplexScalarField.cpp`) and GRTeclyn
(`ComplexScalarField.hpp` + `RadialRecipeMatterDispatch.hpp`). (Historical
single-field mismatch was fixed via `grtresna_independent_scalars`; see also
[Matter model](#foundational-entries-2026-06-10) and
[Matter selector](#matter-selector-scalar-vs-boson-star).)

## Behavior descriptors (the "diversity" axes)

Current descriptor: **`ftl_lifetime`** (`qd_search/descriptors.py`), added in
[v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17):

- **x** — peak trustworthy `f_geo`, ramped `(f_geo − 1e-3)/(2e-1 − 1e-3)`.
- **y** — **FTL-lifetime fraction** (share of frames with a live shortcut).

Back-compat: **`speed_super`** (v14 default); `speed_horizon` retired after the
`theta_plus` centering bug (see [Status](#map-elites-ftl-discovery-status)).

## Scoring model (the "quality" axis)

Fitness is `ftl_first` in `metrics/score/` (CMA-ES may use `robust_ftl`; see
[v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17)).

- **Gauge-invariant FTL dominates — time-averaged since v15; 4D since v18.**
  `operational_ftl_geodesic` (×1000) or `ftl_geo_evolving` when 4D runs; gated on
  `h_quality_ok` and ray reach; mean over trustworthy magnitudes × `structural_persistence`.
- **Dynamical signal next.** `operational_ftl` (×400) + `ftl_persistence` (×300) outweigh
  coordinate `operational_ftl_solved` (×50, shaping only).
- **Persistence-honest health.** `survival = numerical_survival × structural_persistence`.
- **Vetoes / penalties.** Horizon (−500 when corroborated, [v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17)),
  exotic, instability, stationary warp-lens artifacts.

Geodesic ramp: `(f_geo − 1e-3)/(2e-1 − 1e-3)` — full marks need ~20% shortcut
([v9 recalibration](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17)).

Per-component weight table: `grteclyn-wrapper/src/grteclyn_wrapper/metrics/README.md`.

### Plain-English glossary (compact)

Modes: `weighted` (plain sum) and `ftl_first` (validated FTL dominates).

| FTL signals | Role |
|-------------|------|
| `operational_ftl_geodesic` | Frozen per-frame null-ray shortcut; largest weight when 4D off |
| `ftl_geo_evolving` | 4D end-to-end null trace; authoritative when enabled |
| `operational_ftl`, `ftl_persistence` | Evolved coordinate shortcut + final-frame persistence |
| `operational_ftl_solved`, `ftl_precursor`, `channel_progress`, `shift_drive` | Shaping / t=0 hints; gated |

| Health rewards | Role |
|----------------|------|
| `numerical_survival`, `structural_persistence`, `survival` | Integrator + morphology |
| `stability`, `comoving_stability`, `constraint_health`, `lapse_health` | Geometry / solve quality |
| `energy_condition`, `anec_condition`, `tidal_comfort` | Physical energy rules |
| `curvature_activity`, `nontrivial_geometry`, … | Non-flat geometry rewards |

| Penalties | Role |
|-----------|------|
| `exotic_penalty` | Graded 0..−1.6 for negative-energy matter |
| `horizon_penalty` | −500 veto when lapse-collapsed trapped surface corroborated |
| `instability_penalty`, `qei_penalty`, `boundary_penalty` | Geometry churn / bounds |
| `stationary_artifact_penalty` | Static lens pretending to be FTL |

**Non-triviality gate** — health rewards off for flat vacuum.

## Code map (where everything lives)

| Concern | Path |
|---------|------|
| **MAP-Elites** QD loop, archive, descriptors | `grteclyn-wrapper/src/grteclyn_wrapper/search/qd_search/` |
| **Pre-GPU rejection learning** (near-miss pool, shadow archive, Ham/Mom descriptors) | `search/pre_gpu/` |
| Shadow archive persistence | `qd_search/pre_gpu_archive.py` → `pre_gpu_archive.json` |
| **CMA-ES** optimize loop, warm-start | `grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/` |
| Search-space defs (shared QD + CMA-ES) | `search/optimize/spaces.py` |
| FTL champion retention | `search/ftl_retention.py` |
| Scoring (`ftl_first`, `robust_ftl`, `general_ftl`) | `grteclyn-wrapper/src/grteclyn_wrapper/metrics/score/` |
| Metric aggregation | `grteclyn-wrapper/src/grteclyn_wrapper/metrics/aggregation/collector.py` |
| FTL probes | `grteclyn-wrapper/src/grteclyn_wrapper/metrics/probes/ftl/` |
| **4D evolving geodesic** + metric-stack cache | `metrics/probes/ftl/evolving_geodesic.py`, `metric_field.py`, `metric_stack_cache.py` |
| Plotfile → frames + `ftl_timeseries.dat` | `visualisation/process_wave/consume_plotfiles/` |
| **Incremental HQ scoring** | `metrics/aggregation/incremental.py` |
| HQ promotion launcher | `scripts/campaigns/hq/run_batch.sh`, `campaigns/hq/replay_eval.py` |
| Frame → movie stitching | `scripts/plot/make_movies.sh` |
| Matter (evolution, real scalars) | `Source/Matter/GRTresnaIndependentScalars.{hpp,impl.hpp}`, `Examples/RadialRecipe/` |
| Matter (evolution, boson star) | `Source/Matter/ComplexScalarField.hpp`, `RadialRecipeMatterDispatch.hpp` |
| Matter (initial data, real scalars) | `../GRTresna/Examples/ScalarFieldBH/` |
| Matter (initial data, boson star) | `../GRTresna/Examples/BosonStarBH/`, `Source/Matter/ComplexScalarField.cpp` |
| Matter selector / wiring | `grteclyn_wrapper/grtresna/matter_models.py`, `matter_wiring.py` |
| Boson search space (7-D) | `search/optimize/spaces.py` → `grtresna_boson_star_search_space()` |
| Campaign launchers | `scripts/campaigns/qd/run.sh`, `cmaes/run.sh`, `general_ftl/run_all.sh`, **`boson_star/ftl_shell_run.sh`**, `boson_star/run.sh` |

## Building the binaries (GRTresna + GRTeclyn)

**GRTresna** = Chombo + conda-OpenMPI (CPU/MPI). **GRTeclyn** = AMReX + CUDA (GPU).

### One env to set first (every shell)

```bash
export GRTRESNA_ENV=/home/jovyan/.mlspace/envs/grtresna
export SIM_ROOT=/home/jovyan/nachevsky/test/simulation
export CHOMBO_HOME="${SIM_ROOT}/Chombo/lib"
export CONDA_PREFIX="${GRTRESNA_ENV}"
export PATH="${GRTRESNA_ENV}/bin:${PATH}"
export LD_LIBRARY_PATH="${GRTRESNA_ENV}/lib:${LD_LIBRARY_PATH:-}"
```

Shortcut: `source grteclyn-wrapper/scripts/lib/env.sh`. Full recipe:
`grteclyn-wrapper/src/grteclyn_wrapper/README.md`.

### Build GRTresna (initial-data solver, MPI)

**Real scalars (default):**

```bash
cd "${SIM_ROOT}/GRTresna/Examples/ScalarFieldBH"
PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
  make all -j4 CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE
# -> Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex
```

**Boson star (complex scalar):**

```bash
cd "${SIM_ROOT}/GRTresna/Examples/BosonStarBH"
PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
  make all -j4 CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE
# -> Main_BosonStarBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex
```

First time: `cd "${CHOMBO_HOME}" && make lib -j"$(nproc)"`. Header-only edits need force-relink.

### Build GRTeclyn (evolution, GPU)

```bash
cd "${SIM_ROOT}/GRTeclyn/Examples/RadialRecipe"
make COMP=gnu USE_CUDA=TRUE USE_MPI=FALSE CUDA_ARCH=90 -j"$(nproc)"   # H100
make COMP=gnu USE_CUDA=TRUE USE_MPI=TRUE  CUDA_ARCH=90 -j"$(nproc)"   # MPI+CUDA
```

`CUDA_ARCH`: `90` = H100, `80` = A100, `70` = V100.

### Common failures → fixes

| Symptom | Fix |
|---------|-----|
| `mpicxx` / `gfortran` not found | `export PATH="${GRTRESNA_ENV}/bin:${PATH}"` |
| `cannot find -lhdf5` | `export CONDA_PREFIX="${GRTRESNA_ENV}"` before `make` |
| `libmpi.so` / `libhdf5.so` missing at runtime | `export LD_LIBRARY_PATH="${GRTRESNA_ENV}/lib:..."` |
| `CHOMBO_HOME` undefined | `export CHOMBO_HOME="${SIM_ROOT}/Chombo/lib"` |
| Header edit "did nothing" | force-relink |
| `/bin/csh: No such file` | point Chombo at conda `tcsh` (README step 6) |
| `gcc versions later than 12 are not supported` (nvcc) | use gcc≤12, or `-allow-unsupported-compiler`, or existing CUDA binary |

## How to run a campaign

> **Copy-paste commands:** see [Quick start — running campaigns](#quick-start--running-campaigns) at the top of this doc.

### MAP-Elites (QD)

Same as [Stage 0 — MAP-Elites](#stage-0--map-elites-qd). Results: `runs/grtresna_qd/<name>/`.

### CMA-ES refinement after MAP-Elites

| Line | QD source | Objective | Warm-start | Result | Doc |
|------|-----------|-----------|------------|--------|-----|
| **v22 wormhole** | `general_ftl_wormhole_v21` | **`general_ftl`** | QD **063** → CMA-ES **046** (+14.2 pts) | [v22 CMA-ES](#v22-cma-es-wormhole-refinement-general_ftl_wormhole_cmaes_v1-2026-06-18) |
| **v20** | `general_ftl_{wormhole,ring,spin}` | **`general_ftl`** | Ring eval **43** (TBD) | [v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17) |
| **v18 / ftl_4d** | `ftl_4d_v1` | **`ftl_first`** | QD **156** → CMA-ES **144** | [v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17) |
| **v17** | `ftl_discovery_v16` | `robust_ftl` | eval **739** (not king 233) | [v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17) |

**Rule:** CMA-ES must match QD `OBJECTIVE_MODE`, grid, stop time, and 4D profile
(`campaigns/lib/search_common.sh`). Do **not** switch objectives on a warm-started trajectory.
Full launch recipe → [Stage 1 — CMA-ES](#stage-1--cma-es-optimization).

### HQ promotion (full resolution)

Same as [Stage 2 — HQ promotion](#stage-2--hq-promotion). Search at **N=128, L=64, ml=1–2, t=16** →
HQ at **N=256, L=128, ml=3, t=30**.

| Knob | Search | HQ |
|------|--------|-----|
| Grid | 128³, L=64, ml=2 | **256³, L=128, ml=3** |
| Stop time | 16 | **30** |
| 4D mode | `search` (cheap) | `hq` (full stack) |
| Geodesic dirs | **`x y z`** (`GRTECLYN_GEO_DIRECTIONS`) | **`x y z`** (same; wormhole on **z**) |
| Objective | campaign-specific | **`general_ftl`** for v20/v22 wormhole |

Legacy v17 CMA-ES (`robust_ftl`, `ftl_cmaes_v17_robust`): see [v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17).

**Incremental scoring:** when `--evolving-geodesic` is on, only `ftl_geo_evolving` earns
geodesic FTL credit until the end-of-run 4D trace completes — mid-run totals are not comparable
to search finals. See [v10–v20 HQ eval 144](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17).

```bash
tail -f runs/grtresna_promote/*/small_data/score_timeseries.jsonl
bash scripts/plot/make_movies.sh runs/grtresna_promote/<name> --framerate 10
# → writes runs/grtresna_promote/<name>/movies/movie_<field>_<axis>.mp4
```

Concrete candidate lists and results → [campaign log](#campaign-log--runs-analysis).

---

## Campaign log / runs analysis

Reverse-chronological journal. Quick index:

| Campaign / section | Date | Headline |
|--------------------|------|----------|
| [**Static shell post-mortem (`scalar_shell_ftl_rl_v1`)**](#static-shell-post-mortem-scalar_shell_ftl_rl_v1-2026-06-24) | **06-24** | `shell_static=1` zeroed all velocity/omega/shift dims (~5 of 22 dead no-ops). Warp motor structurally off (Pi=0 → S_i=0 → beta=0). Champions are static exotic lenses (maxed `exotic_penalty`, weak 4D). Archive collapsed to **3/64 cells**, stalled 17 iterations. **54% of compute wasted** on non-converging constraint solves. Next: un-pin `shell_static`, activate velocity+spin DOF for moving/rotating warps. |
| [**Paired shell FTL comparison (`boson_shell_ftl_rl_v1` vs `scalar_shell_ftl_rl_v1`)**](#paired-shell-ftl-comparison-boson-vs-scalar-2026-06-23) | **06-23 complete** | **200+200 evals**, matched campaign knobs. Scalar: **92 `gpu_ok`**, **32/92** with `ftl_geo_evolving>0`, champion **eval 166** (**869**, persist 0.76). Boson: **94 `gpu_ok`**, **0 FTL**. **No boson rerun**; promote scalar **eval 166/126** for RL Gate 2. |
| [**Boson shell FTL for RL (`boson_shell_ftl_rl_v1`)**](#boson-shell-ftl-for-rl-boson_shell_ftl_rl_v1-2026-06-23) | **06-23 complete** | Boson arm of paired comparison — **0/94** `f_geo>0`. See [paired results](#paired-shell-ftl-comparison-boson-vs-scalar-2026-06-23). |
| [**Boson star: unpinned QD + frames (v3)**](#boson-star-unpinned-qd--frames-v3-2026-06-18) | **06-18 complete** | First **unpinned 7-D** boson MAP-Elites: **12/12 evals**, **8 `gpu_ok`**, champion **eval 004** (−33.2). Frames on H100; `scalar_activity`/`phi` projections visible. |
| [**Boson star: integration + GPU smoke**](#boson-star-matter-sector-integration--gpu-smoke-2026-06-18) | **06-18 complete** | Matter selector wired; pinned GPU smoke **2/2 `gpu_ok`**. Phase 1 = single centered Gaussian. |
| [**v22 CMA-ES: wormhole refinement**](#v22-cma-es-wormhole-refinement-general_ftl_wormhole_cmaes_v1-2026-06-18) | **06-18 complete** | QD **063** (165.6) → CMA-ES **046** (**179.8**, +14.2). [HQ eval 046](#hq-eval-046-final-results-t30) **complete** — peak **7.57%** 4D @ t≈15.6, **−546** final @ t=30 (horizon kill). [Movies](#hq-eval-046-final-results-t30) in `movies/`. |
| [**v22: Pre-GPU learning + v21 resume**](#v22-pre-gpu-rejection-learning--v21-resume-2026-06-18) | **06-18 complete** | **200 evals** on `general_ftl_wormhole_v21`. Pipelined QD + pre-GPU learning. **Champion eval 191** (161.9, survival 1.00, 18.5% 4D). Eval 063 holds score/FTL records (165.6, 19.3%). [Top 3 + FTL champions](#v22-final-results-top-3--ftl-champions). |
| [**v21: Pipelined QD + GPU tenancy**](#v21-pipelined-qd--gpu-tenancy-tuning-2026-06-17) | **06-17 stopped** | **Pipelined MAP-Elites** (`GpuPool` + `EvalPipeline`). **5 slots/GPU overloaded** H100s at t=16 (~3× slower/evol). **Working config:** 8 GPUs × **1 slot/GPU** + continuous GRTresna; **bottleneck** `MAX_CONCURRENT_GRTRESNA=3` → raise to **5+**. **26 evals:** 11 `gpu_ok`. Cold-start gap → [fixed in v22](#v22-pre-gpu-rejection-learning--v21-resume-2026-06-18). |
| [**v10–v20: pipeline + runs**](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17) | **06-11 → 06-17** | Scoring/geodesic hardening (v7–v16) → **4D probe + HQ** (eval **144** verified **~8%**) → **general_ftl** QD (ring eval **43** **~3.9%**). See glued section. |
| [**Foundational (06-10)**](#foundational-entries-2026-06-10) | 06-10 | Matter model, navigation overhaul, status reset |

---

## Boson shell FTL for RL (`boson_shell_ftl_rl_v1`, 2026-06-23)

**Status:** **complete** — **200/200 evals** (finished **2026-06-23 ~11:27 UTC**). MAP-Elites QD for
closed-loop RL initial-data elites.

**Motivation.** Historical wormhole QD (`general_ftl_wormhole_v*`) used **real scalar** lumps.
RL pump / lump tracker requires **`grtresna_complex_scalar`** (boson shell, same wiring as
[grlab splash](../grlab/README.md)). Splash QD optimizes **`critical_collapse`**, not FTL.
This campaign is the first **boson shell + `general_ftl`** survey aimed at elites with finite
`ftl_geo_evolving` for [RL Tax Man training](../RL/LabJournal.md).

**Run dir:** `runs/grtresna_qd/boson_shell_ftl_rl_v1/`

### Final results (200 evals)

| Metric | Value |
|--------|------:|
| **`gpu_ok`** | **94** (47%) |
| **`postload_rejected`** | **76** (38%; **45%** of gate-tested 170) |
| **`grtresna_rejected`** | 19 |
| **`grtresna_failed`** | 4 (incl. attempt-1 `h5py`; 3 late-run) |
| **`pipeline_interrupted`** | 7 (prior kill/resume) |
| **`f_geo_peak > 0`** | **0 / 94** — no evolved FTL shortcut |
| **`structural_persistence`** (gpu_ok) | mean **0.39**, median **0.34**; **23/94** ≥ 0.5; **9/94** = **1.0** |
| **Best score** | **eval 100** → **21.6** (persist **1.0**, still **f_geo = 0**) |
| **Best persistence** | evals **39, 66, 87, 100, 112, …** at **1.0** (no FTL) |

**Verdict.** Boson shell + exotic wedge can stay **cohesive** (9 configs with perfect persistence
@ t=16) but **did not produce any FTL geometry** in 200 evals. Matter stability alone is
insufficient for RL chassis — need **`f_geo > 0`** plus persistence. See
[paired comparison](#paired-shell-ftl-comparison-boson-vs-scalar-2026-06-23) (scalar arm).

**Champion for archive inspection:** `eval_000100` (score 21.6, persist 1.0, frames in
`eval_000100/frames/`).

### Attempt 1 — `grtresna_failed` (missing `h5py`)

First launch (**~07:36 UTC**) used exotic wedge search and BosonStarBH solves correctly, but
**every scored eval** landed as **`grtresna_failed`** at score **−350** with:

```text
reason: ModuleNotFoundError("No module named 'h5py'")
```

**Root cause.** GRTresna (Fortran/Chombo) writes `InitialDataFinal.3d.hdf5`; the wrapper converts
that to `initial_data.gridinit` in `grtresna/io.py` via **`h5py`**. That package was documented
as a manual `uv pip install h5py` step but was **not** listed in `grteclyn-wrapper/pyproject.toml`.
Earlier wormhole/splash campaigns worked because `h5py` had been installed by hand and survived
in the venv; a later **`uv sync`** on the managed venv (Python 3.14) dropped the undeclared
package. Launch also used **`SKIP_QD_PREFLIGHT_TESTS=1`**, so nothing imported `h5py` before
GPU/GRTresna slots were allocated.

**Symptom vs physics.** Not a boson-matter or exotic-wedge bug — BosonStarBH completed
(`Ham_and_Mom_errors.txt`, HDF5 under `eval_*/grtresna/Outputs/`), then Python died at the
Chombo→`.gridinit` bridge before GPU evolution.

### Fix (2026-06-23)

| Change | File / action |
|--------|----------------|
| Declare dependency | `grteclyn-wrapper/pyproject.toml` — **`h5py>=3.10`** in `[project].dependencies` |
| Install in venv | `cd grteclyn-wrapper && uv sync` |
| Preflight guard | `tests/grtresna/test_grtresna_h5py_dep.py` added to `search_common.sh` preflight list |
| Clean restart | Kill campaign, `rm -rf runs/grtresna_qd/boson_shell_ftl_rl_v1`, relaunch |

**Verification (attempt 2).** After restart, eval dirs show **`initial_data.gridinit`** (HDF5
conversion succeeds). Expect normal mix of `gpu_ok`, `grtresna_rejected`, `postload_rejected` —
not blanket `grtresna_failed` / `h5py`.

**Ops note.** After any `uv sync`, confirm `uv run python -c "import h5py"`. Wormhole-era
campaigns never hit this because `h5py` was an undeclared manual install; it must stay in
`pyproject.toml` for all GRTresna QD paths (scalar and boson).

### Launch (attempt 2)

```bash
cd grteclyn-wrapper
uv sync   # ensures h5py>=3.10 from pyproject.toml

QD_NAME=boson_shell_ftl_rl_v1 \
QD_TARGET_EVALS=200 \
QD_ITERATIONS=30 \
STOP_TIME=16.0 \
PLOT_INTERVAL=320 \
GRTECLYN_FRAMES=1 \
GPU_IDS="0 1 2 3 4 5 6 7" \
GPU_SLOTS_PER_DEVICE=1 \
MAX_CONCURRENT_GRTRESNA=5 \
BATCH_SIZE=8 \
  nohup bash scripts/campaigns/boson_star/ftl_shell_run.sh \
  > ../runs/boson_shell_ftl_rl_v1.launch.log 2>&1 &
```

Use **`SKIP_QD_PREFLIGHT_TESTS=1`** only if pytest is not installed in the venv; the h5py
preflight test is skipped in that case — run `uv run python -c "import h5py"` manually first.

| Knob | Value |
|------|-------|
| Launcher | `scripts/campaigns/boson_star/ftl_shell_run.sh` |
| Matter | `GRTRESNA_MATTER_SECTOR=boson_star`, ansatz **shell**, `grtresna_complex_scalar` |
| Search space | `grtresna_boson_shell_search_space(static=True, allow_exotic=True)` — shell geometry + mass/λ/ω + **exotic wedge** |
| Objective | **`general_ftl`**, descriptor **`ftl_lifetime`**, 4D geodesic **`search`** |
| Grid | N=128, L=64, ml=1 (search defaults) |
| Stop / plots | **t=16**, **`PLOT_INTERVAL=320`** → **6** plotfiles per eval |
| Frames | **`GRTECLYN_FRAMES=1`** |
| Pipeline | v22-style (`USE_PIPELINE=1`, `MAX_CONCURRENT_GRTRESNA=5`, pre-GPU learning) |
| Pins | `grtresna_shell_static=1` only (**no** `scalar_sign` pin) |

### Exotic matter — **enabled (shell wedge search)**

`ftl_shell_run.sh` sets **`GRTRESNA_BOSON_ALLOW_EXOTIC=1`** by default:

- Search dims **`grtresna_shell_exotic_fraction`** and **`grtresna_shell_exotic_phase`** (same idiom as real-scalar wormhole QD)
- Lump exotic flags wired into BosonStarBH paint; exotic-safe GRTresna solver when any lump is exotic
- **`grtresna_scalar_sign` not pinned** — canonical coupling + searchable exotic wedge (not global phantom-only)

Disable exotic wedge search (splash-style canonical-only):

```bash
GRTRESNA_BOSON_ALLOW_EXOTIC=0 PIN_DIMS="grtresna_scalar_sign=1 grtresna_shell_static=1" \
  bash scripts/campaigns/boson_star/ftl_shell_run.sh
```

Global phantom boson (whole-field `scalar_sign=-1`, still no exotic-fraction dims):

```bash
GRTRESNA_BOSON_ALLOW_EXOTIC=0 GRTRESNA_MATTER_COUPLING=exotic \
  PIN_DIMS="grtresna_shell_static=1" bash scripts/campaigns/boson_star/ftl_shell_run.sh
```

### Monitor

```bash
tail -f runs/boson_shell_ftl_rl_v1.launch.log
tail -f runs/grtresna_qd/boson_shell_ftl_rl_v1/trajectory.jsonl
cat runs/grtresna_qd/boson_shell_ftl_rl_v1/ftl_champions.json
```

### Paired shell FTL comparison (boson vs scalar, 2026-06-23)

**Status:** **both complete** — **200/200 evals** each. Same campaign **intent** (static 5-lump
shell, exotic wedge, **`general_ftl`**, t=16, frames, postload **3e-2**, 8×GPU pipeline).
**Only matter sector differs** (`grtresna_complex_scalar` vs `grtresna_independent_scalars`).

**Run dirs:** `runs/grtresna_qd/boson_shell_ftl_rl_v1/` ·
`runs/grtresna_qd/scalar_shell_ftl_rl_v1/`

#### Were these the same conceptual runs?

**Mostly yes on pipeline knobs; not bit-identical on search space.**

| Matched | Boson-only / scalar-only |
|---------|--------------------------|
| Objective **`general_ftl`**, descriptor **`ftl_lifetime`**, 4D geodesic **`search`** | Matter model + GRTresna solver (`BosonStarBH` vs `ScalarFieldBH`) |
| Grid **L=64, N=128, dx=0.5**, **t=16**, **PLOT_INTERVAL=320**, **frames on** | Search bounds: boson amp **0.04–0.12**, mass **0.05–0.35**, **`grtresna_bs_omega`**; scalar amp **0.08–0.16**, mass **0.3–1.5**, **`grtresna_scalar_lambda`**, velocity dims in space (zeroed when **`shell_static=1` pinned**) |
| Pin **`grtresna_shell_static=1`**, exotic wedge ON, **5 lumps** | Boson excludes velocity dims from optimizer; scalar shell space is wider (~23-D vs ~18-D) |
| **`POSTLOAD_MAX_HAM_L2=3e-2`**, pipeline, **`SKIP_QD_PREFLIGHT_TESTS=1`** (scalar launch) | Boson attempt 1 **`h5py`** failure + resume artifacts (`pipeline_interrupted` ×7) |

Fair conclusion: **same experiment design, different matter physics and slightly different
optimizer bounds** — sufficient to test “can boson shell + exotic wedge produce FTL like real
scalar shell under matched scoring?”

#### Side-by-side results (200 evals each)

| Metric | **Boson** | **Scalar** |
|--------|----------:|-----------:|
| **`gpu_ok`** | 94 (47%) | 92 (46%) |
| **`postload_rejected`** | 76 (**45%** of gate-tested) | 26 (**22%**) |
| **`grtresna_rejected`** | 19 (10%) | 74 (**37%**) |
| **`f_geo_peak > 0`** | **0 / 94** | **33 / 92** |
| **`ftl_geo_evolving > 0`** | **0 / 94** | **32 / 92** |
| **`operational_ftl > 0`** | **0 / 94** | **5 / 92** |
| **Persistence mean / median** | 0.39 / 0.34 | **0.74 / 0.87** |
| **Persist = 1.0** | 9 / 94 | **37 / 92** |
| **Score mean / median** | −51 / −54 | −73 / −29 |
| **Best score** | **21.6** (eval **100**, persist 1.0, **no FTL**) | **869.3** (eval **166**, persist 0.76, **ftl_p=0.96**) |
| **Best `f_geo` peak** | 0 | **0.86** (eval **41**, spike only → score **−12**) |

#### Why scores differ so much (not a scorer bug)

Both use the same **`general_ftl`** total (`objectives.py`): FTL terms dominate at
**+1000×`ftl_geo_evolving` + 600×`ftl_persistence` + 200×`operational_ftl`**. Health-only
ceiling (no FTL) is **~O(20–30)** after **`exotic_penalty`** (−112 at −1.6).

| Boson eval 100 (score **21.6**) | Scalar eval 166 (score **869**) |
|--------------------------------|----------------------------------|
| All FTL components **0** | **`ftl_persistence` 0.96** → **+577** |
| Health + persist 1.0 | **`operational_ftl` 0.96** → **+192** |
| | **`ftl_geo_evolving` 0.10** → **+105** |

**Root cause:** boson static exotic shell **never opens a geodesic shortcut** (0/94); real scalar
exotic shell **does** (32/92 sustained). Higher scalar persistence is a **secondary** effect
(heavier searched mass keeps lumps bound); it does **not** explain the score gap — **FTL terms do**.

One scalar spike (eval **41**, peak **86%** `f_geo`) scored **−12** because
**`ftl_lifetime_fraction=0`** — scorer correctly rejects transient gauge artifacts.

#### Verdict and rerun?

| Question | Answer |
|----------|--------|
| **Same conceptual experiment?** | **Yes** on pipeline/objective/grid; **minor** search-bound asymmetry (documented above). |
| **Why boson scores ~20 vs scalar ~870?** | **Zero FTL signal** vs **strong `ftl_persistence` + geodesic credit** — by design of `general_ftl`. |
| **Rerun boson?** | **No** — 200 evals conclusively show **0 FTL**; rerun without new physics/search changes would not change the conclusion. |
| **Rerun scalar?** | **No** for the comparison — scalar arm succeeded. Optional next step: **CMA-ES refine** from **eval 166** or **126** (both persist ≥0.76, finite FTL) for RL chassis, not another blind 200-eval QD. |
| **RL chassis?** | **Scalar shell elite** (e.g. **`eval_000166`**, **`eval_000126`**). Boson arm **rejected** for FTL RL. Run **Gate 2** pump test on scalar elite before PPO. |

**Scalar champions (inspect frames / `initial_data.gridinit`):**

| Eval | Score | Persist | `ftl_persistence` | `ftl_geo_evolving` | Notes |
|------|------:|--------:|------------------:|-------------------:|-------|
| **166** | **869** | 0.76 | **0.96** | 0.10 | Overall score leader |
| **126** | 840 | **1.0** | 0.77 | **0.17** | Best persist + FTL combo |
| **64** | 592 | 0.66 | 0.65 | 0.12 | Mid-run reference |
| **41** | −12 | 0.95 | 0 | 0 | **Anti-pattern** — spike only |

Launchers: `boson_star/ftl_shell_run.sh` · `general_ftl/scalar_shell_ftl_run.sh`

```bash
# Scalar arm (completed 2026-06-23)
SKIP_QD_PREFLIGHT_TESTS=1 QD_NAME=scalar_shell_ftl_rl_v1 QD_TARGET_EVALS=200 \
  GPU_IDS="0 1 2 3 4 5 6 7" GPU_SLOTS_PER_DEVICE=1 \
  bash scripts/campaigns/general_ftl/scalar_shell_ftl_run.sh
```

### Handoff to RL

**Boson arm:** **rejected** — **0/94** with **`f_geo > 0`**. Persistence-only elites (e.g. eval 100)
max out at score **~22**; not usable for **`general_ftl`** Tax Man training.

**Scalar arm:** **promote** **`eval_000166`** or **`eval_000126`** (`gpu_ok`, finite
**`ftl_geo_evolving`**, `initial_data.gridinit` present). Next: **Gate 2** actuation on scalar
elite with **`grtresna_independent_scalars`** pump (wormhole eval 063 is alternate IVP reference,
not this paired shell search).

Interim pump proof remains splash boson **`spacetime_splash_v14_moving/eval_000010`** until Gate 2
passes on the scalar FTL chassis.

### Static shell post-mortem (`scalar_shell_ftl_rl_v1`, 2026-06-24)

Deep analysis of the `scalar_shell_ftl_rl_v1` campaign reveals that **the warp/frame-drag
motor was structurally disabled** throughout the run. The FTL found is a static exotic-matter
lens — the weakest, most artifact-prone mechanism class — and the search stalled hard with
most compute budget wasted.

#### Finding 1 — `shell_static=1` zeroed all dynamical DOF

The launcher pins `grtresna_shell_static=1`
([`scalar_shell_ftl_run.sh`](../../grteclyn-wrapper/scripts/campaigns/general_ftl/scalar_shell_ftl_run.sh)):

```bash
export PIN_DIMS="${PIN_DIMS:-grtresna_shell_static=1}"
```

The shell builder hard-zeroes velocity and rotation when this pin is active
([`config.py`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/config.py)):

```python
if int(round(get_float("grtresna_shell_static", 0.0))) >= 1:
    v_tor = v_pol = v_rad = omega = 0.0
```

The QD optimizer **sampled** these dims (e.g. eval 166: `toroidal_velocity=-0.63`, `omega=0.02`,
`shift_seed=0.45`), but all values were **overwritten to zero** before the GRTresna solve.
Additionally, `shift_seed` is inert when matter is static because
[`paint_shift_seed_on_grid`](../../grteclyn-wrapper/src/grteclyn_wrapper/grtresna/lump_fields.py)
aligns the seed to matter momentum density `S_i`, which is identically zero when `Pi=0`.

**Dead search dimensions — scalar arm (~5 of 22, ~23% of search space):**

| Dimension | Sampled range | Actual effect |
|-----------|--------------|---------------|
| `grtresna_shell_toroidal_velocity` | \[-1.2, 1.2\] | **Zeroed** by `shell_static=1` |
| `grtresna_shell_poloidal_velocity` | \[-0.8, 0.8\] | **Zeroed** |
| `grtresna_shell_radial_velocity` | \[-0.3, 0.3\] | **Zeroed** |
| `grtresna_shell_omega` | \[-0.5, 0.5\] | **Zeroed** |
| `grtresna_shift_seed` | \[-0.6, 0.6\] | **Inert** (no momentum to align to) |

**Boson arm (`boson_shell_ftl_rl_v1`) — same static pin, even fewer DOF:**

The boson launcher ([`ftl_shell_run.sh`](../../grteclyn-wrapper/scripts/campaigns/boson_star/ftl_shell_run.sh))
also pins `grtresna_shell_static=1`, and additionally the boson shell search-space builder
(`grtresna_boson_shell_search_space(static=True)`) **excluded velocity dims entirely** from
the optimizer
([`spaces.py`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/spaces.py)):

```python
if static:
    skip |= {
        "grtresna_shell_toroidal_velocity",
        "grtresna_shell_poloidal_velocity",
        "grtresna_shell_radial_velocity",
        "grtresna_shell_omega",
    }
```

So unlike the scalar arm (which wasted 4 dims sampling velocities that were zeroed),
the boson arm **never even sampled** them — the 19-D boson search space had zero kinematic
DOF. Only `grtresna_shift_seed` remained (inert, as above). The boson campaign's
`base_overrides` confirm: `grtresna_shell_static: 1.0`,
`grtresna_matter_model: grtresna_complex_scalar`.

**Both arms of the paired comparison had the warp motor structurally off.** Neither
campaign tested a moving or rotating matter configuration. The conclusion "boson produces
0 FTL" vs "scalar produces FTL" applies only to the **static exotic-lens** mechanism class —
it says nothing about dynamical/rotating warps, which were never searched. The boson arm
might perform differently with momentum-carrying matter; this remains untested.

#### Finding 2 — warp/frame-drag motor structurally off

The physics chain that creates a moving warp is:

```
velocity/omega → Pi (scalar momentum) → S_i = sign·Pi·∂_iφ → shift β^i (warp channel)
```

With `v_tor = v_pol = v_rad = omega = 0`: **`Pi = 0` → `S_i = 0` → `β^i = 0`**. The initial data
is time-symmetric (zero momentum, zero shift). GRTresna's momentum constraint `M_i = 0` is trivially
satisfied. The matter has no kinetic energy, no frame-dragging, no shift motor — it's a pure
static exotic-matter curvature lens.

#### Finding 3 — champions are static exotic lenses (weak, artifact-prone)

| Eval | Score | `exotic_penalty` | `instability_penalty` | `ftl_geo_evolving` (4D) | Notes |
|------|------:|:----------------:|:---------------------:|:-----------------------:|-------|
| **166** | 869 | **-1.600** (maxed) | **-0.960** | 0.10 | Static exotic lens; score from coordinate proxy |
| **126** | 840 | **-1.600** | **-0.951** | 0.17 | Best 4D but still weak |
| **64** | 592 | **-1.600** | **-0.956** | 0.12 | Same mechanism class |
| **51** | 311 | **-1.600** | **-0.955** | — | `horizon_penalty -1.0` (formed horizon) |
| **38** | 78 | **-1.600** | **-0.957** | — | `stationary_artifact_penalty -1.0` |

All champions share: **maxed exotic penalty** (requires full WEC-violating matter), **near-maxed
instability** (riding the edge of numerical blowup), and **weak authoritative 4D geodesic signal**
(0.10–0.17). Score is carried by coordinate proxies (`operational_ftl` ~0.96), not the
gauge-invariant 4D trace. Some hit the `stationary_artifact_penalty` — the scorer flagging
exactly this class of static lens.

#### Finding 4 — archive collapsed, search stalled, massive compute waste

| Metric | Value | Concern |
|--------|-------|---------|
| **Archive coverage** | **3/64 cells** (0.047) | QD diversity completely collapsed — all top elites in **cell [1,7]** |
| **Stall duration** | **17 iterations** (iter 8→25) at constant coverage, minimal score gain | Diminishing returns after ~80 evals |
| **`grtresna_rejected`** | **74/200** (37%) — 53 at `Ham=100%`, 14 NaN | Proposer sampling far outside feasible region |
| **`postload_rejected`** | **26/200** (13%) | Constraints blow up at higher res after load |
| **Total non-GPU** | **108/200** (**54%**) | Over half the compute budget never reached evolution |
| **Descriptor x** (`f_geo` ramp) | mean **0.028**, max 0.195 | Genuine FTL near-zero for most configs |
| **Descriptor y** (`ftl_lifetime`) | mean **0.348**; 32/92 have y > 0 | Only one basin has any FTL lifetime |

**Convergence curve:** best score −8.7 → 78 → 311 → 592 (iter 8) → 839 (iter 16) → **869** (iter 21) → flat. Zero new cells after iter 8.

#### Finding 5 — RL evolution-control was premised on the wrong foundation

The RL research plan ([`research/RL/research.md`](../RL/research.md)) proposes a matter pump
on `c_Pi`/`c_Pi2` (U(1) rotation) to sustain an FTL throat past t≈21. But:

1. **Actuator–chassis mismatch:** the pump requires `grtresna_complex_scalar` (boson), but the
   only FTL-producing chassis is the real-scalar shell (this campaign).
2. **The FTL here is a static exotic lens, not a dynamical warp.** There is no "throat" to sustain —
   the geometry warps from negative-energy matter at rest. The failure mode is not "throat collapse
   from dispersion" but rather numerical instability of an inherently fragile exotic configuration.
3. **The actuator cannot supply what the geometry needs.** Persistent static-lens FTL requires
   sustained WEC violation, which the scorer penalizes (`exotic_penalty`, `qei_penalty`,
   `horizon_penalty`). A pump adding canonical energy doesn't help; a pump adding exotic energy
   gets clawed back by the scorer.

#### Conclusion — next campaign: dynamical (moving + rotating) shell

**The single highest-leverage change is to un-pin `shell_static`** (`=0`), activating the 5
dormant dimensions and enabling the frame-drag/warp motor (`Pi ≠ 0 → S_i ≠ 0 → β^i ≠ 0`).
This shifts the search from static exotic lenses toward **moving/rotating warps** — the mechanism
class that the [`RotatingWormholeCollapse`](../../Examples/RotatingWormholeCollapse/README.md) example
demonstrates with co-evolving exotic scalar and m=2 quadrupole rotation.

**Concrete changes for a dynamical shell campaign:**

| Change | Detail | Risk |
|--------|--------|------|
| **Un-pin `shell_static`** | `PIN_DIMS=""` or `PIN_DIMS="grtresna_shell_static=0"` | Higher GRTresna rejection (momentum constraint harder); mitigate with `GRTRESNA_MAX_MOM_PCT≥10` |
| **Tighten velocity bounds** (initially) | `toroidal ±0.4`, `radial ±0.15`, `omega ±0.2` | Keeps Mom solve convergent; widen as baseline stabilizes |
| **Add higher angular momentum modes** | Extend `grtresna_shell_mode` upper from 2.0 to 4.0 (m=3, m=4 quadrupole/octupole) | m≥3 modes may need finer grid (ml≥2) |
| **Differential rotation** (future) | Per-lump omega instead of single global `shell_omega` | Requires search-space extension in `spaces.py` |
| **Pre-GPU feasibility** | Strengthen `pre_gpu_learning` or add constraint-residual surrogate to reduce the 54% rejection rate | Training data available from `trajectory.jsonl` + `pre_gpu_archive.json` |
| **Increase GRTresna iterations** | `ITERATIONS=30` → `50` for momentum-carrying configs | Slower per-eval; parallelize with `MAX_CONCURRENT_GRTRESNA≥5` |

**Proposed launcher override:**

```bash
cd grteclyn-wrapper
QD_NAME=scalar_shell_dynamical_v1 \
QD_TARGET_EVALS=200 \
PIN_DIMS="" \
STOP_TIME=16.0 \
PLOT_INTERVAL=320 \
GRTECLYN_FRAMES=1 \
ITERATIONS=50 \
GRTRESNA_MAX_MOM_PCT=10 \
GPU_IDS="0 1 2 3 4 5 6 7" \
GPU_SLOTS_PER_DEVICE=1 \
MAX_CONCURRENT_GRTRESNA=5 \
  nohup bash scripts/campaigns/general_ftl/scalar_shell_ftl_run.sh \
  > ../runs/scalar_shell_dynamical_v1.launch.log 2>&1 &
```

---

## Boson star: unpinned QD + frames (v3, 2026-06-18)

**Status:** **complete** — `boson_star_unpinned_frames_v3` reached **12/12 evals**.
First production-style boson-star MAP-Elites run: **unpinned 7-D search**, **2× H100**,
**frames on**, relaxed Ham/Mom gates, solved-FTL gate **off** for boson sector.

**Run dir:** `runs/grtresna_qd/boson_star_unpinned_frames_v3/`

### Launch

```bash
cd grteclyn-wrapper
QD_NAME=boson_star_unpinned_frames_v3 \
  QD_TARGET_EVALS=12 \
  GPU_IDS="0 1" BATCH_SIZE=2 \
  RANKS=4 ITERATIONS=30 \
  GRTRESNA_MAX_HAM_PCT=10 GRTRESNA_MAX_MOM_PCT=10 \
  STOP_TIME=8 PLOT_INTERVAL=64 \
  GRTECLYN_FRAMES=1 \
  bash scripts/campaigns/boson_star/run.sh
```

| Knob | Value |
|------|-------|
| Search space | **7-D** boson (unpinned) |
| Grid | GRTresna 128³ → evolution 128³, L=64, ml=1 |
| GPUs | 2 × H100, 1 slot/GPU |
| Gates | Ham/Mom ≤ **10%**; **`--no-grtresna-solved-ftl-gate`** (boson has no t=0 FTL precursor) |
| Frames | `GRTECLYN_FRAMES=1`; projections **`phi`** + **`scalar_activity`** |

### Results

| Metric | Value |
|--------|-------|
| Evals | **12** |
| `gpu_ok` | **8** |
| `grtresna_rejected` | 2 |
| `postload_rejected` | 2 |
| Archive | **1 elite**, coverage **1.56%** (cell `[0,0]`) |
| **Champion** | **eval 004**, score **−33.2**, tier `constructed`, survival **1.0** |

**Champion params (eval 004):** `mass=0.276`, `λ=0.032`, `φ_c=0.122`, `width=9.02`, `shift_seed=−0.29`.

All `gpu_ok` evals: **`f_geo=0%`**, **`operational_ftl=0`** — expected for centered canonical
boson stars (no warp geometry). Negative scores dominated by `exotic_penalty` (−1.6) and
`stationary_artifact_penalty` (−1.0), not FTL failure.

### Frames

PNG frames under `eval_*/frames/` (~224 PNGs per `gpu_ok` eval). Use:

| Field | Use |
|-------|-----|
| **`phi_z`**, **`phi_proj_*`** | Real scalar profile — visible centered blob |
| **`scalar_activity_*`** | Combined \|Φ\| + \|Π\| (fixed in consume_plotfiles) |
| **`phi_lump0_*`** | Imaginary part (φ₂) — **~0 for canonical** bosons; flat is correct |

**Frame fixes applied during v3 prep:**

1. Plot vars: GRTeclyn names imaginary components **`phi_lump0`/`Pi_lump0`** (not `phi2`/`Pi2` in plotfiles).
2. `scalar_activity` derived field: combine Re+Im when only `phi_lump0` pair present (was summing lumps only → blank projections).
3. Projection fields: **`phi`** instead of **`phi_lump0`**.

Movies: `bash scripts/plot/make_movies.sh runs/grtresna_qd/boson_star_unpinned_frames_v3/eval_000004 --framerate 10`

### Lessons

- **Solved-FTL gate must stay off** for boson campaigns (`search_common.sh` disables when `GRTRESNA_MATTER_SECTOR=boson_star`).
- Random 7-D samples need **`ITERATIONS=30`** and **`GRTRESNA_MAX_MOM_PCT≥10`** at 128³ for acceptable GPU reach (~67%).
- Boson MAP-Elites with `ftl_lifetime` descriptors collapses to cell **`(0,0)`** until geometry ansatz + FTL-capable configs are searched — Phase 1 validates **pipeline**, not FTL discovery.

---

## Boson star: matter-sector integration + GPU smoke (2026-06-18)

**Status:** **Phase 1 infrastructure complete** — superseded for production validation by
[v3 unpinned QD](#boson-star-unpinned-qd--frames-v3-2026-06-18). Real-scalar lump path
**unchanged** (default `GRTRESNA_MATTER_SECTOR=scalar`).

**Purpose.** Implement [action plan step 3](#iv-strategic-action-plan): add **complex scalar /
U(1) boson star** as an opt-in matter sector so MAP-Elites can search non-dispersing matter
without replacing the existing 23-D shell pipeline.

### What shipped

| Layer | Change |
|-------|--------|
| **Matter selector** | `GRTRESNA_MATTER_SECTOR` × `GRTRESNA_MATTER_COUPLING`; CLI + campaign env vars |
| **GRTresna** | Auto-selects `BosonStarBH` example; `grtresna_complex_scalar` in params + `.matter.json` |
| **GRTeclyn** | `recipe_matter_model=grtresna_complex_scalar`; plot vars `phi_lump0`/`Pi_lump0` (imaginary part) |
| **Search space** | `grtresna_boson_star_search_space()` — **7-D** (mass, λ, φ_c, width, ω, sign, shift_seed) |
| **Campaigns** | `scripts/campaigns/boson_star/{run,cmaes_run,hq_run}.sh` |
| **Tests** | `test_matter_selection.py`, `test_boson_star_ansatz.py`, `test_complex_scalar_wiring.py`, … |

### Validation (2026-06-18)

| Check | Result |
|-------|--------|
| Unit tests (`uv run pytest`) | **19** matter/boson tests pass; **79** in campaign preflight gate |
| GRTresna `BosonStarBH` rebuild | `Main_BosonStarBH3d.*.MPI.ex` rebuilt |
| GRTeclyn CPU smoke | `main3d.gnu.ex` — 1 step, gridinit load OK |
| GRTeclyn GPU smoke | `main3d.gnu.CUDA.ex` on H100 — evolution OK |
| **GPU QD campaign** | `boson_star_gpu_smoke2_*` — **2/2 `gpu_ok`**, 1 MAP-Elite, scores **58.4** / **65.3** |

**Smoke campaign** (pinned default boson params, reduced grid for speed):

```bash
cd grteclyn-wrapper
QD_NAME=boson_star_gpu_smoke \
  QD_TARGET_EVALS=2 GPU_IDS="0" BATCH_SIZE=1 \
  RANKS=4 ITERATIONS=30 \
  GRTRESNA_GRIDINIT_NX=64 GRTRESNA_GRIDINIT_NY=64 GRTRESNA_GRIDINIT_NZ=64 \
  GRTRESNA_EVOLUTION_N_FULL=64 GRTRESNA_EVOLUTION_L_FULL=32 \
  STOP_TIME=0.5 GRTRESNA_MAX_MOM_PCT=15.0 \
  PIN_DIMS="grtresna_scalar_mass=0.1 grtresna_scalar_lambda=0.0 grtresna_bs_phi_c=0.08 grtresna_bs_profile_width=8.0 grtresna_bs_omega=0.0 grtresna_scalar_sign=1.0 grtresna_shift_seed=0.0" \
  SKIP_QD_PREFLIGHT_TESTS=1 \
  bash scripts/campaigns/boson_star/run.sh
```

Run dir: `runs/grtresna_qd/boson_star_gpu_smoke2_20260618T175911Z/`.

Per eval: GRTresna solve → postload gate → **`main3d.gnu.CUDA.ex`** → plotfile consume →
`ftl_first` score. `f_geo=0%` expected — centered boson with no warp geometry.

**Earlier failed smoke** (`boson_star_gpu_smoke_*175738Z`): random 7-D samples + 15 NL iters
→ both **`grtresna_rejected`** (Mom > 5% gate) before GPU. Random boson QD needs more
iterations or looser convergence gates at reduced grid.

### Remaining (post-v3)

| Item | Status |
|------|--------|
| Larger boson QD (80+ evals) | Not run |
| `boson_star + exotic` GPU path | Wired, not smoke-tested on GPU |
| CMA-ES → HQ handoff for boson sector | Launchers exist, not validated |
| Multi-lump boson at shell positions | Future — Phase 1 is single centered Gaussian |
| Boson + shell geometry combined search | Future — matter sector and ansatz are orthogonal today |

---

## v22: Pre-GPU rejection learning + v21 resume (2026-06-18)

**Status:** **complete** — `general_ftl_wormhole_v21` reached **200 evals** (2026-06-18).
Resume from eval 98 after [convergence-fix](#v22-launch-resume) that prevented early stop
before `target_evals`. Same campaign dir; **v22** = code + learning phase + full production run.

**Purpose.** Close the v21 [cold-start gap](#v21-map-elites-cold-start-grtresna-rejections): MAP-Elites
now **learns from graded pre-GPU rejections** while keeping the main FTL archive pure (`gpu_ok`
only). Continue wormhole QD toward **80 evals** with tuned pipeline knobs from v21.

### Pipeline update: pre-GPU rejection learning

| Layer | Role | Learns from |
|-------|------|-------------|
| **Main archive** (`archive.json`) | Real FTL elites | `gpu_ok` only (unchanged) |
| **Near-miss pool** (`search/pre_gpu/near_miss_pool.py`) | Top-K parents by `score` | `grtresna_rejected`, `solved_ftl_rejected`, `postload_rejected` |
| **Shadow archive** (`pre_gpu_archive.json`) | Separate MAP-Elites grid on pre-GPU axes | Same graded rejections |
| **Unified sampler** (`qd_search/sampling.py`) | 60% gpu_ok elite / 15% shadow / 15% near-miss / 5% feasible / 5% random | Renormalizes when sources empty |

**Shadow descriptor axes** (status-aware, not post-GPU FTL axes):

| Status | x-axis | y-axis |
|--------|--------|--------|
| `grtresna_rejected` | Ham quality | Mom quality |
| `solved_ftl_rejected` | convergence (=1) | precursor tilt (solved FTL) |
| `postload_rejected` | convergence / score | postload margin |

**Excluded from mutation parents:** `grtresna_failed`, `gpu_failed`, `pipeline_interrupted`.

**Gate metrics** now flow into `trajectory.jsonl` (`grtresna_convergence`, `postload_gate`).
Legacy v21 records without structured convergence are rebuilt via `Ham=(\d+)%` regex on
`reason` strings.

**CMA-ES (stage 1):** `--warm-start-include-near-miss` seeds from graded rejections in prior
`trajectory.jsonl` (default on with `--grtresna`). **HQ promotion:** unchanged (replay only).

### v22 launch (resume)

```bash
cd grteclyn-wrapper
BRANCH=wormhole PIPELINE_MONITOR=1 \
  QD_RESUME=1 \
  QD_NAME=general_ftl_wormhole_v21 \
  QD_TARGET_EVALS=80 \
  GPU_IDS="0 1 2 3 4 5 6 7" \
  GPU_SLOTS_PER_DEVICE=1 \
  MAX_CONCURRENT_GRTRESNA=5 \
  BATCH_SIZE=8 \
  QD_ITERATIONS=30 \
  SKIP_QD_PREFLIGHT_TESTS=1 \
  bash scripts/campaigns/general_ftl/run_all.sh \
  > ../runs/_logs/general_ftl_wormhole_v21.launch.log 2>&1 &
```

| Knob | v21 final | **v22 resume** |
|------|-----------|----------------|
| `QD_RESUME` | — | **1** |
| `pre_gpu_learning` | — | **true** (auto with `--grtresna`) |
| `MAX_CONCURRENT_GRTRESNA` | 3 | **5** |
| `max_level` | 2 (metadata) | **1** (`search_common.sh` default) |
| `near_miss_pool_size` | — | **32** |

`max_in_flight = 8 gpu_slots + 5 GRTresna = 13`.

### Progress snapshot (resume start)

- **Inherited from v21:** 26 scored evals — 11 `gpu_ok`, 9 `grtresna_rejected`, 3 `postload_rejected`;
  best near-miss **eval 4** score **−195** (Ham=100%, Mom within threshold).
- **On resume:** near-miss pool + shadow archive rebuild from full `trajectory.jsonl` before
  new submits; children should **mutate near best near-misses** instead of uniform random.
- **Reports:** after each batch of 8 completions, `[qd] report` includes
  `near_miss=N shadow_cov=…`; writes `pre_gpu_archive.json`.

**Monitor:**

```bash
tail -f runs/_logs/general_ftl_wormhole_v21.launch.log
tail -f runs/grtresna_qd/general_ftl_wormhole_v21/trajectory.jsonl
grep near_miss runs/_logs/general_ftl_wormhole_v21.log
watch -n5 'nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv'
```

**Outputs:** `runs/grtresna_qd/general_ftl_wormhole_v21/` — adds `pre_gpu_archive.json` alongside
`archive.json` and `trajectory.jsonl`.

### v22 vs v21

| | v21 | **v22** |
|---|-----|---------|
| GRTresna rejections | logged, ignored for sampling | **near-miss pool + shadow archive** |
| Main archive | `gpu_ok` only | `gpu_ok` only (unchanged) |
| `MAX_CONCURRENT_GRTRESNA` | 3 (starved GPUs) | **5** |
| `max_level` | 2 in run metadata | **1** |
| Campaign dir | `general_ftl_wormhole_v21` | **same** (resume) |

### v22 live pipeline validation (2026-06-18, running)

Validated from live trajectory (54 records at time of check, 3 GPU evolutions + 3 GRTresna
CPU solves running concurrently on GPUs 1-3 at 81-87% utilization).

#### Continuous pipeline — no batch barrier

```
Completion order (post-resume, 25 records):
  eval= 7  pipeline_interrupted  [resumed]
  eval=22  pipeline_interrupted  [resumed]
  ...
  eval=37  grtresna_rejected     [CPU ~30s]
  eval=41  grtresna_rejected     [CPU ~30s]
  eval=46  grtresna_rejected     [CPU ~30s]
  eval=45  grtresna_rejected     [CPU ~30s]
  eval=47  grtresna_rejected     [CPU ~30s]
  eval=52  grtresna_rejected     [CPU ~30s]    ← completes BEFORE eval 35
  eval=35  gpu_ok                [GPU ~2min]   ← GPU evolution still running when 52 finished
  eval=36  gpu_ok                [GPU ~2min]
  eval=38  gpu_ok                [GPU ~2min]
  eval=39  gpu_ok                [GPU ~2min]
  eval=49  postload_rejected     [CPU+solved ~45s]
  eval=43  gpu_ok                [GPU ~2min]
  eval=42  gpu_ok                [GPU ~2min]
  eval=53  grtresna_rejected     [CPU ~30s]
  eval=44  gpu_ok                [GPU ~2min]
```

**Key evidence:** eval 52 (CPU rejection, ~30s) completes before eval 35 (GPU evolution, ~2min).
The pipeline submitted eval 52 while eval 35's GPU phase was still running — no batch await.
**16 out-of-order completions** in trajectory confirm high parallelism.

**At validation time:**

| Phase | Active | Evidence |
|-------|--------|----------|
| GPU evolution | eval_000048, _050, _054, _056, _057 | nvidia-smi: GPUs 1-3 at 81-87%, ~8 GB VRAM each |
| GRTresna CPU solve | eval_000058, _059, _060 | 8 MPI ranks each visible in ps |

9 GPU evolutions + 3 CPU solves running simultaneously. Pipeline continuously feeds
GPU slots without waiting for any batch to complete.

#### Pre-GPU learning active in live run

| Metric | Value |
|--------|-------|
| Shadow archive elites | 3 cells occupied |
| `shadow_improved` events | 8 records tagged |
| Structured `grtresna_convergence` | 7 records carry `{ham_pct, mom_pct}` |
| `postload_gate` metrics | 1 record with `{max_hamiltonian_l2: 0.0124}` |

Shadow archive cell placement:

| Cell | Status | Descriptors | Score |
|------|--------|-------------|-------|
| `(7,0)` | postload_rejected | `conv_quality=1.0, tilt=0` | -100 |
| `(0,7)` | grtresna_rejected (Ham=100%, Mom=0%) | `ham_quality=0.05, mom_quality=1.0` | -195 |
| `(0,0)` | grtresna_failed (NaN convergence) | `quality=0, 0` | -350 |

Near-miss pool top candidates: 4× `postload_rejected` (score -100) rank above
8× `grtresna_rejected` (score -195). This correctly prioritizes candidates that
passed more gates (closer to GPU success).

### Schema: continuous pipeline architecture

```
                         ┌─────────────────────────────────────────────────────────┐
                         │              EvalPipeline (no batch barrier)             │
                         │                                                         │
  ┌──────────┐   submit  │   ┌──────────────┐        ┌──────────────────────┐      │
  │  Unified │──────────►│   │ CpuAdmission │  pass  │      GpuPool         │      │
  │  Sampler │           │   │  Controller  │───────►│  lease() → GPU slot  │      │
  │          │◄──────────│   │  (max=5 MPI) │        │  (8 GPUs × 1 slot)   │      │
  └──────────┘ on_complete   └──────┬───────┘        └──────────┬───────────┘      │
       │                            │                           │                   │
       │                            │ reject                    │ gpu_ok/gpu_failed  │
       │                            ▼                           ▼                   │
       │                   ┌─────────────────────────────────────────┐              │
       │                   │            _on_eval_complete            │              │
       │                   │                                         │              │
       │                   │  ┌─── archive_lock ──────────────────┐  │              │
       │                   │  │ if gpu_ok: archive.insert(elite)  │  │              │
       │                   │  │ if graded_rejection:              │  │              │
       │                   │  │   near_miss_pool.consider(record) │  │              │
       │                   │  │   pre_gpu_archive.insert(shadow)  │  │              │
       │                   │  │ trajectory.append(record)         │  │              │
       │                   │  └───────────────────────────────────┘  │              │
       │                   └─────────────────────────────────────────┘              │
       │                                    │                                       │
       └────────────────────────────────────┘                                       │
         immediately submit next candidate                                          │
         (no wait for batch)                                                        │
                         └─────────────────────────────────────────────────────────┘

  Parallelism at any instant:
    CPU phase:  up to MAX_CONCURRENT_GRTRESNA (5) MPI solves
    GPU phase:  up to gpu_slots (8) evolutions
    Total in-flight: up to 13 evals simultaneously
```

### Schema: pre-GPU rejection learning flow

```
  ┌─────────────────────────────────────────────────────────────────────────────┐
  │                        Eval Completion                                       │
  │                                                                             │
  │   status?                                                                   │
  │     │                                                                       │
  │     ├── gpu_ok ─────────────► Main Archive (archive.json)                   │
  │     │                              │                                        │
  │     │                              ▼                                        │
  │     │                    ┌───────────────────┐                              │
  │     │                    │ Unified Sampler   │                              │
  │     │                    │ 60% mutate elite  │◄─── Main Archive             │
  │     │                    │ 15% mutate shadow │◄─── Shadow Archive           │
  │     │                    │ 15% mutate near   │◄─── Near-Miss Pool           │
  │     │                    │  5% feasible box  │◄─── Union bounds             │
  │     │                    │  5% random        │                              │
  │     │                    └────────┬──────────┘                              │
  │     │                             │                                         │
  │     │                             ▼                                         │
  │     │                      Next candidate                                   │
  │     │                                                                       │
  │     ├── grtresna_rejected ──┐                                               │
  │     ├── solved_ftl_rejected ┼──► Near-Miss Pool (top-32 by score)           │
  │     ├── postload_rejected ──┘         │                                     │
  │     │                                 └──► Shadow Archive (pre_gpu_archive) │
  │     │                                      status-aware descriptors:        │
  │     │                                        grtresna: (ham_q, mom_q)       │
  │     │                                        solved_ftl: (1.0, tilt)        │
  │     │                                        postload: (conv_q, margin)     │
  │     │                                                                       │
  │     ├── grtresna_failed ────► trajectory only (no gradient signal)          │
  │     ├── gpu_failed ─────────► trajectory only                               │
  │     └── preflight_rejected ─► trajectory only                               │
  └─────────────────────────────────────────────────────────────────────────────┘

  Score ranking in near-miss pool (higher = closer to GPU success):
    postload_rejected:    -75 to -100  (passed convergence + solved-FTL)
    solved_ftl_rejected:  -75 to -95   (passed convergence)
    grtresna_rejected:    -100 to -350 (failed convergence, graded by Ham/Mom)
    grtresna_failed:      excluded     (solver crash, no useful gradient)

  Cold start (no gpu_ok yet): renormalize → 43% shadow + 43% near-miss + 14% random
```

### Schema: CMA-ES near-miss warm-start (QD→CMA handoff)

```
  Prior QD trajectory.jsonl
        │
        ├── gpu_ok records ──────────► _load_warm_start_vectors (top_k=5)
        │                                     │
        │                                     ▼
        │                              ┌──────────────┐
        └── graded rejections ────────►│ near_miss_   │──► Append to warm vectors
             (--warm-start-include-    │ vectors_from │    (dedup by 8-decimal rounding)
              near-miss, default on)   │ _trajectory  │
                                       │ (near_miss_k │
                                       │  = 4)        │
                                       └──────────────┘
                                              │
                                              ▼
                                    CMA-ES initial population
                                    (gpu_ok first, then near-miss)
                                    First vector = initial mean
                                    Others jittered by sigma
```

### v22 final results: top 3 + FTL champions

**Campaign summary** (`runs/grtresna_qd/general_ftl_wormhole_v21/`):

| Metric | Value |
|--------|-------|
| Evals | **200** (201 trajectory lines incl. 10 `pipeline_interrupted` from prior halt) |
| `gpu_ok` | **103** |
| Archive | **4 elites**, coverage **6.25%**, MAP-Elites best **165.6** (eval 063); **refinement champion eval 191** |
| Throughput (resume 98→200) | ~**0.93 eval/min** (8 GPU × 1 slot, `MAX_CONCURRENT_GRTRESNA=3`) |
| Near-miss pool | 32 (pre-GPU learning active throughout) |

**Post-QD CMA-ES:** seeded from eval **063** → champion **046** at **179.8** (+14.2 pts); see
[v22 CMA-ES section](#v22-cma-es-wormhole-refinement-general_ftl_wormhole_cmaes_v1-2026-06-18).

**Top 3 by score** (final run):

| Eval | Score | `ftl_geo_evolving` | `f_geo_peak` | `operational_ftl` | Survival | Reading |
|------|------:|-------------------:|-------------:|------------------:|---------:|---------|
| **063** | 165.6 | **19.3%** | 4.2% | 0 | 0.94 | Score + `ftl_geo_evolving` record holder (early run) |
| **191** | **161.9** | **18.5%** | 3.8% | 0 | **1.00** | **Champion** — resume-era lead; perfect survival, `f_op_peak` holder, jitter 0.96 |
| **174** | 157.4 | 18.5% | 3.8% | 0 | 1.00 | Stable variant — thick shell (1.40), moderate jitter |

All three: `operational_ftl = 0` but strong **4D null-geodesic** shortcuts — the hallmark of
stationary, non-translating wormholes (see [physical validation](#v22-physical-validation--next-steps)).

Earlier basin members still in archive: **eval 071** (122.8, symmetric/stable), **eval 040** (96.5, compact/massive).

**FTL champions** (`ftl_champions.json` — one eval dir kept per peak metric):

| Metric | Eval | Peak value | Notes |
|--------|------|------------|-------|
| `ftl_geo_evolving` | **063** | **19.3%** | Run-best integrated 4D shortcut |
| `f_op_peak` | **191** | 0.0178 | Best operational precursor |
| `f_geo_evol` (frozen) | **178** | 5.2% | Best static-slice geodesic peak (resume find) |
| `max_local_speed` | **178** | **1.31** | Highest local coordinate speed |
| `superluminal_fraction` | **050** | 77.4% | Largest superluminal volume fraction |
| `ftl_lifetime_fraction` | **008** | 1.00 | Longest FTL persistence |

### v22 physical validation & next steps

**Verdict.** First production run of the pipelined `general_ftl` wormhole campaign in a
**15-D** pinned search space successfully located a **clustered family** of stable,
constraint-clean, non-translating wormholes. Empirical proof that coordinate warp-bubble
metrics and stationary FTL wormhole metrics are physically distinct observables.

#### I. Scoring paradox: `operational_ftl = 0` while `ftl_geo_evolving ≈ 19%`

In a **stationary** wormhole (\(\beta^i \approx 0\)):

1. **Lapse-dominated coordinate speed** — heavy throat matter drops \(\alpha < 1\), so
   \(c_{\rm coord} = \alpha \pm \beta < 1\). Frozen-slice Dijkstra (`operational_ftl`) reads
   subluminal → **0**.
2. **Proper-distance contraction** — bipolar shell geometry warps \(\gamma_{ij}\); proper
   path length through the throat is shorter than coordinate distance. The **4D null-geodesic
   tracer** integrates \(ds^2 = 0\) and captures this → **~19% faster** than flat-space light.

The pipeline has decoupled coordinate gauge artifacts from physical traversable shortcuts.

#### II. Dynamically opening throat

For eval **191**: integrated `ftl_geo_evolving` = **18.5%** but frozen `f_geo_peak` = **3.8%**
at \(t \approx 16\) (eval **063** shows the same pattern at 19.3% / 4.2%). Even with zero
initial velocities, CCZ4 evolution is highly dynamical — negative-energy shell and positive
central regions interact gravitationally. The throat **opens during evolution**; a photon
emitted at \(t=0\) traverses geometry that is actively warping. Only the 4D tracer detects
this time-dependent GR effect.

#### III. Elite basin (cell [1, 7]) — final run

| Eval | Role | Key traits |
|------|------|------------|
| **191** | **Champion** | Resume-era lead — extreme jitter (0.96), dipole +0.57, survival **1.00**, 18.5% 4D shortcut, `f_op_peak` holder |
| **174** | Stable variant | Thick shell (1.40), moderate jitter (0.59) → symmetric, survival 1.00, 18.5% 4D shortcut |
| **063** | Record holder | Early-run peak score (165.6) + `ftl_geo_evolving` record (19.3%); survival 0.94, 52% superluminal volume |

Legacy early-basin members: **071** (symmetric/stable, 13.9% 4D), **040** (compact/massive, survival 0.77).

#### IV. Strategic action plan

| Step | Action | Goal | Status |
|------|--------|------|--------|
| **1** | **CMA-ES refinement** seeded from QD **eval 063** | Local hill-climb in the v22 wormhole basin; push 4D shortcut above QD record while holding survival. | **Done** — [eval 046 @ 179.8](#v22-cma-es-wormhole-refinement-general_ftl_wormhole_cmaes_v1-2026-06-18) (+14.2 vs 063). Stopped gen 6 (48 evals). |
| **2** | **HQ promotion** of CMA-ES **eval 046** → **N=256, L=128, ml=3, t=30** | Does the refined throat survive full resolution and longer time? | **Done** — [HQ final results](#hq-eval-046-final-results-t30): peak **7.57%** 4D @ t≈15.6, channel dies + horizon **−546** @ t=30. Static shell confirmed; structure **44%** persistence at end. |
| **3** | **Complex scalar (boson star) pivot** | Real scalars disperse at \(t \gg 50\). Re-run with complex fields + conserved U(1) charge. | **Phase 1 done** — [smoke + v3 unpinned QD](#boson-star-unpinned-qd--frames-v3-2026-06-18) (12 evals, 8 gpu_ok); larger campaigns + multi-site bosons pending |

**Note:** Action plan originally named eval **191** as refinement seed (resume-era champion, survival 1.00). CMA-ES used **063** instead — the absolute score / `ftl_geo_evolving` record holder — which replayed at **165.6** bit-identically once scoring parity was fixed.

---

## v22 CMA-ES: wormhole refinement (`general_ftl_wormhole_cmaes_v1`, 2026-06-18)

**Status:** **stopped early** (user halt @ gen 6, 48 evals; score plateau ~179.8).
Campaign dir (pruned): `runs/grtresna_cmaes/general_ftl_wormhole_cmaes_v1/` — champion **`eval_000046` only**.

**Purpose.** Stage-1 local refinement after [v22 QD](#v22-final-results-top-3--ftl-champions): hill-climb
within the clustered wormhole basin found by MAP-Elites, using the same 15-D pinned search space,
`OBJECTIVE_MODE=general_ftl`, and pipelined evaluator as QD.

### Seed & launch

Warm-started from QD trajectory `general_ftl_wormhole_v21` with `WARM_START_TOP_K=1` (exact genome
of eval **063**, not resume-era eval 191):

```bash
cd grteclyn-wrapper
RUN_NAME=general_ftl_wormhole_cmaes_v1 \
OBJECTIVE_MODE=general_ftl \
WARM_START_TRAJECTORY="${GRTECLYN_ROOT}/runs/grtresna_qd/general_ftl_wormhole_v21/trajectory.jsonl" \
WARM_START_TOP_K=1 WARM_START_JITTER=0.05 SIGMA0=0.05 \
TARGET_EVALS=150 MAX_GENERATIONS=50 KEEP_TOP_EVAL_DIRS=3 \
GPU_IDS="0 1 2 3 4 5 6 7" GPU_SLOTS_PER_DEVICE=1 MAX_CONCURRENT_GRTRESNA=3 \
PIN_DIMS="$(bash -c 'source scripts/campaigns/lib/general_ftl_pins.sh && ftl_general_ftl_wormhole_pins')" \
  bash scripts/campaigns/cmaes/run.sh
```

| Knob | Value |
|------|-------|
| Search space | 15-D wormhole pins (same as v22 QD) |
| Population | 8 (= GPU count) |
| σ₀ | 0.05 |
| Pipeline | `GpuPool` + `EvalPipeline`, `max_concurrent_grtresna=3` |
| Target | 150 evals (stopped @ gen 6, 48 evals) |

### Scoring parity fix (critical)

First CMA-ES attempt replayed QD eval 063 at **−27** instead of **165.6** — same overrides,
different score. Root cause: CMA-ES path lacked QD's geodesic configuration.

| Issue | QD (v22) | CMA-ES (before fix) | Fix |
|-------|----------|---------------------|-----|
| `GRTECLYN_GEO_DIRECTIONS` | `x y z` (`general_ftl/run_all.sh`) | default **`x` only** | export in `campaigns/lib/search_common.sh` when `OBJECTIVE_MODE=general_ftl` |
| Eval driver | shared `_run_cpu_grtresna_gates` + `_run_gpu_session` | duplicate `_objective` path | CMA-ES GRTresna mode wired to shared `evaluate_overrides` stack |
| Live trajectory | append per eval | batch rewrite per gen | `_track_trajectory` + `live_trajectory_path` |

After fix, warm-start replay (CMA-ES eval **001**) scored **165.608** — **bit-identical** to QD eval **063**
(`ftl_geo_evolving` = 0.192631, `f_geo_evol` = 0.0419, geodesic axis **z**).

HQ promotion received the same fix: `GRTECLYN_GEO_DIRECTIONS=x y z` + `--objective-mode general_ftl`
in `campaigns/lib/promote_common.sh` and `campaigns/hq/run_batch.sh`.

### Improvement over QD top candidate

| Metric | QD eval **063** | CMA-ES eval **046** | Δ |
|--------|----------------:|--------------------:|---|
| **Score** | **165.6** | **179.8** | **+14.2 (+8.6%)** |
| `ftl_geo_evolving` | 19.3% | **20.3%** | +1.0 pp |
| `f_geo_evol` (4D trace) | 4.19% | 4.15% | ~flat |
| `f_geo_peak` (frozen) | 7.07% | 6.77% | −0.3 pp |
| `ftl_lifetime_fraction` | 86% | 86% | — |
| `max_local_speed_peak` | 1.146 | 1.140 | −0.006 |
| Survival | 0.94 | **1.00** (structural 0.998) | improved |
| Geodesic axis | **z** | **z** | same |

CMA-ES did **not** discover a new mechanism class — it tightened the existing v22 wormhole basin:
higher integrated 4D score (`ftl_geo_evolving` × survival weight) with perfect structural persistence.

### Convergence curve (best per generation)

| Gen | Best score | Notes |
|-----|----------:|-------|
| 1 | 174.4 | Warm seed + σ=0.05 neighbors; eval 001 replays 063 @ 165.6 |
| 2 | 174.8 | |
| 3 | 178.0 | |
| 4 | 179.1 | |
| 5 | 179.6 | Population mean ≈ 178 — basin collapsed |
| 6 | **179.8** | **Champion eval 046**; last 4 gen bests spread only +1.8 |

Stopped manually when further gens were unlikely to beat 179.8 (σ=0.05 jitter around a single optimum).

### v22 CMA-ES HQ promotion

One HQ replay suffices — top CMA-ES evals 039/045/046 were σ-neighbors in the same basin (scores
179.6–179.8), not distinct mechanisms.

```bash
cd grteclyn-wrapper
SOURCE_RUN="${GRTECLYN_ROOT}/runs/grtresna_cmaes/general_ftl_wormhole_cmaes_v1" \
NAME_PREFIX=general_ftl_wormhole_cmaes_v1 \
OBJECTIVE_MODE=general_ftl \
CANDIDATES="46 0" \
  bash scripts/campaigns/hq/run_batch.sh
```

Output: `runs/grtresna_promote/general_ftl_wormhole_cmaes_v1_hq_eval000046/` (L=128, N=256, ml=3, t=30,
4D mode **`hq`**, dirs **`x y z`**, frames **on**, incremental `score_timeseries.jsonl`).

### HQ eval 046: final results (t=30)

**Status:** **complete** (2026-06-18). Run dir:
`runs/grtresna_promote/general_ftl_wormhole_cmaes_v1_hq_eval000046/`.

Fresh GRTresna solve converged (Ham **0.00026%**, Mom **0%**). GPU evolution finished **t=30**
(`simulation_exit_code=0`, ~105 min wall). Incremental consumer + frame rendering complete.

#### Score: search vs HQ (not directly comparable)

| Stage | Eval | Score (`general_ftl`) | Grid / time | 4D profile |
|-------|------|----------------------|-------------|------------|
| **CMA-ES search** | **046** | **+179.8** | N=128, L=64, ml=1–2, **t=16** | `search` (3-ray stack) |
| **HQ promotion** | **046** | **−546.3** @ t=30 | N=256, L=128, ml=3, **t=30** | `hq` (full verify) |

HQ totals come from `small_data/score_timeseries.jsonl` (incremental `general_ftl` with horizon
penalty). Recomputed end-of-run metrics agree within ~10 pts (**−556** from full metric stack).

**Do not rank HQ −546 against search +179.8** — different resolution, stop time, geodesic mode, and
HQ applies cumulative horizon / exotic penalties over **t=30** that search never sees at **t=16**.

#### FTL timeline (HQ incremental)

| Phase | Sim time | Score | `f_geo` / peak | Notes |
|-------|----------|-------|----------------|-------|
| Throat opening | t ≈ 0–8 | −12 → −45 | peak builds to **~7.2%** | Best incremental window |
| **Peak 4D shortcut** | **t ≈ 15.6** | **−47.1** | **`f_geo_peak` ≈ 7.57%** | Best score + best geodesic peak |
| Horizon trigger | **t ≈ 21.1** | **−520** | `horizon_penalty → −1` | Score cliff (same pattern as [v16 HQ elites](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17)) |
| Channel death | t ≈ 27–30 | −538 → **−546** | instant **`f_geo` → 0%** | 4D channel gone by end |
| **Final @ t=30** | **30.0** | **−546.3** | peak **7.57%**, instant **0%** | `max_local_speed ≈ 1.04` |

End-of-run **`small_data/evolving_geodesic.json`**: **`f_geo = 0.0`** (end-to-end HQ trace on **x** —
no surviving shortcut at t=30). Search-stage eval 046 had **`ftl_geo_evolving ≈ 20.3%`** at t=16 on
the cheaper `search` profile; HQ **does not verify** that headline at full resolution / longer time.

#### Structural / penalty breakdown @ t=30

| Metric | HQ @ t=30 | CMA-ES search @ t=16 |
|--------|-----------|----------------------|
| `structural_persistence` | **0.44** | **0.998** |
| `horizon_penalty` | **−1.0** | **0.0** |
| `exotic_penalty` | **−1.6** (saturated) | **−1.57** |
| `max_local_speed` | **1.04** | **1.14** (evolved peak) |
| `ftl_geo_evolving` | **0.0** | **0.203** |

**Verdict:** The v22 static-matter wormhole **opens a real 4D throat mid-run** (peak **7.57%** —
comparable in order to [eval 144 HQ](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17) **7.96%**)
but **does not persist** to t=30 at ml=3 / 256³. Horizon penalty kills the leaderboard total; this
genome is a **mid-run FTL demonstrator**, not a t=30 survivor.

#### Movies

15 field movies (126 frames each, 10 fps) under:

`runs/grtresna_promote/general_ftl_wormhole_cmaes_v1_hq_eval000046/movies/`

Key diagnostics: `movie_local_speed_z.mp4`, `movie_chi_z.mp4`, `movie_shift1_z.mp4`,
`movie_lump_activity_proj_z.mp4`.

```bash
cd grteclyn-wrapper
bash scripts/plot/make_movies.sh \
  ../runs/grtresna_promote/general_ftl_wormhole_cmaes_v1_hq_eval000046 \
  --framerate 10
# → writes .../movies/movie_<field>_<axis>.mp4
```

#### Artifacts

| File | Content |
|------|---------|
| `small_data/score_timeseries.jsonl` | Incremental score / components vs sim time |
| `small_data/ftl_timeseries.dat` | FTL probe time series |
| `small_data/evolving_geodesic.json` | Final HQ 4D trace (`f_geo=0`) |
| `metadata.json` | GRTresna convergence, solved-geometry FTL @ t=0 grid |
| `movies/*.mp4` | Frame movies (15 fields) |

Inspect:

```bash
tail -1 runs/grtresna_promote/general_ftl_wormhole_cmaes_v1_hq_eval000046/small_data/score_timeseries.jsonl | python3 -m json.tool
```

### HQ eval 046: matter dynamics

**Status:** **complete** (2026-06-18) — see [final results](#hq-eval-046-final-results-t30) for scores
and FTL timeline. Matter class unchanged from search: **static shell**, no lump locomotion. This is
**not** a moving-matter / frame-dragging wormhole — it is the [stationary wormhole basin](#v22-physical-validation--next-steps) where FTL comes from **geometry + 4D null trace**, not from lumps that translate or spin.

#### Static matter configuration (no lump locomotion)

| Override | HQ value | Effect |
|----------|----------|--------|
| `grtresna_shell_static` | **1.0** | Forces all lump velocities and `omega` to **zero** (`build_grtresna_config`) |
| `grtresna_shell_*_velocity` | **0** | No toroidal / poloidal / radial matter currents |
| `grtresna_shell_omega` | **0** | No spin |
| `grtresna_shell_lumps` | **5** | Shell ansatz — bipolar ±z lobes + throat (fixed centres from GRTresna solve) |
| `grtresna_shift_seed` | **0.55** | Initial shift / frame-drag **seed** (geometry channel, not bulk matter motion) |
| `grtresna_shell_radial_jitter` | **0.82** | Asymmetric shell placement (same basin as QD/CMA-ES) |
| `grtresna_shell_dipole_amp` | **+0.59** | Dipole asymmetry on the static shell |

Every resolved lump in `grtresna/params.txt` has `lump*_velocity = 0 0 0` and `lump*_omega = 0`.
Search-stage eval **046** already reported `comoving.stationary: true` and `scalar_pi_range: 0.0`
at t=16 — HQ inherits that physics class.

#### What evolves vs what is frozen

| Frozen during evolution | Evolves during evolution |
|-------------------------|--------------------------|
| Lump **centres** (shell geometry from elliptic solve) | **Spacetime** — \(\chi\), lapse, throat opening ([§ II](#ii-dynamically-opening-throat)) |
| Bulk linear / angular **momentum** (\(\Pi \approx 0\)) | **Scalar amplitudes** at fixed `phi_lump*` / `Pi_lump*` slots |
| Translating / rotating cloud dynamics | **Shift** seeded at t=0; metric shear through the throat |
| Moving-matter FTL channel | **4D null-geodesic shortcut** on **z** (`ftl_geo_evolving`; HQ verify mode) |

The HQ frame fields `lump_activity`, `phi_lump_sum`, and `Pi_lump_sum` are **diagnostic overlays**,
not evidence of lump translation:

- **`lump_activity`** = \(\sum_k \sqrt{\phi_k^2 + \Pi_k^2}\) summed over per-lump fields at **fixed**
  grid locations — shows breathing / redistribution, not bulk drift.
- **`Pi_lump_sum`** stays flat (consumer z-lims \(\sim 10^{-13}\)) — no momentum transport.
- **`phi_lump_sum`** shows modest amplitude variation (\(\lesssim 0.01\) early in the run) while the
  throat geometry evolves in `chi` / `local_speed` frames.

Do **not** interpret `phi_lump_sum_z` movies as “lumps moving through the box”; interpret them as
**field activity on a static shell** while the [dynamically opening throat](#ii-dynamically-opening-throat)
develops.

#### HQ answers (t=30)

Search stopped at **t=16**, N=128. HQ results for the static-matter throat:

1. **Persistence to t=30** — **partial failure.** `structural_persistence` falls from **0.998** (search)
   to **0.44** (HQ @ t=30). Throat geometry evolves through mid-run ([movies](#hq-eval-046-final-results-t30))
   but the 4D channel dies before the end.
2. **Score vs CMA-ES 179.8** — **does not hold.** Final HQ incremental total **−546** (horizon kill @
   t≈21). Peak incremental score **−47** @ t≈15.6 with **`f_geo_peak` 7.57%** — use that window for
   physical comparison, not the t=30 total (same lesson as [v16 HQ](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17)).
3. **Stationary class** — **confirmed.** `grtresna_shell_static=1`, all lump velocities/omega **0**;
   `Pi_lump_sum` flat; no centre drift. Mechanism remains **geometry + 4D null trace**, not moving matter.

Moving-matter tests (`grtresna_shell_static=0`, non-zero shell velocities) require a **new** search
genome; this HQ run is explicitly **not** that experiment.

---

## v21: Pipelined QD + GPU tenancy tuning (2026-06-17)

**Status:** **stopped early** (user halt at ~30 evals; **26** trajectory lines saved).
Campaign dir preserved: `runs/grtresna_qd/general_ftl_wormhole_v21/`.

**Purpose.** Validate the pipelined MAP-Elites evaluator (`search/gpu_pool.py`,
`search/eval_pipeline.py`) on production grid/time: GRTresna CPU solves overlap with
GPU `main3d` evolution across **8× H100**, wormhole-only pins (same as v20 wormhole).

### Architecture (unchanged from design)

| Component | Role |
|-----------|------|
| `GpuPool` | `total_slots = len(gpu_ids) × gpu_slots_per_device`; blocking lease before GPU phase |
| `EvalPipeline` | Cross-batch queue; no batch barrier; one completion → one new submit |
| `CpuAdmissionController` | Caps concurrent GRTresna MPI solves (`MAX_CONCURRENT_GRTRESNA`) |
| QD driver | Pipelined MAP-Elites (`use_pipeline=True` default) |

`max_in_flight = gpu_slots + max_concurrent_grtresna` (e.g. 8 + 3 = 11 with final config).

### Run history (three launches, same campaign name)

#### Attempt 1 — multi-tenant GPU overload (failed throughput)

```bash
GPU_SLOTS_PER_DEVICE=5   # 8×5 = 40 GPU leases
MAX_CONCURRENT_GRTRESNA=6
```

**Result:** Multiple `main3d` sessions stacked on the same H100 (GpuPool fills GPU 0
first). At `stop_time=16` each solo evolution peaks ~**44 GB** VRAM; **5 concurrent**
on one GPU caused severe slowdown — wall time per eval ~**3×** vs solo (~35 min vs
~12 min). **Do not use 5 slots/GPU at t=16** on 80 GB cards with this gridinit.

#### Attempt 2 — 1 slot/GPU, high GRTresna concurrency

```bash
GPU_SLOTS_PER_DEVICE=1   # 8 GPU leases total
MAX_CONCURRENT_GRTRESNA=6
```

**Result:** GPU memory healthy (~one evolution per device), but GRTresna admission
still contended with long GPU phases.

#### Attempt 3 — production config (successful pipeline shape)

```bash
GPU_SLOTS_PER_DEVICE=1
MAX_CONCURRENT_GRTRESNA=3
BATCH_SIZE=8
QD_TARGET_EVALS=80
```

**Result:** **8 GPUs each running one evolution**; **up to 3 GRTresna solves in flight**
continuously feeding the pipeline. First `gpu_ok` after ~5 min wall; archive began
filling. Typical `gpu_ok` eval ~5–8 min; one AMR outlier (`eval_000007`, `max_level=2`)
hit ~**39 GB** VRAM and ~30 min wall.

**Trajectory (26 logged):** 11 `gpu_ok`, 9 `grtresna_rejected`, 3 `postload_rejected`,
3 `grtresna_failed`. Best scored `gpu_ok`: eval **19** score **+33.7**. Best GRTresna
near-miss: score **−195** (Ham=100%, Mom within threshold) — logged but **not used
for sampling** (see below).

### Issue: GPU idle gaps — GRTresna batch too small

With `gpu_slots=8` and `MAX_CONCURRENT_GRTRESNA=3`, the pipeline often had **3–5 GPUs
idle** for multi-minute stretches: GPU work finished faster than new gridinits arrived.
The evaluator submits one new candidate per completion (no batch barrier), so throughput
is limited by how fast GRTresna can produce passing initial data.

**Fix:** **Increase `MAX_CONCURRENT_GRTRESNA`** (GRTresna batch / CPU admission cap) to
**5–6** so more solves run in parallel while GPUs evolve — target `max_in_flight ≈ 13–14`
on 8 GPUs. Do **not** compensate by raising `GPU_SLOTS_PER_DEVICE` at t=16.

| Knob | Attempt 1 (bad) | Attempt 3 (ran) | **Recommended relaunch** |
|------|-----------------|-----------------|--------------------------|
| `GPU_SLOTS_PER_DEVICE` | 5 | 1 | **1** |
| `MAX_CONCURRENT_GRTRESNA` | 6 | 3 | **5–6** |
| `gpu_slots` | 40 | 8 | **8** |
| GPU util @ t=16 | overloaded, slow | healthy, GRTresna-starved | balanced |

### VRAM sizing (eval_000005 gridinit replays, H100 80 GB)

Benchmark (`scripts/benchmarks/gpu_gridinit_load.sh`) on real 128³ wormhole gridinit —
not bare `sweep` (~5 GB underestimate).

| Concurrent `main3d` | `stop_time` | Peak VRAM | Notes |
|--------------------|-------------|-----------|-------|
| 1 | 4 | ~9.7 GB | matches QD evaluate path |
| 2 | 4 | ~17.4 GB | ~8.7 GB/eval, linear |
| 5 | 4 | ~52 GB | fits 80 GB at short t only |
| 1 | **16** | **~44 GB** | AMR; **1 slot/GPU** at t=16 |
| 5 | **16** | — | **avoid** — contention, not OOM |

### v21 MAP-Elites cold start (GRTresna rejections)

**Bug / design gap:** MAP-Elites only inserts `gpu_ok` into `archive.json`. GRTresna
rejections are written to `trajectory.jsonl` with **graded fitness** (CMA-ES uses this;
QD does not). Until the first `gpu_ok`, `_sample_next_candidate` draws **uniform random**
vectors — ~40% of v21 evals were pre-GPU rejects with no learning signal.

Rejections also bin to descriptor cell `(0,0)` because FTL axes need GPU evolution.

**Fixed in [v22](#v22-pre-gpu-rejection-learning--v21-resume-2026-06-18):** near-miss parent pool
+ shadow pre-GPU archive (`search/pre_gpu/`). Main `archive.json` stays `gpu_ok` only; graded
rejections guide sampling via `NearMissPool` and `pre_gpu_archive.json`.

### Recommended relaunch

See [v22 resume command](#v22-launch-resume). v21 stopped campaign preserved at
`runs/grtresna_qd/general_ftl_wormhole_v21/`.

**Monitor:**

```bash
tail -f runs/_logs/general_ftl_wormhole_v21.launch.log
tail -f runs/grtresna_qd/general_ftl_wormhole_v21/trajectory.jsonl
watch -n5 'nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv'
```

**Outputs:** `runs/grtresna_qd/general_ftl_wormhole_v21/` — `general_ftl` objective,
`ftl_lifetime` descriptors, `GRTECLYN_GEO_DIRECTIONS=x y z`.

### v21 vs v20

| | v20 | **v21** |
|---|-----|---------|
| Evaluator | batch barrier | **pipelined** CPU→GPU, no batch wait |
| GPU tenancy | 1 evol / GPU | **1 evol / GPU** (after tuning; not 5) |
| GRTresna | serialized per wave | **continuous** solves (`MAX_CONCURRENT_GRTRESNA`) |
| Campaign | 3-class parallel | **wormhole-only** pipeline validation |
| Target | 30 iter × 3 classes | **80 evals** single class |

---

## v10–v20: pipeline evolution & runs (2026-06-11 — 2026-06-17)

Single chronicle for scoring/metric/probe campaigns **v7–v20** (named `ftl_discovery_v*` through
`general_ftl_*`). **Pipeline changes** and **run outcomes** only — no per-version chapters.

**Main pipeline improvements**

| Era | Date | Pipeline change | Code / config |
|-----|------|-----------------|---------------|
| **v4** | 06-11 | Persistence-honest `survival`; FTL shaping × `structural_persistence`; mass search + velocity caps | `metrics/score/`, 18-D space, `STOP_TIME` 8 |
| **v7→v8** | 06-11 | Null-ray **forward** launch (`future_null_cov`); relative H-drift gate `H_REL_TOL=1e-2` | `metrics/probes/ftl/geodesic.py` |
| **v9** | 06-11 | Geodesic ramp target 5%→**20%**; weight ×1500→×1000; rebalance `ftl_first` (coordinate shaping down) | `metrics/score/ftl.py` |
| **v10** | 06-11 | HQ killed v9 shortcuts → `STOP_TIME` 8→**16**; searchable `shell_static` toggle | `campaigns/qd/run.sh` |
| **v11** | 06-12 | Geodesic reward × `structural_persistence`; exotic/energy weights **100× / 40×** | `metrics/score/ftl.py` |
| **v13–v14** | 06-12 | λφ⁴; matter layouts 0–4; per-lump profile + cloud; geodesic contradiction gate; **23-D** space | GRTresna + `spaces.py` |
| **v15** | 06-13 | Per-frame `ftl_timeseries.dat`; time-mean geodesic; descriptor **`ftl_lifetime`** | `consume_plotfiles/`, `qd_search/descriptors.py` |
| **v16** | 06-13 | FTL champion retention (`ftl_retention.jsonl`); horizon penalty needs **lapse corroboration** | `search/ftl_retention.py`, `collapse.py` |
| **v17** | 06-14 | CMA-ES **`robust_ftl`** warm-start from healthy QD survivors (not raw king) | `search/optimize/`, `ftl_cmaes_v17_robust` |
| **HQ batch** | 06-15 | Incremental `score_timeseries.jsonl`; HQ ml=3 `regrid_interval` fix | `metrics/aggregation/incremental.py` |
| **4D probe** | 06-15→16 | End-to-end null trace through metric stack; `ftl_geo_evolving` headline; search vs HQ profiles | `evolving_geodesic.py`, `evolving_geodesic_options.py` |
| **v18 / ftl_4d** | 06-16→17 | Full **QD → CMA-ES → HQ** on 4D metric; matched `ftl_first` objective | `ftl_4d_v1`, `ftl_4d_cmaes_v1` |
| **v20** | 06-17 | **`general_ftl`** objective; `--pin-dimension`; `GRTECLYN_GEO_DIRECTIONS=x y z` | `objectives.py`, `general_ftl/run_all.sh` |

**4D frozen vs evolving (design lesson):** eval **086** smoke 4D **1.42%** vs frozen **5.75%**; HQ t=30
4D **0%** (negative control); eval **144** HQ 4D **7.96%** verified (positive control). Frozen
mid-run peaks overstate dynamic warps — **4D trace is authoritative** when enabled.

**HQ lesson (v16 elites):** trust **incremental peak** (~t≈9–12), not final t=30 — horizon −500 kills
static-lump winners; eval **177** (dynamic, CMA-ES) only finisher positive (+67).

**Alcubierre control (06-12):** probes validated ~**32%** f_geo @ 129³; QD H-gate fixed with
`GEO_REFINE_N=129` re-probe.

**Next directions (06-15):** persistent standing channel vs transport worldtube; exotic-energy Pareto;
boson-star matter; CMA-ME; analytic metric extraction.

**Campaign runs (v7–v20)**

| Run | Date | Evals / gpu_ok | Best eval | Score / f_geo | Headline | Run dir |
|-----|------|----------------|-----------|---------------|----------|---------|
| `ftl_discovery_v7` | 06-11 | 88 / 53 | 71 | 606 / geo=0 | Coordinate precursors; geodesic blind (fixed v8) | `runs/grtresna_qd/ftl_discovery_v7/` |
| `ftl_discovery_v9` | 06-11 | 88 / 54 | 40 | 266 / **2.6%** | First certified geodesic shortcuts | `runs/grtresna_qd/ftl_discovery_v9/` |
| `ftl_discovery_v10` | 06-11 | — | — | — | HQ rejected all v9 promotes; longer QD window | — |
| `ftl_discovery_v11` | 06-12 | 400 | — | — | Persistence-gated geodesic + physicality pressure | `runs/grtresna_qd/ftl_discovery_v11/` |
| `ftl_discovery_v13` | 06-12 | 278 | — | — | λφ⁴ + layouts; zero geodesic until gate | `runs/grtresna_qd/ftl_discovery_v13/` |
| `ftl_discovery_v14` | 06-12 | 504 / 351 | **231** | 551 / **5.30%** | Ring layout dominates; Alcubierre ~6× higher f_geo but metric-first | `runs/grtresna_qd/ftl_discovery_v14/` |
| `ftl_discovery_v15` | 06-13 | — | 231 | peak **7.43%** @ t=9.6 | Time-resolved scoring validated | `runs/grtresna_qd/ftl_discovery_v15/` |
| `ftl_discovery_v16` | 06-13 | ~971 | **233** | 652 / 5.88% peak | FTL retention; horizon fix; plateau | `runs/grtresna_qd/ftl_discovery_v16/` |
| `ftl_cmaes_v17_robust` | 06-14→15 | 200 / 163 | **177** | 312 / **5.65%**, timeavg 16.3% | CMA-ES +0.26 pp f_geo vs seed 739 | `runs/grtresna_cmaes/ftl_cmaes_v17_robust/` |
| HQ v16+v17 | 06-15 | 4 promotes | **233** incr. | **749** @ t≈11.8 | Only **177** final +67; 233/446/676 horizon-killed | `runs/grtresna_promote/l128n256t30_*/` |
| `ftl_max_speed_no_penalty_v1` | 06-15 | 200 / 100 | **92** | +27.5 / **27.5%** timeavg | Side survey; max **1.58 c** coord (eval 70); not v16-comparable | `runs/grtresna_qd/ftl_max_speed/ftl_max_speed_no_penalty_v1/` |
| `ftl_4d_v1` | 06-16 | 192 / 105 | **156** | **508** / `ftl_geo_evol` 0.346 | First 4D-in-loop QD | `runs/grtresna_qd/ftl_4d/ftl_4d_v1/` |
| `ftl_4d_cmaes_v1` | 06-16 | 144 / 140 | **144** | **596** / 0.395 | +88 pts vs QD 156 | `runs/grtresna_cmaes/ftl_4d_cmaes_v1/` |
| HQ eval **144** | 06-17 | 1 | **144** | **283** / 4D **7.96%** | **Verified** gauge-invariant shortcut (5/5 rays) | `runs/grtresna_promote/l128n256t30_ftl_4d_cmaes_qd_eval000144/` |
| `general_ftl_{wormhole,ring,spin}` | 06-17 | 172 / 130 | ring **43** | **196** / search **~3.9%** | Ring **20** 4D hits; wormhole **0**; stopped early, top-3 kept | `runs/grtresna_qd/general_ftl_*/` |

**Eval 177 physics (HQ, 06-15):** null transit **~5.9%** early (`≈1.063 c` signal); matter sub-luminal;
exotic total **~5–24× < Alcubierre**, per-shortcut energy **≈ comparable**. Not transport — pulsing lens.

**v20 launch:**

```bash
cd grteclyn-wrapper
MODE=par QD_ITERATIONS=30 \
  bash scripts/campaigns/general_ftl/run_all.sh \
  > ../runs/_logs/general_ftl_par.launch.log 2>&1 &
```

**v18 QD → CMA-ES → HQ (representative):**

```bash
# QD
RUNS_DIR="${GRTECLYN_ROOT}/runs/grtresna_qd/ftl_4d" QD_NAME=ftl_4d_v1 QD_TARGET_EVALS=200 \
  bash scripts/campaigns/qd/run.sh
# CMA-ES
RUN_NAME=ftl_4d_cmaes_v1 WARM_START_TRAJECTORY=.../ftl_4d_v1/trajectory.jsonl \
  WARM_START_TOP_K=5 SIGMA0=0.05 MAX_GENERATIONS=18 bash scripts/campaigns/cmaes/run.sh
# HQ
SOURCE_RUN=.../ftl_4d_cmaes_v1 CANDIDATES="144 0" bash scripts/campaigns/hq/run_batch.sh
```

---

## Foundational entries (2026-06-10)

### Matter model — reference

Campaign evolves **N independent massive real scalar fields** ("lumps") via
`grtresna_independent_scalars` — 5 lumps, `recipe_scalar_mass` searched, per-lump exotic flags.

| Side | Key paths |
|------|-----------|
| GRTeclyn | `RadialRecipeMatterDispatch.hpp`, `GRTresnaIndependentScalars.{hpp,impl.hpp}`, `StateVariables.hpp` |
| GRTresna | `MatterParams.hpp` (lump_t), `MyMatterFunctions.cpp` |

Potential: `V = ½ m² (Σφ_k)²` (+ λφ⁴ when searched). Lumps interact via shared gravity + mass term;
O(1) boosts + light mass → fly-away (fixed v4 mass search + velocity caps).

**Roadmap (done in bold):** **(1)** search mass + cap boosts; **(2)** λφ⁴; **(3)** complex scalar /
boson star (Phase 1 — see [Boson star integration](#boson-star-matter-sector-integration--gpu-smoke-2026-06-18));
(4) per-lump independent mass; (5) multi-site boson stars at shell lump centres.

### Navigation overhaul (2026-06-10)

1. `speed_horizon` → **`speed_super`** descriptor.
2. ~82% pre-GPU rejection → tightened bounds, feasible-box sampling, harder GRTresna solve.
   `ftl_discovery_nav`: GPU-reach ~40%.

**Scoring fix (after ~90 evals):** stationary warp-lens artifacts — reliability-gate geodesic;
zero shaping when stationary + no dynamical FTL. Eval 83: 1164→−247.

### MAP-Elites FTL Discovery Status

Status: **reset**. `theta_plus` measured from origin not `grid_center` → false horizon penalty.
Fixed in `RadialRecipeLevel.cpp`; `ftl_discovery_postfix` confirmed `horizon_penalty=0`.
~93% pre-GPU rejection exposed navigation defects fixed above.
