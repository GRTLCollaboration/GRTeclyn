# Geometry-First Matter Projection

## Overview

The geometry-first approach inverts the usual numerical relativity workflow.
Instead of specifying matter and letting the solver determine the spacetime,
we **define a target geometry**, algebraically reconstruct the matter
distribution required to support it, fit scalar-field lumps to that
distribution, and then solve the constraints with GRTresna to verify the
geometry is actually realised.

The pipeline lives in `scripts/search/project_geometry_motif.py` and the
supporting modules under `src/grteclyn_wrapper/projection/` and
`src/grteclyn_wrapper/grtresna/fit/`.

---

## Pipeline Stages

```
Geometry-first episode          Motif scout               Matter fitting
(overrides, constraints)  -->  extract_motif_from_episode  -->  fit_matter_from_motif
                                       |                            |
                                  motif.json              fitted_matter.json
                                       |                            |
                                       v                            v
                            +----------+----------+    +-----------+-----------+
                            |  rho_req support    |    |  Gaussian lump params  |
                            |  regions, FTL flags |    |  (amp, width, center,  |
                            |  exotic_needed, etc |    |   velocity, omega)     |
                            +---------------------+    +-----------------------+
                                                                  |
                                                                  v
                                                        GRTresna solve
                                                         (constraint solver)
                                                                  |
                                                                  v
                                                        .gridinit (initial data)
                                                                  |
                                          +-----------------------+-----------------------+
                                          |                                               |
                                    Mismatch check                              Preservation check
                                    (chi_l2, beta_l2,                      (compare_motif_preservation:
                                     convergence penalty)                   f_op, polarity, support)
                                          |                                               |
                                          +-----------------------+-----------------------+
                                                                  |
                                                          Pass / Reject / Promote
```

### 1. Motif scout (`extract_motif_from_episode`)

Reads a geometry-first episode (RadialRecipe overrides + constraint norms)
and extracts a `GeometryMotif`:

- **Conformal factor** χ(r) and shift β(r) profiles from the recipe basis
- **Required matter** ρ_req from the Hamiltonian constraint
  (`rho = ±(R + ⅔K²) / 16π`, sign depends on phantom flag)
- **Support regions**: contiguous radial segments where |ρ_req| exceeds a
  threshold, each tagged as exotic (negative density) or canonical
- **FTL metrics**: f_op, f_null, f_portal, f_throat, f_asymmetry
- **Momentum target**: inferred from β asymmetry
- **Exotic flag**: `exotic_needed = True` if any support region has negative
  ρ or the constraint file reports `min_rho_required < 0`

Output: `motif.json`

### 2. Matter fitting (`fit_matter_from_motif`)

Maps support regions to Gaussian scalar lumps:

- One lump per support region (up to `max_lumps`)
- **Ring splitting**: when fewer support regions than `max_lumps`, the
  dominant region is split into additional lumps arranged in a ring in the
  xy-plane, giving the solver more degrees of freedom
- Each lump carries: amplitude, width, center (x,y,z), velocity (vx,vy,vz),
  omega, mode, exotic flag
- Exotic lumps require `maximal_slicing = True` in the solver

Output: `fitted_matter.json`

### 3. GRTresna solve

Solves the Hamiltonian and momentum constraints for the fitted matter:

- Lichnerowicz conformal-decomposition solver
- Non-linear (Picard) iteration with adaptive exit/stall tolerances
- For exotic matter: `max_NL_iterations = 200`, `nl_stall_tolerance = 0.005`
  (defaults are too aggressive and exit at ~96% Ham)
- Output: `initial_data.gridinit` + convergence diagnostics

### 4. Mismatch check (`compute_mismatch`)

Compares the solved geometry against the motif target:

- **chi_l2**: L2 norm of χ_solved − χ_target (resampled to common grid)
- **beta_l2**: L2 norm of β_solved − β_target
- **exotic_penalty**: penalty proportional to total exotic matter amplitude
- **convergence_penalty**: tanh-saturated Ham residual (see Phase 1a below)

### 5. Preservation check (`compare_motif_preservation`)

Higher-level geometric fidelity check:

- **f_op retention**: does the operational FTL feature survive?
- **polarity retention**: does the sign structure of χ match?
- **support_localized**: is the matter concentrated in the target region?
- **shift_alignment** / **momentum_alignment**: does the shift match?
- **retention_score**: weighted composite (threshold: 0.5 to promote)

---

## Iterative Adjustment Loop (CMA-ES)

The one-shot pipeline (stages 1–5) often produces a geometry that doesn't
match the target well — the initial matter guess is just a rough fit. The
`--iterate N` flag wraps stages 3–5 in a CMA-ES optimisation loop that
adjusts lump parameters to minimise the mismatch.

### Optimiser architecture

```
                    Amplitude pre-conditioning (Phase 1b)
                              |
                              v
                    CMA-ES (cma.CMAEvolutionStrategy)
                       /          \
                      /            \
              ask() -> solutions    tell(fitnesses)
                    |                 ^
                    v                 |
              _evaluate_candidate ----+
              (parallel via ThreadPoolExecutor)
                    |
        +-----------+-----------+
        |                       |
   _solve_with_fallback     compute_mismatch
   (Phase 1c ladder)        (Phase 1a two-phase fitness)
```

### Phase 1 improvements (implemented)

#### 1a. Two-phase fitness shaping (`mismatch.py`)

**Problem**: The original fitness was `chi_l2 + beta_l2 + convergence_penalty`
with the convergence penalty being the raw Ham percentage. For exotic matter
where Ham ~96%, this term (96) dominated the geometry mismatch (~0.05) by
three orders of magnitude, so CMA-ES had no geometry signal to follow.

**Solution**: A two-phase fitness function:

- **Feasibility phase** (Ham > 5%): fitness is dominated by a
  tanh-saturated convergence penalty:
  `convergence_penalty = 50 * tanh(Ham% / 20)`
  This caps the penalty at ~50 and gives CMA-ES a smooth gradient toward
  feasibility without the geometry term being drowned out.
- **Geometry phase** (Ham ≤ 5%): fitness is dominated by the L2 mismatch
  terms (chi_l2, beta_l2) plus a small exotic penalty. The convergence
  penalty is now negligible.

Constants: `FEASIBILITY_THRESHOLD = 5.0%`, `FEASIBILITY_WEIGHT = 50.0`,
`CONV_TANH_SCALE = 20.0`.

#### 1b. Amplitude pre-conditioning (`iterate.py`)

**Problem**: CMA-ES wasted evaluations finding a feasible amplitude — the
initial amplitude was too aggressive, causing solver crashes or Ham ~96%.

**Solution**: Before CMA-ES, run a bisection on a global amplitude scale
factor (1.0 → 0.5 → 0.25 → ... → 0.0625). For each scale, run GRTresna and
check if Ham < 10%. Use the largest feasible amplitude as the CMA-ES
starting point.

- `PRECOND_HAM_THRESHOLD = 10.0%`
- `PRECOND_MAX_STEPS = 5` (halving each time)
- `PRECOND_MIN_SCALE = 0.0625` (don't go below 1/16)

If no feasible amplitude is found, CMA-ES starts from the original — the
two-phase fitness (1a) will still try to find feasibility.

#### 1c. Solver fallback ladder (`iterate.py`)

**Problem**: Some CMA-ES candidates crash GRTresna or produce very high
residuals, wasting the evaluation.

**Solution**: On crash or Ham > 50%, retry once with safer settings:

- Relaxation factor: `min(current, 0.4)` (more conservative)
- Max iterations: `max(current, 150)` (more time to converge)
- Amplitude: `× 0.7` (reduce source strength)

If the fallback succeeds, use its result. If it also fails, return the gate
fitness (penalty).

#### 1d. Tighter bounds from motif support (`iterate.py`)

**Problem**: Default bounds (±64 for lump centers) span the entire grid box.
CMA-ES wastes time exploring empty regions.

**Solution**: Compute per-lump bounds from motif support regions:

- Center bounds: support-region centroid ± `max(2×width, 8)` (capped at
  `0.4 × grtresna_L`)
- Amplitude, width, velocity bounds unchanged

When no support regions exist, falls back to default ±64.

#### Multi-lump ring splitting (`fit/motif.py`)

**Problem**: A single support region produces a single lump, limiting the
solver's degrees of freedom for strong-curvature targets.

**Solution**: When `len(support_regions) < max_lumps`, split the dominant
region into `max_lumps` sub-lumps arranged in a ring in the xy-plane:

- Ring radius: `max(dominant.width × 0.8, 2.0)`
- Sub-lump amplitude: `original_amp / max_lumps`
- Sub-lump width: `dominant.width × 0.7` (clamped to [MIN, MAX])

#### Solver tuning for exotic matter

**Problem**: The default GRTresna solver exits at `NL_stall_tolerance = 2%`,
which stops at ~96% Ham for exotic configurations — the solver thinks it's
stalled when it's actually slowly converging.

**Solution**: For the iteration loop, use `max_NL_iterations = 200` and
`nl_stall_tolerance = 0.005` (0.5%). This gives the solver enough time to
push Ham down past the stall floor.

---

## Results

### Test configuration

- **Grid**: n=32, L=64 (GRTresna), L=16 (FTL integration)
- **Optimizer**: 40 evals, popsize=8, 5 generations, 4 concurrent solves
- **Matter**: 5 Gaussian lumps (ring-split from 1 support region)

### Exotic matter (phantom, χ_coeff_0 = −0.2, β_coeff_0 = 0.35)

| Metric | Value |
|--------|-------|
| Best fitness | 46.13 |
| chi_l2 | 0.055 |
| beta_l2 | 0.086 |
| convergence_penalty | ~92 (tanh) |
| Ham residual | ~96% (solver stalled) |
| Failed solves | 4 / 40 |
| Preservation score | 0.089 |
| Converged | No |

**Diagnosis**: GRTresna's Lichnerowicz solver cannot converge for this
exotic configuration — Ham plateaus at ~96% regardless of amplitude. The
physics bottleneck is the matter ansatz: a single Gaussian lump (even
ring-split into 5) cannot support the target curvature within the
convergence basin. The two-phase fitness correctly identifies this
(infeasibility phase) but CMA-ES has no path to feasibility.

### Canonical matter (no phantom, χ_coeff_0 = −0.08, β_coeff_0 = 0.0)

| Metric | Value |
|--------|-------|
| Best fitness | **0.026** |
| chi_l2 | 0.023 |
| beta_l2 | 0.000 |
| convergence_penalty | 0.007 |
| Ham residual | **0.61%** |
| Failed solves | **0 / 40** |
| Preservation score | **0.483** |
| Converged | No (threshold 0.5) |
| Precondition | Original amplitude feasible |

**Diagnosis**: The optimizer works end-to-end. Pre-conditioning confirmed
feasibility immediately (Ham=0.61%). CMA-ES found the best geometry match
in generation 1 (fitness 0.026). The preservation score (0.483) is just
below the 0.5 promotion threshold — the 5-lump ring spreads matter wider
than the single-region target, causing `support_localized = False`.

### Comparison

| Metric | Exotic | Canonical | Improvement |
|--------|--------|-----------|-------------|
| Best fitness | 46.13 | 0.026 | 1800× |
| Ham residual | 96% | 0.61% | 157× |
| Failed solves | 4/40 | 0/40 | — |
| Preservation | 0.089 | 0.483 | 5.4× |

---

## How to Run

### One-shot projection (no iteration)

```bash
.venv/bin/python scripts/search/project_geometry_motif.py \
  /path/to/geometry_first_episode \
  --out-dir /path/to/output \
  --mode solve-only \
  --max-lumps 5 \
  --gridinit-n 32 \
  --grtresna-L 64.0 \
  --mpi-ranks 4 \
  --ftl-L 16.0
```

### Iterative projection (CMA-ES loop)

```bash
.venv/bin/python scripts/search/project_geometry_motif.py \
  /path/to/geometry_first_episode \
  --out-dir /path/to/output \
  --mode solve-only \
  --iterate 40 \
  --iterate-popsize 8 \
  --max-concurrent-grtresna 4 \
  --max-lumps 5 \
  --gridinit-n 32 \
  --grtresna-L 64.0 \
  --mpi-ranks 4 \
  --ftl-L 16.0
```

### Non-exotic (canonical) matter

Add `--no-phantom` to use positive-energy scalar fields. The episode's
`constraint_norms.dat` must have `min_rho_required > 0` (column 4) for
`exotic_needed` to be `False`.

```bash
.venv/bin/python scripts/search/project_geometry_motif.py \
  /path/to/episode \
  --out-dir /path/to/output \
  --mode solve-only \
  --no-phantom \
  --iterate 40 \
  --max-lumps 5 \
  --gridinit-n 32 \
  --grtresna-L 64.0 \
  --mpi-ranks 4 \
  --ftl-L 16.0
```

### Full solve-and-evolve

```bash
.venv/bin/python scripts/search/project_geometry_motif.py \
  /path/to/episode \
  --out-dir /path/to/output \
  --mode solve-and-evolve \
  --iterate 40 \
  --max-lumps 5 \
  --gridinit-n 32 \
  --grtresna-L 64.0 \
  --mpi-ranks 4 \
  --ftl-L 16.0 \
  --stop-time 100 \
  --plot-interval 10
```

### Reusing existing motif/fitted matter

```bash
.venv/bin/python scripts/search/project_geometry_motif.py \
  /path/to/episode \
  --out-dir /path/to/output \
  --mode solve-only \
  --motif-json /path/to/existing/motif.json \
  --fitted-matter-json /path/to/existing/fitted_matter.json \
  --iterate 40
```

### Key CLI flags

| Flag | Default | Description |
|------|---------|-------------|
| `--mode` | — | `fit-only`, `solve-only`, `solve-and-evolve` |
| `--max-lumps` | 3 | Max Gaussian lumps (ring-split if fewer support regions) |
| `--iterate N` | 0 (off) | Max CMA-ES evaluations (0 = one-shot) |
| `--iterate-popsize` | 8 | CMA-ES population size |
| `--iterate-sigma0` | 0.2 | CMA-ES initial step size |
| `--max-concurrent-grtresna` | 6 | Parallel solver calls during iteration |
| `--phantom` / `--no-phantom` | phantom | Controls sign of ρ_req reconstruction |
| `--promote-retention-min` | 0.5 | Preservation score threshold to promote candidate |
| `--gridinit-n` | 64 | Grid resolution per axis |
| `--grtresna-L` | 128.0 | GRTresna box half-size |
| `--ftl-L` | 16.0 | FTL metric integration radius |
| `--mpi-ranks` | 8 | MPI ranks for GRTresna |

### Output artifacts

```
out_dir/
├── motif.json                      # Extracted geometry motif
├── fitted_matter.json              # Fitted lump parameters
├── momentum_target.json            # Inferred momentum target
├── initial_data.gridinit           # Solved initial data
├── projection_report_preservation.json   # Preservation check result
├── iterate/                        # (only with --iterate)
│   ├── iteration_log.jsonl         # Per-eval fitness log
│   ├── iterate_summary.json        # Final summary
│   ├── best_fitted_matter.json     # Best matter config
│   ├── best/initial_data.gridinit  # Best solved initial data
│   ├── precondition/               # Amplitude bisection traces (cleaned up)
│   └── eval_XXXXXX/                # Per-eval GRTresna work dirs (pruned)
└── grtresna/                       # One-shot GRTresna work dir
```

---

## Physics Caveats

### Expressivity bound

GRTresna cannot pin arbitrary geometries with a limited set of Gaussian
lumps. The Lichnerowicz equation is a non-linear elliptic PDE, and the
Gaussian ansatz spans only a small subspace of possible matter
distributions. Strong-curvature targets (large |χ_coeff|) may be
unreachable regardless of optimisation.

### Matching initial data ≠ matching spacetime

Even a perfect t=0 match between the solved geometry and the target motif
does not guarantee that the spacetime *evolves* to maintain that geometry.
The constraint-satisfying initial data may diverge from the target under
Einstein evolution. The `solve-and-evolve` mode addresses this by running
a short GRTeclyn evolution and checking post-load constraints.

### Underdetermination

Many matter configurations produce the same constraint-satisfying
geometry. This can be leveraged: the exotic penalty in the fitness
function prefers less exotic matter among configurations that produce
similar geometry mismatch.

### Exotic matter convergence

The Lichnerowicz solver struggles with negative-energy (phantom) scalar
fields. The Hamiltonian residual plateaus at ~96% for aggressive exotic
configurations, even with 200 iterations and tight stall tolerance. This
is a fundamental limitation of the Picard iteration for this class of
matter, not a code bug.

---

## Next Steps

### Phase 2: Richer matter basis

The current Gaussian lump ansatz is too restrictive. A multi-sector
encoder would allow:

- **Shells**: hollow spherical matter distributions (thin-shell theorems)
- **Boson stars**: complex scalar field solitons (already supported by
  GRTresna's `recipe_matter_model`)
- **Q-balls**: non-topological solitons with time-harmonic ansatz
- **Multi-component**: mix of canonical + exotic sectors

This would expand the expressivity of the matter ansatz and allow the
optimizer to find configurations that the Gaussian basis cannot reach.

### Phase 3: Structural outer loop (MAP-Elites)

Instead of optimising lump *parameters* for a fixed *structure* (number
and type of lumps), use MAP-Elites to explore different matter structures:

- **Outer loop**: MAP-Elites over matter structure (n_lumps, lump types,
  shell vs. blob, exotic fraction)
- **Inner loop**: CMA-ES over continuous parameters (amplitudes, widths,
  centers) for each structure

This separates the discrete structural search from the continuous
parameter optimisation and maintains an archive of diverse solutions.

### Phase 4: Better geometry targets (implemented)

#### 4a. 2D slice mismatch

The mismatch is now computed on the full xz-plane slice (y = domain centre),
not just the 1D midline. This captures angular structure that the radial
profile misses — critical for ring-distributed lumps and non-spherical
matter configurations.

- **Target**: the spherically symmetric RecipeBasis is evaluated on
  r = sqrt(x² + z²), producing 2D chi and beta fields. Radial beta is
  decomposed into beta_x and beta_z components.
- **Solved**: extracted from the gridinit using the existing
  `build_xz_slice_from_gridinit` infrastructure, then bilinearly resampled
  onto the target grid via `scipy.RegularGridInterpolator`.
- **Weights**: `W_CHI_2D = 2.0`, `W_BETA_2D = 1.5` (higher than 1D weights
  since the 2D slice carries more information).

#### 4b. K_ij matching

The extrinsic curvature tensor is now included in the mismatch:

- **K (trace)**: extracted from the gridinit `K` component. Target is K=0
  (maximal slicing). Weight: `W_KIJ = 0.5`.
- **A_ij (traceless)**: extracted from `A11..A33` components, combined into
  a scalar proxy |A_ij|²/2. Target is A_ij=0 (conformally flat). Weight:
  `W_AIJ = 0.3`.

This provides a stronger constraint on the initial data beyond just the
spatial geometry (chi, beta).

#### 4c. Feasibility pre-check

Before running GRTresna, `feasibility_precheck()` estimates whether the
target geometry is likely to converge, based on the peak |ρ_req| from the
motif's support regions:

| Risk level | ρ_peak range | Action |
|------------|-------------|--------|
| safe | < 0.5 | Proceed normally |
| marginal | 0.5 – 2.0 | Proceed with warning (expect slower convergence) |
| hard | > 2.0 | Log warning (likely infeasible for Gaussian ansatz) |

The pre-check does not skip solves — even "hard" targets are attempted
because the actual convergence depends on the matter ansatz and solver
settings. But it provides early diagnostics and could be used to skip
obviously infeasible targets in batch runs.

#### Impact on test results

Adding 2D + K_ij terms improved the optimizer's behavior:

| Metric | 1D only (v2) | 2D + K_ij (v3) |
|--------|-------------|----------------|
| Best fitness | 0.026 | 0.047 |
| Preservation | 0.483 | **0.490** |
| Improvement over generations | plateaued gen 1 | improved through gen 4 |
| chi_1d | 0.023 | 0.016 |
| chi_2d | — | 0.006 |
| K_ij | — | 0.004 |
| A_ij | — | 0.007 |

The fitness is higher in absolute terms (more terms), but the optimizer
keeps improving longer and the preservation score is closer to the 0.5
threshold. The 2D chi mismatch (0.006) is smaller than 1D (0.016) because
the ring-distributed lumps produce a more uniform 2D field.

### Phase 5: Close the loop with evolution

Promote converged candidates (preservation score ≥ threshold) to full
GPU evolutions in GRTeclyn:

- Run a short evolution (t = 0–50M) and check if the geometry persists
- Compare the evolved geometry against the motif target at multiple
  time slices
- Use the evolution mismatch as a secondary fitness signal for the
  optimizer

---

## MAP-Elites Campaign (Phase 3 — implemented)

The CMA-ES loop above uses a coarse grid (N=32, RANKS=4) and a Gaussian
lump ansatz. The **MAP-Elites campaign** replaces it with the same
infrastructure as the QD FTL search: shell ansatz, N=128 grid, RANKS=8,
and a quality-diversity archive that explores the geometry-mismatch
landscape.

### How to run

```bash
# 1. Generate a motif.json from a geometry-first episode
.venv/bin/python scripts/search/project_geometry_motif.py \
  /path/to/episode --out-dir /tmp/motif_gen --mode fit-only

# 2. Run the MAP-Elites campaign
MOTIF_JSON=/tmp/motif_gen/motif.json \
  GPU_IDS="0 1 2 3 4 5 6 7" \
  QD_TARGET_EVALS=200 \
  bash scripts/campaigns/geometry_first/run.sh
```

Key environment variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `MOTIF_JSON` | required | Path to target motif.json |
| `GPU_IDS` | 0–7 | GPU list for parallel GRTresna solves |
| `QD_TARGET_EVALS` | — | Total evals (e.g. 200) |
| `RANKS` | 8 | MPI ranks per GRTresna solve |
| `BINS` | 8 | MAP-Elites grid (8×8 = 64 cells) |
| `SHELL_PROFILE` | compact | Shell ansatz profile preset |
| `LUMPS` | 5 | Number of lumps in the shell distribution |

### CMA-ES vs MAP-Elites

| | CMA-ES (`--iterate`) | MAP-Elites (campaign) |
|---|---|---|
| **Search** | Gradient-free, single trajectory | Quality-diversity archive (64 cells) |
| **Ansatz** | Gaussian lumps (45D) | Shell — Fibonacci sphere (23D) |
| **Grid** | N=32×32×16, L=64 | N=64+AMR L3, gridinit 128³, L=128 |
| **RANKS** | 4 | 8 |
| **Cell size** | 2.0 | 0.5 (4× finer) |
| **Solve time** | ~15s | ~30–60s (higher resolution) |
| **Diversity** | Single best solution | Archive of diverse solutions |
| **Descriptors** | None | chi match × convergence quality |
| **Score** | Minimise mismatch | Maximise (100 − mismatch) |
| **Refinement** | Built-in | Follow with `campaigns/cmaes/run.sh` |

### Is it better?

**Yes, in three ways:**

1. **Resolution**: N=128 vs N=32 means the solver can resolve curvature
   features 4× smaller. The CMA-ES runs couldn't pass `support_localized`
   because the coarse grid smoothed χ to ~0.06 deviation — at N=128 the
   same target produces sharper features.

2. **Ansatz**: The shell ansatz distributes lumps over a full sphere
   (Fibonacci lattice) with toroidal/poloidal currents, dipole/quadrupole
   modulation, and exotic-sector wedges. This is a richer matter basis
   than the Gaussian ring, giving the solver more expressive power.

3. **Diversity**: MAP-Elites maintains an archive of solutions across
   the chi-match × convergence-quality plane. Instead of a single
   optimum, you get a family of candidates — some with perfect chi
   match, some with better convergence, some with different matter
   layouts. CMA-ES converges to one point and can't recover from a
   bad basin.

**Trade-off**: each eval is ~2–4× slower (higher resolution, more MPI
ranks). A 200-eval campaign takes ~1.5–3 hours vs ~4 minutes for the
CMA-ES loop. But the quality and diversity of results justify the cost.

**Test result** (4 evals, RANKS=4): best score 99.91 (mismatch=0.089),
chi_l2=0.040, chi_2d_l2=0.015, Ham=0.75%. Two elites in two cells —
diversity working from the first batch.

### Next steps

1. **Run a full 200-eval campaign** on a strong target motif to fill the
   archive and identify the best cell.
2. **CMA-ES refinement**: seed `campaigns/cmaes/run.sh` with the best
   MAP-Elites elite to fine-tune continuous parameters.
3. **Evolution check**: promote the best candidate to a GRTeclyn
   evolution (HQ replay) to verify the geometry persists in time.
4. **Exotic targets**: test with phantom matter motifs once the canonical
   pipeline produces passing preservation scores.
5. **Fix unrelated QD-resume test** (`tests/search/test_proposed_extensions.py::
   test_qd_search_resume_continues_eval_counter`) — pre-existing, archive
   records 0 elites on resume; not a geometry-first regression.

---

## Warp-Drive (Alcubierre) Matter Reconstruction — Implemented

### Can we reconstruct the matter required for an Alcubierre warp drive?

**Yes — the exact stress-energy is reconstructed and cross-checked.** The
pipeline now supports Alcubierre warp-bubble targets end-to-end: motif
generation, toroidal exotic matter fitting, warp-target mismatch, and a
warp-factory cross-check that validates the reconstructed `rho` against the
exact `T = G/8 pi` from the analytic 4-metric.

The honest technical picture:

- The infrastructure to *specify and analyse* an Alcubierre metric already
  existed: `metrics/probes/warpfactory.py` builds the analytic 4-metric
  (`alcubierre_metric`, `_alcubierre_shape`) and computes the exact
  stress-energy `T_{mu nu} = G_{mu nu} / 8 pi` plus energy-condition
  diagnostics. `initial_data/seeds.py::alcubierre_warp` already emits a
  radial `beta^x(r) = -v f(r)` recipe seed.
- The geometry-first *reconstruction* pipeline was built around
  **spherically-symmetric conformal geometry** (`RecipeBasis` chi(r),
  beta(r)) with a **maximal-slicing / conformally-flat target**
  (K -> 0, A_ij -> 0). Alcubierre is the opposite regime: on a t=const
  slice the spatial metric is **flat (chi = 1)** and *all* the physics —
  including the negative-energy toroidal ring — lives in the **extrinsic
  curvature K_ij**, sourced by an **axial shift** beta^x = -v f(r_s).
- The old `rho_req` reconstruction (from the conformal Hamiltonian
  constraint) returns ~0 for a flat-chi Alcubierre slice: the exotic energy
  is hidden in K_ij, not in chi. The new `motif_from_alcubierre` derives
  `rho_req` from the analytic Eulerian energy density instead.

### Physics target

For a bubble travelling along x at speed v with shape function f(r_s):

- Spatial metric: gamma_ij = delta_ij (flat), lapse alpha = 1.
- Shift: beta^i = (-v f(r_s), 0, 0), r_s = sqrt((x - x_s)^2 + y^2 + z^2).
- Eulerian energy density (Alcubierre 1994):
  rho = -(v^2 / 32 pi) * (rho_cyl^2 / r_s^2) * (df/dr_s)^2 <= 0,
  a **toroidal ring** of negative energy in the plane perpendicular to
  travel, peaked at the bubble wall (r_s ~ bubble_radius).
- Extrinsic curvature K_ij is non-zero (A_ij != 0 from the shift gradient),
  K (trace) = 0.

### Implementation (5 steps, all complete)

1. **Alcubierre motif source** (`initial_data/motif.py`):
   `motif_from_alcubierre(velocity, bubble_radius, sigma, ...)` produces a
   `GeometryMotif` with flat chi, an axial shift (reusing the
   `alcubierre_warp` seed's beta fit), and `rho_req` derived from the
   analytic Eulerian energy density (not the conformal Hamiltonian path).
   Also exports `alcubierre_shape_function` and `_alcubierre_eulerian_rho`.

2. **Toroidal support regions** (`initial_data/motif.py`):
   `_toroidal_support_regions` extracts a ring-shaped exotic support region
   at cylindrical radius ~ bubble_radius in the plane perpendicular to the
   transport axis, tagged `exotic=True`.

3. **Toroidal exotic ansatz** (`grtresna/fit/motif.py`):
   `fit_matter_from_motif` now detects the `toroidal` momentum template and
   distributes lumps in a **yz-plane ring** (perpendicular to the x-transport
   axis) centred on the transport axis, at cylindrical radius = bubble
   radius. All lumps are exotic with axial velocity and toroidal omega.

4. **Warp-target mismatch** (`projection/mismatch.py`):
   `_target_2d_slice_warp` evaluates flat chi (1.0), axial beta_x = -v f(r_s),
   K=0, and **non-zero A_ij** from the shift gradient. `compute_mismatch`
   detects warp motifs via `_is_warp_motif` and uses these targets instead
   of the spherical-recipe defaults (which incorrectly pushed K_ij -> 0).

5. **CLI + warp-factory cross-check** (`scripts/search/project_geometry_motif.py`):
   `--target alcubierre` with `--warp-velocity`, `--warp-bubble-radius`,
   `--warp-sigma` generates the motif and runs `warp_factory_cross_check`,
   which compares the reconstructed `rho` against the exact `T = G/8 pi`
   from `warpfactory.alcubierre_metric` and reports NEC/WEC violation
   fractions and the exotic energy budget.

### How to run

```bash
# Generate motif + fitted matter + warp-factory cross-check (no solve)
.venv/bin/python scripts/search/project_geometry_motif.py \
  --target alcubierre \
  --warp-velocity 0.5 \
  --warp-bubble-radius 2.0 \
  --warp-sigma 2.0 \
  --out-dir /tmp/alcubierre_test \
  --mode fit-only \
  --ftl-L 8.0

# Full solve (will hit the exotic convergence wall — see feasibility bound)
.venv/bin/python scripts/search/project_geometry_motif.py \
  --target alcubierre \
  --warp-velocity 0.5 \
  --warp-bubble-radius 2.0 \
  --warp-sigma 2.0 \
  --out-dir /tmp/alcubierre_solve \
  --mode solve-only \
  --max-lumps 5 \
  --gridinit-n 32 \
  --grtresna-L 64.0 \
  --mpi-ranks 4 \
  --ftl-L 8.0
```

### Results (v=0.5, R=2.0, sigma=2.0)

#### Reconstruction (analytic, cross-checked)

| Metric | Value |
|--------|-------|
| rho_l2 (motif vs T=G/8pi) | **0.000027** |
| rho_min (motif) | -0.00213 |
| rho_min (warp-factory) | -0.00224 |
| NEC violation fraction | 33.9% |
| WEC violation fraction | 40.8% |
| Exotic energy | 0.0565 (geometric units) |
| Support regions | 1 (toroidal, exotic) |
| Fitted lumps | 5 (toroidal ring, exotic, axial v=0.6) |

The reconstructed `rho` matches the exact `T = G/8 pi` to within finite-
difference noise (L2 = 2.7e-5), confirming the analytic reconstruction is
correct. The NEC is violated at 33.9% of grid points — exotic matter is
required, as expected for an Alcubierre warp bubble.

#### GRTresna solve (measured, n=32, L=64, 4 MPI ranks, Q-ball lumps)

After upgrading the toroidal lumps from plain Gaussian wave packets
(scalar_mass=0, no binding) to **Q-ball solitons** (compact preset: m=1,
λ=640, μ=85333, ω=0.8, ODE profile, equilibrium amplitude) — the same
couplings the QD MAP-Elites campaign uses — and switching to a full-box
domain (lo_boundary=(0,0,0), N=(32,32,32), lumps centred at (L/2,L/2,L/2)):

| Metric | Plain Gaussian | Q-ball soliton | QD champion (eval 118) |
|--------|---------------:|---------------:|-----------------------:|
| Ham residual | 95.83% | **13.5%** (still dropping at iter 43) | 1.1% (evolved) |
| Mom residual | 0.01% | **66.86%** (stuck) | 0.1% (evolved) |
| shift_max | 0.0 | **0.0** | **0.0** |
| chi range | 0.997–1.000 | 0.996–1.010 | 0.954–1.335 |
| scalar phi max | 0.014 | 0.018 | 0.213 |

**Root cause — GRTresna does not solve for the shift from scalar matter.**

The critical finding: GRTresna's constraint solver is a **Lichnerowicz
solver for the conformal factor (chi/psi) only**. It does NOT solve the
momentum constraint for the shift vector (Beta^i) from scalar field
sources. The shift is always **exactly zero** — confirmed in both the
Alcubierre solve AND the QD MAP-Elites champion (eval 118: beta_mean=0.0,
stationary=True, shift1/2/3 all exactly 0.0).

The QD campaign achieves FTL-like effects through **chi curvature alone**
(conformal geometry, not shift) — the exotic Q-ball lumps curve the
conformal factor, creating a gauge-invariant null-geodesic shortcut without
any shift. This is a fundamentally different mechanism from the Alcubierre
warp drive, which requires `beta^x = -v·f(r_s)` — the shift IS the warp
bubble.

**Why the momentum constraint is stuck at 66.86%**: The `use_compact_Vi_ansatz`
flag enables a momentum constraint solve, but the ansatz appears designed
for BH binary initial data (where point-source momenta drive the shift),
not for distributed scalar field sources. The scalar field's momentum
density (S_i from Pi and grad_phi) cannot be represented by the compact Vi
ansatz, so the momentum constraint residual plateaus.

**Implication**: The current GRTresna solver **cannot reconstruct an
Alcubierre warp bubble**, regardless of the matter configuration. This is
not an exotic matter convergence problem (Ham converges fine with Q-ball
lumps). It is a fundamental limitation: the solver does not produce shifts
from scalar field matter. Reconstructing an Alcubierre warp bubble would
require either:
  1. A momentum constraint solver that can solve for Beta^i from scalar
     field S_i (not just the compact Vi ansatz), or
  2. Painting the Alcubierre shift directly as a recipe seed (bypassing the
     momentum constraint solve), then solving only the Hamiltonian
     constraint for chi given the fixed shift + matter.

### Feasibility bound (confirmed by measured solve)

- **Exact reconstruction** of the required stress-energy is analytic and
  reliable (`T = G/8 pi` from the warp-factory metric, cross-checked to
  L2 = 2.7e-5): this quantifies exactly how much negative energy the bubble
  needs (exotic energy = 0.0565 geometric units for v=0.5, R=2.0).
- **Q-ball lumps converge the Hamiltonian constraint** (Ham → 13.5% and
  still dropping at iteration 43). The exotic matter is NOT the problem —
  the same Q-ball couplings the QD campaign uses work fine for Ham.
- **GRTresna does not solve for the shift from scalar matter.** The shift
  is exactly 0.0 in both the Alcubierre solve and the QD champion (eval
  118). The momentum constraint is stuck at 66.86% because the compact Vi
  ansatz cannot represent distributed scalar field momentum sources.
- **The Alcubierre warp drive requires a non-zero shift** (beta^x = -v·f).
  Without it, there is no warp bubble. The QD campaign's FTL effects come
  from chi curvature, a fundamentally different mechanism.
- **Root cause in GRTresna C++**: `set_output_data` in
  `GRTresna/Source/Tools/WriteOutput.H` never writes `c_shift1/2/3` — they
  stay at 0.0 from the initial `setVal(0.0)` loop. The CTTK method
  (`Source/Methods/CTTK.impl.hpp`) *does* solve the momentum constraint for
  V_i (with source `8πG·S_i`), but the output stage discards V_i and never
  converts it to a shift. This is why the shift is exactly 0.0 in every
  GRTresna solve, including the QD champion.

### Next steps (in order of increasing effort)

#### Step 1: Post-process — re-paint the Alcubierre shift onto the solved gridinit

**Effort: Small (Python only, no C++ changes)**

After GRTresna solves for chi and A_ij (with shift=0), paint the Alcubierre
shift `beta^x = -v·f(r_s)` back onto the gridinit in Python.  This gives:
- chi consistent with the exotic Q-ball matter (from the solve)
- The prescribed Alcubierre shift (painted post-solve)

The Hamiltonian constraint won't be exactly satisfied (chi was solved
without the shift's contribution to the extrinsic curvature), but it should
be close, and the GRTeclyn evolution will adjust.  This tests whether the
exotic toroidal ring + prescribed shift produces a stable warp bubble.

Implementation:
1. In `project_geometry_motif.py`, after the GRTresna solve, read the
   gridinit and overwrite `shift1` with `-v * f(r_s)` centred at the domain
   centre, using `alcubierre_shape_function`.
2. Re-run the preservation check with the painted shift.
3. Evolve with GRTeclyn to see if the warp bubble survives.

#### Step 2: Fix GRTresna's `set_output_data` to convert V_i to shift

**Effort: Medium (one C++ function change in `WriteOutput.H`)**

The CTTK method already solves the momentum constraint for V_i.  The shift
is `Beta^i = (2/7) * psi^6 * V^i` (B&S eq. 3.38, depending on convention).
Add this conversion to `set_output_data` so the solved shift is written to
the gridinit.  This would let GRTresna produce a shift from the scalar
field's momentum density — the "correct" fix.

Implementation:
1. In `GRTresna/Source/Tools/WriteOutput.H`, `set_output_data`, add:
   ```cpp
   // Convert momentum constraint solution V_i to shift Beta^i
   Real psi = multigrid_vars_box(iv, c_psi_reg) + psi_bh;
   Real psi6 = pow(psi, 6.0);
   grchombo_vars_box(iv, c_shift1) = (2.0/7.0) * psi6 * multigrid_vars_box(iv, c_V1_0);
   grchombo_vars_box(iv, c_shift2) = (2.0/7.0) * psi6 * multigrid_vars_box(iv, c_V2_0);
   grchombo_vars_box(iv, c_shift3) = (2.0/7.0) * psi6 * multigrid_vars_box(iv, c_V3_0);
   ```
2. Re-run the Alcubierre solve and check if the shift is non-zero.
3. If the momentum constraint is still stuck at ~67%, the compact Vi ansatz
   may need to be replaced with the non-compact (periodic) ansatz
   (`use_compact_Vi_ansatz = 0`).

#### Step 3: Two-step solve (shift-prescribed Hamiltonian solve)

**Effort: Medium (Python pipeline + possibly C++ params)**

1. Paint the Alcubierre shift via the `alcubierre_warp` recipe seed (which
   already fits `beta^x = -v·f(r)` as `recipe_beta_coeff_*` overrides).
2. Run GRTresna with the shift fixed as a background, solving only for chi
   given the shift + exotic matter.  This requires telling GRTresna to keep
   the shift fixed — either a new param (`fix_shift = 1`) or a CTTK
   modification that skips the momentum constraint solve and uses the
   painted shift for the A_ij source.
3. The Hamiltonian constraint with a fixed shift sources A_ij through the
   extrinsic curvature, so chi will adjust to support both the matter and
   the shift.  This is the physically correct reconstruction.

#### Step 4: Full momentum constraint solve from scalar S_i

**Effort: Large (C++ solver work)**

The CTTK RHS already has `8πG·S_i` as the source for V_i.  Even with the
Step 2 output fix, the momentum constraint may not converge for distributed
scalar sources (Mom stuck at 66.86% in the current run).  This would
require:
1. Investigating why the momentum constraint stalls — is it the compact Vi
   ansatz, the multigrid solver tolerance, or the scalar field S_i being
   too weak/delocalised?
2. Potentially switching to `use_compact_Vi_ansatz = 0` (non-compact /
   periodic ansatz, B&S eq. B.3) for distributed sources.
3. Increasing `max_NL_iterations` and relaxing the stall tolerance for the
   momentum solve.

### Recommended order

Start with **Step 1** (post-process shift painting) — it's the fastest way
to get a testable gridinit with both the exotic matter and the Alcubierre
shift.  If the GRTeclyn evolution shows the warp bubble is stable, that's
the deliverable.  If not, proceed to **Step 2** (fix the V_i → shift
conversion in GRTresna's output) to let the solver produce the shift
self-consistently from the scalar field momentum.

### Verification

- 30 new tests in `tests/projection/test_alcubierre.py`: shape function,
  motif generation, toroidal fitting, warp-target mismatch, and warp-
  factory cross-check. All pass.
- `pytest tests/projection -q`: 79 passed (49 original + 30 new).
- `pytest tests/search -q`: 99 passed, 1 pre-existing unrelated failure
  (`test_qd_search_resume_continues_eval_counter`).

---

## File Reference

| File | Purpose |
|------|---------|
| `scripts/search/project_geometry_motif.py` | CLI entry point (CMA-ES loop) |
| `scripts/campaigns/geometry_first/run.sh` | MAP-Elites campaign launcher |
| `src/grteclyn_wrapper/initial_data/motif.py` | Motif extraction from episodes |
| `src/grteclyn_wrapper/grtresna/fit/motif.py` | Matter fitting (lumps, ring splitting) |
| `src/grteclyn_wrapper/projection/iterate.py` | CMA-ES iteration loop |
| `src/grteclyn_wrapper/projection/mismatch.py` | Two-phase fitness + 2D/K_ij mismatch + feasibility pre-check |
| `src/grteclyn_wrapper/projection/motif_preservation.py` | Preservation check |
| `src/grteclyn_wrapper/search/qd_search/driver.py` | MAP-Elites driver (`geometry_first` mode) |
| `src/grteclyn_wrapper/search/qd_search/descriptors.py` | `geometry_first` behavior descriptors |
| `tests/projection/test_iterate.py` | Tests (37 passing) |
| `tests/projection/test_alcubierre.py` | Alcubierre warp-drive tests (30 passing) |
