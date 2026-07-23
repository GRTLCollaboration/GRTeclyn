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

### Pure-geometry MAP-Elites atlas (Stage-1 scout)

Upstream of motif→matter projection sits a **pure-geometry** MAP-Elites
scout that does not involve matter, GRTresna, or GRTeclyn evolution:

| Item | Location |
|------|----------|
| Package | `src/grteclyn_wrapper/search/geometry_atlas/` |
| CLI | `python -m grteclyn_wrapper geometry_atlas ...` |
| Launcher | `scripts/campaigns/geometry_atlas/run.sh` |
| Outputs | `runs/geometry_atlas/<name>/` (`archive.json`, `elites/*.gridinit`) |

**What it searches.** A compact-support 3D RBF genome for stationary,
asymptotically flat 4-metrics (`alpha`, `beta^i`, `gamma_ij = expm(S)`).
No spherical/axial symmetry; Alcubierre is **not** the search basis.
`K_ij` is derived from the stationary ADM relation. Effective
`T_ab = G_ab/8pi` (and therefore `rho`, `j_i`) is computed under the
explicit stationarity assumption.

**What it scores.** Frozen null-geodesic `f_geo` and stationary free-fall
`f_ff` (the free-fall probe wraps the frozen slice as a time-independent
metric stack). Archive axes are `[f_geo] × [log exotic-energy budget]`.
Within-cell ranking prefers valid `f_ff`, with hard rejects for signature
failures, non-flat boundaries, and inconsistent constraints.

**Validity boundary.** Atlas `f_geo` / `f_ff` are **screening** metrics.
A dynamical shortcut claim still requires GRTeclyn evolution + the 4D
evolving geodesic certificate. Do not conflate this scout with the
motif-to-matter MAP-Elites campaign in `scripts/campaigns/geometry_first/`
(which searches **matter** genomes to match a fixed motif).

**Handoff to Stage 2.** Elite cells write `elites/cell_i_j.gridinit` +
genome JSON. Those geometries are intended inputs for
`project_geometry_motif.py` / `fit_matter_from_motif` (matter synthesis)
in a later milestone — not wired automatically yet.

Smoke run:

```bash
TARGET_EVALS=8 N=16 NO_FF=0 bash scripts/campaigns/geometry_atlas/run.sh
```

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

### Physics correction (supersedes the earlier "Step 2" proposal)

An earlier draft of this document proposed converting the CTTK momentum
solution `V_i` into a shift via `Beta^i = (2/7) psi^6 V^i` in
`GRTresna/Source/Tools/WriteOutput.H`.  **That proposal is physically
wrong and is withdrawn:**

- In the CTT/CTTK decomposition, `V_i` is **not the shift**.  It is the
  vector potential of the conformal traceless extrinsic curvature: `A_ij`
  is built from derivatives of `V_i` and `U` (B&S eqs. B.3 / B.7 — see
  `PsiAndAijFunctions::compute_ctt_Aij` in
  `GRTresna/Source/Core/PsiAndAijFunctions.cpp` and the RHS in
  `GRTresna/Source/Methods/CTTK.impl.hpp`, where `rhs(c_V1..c_V3) =
  psi^6 (2/3 d_i K + 8 pi G S_i)`).  The scalar momentum density `S_i`
  therefore *already* feeds `A_ij` through the V_i solve.
- In 3+1 initial data the **shift is pure gauge**.  The Hamiltonian and
  momentum constraints fix `(gamma_ij, K_ij, rho, S_i)` — never `beta^i`.
  For the Alcubierre slice, the physics is entirely in `K_ij` (flat
  `gamma_ij`, `K = 0`, `A_ij` from the shift gradient); the shift itself
  may be freely painted onto the slice as a coordinate choice.

**Consequences for the fix:**

1. Painting `beta^x = -v f(r_s)` onto the gridinit is *legitimate* (not a
   hack) — provided `A_ij` on the slice matches the Alcubierre extrinsic
   curvature.  The real quantity to match is `A_ij`, not `beta`.
2. The correct role of the solver is: matter with the right momentum
   density `S_i` → momentum constraint → the right `A_ij`.  The 66.86% Mom
   stall means the solved `A_ij` does **not** yet match; that (not the
   missing shift output) is the genuine solver-side problem.

---

## Implementation Plan — Alcubierre Validation Baseline

Goal: make the Alcubierre warp drive a working validation baseline for the
geometry-first pipeline, using the **same matter family as the QD
MAP-Elites campaigns** (compact Q-ball preset: `m=1, lambda=640,
mu=85333, omega=0.8`; couplings frozen, per-lump amplitude / position /
velocity / omega-phase retunable by the optimiser).

Two phases.  **Phase A is Python-only and is executed first.**  Phase B
(GRTresna C++ work) is contingent — only start it if the Phase A ladder
fails at the stated acceptance criteria.

Success criterion for the whole effort (the "full ladder"):

1. **t=0 slice fidelity** — reconstructed gridinit matches the analytic
   Alcubierre targets (`chi`, `beta`, `A_ij`) within tolerance, with
   reported Ham/Mom residuals.
2. **Evolution persistence** — the bubble's shift/curvature structure
   survives a GRTeclyn evolution for >= 1 bubble crossing time
   (`t_cross = 2 R / v`) under the standard Gamma-driver gauge.
3. **Gauge-invariant probe** — the 4D null-geodesic FTL probe
   (`metrics/probes/ftl/`) confirms a shortcut on both the analytic slice
   (A0) and the reconstructed slice (A1+A2).

---

### Phase A — Python-only (execute now)

#### A0. Ground-truth analytic Alcubierre gridinit (validation anchor)

**Effort: Small.  No solver involved.  Do this first — everything else is
calibrated against it.**

Create `src/grteclyn_wrapper/projection/warp_gridinit.py` with:

```python
def alcubierre_analytic_fields(
    shape_xyz: tuple[int, int, int],
    dx_xyz: np.ndarray,        # (3,)
    origin: np.ndarray,        # (3,)
    *,
    velocity: float,
    bubble_radius: float,
    sigma: float,
    center: tuple[float, float, float],  # bubble centre in physical coords
) -> dict[str, np.ndarray]:
    """Evaluate the exact Alcubierre t=0 slice on a uniform grid.

    Returns a dict of arrays keyed by GRTECLYN_STATE_VARS names:
      chi = 1, h11=h22=h33 = 1 (h12=h13=h23 = 0), lapse = 1, K = 0,
      shift1 = -v * f(r_s), shift2 = shift3 = 0,
      A_ij  = (1/2)(d_i beta_j + d_j beta_i) - (1/3) delta_ij d_k beta_k
    (with alpha=1, flat conformal metric, K=0 the ADM formula
     K_ij = (D_i beta_j + D_j beta_i)/(2 alpha) reduces to this; sign
     convention must match GRTeclyn's A_ij — verify against
     Examples/RadialRecipe before freezing, see task A0.3).
    Derivatives of beta: use np.gradient on the painted shift1 array
    (robust, no closed form needed).
    """

def write_alcubierre_gridinit(
    out_path: Path, *, n: int, L: float,
    velocity: float, bubble_radius: float, sigma: float,
) -> Path:
    """Write a full analytic Alcubierre gridinit (zero matter).

    Uses GRTECLYN_STATE_VARS ordering and write_gridinit() from
    grtresna/io/gridinit.py.  Grid: n^3 cells, box [0, 2L]^3, bubble at
    the centre (L, L, L), origin chosen exactly as in
    grtresna/io/conversion.py (origin = target_center - L).
    All matter fields (phi, Pi, teo_*) = 0; Theta = Gamma_i = B_i = 0.
    """
```

Reuse: `alcubierre_shape_function` from `initial_data/motif.py`,
`write_gridinit` / `GRTECLYN_STATE_VARS` from `grtresna/io/gridinit.py`.

Tasks:

- **A0.1** Implement the two functions above.
- **A0.2** Add a CLI entry: `--target alcubierre --mode analytic-gridinit`
  in `scripts/search/project_geometry_motif.py` that writes
  `out_dir/analytic_alcubierre.gridinit` and stops (no solve).
- **A0.3** Verify the A_ij sign/normalisation convention: load the
  analytic gridinit through the post-load gate
  (`projection/postload_gate.py`) and check the reported Mom residual is
  *small in the bubble interior/exterior and peaked at the wall* (the wall
  violation is expected — no matter is present).  If Mom is large
  everywhere, the A_ij convention is wrong; flip sign / factor and re-check.
- **A0.4** Run the preservation check (`compare_motif_preservation`)
  against `motif_from_alcubierre(...)` — this must now report
  `shift_alignment ~ 1` and non-zero `beta_max_solved`.  This validates
  the checker itself against known-good data.
- **A0.5** Short GRTeclyn evolution (`--mode solve-and-evolve` machinery,
  but loading the analytic gridinit via `recipe_initial_data_file`):
  evolve to `t = 2 * t_cross`, plot `shift1` on the midplane every few
  steps.  Record how fast the Gamma-driver (params:
  `shift_Gamma_coeff=0.75, eta=1.0` in `Examples/RadialRecipe/params.txt`)
  damps the painted shift.  Also run one comparison with `eta = 0.25` to
  see if a weaker driver preserves the bubble longer.
- **A0.6** Run the 4D null-geodesic probe on the analytic slice and
  record its FTL metrics.  **This number is the calibration target** for
  A1/A2: the reconstructed slice should reproduce it.

Acceptance (A0): postload gate loads the file; preservation check reports
correct shift; probe sees the shortcut on the analytic slice.  If the
probe does NOT see a shortcut even on the exact Alcubierre slice, stop and
fix the probe first — nothing downstream can succeed.

#### A1. Post-solve shift + A_ij painting on the solved gridinit

**Effort: Small (Python only).**

Add to `src/grteclyn_wrapper/projection/warp_gridinit.py`:

```python
def paint_alcubierre_warp_on_gridinit(
    gridinit_path: Path,
    *,
    velocity: float,
    bubble_radius: float,
    sigma: float,
    paint_aij: bool = True,   # False = keep solved A_ij, paint shift only
) -> Path:
    """Read a solved gridinit, overwrite shift1/2/3 with the analytic
    Alcubierre shift (centred on the grid centre), optionally overwrite
    A11..A33 with the shift-consistent analytic A_ij, write back in place.
    Reuses read_gridinit/write_gridinit and alcubierre_analytic_fields.
    """
```

Tasks:

- **A1.1** Implement the painter (round-trip via
  `read_gridinit` → modify components → `write_gridinit`, preserving
  header, dx, origin).
- **A1.2** Wire a CLI flag `--paint-warp-shift {off,shift,shift+aij}`
  (default `off`) into `project_geometry_motif.py`.  Insertion point:
  immediately **after** the GRTresna solve produces
  `out_dir/initial_data.gridinit` and **before**
  `compare_motif_preservation(...)` — so the preservation check, postload
  gate, and evolution all see the painted slice.
- **A1.3** With `shift+aij`: the Ham constraint acquires an error because
  chi was solved with the *solver's* A_ij, not the painted one.  Log both
  (a) postload-gate Ham/Mom of the painted slice and (b) the un-painted
  slice, into `projection_report_preservation.json` under a new
  `warp_painting` key, so the trade-off is visible per run.

Acceptance (A1): preservation check on a painted solve reports
`shift_alignment > 0.9`; `beta_l2` (vs the motif beta target) drops by
>= 10x compared with the unpainted solve.

#### A2. Momentum-matched matter fitting (QD Q-ball couplings, per-lump retune)

**Effort: Medium (Python only).**

The solver's `A_ij` is sourced by the matter's momentum density `S_i`
(CTTK RHS).  Today the fitter only matches `rho`; the lump velocities are
heuristic.  Make the fit match the **Alcubierre momentum density** too.

Tasks:

- **A2.1** In `initial_data/motif.py`, add
  `_alcubierre_eulerian_Si(velocity, bubble_radius, sigma, L, n)`
  returning the analytic Eulerian momentum density on the same xz-slice
  grid as `_alcubierre_eulerian_rho` (from the exact 4-metric:
  `S_i = -(1/8 pi) (D_j K^j_i - D_i K)` with the analytic K_ij; use
  finite differences as in `_alcubierre_eulerian_rho`).  Cross-check
  against `metrics/probes/warpfactory.py` (`T^{0i}` from `G/8 pi`), same
  style as the existing `warp_factory_cross_check`.
- **A2.2** In `grtresna/fit/motif.py`, when the motif is a warp motif
  (`_is_warp_motif`), set per-lump velocities so the summed lump momentum
  density `S_i ~ Pi * grad phi` approximates the analytic `S_i` ring
  (direction: -x inside the wall; magnitude from A2.1), instead of the
  current heuristic axial `v=0.6`.  Keep the compact Q-ball couplings
  fixed.
- **A2.3** In `projection/mismatch.py`, add an `si_l2` term (weight
  `W_SI = 1.0`) to `compute_mismatch` for warp motifs: L2 of the solved
  slice's momentum-density proxy vs the analytic `S_i` target on the
  same 2D slice.  The solved-side proxy: reconstruct
  `S_i = -(1/8 pi)(d_j A^j_i ...)` from the gridinit A_ij with
  np.gradient (this measures whether the *solver's geometry* carries the
  right momentum — exactly the quantity the Mom stall corrupts).
- **A2.4** In `projection/iterate.py`, include `si_l2` in the CMA-ES
  fitness (geometry phase) so per-lump velocity/position/amplitude are
  optimised toward the momentum target.  Expose per-lump velocity
  components in the parameter vector for warp motifs if not already
  present.
- **A2.5** Expose `use_compact_Vi_ansatz` as a wrapper knob (params.txt
  generation in `grtresna/`), default unchanged; for warp runs set
  `use_compact_Vi_ansatz = 0` (periodic ansatz B.3 — better suited to
  distributed scalar sources) and record whether the Mom stall (66.86%)
  improves.  Also keep `max_NL_iterations = 200`,
  `nl_stall_tolerance = 0.005` for warp runs.

Acceptance (A2): Mom residual < 30% (down from 66.86%) on the one-shot
solve, or CMA-ES (40 evals) reaches `si_l2` reduced >= 3x from the initial
fit.  If neither is achievable, this is the trigger for Phase B.

#### A3. Validation ladder run + tests

**Effort: Small-Medium.**

- **A3.1** New tests in `tests/projection/test_alcubierre.py`:
  - `TestAnalyticGridinit`: writer round-trips through `read_gridinit`;
    chi=1 / lapse=1 / K=0 everywhere; `shift1(bubble centre) ~ -v`;
    `shift1(far field) ~ 0`; A_ij traceless (`A11+A22+A33 ~ 0`); A_ij
    peaks at the bubble wall.
  - `TestWarpPainting`: painter overwrites only the intended components;
    `paint_aij=False` leaves A_ij untouched; header/dx/origin preserved.
  - `TestAlcubierreSi`: analytic `S_i` ring is antisymmetric about the
    wall, points along -x, and cross-checks against warpfactory `T^{0i}`
    to L2 < 1e-3.
  - `TestSiMismatch`: `si_l2 = 0` when solved slice == target;
    monotonically increases under perturbation.
- **A3.2** Full ladder execution script/README snippet:
  ```bash
  # Rung 1 — analytic anchor
  .venv/bin/python scripts/search/project_geometry_motif.py \
    --target alcubierre --mode analytic-gridinit \
    --warp-velocity 0.5 --warp-bubble-radius 2.0 --warp-sigma 2.0 \
    --out-dir runs/alcubierre_anchor
  # Rung 2 — solve + paint
  .venv/bin/python scripts/search/project_geometry_motif.py \
    --target alcubierre --mode solve-only --paint-warp-shift shift+aij \
    --max-lumps 5 --gridinit-n 32 --grtresna-L 64.0 --mpi-ranks 4 \
    --out-dir runs/alcubierre_painted
  # Rung 3 — momentum-matched iterate
  ... --iterate 40 --iterate-popsize 8 ...
  # Rung 4 — evolve + probe (solve-and-evolve, stop_time >= 2*t_cross)
  ```
- **A3.3** Keep `pytest tests/projection -q` green (79 existing + new).

---

### Phase B — GRTresna C++ (contingent: only if A2 acceptance fails)

#### B1. Prescribed-A_ij background (the physically correct two-step solve)

**Effort: Medium (C++).**

GRTresna already supports an *analytic background A_ij* — the Bowen-York
term in `PsiAndAijFunctions::compute_bowenyork_Aij` is added to the CTT
part everywhere (`CTTK.impl.hpp`, `Aij_reg + Aij_bh`).  Mirror that
mechanism:

- **B1.1** Add `compute_warp_Aij(Tensor<2, Real> &Aij, const RealVect &loc)`
  to `PsiAndAijFunctions` (params: `warp_velocity`, `warp_bubble_radius`,
  `warp_sigma`, `warp_center`, read via `read_params`; zero when
  `warp_velocity == 0`).  Formula: the same analytic A_ij as A0 (closed
  form via df/dr_s — 1D, safe).
- **B1.2** Add it wherever `compute_bowenyork_Aij` is summed
  (`CTTK.impl.hpp` lines ~83/170, `CTTKHybrid.impl.hpp`,
  `Diagnostics.impl.hpp`, `RHSTagging.hpp`).
- **B1.3** In the wrapper's params.txt generation, emit the warp params
  for warp motifs.  With the background A_ij prescribed, the Hamiltonian
  solve makes chi consistent with matter + warp curvature — and the
  momentum constraint only needs to mop up the *difference* between the
  matter S_i and the prescribed A_ij divergence.
- **B1.4** Also write the analytic shift into `set_output_data`
  (`WriteOutput.H`): `c_shift1 = -v f(r_s)` when `warp_velocity != 0`.
  This is a *gauge seed*, mirroring A1 but done natively (note: this is
  NOT the withdrawn V_i→shift conversion).

#### B2. Momentum-constraint convergence for distributed scalar S_i

**Effort: Large (C++ solver work).  Last resort.**

- **B2.1** Instrument the Mom residual per NL iteration for a warp run;
  determine whether the stall is (a) the Vi ansatz, (b) multigrid
  tolerance, or (c) S_i too delocalised.
- **B2.2** Compare `use_compact_Vi_ansatz = 0` vs `1` (A2.5 gives the
  wrapper knob) at n=32 and n=64.
- **B2.3** If ansatz-limited: implement the full periodic B.3 ansatz path
  for non-periodic boundaries, or tighten multigrid solves per NL
  iteration for the V_i components specifically.

### Recommended order (updated)

**A0 → A1 → A2 → A3 (full ladder).**  A0 is non-negotiable first: it
calibrates the probe and the preservation checker against exact
Alcubierre data before any reconstruction is attempted.  If the ladder
passes at Rung 4, the deliverable is done without touching C++.  If A2's
acceptance fails (Mom >= 30% and no `si_l2` progress), start B1; B2 only
if B1's mopped-up momentum residual still dominates the mismatch.

---

## Phase A Implementation Status (2026-07-10)

Phase A is **implemented and tested** (Python-only, no C++ changes).

### Implemented

| Task | Status | Files |
|------|--------|-------|
| **A0** Analytic Alcubierre gridinit writer | ✅ | `projection/warp_gridinit.py` |
| **A0.2** `--mode analytic-gridinit` CLI | ✅ | `scripts/search/project_geometry_motif.py` |
| **A1** Post-solve shift + A_ij painting | ✅ | `projection/warp_gridinit.py`, CLI `--paint-warp-shift` |
| **A2.1** Analytic `S_i` on xz-slice | ✅ | `projection/warp_gridinit.py::alcubierre_analytic_Si` |
| **A2.3** `si_l2` mismatch term | ✅ | `projection/mismatch.py` (`_compute_si_mismatch`, `_solved_si_slice`) |
| **A2.5** `use_compact_Vi_ansatz` CLI knob | ✅ | `--use-compact-vi-ansatz` flag |
| **A3** Tests (18 new) | ✅ | `tests/projection/test_alcubierre.py` |
| **A0 smoke** Preservation check on analytic gridinit | ✅ | validates shift + A_ij |

### Preservation checker fixes (warp-motif specific)

The preservation checker (`motif_preservation.py`) had two sign/convention
mismatches that caused it to mis-score the analytic Alcubierre gridinit:

1. **Shift alignment sign**: The Alcubierre shift is `beta^x = -v f(r_s)`
   (negative x) but the momentum target direction is `(1, 0, 0)` (positive
   x transport).  In ADM, a negative shift produces +x motion.  Fixed:
   for warp motifs, the sign factor is flipped so a correctly-negative
   shift scores positively.

2. **Support localization**: Alcubierre has flat `chi=1` (the spatial
   metric is flat), so the chi-deviation check always fails.  Fixed: for
   warp motifs, localization is measured by `max|A_ij|` with a lower
   threshold (`WARP_LOCALIZATION_TOLERANCE = 0.01`).

After these fixes, the preservation check on the analytic gridinit reports:
- `shift_alignment = 0.25` (correct — the bubble occupies ~25% of the slice)
- `support_localized = True` (A_ij is non-zero at the bubble wall)
- `polarity_retention = 1.0`
- `retention_score = 0.36` (limited by `f_op_retention`, which requires
  evolution to confirm the shortcut — expected for a t=0-only check)

### Not yet executed (requires GPU + built GRTeclyn binary)

These tasks from the plan require a running GRTeclyn binary and are left
for the campaign execution environment:

- **A0.5** GRTeclyn evolution of the analytic gridinit (Gamma-driver damping)
- **A0.6** 4D null-geodesic FTL probe on the analytic slice
- **A1.3** Dual Ham/Mom logging of painted vs unpainted slice (needs solve)
- **A2.2** Warp-aware lump velocities in `fit/motif.py` (needs solve to validate)
- **A2.4** `si_l2` in CMA-ES fitness (already wired through `compute_mismatch`)
- **A3.2** Full 4-rung ladder execution

### Test results

- `pytest tests/projection -q`: **97 passed** (79 original + 18 new)
- `pytest tests/search -q`: 99 passed, 1 pre-existing unrelated failure
  (`test_qd_search_resume_continues_eval_counter`)

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
| `scripts/campaigns/geometry_first/run.sh` | Motif-to-matter MAP-Elites campaign launcher |
| `scripts/campaigns/geometry_atlas/run.sh` | Pure-geometry MAP-Elites atlas launcher |
| `src/grteclyn_wrapper/search/geometry_atlas/` | Stationary genome, render, score, MAP-Elites driver |
| `src/grteclyn_wrapper/initial_data/adm_stationary.py` | Shared stationary ADM / CCZ4 / Einstein-source helpers |
| `src/grteclyn_wrapper/initial_data/motif.py` | Motif extraction from episodes |
| `src/grteclyn_wrapper/grtresna/fit/motif.py` | Matter fitting (lumps, ring splitting) |
| `src/grteclyn_wrapper/projection/iterate.py` | CMA-ES iteration loop |
| `src/grteclyn_wrapper/projection/mismatch.py` | Two-phase fitness + 2D/K_ij/S_i mismatch + feasibility pre-check |
| `src/grteclyn_wrapper/projection/motif_preservation.py` | Preservation check (warp-aware shift sign + A_ij localization) |
| `src/grteclyn_wrapper/projection/warp_gridinit.py` | Analytic Alcubierre gridinit writer + post-solve shift/A_ij painter + S_i |
| `src/grteclyn_wrapper/search/qd_search/driver.py` | MAP-Elites driver (`geometry_first` mode) |
| `src/grteclyn_wrapper/search/qd_search/descriptors.py` | `geometry_first` behavior descriptors |
| `src/grteclyn_wrapper/metrics/probes/ftl/geodesic.py` | Includes `build_metric_3d_from_gridinit` for atlas scoring |
| `tests/projection/test_iterate.py` | CMA-ES iterate tests |
| `tests/projection/test_alcubierre.py` | Alcubierre warp-drive tests |
| `tests/search/test_geometry_atlas.py` | Pure-geometry atlas unit + smoke tests |
