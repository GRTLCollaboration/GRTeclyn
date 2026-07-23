# Geometry-First Lab Journal

## How it works (overview)

A three-part loop that never runs a physics solver — it *paints* candidate
spacetimes and measures them.

**1. Genome → geometry.** A flat coefficient vector in two parts:
* **RBF bumps** — control centers on a sphere, each adding a smooth Gaussian to
  lapse / shift / spatial-metric (free-form "clay").
* **Analytic tail** — named coherent topologies the clay cannot form: Alcubierre
  shift tube, compression tunnel, conformal lens, wormhole throat (`ansatz.py`).

**2. Render → legal metric** (`render.py`). Two constructions keep every
candidate physical:
* `gamma = exp(S)` ⇒ the spatial metric is *always* positive-definite (no
  signature flips).
* a **compact envelope** forces the deformation to zero at a fixed radius ⇒
  asymptotically flat by construction (cannot cheat via the box boundary).

**3. Probe → score** (`score.py`). A frozen null-geodesic probe fires light rays
and measures how much faster than flat-space light they arrive → `f_geo`. The
required negative (exotic) energy is computed analytically. Each candidate is
filed into a **MAP-Elites archive** indexed by *family* (shift-fraction) ×
*cost* (log exotic energy); each cell keeps its best `f_geo`.

**Seeds.** Three entry points are injected at start (`driver.py`): exact
Minkowski (`zero_genome`), a hand-written Alcubierre bubble
(`seed_alcubierre_genome`), and the trusted CMA champion JSON (`--seed-genome`).
Seeds are starting points, not answers.

**Evolution on top of seeds.** Each round MAP-Elites picks archive elites,
applies bounded Gaussian mutation (`mutate_genome`), and injects a fraction of
fresh random genomes — genuine exploration beyond the seeds.

**Anti-"wibbly-wobbly" guards** (stop the optimizer gaming numerical junk):
* **Hard rejects** (`_rank_score`, score `-1e9`): non-positive lapse,
  non-SPD metric, non-flat boundary (>5% α/β at edge), or inconsistent geometry
  (`ham_l2 > 5`).
* **Null-Hamiltonian drift check**: under-resolved / sharp geometries lose the
  null condition; those results are flagged untrusted and their `f_geo` **halved**
  (exploiting integrator error is penalised, not rewarded).
* **Resolution honesty**: run at `h=0.5`; champion `f_geo` confirmed stable on an
  `N=48/64/96` ladder (not a grid artifact).
* **Exotic-energy axis**: geometries demanding absurd negative energy land in a
  separate expensive cell instead of dominating.

**Honest caveats.** The drift penalty *softens* (halves) rather than zeros, and
`f_geo` is a **frozen-slice screen**, not a verdict. This is deliberate: Stage-1
is a cheap breadth scout. Two harder physical filters sit downstream — **Stage-2**
(GRTresna constraint solve: can matter source it?) and **Stage-3** (GRTeclyn
evolution: does it survive in time?) — which eliminate any geometry that merely
games the frozen probe.

---

## 2026-07-23 — Pure-geometry MAP-Elites atlas MVP + first wide hunt

**Motivation.** Matter-first search (paper pipeline) only finds geometries that the Q-ball ansatz can source. Stage-1 inverse design flips this: scout broad stationary metrics first, score frozen shortcut strength vs exotic cost, then later ask matter/solvers how to build the elites.

**What was done.** Implemented CPU `geometry_atlas` (RBF genome → ADM/CCZ4 gridinit → frozen `f_geo`, optional `f_ff`, MAP-Elites archive). No GRTresna, no GRTeclyn. Ran first wide hunt:

| Knob | Value |
|------|-------|
| Run | `runs/geometry_atlas/geometry_atlas_wide_n32_e200_20260723T015531Z/` |
| Evals | 200 (`NO_FF=1`) |
| Grid | \(N=32\), \(L=64\), bins \(8\times8\) |

**Results.**

| | |
|--|--|
| Best score | **8.34** (\(\approx 100\times f_{\rm geo}\); FF off) |
| Best frozen \(f_{\rm geo}\) | **0.083** (~8.3%) |
| Elites / coverage | 3 / 4.7% |
| Best cell | `[0,5]` (low shortcut bin × mid exotic cost) |

Screening only — not comparable to paper \(S\sim852\) / evolved \(f_{\rm geo}=20.2\%\) (cand. 146). Next: free-fall ranking pass and/or Stage-2 matter infill on elites.

## 2026-07-23 — Probe calibration + Alcubierre/\(K_{ij}\) genome + CMA-ES

**Motivation.** First wide hunt stalled at frozen \(f_{\rm geo}=8.3\%\) because endpoints spanned the full box, the genome lacked Alcubierre controls, and MAP-Elites was exploring a poorly conditioned space without a hill-climb.

**What was done.**
1. Localised probe endpoints about the support (`probe_half_length`).
2. Calibration anchors: analytic Alcubierre + cand.146 `initial_data.gridinit`.
3. Genome upgrade: trailing Alcubierre controls + per-center \(K_{ij}\) RBF modes.
4. Focused CMA-ES refine mode (`MODE=cmaes`).

**Calibration** (`runs/geometry_atlas/geometry_atlas_calibration_20260723/`):

| Case | frozen \(f_{\rm geo}\) |
|------|------------------------|
| Alcubierre localised | **0.220** |
| Alcubierre full-box | 0.140 |
| Cand.146 t=0 IVP (subsampled) | 0.019 |

Probe OK. Localisation gain \(+8\) pp. Cand.146 t=0 is *not* the paper frozen-peak \(29.45\%\) (evolved snapshot).

**CMA-ES** (`geometry_cmaes_n32_e24_fixalc_20260723/`, \(N=32\), \(L=64\), 24 evals, `NO_FF=1`): best frozen \(f_{\rm geo}=0.133\) (seed Alcubierre; CMA did not beat seed in 24 steps). Beats the earlier RBF-only hunt (0.083); still below the denser Alcubierre calibration grid.

## 2026-07-23 — Maximise frozen \(f_{\rm geo}\) (Alcubierre-only CMA)

**Motivation.** Prior CMA stalled because (i) \(v\gtrsim1.2\) null-H drift zeroed \(f_{\rm geo}\), and (ii) best-direction picked flat transverse axes over the warp axis.

**Fixes.** Rank axes by raw \(f_{\rm geo}\); soft-penalise h-drift; CMA objective = raw \(f_{\rm geo}\) with mild drift cost; `--alc-only` 4-D search.

**Run.** `geometry_cmaes_maxfgeo_n48_L32_e100_20260723/` — \(N=48\), \(L=32\), 100 evals, `ALC_ONLY=1`, \(v_{\max}=2\).

| | frozen \(f_{\rm geo}\) | notes |
|--|--|--|
| Best raw | **0.436** (eval 63) | \(v\approx1.96\), \(R\approx8\), h-drift high |
| Best trusted | **0.413** (eval 70) | h-quality OK |
| Saved best (obj) | 0.217 stored / ~0.43 raw | `best.gridinit` |

Beats paper cand.146 frozen peak \(29.45\%\) under this stationary Alcubierre probe (screening only; not evolved 4D). Artifacts: `best_trusted.gridinit`, `best_trusted_genome.json`.

## 2026-07-23 — Breadth MAP-Elites at \(h=0.5\) (morphology axes)

**Motivation.** Goal is not Alcubierre tuning alone: fill an atlas of geometry *families* that maximise FTL metrics. Replaced \(f_{\rm geo}\) as a descriptor with shift-fraction morphology; warm-started from the trusted CMA elite.

**Run.** `geometry_atlas_breadth_n64_L32_e220_20260723/` — \(N=64\), \(L=32\) (\(h=0.5\)), 220 evals, \(n_{\rm centers}=3\), bins \(8\times8\), `NO_FF=1`, seed `cma_maxfgeo_trusted.json`.

| | |
|--|--|
| Coverage | **17 / 64** (26.6%) |
| QD-score \(\sum f_{\rm geo}\) | **1.805** |
| Best frozen \(f_{\rm geo}\) | **0.448** (cell `[7,6]`, shift≈1) |
| Mid-morphology best | 0.116 (cell `[5,5]`, shift≈0.63) |

**Resolution check.** Trusted CMA seed re-scored at \(N=48/64/96\) (\(L=32\), \(h=0.67/0.50/0.33\)): \(f_{\rm geo}\approx0.41\) stable to \(\pm0.02\) — resolution-converged, not a coarse-grid artifact.

**What the search converged to.** The top-20 evals live almost entirely in cell `[7,6]` (16/20) and `[7,5]` (4/20): shift-fraction \(\approx1\), mid–high exotic cost. Best-so-far jumped to \(0.43\) at the warm-started seed (evals 0–6) and only crept to \(0.448\) over 220 evals via local mutation of the Alcubierre channel. Lower-shift bins (0–6) were reached but plateaued at \(f_{\rm geo}\lesssim0.12\). **Conclusion: under the frozen null-geodesic probe, shortcut strength is monopolised by the shift channel (Alcubierre-like tubes); curvature/lapse-only "conformal lens" families are genuinely weaker, not merely under-sampled.**

**Corrections to earlier entries (coarse-resolution).** The first wide hunt (\(N=32\), \(L=64\) → \(h=2.0\)) and its \(f_{\rm geo}=0.083\) were reported as a screening ceiling. That number was **not physical**: at \(h=2.0\) the warp wall is unresolved and the full-box endpoints dilute the ray, so it under-reads the score badly. Corrected view: the same geometry family scores \(f_{\rm geo}\approx0.41\)–\(0.45\) once probed at \(h\le0.5\) with localised endpoints. The earlier claim of a "resolution artifact making finer grids score *lower*" was also wrong — it came from a run where the warm-start genome was silently skipped and a weaker default seed was scored; re-scoring the correct genome shows \(f_{\rm geo}\) is flat across \(h=0.67\to0.33\). Treat pre-`h=0.5` \(f_{\rm geo}\) figures in this journal as lower bounds only.

**Next steps.**
1. **Stage-2 matter infill** on the 3 shift-heavy elites (cells `[7,5/6/7]`) — hand off to `project_geometry_motif` + GRTresna constraint solve; this is the real test of whether the geometry is sourceable.
2. **Widen the genome away from shift** — the atlas is telling us curvature-only families are weak; either accept that (report it as a finding) or add richer non-shift ansätze (e.g. wormhole-throat / conformal-well modes) before spending more QD budget.
3. **Free-fall ranking pass** (`NO_FF=0`) on the current elites to check the shortcut survives a physical observer, not just the frozen null probe.
4. **Longer/finer QD only if** step 2 adds expressive non-shift modes; otherwise QD budget is wasted re-discovering the same `[7,*]` tube.

## 2026-07-23 — Beyond Alcubierre: analytic topology registry

**Motivation.** The atlas kept re-discovering the shift tube because the only *coherent* analytic mode was Alcubierre; the RBF `S`-bumps are isotropic Gaussians that cannot organise into a tunnel or throat. To test whether pure-curvature shortcuts exist, add explicit non-shift topologies.

**What was done.** Refactored the trailing Alcubierre block into an open/closed **analytic ansatz registry** (`ansatz.py`): each topology is a pure function `(points) -> (dalpha, dbeta, dS)` composed into the genome, with `gamma = expm(S)` guaranteeing SPD + asymptotic flatness. Added three shift-free families:

| Ansatz | Mechanism | Params |
|--|--|--|
| Compression tunnel | anisotropic `gamma` contraction along an axis | strength, width, axis |
| Conformal lens | isotropic `psi^4` well in a spherical shell | strength, radius, sigma |
| Wormhole throat | radial (`r̂⊗r̂`) contraction shell | strength, throat_radius, sigma |

Genome tail grew 4 → 13 params; legacy Alcubierre-only elites migrate transparently (`_migrate_coeffs`). CLI ablation flags `--no-tunnel/--no-lens/--no-throat`; CMA `--alc-only` now frees the whole analytic tail. 19/19 unit tests pass.

**Single-parameter screening** (`N=48`, `L=32`, localised probe), each topology alone with **zero shift**:

| Ansatz | frozen \(f_{\rm geo}\) | shift frac |
|--|--|--|
| Tunnel | **0.094** | 0.000 |
| Lens | 0.055 | 0.000 |
| Throat | 0.035 | 0.000 |

Curvature-only shortcuts are real (the tunnel already rivals the prior best mid-morphology elite, \(0.12\), at a single untuned setting). **Corrects the earlier "curvature families are genuinely weak" conclusion — they were unreachable, not weak.**

**Run launched.** `geometry_atlas_topologies_n64_L32_e260_20260723/` — all topologies on, \(N=64\), \(L=32\), 260 evals, warm-started from the trusted CMA elite. Results pending.
