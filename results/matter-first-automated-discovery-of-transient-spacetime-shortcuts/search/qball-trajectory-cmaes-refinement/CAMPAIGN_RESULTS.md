# Campaign results: `qball_traj_fgeo_depth_cmaes_v1` — CMA-ES depth refinement

**Run:** `runs/grtresna_cmaes/qball_traj_fgeo_depth_cmaes_v1`
**Window:** 2026-08-06 03:58 → 21:25 UTC (17 h 27 m wall clock, single node, single rank)
**Budget:** 200 evaluations · 13 generations × population 16 (final generation truncated to 8 by the eval cap)
**Objective:** `f_geo_depth` — raw uncapped evolving-geodesic path saving, identical to the preceding QD campaign
**Warm start:** best row of the aborted first attempt's generation 1 (`f_geo` = 0.4628)

**Best result: eval 195 — `f_geo` = 0.4838 (48.38 % path saving), score 4862.61.**
That is **+2.5 percentage points** over the MAP-Elites champion (45.90 %) and **+2.1 points** over the
warm-start seed it began from. The improvement was **still accelerating when the eval budget ran out** —
the single best result of the whole campaign came from the final, half-sized generation.

Two caveats carried through the whole document and stated once here: (a) roughly **four fifths of the
gain is clean and one fifth coincides with a jump to a different initial-data solution branch that
needs 40 % more exotic matter** (§3.3); (b) the depth measurement remains **clipped by the emission
window** — 48.38 % is a lower bound, not the recipe's peak (§3.4).

---

## 1. What this campaign was for

The preceding MAP-Elites campaign (`qball_traj_fgeo_depth_v1`, 200 evals) established that its own
record recipe was inherited unchanged from its warm start: 200 broad map evaluations never improved on
it, and the genome distance from seed to champion was exactly 0.000. That is the signature of a search
whose step sizes are too large for the local structure — it can find the basin but cannot climb inside
it.

This campaign is the direct follow-up: hand the same objective and the same physics to a **local,
directed optimizer** (CMA-ES, initial step size 0.05 of each parameter's range) and let it climb.

Search space, physics, gates and scoring are byte-for-byte the QD setup — 39 free dimensions describing
five Q-ball lump trajectories (per-lump orbital radius, phase, two tilt angles, rotation rate, radial
velocity, exotic-matter dial; plus global breathing and z-bobbing amplitude/frequency), bicomplex scalar
matter, five lumps, retrograde-only orbits.

Everything is measured the same way as before: a light ray is traced end-to-end through the **evolving**
spacetime between antipodal points on the r = 8 shell, and the score is how much earlier it arrives than
the same ray in flat space (14.4 time units). Honesty gates are unchanged — the complete ray bundle must
arrive (3 of 3), the numerical constraint drift must stay small, and the initial data must satisfy the
Hamiltonian constraint below the acceptance threshold of 0.03 or the candidate is rejected before it
ever reaches a GPU.

---

## 2. Headline numbers

| Quantity | Value |
|---|---|
| Evaluations | 200 / 200 (13 generations, last one truncated to 8 by the eval cap) |
| Wall clock | 17 h 27 m |
| **Best depth `f_geo`** | **0.483815** (eval 195) — 48.38 % path saving |
| Best score | 4862.61 |
| Best arrival time | 25.433 vs 14.400 flat-space reference, ray emitted at t = 18 |
| Warm-start seed | 0.462849 (score 4648.41) |
| Previous QD champion | 0.459 (score 4609.29) |
| **Gain over seed** | **+2.10 points** of path saving |
| **Gain over QD champion** | **+2.5 points** |
| Best on the clean branch (§3.3) | 0.479694 (eval 185) — +1.68 over the seed with no branch change |
| Clean evaluations | 175 / 200 |
| Crashed during evolution | 14 |
| Rejected at the constraint gate before GPU | 11 |
| **Solver timeouts** | **0** |
| Retained eval directories | 12 (top 10 by score + 2 diagnostics) |

### Generation-by-generation

| Gen | Best depth | Median depth (in basin) | Clean | Crash | Gate reject |
|---|---|---|---|---|---|
| 1 | 0.4661 | 0.4628 | 11 | 2 | 3 |
| 2 | 0.4669 | 0.4636 | 14 | 1 | 1 |
| 3 | 0.4680 | 0.4653 | 14 | 1 | 1 |
| 4 | 0.4716 | 0.4669 | 16 | 0 | 0 |
| 5 | 0.4715 | 0.4682 | 14 | 1 | 1 |
| 6 | 0.4716 | 0.4707 | 13 | 2 | 1 |
| 7 | 0.4736 | 0.4716 | 15 | 0 | 1 |
| 8 | 0.4763 | 0.4737 | 14 | 2 | 0 |
| 9 | 0.4772 | 0.4745 | 15 | 0 | 1 |
| 10 | 0.4777 | 0.4762 | 14 | 1 | 1 |
| 11 | 0.4799 | 0.4777 | 14 | 1 | 1 |
| 12 | 0.4809 | 0.4787 | 14 | 2 | 0 |
| **13** | **0.4838** | **0.4810** | 7 | 1 | 0 |

Full per-evaluation and per-generation tables: [`analysis/per_eval.csv`](analysis/per_eval.csv),
[`analysis/per_generation.csv`](analysis/per_generation.csv).

---

## 3. Findings

### 3.1 The search did not converge — it ran out of budget

The **median** in-basin result rose monotonically in all thirteen generations, 0.4628 → 0.4810, without
a single reversal. The best rose in twelve of thirteen (generation 5 was flat at the fourth decimal).
The three deepest results of the entire campaign — 0.4838, 0.4829, 0.4822 — are evals 195, 196 and 197,
all from the last generation, which only got 8 of its 16 slots before the eval cap stopped the run.

Two independent signals say the same thing. A converged optimizer shows the best flattening while the
median catches up to it; here **both curves climbed together at a roughly constant rate to the last
measurement**, with the gap between them (~0.003) essentially unchanged from generation 6 onward. The
population was translating, not contracting.

The practical reading: 200 evaluations was a budget decision, not a physical limit. Extending this run
would be expected to continue producing gains at roughly 0.1–0.4 points per generation. Where that stops
is not determined by this data.

### 3.2 The landscape is a very narrow ridge — which explains the QD result

The winning genome sits a normalized distance of just **0.0396** from the seed across all 39 dimensions.
The single largest parameter change is **2.1 % of that parameter's allowed range**; 38 of 39 dimensions
moved less than 2 %. That tiny, coordinated displacement bought +2.1 points of path saving.

| Largest parameter moves, seed → winner | From | To | % of range |
|---|---|---|---|
| `trajectory_lump4_tilt_phi` | 3.7854 | 3.6520 | −2.1 % |
| `trajectory_lump1_phase0` | 4.2143 | 4.2994 | +1.4 % |
| `trajectory_lump4_phase0` | 5.3241 | 5.2454 | −1.3 % |
| `trajectory_lump3_tilt_phi` | 1.3885 | 1.3219 | −1.1 % |
| `trajectory_lump0_R0` | 5.9043 | 5.8391 | −1.0 % |
| `trajectory_lump4_R0` | 7.1113 | 7.1700 | +0.9 % |

This retroactively explains the MAP-Elites null result rather than contradicting it. The productive
direction is a narrow ridge in a 39-dimensional space; a broad sampler that perturbs parameters by tens
of percent steps off it every time and lands back at or below the seed. It is not that the map search
was run badly — the two methods were doing different jobs, and this is the job the local optimizer is
for. The pairing (map search finds the basin, local search climbs it) worked exactly as intended.

The corollary is a warning: because the ridge is this narrow, **the result is sensitive to anything that
perturbs the geometry at the sub-percent level** — resolution, domain size, interpolation. That is a
direct argument for validating the winner at higher resolution before treating 48.38 % as established
(§5, §8).

### 3.3 The gain splits into a clean part and a branch-change part

This is the most important caveat in the campaign, and it needs stating precisely.

Across the 147 in-basin evaluations, two diagnostics of the initial data — the Hamiltonian constraint
residual and the total amount of negative energy density the configuration requires — do not vary
smoothly. They fall into **two tight clusters separated by a wide gap**:

| Branch | Constraint residual | Exotic-matter integral | Count | Depth range | Median depth |
|---|---|---|---|---|---|
| **A** | ≈ 0.0162 | 30.10 – 30.20 | 106 | 0.4589 – 0.4797 | 0.4684 |
| **B** | ≈ 0.0286 | 42.32 – 42.41 | 41 | 0.4736 – 0.4838 | 0.4776 |

Within each branch the exotic-matter integral varies by about 0.3 %; between branches it jumps 40 %.
Branch B also sits at a constraint residual 76 % higher than branch A — **0.0286 against an acceptance
gate of 0.0300**, i.e. 95 % of the way to rejection.

The two branches are not two different recipes. The nearest cross-branch pair, eval 185 (branch A) and
eval 195 (branch B), are separated by a normalized genome distance of **0.0206**, with no single
parameter differing by more than 1.1 % of its range — yet they land on opposite sides of the gap. The
initial-data solver is finding a discretely different solution for near-identical inputs. Inspecting the
resulting lump placements confirms it: eval 195's configuration puts two lumps far off the equatorial
plane (z ≈ −5.8 and −3.6) that eval 185's does not have, despite the orbital parameters being nearly
identical — a bifurcation, most plausibly driven by the strong z-bobbing (amplitude 2.57) making lump
positions a stiff function of orbital phase.

Branch B first appears at eval 122 (generation 8) as a single member and then takes over the population:
1 → 5 → 8 → 10 → 11 → 6 of 6 by the final generation. The optimizer found it and preferred it, because
it is genuinely deeper.

**Decomposition of the +2.10 point gain:**

- **+1.68 points are clean.** Branch A alone climbed 0.4628 → 0.4797 across generations 1–12 with its
  constraint residual pinned at 0.0162 and its exotic-matter requirement pinned at ~30.15. This part of
  the gain is bought with nothing but better lump choreography.
- **+0.41 points came with the branch switch.** Within the generations where both branches coexisted
  (9–12), branch B beat the best branch-A member of the same generation by +0.14 to +0.42 points — a
  real but modest edge, purchased with 40 % more exotic matter and a constraint residual pressed against
  the acceptance gate.

For any claim that needs to be conservative, **the defensible refined number is 47.97 % (eval 185,
branch A)**, not 48.38 %. The 48.38 % result is not disqualified — it passes every gate the campaign
enforces, with a full 3-of-3 ray bundle and a numerical drift of 0.0011 — but it is the deepest point of
a branch whose initial data is measurably worse conditioned, and it should not be promoted without a
higher-resolution check that its constraint residual does not drift past 0.03 (§8). Both genomes are
preserved in [`genomes/`](genomes/) for exactly this reason.

A note on what this is *not*: the branch-B configurations were not sneaking past the gate. The gate did
its job — it rejected 11 candidates outright, and it admitted these. The concern is narrower: the
optimizer is now operating close enough to the acceptance boundary that the boundary, rather than the
physics, could start shaping the result.

### 3.4 Depth is still clipped by the measurement window

The winner's emission sweep is **monotonically increasing to the very last launch**:

| Emission time | 0 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | **18** |
|---|---|---|---|---|---|---|---|---|---|---|
| Path saving | 0.213 | 0.243 | 0.273 | 0.305 | 0.339 | 0.372 | 0.404 | 0.433 | 0.460 | **0.484** |

Ten launches at intervals of 2, the last at t = 18, and the deepest result is the last one. There is no
turnover anywhere in the sweep — every extra unit of time the geometry is allowed to develop before the
ray is emitted buys roughly 0.015 more path saving, with no sign of saturation.

This reproduces the QD campaign's finding at the new depth: **48.38 % is a lower bound set by where the
measurement stopped, not by where the physics stopped.** The recipe was never asked what it can do after
t = 18.

The gradient is steady enough (≈ 0.0135 per unit time over the last four launches) that extrapolation is
tempting; it should be resisted, because a turnover has never been observed and its absence over
t ∈ [0, 18] says nothing about t > 18. The correct response is to measure it, which is what the planned
t = 64 promotion run does.

### 3.5 Recorded diagnostics behaved as designed

`f_geo_depth` scores path saving and nothing else. Matter retention, confinement, coordinate speeds and
persistence are **logged as diagnostics and carry no weight in the objective** — they are recorded so
that cross-campaign comparisons are possible, not optimized toward or away from.

Their values are unremarkable and consistent between the seed and the winner:

| Diagnostic | Seed | Winner |
|---|---|---|
| Final confined fraction | 0.180 | 0.171 |
| RMS radius spread ratio | 3.33 | 3.31 |
| Peak local coordinate speed (min over samples) | 1.895 | 1.897 |
| Peak shift magnitude | 0.161 | 0.168 |
| Numerical constraint drift on the ray | 0.00103 | 0.00108 |

Comparing the twenty deepest results against the other 127 in-basin results, every one of these
diagnostics is statistically indistinguishable between the groups — the only quantities that separate
the deep set are the two branch markers of §3.3. In other words, within this recipe family, depth is not
trading against any recorded diagnostic; it is trading against exotic-matter content and initial-data
conditioning specifically.

One methodological caution for anyone reading `analysis/per_eval.csv`: the correlations between depth
and the physics diagnostics look strong (|r| = 0.5–0.9) when computed over all clean rows, but that is
an artifact of the failed and degenerate runs at the bottom of the range dragging every metric down
together. Restricted to the in-basin population where the comparison is meaningful, the apparent
relationships are the §3.3 branch structure and little else.

### 3.6 Failure taxonomy, and the infrastructure fixes that held

| Outcome | Count | Note |
|---|---|---|
| Clean evaluation | 175 | |
| Crashed during GPU evolution | 14 | Scattered across generations, no clustering; normal exploration cost |
| Rejected at the constraint gate | 11 | Offspring stepping over the 0.03 acceptance threshold |
| **Solver timeout** | **0** | |

Two infrastructure problems were fixed immediately before this run, and both held for its full duration.

**Solver timeout.** The first attempt at this campaign lost an entire generation — six consecutive
evaluations scoring −350 — when the initial-data solver exceeded its 900-second limit under six-way
concurrency. Raising the limit to 2400 s eliminated the failure class completely: zero timeouts in 200
evaluations.

**Pipeline starvation.** The first attempt also ran with population 4 on 4 GPUs. Because CMA-ES cannot
produce the next generation until the current one is fully evaluated, any candidate that failed fast or
was gate-rejected left a GPU idle for the rest of the generation — measured at roughly 50 % idle. Raising
the population to 16 (4 × GPU slots) restored streaming, and GPU occupancy through this run was
consistent with the MAP-Elites campaigns. This is now the documented default in both launchers and the
wrapper README: **population must be at least 4 × GPU slots and never equal to the GPU count.**

The gate-rejection rate of 11/200 (5.5 %) is low enough that no change to the initial step size is
warranted. Should it rise in a follow-up run, the correct response is to shrink the step size
(0.05 → 0.03), not to loosen the 0.03 acceptance gate — loosening it would break score comparability
with every campaign recorded so far, and given §3.3 it would specifically license the optimizer to buy
depth with worse-conditioned initial data.

---

## 4. Winner configuration (eval 195)

Full genome: [`genomes/winner_eval195_familyB.json`](genomes/winner_eval195_familyB.json).
Conservative branch-A alternative: [`genomes/best_clean_branch_eval185_familyA.json`](genomes/best_clean_branch_eval185_familyA.json).
Run directory extract: [`run/eval_000195/`](run/eval_000195/).

| Property | Value |
|---|---|
| Path saving `f_geo` | 0.483815 |
| Score | 4862.61 |
| Ray emitted / arrives | t = 18.0 → 25.433 (flat-space reference 14.400) |
| Ray bundle | 3 of 3 arrived, 0 captured |
| Numerical drift on the ray | 0.00108 (well inside tolerance) |
| Best probe axis | y |
| Hamiltonian residual (max / final) | 0.02858 / 0.01314 — branch B, 95 % of the 0.03 gate |
| Momentum residual (max / final) | 0.00147 / 0.00106 |
| Exotic-matter integral | 42.32 |
| Most negative energy density required | −0.03289 |
| Lump count / matter model | 5 / bicomplex scalar |
| Exotic dials (lumps 0–4) | 0.870 / 0.512 / 0.687 / 0.613 / 0.986 |
| Orbital radii (lumps 0–4) | 5.839 / 2.802 / 2.188 / 5.589 / 7.170 |
| z-bobbing amplitude / frequency | 2.568 / 0.810 |
| Breathing amplitude | 0.029 (essentially none) |

The physical picture is unchanged from the QD champion and worth restating: five strongly exotic-leaning
lumps on tilted retrograde orbits, bobbing hard out of the equatorial plane at a frequency close to the
orbital rate, with almost no radial breathing. The refinement did not change the mechanism; it tuned the
relative phases and tilts of the five lumps so that their contributions to the ray's path line up
slightly better.

---

## 5. Caveats

1. **Compare depths, not scores, across campaigns.** Scores are objective-mode specific. `f_geo`
   is the transferable number.
2. **48.38 % is a lower bound**, clipped by the emission window (§3.4). The recipe was never measured
   past t = 18.
3. **The conservative refined figure is 47.97 %** (§3.3). The extra 0.41 points came with a jump to an
   initial-data branch requiring 40 % more exotic matter, at 95 % of the constraint acceptance gate.
4. **This is a single-resolution result.** Everything here was measured at domain L = 64, N = 128,
   one refinement level, to t = 32. Given how narrow the productive ridge is (§3.2), a resolution study
   is required before the number is treated as converged. This is the primary purpose of the planned
   promotion run.
5. **The search was still improving when it stopped** (§3.1). 48.38 % is this budget's best, not this
   recipe family's ceiling.
6. **Coordinate-dependent quantities in the diagnostics** (local speeds, shift magnitudes) are gauge
   artifacts and are not evidence of anything physical on their own. They are recorded for continuity
   with earlier campaigns.
7. Retention, confinement and persistence are **recorded diagnostics with no weight in this objective**
   (§3.5). They should not be read as targets this campaign pursued or failed.

---

## 6. What's in this pack

```
qball-trajectory-cmaes-refinement/
├── CAMPAIGN_RESULTS.md          this document
├── README.md                    layout and headline
├── analysis/
│   ├── per_eval.csv             all 200 rows: generation, status, depth, constraints, branch label
│   └── per_generation.csv       13 rows: best/median depth, branch counts, failure counts
├── genomes/
│   ├── winner_eval195_familyB.json              the 48.38 % genome (39 dims + measurements)
│   ├── best_clean_branch_eval185_familyA.json   the conservative 47.97 % genome
│   └── warm_start_seed.json                     the 46.28 % starting point
├── logs/
│   └── cmaes_progress_evals_1-200.log           per-eval and per-generation progress lines
└── run/                         campaign extract, heavy binaries removed
    ├── trajectory.jsonl         all 200 scoreboard rows with full genomes
    ├── metadata.json            campaign configuration and search-space bounds
    ├── warm_start_gen1_seed.jsonl   the warm-start trajectory this run began from
    ├── ftl_champions.json, ftl_retention.jsonl
    └── eval_000193 … eval_000198    top-6 eval dirs: params.txt, score.json,
                                     small_data/ (emission sweep, timeseries),
                                     grtresna/ (solver residual history), run.log
```

Excluded on purpose, as with every pack in `results/`: `initial_data.gridinit` (531 MB per eval),
`small_data/metric_stack/` (16 MB per eval), plotfiles and checkpoints. Those stay in the gitignored
`runs/` tree on this machine. Pack size: 27 MB.

---

## 7. Run commands & properties

### Launcher

[`grteclyn-wrapper/scripts/campaigns/qball_trajectory/cmaes_run.sh`](../../../../grteclyn-wrapper/scripts/campaigns/qball_trajectory/cmaes_run.sh)
→ [`scripts/campaigns/cmaes/run.sh`](../../../../grteclyn-wrapper/scripts/campaigns/cmaes/run.sh)
→ `grteclyn_wrapper … optimize`.

The exact detached launch used for this run:

```bash
cd grteclyn-wrapper
setsid nohup env \
  RUN_NAME=qball_traj_fgeo_depth_cmaes_v1 \
  OBJECTIVE_MODE=f_geo_depth \
  WARM_START_TRAJECTORY=<runs>/grtresna_cmaes/qball_traj_fgeo_depth_cmaes_v1_gen1_seed.jsonl \
  WARM_START_TOP_K=1 WARM_START_JITTER=0.05 \
  POPULATION=16 MAX_GENERATIONS=13 TARGET_EVALS=200 \
  SIGMA0=0.05 SEED=7 \
  GPU_IDS="0 1 2 3" MAX_CONCURRENT_GRTRESNA=6 \
  GRTRESNA_TIMEOUT=2400 GRTRESNA_RANKS=1 \
  STOP_TIME=32 \
  GRTECLYN_GEO_EMIT_INTERVAL=2 GRTECLYN_GEO_MAX_EMISSIONS=10 \
  POSTLOAD_MAX_HAM_L2=3e-2 POSTLOAD_MAX_MOM_L2=1e-2 \
  bash scripts/campaigns/qball_trajectory/cmaes_run.sh \
  > <runs>/grtresna_cmaes/qball_traj_fgeo_depth_cmaes_v1_launch.log 2>&1 </dev/null &
```

### Settings that define this run

| Setting | Value | Why it matters |
|---|---|---|
| Objective | `f_geo_depth` | Raw uncapped path saving; identical to the QD campaign so depths compare |
| Population | 16 | 4 × GPU slots. Never set equal to the GPU count — see §3.6 |
| Generations / eval cap | 13 / 200 | The cap truncated generation 13 to 8 evaluations |
| Initial step size | 0.05 | Fraction of each parameter's range |
| Random seed | 7 | |
| Evolution grid | L = 64, N = 128, 1 refinement level, t = 32 | |
| Emission sweep | interval 2, 10 launches (t = 0 … 18) | The clipping in §3.4 comes from here |
| Probe directions | x, y, z (winner used y) | |
| Constraint gate | Hamiltonian 0.03, momentum 0.01 | Rejects before GPU; 11 rejections |
| Solver timeout | 2400 s | Raised from 900 s; zero timeouts resulted |
| MPI ranks | 1 | MPI is broken on this node — always single-rank |
| Concurrent solves | 6 | |

### Monitoring, resuming, stopping

```bash
# progress
grep -a "^\[optimize\] \(eval\|gen\)" <runs>/grtresna_cmaes/<name>_launch.log | tail -20

# resume a killed run (checkpoint written after every generation)
RESUME=1 RUN_NAME=<same name> bash scripts/campaigns/qball_trajectory/cmaes_run.sh

# stop — orchestrator first, then verify
bash scripts/campaigns/stop_campaign.sh <name>
pgrep -af "grteclyn_wrapper|grtresna"     # must come back empty
```

The resume capability (`--resume` / `RESUME=1`) was added to the driver during this run, in response to
the first attempt having to be restarted from scratch after the solver-timeout failure. It checkpoints
optimizer state, counters and the random stream after every generation and writes atomically, so an
infrastructure failure no longer costs learned search directions. **This run predates the feature and
therefore has no checkpoint file** — it completed, so none was needed.

### Replaying the winner

```bash
python scripts/campaigns/hq/replay_eval.py \
  --trajectory results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-cmaes-refinement/run/trajectory.jsonl \
  --eval 195 --objective-mode f_geo_depth
```

---

## 8. Recommended next steps

1. **Promote the winner at higher resolution and a longer window.** This is the already-agreed next
   action and it answers three open questions at once: whether the emission sweep turns over past t = 18
   (§3.4), whether 48.38 % survives a finer grid given how narrow the ridge is (§3.2), and whether
   branch B's constraint residual stays inside the gate at higher resolution (§3.3). Target: domain
   L = 128, resolution N = 384, single rank, run to t = 64, dense emission sweep, frames enabled, movies
   afterwards.
2. **Promote eval 185 alongside it.** The two genomes differ by 2 % of normalized parameter space but
   sit on different initial-data branches. Running both at high resolution is the cleanest way to
   establish whether the branch-B advantage is physical or an artifact of the coarse-grid solve — and it
   costs one extra run.
3. **Extend the refinement.** The search was still climbing at a steady rate when the budget ran out
   (§3.1). A continuation of 200 more evaluations, warm-started from eval 195 with the step size reduced
   to 0.03, is the highest-expected-value cheap follow-up. With the resume feature now in place it can
   also be run incrementally.
4. **Watch the constraint gate, do not move it.** If gate rejections climb above ~15 % in the
   continuation, shrink the step size. Loosening the 0.03 threshold would break comparability with every
   prior campaign and would let depth be bought with worse-conditioned initial data (§3.6).
