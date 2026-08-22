# The paper campaign — plan and run-book

**Goal.** Produce every number, figure and table in
`research/bondi_dipole/bondi_dipole.tex` from runs on the corrected code path,
so that nothing in the paper rests on the old campaign — which ran with the
double-subtracted potential (canonical stars were collapsing) and the displaced
initial-data metric (stars born off the centre of their own wells). The old
paper's "gap closing", "contact contamination" and "halo bias" storyline was
built on those artefacts; the corrected result is cleaner and stronger: **both
stars accelerate together at a = GM/d², the separation stays put, total
momentum stays zero.**

**Method.** One cell at a time. Each phase has a launch command, a cost, and a
pass gate; nothing in the next phase starts until the gate is green. No
pipeline script — the checking between steps is the point.

The seven artefact rules at the top of `grteclyn-wrapper/README.md` are the
constitution of this campaign; this plan just applies them cell by cell.

---

## 0. Decisions taken up front (so they are not re-litigated per cell)

### Mesh refinement: `max_level = 0` everywhere. Not negotiable, and here is why

1. **The tagger would never fire anyway.** Refinement triggers at a threshold
   of 0.02; the corrected runs peak at |chi − 1| ≈ 0.005. A `max_level=1` cell
   would evolve on level 0 for its entire life, paying AMR bookkeeping for
   nothing. (The old d=8 data *did* tag at t ≈ 47.5 — because the potential bug
   was collapsing the star to chi ≈ 0.43. That physics no longer exists.)
2. **The convergence ladder must mean exactly one thing** — the cell size. A
   tagging criterion re-evaluated per rung refines a different region on each
   rung, and the ladder stops measuring convergence.
3. **The initial-data fix is proven at matched uniform spacing.** The whole
   artefact was the metric arriving displaced when grids disagree; refined
   levels would put interpolation back into the transfer path, and
   `check_gridinit_alignment.py` only certifies level 0.
4. **Memory.** A refined N=128 cell measured ~35–41 GB against 6 GB unigrid;
   refined N=256 would not fit an 80 GB card. Unigrid never comes close:
   6 GB (N=128) → 20 GB (N=192) → 48 GB (N=256), all on one card.
5. **The known RadialRecipe MPI+CUDA segfault is in the AMR fill path** —
   staying unigrid is also what keeps the multi-GPU option alive.

If a referee asks for an AMR cross-check: one `max_level=1` d=10 cell, same
threshold, is cheap and is predicted to be *identical* because the tagger never
fires. Keep it in the back pocket; do not spend it now.

### MPI and the GPUs

- **Elliptic solves: multi-rank MPI, verified in production.** Every archived
  cell already solved on 32 ranks (~20 min at 256³). The node has 128 CPUs:
  up to three 32-rank solves may overlap; do not run four.
- **Evolutions: one GPU per run is the default.** No cell in this matrix needs
  more than 48 of a card's 80 GB, so multi-GPU is a wall-clock option, not an
  OOM necessity. The four cards are best spent running four *cells* at once.
- **Multi-GPU evolution is unverified on RadialRecipe** (the AMR segfault is
  documented; the unigrid path has simply never been tried). Phase 0 smoke-tests
  it. If it passes, the two N=256 rungs may take two cards each
  (`BONDI_EVO_RANKS=2 BONDI_GPU="0,1"`) to halve the longest runs; if it fails,
  nothing in the plan changes.

### Time window

Every quantitative claim is fitted on **t ≤ 200**, the window all archived
cells share (fit from t ≥ 5; before that the gauge is settling). The stability
survey stays at its archived t = 120. Anything running longer is illustration,
never a fit window.

### Movies

Every new cell launches with `GRTECLYN_FRAMES_CACHE_SLICES=1`, so after the
campaign every series can be re-rendered on one fixed colour scale:

```bash
grteclyn-wrapper/scripts/plot/rerender_frames.py <episode>/frames --movies
```

Plotfiles are still consumed and deleted on the fly — the cache keeps only the
2-D slice behind each frame (~1.4 GB per pair cell; delete
`frames/_slice_cache/` once the movies are final).

### Naming and staging

Cell names carry the full grid: `<what>_<dN>_L<box>_N<cells>_lev0`. New runs
land in `runs/bondi/staging/<cell>` and are **moved into
`runs/bondi/runaway_paper/` only after their pass gate is green**, so the
archive never contains an unchecked run. After a cell's alignment and t=0
gates pass, delete its `initial_data.gridinit` (0.5–4.4 GB; regenerable from
the cell's own `launch.sh`).

---

## 1. What already exists and is reused as-is

Paper-ready, in `runs/bondi/runaway_paper/` (details in its README):

| cells | feeds |
|---|---|
| `runaway_pair_d{08,10,12,16}_L64_N128_lev0`, t=200 | separation scan, a·d² = GM table, point-mass limit, headline figures |
| `control_lone_canonical…` / `control_lone_phantom…` (off-centre, x=37), t=200 | single-star null rows of the run matrix |
| `stability/canonical_w{075,080,085,090}…`, t=120 | stability survey section |

Measured anchors these provide: lone-star drift ≤ 1.8e-03 over t=200; four-point
scan gives exponent −2.051 (exact −2) and a·d² → 0.014350; d=10 acceleration
constant across four disjoint fit windows.

**What is missing** is everything below: the same-sign pair nulls, the mirror,
the resolution ladder, the mass-scaling lever, and the wave-zone box.

---

## 2. The phases

Baseline launch environment (identical to the archived d=10 cell; every phase
edits only the lines it names). One block = one cell = one command:

```bash
GRTECLYN_FRAMES_CACHE_SLICES=1 \
BONDI_GPU=<card> BONDI_STOP_TIME=200 \
BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 \
BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="runs/bondi/staging/<cell>" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh
```

The solve grid always obeys the alignment law
`GRTRESNA_N = NFULL × (GRTRESNA_DOMAIN_L / LFULL)` with no solve refinement —
that is the fix that killed the drift artefact, and it is re-verified per cell.

### The universal pass gate (every cell, before its GPU time is trusted)

1. **Alignment** — `research/bondi_dipole/check_gridinit_alignment.py` on the
   fresh `initial_data.gridinit`: metric-minus-matter centroid offset `0.0000`.
2. **Solve exit** — `metadata.json → grtresna_convergence`: Ham% and Mom% at
   or below the cell's tolerance (the d=10 archive landed at 0.00086 / 0.00080).
   Remember the gate binds the phantom side, not the canonical — the phantom
   number is the one to watch.
3. **Matter first** — t=0 row of `small_data/sector_barycenters.dat` matches
   the archived d=10 cell's t=0 row (totals and rms per sector). A dissolved or
   half-painted star makes every later geometry number meaningless.
4. **During the run** — sector total weight flat; scratch dir not growing
   (consumer keeping up).

### Phase 0 — preflight (minutes, no science GPU time)

- [ ] Binary current: rebuild if any source changed since `runaway_paper` was
      produced; `git status` clean, note the commit.
- [ ] `BONDI_DRYRUN=1` on the Phase 1 PP command: printed grid, omegas, solve
      N, tolerance all as intended.
- [ ] Disk: ≥ 100 GB free for staging.
- [ ] **MPI evolution smoke test** (the only new machinery question):

```bash
BONDI_EVO_RANKS=2 BONDI_GPU="0,1" BONDI_STOP_TIME=1 \
BONDI_NFULL=64 BONDI_GRTRESNA_N=128 BONDI_NL_TOL=0.1 BONDI_NL_STALL_TOL=0.002 \
BONDI_GRTRESNA_RANKS=8 BONDI_SCRUTINY=0 \
BONDI_RUNS_DIR="runs/bondi/staging/smoke_mpi_evo" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh
```

  Pass: exit 0, no segfault at the first RK4 advance. Records whether
  `BONDI_EVO_RANKS=2` is available for the N=256 rungs. Either outcome is fine.
  Delete the smoke cell afterwards.

### Phase 1 — the sign matrix (three cells, N=128, ~1.1 h each on the GPU)

The falsifiable core of the paper: only the mixed pair may move.

| cell | change vs baseline | must show |
|---|---|---|
| `control_pair_pp_d10_L64_N128_lev0` | `BONDI_S0=0 BONDI_S1=0`, drop `BONDI_S1_OMEGA` | gap **shrinks** (attraction), pair barycentre still (≤ 5e-03 at t=200) |
| `control_pair_mm_d10_L64_N128_lev0` | `BONDI_S0=1 BONDI_S1=1`, `BONDI_S0_OMEGA=0.7603`, drop `BONDI_S1_OMEGA` | gap **grows** (two negative masses recede), barycentre still |
| `control_mirror_mp_d10_L64_N128_lev0` | `BONDI_S0=1 BONDI_S1=0`, `BONDI_S0_OMEGA=0.7603 BONDI_S1_OMEGA=0.75` | acceleration equal and **opposite** to the archived d=10 cell within a few % |

All three can run at once (three cards, three staggered 32-rank solves — start
them ~5 min apart). Gate for the phase: the two same-sign barycentres sit at
the lone-star noise floor while their gaps move in the *predicted opposite
directions* — attraction where both masses are positive, recession where both
are negative — and the mirror reverses the runaway. Any same-sign cell whose
barycentre drifts like the mixed pair kills the claim; stop and debug, do not
proceed.

### Phase 2 — the resolution ladder (the expensive phase, ~45 GPU-hours)

Headline cell d=10 plus the PP null, each at N=192 and N=256. The tolerance
scales as dx⁴ — this is the one-input fix the old paper could only promise:
every rung previously exited the solve at the *same* residual, so the ladder
had no order to measure. Now the initial violation must fall ~16× per
resolution doubling, and the evolution's own convergence becomes visible.

| rung | grid changes | solve changes | GPU cost (t=200) |
|---|---|---|---|
| N=192 | `BONDI_NFULL=192 BONDI_PLOT_INTERVAL=120` | `BONDI_GRTRESNA_N=384 BONDI_NL_TOL=0.000395 BONDI_NL_STALL_TOL=0.0000079` | ~5.5 h |
| N=256 | `BONDI_NFULL=256 BONDI_PLOT_INTERVAL=160` | `BONDI_GRTRESNA_N=512 BONDI_NL_TOL=0.000125 BONDI_NL_STALL_TOL=0.0000025 BONDI_GRTRESNA_TIMEOUT=43200` | ~17 h (or 2 cards if Phase 0 passed) |

Cells: `runaway_pair_d10_L64_N{192,256}_lev0` (baseline env otherwise) and
`control_pair_pp_d10_L64_N{192,256}_lev0` (Phase 1 PP env otherwise).

Order: run the two N=192 cells first (~5.5 h, overnight is generous), check,
then commit to the two N=256 cells. At N=256 watch the consumer lag — if the
scratch dir grows without bound, stop, double `BONDI_PLOT_INTERVAL`, relaunch.

Gate: (a) solve exits *at* the tighter tolerances — this is new territory for
the solver, and if it stalls above the gate the ladder is off until that is
understood; (b) fitted acceleration a agrees across N=128/192/256 — the spread
is the paper's error bar and ~5% or less is the expectation; (c) the PP
barycentre residual *shrinks* rung to rung; (d) t=0 Ham falls ~16× per rung.
Remember the standing trap: refining the evolution grid alone makes t=0
*worse* — it is the tolerance scaling that pays for the ladder, which is why
the two are locked together in the table above.

### Phase 3 — gravity scales with the source (one cell, ~1.5 h)

The cleanest differential test, and one the grid bug never touched: change the
phantom's mass by ~20%, nothing else. Each star is accelerated by the *other's*
mass, so the canonical star's acceleration must drop by exactly the mass
ratio while the phantom's own acceleration must not move.

First compute the exact retuned mass (no GPU): extend the omega list in
`results/bondi-dipole-runaway/analysis/star_family_scan.py` to cover
0.775–0.80 and pick the phantom omega nearest |M−| ≈ 0.0115 (≈ 0.8 × matched).
The branch anchors: omega 0.750 → −0.015131, 0.775 → −0.013226, 0.7603 →
−0.014350 (the matched point).

Cell `massscale_pair_d10_w<omega>_L64_N128_lev0`: baseline with
`BONDI_S1_OMEGA=<omega>`. Gate: canonical a ratio (this cell / archived d=10) =
mass ratio within a few %; phantom a unchanged within a few %; total momentum
now sums to the expected *non*-zero rate, which is itself a check.

### Phase 4 — the wave zone (one cell, ~9 h)

Doubled box, same cell size, so the light-crossing distance grows without
touching resolution — the wave-zone question is distance, not dx.

Cell `wavezone_pair_d10_L128_N256_lev0`:

```
BONDI_NFULL=256 BONDI_LFULL=128 BONDI_RADII="16 24 32 40" \
BONDI_PSI4_HIGHER_L=1 \
BONDI_GRTRESNA_DOMAIN_L=256 BONDI_GRTRESNA_N=512 BONDI_GRTRESNA_TIMEOUT=43200
```

(the solve box doubles with the evolution box, keeping both the spacing match
and the far outer boundary). 48 GB, one card. Gate: r·psi4 (l=2) amplitude
constant across the four shells; features arrive at retarded time between
shells; the l=2 signal at R=16 consistent with the L=64 runs' extraction.
The dipole-radiation subsection stands or falls here — for a pair with zero
total momentum the mass dipole is static, so the prediction is *no* l=1-type
growth; whatever is measured is reported.

### Phase 5 — optional garnish (only if the paper wants the figures)

- `longrun_pair_d10_t400_L64_N128_lev0` — baseline with `BONDI_STOP_TIME=400`,
  ~2.2 h: the sustained-acceleration money plot (velocity still growing
  linearly at t=400). Illustration only; fits stay on t ≤ 200.
- `control_pair_mm_d10_L64_N{192,256}_lev0` — completes the null-ladder figure
  with MM alongside PP. Same cost as the PP rungs; skip unless the convergence
  figure looks thin without it.

### Phase 6 — analysis, movies, packing (no GPU)

1. Update the hardcoded cell lists in `results/bondi-dipole-runaway/analysis/`
   (`separation_scaling.py`, `convergence_check.py`, `momentum_balance.py`,
   `newtonian_reference.py`, `make_tables.py`, …) from the old `convA_*` names
   to the `runaway_paper` names, then regenerate every table.
2. Re-render every movie on its fixed colour scale (`rerender_frames.py
   <episode>/frames --movies`), then delete each `frames/_slice_cache/`.
3. `research/bondi_dipole/pack_runaway.sh` to refresh the git-tracked extract
   pack in `results/bondi-dipole-runaway/campaign/`; update both READMEs with
   the new cells.
4. Rewrite the results sections of `bondi_dipole.tex` from the new tables —
   the old text's gap-closing/contact-contamination narrative does not survive
   the fix and must not be patched around.

---

## 3. Budget

| phase | cells | GPU-hours | wall (4 cards) |
|---|---|---|---|
| 1 sign matrix | 3 | ~3.5 | ~1.5 h |
| 2 ladder | 4 | ~45 | ~1.5 days |
| 3 mass scale | 1 | ~1.5 | with phase 2 |
| 4 wave zone | 1 | ~9 | with phase 2 |
| 5 optional | 0–3 | 0–25 | — |
| **total (required)** | **9** | **~59** | **~2–3 days** |

Solves add ~20 min (256³) to ~4 h (512³) of CPU time per cell, overlapping GPU
work on other cells. The entire required matrix is under three days — the cost
of *not* checking between phases was three weeks of artefact archaeology, which
is the ratio to remember when tempted to launch everything at once.
