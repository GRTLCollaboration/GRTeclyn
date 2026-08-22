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

## The scoreboard — every run the paper needs

Tick a box only when the cell has passed its gate and been moved into
`runs/bondi/runaway_paper/`. Details for each are in the phase it belongs to.

**Already in the archive (corrected path, gates passed):**

- [x] `runaway_pair_d08_L64_N128_lev0` — separation scan point, strong-signal pictures
- [x] `runaway_pair_d10_L64_N128_lev0` — the headline cell; base rung of the ladder
- [x] `runaway_pair_d12_L64_N128_lev0` — separation scan point
- [x] `runaway_pair_d16_L64_N128_lev0` — separation scan point
- [x] `control_lone_canonical_L64_N128_lev0` — single-star null (box centre)
- [x] `control_lone_phantom_L64_N128_lev0` — single-star null (off-centre, the sharper test)
- [x] `stability/canonical_w{075,080,085,090}_L64_N128_lev0` — stability survey, t=120

**Required — the paper is not submittable without these:**

- [ ] `smoke_mpi_evo` — phase 0, minutes; answers whether 2-GPU evolution works (deleted after)
- [ ] `control_pair_pp_d10_L64_N128_lev0` — phase 1, ~1.1 h; two canonical stars: gap shrinks, barycentre still
- [ ] `control_pair_mm_d10_L64_N128_lev0` — phase 1, ~1.1 h; two phantom stars: gap grows, barycentre still
- [ ] `control_mirror_mp_d10_L64_N128_lev0` — phase 1, ~1.1 h; sectors swapped: runaway reverses
- [ ] `runaway_pair_d10_L64_N192_lev0` — phase 2, ~5.5 h; middle rung of the ladder
- [ ] `runaway_pair_d10_L64_N256_lev0` — phase 2, ~17 h; finest rung of the ladder
- [ ] `control_pair_pp_d10_L64_N192_lev0` — phase 2, ~5.5 h; null residual must shrink with the grid
- [ ] `control_pair_pp_d10_L64_N256_lev0` — phase 2, ~17 h; null residual, finest rung
- [ ] `massscale_pair_d10_w<omega>_L64_N128_lev0` — phase 3, ~1.5 h; lighter phantom: pull scales with the source
- [ ] `wavezone_pair_d10_L128_N256_lev0` — phase 4, ~9 h; doubled box, four extraction shells

**Optional — only if the paper wants the figure:**

- [ ] `longrun_pair_d10_t400_L64_N128_lev0` — ~2.2 h; sustained-acceleration money plot
- [ ] `control_pair_mm_d10_L64_N192_lev0` — ~5.5 h; MM alongside PP in the null ladder
- [ ] `control_pair_mm_d10_L64_N256_lev0` — ~17 h; MM null, finest rung
- [ ] `amrcheck_pair_d10_L64_N128_lev1` — ~1.5 h; referee-proofing only (predicted identical to lev0)

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

All four cards (0–3) stay busy throughout; the gates order the *families*,
not the queue:

| wave | GPU 0 | GPU 1 | GPU 2 | GPU 3 |
|---|---|---|---|---|
| smoke | 2-rank smoke test (cards 0,1) | ← | — | — |
| A | PP d10 N128 | MM d10 N128 | MP mirror N128 | mass-scale (after the CPU family scan) |
| B — phase-1 gate green | PM N192 | PP N192 | wave-zone L128 | long-run t400 (optional) |
| C — N=192 gate green | PM N256 | (+0 if 2-rank passed) | PP N256 | (+2 if 2-rank passed) |

### Time window

Every quantitative claim is fitted on **t ≤ 200**, the window all archived
cells share (fit from t ≥ 5; before that the gauge is settling). The stability
survey stays at its archived t = 120. Anything running longer is illustration,
never a fit window.

### Frames and movies — only where a figure needs them

Drawing frames is the dominant per-plotfile CPU cost, and most cells' product
is numbers, not pictures. The matrix default is therefore `GRTECLYN_FRAMES=0`
(metrics-only consumer; plotfiles are still deleted on the fly — that gate
does not depend on frames). Frames are switched on only for the cells a
figure or movie actually comes from: the P/M singles (archived, already have
frames) and the finest-grid PM and PP cells. Those launch with
`GRTECLYN_FRAMES_CACHE_SLICES=1`, so afterwards every kept series can be
re-rendered on one fixed colour scale:

```bash
grteclyn-wrapper/scripts/plot/rerender_frames.py <episode>/frames --movies
```

Plotfiles are still consumed and deleted on the fly — the cache keeps only the
2-D slice behind each frame (~1.4 GB per pair cell; delete
`frames/_slice_cache/` once the movies are final).

### Where the GW shells may sit — outside the stars, inside the sponge

Measured at t=0 on the archived d=10 cell: star cores at ±5 from the box
centre with rms radius 4.34, so ~90% of the matter lives inside r ≈ 11.5 and
the tails reach further. The template's inner extraction shell at r = 8 passes
straight through both stars — it is a near-zone monitor, never a GW
measurement. And the in-code Weyl extraction has been OFF in every cell so
far: all psi4 to date is the consumer's coarse per-plotfile shell sampling.

| box | matter ends | sponge (inner→outer) | shells that mean radiation |
|---|---|---|---|
| L=64 | ~11.5 + tails | 24→32 | r = 16 only — quote nothing else |
| L=128 wave zone | ~11.5 + tails | **48→60, must be moved out** | 16, 24, 32, 40 |

The sponge radii are absolute numbers, not box-scaled: left at their 24→32
defaults, the doubled box would run its dissipation band straight through
three of the four planned shells. Phase 4 sets them explicitly, and switches
the dense in-code extraction on with the new `BONDI_EXTRACTION_RADII` knob.

### Naming and staging

Cell names follow the archive's pattern exactly —
`runaway_pair_d08_L64_N128_lev0` reads as: what the run is, star separation,
box side, cells per side, refinement depth. Every new cell keeps that shape:
`<what>_<dN>_L<box>_N<cells>_lev0`, with extra qualifiers (a retuned omega, a
longer stop time) slotted in before the grid part, e.g.
`massscale_pair_d10_w0790_L64_N128_lev0`, `longrun_pair_d10_t400_L64_N128_lev0`.
A name alone must tell you what the folder holds and on which grid. New runs
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
GRTECLYN_FRAMES=0 \
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
`control_pair_pp_d10_L64_N{192,256}_lev0` (Phase 1 PP env otherwise). The two
N=256 cells are the figure cells: they replace `GRTECLYN_FRAMES=0` with
`GRTECLYN_FRAMES_CACHE_SLICES=1`; every other new cell in the matrix stays
frameless.

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
BONDI_EXTRACTION_RADII="16 24 32 40" BONDI_PSI4_HIGHER_L=1 \
BONDI_SPONGE_INNER=48 BONDI_SPONGE_OUTER=60 \
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
  linearly at t=400). It exists for the picture, so it launches with frames
  and the slice cache on. Illustration only; fits stay on t ≤ 200.
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

---

## 4. Outlook — how far toward the speed of light can the dipole be pushed?

Not all the way: a constantly self-accelerating pair follows hyperbolic
motion, approaching c asymptotically, never reaching it. The real question is
what fraction of c fits in a computable box. With the measured d=8
acceleration (2.32e-4) and the pair racing down a long thin box (transverse
side kept at 64, cell 0.5, pair launched near the −x face):

| target | run length t | distance | box length | memory | GPU-days |
|---|---|---|---|---|---|
| 0.10 c | 430 | 22 | 160 | 15 GB | 0.2 |
| 0.30 c | 1,360 | 210 | 352 | 33 GB | 1.7 |
| 0.50 c | 2,500 | 670 | 800 | 75 GB | 7 |
| 0.70 c | 4,200 | 1,700 | 1,856 | 174 GB — 3+ cards | 28 |
| 0.90 c | 8,900 | 5,600 | 5,728 | 537 GB — over the node | 180 |

0.3c costs under two days on one card and already carries a 4.8% deviation
from the Newtonian v = a·t — a measurable test of relativistic
self-acceleration. 0.5c (15% deviation) is a week and just fits one card.
Beyond that the box, not the time, is the wall: the distance needed to go
relativistic is d²/GM ≈ 900 star radii, which is why literal c is petabytes.

Ways to cut the cost, strongest first:

1. **Heavier stars — cost falls as 1/M².** Time, distance, box and step count
   all scale as 1/M, so GPU time to a fixed speed scales as 1/M². The
   ω=0.75 star (M = 0.01435) sits far below the stable branch's mass maximum
   (the branch runs from ω ≈ 0.67 up). Run the family scan (CPU-only) over
   ω = 0.67–0.75 and take the heaviest stable canonical that has a
   mass-matched phantom partner: 3× the mass makes every row above 9×
   cheaper — 0.5c in under a day.
2. **A recentring (comoving) box.** The box is long only because the pair
   moves. Shifting every field back by a whole number of cells at intervals
   is an exact copy — no interpolation, the same reasoning that fixed the
   initial data — with the sponge eating the seam at the trailing face. That
   caps the box at ~L=128 forever and makes cost linear in run time: 0.9c in
   a small box becomes days, not half a year. New default-off feature; needs
   one validation cell (recentred vs not, over a window both can reach).
3. **Both together**: with 3× mass, 0.9c needs only t ≈ 3,000 in a small
   static box — even paying for the finer grid the Lorentz-contracted star
   eventually demands (γ = 2.3 at 0.9c halves the star's grid footprint;
   plan on dx = 0.25 past ~0.6c), it stays in single-digit GPU-days.
4. Smaller levers: separation below 8 buys ~30% but enters envelope overlap;
   an initial boost via the solver's boost knob skips the slow early phase,
   but that code path has never been validated — it gets a validation cell
   before it gets believed.

None of this is in the required matrix; it is the natural follow-up paper.
