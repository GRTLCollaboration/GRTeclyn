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

The eight artefact rules at the top of `grteclyn-wrapper/README.md` are the
constitution of this campaign; this plan just applies them cell by cell.

---

## The scoreboard — every run the paper needs

Tick a box only when the cell has passed its gate and been moved into
`runs/bondi/staging/archive/`. Details for each are in the phase it belongs to.

**Already in the archive (corrected path, gates passed):**

- [x] `runaway_pair_d08_L64_N128_lev0` — separation scan point, strong-signal pictures
- [x] `runaway_pair_d10_L64_N128_lev0` — the headline cell; base rung of the ladder
- [x] `runaway_pair_d12_L64_N128_lev0` — separation scan point
- [x] `runaway_pair_d16_L64_N128_lev0` — separation scan point
- [x] `control_lone_canonical_L64_N128_lev0` — single-star null (box centre)
- [x] `control_lone_phantom_L64_N128_lev0` — single-star null (off-centre, the sharper test)
- [x] `stability/canonical_w{075,080,085,090}_L64_N128_lev0` — stability survey, t=120

**Required — the paper is not submittable without these:**

- [x] `smoke_mpi_evo` — phase 0, **passed 2026-08-22**: 50/50 steps on cards 0+1, no segfault, exit 0. Two-GPU evolution is available at `max_level=0`. Cell deleted.
- [x] `control_pair_pp_d10_L64_N128_lev0` — phase 1, **evolved to t=200 on 2026-08-22**; barycentre 32.00073, drift +0.00073, min χ = 0.97947. Gap unmeasurable (see phase 1)
- [x] `control_pair_mm_d10_L64_N128_lev0` — phase 1, **evolved to t=200 on 2026-08-22**; barycentre 31.99972, drift −0.00028, min χ = 1.00000 (exactly flat). Gap unmeasurable (see phase 1)
- [x] `control_mirror_mp_d10_L64_N128_lev0` — phase 1, **evolved to t=200 on 2026-08-22**; runaway reverses exactly: displacement and acceleration both −1.00002× the archived cell. Frameless by decision
- [ ] `runaway_pair_d10_L64_N192_lev0` — phase 2, ~5.5 h; middle rung of the ladder
- [ ] `runaway_pair_d10_L64_N256_lev0` — phase 2, ~17 h; finest rung of the ladder — **frames + slice cache** (headline movie)
- [ ] `control_pair_pp_d10_L64_N192_lev0` — phase 2, ~5.5 h; **solve gate held flat at
  0.002 (not the dx⁴ value) — the same-sign floor, see Phase 2**; the *evolution* null
  residual must still shrink with the grid, the t=0 violation will not
- [ ] `control_pair_pp_d10_L64_N256_lev0` — phase 2, ~17 h; null residual, finest rung — **frames + slice cache** (the null's movie)
- [x] `massscale_pair_d10_w0804_L64_N128_lev0` — phase 3, **evolved to t=200 on 2026-08-22**; lighter phantom (M = −0.011472, 79.95% of matched). Drift +2.6732, fitted a = 1.3365e-04 = 0.928× the matched cell's 1.4404e-04, so the pull does scale with the source. The pair no longer moves rigidly — the two stars pull each other unequally and the separation closes 10.000 → 9.408 — so the exact ratio needs the d-corrected fit, not the raw number
- [ ] `wavezone_pair_d10_L128_N256_lev0` — phase 4, ~9 h; doubled box, four extraction shells

**Optional — only if the paper wants the figure:**

- [ ] `longrun_pair_d10_t400_L64_N128_lev0` — ~2.2 h; sustained-acceleration money plot — **frames + slice cache**
- [ ] `control_pair_mm_d10_L64_N192_lev0` — ~5.5 h; MM alongside PP in the null ladder
- [ ] `control_pair_mm_d10_L64_N256_lev0` — ~17 h; MM null, finest rung
- [ ] `amrcheck_pair_d10_L64_N128_lev1` — ~1.5 h; referee-proofing only (predicted identical to lev0)
- [ ] `chase_pair_d08_v03c_Lx352_L64_N128_lev0` — ~1.7 GPU-days; ride the runaway to 0.3c in a long box (section 4). **Superseded on paper by the recentring box — ~7.5 h for the same physics; see section 5 for the build plan.** Follow-up paper material either way

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
does not depend on frames). Frames are switched on only for the cells a figure
or movie actually comes from, and those launch with
`GRTECLYN_FRAMES_CACHE_SLICES=1` so every kept series can be re-rendered on one
fixed colour scale afterwards.

**The cells that draw, and the cells that do not** (settled 2026-08-22):

| draws frames | why |
|---|---|
| `runaway_pair_d10_L64_N256_lev0` | the headline movie, at the best resolution the campaign has |
| `control_pair_pp_d10_L64_N256_lev0` | the null beside it, same grid and same colour scale |
| `longrun_pair_d10_t400_L64_N128_lev0` | the sustained-acceleration money plot; it exists for the picture |
| P/M singles | archived, already have frames |

| frameless | why |
|---|---|
| `control_mirror_mp_d10_L64_N128_lev0` | **decided while it was already running.** Its product is the reversed acceleration, a number, and it is the mirror of a cell that already has a movie |
| `control_pair_pp/mm_d10_L64_N128_lev0`, `massscale_…` | numbers only |
| both N=192 rungs, `wavezone_…` | ladder and extraction cells; nothing is read off a picture |

**Frames are ON unless a cell switches them off.** `consumer_frames_enabled()`
returns true for anything that is not literally `GRTECLYN_FRAMES=0/off/no/false`
— an unset variable draws frames. So a frameless cell needs the flag written out
explicitly, while a frames cell works either way. Every frames cell in this
campaign nevertheless states `GRTECLYN_FRAMES=1` beside its slice cache, because
"frames on" should be visible in the launcher rather than inferred from a
default: verified across all ten launchers on 2026-08-22 before wave B started.

**This cannot be revisited after the fact** (README rule 6): the plotfiles a
frame is drawn from are deleted as they are consumed, so a frameless cell can
never be re-rendered — it can only be re-run. A movie of the mirror cell would
therefore cost a fresh ~1.1 h cell, which is cheap if the paper turns out to
want it and wasted if it does not.

The re-render, once a frame-bearing cell lands:

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
`runs/bondi/staging/archive/` only after their pass gate is green**, so the
archive never contains an unchecked run. After a cell's alignment and t=0
gates pass, delete its `initial_data.gridinit` (0.5–4.4 GB; regenerable from
the cell's own `launch.sh`).

### Running order — launched by hand, one solve at a time

**Cells are started manually, one command each. Do not write a queue script
that launches them by itself.** An orchestrator outlives the session that made
it: on 2026-08-22 a leftover queue from an earlier session was found orphaned
to init, still launching cells on its own, two 32-rank solves at once, while a
replacement queue was being started. Nothing in the run directories shows that
is happening. Before any launch, look for leftovers and kill them first:

```bash
python3 grteclyn-wrapper/scripts/ops/sweep_ranks.py            # what is alive
python3 grteclyn-wrapper/scripts/ops/sweep_ranks.py --kill solves
```

`sweep_ranks.py` walks `/proc` directly because `ps`, `top`, `pgrep` and `free`
are all broken on this node. Kill any orchestrator **first**, then the workers,
then re-run the bare command to verify — killing workers alone just makes an
orchestrator advance to the next cell.

Four cards, but only **one 32-rank elliptic solve may run at any moment**
(README rule 10: a second solve starves every evolution in flight down to a
fifth speed). So each launch waits for the previous cell's solve to hand over
to the GPU — that wait is a decision to take, not something to automate.

Wave B, in the intended order with its pinned card:

| # | cell | card | ~GPU-hours |
|---|---|---|---|
| 1 | `runaway_pair_d10_L64_N192_lev0` | 1 | 5.5 |
| 2 | `wavezone_pair_d10_L128_N256_lev0` | 2 | 9 |
| 3 | `control_pair_pp_d10_L64_N192_lev0` | 3 | 5.5 |
| 4 | `runaway_pair_d10_L64_N256_lev0` | 0 | 17, frames |
| 5 | `control_pair_pp_d10_L64_N256_lev0` | 1 | 17, frames |
| 6 | `longrun_pair_d10_t400_L64_N128_lev0` | 3 | 2.2, frames |

A relaunch needs the cell's half-built episode directory deleted first, or the
wrapper refuses to start with "already exists" and exits within a second. Check
for that before assuming a launch took.

Wave B went up 2026-08-22 03:2x–03:32 on all four cards. The last two wait for
a card and are launched by hand when one frees:

| when this finishes | the card it frees | launch this |
|---|---|---|
| `runaway_pair_d10_L64_N192_lev0` | 1 | `control_pair_pp_d10_L64_N256_lev0` |
| `control_pair_pp_d10_L64_N192_lev0` | 3 | `longrun_pair_d10_t400_L64_N128_lev0` |

```bash
# card 1, after runaway N192 lands  (~4 h solve + ~17.6 h evolve, frames on)
bash runs/bondi/staging/control_pair_pp_d10_L64_N256_lev0/launch.sh

# card 3, after PP N192 lands       (~0.4 h solve + ~2.2 h evolve, frames on)
bash runs/bondi/staging/longrun_pair_d10_t400_L64_N128_lev0/launch.sh
```

Run them from the repository root — the launchers take that path as given, and
a stale working directory is the one way they fail instantly with
"No such file or directory". Preflight both first, and check no solve is
already running before starting a second one.

**Projected finish of the whole stack: Sun 23 Aug, 08:00 at best, ~15:00 with
contention margin.** The critical path is `control_pair_pp_d10_L64_N256_lev0`,
which cannot start until a card frees at ~Sat 10:20 and then needs 21.6 h of
its own; everything else is done by Sun 01:10. The only way to pull the stack
in is to give that cell a card sooner.

---

## 1. What already exists and is reused as-is

Paper-ready, in `runs/bondi/staging/archive/` (details in its README):

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
2. **Solve exit** — the tolerance is a request, not a guarantee (rule 8): the
   solver can also leave by stalling or by running out of its 50 iterations,
   and it says "converged" either way. Read the verdict, not the label:

   ```bash
   grteclyn-wrapper/.venv/bin/python research/bondi_dipole/check_solve_exit.py \
       runs/bondi/staging/<cell>
   ```

   Exit 0 or the cell does not count. The d=10 archive landed at 0.00086 /
   0.00080 against 0.002 — only ~2x of headroom, so this is not a formality.
   Remember the gate binds the phantom side, not the canonical.
3. **Matter first** — t=0 row of `small_data/sector_barycenters.dat` matches
   the archived d=10 cell's t=0 row (totals and rms per sector). A dissolved or
   half-painted star makes every later geometry number meaningless.
4. **During the run** — sector total weight flat; scratch dir not growing
   (consumer keeping up).
5. **After the run** — write the measured wall time into the ledger in section
   3. Every episode records it, so this is a lookup, not a stopwatch:

   ```bash
   python3 -c "import json,sys;print(json.load(open(sys.argv[1]))['simulation_elapsed_seconds']/3600)" \
       runs/bondi/staging/<cell>/*/metadata.json
   ```

   The estimates in this plan came from the archive; the ledger is what the
   next campaign should cost its runs from.

### Phase 0 — preflight (minutes, no science GPU time)

- [ ] Binary current: rebuild if any source changed since the archive was
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

  **Result (2026-08-22): passed.** The solve converged in 8 nonlinear
  iterations (Ham 0.082%, Mom 0.086%), the evolution ran all 50 steps split
  evenly across the two cards (identical FAB footprint per rank), and the run
  exited 0 with nothing in the log resembling a fault. So the AMR segfault is
  confined to the refined fill path, exactly as suspected, and the two N=256
  rungs may use `BONDI_EVO_RANKS=2 BONDI_GPU="0,1"` to halve their wall time.

### Phase 1 — the sign matrix (three cells, N=128, ~1.1 h each on the GPU)

### Result — the mirror reverses, and the same-sign pairs do not move

Re-measured 2026-08-22 from `sector_barycenters.dat`, **after both cells
finished and their consumers drained**. "Drift" is the mean of the two sector
barycentres relative to birth; the acceleration is the quadratic coefficient of
a fit over the full shared window `t ≤ 200`. An earlier version of this table
quoted a mid-run snapshot read off the core tracker at `t = 147` and is
superseded — the numbers below are the ones that reproduce from the archived
files.

| | archive (canonical left, phantom right) | mirror (sectors swapped) | ratio |
|---|---|---|---|
| drift at t=200 | **+2.88147** | **−2.88153** | **−1.000022** |
| acceleration (fit, t ≤ 200) | **+1.44877e−04** | **−1.44878e−04** | **−1.000010** |
| worst disagreement over the whole run | — | — | **2.2e−05** relative |

Swapping which star carries positive mass and which carries negative mass
flips the direction of travel and changes nothing else — same speed, same
acceleration, same wobble, agreeing to two parts in a hundred thousand. That is the
point of the control: a drift produced by the grid, the boundary or the solver
would have kept pointing the same way when the physics was swapped. The
direction follows the matter, not the machine.

The same-sign nulls, same date, same files, same window:

| | MM (two phantoms) | PP (two canonicals) | mirror (mixed) |
|---|---|---|---|
| drift of the pair, birth → t=200 | **−0.00028** | **+0.00073** | −2.88153 |
| sector barycentre at t=200 | 31.99972 | 32.00073 | moves 2.9 |
| peak field amplitude, birth → t=200 | 0.02443 → 0.02435 | 0.02457 → 0.02419 | 0.02452 → 0.02443 |
| min χ over the run (1 = flat, 0 = horizon) | **1.00000** | 0.97947 | 0.98889 |
| sector field weight, birth → t=200 | 8.06 → 31.53 (3.9×) | 7.83 → 32.37 (4.1×) | 3.91 → 4.73 (1.2×) |

**Four** orders of magnitude separate the mixed pair from either same-sign pair
— 2.88 against 3e−04 and 7e−04. MM sits at exactly flat geometry and does not
move at all: two negative masses cancel each other's pull to nothing and there
is no dipole to drive. Nothing in the phase is anywhere near collapse — the
lowest χ on the board is 0.979, and a horizon needs it near zero.

The last row is the one to keep an eye on. Both same-sign cells **quadruple**
their sector field weight by t = 200 while the mixed pair grows only 1.2×.
That is the same-sign halo spreading, already documented — the cores stay put
and steady, the envelope does not. It is why the same-sign nulls are quoted on
the barycentre and the core, and why nothing from these two cells is quoted
past t ≈ 60 without saying so.

The phase gate is met on the barycentre and the mirror. The gap half of it is
not measurable with the present tracker; see below.

**Measured 2026-08-22, and it is the documented behaviour**
(`MatterDebugg.md`): the two same-sign cells grow a large halo at late times —
tracked activity up ~6x, rms radius reaching the domain, confined fraction
falling to ~0.10 — while their cores survive intact (peak amplitude steady to
within 3%). This is not a fault. In a same-sign cell both lumps live in *one*
field, so the potential's cross terms give them a direct scalar interaction on
top of gravity; the mixed pair's two sectors share no potential coupling and
meet only through gravity, which is exactly why it is the clean case. The
sponge cannot help: it is zero inside its inner radius, where essentially all
of the halo sits, and it is Kreiss-Oliger dissipation tuned to grid-scale
noise rather than an absorbing layer. **Do not quote a same-sign cell past
t ~ 60.**

The barycentre null is unaffected and is the number this phase exists for: at
t=80 both same-sign pairs sit within 3e-04 of their birth position, six times
below the lone-star noise floor and 1500x below the mixed pair's 0.467 over
the same window.

One measurement this phase does *not* deliver: the **gap** between two
same-sign stars. Both live in one sector, so the core tracker locks onto the
interference peak at their midpoint and reports the separation as `nan` on
every row. The attraction/recession half of the gate needs a two-lump-aware
tracker, which does not exist yet — the barycentre half stands on its own.


The falsifiable core of the paper: only the mixed pair may move.

| cell | change vs baseline | must show |
|---|---|---|
| `control_pair_pp_d10_L64_N128_lev0` | `BONDI_S0=0 BONDI_S1=0`, drop `BONDI_S1_OMEGA`, **`BONDI_GRTRESNA_MAXIMAL_SLICING=1`** (see below) | gap **shrinks** (attraction), pair barycentre still (≤ 5e-03 at t=200) |
| `control_pair_mm_d10_L64_N128_lev0` | `BONDI_S0=1 BONDI_S1=1`, `BONDI_S0_OMEGA=0.7603`, drop `BONDI_S1_OMEGA` | gap **grows** (two negative masses recede), barycentre still |
| `control_mirror_mp_d10_L64_N128_lev0` | `BONDI_S0=1 BONDI_S1=0`, `BONDI_S0_OMEGA=0.7603 BONDI_S1_OMEGA=0.75` | acceleration equal and **opposite** to the archived d=10 cell within a few % |

**The all-canonical cell needs one extra knob, and it is the easiest thing in
this campaign to get wrong.** The solver turns maximal slicing on *by itself*
whenever any lump carries negative energy — the CTTK ansatz
`K = sign·sqrt(24πGρ)` is imaginary for `ρ < 0`, so it has no choice. Every
phantom-bearing cell therefore starts from `K = 0`. PP is the only cell in the
whole matrix with no phantom star, so left alone it starts on an
already-collapsing slice and is not comparable to the very cell it is the null
for. Forcing the flag also picks up the rest of the matched path (psi
relaxation 0.6, psi floor 0.1, arithmetic coefficient averaging), so the two
cells then differ only in the sign of the mass.

This was caught on 2026-08-22 by comparing the solve records of the four wave-A
cells: PP was the only one showing `maximal_slicing=0`, and its residuals
oscillated over four orders of magnitude between iterations instead of falling
by a clean factor of ~2.5 like every other cell. The run was stopped before it
reached the GPU and relaunched with the flag. **Every same-sign-canonical cell
in this plan carries it — including the PP ladder rungs in phase 2.**

**The archive was audited for the same fault and is clean** — nothing there
needs re-running. The lone-canonical cell, the only other all-canonical cell
ever produced, set the flag explicitly in its own `launch.sh`. The stability
survey reused pre-solved `K = 0` data. And the decisive check is measured
rather than inferred: the trace of the extrinsic curvature at the innermost
shell reads exactly `0.00000` at birth in every one of the ten archived cells,
which is only true of a flat start. A cell born on the CTTK slice would carry
`K` of order `1e-01` there and show the lapse swinging within the first few
time units.

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

The tolerances below are the **mixed-cell** values. The same-sign rungs cannot
use them — see "The same-sign rungs keep the N=128 gate" below.

| rung | grid changes | solve changes | GPU cost (t=200) |
|---|---|---|---|
| N=192 | `BONDI_NFULL=192 BONDI_PLOT_INTERVAL=120` | `BONDI_GRTRESNA_N=384 BONDI_NL_TOL=0.000395 BONDI_NL_STALL_TOL=0.0000079` | ~5.5 h |
| N=256 | `BONDI_NFULL=256 BONDI_PLOT_INTERVAL=160` | `BONDI_GRTRESNA_N=512 BONDI_NL_TOL=0.000125 BONDI_NL_STALL_TOL=0.0000025 BONDI_GRTRESNA_TIMEOUT=43200` | ~17 h (or 2 cards if Phase 0 passed) |

Cells: `runaway_pair_d10_L64_N{192,256}_lev0` (baseline env otherwise) and
`control_pair_pp_d10_L64_N{192,256}_lev0` (Phase 1 PP env otherwise — which
means they carry `BONDI_GRTRESNA_MAXIMAL_SLICING=1` too; a PP rung without it
is not a null for anything). The two
N=256 cells are the figure cells: they replace `GRTECLYN_FRAMES=0` with
`GRTECLYN_FRAMES_CACHE_SLICES=1`; every other new cell in the matrix stays
frameless.

Order: run the two N=192 cells first (~5.5 h, overnight is generous), check,
then commit to the two N=256 cells. At N=256 watch the consumer lag — if the
scratch dir grows without bound, stop, double `BONDI_PLOT_INTERVAL`, relaunch.

Gate: (a) the **mixed** rungs exit *at* the tighter tolerances (measured
2026-08-22: N=192 converged at step 14, Ham 3.52e-04 %, still improving 58% a
step — the dx⁴ ladder is sound for these cells); (b) fitted acceleration a
agrees across N=128/192/256 — the spread is the paper's error bar and ~5% or
less is the expectation; (c) the PP barycentre residual *shrinks* rung to rung,
but any residual that does *not* shrink must be checked against the same-sign
floor below before it is called physical; (d) t=0 Ham falls ~16× per rung —
**mixed rungs only**; the PP rungs are flat across resolution by construction.
Remember the standing trap: refining the evolution grid alone makes t=0
*worse* — it is the tolerance scaling that pays for the ladder, which is why
the two are locked together in the table above.

#### The same-sign rungs keep the N=128 gate — measured 2026-08-22

The dx⁴ ladder is right for the mixed cells and **unreachable** for the
same-sign ones. `control_pair_pp_d10_L64_N192_lev0` was launched at
`BONDI_NL_TOL=0.000395`, ran 16 nonlinear steps, and stopped improving:

| step | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|
| Ham [%] | 9.55e-04 | 6.27e-04 | 5.57e-04 | 5.46e-04 | 5.44e-04 |
| improvement on the step | 53% | 34% | 11% | 2.1% | 0.4% |

It floors at 5.4e-04 %, *above* the 3.95e-04 % it was asked for. The cell was
killed at step 16; left alone it would have ground through all 50 steps and
exited by `max_NL_iterations` — the door that does not count — because the
stall detector requires **both** residuals to flatten and Mom was still falling
59% a step (7.02e-05 at step 15, 2.86e-05 at step 16). Note the trap: a stall
tolerance tuned to the Hamiltonian will never fire while the momentum residual
is healthy.

**Cause, confirmed in the solver source.** `BoundaryConditions::fill_constraint_box`
fills the psi correction on the outer boundary with a hard constant 0.0 (the
comment reads "zero for psi"), and `CTTKHybrid::initialise_multigrid_vars` sets
`c_psi_reg` to the constant `regularised_part_psi`. Together these pin the
conformal factor at the solve boundary to exactly 1 — flat space — for the whole
solve. That is the correct condition only for a configuration with zero net ADM
mass. A same-sign pair carries 2M, whose true conformal factor at the boundary is
1 + M/(2R); the solve is structurally unable to represent it. The resulting error
is set by the box, not by the grid, so refinement cannot remove it.

The floor tracks net mass, not resolution — this is the discriminating evidence:

| cell | net ADM mass | lowest Ham reached | still improving? |
|---|---|---|---|
| `runaway_pair_*` (mixed) | 0 | 3.52e-04 % | yes, 58%/step |
| `control_mirror_mp_*` (mixed) | 0 | 8.62e-04 % | yes, 59%/step |
| `massscale_pair_*` (light phantom) | +0.0029 | 8.58e-04 % | yes, 59%/step |
| `control_pair_pp_*` N=192 | +0.0287 | 5.44e-04 % | **no — floored** |

Ten times less leftover mass and no floor appears at all; the same leftover mass
at two resolutions gives the same floor. Note also that `control_pair_pp_*` and
`control_pair_mm_*` at N=128 passed only because their 0.002 gate caught them
(1.27e-03 and 1.33e-03) before the floor was revealed — one step earlier both
were still above it. The N=128 pass was luck, not headroom.

**Fixes considered and rejected.** Moving the boundary out shrinks the error as
1/R while the solve cost grows as R³: merely reaching the dx⁴ gate at N=192 needs
R×1.4 (2.6× cost, zero margin), and a factor-3 margin needs R×2.7 (~20× cost, a
~630 GB intermediate). Replacing the hard zero with a 1/r falloff for the
correction is the physically correct fix and costs nothing at runtime, but it is
new code in the shared elliptic solver and would have to be revalidated before it
could produce paper data.

**Decision: accept the floor and say so.** The same-sign rungs run at a flat
`BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004` at *every* resolution — the same
gate the N=128 rung already used — because their achievable accuracy is set by
the solve boundary and does not improve with cell size. All three PP rungs
therefore begin from initial data of the same quality (~1e-03 % Ham). That is the
honest description of what this solver delivers for a net-mass configuration, and
it is what the paper should state.

What this costs, and what it does not:

| claim | affected? |
|---|---|
| PP/MM as null controls (no self-acceleration) | **no** — the boundary error is spherically symmetric and cannot manufacture a dipole |
| the headline runaway ladder at N=128/192/256 | **no** — mixed pairs have zero net mass, the flat boundary is *correct* for them, and they show no floor |
| PP initial violation falling with resolution | **yes** — it does not fall; do not quote a convergence order for the PP rungs |
| PP evolution convergence rung to rung | partly — the rungs differ in evolution grid only, initial data held at a common quality, which is the cleanest reading available here |

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

**Done 2026-08-22 (CPU, no GPU): the frequency is omega = 0.8040**, giving
M = −0.011472, which is 79.95% of the matched phantom mass. Neighbours for
reference: 0.7900 → −0.012267, 0.8000 → −0.011692, 0.8100 → −0.011162.

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
   to the archive cell names, then regenerate every table.
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

### Measured wall time — the ledger

Filled in from each cell's own record after it lands (gate step 5), so the next
campaign costs its runs from measurements rather than from the estimates above.
The last column is the useful one: GPU-hours per 1000 units of evolution time,
which is what actually transfers between cells.

| cell | N | t | cards | GPU-hours | h / 1000 t |
|---|---|---|---|---|---|
| `control_lone_canonical_L64_N128_lev0` | 128 | 200 | 1 | **1.09** | 5.5 |
| `control_lone_phantom_L64_N128_lev0` | 128 | 200 | 1 | **1.10** | 5.5 |
| `runaway_pair_d08_L64_N128_lev0` | 128 | 200 | 1 | **1.10** | 5.5 |
| `runaway_pair_d10_L64_N128_lev0` | 128 | 200 | 1 | **1.09** | 5.5 |
| `runaway_pair_d12_L64_N128_lev0` | 128 | 200 | 1 | **1.09** | 5.4 |
| `runaway_pair_d16_L64_N128_lev0` | 128 | 200 | 1 | **1.10** | 5.5 |
| `canonical_w075_L64_N128_lev0` | 128 | 120 | 1 | **0.66** | 5.5 |
| `canonical_w080_L64_N128_lev0` | 128 | 120 | 1 | **0.65** | 5.4 |
| `canonical_w085_L64_N128_lev0` | 128 | 120 | 1 | **0.66** | 5.5 |
| `canonical_w090_L64_N128_lev0` | 128 | 120 | 1 | **0.66** | 5.5 |
| `control_pair_mm_d10_L64_N128_lev0` | 128 | 200 | 4 | *1.41* | *7.0* ⚠ |
| `control_mirror_mp_d10_L64_N128_lev0` | 128 | 200 | 4 | *1.41* | *7.0* ⚠ |

⚠ **These two rows are contaminated and must not be used for costing.** They
ran through the window in which 133 orphaned solver ranks were still on the
machine, which cut the evolution rate from ~183 to 24–56 units of t per hour
until the ranks were swept. Their physics is unaffected — the mirror reproduces
the archive to two parts in a hundred thousand, and every timestep is present at
uniform spacing — but the clock they were measured on was wrong, and the true
cost of these cells is the 5.5 h/1000 t measured on the clean archive. The same
caveat applies to the PP and mass-scale cells of the same wave.

At `N = 128`, `L = 64`, no refinement, the cost is **5.5 GPU-hours per 1000
units of t** and it does not vary by more than 2% across ten cells — separation,
sector signs and star frequency all cost the same. So the run length is the only
thing that sets the bill at this grid, and the estimates for the new N=128 cells
(~1.1 h at t=200) are measurements, not guesses.

Rows for the cells this campaign adds go in as they finish; the resolution
rungs are the interesting ones, since the expected scaling is `N⁴` (three
dimensions plus the shorter timestep) — that is what predicts ~5.5 h at N=192
and ~17 h at N=256, and the ledger is where that prediction gets tested.

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

---

## 5. The recentring box ("running road") — realisation plan

This is outlook item 2 above, worked out far enough to build. The chase cell is
long only because the pair travels; if the pair is periodically carried back to
the middle of a small box, the box never has to grow. Cost stops scaling with
the target speed and becomes linear in run time.

### What moves is the data, not the box

The domain stays exactly what every mixed cell already uses: 64 on a side, 128
cells, centre at 32. Nothing about the grid, the extraction radii or the sponge
changes. When the pair's centroid has drifted 8 units toward the front wall,
every field is copied back by **exactly 16 cells**, which puts the pair back at
the centre, and 8 is added to a running odometer:

```
true displacement = position on the grid + odometer
```

16 cells is a whole number of cells, so the copy is a pure relabelling of which
cell holds which value — **no interpolation, exact to the last bit**. This is
the same reasoning that fixed the initial-data alignment: shifts that land on
the mesh are free, shifts that land between cells are not, so only the former
are ever taken.

Reaching 0.3c means travelling 210 units, so about **26 shifts over the whole
run**. Each one moves roughly 600 MB inside the card and takes milliseconds.
Early on the pair barely moves and shifts are rare; near 0.3c one falls every
27 units of time.

### Why not move the coordinates instead

The obvious alternative is to leave the data alone and let the coordinates chase
the pair — a comoving gauge. It is rejected on presentation grounds, not
technical ones. The entire claim of the paper is *the pair accelerates*.
Measuring that acceleration in a gauge purpose-built to follow the pair is a
much harder result to defend, and invites exactly the objection the whole
control matrix exists to kill. Translating the data leaves the evolution
equations, the gauge and the boundary untouched; only the labels change.

### The two seams, and why neither one reaches the stars

| face | what happens | distance from the pair | sits in |
|---|---|---|---|
| front (leading) | 16 cell-layers have no source and must be invented | 32 | sponge band |
| back (trailing) | 16 cell-layers of already-departed wake are discarded | 32 | sponge band |

The front sliver is filled from the outermost surviving layer using the same
1/r falloff the outer boundary condition already assumes, so the fabricated
values are continuous with what the boundary would have produced anyway. The
back sliver is thrown away; it holds radiation that has already left and cannot
travel back inward to reach the pair before the next shift.

**Nothing inside radius 24 is ever fabricated or discarded.** Both seams live
entirely in the outer shell.

### The sponge is what makes this legal

The sponge already exists (`Source/Grids/SpongeZone.hpp`), is already switched
on in every mixed cell, and already covers the band from radius 24 to 32 — it
is the thing that bought evolution past t = 60 in the first place. It is
anchored at fixed radius about the box centre, and after every shift the pair is
back at that centre, so **the sponge's geometric relationship to the pair is
restored exactly, every time**. Without it the fresh seam would ring and the
ringing would eventually turn around; with it the seam is absorbed before it
can.

The important point for a referee: relative to the fixed 64-box already used for
every mixed cell, the recentring box adds **no new approximation**. That box
already truncates the domain at radius 32 and already sponges everything outside
24. The recentring box keeps that same truncation centred on the pair instead of
letting the pair drift into one wall — which is strictly the safer of the two,
because the pair never approaches the boundary at all.

### What this cell cannot answer

The wake is discarded, so this run yields **no radiated momentum and no wave
zone**. It answers how fast the pair goes and nothing else. The wave-zone cell
in phase 4 remains the only source for the radiation field, and stays in the
matrix.

### Deciding when to shift

After each step, one reduction over the grid gives the mass-weighted centroid of
the pair. The same reduction pattern is already in the binary and already proven
on the cards — it is what steers the pump's spotlight
(`Examples/RadialRecipe/RLLumpTracker.hpp`) — but the recentring module gets its
own copy rather than reaching into the reinforcement-learning path, so the two
stay uncoupled.

Threshold: shift when the centroid is more than 8 units from the box centre.
Large enough that shifts are rare, small enough that the pair stays deep inside
the sponge's inner radius at all times.

### Checkpoints

Every cell in this campaign currently runs with checkpointing switched off,
which is right for a 1.1-hour cell and wrong for this one — it is between 8 and
50 hours. Switch it on. A checkpoint is about 1.2 GB at 128 cells; keeping a
rolling pair is nothing against the 3.5 TB free.

One piece of genuinely new state has to go into the checkpoint: **the odometer**.
Without it a restart resumes with the trajectory silently reset to zero, and the
resulting plot looks plausible while being wrong — the worst possible failure
mode. It is written and read back explicitly, and the restart is verified by
checking that the first line after a restart continues the last line before it.

### Code shape

Obeys the module rule: new behaviour lives in its own file, behind its own
switch, writing its own output.

| piece | where | note |
|---|---|---|
| the shift and the odometer | new `Source/Grids/GridTreadmill.hpp` | self-contained |
| the switch | one parameter, **off by default** | every archived cell is unaffected |
| the log | new `small_data/treadmill.dat` — step, time, cells shifted, cumulative offset | no existing column contract widened |
| the call | the post-timestep hook already present in the evolution level, coarsest level only | |

The shift itself is about fifteen lines: copy the state into a scratch buffer,
relabel the scratch buffer's boxes 16 cells back in x, copy it back by box
intersection, then fill the front sliver. The framework provides both halves of
that directly.

**Restriction:** single-level only. Every cell in this campaign is single-level
by decision, so this is not a limitation in practice, but the module refuses to
start with a hard error if refinement is on rather than silently doing something
wrong on a refined grid.

### Validation — three cells, and no result is believed before them

1. **Forced-shift null (~1.1 h).** The headline cell, shifting on a fixed
   schedule with zero net displacement: 16 cells forward, then 16 cells back.
   The trajectory must come out unchanged to round-off. This isolates seam
   damage from physics — if anything is wrong with the copy or the fill, it
   shows here with nothing else to hide behind.
2. **Recentred against the archive (~1.1 h).** The headline cell with
   recentring armed and the threshold dropped to 2 units, so several shifts
   land inside t ≤ 200. The trajectory must reproduce the archived cell.
   *Pass gate:* displacement and acceleration agree with the archive to the
   same tolerance the mirror control met — 1.4e-04 relative.
3. **Box-size check, once the pair is past 0.3c (~15 h).** Recentring at 96 on
   a side and 192 cells (same cell size) against 64/128 over the same window.
   This is the one that actually matters at speed: the outer boundary assumes a
   static centre, and a pair moving at a third of light speed is not static
   relative to the grid. Cells 1 and 2 cannot detect this, because the pair is
   barely moving by t = 200.

### What it buys

| target | run length t | long box, as planned | recentring box, 64 / 128 |
|---|---|---|---|
| 0.10 c | 430 | 0.2 GPU-days | **2.4 h** |
| 0.30 c | 1,360 | 1.7 GPU-days | **7.5 h** |
| 0.50 c | 2,500 | 7 GPU-days | **14 h** |
| 0.70 c | 4,200 | 28 GPU-days | **23 h** † |
| 0.90 c | 8,900 | 180 GPU-days | **49 h** † |

† past roughly 0.6c the Lorentz-contracted star needs half the cell size, which
multiplies these by about 16 — still days rather than months.

Memory at the 0.3c target falls from 33 GB to 6 GB, so the run fits on one card
with room to spare instead of filling it.

The headline: **the chase cell drops from 1.7 GPU-days to about 7.5 GPU-hours**,
and — unlike the long box — its cost no longer grows with the target speed.

### The cell

`chase_pair_d08_v03c_recentred_L64_N128_lev0` replaces
`chase_pair_d08_v03c_Lx352_L64_N128_lev0` in the scoreboard once validation
cells 1 and 2 pass. The long box stays on the books as the fallback if they do
not, and as the cross-check if the follow-up paper wants one.
