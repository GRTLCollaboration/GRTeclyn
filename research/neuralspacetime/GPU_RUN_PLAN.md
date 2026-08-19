# GPU run plan — closing the gap between the post-bugfix campaigns and the paper

Status 2026-08-18. Branch `feature/interstellar`. Companion to the article
at [`article/research.tex`](article/research.tex) (the live manuscript;
`researchnew.tex` is a legacy draft), [`Debug.md`](Debug.md),
[`DebugPreGPU.md`](DebugPreGPU.md), and the campaign packs under
[`results/`](../../results/).

This is the end-to-end queue for the remaining GPU work, in dependency order,
with the decision each block exists to make, and (§12) the exact launch
commands so each slot can fire the moment the GPUs free up. Rule inherited
from both debug logs: **a number goes into the paper only after it has a
production-tier run, a stated acceptance criterion fixed before the run, and
a pack under `results/` that regenerates it.**

---

## 0. Scoreboard — what the paper needs, and where each piece stands

*Last updated 2026-08-19 04:15. Update this block whenever a run lands; the
numbered gaps map to §1 "What the paper still lacks".*

> **Capacity changed 2026-08-19: we have 3 GPUs, not 4.** GPU 3 is reserved
> for other users — never schedule it (§2). Every "4 GPUs" and every
> wall-clock estimate further down this document predates that and runs ~4/3
> long.

| # | Goal | Status |
|---|---|---|
| 1 | In-scope (gated) CMA-ES refinement stage | ✅ **done** — `qball-trajectory-fgeo-max-refinement`, champion eval 193 |
| 2 | Production replay of a post-fix champion | ⚠️ **partly** — both depth exhibits ran; the *gated headline* died mid-run |
| 3 | Convergence + domain matrix for the new result | ⬜ not started (Phase 3, gated on #2) |
| 4 | Depth numbers are lower bounds → measure the true peak | 🔄 **in progress** — first t = 64 runs finishing now |
| 5 | Mechanism for the new champions (ablation + free-fall) | ⬜ not started |
| 6 | Post-fix canonical-only control under the gated objective | ⬜ not started |
| 7 | Dead numbers in the draft | ✅ resolved by exclusion — candidate 146 left the paper (§12.9) |
| 8 | Efficiency accounting / exotic-matter frontier | ⬜ not started (Phase 4) |

**Running now**

| run | source | where | state |
|---|---|---|---|
| DEPTH-X | depth eval 195 | GPU 1 | t ≈ 60 of 64, healthy — finishing ~04:05 |
| DEPTH-A | depth eval 185 | GPU 2 | t ≈ 60 of 64, healthy — finishing ~04:05 |

**Done this cycle**

- **Resolution ceiling measured, not extrapolated.** N = 240 → 49.8 GB and
  29.4 code units/h; N = 256 → 62 GB at start and ~20.5 units/h; N = 288 →
  OOM. Per §5 the intermediate rung therefore **upgrades from 224 to 240**
  (applies to FMAX-RI and DEPTH-RI).
- **MPI re-verified and the §2 constraint is now wrong** — see below.

**Failed this cycle**

- **FMAX-RM (gated headline, eval 193) — Arena OOM at t = 36.6 of 64.** The
  Arena grew 62 → 77 GB as the matter dispersed and more cells were tagged;
  it aborted asking for 17 MiB more. `checkpoint_interval` was `-1`, so there
  was no restart point and the run was lost outright. Artifacts deleted
  2026-08-19. **This is the blocker: the Phase-2 exit gate freezes the FMAX
  champion, so Phase 3 cannot open until FMAX-RM is re-run** — with
  checkpoints on, and at a grid/refinement setting that leaves headroom (see
  the wrapper README's "Arena OOM part-way through a long AMR run").
- N = 288 memory probe — OOM, expected and recorded (§12.4).

**Corrections pending against this document**

- §2 says *"MPI is broken; never raise RANKS."* **That is no longer true.**
  On 2026-08-19 `mpirun` was verified working on the current node at 1 and 4
  ranks with correct rank ids; the old failure belonged to a node the pod has
  since left. What remains unproven is (a) GRTresna multi-rank solves and (b)
  RadialRecipe MPI+CUDA, which last crashed on the *first* AMR advance in
  July — on a different node, with code that has changed since, and untested
  since. Single-GPU AMR is unaffected and is what runs today. §340's "no MPI
  repair unless the matrix fails its acceptance band" trigger has now fired.
- Missing manifests: the **depth mini-ladder has none**, and the fgeo_max
  Phase-3 manifest lacks **FMAX-PF** and the free-fall companions that §6
  requires.

---

## 1. Where we stand, and the scope ruling

All three post-bugfix search campaigns are done and packed. Everything they
measured is **search tier** (L = 64, N = 128, max_level 1, t ≤ 32) — nothing
is production-tier yet.

| campaign (pack) | objective | headline |
|---|---|---|
| [`qball-trajectory-evolving-geodesic-shortcut-search`](../../results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-evolving-geodesic-shortcut-search/CAMPAIGN_RESULTS.md) | `f_geo_max` — **persistence-gated** (the paper's agenda) | eval 322: 24 % depth with 25 % retention; ≥ 2 phantom lumps required |
| [`qball-trajectory-geodesic-depth-search`](../../results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-geodesic-depth-search/CAMPAIGN_RESULTS.md) | `f_geo_depth` — raw depth, no gating | 45.9 %, window-clipped; QD could not beat its seed |
| [`qball-trajectory-cmaes-refinement`](../../results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-cmaes-refinement/CAMPAIGN_RESULTS.md) | `f_geo_depth` | 48.38 % (eval 195) / 47.97 % clean-branch (eval 185), window-clipped |
| [`bondi-dipole-runaway`](../../results/bondi-dipole-runaway/README.md) | — | dressed ultraweak stars, stable solitons; separate paper; phase 2 on the GPUs now |

**Scope ruling.** The paper's contribution is the *gated* search — transport
credited only while the source holds together. That is the `f_geo_max`
lineage. The pure-depth lineage answers a different question ("how deep can
this ansatz cut, ignoring survival") and stays in the paper only as the
strength end of the strength-vs-persistence Pareto frontier — the same role
the rejected evals 99/136 play in the current draft — **not** as the
headline. (Decision 2026-08-18: both lineages *do* ship — the gated champion
as the headline with the full matrix, the depth champion as the second
result with a mini resolution ladder, §6 — because the two results explain
each other: depth proves the search saturates the known kinematic ceiling,
the gated run proves persistence selects a genuinely different object.)
Consequence: the CMA-ES refinement that exists
(`qball_traj_fgeo_depth_cmaes_v1`) refined the *out-of-scope* objective. The
paper's own ladder — MAP-Elites (`fgeo_v1`) → CMA-ES **under the same gated
rule** → HQ → matrix — is missing its refinement stage. That run comes first.

### What the paper still lacks

1. **The in-scope CMA-ES stage.** No CMA-ES has ever refined the
   persistence-gated champion (eval 322). The depth campaign proved the
   landscape is a narrow ridge that CMA-ES climbs well; the same treatment
   under `f_geo_max` is the missing rung of the paper's ladder.
2. **No production replay of any post-fix champion.** Nothing has been
   measured above N = 128 / max_level 1.
3. **No convergence or domain matrix for the new result** (the matrix in
   `article/research.tex` validates candidate 146 only).
4. **Depth numbers are lower bounds.** Every emission sweep still rises at
   its last allowed launch. The true peak has never been measured.
5. **Mechanism unknown for the new champions.** The depth exhibit's 48 %
   sits at the shift-tube frozen ceiling (48.5 %), far above the curvature
   branch (21.8 %). The ablation + free-fall metric must be re-run before
   any rewrite — the sourceability-barrier narrative may move.
6. **No post-fix canonical-only control** under the current pipeline and the
   gated objective.
7. **Dead numbers still in the draft** (per `Debug.md` PART A and
   `DebugPreGPU.md` PART E): E_pump 1.6e-17, "quasi-static scaffold"
   provenance, "5/5 rays at every launch" at t_emit = 12, "corroborated
   trapped surface at t≈27", 72 %→11 % confinement, absolute constraint
   norms from the diluted level-0 columns. **Resolved by exclusion
   (2026-08-18): candidate 146 leaves the paper entirely** — the rebuilt
   draft quotes only post-fix runs, and its heavy run tree was deleted
   (§12.9). The audit pack under `results/` stays in git as the record of
   the superseded iteration.
8. **No efficiency accounting.** The paper never answers the question
   referees respect most — *how much shortcut does a unit of exotic matter
   buy?* — even though both quantities (f_geo and the negative-energy
   integral) are measured per eval, and the draft's own
   quantum-inequality/Casimir discussion begs for exactly this frontier.
   Closed by the frontier search in Phase 4 (adopted 2026-08-18 from
   `research/nextsteps.md` item F3, together with the free-fall-certificate
   upgrade F7, §7). The rest of `nextsteps.md` belongs to other papers
   (wormhole arc, gw_beam/reward-hacking, Bondi) and stays out of this
   queue.

---

## 2. Standing constraints (apply to every run below)

- **3 usable GPUs — 0, 1 and 2 only — ~81 GB each, one node.** Changed
  2026-08-19: **GPU 3 belongs to other users and must never be scheduled.**
  Never pass `--gpu 3` or include 3 in a `--gpus` list, and cap every queue
  and fan-out at three concurrent slots. Anywhere below that says "4 GPUs"
  predates this and is wrong: wall-clock estimates stretch by ~4/3, and the
  CMA-ES population floor (next bullet) becomes **≥ 12**, not 16.
- **Single rank, for now.** The old RF (N = 384, 3 ranks) and DL (N = 320,
  2 ranks) **cannot currently be
  reproduced** — the proven single-GPU ceiling at max_level 3 is N = 256
  (BCMA-RM), and N = 384 **OOMs outright on one GPU** (confirmed 2026-08-18).
  The measured budget (bondi node): **8.8 GB at N = 128 scaling as N³** →
  ~30 GB at 192, ~47 GB at 224, ~70 GB at 256, **~100 GB at 288 — over the
  card even before AMR overhead**. So the ladder plans on 192/224/256 and
  treats anything above 256 as a probe-only curiosity (§12.4), never a
  dependency.
- **CMA-ES population ≥ 4 × GPU slots, never = GPU count** (the generation
  barrier starves the pipeline otherwise).
- **Composite constraint norms are mandatory.** Every paper-tier run must be
  on the current binary (17-column `constraint_norms.dat`); every quoted
  constraint figure is `L2_Ham_amr` / `Linf_Ham_amr` / refined-region, never
  the diluted cols 2–3.
- **Runs to t ≥ 40 need a denser metric stack**
  (`GRTECLYN_METRIC_STACK_N_SPACE=257`) and must pass the `cache_fidelity`
  gate on `min_chi`. Never `--force` a failed fidelity gate.
- **Emission-sweep honesty:** report per-launch bundle completeness; a
  launch without a full bundle is "unmeasured", never "zero".
- Stop detached campaigns only via `stop_campaign.sh` (orchestrator first),
  then verify with `pgrep` — killing workers only advances the queue.
- Packs scrub machine identity at pack time; grep before committing.
- The scorer changed 2026-08-05 (depth cap removed from all modes): new
  scores are not comparable to capped-era scores. Compare depths and
  retentions across campaigns, never scores.

---

## 3. Phase 0 — let bondi phase 2 finish (running now)

The `convA_mm_n{128,192,256}` and `convA_pm_sep12_n{128,192,256}` ladders
occupy the GPUs. Nothing here preempts them. On completion re-run
`research/bondi_dipole/pack_campaign.sh` and fire Phase 1 (chain command in
§12.2). No further bondi work is owed by *this* paper.

---

## 4. Phase 1 — the in-scope CMA-ES refinement (4 GPUs, ~1 day)

CMA-ES under **`f_geo_max`** (persistence-gated, now uncapped-linear ×
retention), warm-started from the `fgeo_v1` scoreboard — top row is the
eval-322 basin. Same search grid, stop time and gates as the QD stage, so
the ladder reads exactly like the paper's: quality-diversity finds the
basin, covariance adaptation climbs it.

- 200 evals, population 16 (4 × GPU slots), σ₀ = 0.05, checkpointed
  (resumable) after every generation.
- Expectation from the depth twin: the ridge is narrow (sub-percent moves),
  so gains are plausible in *both* depth and retention; if nothing improves
  on eval 322, eval 322 itself is promoted and the stage still closes the
  ladder honestly.
- Watch the gate-rejection rate; if it climbs past ~15 %, shrink σ₀ to 0.03.
  Never loosen the 0.03 Hamiltonian gate.

**Output:** the paper's champion genome (call it FMAX-champ), frozen before
anything else runs.

---

## 5. Phase 2 — promotion runs (4 GPUs, ~1 day)

| run | source | config | GPU |
|---|---|---|---|
| FMAX-RM (headline HQ) | FMAX-champ | L = 128, N = 256, max_level 3, regrid 0.02, **t = 64**, emission sweep t_emit = 0…48 step 2, 5 rays, 3 axes, dense metric stack | 0 |
| DEPTH-X (Pareto exhibit) | depth eval 195 (pack) | same config, objective `f_geo_depth` | 1 |
| DEPTH-A (fallback exhibit) | depth eval 185 (pack) | same — the clean-branch genome, so the branch question answers itself in one phase | 2 |
| memory probes | FMAX-champ | N = 288 and N = 240 at max_level 3, t = 2, watch `nvidia-smi` — 288 is expected to fail (for the record); if 240 fits it upgrades the intermediate rung from 224 | 3 |

(Candidate-146 re-measure runs were dropped 2026-08-18: the paper is rebuilt
on post-fix provenance only — §12.9. The GPU column is the nominal
assignment; in practice every run goes through the §12.10 queue and takes
whichever GPU frees first, with its constraint solve prestaged on CPU.)

**Decisions this phase makes:**

1. **Where the gated champion's advantage peaks** (t = 64 window ends the
   "lower bound" era).
2. **Mechanism branch** — ablation (shift-deleted / lapse-flattened /
   spatial-only) + free-fall f_ff on both FMAX-RM and DEPTH-X. If DEPTH-X is
   shift-dominated, the Pareto/ceiling section gains its strongest exhibit;
   if FMAX-champ is, the sourceability narrative itself is rewritten.
3. **Production-solve conditioning** — the champion's constraint residual at
   the N = 256 solve, quoted from composite norms. The DEPTH-X / DEPTH-A
   pair settles the branch question directly: if branch B degrades past the
   0.03 gate at the production solve, DEPTH-A's clean-branch 47.97 % is the
   exhibit.
4. **How the runs end** — collapse diagnostics over the full window; collapse
   called on lapse/chi/max|K| only, no AH-finder language.
5. **Ladder ceiling** — the probes fix RF′ and DL grids in the Phase-3
   manifest (320 if it fits, else 288, else the unigrid fallback).

**Exit gate:** freeze FMAX-champ + emission protocol; fill in the Phase-3
manifest grids. Nothing in Phase 3 launches before this.

---

## 6. Phase 3 — convergence and domain matrix on the frozen champion (~2 days)

Via a cloned promotion campaign (`promote/fgeo_max_cmaes_v1`, §12.5), same
framework that ran the candidate-146 matrix:

| cell | L | N | h | note |
|---|---|---|---|---|
| FMAX-RC | 128 | 192 | 0.667 | coarse |
| FMAX-RM | 128 | 256 | 0.500 | reference — reuse Phase 2 |
| FMAX-RI | 128 | 224 | 0.571 | intermediate rung — the ladder is **192/224/256**; N ≥ 288 does not fit one card (§2 memory table) and N = 384 OOMs, so 256 is the finest grid and the paper says so plainly |
| FMAX-DS | 96 | 192 | 0.500 | domain ladder |
| FMAX-DS2 | 112 | 224 | 0.500 | second domain rung — replaces the old DL: L = 160/N = 320 needs ~137 GB, impossible single-GPU; the domain series is L = 96/112/128 at fixed h |
| FMAX-PF | 128 | 256 | 0.500 | pump-free twin (t_pump = 0), t = 64 |
| freefall companions | — | — | — | `manifest_freefall.json` cells on the stored stacks |
| DEPTH-RC | 128 | 192 | 0.667 | **depth mini-ladder** — result 2 gets its own error bar |
| DEPTH-RI | 128 | 224 | 0.571 | same intermediate rung as FMAX |
| DEPTH-RM | 128 | 256 | 0.500 | = the Phase-2 DEPTH-X run, reused as the finest rung |

The depth mini-ladder (commands in §12.5) is deliberately small: its job is
to turn the ~48 % from a search-tier window-clipped number into a defensible
supporting result with a spread, not to be a second precision measurement.
Quote its grid-to-grid spread exactly like the bondi tables.

**Pre-registered acceptance (fix before launching):** fine–medium relative
difference ≤ 10 % on peak f_geo; f_ff ladder spread ≤ 10 %; full ray bundles
at the quoted launch; quote grid-to-grid **spread** as the error bar (no
order fit on f_geo — measure observed order on the composite constraint
norms instead). If AMR regridding noise spoils the ladder, fall back to the
bondi-proven unigrid methodology (max_level 0, three cell sizes,
`convergence_check.py`-style spread tables); unigrid ceiling on this node is
N = 256 (~70 GB).

---

## 7. Phase 4 — controls, frontier search, and close-outs

1. **Canonical-only control (4 GPUs, ~1 day).** MAP-Elites under the *same
   gated objective* (`f_geo_max`), current pipeline, every per-lump exotic
   dial pinned to 0, 200 evals, same gates and grid. Replaces the pre-fix
   canonical bound; keeps the "phantom sector required" claim commensurable
   with the headline. Verify on the first evals that
   `recipe_scalar_field_signs` is all +1.
2. **Efficiency-frontier search (FRONTIER-1, 4 GPUs, ~1 day; adopted from
   `nextsteps.md` F3).** MAP-Elites with the *same* gated objective but
   descriptor axes **(f_geo, exotic-energy budget)** — the Pareto map of
   "how much shortcut per unit negative energy, under full evolution". Both
   axes are already measured per eval; the only prep is registering a new
   descriptor mode (§12.2, a small default-off extension — none of the
   existing modes pairs the two). 200 evals, search tier, seeded from both
   lineages' champions so the frontier's two known corners are occupied
   from generation zero. Why it earns a GPU-day: it is the figure that
   unifies the paper — the headline and the depth champion become two
   labelled points on one measured frontier — and it converts the draft's
   quantum-inequality/Casimir discussion from prose into a plot. Command in
   §12.6.
3. **CPU-only close-outs (no GPU):** ablation + endpoint-gauge +
   proper-time tables for the new champion from stored stacks;
   `cache_fidelity` report attached to every quoted f_geo; calibrate the
   post-load gate's composite thresholds from the accumulated composite
   columns (closes `DebugPreGPU.md` PG-3); **free-fall certificate curve
   (`nextsteps.md` F7):** the free-fall observer probe is already built and
   env-gated — §12.1 now switches it on for every paper-tier replay, so
   each matrix cell records one certificate for free — and here its
   emission time is swept over the same grid as the geodesic sweep on the
   stored headline stack, giving f_ff(t_emit) alongside f_geo(t_emit). The
   free-faller number is the paper's most gauge-honest observable; this
   turns it from a single point into its strongest curve at zero GPU cost.

---

## 8. Phase 5 — the new-matter question (dressed stars)

**Recommendation: do not rebase this paper on the new matter model.** Ship
on the Q-ball ansatz it actually searched; run the cheap pilot so the
follow-up decision is made on data.

Why not rebase: the dressed-star model is validated for two *static* stars
(ω = 0.55 ultraweak rung); five boosted stars on tilted orbits is untested
physics, and folding it in means a new search plus a new validation stack —
months, for a paper that is otherwise one refinement + one matrix from done.
The paper already reports dispersal honestly as *the* limitation and the
bondi pack documents the cure; "stable-matter trajectory search" is the
natural flagship follow-up.

Why pilot anyway: the champions retain 14–25 % of their matter and referees
will ask whether the result is matter-model-limited. The plumbing exists
today — the bondi phase-2 reruns drive dressed stars through
`replay_eval.py` (`grtresna_bs_selfgrav=1`, ultraweak λ/μ, per-lump
`bs_omega`, per-sector profile tables).

| pilot | what | cost | keep/kill rule |
|---|---|---|---|
| NM-1 | replay FMAX-champ's 5-lump trajectory with dressed ultraweak stars in place of flat-space Q-ball profiles; search grid, t = 32; pump **off** (its target amplitude 0.075 is 4× the ultraweak star's φ_c = 0.0197 — pumping toward it would inject, not stabilise) | 1 GPU, hours | if the solve won't converge or f_geo ≈ 0, stop here |
| NM-2 | only if NM-1 shows f_geo > 0: 100-eval warm-started search | 4 GPUs, ~0.5 day | **≥ 3× retention at t = 32 at comparable depth → flagship follow-up campaign/paper**; outlook mention in this paper either way |

Solve tolerances for any NM run: the tightened bondi values
(`NL_exit_tolerance = 0.1 %`, `NL_stall_tolerance = 0.002`) — the loose gate
demonstrably sheds a χ ring that would contaminate the corridor. Main open
physics risk: boosted painting of orbiting dressed stars (bondi validated
at-rest only) — that is exactly what the pilot is for.

---

## 9. Phase 6 — paper integration gates

Before the rewrite is called done, walk both debug logs' claim tables:

- `Debug.md` PART A §1 ("MUST NOT say") — every row either re-measured in
  Phases 2–4 or absent from the text.
- `DebugPreGPU.md` PART E — resolved by exclusion: candidate 146 does not
  appear in the rebuilt paper; the post-fix campaigns' clean provenance
  carries the discovery narrative end to end.
- The depth lineage appears as **result 2** (search-tier tables + the
  DEPTH-X production point + the §6 mini-ladder spread), explicitly
  labelled as outside the gated objective.
- The frontier figure (FRONTIER-1 archive with FMAX-champ and DEPTH-X as
  labelled points) appears in the discussion, wired to the
  quantum-inequality/Casimir section; the f_ff(t_emit) curve appears next
  to the f_geo sweep wherever the headline number is quoted.
- Every quoted number names its source run and pack; every pack regenerates
  its tables by script. New packs: `qball-trajectory-fgeo-max-refinement`
  (Phase 1), `qball-trajectory-hq-promotion` (Phases 2–3), the canonical
  control and the frontier campaign (Phase 4).

---

## 10. Schedule sketch (4 GPUs, dependency order)

| slot | work | GPUs | est. wall |
|---|---|---|---|
| now | Phase 0: bondi mm + sep12 ladders (running) | 4 | ~1 day |
| 1 | Phase 1: CMA-ES refinement under `f_geo_max` | 4 | ~1 day |
| 2 | Phase 2: FMAX-RM + DEPTH-X + DEPTH-A + probes | 4 | ~1 day |
| 3 | analysis: freeze champion, fix manifest grids | — | hours |
| 4 | Phase 3: RC, RF, DS, DL, PF (+freefall companions) + depth mini-ladder | 4 | ~2–2.5 days |
| 5 | Phase 4: canonical-only control | 4 | ~1 day |
| 6 | Phase 4: FRONTIER-1 efficiency-frontier search | 4 | ~1 day |
| 7 | Phase 5: NM-1 pilot (NM-2 only on green light); f_ff curve + close-outs on CPU in parallel | 1–4 | 0.5–1 day |

Total ≈ 8–9 GPU-days after the bondi ladders clear — and less wall-clock
than that suggests, because the §12.10 queue prestages every constraint
solve on CPU and re-feeds freed GPUs continuously (the old chain pattern
idled each card up to ~70 min per cell during its solve). If the program
must be cut, FRONTIER-1 is the first thing to drop and the depth mini-ladder
the second — both improve the paper, neither blocks it.

---

## 11. Explicit non-goals

- No bigger boxes (reach is not the limitation; 512³ does not fit this node
  single-rank).
- No raising of the 5-lump cap (`DebugPreGPU.md` PG-9 part 3 — physics/cost
  decision, deferred).
- No MPI repair unless the matrix *fails* its acceptance band and a true
  N = 384 run becomes load-bearing.
- No GW-beaming campaign work (separate paper, separate queue).
- No full dressed-star search in this paper's queue (pilot only, §8).
- No shift-dominance/sourceability search (`nextsteps.md` F1/AT2) in this
  queue — it is its own campaign and likely its own paper; Phase 2's
  ablation already places both champions on the mechanism map, which is all
  this paper needs.
- No reward-hacking methods write-up here (`nextsteps.md` G5+N9) — it is
  publishable as-is from existing material, zero GPU, separate venue;
  harvest after this queue ships.

---

## 12. Exact launch commands

Everything runs from `grteclyn-wrapper/` unless stated. Conventions follow
[`results/bondi-dipole-runaway/LAUNCH.md`](../../results/bondi-dipole-runaway/LAUNCH.md)
and the wrapper [`README.md`](../../grteclyn-wrapper/README.md).

§12.3–12.7 give the per-run commands; **§12.10 wraps them in the work queue**
so the whole program runs hands-off — GPUs re-fed continuously, constraint
solves prestaged on CPU, nothing babysat. Prefer the queue; the raw commands
remain the reference for what each job does.

**Run tree (restructured 2026-08-18, see `runs/neuralspacetime/README.md`):**
everything this paper produces lives under one root —

```
runs/neuralspacetime/
  search/map_elites/   MAP-Elites campaigns (fgeo_v1, fgeo_depth_v1, pump_v2, canonical control)
  search/cma_es/       CMA-ES refinements (fgeo_depth_cmaes_v1, fgeo_max_cmaes_v1)
  hq/                  production replays + matrix; sources/ = frozen champions;
                       _cache/ = candidate-146 metric stacks (RM/RC/RF/DS/PF)
  experiments/         matter_model_ab/ + pump/ (pre-campaign A/B debugging)
  atlas/               geometry-first atlas campaigns
  pilots/              dressed-star pilots (NM-1/NM-2)
  _logs/               monitor logs, gpu samples, .done markers
```

The legacy roots (`runs/grtresna_qd`, `runs/grtresna_cmaes`,
`runs/grtresna_promote`, `runs/_logs`) are retired; the shared engines
(`qd/run.sh`, `cmaes/run.sh`, `lib/promote_common.sh`,
`lib/pipeline_monitor.sh`, `hq/replay_eval.py`, `stop_campaign.sh`) default
to the new tree. Dormant campaign families (rl, gw_beam, splash,
general_ftl) still name the old roots in their own launchers — pass
`RUNS_DIR` or patch them when those queues wake up. `runs/bondi_rerun` and
`runs/rotating_wormhole` are separate projects, untouched.

### 12.1 Conventions that bite

- **Detached launches spell out `/usr/bin/env`.** On this node other users'
  `PATH` entries shadow `env` with a sourceable snippet that exits 0 without
  running anything — a bare `env` launch writes an empty log and looks
  merely slow. (Measured 2026-08-17; see LAUNCH.md.)
- **Verify every detached launch with `pgrep -f <run name>`** before
  trusting its log.
- **Stop** only with
  `bash scripts/campaigns/stop_campaign.sh <run dir>` then confirm
  `pgrep -af "grteclyn_wrapper|grtresna|main3d"` is empty of that run.
- **Pump convention: the pump runs for the ENTIRE simulation, never stopped
  mid-run** (`RL_PUMP_STOP_TIME=-1`). A `params.txt` *without* an
  `rl_pump_stop_time` line already means this — the wiring writes the key
  only for values ≥ 0 and the evolution default is −1 = never stop. Because
  a negative value also erases the scorer's fallback emit floor, every pin
  of `RL_PUMP_STOP_TIME=-1` must be paired with
  `GEODESIC_EMIT_MIN_TIME=4`. The only deliberate pump-off runs are the two
  controls: the `manifest_pumpfree.json` twin (t_pump = 0, §12.5) and the
  NM-1 dressed-star pilot (`rl_pump_stop_time=0`, §12.7). Nothing else may
  set a non-negative value.
- Common env for every paper-tier (t ≥ 40) replay:
  `GRTECLYN_GEO_EMIT_INTERVAL=2 GRTECLYN_GEO_MAX_EMISSIONS=25
  GRTECLYN_METRIC_STACK_N_SPACE=257
  GRTECLYN_FREEFALL_OBSERVER_TIMING=1`. The last flag turns on the built-in
  free-fall observer certificate (default emission τ = 4.0, override with
  `GRTECLYN_FREEFALL_EMISSION_TAU`) — one
  `small_data/freefall_observer_timing.json` per replay, no extra cost. The
  full f_ff(t_emit) curve is a Phase-4 CPU close-out: re-run the probe over
  the emission grid on the stored stack.

### 12.2 Prepare now (no GPU needed) + chain off bondi

```bash
cd grteclyn-wrapper

# 1. DONE — the promotion campaign clone EXISTS with every pin baked in:
#    scripts/campaigns/promote/fgeo_max_cmaes_v1/
#    campaign.env.sh pins (each with its rationale in a comment there):
#      OBJECTIVE_MODE=f_geo_max, SCORE_EXOTIC_PENALTY_WEIGHT=0,
#      RL_PUMP_STOP_TIME=-1 + GEODESIC_EMIT_MIN_TIME=4 (§12.1),
#      GRTRESNA_RANKS=1, GRTRESNA_MAX_HAM_PCT=5 / MOM_PCT=5,
#      GRTECLYN_METRIC_STACK_N_SPACE=257, GRTECLYN_FREEFALL_OBSERVER_TIMING=1
#    manifest.json / manifest_freefall.json carry NO pump key (absent key =
#    pump on for the entire run); manifest_pumpfree.json alone keeps
#    rl_pump_stop_time=0 and declares "pump_off_control": true — the
#    deliberate control (§12.1).
#
#    The pump rule is LAUNCHER-ENFORCED now, not a checklist item:
#    - promote/lib/run_matrix.sh runs lib/validate_pump_convention.py on
#      every launch and --list; a manifest or env that stops the pump
#      without "pump_off_control": true is refused (exit 3), and pump-on
#      manifests require GEODESIC_EMIT_MIN_TIME in the env.
#    - qball_trajectory/cmaes_run.sh requires an explicit RL_PUMP_STOP_TIME
#      (no silent default; exit 2) and refuses a negative value without an
#      explicit floor.
#    Spec: grteclyn-wrapper README "Pump convention (enforced at launch)".
#    Tests: tests/scripts/test_pump_convention.py (11 tests).

# 2. Validate without GPUs (the pump validator runs inside --list too;
#    smoke-verified: all three manifests list clean, the retired bicomplex
#    template manifest is refused with exit 3):
DRY_RUN=1 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh --list

# 2b. Register the frontier descriptor mode (code task, no GPU — needed
#     before the Phase-4 FRONTIER-1 launch, §12.6). Add a mode named
#     `exotic_frontier` to the QD descriptor module + the CLI choices list:
#     axis 1 = evolving f_geo, axis 2 = log10 of the exotic-energy integral.
#     SIGN CONVENTION (two opposite ones live in the repo — do not mix):
#     axis 2 MUST use `integral_negative_rho`, the POSITIVE magnitude from
#     constraint_norms col 6 — template: search/geometry_atlas/score.py
#     (log10(max(value, floor)), the paper's E_- convention). Do NOT wire it
#     to warpfactory's `total_negative_energy`, which is SIGNED (<= 0):
#     log10 of it is NaN for every candidate, and since the PG-8 fix maps
#     NaN to bin 0 instead of raising, the whole archive would silently
#     collapse into one column and look like a converged frontier.
#     Own extension, default-off, unit test alongside the existing
#     descriptor tests (include a negative-input case); smoke it with
#     QD_TARGET_EVALS=2 before the real launch.

# 3. Chain Phase 1 to fire when the bondi ladders finish (detached).
#    NOTE: superseded by the §12.10 queue for everything replay-shaped —
#    keep this until-loop only for the CMA-ES stage itself, which manages
#    its own 4-GPU pipeline and is not a queue job:
setsid nohup bash -c 'until ! pgrep -f "bondi_sg_pair" >/dev/null; do sleep 600; done; \
  bash scripts/campaigns/qball_trajectory/cmaes_run.sh' \
  > ../runs/neuralspacetime/search/cma_es/fgeo_max_chain.log 2>&1 < /dev/null & disown
# (env block for cmaes_run.sh as in §12.3 — or reuse
#  scripts/campaigns/bondi_dipole/chain_on_gpu_free.sh as the template.)
```

### 12.3 Phase 1 — CMA-ES under `f_geo_max`

```bash
cd grteclyn-wrapper
setsid nohup /usr/bin/env \
  RUN_NAME=qball_traj_fgeo_max_cmaes_v1 \
  OBJECTIVE_MODE=f_geo_max \
  WARM_START_TRAJECTORY="$PWD/../runs/neuralspacetime/search/map_elites/qball_traj_fgeo_v1/trajectory.jsonl" \
  WARM_START_TOP_K=1 WARM_START_JITTER=0.05 \
  POPULATION=16 MAX_GENERATIONS=13 TARGET_EVALS=200 \
  SIGMA0=0.05 SEED=11 \
  GPU_IDS="0 1 2 3" MAX_CONCURRENT_GRTRESNA=6 \
  GRTRESNA_TIMEOUT=2400 RANKS=1 \
  STOP_TIME=26 \
  RL_PUMP_STOP_TIME=-1 GEODESIC_EMIT_MIN_TIME=4 \
  SCORE_EXOTIC_PENALTY_WEIGHT=0 \
  POSTLOAD_MAX_HAM_L2=3e-2 POSTLOAD_MAX_MOM_L2=1e-2 \
  bash scripts/campaigns/qball_trajectory/cmaes_run.sh \
  > ../runs/neuralspacetime/search/cma_es/qball_traj_fgeo_max_cmaes_v1_launch.log 2>&1 < /dev/null & disown
```

Notes: `STOP_TIME=26` mirrors the `fgeo_v1` physics window so the refinement
is stage-comparable to its QD parent (the depth twin used 32 — do not mix).
`WARM_START_TOP_K=1` centres on the eval-322 basin. `SEED=11` decorrelates
from the seed-7 noise every earlier campaign replayed.

Four pins exist because `cmaes_run.sh` carries bicomplex-era defaults that
do **not** match the `fgeo_v1` parent (run_fgeo.sh → run.sh) this campaign
warm-starts from — without them the refinement optimizes different physics
than it inherits:

- `RANKS=1` — the search path reads `RANKS` (`lib/search_common.sh:27`,
  default **8** → mpirun segfault). `GRTRESNA_RANKS` is the *promote*-path
  spelling (`lib/promote_common.sh:18`); setting the wrong one here is
  silent.
- `RL_PUMP_STOP_TIME=-1` — cmaes_run.sh defaults to pump-off-at-4; the
  parent ran the pump for the full window (-1).
- `GEODESIC_EMIT_MIN_TIME=4` — with the pump at -1 the fallback ray-launch
  floor (the score code reads `GEODESIC_EMIT_MIN_TIME`, else a
  *non-negative* `RL_PUMP_STOP_TIME`) would vanish; the parent had the
  floor at 4 via run.sh's default. Restate it explicitly.
- `SCORE_EXOTIC_PENALTY_WEIGHT=0` — cmaes_run.sh defaults 0.2; the parent's
  f_geo_max scoring has no exotic penalty (run_fgeo.sh pins 0).

Monitor / resume / stop (note: resume must repeat the same env block):

```bash
grep -a "^\[optimize\] \(eval\|gen\)" ../runs/neuralspacetime/search/cma_es/qball_traj_fgeo_max_cmaes_v1_launch.log | tail
RESUME=1 RUN_NAME=qball_traj_fgeo_max_cmaes_v1 \
  OBJECTIVE_MODE=f_geo_max RANKS=1 STOP_TIME=26 \
  RL_PUMP_STOP_TIME=-1 GEODESIC_EMIT_MIN_TIME=4 SCORE_EXOTIC_PENALTY_WEIGHT=0 \
  bash scripts/campaigns/qball_trajectory/cmaes_run.sh
bash scripts/campaigns/stop_campaign.sh ../runs/neuralspacetime/search/cma_es/qball_traj_fgeo_max_cmaes_v1
```

### 12.4 Phase 2 — promotion, exhibit, re-measure, probes

Headline HQ through the promotion framework (freeze picks the campaign best;
pass `SOURCE_EVAL_ID=<n>` to override):

```bash
cd grteclyn-wrapper
bash scripts/campaigns/promote/fgeo_max_cmaes_v1/freeze.sh
GPU_ID=0 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-RM
```

Depth-lineage Pareto exhibit, replayed straight from the pack (same pattern
as the live bondi reruns):

```bash
setsid nohup /usr/bin/env \
  GRTECLYN_GEO_EMIT_INTERVAL=2 GRTECLYN_GEO_MAX_EMISSIONS=25 \
  GRTECLYN_METRIC_STACK_N_SPACE=257 GRTECLYN_FREEFALL_OBSERVER_TIMING=1 \
  GEODESIC_EMIT_MIN_TIME=4 \
  .venv/bin/python scripts/campaigns/hq/replay_eval.py \
  ../results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-cmaes-refinement/run/eval_000195 \
  --name depth195_hq_L128_N256_t64 --runs-dir ../runs/neuralspacetime/hq \
  --gpu 1 --n-full 256 --l-full 128 --max-level 3 --regrid-threshold 0.02 \
  --stop-time 64 --plot-interval 72 \
  --objective-mode f_geo_depth --evolving-geodesic \
  --grtresna-ranks 1 --grtresna-iterations 50 --grtresna-timeout 7200 \
  --consumer-radii 12 18 24 \
  > ../runs/neuralspacetime/hq/depth195_hq.launch.log 2>&1 < /dev/null & disown
```

DEPTH-A, the clean-branch fallback exhibit, identically configured (the
branch-B conditioning question then answers itself inside Phase 2). Eval
185's directory was pruned from both the live tree and the pack — only its
genome ships (`genomes/best_clean_branch_eval185_familyA.json`). Its
`overrides` block holds **only the 39 searched dials**; the 30 fixed keys
(matter model, μ, λ, mass, ω, lump count, well depths, grid) are campaign
constants that never entered the genome. Copying the genome alone as
`metadata.json` would replay the wrong object — materialize the stub by
overlaying the 185 dials on a sibling eval's full metadata (195's fixed
keys are the same campaign constants):

```bash
# (Already materialized during the pre-GPU smoke run — 69 overrides, loads
# through replay_eval's own reader, no pump key. Re-running is idempotent.)
mkdir -p ../runs/neuralspacetime/hq/sources/depth_eval185_stub
.venv/bin/python - <<'PY'
import json, pathlib
pack = pathlib.Path('../results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-cmaes-refinement')
base = json.loads((pack / 'run/eval_000195/metadata.json').read_text())
genome = json.loads((pack / 'genomes/best_clean_branch_eval185_familyA.json').read_text())
assert set(genome['overrides']) <= set(base['overrides'])
base['overrides'].update(genome['overrides'])
base['stub_note'] = ('eval-185 familyA searched dials overlaid on eval-195 '
                     'fixed keys (185 dir pruned from pack; genome = 39 '
                     'searched dials only)')
out = pathlib.Path('../runs/neuralspacetime/hq/sources/depth_eval185_stub/metadata.json')
out.write_text(json.dumps(base, indent=2) + '\n')
print('wrote', out, '-', len(base['overrides']), 'overrides')
PY

setsid nohup /usr/bin/env \
  GRTECLYN_GEO_EMIT_INTERVAL=2 GRTECLYN_GEO_MAX_EMISSIONS=25 \
  GRTECLYN_METRIC_STACK_N_SPACE=257 GRTECLYN_FREEFALL_OBSERVER_TIMING=1 \
  GEODESIC_EMIT_MIN_TIME=4 \
  .venv/bin/python scripts/campaigns/hq/replay_eval.py \
  ../runs/neuralspacetime/hq/sources/depth_eval185_stub \
  --name depth185_hq_L128_N256_t64 --runs-dir ../runs/neuralspacetime/hq \
  --gpu 2 --n-full 256 --l-full 128 --max-level 3 --regrid-threshold 0.02 \
  --stop-time 64 --plot-interval 72 \
  --objective-mode f_geo_depth --evolving-geodesic \
  --grtresna-ranks 1 --grtresna-iterations 50 --grtresna-timeout 7200 \
  --consumer-radii 12 18 24 \
  > ../runs/neuralspacetime/hq/depth185_hq.launch.log 2>&1 < /dev/null & disown
```

Memory probes for the Phase-3 intermediate rung (kill after t = 2; ~10 min
each — 288 is expected to OOM per the §2 memory table; run it once for the
record, and if 240 fits use 240 instead of 224):

```bash
for N in 288 240; do
  .venv/bin/python scripts/campaigns/hq/replay_eval.py \
    "$(ls -d ../runs/neuralspacetime/hq/sources/qball_traj_fgeo_max_cmaes_v1/eval_* | head -1)" \
    --name mem_probe_N${N} --runs-dir ../runs/neuralspacetime/hq/probes \
    --gpu 3 --n-full ${N} --l-full 128 --max-level 3 --regrid-threshold 0.02 \
    --stop-time 2 --objective-mode f_geo_max \
    --grtresna-ranks 1 --grtresna-timeout 7200
done
# watch in a second shell: nvidia-smi --query-gpu=index,memory.used --format=csv -l 10
```

Largest N that stays comfortably under ~75 GB becomes the FMAX-RI grid in
the manifest (224 is the safe default; nothing above 256 enters the queue).

### 12.5 Phase 3 — the matrix

```bash
cd grteclyn-wrapper
GPU_ID=0 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-RC
GPU_ID=1 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-RI
GPU_ID=2 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-DS
GPU_ID=3 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-DS2

# pump-free twin + freefall companions, same runner, different manifest.
# MANIFEST must be ABSOLUTE: every campaign shell cd's to GRTECLYN_ROOT
# (scripts/lib/env.sh:88), so wrapper-relative paths stop resolving:
MANIFEST="$PWD/scripts/campaigns/promote/fgeo_max_cmaes_v1/manifest_pumpfree.json" \
  GPU_ID=0 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh --list   # then launch its cell id
MANIFEST="$PWD/scripts/campaigns/promote/fgeo_max_cmaes_v1/manifest_freefall.json" \
  GPU_ID=1 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh --list   # then launch its cell ids
```

`DRY_RUN=1` any cell first; launch logs land in the campaign's
`VALIDATION_LAUNCH_LOG_DIR`.

Depth mini-ladder (result 2's error bar) — same pattern as the §12.4
DEPTH-X command, only the grid changes; slot each onto whichever GPU frees
first. DEPTH-RM is the Phase-2 run, reused as the middle rung:

```bash
for CELL in "depth195_hq_L128_N192_t64 192 2" "depth195_hq_L128_N224_t64 224 3"; do
  set -- $CELL   # $1 = name, $2 = N (224 → 240 if the §12.4 probe passed), $3 = GPU
  setsid nohup /usr/bin/env \
    GRTECLYN_GEO_EMIT_INTERVAL=2 GRTECLYN_GEO_MAX_EMISSIONS=25 \
    GRTECLYN_METRIC_STACK_N_SPACE=257 GRTECLYN_FREEFALL_OBSERVER_TIMING=1 \
    GEODESIC_EMIT_MIN_TIME=4 \
    .venv/bin/python scripts/campaigns/hq/replay_eval.py \
    ../results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-cmaes-refinement/run/eval_000195 \
    --name "$1" --runs-dir ../runs/neuralspacetime/hq \
    --gpu "$3" --n-full "$2" --l-full 128 --max-level 3 --regrid-threshold 0.02 \
    --stop-time 64 --plot-interval 72 \
    --objective-mode f_geo_depth --evolving-geodesic \
    --grtresna-ranks 1 --grtresna-iterations 50 --grtresna-timeout 7200 \
    --consumer-radii 12 18 24 \
    > "../runs/neuralspacetime/hq/$1.launch.log" 2>&1 < /dev/null & disown
done
```

### 12.6 Phase 4 — canonical-only control

`PIN_DIMS` **replaces** the launcher default, so the physics pins must be
restated alongside the five exotic pins:

```bash
cd grteclyn-wrapper
setsid nohup /usr/bin/env \
  QD_NAME=qball_traj_fgeo_canonical_v1 \
  QD_TARGET_EVALS=200 SEED=12 \
  PIN_DIMS="grtresna_scalar_mass=1.0 grtresna_scalar_lambda=640 grtresna_bs_omega=0.8 \
trajectory_lump0_well_depth=0.15 trajectory_lump1_well_depth=0.15 trajectory_lump2_well_depth=0.15 \
trajectory_lump3_well_depth=0.15 trajectory_lump4_well_depth=0.15 trajectory_well_width=1.667 \
trajectory_lump0_exotic=0 trajectory_lump1_exotic=0 trajectory_lump2_exotic=0 \
trajectory_lump3_exotic=0 trajectory_lump4_exotic=0" \
  GPU_IDS="0 1 2 3" RANKS=1 \
  bash scripts/campaigns/qball_trajectory/run_fgeo.sh \
  > ../runs/neuralspacetime/search/map_elites/qball_traj_fgeo_canonical_v1_launch.log 2>&1 < /dev/null & disown
```

The default warm start (prior elites) is kept deliberately: under the pins
they decode as the *same choreographies with the phantom sector removed* —
the strongest possible matched control. **Verify on the first completed
evals** that every `params.txt` carries `recipe_scalar_field_signs = 1 1 1 1 1`
and the exotic integral is 0 — and run it as a command, not an eyeball
(this is the ONLY guard: the automated sign rail is a no-op for the
bicomplex model — `sign_consistency.py` compares the requested signs
against themselves and can never fail there):

```bash
grep -H "recipe_scalar_field_signs" \
  ../runs/neuralspacetime/search/map_elites/qball_traj_fgeo_canonical_v1/eval_*/params.txt \
  | grep -v "= 1 1 1 1 1" && echo "SIGN PIN FAILED — stop the campaign" || echo signs-ok
```

Abort and fix the pins if it prints anything but `signs-ok`.

**FRONTIER-1 — efficiency-frontier search** (fires after the control;
requires the `exotic_frontier` descriptor mode from §12.2 item 2b). Same
gated objective, new descriptor axes; the seed list drops both lineages'
champions onto the frontier's known corners:

```bash
cd grteclyn-wrapper
# preflight: the mode must exist BEFORE detaching — an unregistered
# descriptor dies at argparse inside setsid/nohup, silently, into a log
# nobody is watching:
grep -q exotic_frontier src/grteclyn_wrapper/cli/parser.py \
  || { echo "exotic_frontier not registered — do §12.2 item 2b first"; false; }

setsid nohup /usr/bin/env \
  QD_NAME=qball_traj_fgeo_frontier_v1 \
  OBJECTIVE_MODE=f_geo_max \
  DESCRIPTOR_MODE=exotic_frontier \
  QD_TARGET_EVALS=200 SEED=13 \
  SEED_EVAL_DIRS="$(ls -d ../runs/neuralspacetime/hq/sources/qball_traj_fgeo_max_cmaes_v1/eval_* 2>/dev/null | head -1) \
../runs/neuralspacetime/search/map_elites/qball_traj_fgeo_v1/eval_000322 \
../results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-cmaes-refinement/run/eval_000195 \
../runs/neuralspacetime/hq/sources/depth_eval185_stub" \
  GPU_IDS="0 1 2 3" RANKS=1 \
  bash scripts/campaigns/qball_trajectory/run_fgeo.sh \
  > ../runs/neuralspacetime/search/map_elites/qball_traj_fgeo_frontier_v1_launch.log 2>&1 < /dev/null & disown
```

Setting `SEED_EVAL_DIRS` explicitly replaces `run_fgeo.sh`'s pump_v2
default — intended: four champion seeds anchor the corners and MAP-Elites
fills the frontier between them. The archive plot (f_geo vs log exotic
budget) drops straight into the discussion section next to the
quantum-inequality bound.


### 12.7 Phase 5 — dressed-star pilot NM-1

Working reference for the override block:
`scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh` (running these exact
overrides in production right now). Pilot = FMAX-champ's trajectory, matter
swapped to self-gravitating ultraweak stars, pump off, tightened solve:

```bash
cd grteclyn-wrapper
mkdir -p ../runs/neuralspacetime/pilots
setsid nohup /usr/bin/env \
  GEODESIC_EMIT_MIN_TIME=4 \
  .venv/bin/python scripts/campaigns/hq/replay_eval.py \
  "$(ls -d ../runs/neuralspacetime/hq/sources/qball_traj_fgeo_max_cmaes_v1/eval_* | head -1)" \
  --name fmax_dressed_pilot --runs-dir ../runs/neuralspacetime/pilots \
  --gpu 0 --n-full 128 --l-full 64 --max-level 1 --regrid-threshold 0.02 \
  --stop-time 32 --objective-mode f_geo_max --evolving-geodesic \
  --grtresna-ranks 1 --grtresna-iterations 50 \
  --grtresna-nl-exit-tolerance 0.1 --grtresna-nl-stall-tolerance 0.002 \
  --grtresna-timeout 7200 \
  --extra-override grtresna_scalar_lambda=10240 \
  --extra-override grtresna_scalar_mu=21845333 \
  --extra-override grtresna_bs_omega=0.55 \
  --extra-override grtresna_bs_selfgrav=1 \
  --extra-override rl_pump_stop_time=0 \
  > ../runs/neuralspacetime/pilots/fmax_dressed_pilot.launch.log 2>&1 < /dev/null & disown
```

Leave the champion's per-lump trajectory and exotic dials untouched (they
come from its params); do **not** set `grtresna_boost_lumps=0` — the orbits
need their boosts, and whether boosted dressed-star painting behaves is the
pilot's actual question. Sanity-check the t = 0 row of
`small_data/sector_barycenters.dat` against the star fingerprints in
LAUNCH.md §4 before believing anything downstream.

### 12.8 After each phase

```bash
# cache-fidelity check on every paper-tier (t >= 40) run of the phase.
# The cache_fidelity gate lives ONLY in this script — NOTHING on the launch
# path enforces it (replay_eval, run_batch and the promote framework never
# call it). It refuses to score an unfaithful metric stack; §2 forbids
# --force. A refusal here means the stack does not match the spacetime that
# evolved — the exact bug-#9 shape that once faked an FTL shortcut:
GEODESIC_EMIT_MIN_TIME=4 .venv/bin/python \
  scripts/campaigns/rl/score_evolving_geodesic.py <each t>=40 run dir>

# pack the finished campaign (create the pack script alongside the others,
# scrubbing machine paths like research/bondi_dipole/pack_*.sh do), then:
grep -rnE "/(home|users)/|$(whoami)" results/<new-pack>/ && echo "SCRUB FAILED" || echo clean
```

### 12.9 Candidate-146 exclusion and the `hq/` deletion (2026-08-18)

Decision: candidate 146 was discovered by the pre-fix pipeline (broken
genome round-trip, PG-1/PG-9), so the rebuilt paper excludes it entirely and
quotes only post-fix runs — no per-claim "this number survives the bug"
caveats. Consequences, all applied:

- `runs/neuralspacetime/hq/` (49 GB: the 146 metric-stack cache, freefall
  companions, frozen sources) was **deleted**; the dir is empty and ready
  for the FMAX promotion runs. The light record of the superseded iteration
  — validation mirrors, provenance, figures — remains in git under
  `results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/`
  (18 MB, 282 tracked files).
- The BCMA-RM / PF-RM re-measure runs were dropped from Phases 2 and 4; the
  freed GPU slot runs DEPTH-A instead.
- `scripts/campaigns/promote/bicomplex_cmaes_v1/` is retired as a live
  campaign and kept only as the clone template for `fgeo_max_cmaes_v1`
  (§12.2).

### 12.10 The queue — press and forget (2026-08-18)

Two pipeline problems, both fixed and tested:

1. **GPUs were chained pairwise, not queued.** Each successor watched one
   predecessor; nothing generic handed the next job to whichever GPU freed
   first. New: `scripts/campaigns/lib/gpu_queue.sh` — one detached runner
   per pool, one worker per slot, jobs claimed atomically from `pending/` in
   name order, failures isolated in `failed/`, follow-up jobs picked up
   whenever they appear. Tests: `tests/scripts/test_gpu_queue.py`.
2. **The constraint solve blocked the GPU.** Every production cell runs a
   CPU-only GRTresna solve (up to ~70 min at N = 256) before touching its
   card — the measured bondi bottleneck. New: `--solve-only` on
   `replay_eval.py` (and `SOLVE_ONLY=1` through the promote framework's
   `run_batch.sh`) runs just the solve + CPU gates and writes the
   `initial_data.gridinit`; the evolve job then reuses it via the
   already-existing `--gridinit` / `GRIDINIT=` path and starts on the GPU
   instantly. Tests: `tests/search/test_evaluation_phases.py`
   (solve-only cases).

**Start the runners once — today is fine.** The GPU workers gate on
`nvidia-smi` used memory, so they wait politely while bondi still owns a
card and take each device over the moment it clears (this supersedes the
§12.2 chain-off-bondi until-loop):

```bash
cd grteclyn-wrapper
QROOT="$PWD/../runs/neuralspacetime/_queue"
mkdir -p "$QROOT/gpu/pending" "$QROOT/solve/pending"

setsid nohup /usr/bin/env QUEUE_GPU_MEM_MAX_MB=2000 \
  bash scripts/campaigns/lib/gpu_queue.sh "$QROOT/gpu" 0 1 2 3 \
  > "$QROOT/gpu_runner.log" 2>&1 < /dev/null & disown
# CPU solve pool: 2 slots while bondi still needs the cores, raise to 4 after
setsid nohup /usr/bin/env \
  bash scripts/campaigns/lib/gpu_queue.sh "$QROOT/solve" s1 s2 \
  > "$QROOT/solve_runner.log" 2>&1 < /dev/null & disown
pgrep -af gpu_queue.sh   # verify both runners per §12.1
```

**The prestage pattern** (worked example: DEPTH-X). The evolve job is staged
next to the queue; the solve job releases it only when its gridinit is
ready, so the GPU pool never sees a job it cannot start instantly. The depth
exhibits' solves can be enqueued **today** — their genomes exist; FMAX cells
join after the Phase-1 freeze:

```bash
cat > "$QROOT/staged_200_depth195_evolve.job" <<EOF
cd "$PWD"
export GRTECLYN_GEO_EMIT_INTERVAL=2 GRTECLYN_GEO_MAX_EMISSIONS=25
export GRTECLYN_METRIC_STACK_N_SPACE=257 GRTECLYN_FREEFALL_OBSERVER_TIMING=1
export GEODESIC_EMIT_MIN_TIME=4
.venv/bin/python scripts/campaigns/hq/replay_eval.py \
  ../results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-cmaes-refinement/run/eval_000195 \
  --name depth195_hq_L128_N256_t64 --runs-dir ../runs/neuralspacetime/hq \
  --gpu "\$QUEUE_GPU" --n-full 256 --l-full 128 --max-level 3 \
  --regrid-threshold 0.02 --stop-time 64 --plot-interval 72 \
  --objective-mode f_geo_depth --evolving-geodesic \
  --gridinit ../runs/neuralspacetime/hq/presolve/depth195_presolve_N256/initial_data.gridinit \
  --consumer-radii 12 18 24
EOF

cat > "$QROOT/solve/pending/100_depth195_solve.job" <<EOF
cd "$PWD"
.venv/bin/python scripts/campaigns/hq/replay_eval.py \
  ../results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-cmaes-refinement/run/eval_000195 \
  --name depth195_presolve_N256 --runs-dir ../runs/neuralspacetime/hq/presolve \
  --n-full 256 --l-full 128 --max-level 3 \
  --objective-mode f_geo_depth \
  --grtresna-ranks 1 --grtresna-iterations 50 --grtresna-timeout 7200 \
  --solve-only
mv "$QROOT/staged_200_depth195_evolve.job" "$QROOT/gpu/pending/"
EOF
```

Repeat the pair for DEPTH-A (via its §12.4 stub), the depth mini-ladder
(§12.5 grids), NM-1, and — after the freeze — every FMAX matrix cell
through the promote framework:

```bash
# CPU prestage of a matrix cell (solve job body):
RUNS_DIR=../runs/neuralspacetime/hq/presolve SOLVE_ONLY=1 FOREGROUND=1 GPU_ID=0 \
  bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-RC
# matching evolve job body (glob resolves the presolved episode by prefix):
GRIDINIT="$(ls ../runs/neuralspacetime/hq/presolve/*RC*_hq_eval*/initial_data.gridinit)" \
  FOREGROUND=1 GPU_ID="$QUEUE_GPU" \
  bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-RC
```

**Rules that keep it hands-off:**

- Jobs run in the **foreground** of their worker (`FOREGROUND=1` for
  promote cells; never `setsid`/`nohup` inside a job) — the runner is the
  only thing detached.
- Enqueue evolve jobs only via their solve job's `mv` (or directly when no
  solve is needed, e.g. memory probes) — the GPU pool must never wait on a
  solve.
- On-the-fly plotfile extraction is already configured per tier
  (`CONSUMER_JOBS=2` for HQ AMR cells — the NFS-validated sweet spot; the
  bondi-style unigrid cells validated 8) and the drain is off for promote
  cells, so the post-GPU CPU tail per job is minutes, not the bondi hour.
- **Stop dispatch:** `touch "$QROOT/gpu/STOP"` (running cells finish; same
  for the solve pool). **Stop a running cell:** `stop_campaign.sh` on its
  run dir as always — the queue marks that job failed and moves on. Requeue
  with `mv failed/X pending/X`.
- One runner per pool (`runner.lock` enforces it); after a node reboot just
  restart the runners — stale claims are requeued automatically.
