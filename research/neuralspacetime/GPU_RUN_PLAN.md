# GPU run plan — closing the gap between the post-bugfix campaigns and the paper

Status 2026-08-18. Branch `feature/interstellar`. Companion to
[`researchnew.tex`](researchnew.tex), [`Debug.md`](Debug.md),
[`DebugPreGPU.md`](DebugPreGPU.md), and the campaign packs under
[`results/`](../../results/).

This is the end-to-end queue for the remaining GPU work, in dependency order,
with the decision each block exists to make, and (§12) the exact launch
commands so each slot can fire the moment the GPUs free up. Rule inherited
from both debug logs: **a number goes into the paper only after it has a
production-tier run, a stated acceptance criterion fixed before the run, and
a pack under `results/` that regenerates it.**

---

## 1. Where we stand, and the scope ruling

All three post-bugfix search campaigns are done and packed. Everything they
measured is **search tier** (L = 64, N = 128, max_level 1, t ≤ 32) — nothing
is production-tier yet.

| campaign (pack) | objective | headline |
|---|---|---|
| [`qball-trajectory-evolving-geodesic-shortcut-search`](../../results/qball-trajectory-evolving-geodesic-shortcut-search/CAMPAIGN_RESULTS.md) | `f_geo_max` — **persistence-gated** (the paper's agenda) | eval 322: 24 % depth with 25 % retention; ≥ 2 phantom lumps required |
| [`qball-trajectory-geodesic-depth-search`](../../results/qball-trajectory-geodesic-depth-search/CAMPAIGN_RESULTS.md) | `f_geo_depth` — raw depth, no gating | 45.9 %, window-clipped; QD could not beat its seed |
| [`qball-trajectory-cmaes-refinement`](../../results/qball-trajectory-cmaes-refinement/CAMPAIGN_RESULTS.md) | `f_geo_depth` | 48.38 % (eval 195) / 47.97 % clean-branch (eval 185), window-clipped |
| [`bondi-dipole-runaway`](../../results/bondi-dipole-runaway/README.md) | — | dressed ultraweak stars, stable solitons; separate paper; phase 2 on the GPUs now |

**Scope ruling.** The paper's contribution is the *gated* search — transport
credited only while the source holds together. That is the `f_geo_max`
lineage. The pure-depth lineage answers a different question ("how deep can
this ansatz cut, ignoring survival") and stays in the paper only as the
strength end of the strength-vs-persistence Pareto frontier — the same role
the rejected evals 99/136 play in the current draft — **not** as the
headline. Consequence: the CMA-ES refinement that exists
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
   `researchnew.tex` validates candidate 146 only).
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

---

## 2. Standing constraints (apply to every run below)

- **4 GPUs, ~81 GB each, one node, single rank.** MPI is broken; never raise
  RANKS. The old RF (N = 384, 3 ranks) and DL (N = 320, 2 ranks) **cannot be
  reproduced** — the proven single-GPU ceiling at max_level 3 is N = 256
  (BCMA-RM). Anything above 256 gets a short memory probe (§12.4) before it
  enters the queue.
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
| memory probes | FMAX-champ | N = 320 and N = 288 at max_level 3, t = 2, watch `nvidia-smi` | 3 |

(Candidate-146 re-measure runs were dropped 2026-08-18: the paper is rebuilt
on post-fix provenance only — §12.9.)

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
| FMAX-RF | 128 | 288–320 | 0.44–0.40 | fine — as the probe allows; say plainly in the paper that N = 384 is not reachable single-rank |
| FMAX-DS | 96 | 192 | 0.500 | domain ladder |
| FMAX-DL | 160 | 320 | 0.500 | domain ladder — only if the N = 320 probe passed; else drop and rely on DS + the search-tier L = 64 vs 128 comparison |
| FMAX-PF | 128 | 256 | 0.500 | pump-free twin (t_pump = 0), t = 64 |
| freefall companions | — | — | — | `manifest_freefall.json` cells on the stored stacks |

**Pre-registered acceptance (fix before launching):** fine–medium relative
difference ≤ 10 % on peak f_geo; f_ff ladder spread ≤ 10 %; full ray bundles
at the quoted launch; quote grid-to-grid **spread** as the error bar (no
order fit on f_geo — measure observed order on the composite constraint
norms instead). If AMR regridding noise spoils the ladder, fall back to the
bondi-proven unigrid methodology (max_level 0, three cell sizes,
`convergence_check.py`-style spread tables); unigrid ceiling on this node is
N = 256 (~70 GB).

---

## 7. Phase 4 — controls and re-measurements

1. **Canonical-only control (4 GPUs, ~1 day).** MAP-Elites under the *same
   gated objective* (`f_geo_max`), current pipeline, every per-lump exotic
   dial pinned to 0, 200 evals, same gates and grid. Replaces the pre-fix
   canonical bound; keeps the "phantom sector required" claim commensurable
   with the headline. Verify on the first evals that
   `recipe_scalar_field_signs` is all +1.
2. **CPU-only close-outs (no GPU):** ablation + endpoint-gauge +
   proper-time tables for the new champion from stored stacks;
   `cache_fidelity` report attached to every quoted f_geo; calibrate the
   post-load gate's composite thresholds from the accumulated composite
   columns (closes `DebugPreGPU.md` PG-3).

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
- The depth lineage appears as the Pareto strength-end exhibit (search-tier
  tables + the DEPTH-X production point), explicitly labelled as outside the
  gated objective.
- Every quoted number names its source run and pack; every pack regenerates
  its tables by script. New packs: `qball-trajectory-fgeo-max-refinement`
  (Phase 1), `qball-trajectory-hq-promotion` (Phases 2–3), the canonical
  control (Phase 4).

---

## 10. Schedule sketch (4 GPUs, dependency order)

| slot | work | GPUs | est. wall |
|---|---|---|---|
| now | Phase 0: bondi mm + sep12 ladders (running) | 4 | ~1 day |
| 1 | Phase 1: CMA-ES refinement under `f_geo_max` | 4 | ~1 day |
| 2 | Phase 2: FMAX-RM + DEPTH-X + DEPTH-A + probes | 4 | ~1 day |
| 3 | analysis: freeze champion, fix manifest grids | — | hours |
| 4 | Phase 3: RC, RF, DS, DL, PF (+freefall companions) | 4 | ~1.5–2 days |
| 5 | Phase 4: canonical-only control | 4 | ~1 day |
| 6 | Phase 5: NM-1 pilot (NM-2 only on green light) | 1–4 | 0.5–1 day |

Total ≈ 6–7 GPU-days after the bondi ladders clear. Chain each slot off the
previous (§12.2) rather than babysitting launches.

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

---

## 12. Exact launch commands

Everything runs from `grteclyn-wrapper/` unless stated. Conventions follow
[`results/bondi-dipole-runaway/LAUNCH.md`](../../results/bondi-dipole-runaway/LAUNCH.md)
and the wrapper [`README.md`](../../grteclyn-wrapper/README.md).

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
- Common env for every paper-tier (t ≥ 40) replay:
  `GRTECLYN_GEO_EMIT_INTERVAL=2 GRTECLYN_GEO_MAX_EMISSIONS=25
  GRTECLYN_METRIC_STACK_N_SPACE=257`.

### 12.2 Prepare now (no GPU needed) + chain off bondi

```bash
cd grteclyn-wrapper

# 1. Clone the promotion campaign for Phase 3 and edit it:
cp -r scripts/campaigns/promote/bicomplex_cmaes_v1 scripts/campaigns/promote/fgeo_max_cmaes_v1
#    campaign.env.sh : CAMPAIGN_NAME=fgeo_max_cmaes_v1
#                      LIVE_RUN=${GRTECLYN_ROOT}/runs/neuralspacetime/search/cma_es/qball_traj_fgeo_max_cmaes_v1
#                      FREEZE_ROOT=${GRTECLYN_ROOT}/runs/neuralspacetime/hq/sources/qball_traj_fgeo_max_cmaes_v1
#                      OBJECTIVE_MODE=f_geo_max
#                      RL_PUMP_STOP_TIME = the search value (pump on for the
#                      full run — copy `rl_pump_stop_time` from the champion's
#                      params.txt, do NOT keep the bicomplex default 4)
#    manifest.json   : ids FMAX-*, objective_mode f_geo_max,
#                      defaults stop_time 64, geo_emit_interval 2, geo_max_emissions 25,
#                      evolution_mpi_ranks 1 on EVERY cell (multi-rank is dead),
#                      ladder N 192/256/<probe result>, DS L96/N192, DL L160/N320 (probe-gated)
#    manifest_pumpfree.json / manifest_freefall.json : same edits, t_pump=0 twin + observer cells

# 2. Validate without GPUs:
DRY_RUN=1 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh --list

# 3. Chain Phase 1 to fire when the bondi ladders finish (detached):
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
  GRTRESNA_TIMEOUT=2400 GRTRESNA_RANKS=1 \
  STOP_TIME=26 \
  POSTLOAD_MAX_HAM_L2=3e-2 POSTLOAD_MAX_MOM_L2=1e-2 \
  bash scripts/campaigns/qball_trajectory/cmaes_run.sh \
  > ../runs/neuralspacetime/search/cma_es/qball_traj_fgeo_max_cmaes_v1_launch.log 2>&1 < /dev/null & disown
```

Notes: `STOP_TIME=26` mirrors the `fgeo_v1` physics window so the refinement
is stage-comparable to its QD parent (the depth twin used 32 — do not mix).
`WARM_START_TOP_K=1` centres on the eval-322 basin. `SEED=11` decorrelates
from the seed-7 noise every earlier campaign replayed. Monitor / resume /
stop:

```bash
grep -a "^\[optimize\] \(eval\|gen\)" ../runs/neuralspacetime/search/cma_es/qball_traj_fgeo_max_cmaes_v1_launch.log | tail
RESUME=1 RUN_NAME=qball_traj_fgeo_max_cmaes_v1 bash scripts/campaigns/qball_trajectory/cmaes_run.sh
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
  GRTECLYN_METRIC_STACK_N_SPACE=257 \
  .venv/bin/python scripts/campaigns/hq/replay_eval.py \
  ../results/qball-trajectory-cmaes-refinement/run/eval_000195 \
  --name depth195_hq_L128_N256_t64 --runs-dir ../runs/neuralspacetime/hq \
  --gpu 1 --n-full 256 --l-full 128 --max-level 3 --regrid-threshold 0.02 \
  --stop-time 64 --plot-interval 72 \
  --objective-mode f_geo_depth --evolving-geodesic \
  --grtresna-ranks 1 --grtresna-iterations 50 --grtresna-timeout 7200 \
  --consumer-radii 12 18 24 \
  > ../runs/neuralspacetime/hq/depth195_hq.launch.log 2>&1 < /dev/null & disown
```

DEPTH-A, the clean-branch fallback exhibit, identically configured (the
branch-B conditioning question then answers itself inside Phase 2):

```bash
setsid nohup /usr/bin/env \
  GRTECLYN_GEO_EMIT_INTERVAL=2 GRTECLYN_GEO_MAX_EMISSIONS=25 \
  GRTECLYN_METRIC_STACK_N_SPACE=257 \
  .venv/bin/python scripts/campaigns/hq/replay_eval.py \
  ../results/qball-trajectory-cmaes-refinement/run/eval_000185 \
  --name depth185_hq_L128_N256_t64 --runs-dir ../runs/neuralspacetime/hq \
  --gpu 2 --n-full 256 --l-full 128 --max-level 3 --regrid-threshold 0.02 \
  --stop-time 64 --plot-interval 72 \
  --objective-mode f_geo_depth --evolving-geodesic \
  --grtresna-ranks 1 --grtresna-iterations 50 --grtresna-timeout 7200 \
  --consumer-radii 12 18 24 \
  > ../runs/neuralspacetime/hq/depth185_hq.launch.log 2>&1 < /dev/null & disown
```

Memory probes for the Phase-3 fine legs (kill after t = 2; ~10 min each):

```bash
for N in 320 288; do
  .venv/bin/python scripts/campaigns/hq/replay_eval.py \
    "$(ls -d ../runs/neuralspacetime/hq/sources/qball_traj_fgeo_max_cmaes_v1/eval_* | head -1)" \
    --name mem_probe_N${N} --runs-dir ../runs/neuralspacetime/hq/probes \
    --gpu 3 --n-full ${N} --l-full 128 --max-level 3 --regrid-threshold 0.02 \
    --stop-time 2 --objective-mode f_geo_max \
    --grtresna-ranks 1 --grtresna-timeout 7200
done
# watch in a second shell: nvidia-smi --query-gpu=index,memory.used --format=csv -l 10
```

Largest N that stays comfortably under ~75 GB becomes FMAX-RF (and gates
FMAX-DL) in the manifest.

### 12.5 Phase 3 — the matrix

```bash
cd grteclyn-wrapper
GPU_ID=0 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-RC
GPU_ID=1 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-RF
GPU_ID=2 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-DS
GPU_ID=3 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-DL   # only if probe passed

# pump-free twin + freefall companions, same runner, different manifest:
MANIFEST=scripts/campaigns/promote/fgeo_max_cmaes_v1/manifest_pumpfree.json \
  GPU_ID=0 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh --list   # then launch its cell id
MANIFEST=scripts/campaigns/promote/fgeo_max_cmaes_v1/manifest_freefall.json \
  GPU_ID=1 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh --list   # then launch its cell ids
```

`DRY_RUN=1` any cell first; launch logs land in the campaign's
`VALIDATION_LAUNCH_LOG_DIR`.

### 12.6 Phase 4 — canonical-only control

`PIN_DIMS` **replaces** the launcher default, so the physics pins must be
restated alongside the five exotic pins:

```bash
cd grteclyn-wrapper
setsid nohup /usr/bin/env \
  QD_NAME=qball_traj_fgeo_canonical_v1 \
  QD_TARGET_EVALS=200 \
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
evals** that `params.txt` carries `recipe_scalar_field_signs = 1 1 1 1 1`
and the exotic integral is 0; abort and fix the pins if not.


### 12.7 Phase 5 — dressed-star pilot NM-1

Working reference for the override block:
`scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh` (running these exact
overrides in production right now). Pilot = FMAX-champ's trajectory, matter
swapped to self-gravitating ultraweak stars, pump off, tightened solve:

```bash
cd grteclyn-wrapper
mkdir -p ../runs/neuralspacetime/pilots
setsid nohup /usr/bin/env \
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
