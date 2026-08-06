# Campaign results: `qball_traj_fgeo_depth_v1` — pure geodesic-depth hunt

**Status:** complete — 200/200 evals, converged (`validation.json`: score stalled, front stability 1.0).
**Launched:** 2026-08-05 00:29 UTC · **Finished:** 2026-08-05 10:43 UTC (~10 h wall on GPUs 0–3).
**Best score:** 4609.29 (eval 7) — a **45.9 % evolving-geodesic path saving**, the deepest honest
shortcut measured in any campaign so far (v1's reward saturated at 20 %; its champion measured 24 %).
**The catch:** eval 7's genome is **bit-identical to the warm-start seed** (v1 eval 370). Two hundred
evaluations of mutation around the seven v1 elites never produced anything deeper than what they
inherited.

---

## 1. What this campaign searched for

Same 39-D MAP-Elites search as `qball_traj_fgeo_v1` (5 Q-ball lumps × orbit/tilt/spiral/exotic
dials + shared breathing/bobbing), with one deliberate change of objective:

- **Objective `f_geo_depth`:** score ≈ 10000 × raw `f_geo` — the fractional path saving of a light
  ray traced end-to-end through the *evolving* metric, best over the emission sweep
  (emit floor t ≥ 4). **No 20 % saturation, no multiplication by matter retention.** The v1
  campaign optimized retention because depth clipped at 20 %; this campaign pays depth alone.
  Honesty gates stay upstream: depth counts only if the ray bundle actually arrives
  (n_reached = n_rays) with small Hamiltonian drift, plus the usual post-load constraint gate.
- **Descriptors:** `ftl_peak_strength` × `ftl_lifetime_fraction`, 8 × 8 bins.
- **Warm start:** the seven kept v1 elites (evals 2, 141, 179, 226, 272, 322, 370).
- Exotic-matter penalty recorded but **weight 0**; pump-energy tax 40; dispersion gate 1.0 on the
  operational (coordinate) terms only.
- Physics window unchanged: stop_time 32, N=128 / L=64, max_level 1, pump on for the full run,
  measurement shell r = 8.

## 2. Headline numbers

| Metric | Value |
|---|---|
| Evals | 200 = 138 scored clean + 36 GPU crashes (exit 6) + 26 post-load rejections |
| Best score | **4609.29** — eval 7, set within the first two batches, never beaten in the remaining ~193 |
| Best depth `f_geo` | **0.459** (eval 7 = v1 seed 370 verbatim) |
| Best *new-genome* depth | 0.433 (eval 26, a real mutation of v1/226) |
| Archive | 9 elites, 14.1 % coverage; tier mix 6 operational + 3 constructed |
| Validation | 6 survivors, min tier "operational", converged = true |
| Per-metric champions | depth 0.459 (e7) · gated depth 0.350 (e14) · f_op 0.481 (e80) · superluminal fraction 1.0 (e5) · max local speed 4.2 × 10¹⁰ (e63 — coordinate artifact, see §3.5) |

## 3. Findings

### 3.1 The seed was already the local optimum
The all-time record is the inherited recipe itself: eval 7's genome has distance **0.000** to v1
eval 370. Nearby jitters re-scored almost as well (eval 5 ≈ v1/272 at distance 0.47 → depth 0.344),
and the best genuinely new genome (eval 26, distance 3.3 from v1/226) reached 0.433 — 94 % of the
record. Ten distinct top genomes cluster at depth 0.38–0.46. Conclusion: the v1 elite basin has a
depth ceiling of ~0.43–0.46, and 200 evals of local mutation could not climb out of it. A future
depth hunt needs bigger mutation steps or cold-started cells, not more of the same.

### 3.2 Depth is window-limited, not physics-limited
The champion's emission sweep is **monotonically rising to the last allowed emission**:

| t_emit | 0 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | **18** |
|---|---|---|---|---|---|---|---|---|---|---|
| f_geo | .202 | .230 | .259 | .291 | .323 | .353 | .383 | .412 | .437 | **.459** |

A ray launched at t = 18 crosses in 7.79 time units vs 14.4 flat — 54 % of the flat crossing time —
and every later launch would presumably do better, but stop_time = 32 ends the measurement. The
bundle is honest: 3/3 rays arrive, max Hamiltonian drift 1.0 × 10⁻³. **The single cheapest next
experiment is re-running the champion recipe with stop_time ~64** to find where depth actually
peaks.

### 3.3 Retention: recorded, not rewarded — behaved as designed
This objective deliberately pays raw depth only; matter retention is logged as a diagnostic, not a
score term. For the record: the deep elites hold 14–19 % of their lump matter at t = 32 (champion
18 %), and the depth channel stays open the whole run (lifetime fraction 1.0). This confirms the
design premise — v1's scores were retention-dominated, and freeing depth from that multiplier is
what let the measurement reach 45.9 %. These numbers are reference data for cross-campaign
comparisons, not a failure mode of this campaign.

### 3.4 The archive collapsed onto the max-lifetime row
8 of 9 elites sit in the top lifetime bin (y = 7, channel open ≥ the full window), spread across
peak-strength bins x = 0…7; the ninth is the degenerate no-FTL cell [0,0]. Depth and lifetime
co-select under this objective — short-lived deep channels never appeared. Coverage stalled at
14.1 % by mid-campaign.

### 3.5 Coordinate metrics remain artifacts
The `max_local_speed` champion (eval 63) logged 4.2 × 10¹⁰ c with perfect confinement fraction 1.0
but structural persistence 10⁻¹⁹¹ — a gauge blow-up, not transport. This is exactly why coordinate
channels are kept out of the objective and only logged as champions/descriptors.

### 3.6 Failure taxonomy
36/200 (18 %) died as GPU crashes (exit code 6) and 26/200 (13 %) were rejected by the post-load
constraint gate ("post-load constraints exceed threshold") — both concentrated in
large-breathing / large-z corners of the space. Roughly a third of the budget bought no signal;
tightening the sampler priors around known-viable amplitude ranges would buy back ~60 evals.

## 4. Champion configuration (eval 7 ≡ v1 eval 370)

Full parameters: [`run/eval_000007/params.txt`](run/eval_000007/params.txt); genome + all 39
overrides in [`run/trajectory.jsonl`](run/trajectory.jsonl) (eval 7). Character of the recipe:

- 5 lumps, **all exotic-leaning** (exotic dials 0.87, 0.52, 0.69, 0.61, 0.99), radii 5.9 / 2.8 /
  2.2 / 5.5 / 7.1.
- Strong vertical bobbing (z-amplitude 2.59, ω_z 0.80), almost **no breathing** (0.037), slow
  retrograde orbital rotations.
- Geodesic measurement: probe axis y, 3/3 rays arrive, best emission t = 18 → arrival t = 25.79
  (flat crossing 14.4).
- Score components: depth term ≈ 4589 of the 4609 total; retention 0.18; ANEC condition 1.0;
  constraint growth 0.936; exotic penalty −1.6 recorded at weight 0.

## 5. Caveats for the paper

- **Depth is measurement-window-clipped** (§3.2) — 45.9 % is a lower bound on this recipe's peak.
- **Retention is outside this objective** (§3.3): the 14–19 % end-of-run retention of the deep
  elites is a recorded diagnostic. Scores here are not comparable to v1 scores, which multiplied
  by retention — compare depths, not scores, across campaigns.
- The exotic penalty was weight-0: all deep recipes are heavily exotic; nothing here says a
  canonical-matter recipe can do this.
- Champions on coordinate metrics (local speed, superluminal fraction) are gauge artifacts kept
  for the record only.

## 6. What's in this pack

The full campaign directory is `runs/grtresna_qd/qball_traj_fgeo_depth_v1` on the produce node;
heavy binaries (531 MB `initial_data.gridinit` per eval, 16 MB metric stacks, 312 MB launch log,
pruned eval dirs) stay there. This pack keeps every number the analysis above uses.
Note: the campaign's keep-top-3 pruning deleted the survivor-front episode dirs (evals 49, 89,
109, 179) before finish; their scores and genomes survive in `trajectory.jsonl`.

| Path | Content |
|---|---|
| `run/trajectory.jsonl` | all 200 scoreboard rows: score, components, descriptors, full genome/overrides |
| `run/archive.json`, `run/validation.json`, `run/pre_gpu_archive.json` | final MAP-Elites state + convergence report + surrogate archive |
| `run/ftl_champions.json`, `run/ftl_retention.jsonl` | per-metric champions and crowning events |
| `run/metadata.json` | campaign config: search space, base overrides, warm-start, gates |
| `run/eval_000005 / 7 / 14 / 26 / 63 / 80` | kept champion dirs: `params.txt`, `score.json`, diagnostics `data/`, `small_data/` (incl. `evolving_geodesic.json` with the emission sweeps), per-eval `run.log` |
| `logs/qd_progress_evals_1-200.log` | `[qd]` extract of the launch log: batch reports, prunes, record announcements |

**Regenerator:** `grteclyn-wrapper/scripts/campaigns/qball_trajectory/run_fgeo_depth.sh`.

## 7. Run commands & properties

All commands run from `grteclyn-wrapper/` on the compute node (this machine — nothing to copy
anywhere). MPI is broken cluster-wide: everything runs single-rank (`RANKS=1`), never raise it.

**Launch** (detached, as actually run; output → `runs/grtresna_qd/<name>_launch.log`):

```bash
cd grteclyn-wrapper
setsid nohup bash scripts/campaigns/qball_trajectory/run_fgeo_depth.sh \
  > ../runs/grtresna_qd/qball_traj_fgeo_depth_v1_launch.log 2>&1 < /dev/null &
```

**Resume after an interruption** (same command with `QD_RESUME=1`):

```bash
QD_RESUME=1 setsid nohup bash scripts/campaigns/qball_trajectory/run_fgeo_depth.sh \
  > ../runs/grtresna_qd/qball_traj_fgeo_depth_v1_launch.log 2>&1 < /dev/null &
```

**Stop** (the one sanctioned way — kills the orchestrator first, then workers, then verifies):

```bash
bash scripts/campaigns/stop_campaign.sh ../runs/grtresna_qd/qball_traj_fgeo_depth_v1
pgrep -af qball_traj_fgeo_depth_v1   # must come back empty
```

**Monitor progress** (batch reports and record announcements):

```bash
grep '\[qd\]' ../runs/_logs/qball_traj_fgeo_depth_v1.log | tail -5
```

**Properties the launcher pins** (everything else inherits `run.sh` defaults):

| Setting | Value | Meaning |
|---|---|---|
| `QD_NAME` | `qball_traj_fgeo_depth_v1` | run dir `runs/grtresna_qd/<name>` |
| `OBJECTIVE_MODE` | `f_geo_depth` | raw depth is the only first-order reward (1 % saving = 100 pts), uncapped, not survival-multiplied |
| `SCORE_EXOTIC_PENALTY_WEIGHT` | `0` | phantom lumps are free fuel |
| `STOP_TIME` | `32` | up from v1's 26 — gives the implosion room |
| `GRTECLYN_GEO_MAX_EMISSIONS` | `10` | emission sweep t = 0…18 (up from v1's 7 launches) |
| `SEED_EVAL_DIRS` | the 7 kept `qball_traj_fgeo_v1` eval dirs | warm start, re-scored under this objective |
| `GPU_IDS` / `RANKS` | `0 1 2 3` / `1` | 4 GPU slots, single-rank solves |
| From `run.sh` | `DESCRIPTOR_MODE=ftl_lifetime`, `GEODESIC_EMIT_MIN_TIME=4`, `SCORE_FTL_DISPERSION_GATE=1.0`, `SCORE_PUMP_ENERGY_WEIGHT=40`, 200 evals, batch 4, 8×8 bins | |

**Replay the champion** (next-step 1 — same recipe, doubled window; adjust `--gpu` to a free slot):

```bash
PYTHONPATH=src .venv/bin/python scripts/campaigns/hq/replay_eval.py \
  ../runs/grtresna_qd/qball_traj_fgeo_depth_v1/eval_000007 \
  --name depth_champ_t64 --runs-dir ../runs/depth_champion_longrun \
  --gpu 1 --stop-time 64 --grtresna-ranks 1
```

## 8. Recommended next steps

1. **Champion long-run:** replay eval 7's recipe at stop_time 64 with a dense emission sweep —
   locate the true depth peak (§3.2). One GPU, a few hours.
2. **Escape the basin:** larger mutation scale or a fraction of cold-random individuals; the
   current search never improved on its own seed (§3.1).
