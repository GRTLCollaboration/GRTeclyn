# FTL discovery campaigns (three-stage pipeline)

Matter-first FTL search in three **separate** stages. Each stage has its own launcher, resolution, and output directory. Stages share physics conventions but **not** the same grid or stop time.

```
Stage 0  QD (MAP-Elites)     scripts/campaigns/qd/run.sh
              │
              ▼  trajectory.jsonl + eval_* dirs
Stage 1  CMA-ES refine       scripts/campaigns/cmaes/run.sh   (optional)
              │
              ▼  warm-start from stage 0/1 trajectory
Stage 2  HQ promotion        scripts/campaigns/hq/run_batch.sh
              │
              ▼  runs/grtresna_promote/{campaign}_hq_eval*/
```

## Stage comparison

| | **QD** | **CMA-ES** | **HQ** |
|---|--------|------------|--------|
| Script | `qd/run.sh` | `cmaes/run.sh` | `hq/run_batch.sh` |
| Algorithm | MAP-Elites 8×8 archive | CMA-ES hill-climb | Replay + fresh GRTresna |
| Grid | N=128, L=64, ml=**2** | **same as QD** | N=**256**, L=**128**, ml=**3** |
| Stop time | **16** | **16** | **30** |
| Frames | off | off | **on** |
| Plotfile consumer | consume + delete, keep 3 | same | same (via `evaluation.py`) |
| 4D geodesic | `search` (fast) | `search` (fast) | `hq` (full verify) |
| Objective | `ftl_first` | `ftl_first` (default) | `ftl_first` (via replay) |
| Output | `runs/grtresna_qd/…` | `runs/grtresna_cmaes/…` | `runs/grtresna_promote/…` |
| Shared config | `lib/search_common.sh` | `lib/search_common.sh` | `lib/promote_common.sh` |

QD and CMA-ES **must** stay aligned (same `search_common.sh`) so warm-started CMA-ES scores are comparable to the QD trajectory.

HQ is **intentionally different**: higher resolution and longer time stress-test whether shortcuts survive refinement.

---

## Stage 0 — MAP-Elites (QD)

```bash
cd grteclyn-wrapper
QD_NAME=ftl_4d_v1 QD_TARGET_EVALS=200 \
  RUNS_DIR="${GRTECLYN_ROOT}/runs/grtresna_qd/ftl_4d" \
  nohup bash scripts/campaigns/qd/run.sh \
  > "${GRTECLYN_ROOT}/runs/grtresna_qd/ftl_4d/launch.log" 2>&1 &
```

Key env vars: `QD_NAME`, `QD_TARGET_EVALS`, `QD_RESUME=1`, `DESCRIPTOR_MODE=ftl_lifetime`, `GPU_IDS`, `SKIP_QD_PREFLIGHT_TESTS=1`.

Outputs: `trajectory.jsonl`, `archive.json`, `ftl_champions.json`, retained `eval_*/`.

---

## Stage 1 — CMA-ES (elite improvement)

Warm-start from a finished QD (or prior CMA-ES) trajectory:

```bash
cd grteclyn-wrapper
RUN_NAME=ftl_4d_cmaes_v1 \
WARM_START_TRAJECTORY="${GRTECLYN_ROOT}/runs/grtresna_qd/ftl_4d/ftl_4d_v1/trajectory.jsonl" \
WARM_START_TOP_K=8 SIGMA0=0.08 MAX_GENERATIONS=25 \
  nohup bash scripts/campaigns/cmaes/run.sh \
  > "${GRTECLYN_ROOT}/runs/grtresna_cmaes/ftl_4d_cmaes_v1/launch.log" 2>&1 &
```

Use `OBJECTIVE_MODE=robust_ftl` only if you deliberately want the v17-style persistence/exotic rebalancing (scores will **not** match QD `ftl_first` totals).

---

## Stage 2 — HQ promotion

Replays selected eval genomes at full resolution. **Not** the same run as QD/CMA-ES — launches separate jobs via `hq/replay_eval.py`.

```bash
cd grteclyn-wrapper
SOURCE_RUN="${GRTECLYN_ROOT}/runs/grtresna_qd/ftl_4d/ftl_4d_v1" \
CANDIDATES="156 0 142 1" \
NAME_PREFIX=ftl_4d \
  bash scripts/campaigns/hq/run_batch.sh
```

Candidate selection (one of):

- `CANDIDATES="eval gpu eval gpu …"` — explicit pairs
- `CANDIDATES_FILE=path/to/pairs.txt`
- `TOP_K=3 MIN_FTL_GEO_EVOL=0.03` — auto-pick top scores with 4D geo floor

`SOURCE_RUN` may point to a QD or CMA-ES campaign dir (must contain `eval_XXXXXX/` + `trajectory.jsonl` for `TOP_K`).

Run folders are named `{NAME_PREFIX}_hq_eval{eval_id}` (e.g. `ftl_4d_cmaes_hq_eval000144`), matching QD/CMA-ES campaign slugs instead of embedding `l128n256t30`. Set `INCLUDE_RESOLUTION_IN_NAME=1` for legacy `l128n256t30_{prefix}_hq_eval*` names.

Single-eval replay (manual):

```bash
uv run python scripts/campaigns/hq/replay_eval.py \
  runs/grtresna_qd/ftl_4d/ftl_4d_v1/eval_000156 \
  --name ftl_4d_hq_eval000156 \
  --gpu 0 --evolving-geodesic
```

Monitor incremental score: `tail -f runs/grtresna_promote/*/small_data/score_timeseries.jsonl`

---

## Shared libraries

| File | Used by |
|------|---------|
| `lib/bootstrap.sh` | QD + CMA-ES path setup |
| `lib/search_common.sh` | QD + CMA-ES grid, gates, 4D search profile, pytest preflight |
| `lib/promote_common.sh` | HQ batch + replay defaults |

---

## Typical `ftl_4d` workflow

1. **QD** until `QD_TARGET_EVALS` — note best eval + `f_geo_evol` in `trajectory.jsonl`
2. **CMA-ES** warm-started from that trajectory — same scores, local refinement
3. **HQ** promote top 1–3 from QD and/or CMA-ES — falsification at t=30, full 4D HQ trace

Analysis helpers (still under `scripts/search/`): `report_campaign_ftl.py`, `validate_tiers.py`.

---

## Smoke test (all stages)

Dry-run + preflight without launching full GPU jobs:

```bash
cd grteclyn-wrapper
bash scripts/campaigns/smoke_test.sh
```

Covers: pytest preflight (60 tests), QD dry-run, CMA-ES dry-run + warm-start from `ftl_4d_v1`, HQ batch dry-run + `TOP_K` auto-pick, `replay_eval.py --help`.
