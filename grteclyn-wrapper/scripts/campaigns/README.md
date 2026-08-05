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
| Plotfile scratch | node-local `/tmp` | node-local `/tmp` | node-local `/tmp` |
| 4D geodesic | `search` (fast) | `search` (fast) | `hq` (full verify) |
| Objective | `ftl_first` | `ftl_first` (default) | `ftl_first` (via replay) |
| Output | `runs/grtresna_qd/…` | `runs/grtresna_cmaes/…` | `runs/grtresna_promote/…` |
| Shared config | `lib/search_common.sh` | `lib/search_common.sh` | `lib/promote_common.sh` |

QD and CMA-ES **must** stay aligned (same `search_common.sh`) so warm-started CMA-ES scores are comparable to the QD trajectory.

HQ is **intentionally different**: higher resolution and longer time stress-test whether shortcuts survive refinement.

---

## Stopping a campaign — one tool, one way

```bash
# Preview what would be killed (touches nothing):
bash scripts/campaigns/stop_campaign.sh --dry-run <runs_dir | campaign_name>
# Stop for real, with verification and escalation:
bash scripts/campaigns/stop_campaign.sh <runs_dir | campaign_name> [...]
```

Works for **every** campaign type (QD, CMA-ES, HQ replays, Bondi matrices,
one-off ladders). Names resolve under `runs/grtresna_qd/<name>` then
`runs/<name>`. There are deliberately **no per-campaign stop scripts**.

Do NOT stop a campaign by hand with the pid captured at launch or by
pattern-killing workers. Three failure modes, all hit in practice
(2026-08-05, `bondi_dipole_v1` post-mortem — full details in the
`stop_campaign.sh` header):

1. `$!` after `setsid nohup ... &` is the **dead setsid parent**, not the
   session-leader launcher — killing it (or its group) hits nothing.
2. Killing workers (evolution binary, consumer) looks like a **finished step**
   to the orchestrator, which then launches the next eval — the campaign
   appears to refuse to die.
3. GRTresna solvers detach into their **own session/pgid** and survive
   group-kills of the launcher.

The tool therefore freezes the queue first (orchestrators + their shell
ancestors, found via `launcher.pid`, driver argv and a parent-walk), then
sweeps workers by runs dir and scratch path, then **verifies with pgrep and
escalates to SIGKILL** until nothing survives.

**Launcher contract** — every launcher registers its true PID once the runs
dir is known (one line; `exec` into a generic runner preserves it):

```bash
source "${SCRIPT_DIR}/../lib/launcher_common.sh"
campaign_register_launcher "${RUNS_DIR}"
```

Discovery covers unregistered launchers too, but the pid file is the fast,
unambiguous path. Add the call to any new launcher.

---

## Plotfile scratch is node-local — mandatory for every stage

**Plotfiles must never be written to NFS.** They are write-once, read-once,
delete-immediately transients. `amr.plot_file` / `amr.check_file` are
independent of `output_path`, so send the heavy data to node-local NVMe and
keep only results on the shared filesystem.

```
GPU sim ──plotfiles──▶ /tmp/<scratch>/<run>/RadialRecipePlt*   (local NVMe)
   │                        │ consumer extracts (~16 s), then deletes
   │                        ▼
   └──.dat, KB/s──▶ <output_path>/data/  +  small_data/         (NFS, ~12 MB)
```

Every params.txt carries:

```
output_path    = "<RUNS_DIR>/<run>"                        # NFS: .dat + small_data
amr.plot_file  = "/tmp/grteclyn_scratch/<run>/RadialRecipePlt"
amr.check_file = "/tmp/grteclyn_scratch/<run>/RadialRecipeChk"
```

and the consumer runs with `--data /tmp/grteclyn_scratch/<run>`
`--out <RUNS_DIR>/<run>/small_data`.

**This is now the default, for every stage.** `core/scratch.py` applies it in
`episode_path_overrides`, so QD, CMA-ES, HQ and the post-load gate all get it
without a launcher change; the consumer, the scoring plotfile lookup and the
cleanup paths follow the same mapping. `GRTECLYN_SCRATCH=/other/path` moves the
root, `GRTECLYN_SCRATCH=0` restores the old episode-directory layout, and an
unwritable root falls back to the episode directory with a warning rather than
killing the campaign. A launcher that sets `amr.plot_file` explicitly (the RL
gate scripts, the pump queues) still wins — an explicit override is never
overwritten.

**Why.** A 256³ ml=3 plotfile is ~3.2 GB written and read back every ~288 s
per run. On NFS that capped concurrency at **2 runs** — consumers blocked in
NFS I/O (`D` state), backlogs grew, plotfiles accumulated. On node-local NVMe
the same node sustains **6 concurrent HQ runs**; NFS traffic drops from
~130 MB/s to KB/s.

Measured, mode-0 pump ladder, 2026-07-28 — 6 × (N=256, L=128, t=30) on GPUs 0–5:

| | |
|---|---|
| Local scratch per run | 8.8 GB (3 plotfiles resident) |
| NFS run dir per run | ~12 MB (~30 MB at t=30; `metric_stack` dominates) |
| Extraction | 15.7 s vs 288 s cadence (18× headroom) |
| Plotfiles leaked to NFS | 0 |
| Overlay usage | 53 GB of 1.1 TB free |

### Rules

1. **Cloning a baseline `params.txt`** — `amr.plot_file` / `amr.check_file` are
   absolute paths into the *source* run dir. Strip and re-emit them, or the
   clone writes into the baseline and a `--delete` consumer prunes the
   baseline's own plotfiles. Dry-run and assert each of `output_path`,
   `amr.plot_file`, `amr.check_file` occurs exactly once.
2. **Purge scratch only when everything was extracted** — scratch is transient,
   so an unextracted plotfile deleted from `/tmp` is lost (on NFS it would have
   survived). Compare resident plotfiles against
   `small_data/consume_state.json`; if any are missing, keep the directory and
   log it.
3. **Drain for 600 s** after the simulation exits before killing the consumer.
   120–180 s truncated confinement data in three fast-ladder runs.
4. **Budget the disk** — `n_runs × (keep_last + jobs) × plotfile_size`. Check
   `df -h /tmp`. Checkpoints are **not** pruned by the consumer; with
   `checkpoint_interval > 0` they accumulate on scratch.
5. **Never run an external deletion loop** alongside the consumer — see the
   [race-condition post-mortem](../../README.md#plotfile-pruning-failure-mode-and-fix).

### Keep writes inside your own directories

Call the venv python directly instead of `uv run` (which writes to
`~/.cache/uv`), and pin library caches into scratch:

```python
cmd = [str(WRAPPER / ".venv/bin/python"), "-m",
       "grteclyn_wrapper.visualisation.process_wave.consume_plotfiles", ...]
env["XDG_CACHE_HOME"]      = f"{SCRATCH}/_cache"
env["UV_CACHE_DIR"]        = f"{SCRATCH}/_cache/uv"
env["MPLCONFIGDIR"]        = f"{SCRATCH}/_cache/mpl"
env["TMPDIR"]              = f"{SCRATCH}/_cache/tmp"
env["PYTHONPYCACHEPREFIX"] = f"{SCRATCH}/_cache/pyc"
```

On the shared cluster `~/.local/bin` and `~/.local/lib` are owned by
`nobody:nogroup` and are not writable — admin policy is write only to your own
directories. The complete set of write targets for a campaign is
`/tmp/<scratch>/` plus the NFS run directory. Nothing else.

`core/scratch.py` pins the same variables for every child process the wrapper
spawns (`cache_env`), filling in only the ones the operator has not already
set. Reference launcher implementation, predating the central one:
[`rl/pump_ladder_queue.py`](rl/pump_ladder_queue.py) (mode-0 pump ladder).

---

### End-to-end orchestrator (QD → CMA-ES → HQ)

Single script under one campaign root (`runs/campaigns/<CAMPAIGN_NAME>/`):

```bash
cd grteclyn-wrapper
CAMPAIGN_NAME=general_ftl_wormhole_v22 \
  bash scripts/campaigns/run_full_campaign.sh
```

Stages (sequential, blocking):

1. **QD** → `qd/` (200 evals, keep top 3 + FTL champions)
2. **CMA-ES** → `cmaes/` (warm-start **top-1** QD eval, keep top 3)
3. **HQ** → `promote/` (replay **top-1** CMA-ES eval, foreground)

Resume / partial:

```bash
RESUME=1 STAGE=cmaes CAMPAIGN_NAME=... bash scripts/campaigns/run_full_campaign.sh
```

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
| `lib/launcher_common.sh` | every launcher — writes `launcher.pid` for `stop_campaign.sh` |
| `lib/search_common.sh` | QD + CMA-ES grid, gates, 4D search profile, pytest preflight |
| `lib/general_ftl_pins.sh` | v20 wormhole / ring / spin `--pin-dimension` bundles |
| `lib/pipeline_monitor.sh` | Optional `PIPELINE_MONITOR=1` GPU/pipeline sampling |
| `lib/promote_common.sh` | HQ batch + replay defaults |

---

## v20 `general_ftl` (wormhole / ring / spin)

Orchestrator only — still launches **`qd/run.sh`** (stage 0), never a parallel code path.

```bash
cd grteclyn-wrapper
MODE=par QD_ITERATIONS=30 bash scripts/campaigns/general_ftl/run_all.sh
```

Single branch:

```bash
BRANCH=wormhole GPU_IDS="0 1 2 3" QD_TARGET_EVALS=80 bash scripts/campaigns/general_ftl/run_all.sh
```

Multi-slot pipeline test (8 evals, 4 GPUs, 2 concurrent evolutions per GPU):

```bash
BRANCH=wormhole PIPELINE_MONITOR=1 \
  QD_TARGET_EVALS=8 GPU_IDS="0 1 2 3" GPU_SLOTS_PER_DEVICE=2 BATCH_SIZE=8 \
  STOP_TIME=4.0 PLOT_INTERVAL=40 QD_ITERATIONS=4 SKIP_QD_PREFLIGHT_TESTS=1 \
  bash scripts/campaigns/general_ftl/run_all.sh
```

**v21** — wormhole relaunch, 5 concurrent GPU evolutions per device (pipelined MAP-Elites):

```bash
BRANCH=wormhole PIPELINE_MONITOR=1 \
  QD_NAME=general_ftl_wormhole_v21 QD_TARGET_EVALS=80 \
  GPU_IDS="0 1 2 3 4 5 6 7" GPU_SLOTS_PER_DEVICE=5 BATCH_SIZE=40 \
  bash scripts/campaigns/general_ftl/run_all.sh
```

See [MapElites.md v21](../../../../research/neuralspacetime/MapElites.md#v21-multi-slot-gpu-pipeline-5-evols-per-gpu-2026-06-17) for VRAM sizing notes.

`PIPELINE_MONITOR=1` writes `runs/_logs/${QD_NAME}.pipeline_monitor.csv` and `.pipeline_summary.txt`.
Preflight smokes stay in `scripts/campaigns/smoke_test.sh` — not under `general_ftl/`.

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
