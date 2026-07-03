# RL Campaign — Lab Journal (v1)

Reverse-chronological log for the closed-loop RadialRecipe RL stack. Newest first.

---

## 2026-07-03 — `gw_beam_qd100_v4` MAP-Elites (100 evals, COMPLETE)

First GW-beam QD run with **v4 hard gates** (Ψ₄ truncation at spike, `gw_health_multiplier`, archive tier gate). Fresh campaign; v2/v3 run dirs removed.

### Config

| Setting | Value |
|---------|-------|
| Objective / descriptor | `gw_beam` |
| Box | L=64, N=128, t_stop=24 |
| Ψ₄ extraction | **R=20** (was R=12 in v3) |
| Matter | 5 canonical Q-balls, trajectory ansatz |
| Search clamp | `trajectory_v_rad_inward_floor=-0.1` |
| Retention | top-3 by score + in-flight pipeline (~25–30 dirs on disk) |

Run dir: `runs/grtresna_qd/gw_beam_qd100_v4/`

### Outcome

| Metric | Value |
|--------|-------|
| Evals completed | 100/100 |
| Healthy (score > 0) | **22** |
| Collapse vetoed (`mult=0`, ~−116) | **77/100** (77%) |
| Archive elites | **5** (rejected tiers excluded) |
| Best score | **eval 88 @ 3.09** |

### Top 5

| Eval | Score | Tier | beam_ratio | Notes |
|------|-------|------|------------|-------|
| 88 | 3.09 | constructed | 0.14 | Campaign best; full t=24, max Ham≈0.08 |
| 61 | 2.82 | constructed | **0.28** | Best **directional** GW in top tier |
| 52 | 2.78 | constructed | 0.27 | |
| 55 | 2.71 | constructed | 0.25 | |
| 49 | 2.70 | constructed | 0.18 | |

### Findings

1. **Hard gates worked.** v3 collapse exploit (eval 51 @ 336) is gone — collapse runs score ~**−116** (`mult=0`, spike penalty −112.5), not +300.
2. **No strong GW emitter found.** Best scores are modest (~2–3), dominated by survival/health bonuses; mean Ψ₄ power remains ~10⁻⁴–10⁻³. Eval 88 wins overall score but eval **61 has higher beam ratio** (28% vs 14%).
3. **Search still collapse-heavy early.** 77% of evals hit constraint spike before t≈18 (R=20 light-travel gate). Optimizer is probing aggressive configs; 22 stable survivors is a reasonable yield.
4. **Archive is trustworthy.** Only `tier ≥ constructed` runs occupy cells; collapse modes stay in trajectory/near-miss pool only.

### Next

- Inspect eval 61/88 Ψ₄ + constraint plots before trusting any as physical emitters.
- Consider tighter inward-motion limits, longer stop_time, or CMA-ES warm-start from eval 61 (beam) or 88 (score).

---

## Current state (2026-06-23)

### Gates passed

| Gate | What it proves | How to run | Artifact |
|------|----------------|------------|----------|
| **0A** | IVP neutrality: `rl_enabled=0` matches `rl_enabled=1` + neutral agent (rel err 1e-10 on L2_Ham / min_lapse) | `bash grteclyn-wrapper/scripts/campaigns/rl/gate0a_run.sh` | `runs/rl_gate0a/spacetime_splash_v14_eval010` |
| **0B** | ZMQ bridge smoke: dummy agent, no hang, obs_dim=14 | `GATE0B_STOP_TIME=1.0 bash …/gate0b_run.sh` | `runs/rl_gate0b/spacetime_splash_v14_eval010` |
| **1** | Tax Man T1–T4 + live `SpacetimeFtlEnv` (plot consumer, score_timeseries, clean exit) | `bash …/gate1_run.sh` | `runs/rl_gate1/spacetime_splash_v14_eval010` |
| **2** | Kamikaze actuation proof: multi-lump pump steers field, reward responds, no NaN | `bash …/gate2_run.sh` | `runs/rl_gate2/gate2_kamikaze` |

**Not yet run:** PPO dry run on FTL-producing elite.

### Chassis

- **Elite IC:** `runs/grtresna_qd/spacetime_splash_v14_moving/eval_000010` (boson `grtresna_complex_scalar`, gridinit present)
- **Binary:** `Examples/RadialRecipe/main3d.gnu.CUDA.ex` (`USE_RL=TRUE`, links local libzmq)
- **Build:** `make USE_RL=TRUE USE_CUDA=TRUE USE_MPI=FALSE CUDA_ARCH=90 NVCC_CCBIN=/usr/bin/g++`
- **ZMQ:** `GRTeclyn/local/zeromq/` (no sudo); `LD_LIBRARY_PATH` must include it at runtime

### Observation / action layout (N lumps)

| | Dim | Content |
|---|-----|---------|
| **Obs** | 6 + 8N | global: min_χ, min_lapse, max\|K\|, L2_Ham, L2_Mom, t; per lump: x,y,z, size, mass, peak, min_lapse, min_χ |
| **Act** | 3N + 2 | per lump: amp, freq, phase (raw ∈ [-1,1]); gauge: lapse_advec, shift_Γ |

C++ (`RLActionApplier`) is the sole raw→physical mapper; Python only clips to [-1, 1].

---

## 2026-06-23 — Gate 2 PASS: multi-lump actuation proof

### Results (`gate2_run.sh` / `gate2_validate.py`)

| Check | Result |
|-------|--------|
| T0 multi-lump dims | PASS — obs=30, act=11, num_lumps=3 |
| T1 episodes ran | PASS — neutral=12 steps, kamikaze=12 steps, no crash |
| T2 pump perturbs field | PASS — L2_Ham diff=3.64e-3 (kamikaze 0.0083 vs neutral 0.0046) |
| T3 reward responds | PASS — total reward diff=5.95 (kamikaze penalized more) |
| T4 governor | PASS — kamikaze L2_Ham grew faster (fence not triggered in short 0.5s run) |
| T5 no NaN | PASS |

**Key findings:**
- 3 pump sites at (8,8,8), (6,8,8), (10,8,8) — multi-site spatial targeting works
- EMA ramp visible: kamikaze L2_Ham reverses decay at step 5, grows exponentially by step 11
- Lump seed bug fixed: was at (0,0,0) (domain corner, envelope=0.25%), moved to grid center
- GPU build confirmed: `main3d.gnu.CUDA.ex` with RLBridge fix, 2211 code units/h on H100
- `f_geo = 0` everywhere — splash IC has no FTL geometry (expected; need FTL elite)

### Reward signal finding

`score_timeseries.jsonl` flows incrementally with per-frame `f_geo`, `operational_ftl_geodesic`,
`max_local_speed`, `shift_drive`. All zero because splash IC does not produce FTL.
With FTL-producing IC these become the mid-episode dense reward signal.
`ftl_geo_evolving` (4D null trace) only resolves end-of-run — serves as terminal audit bonus.

---

## 2026-06-23 — Gates 0A→1 complete; ZMQ hardened

### What was built (commits `3a505c7` → `d64e0e6`)

**C++ (GRTeclynCore/RL + RadialRecipe)**

- Matter pump on `ComplexScalarField` (multi-site, tanh governor, per-lump amp/freq/phase)
- Per-lump tracker + observation collector; `RLBridge` ZMQ REP
- Handshake: REP `recv(hello)` before first `send(obs)` (`0a26d64`)
- **EFSM fix:** track `m_awaiting_reply`; on action timeout consume late reply before next `send(obs)`; terminate after 3 consecutive timeouts (`d64e0e6`)
- Build fixes: `FillPatch` qualification, nvlink ODR, GPU governor precompute, `-lzmq` rpath, `NVCC_CCBIN=/usr/bin/g++`

**Python (`grteclyn_wrapper/rl/`)**

- `SpacetimeFtlEnv` — Gymnasium env, Tax Man dense reward + fences + end-of-episode audit
- `SubprocessEpisodeLauncher` — direct exe (no mpirun), plot consumer, conda libstdc++ stripped from sim env
- `dummy_agent.py`, `zmq_client.py`, `train_motor.py` (SB3 PPO stub), `params.py`, tests under `tests/rl/`
- Gate scripts: `gate0a_run.sh`, `gate0a_compare.py`, `gate0b_run.sh`, `gate1_run.sh`, `gate1_validate.py`

### What was tested

| Layer | Result |
|-------|--------|
| Python unit tests | 449 passed, 4 skipped (`uv run pytest`) |
| Gate 0A | PASS — t=2 and t=16, zero rel err on tracked scalars |
| Gate 0B | PASS — multiple ZMQ exchanges, obs_dim=14 |
| Gate 1 static | PASS — T1 γ=0.999; T2 audit clip −2000; T3 drain-before-audit; T4 WEC None-guard |
| Gate 1 live | PASS — 8 env steps, `score_timeseries.jsonl` written, sim finalizes cleanly |
| Prior blocker | RESOLVED — ZMQ EFSM (`Operation cannot be accomplished in current state`) from REP desync after action timeout |

Gate 1 uses a **small grid** override (N=32, L=16, max_level=0, rl_coarse_step_interval=4, rl_zmq_timeout_ms=300000) for fast smoke; production training will use full elite grid params.

---

## Next steps (priority order)

### 1. Switch to FTL-producing elite initial data

Current splash boson star (`spacetime_splash_v14_moving/eval_000010`) is a collapse/oscillation
IC — it produces **zero FTL geometry** (`max_local_speed ≈ 0.997`, `f_geo = 0` at all frames).
All per-frame FTL proxies (`f_geo`, `operational_ftl_geodesic`, `shift_drive`) are identically zero.

**Needed:** promote from `boson_shell_ftl_rl_v1` QD (or reuse a proven `general_ftl_wormhole` elite)
that produces non-zero `f_geo` per-frame. The agent's pump should then steer the field to
*maintain or extend* existing FTL geometry.

### 2. Reward signal architecture (resolved understanding)

| Signal | When available | Role in training |
|--------|---------------|------------------|
| `f_geo` / `operational_ftl_geodesic` | **Per-frame** (incremental from plot consumer) | Dense reward via `compute_dense_reward` — `1000 × ftl_term` |
| `max_local_speed` | Per-frame | Shaping proxy (coordinate speed) |
| `shift_drive` | Per-frame | Shaping proxy (shift contribution) |
| `ftl_geo_evolving` (4D null trace) | **End-of-run** (recomputes over all saved frames) | Terminal audit bonus |
| `-500 × L2_Ham`, `-50 × lapse_penalty` | Per-step (from C++ obs directly) | Defensive dense shaping |

Key insight: with the RIGHT initial data, `f_geo` flows per-frame and gives the agent
a positive mid-episode gradient. It is only zero now because splash IC has no FTL.
The 4D trace (`ftl_geo_evolving`) is the authoritative end-of-run validation, not
the training signal.

### 3. `train_motor.py` dry run on FTL elite

Once FTL-producing IC is available: short PPO (50k steps, 1 GPU) — verify:
- `f_geo` is non-zero in dense reward
- Gradients flow, VecNormalize works
- Return improves over random baseline
- Agent learns to pump in a way that maintains/extends `f_geo`

### 4. Boson-star `general_ftl` chassis QD

Generate elite ICs tuned for RL control (not QD reuse).

**In progress:** MAP-Elites QD `boson_shell_ftl_rl_v1` — boson shell + `general_ftl` (see below).

---

## Bosonic matter vs scalar FTL QD

Historical wormhole MAP-Elites (`general_ftl_wormhole_v*`) used **real scalar lumps**
(`GRTRESNA_MATTER_SECTOR=scalar`). RL needs **`grtresna_complex_scalar`** (U(1) boson).

| Campaign | Launcher | Matter | Objective | Use |
|----------|----------|--------|-----------|-----|
| Wormhole FTL | `general_ftl/run_all.sh` | real scalar shell | `general_ftl` | warp motors (not RL pump) |
| Splash / collapse | `splash/run.sh` | **boson shell** | `critical_collapse` | geometry splash ([grlab](../grlab/README.md)) |
| **Boson FTL (RL)** | **`boson_star/ftl_shell_run.sh`** | **boson shell** | **`general_ftl`** | RL chassis + `f_geo` elites |
| Centered boson FTL | `boson_star/run.sh` | centered 7-D boson | `ftl_first` | pipeline smoke ([MapElites](../neuralspacetime/MapElites.md)) |

Interim RL gates reuse splash elite `spacetime_splash_v14_moving/eval_000010` (collapse QD, not FTL).
Production RL should promote from **`runs/grtresna_qd/boson_shell_ftl_rl_v1/`** once QD completes.

### Launch — boson shell FTL for RL (200 evals, frames, t=16, 6 plots)

```bash
cd grteclyn-wrapper
QD_NAME=boson_shell_ftl_rl_v1 \
QD_TARGET_EVALS=200 \
QD_ITERATIONS=30 \
STOP_TIME=16.0 \
PLOT_INTERVAL=320 \
GRTECLYN_FRAMES=1 \
GPU_IDS="0 1 2 3 4 5 6 7" \
GPU_SLOTS_PER_DEVICE=1 \
MAX_CONCURRENT_GRTRESNA=5 \
BATCH_SIZE=8 \
  nohup bash scripts/campaigns/boson_star/ftl_shell_run.sh \
  > ../runs/boson_shell_ftl_rl_v1.launch.log 2>&1 &
```

Monitor: `tail -f runs/grtresna_qd/boson_shell_ftl_rl_v1/trajectory.jsonl`

---

## Quick reference

```bash
# Build
cd Examples/RadialRecipe
make USE_RL=TRUE USE_CUDA=TRUE USE_MPI=FALSE CUDA_ARCH=90 NVCC_CCBIN=/usr/bin/g++ -j$(nproc)

# Gates
bash grteclyn-wrapper/scripts/campaigns/rl/gate0a_run.sh
GATE0B_STOP_TIME=1.0 bash grteclyn-wrapper/scripts/campaigns/rl/gate0b_run.sh
bash grteclyn-wrapper/scripts/campaigns/rl/gate1_run.sh

# Tests
cd grteclyn-wrapper && uv sync --extra rl && uv run pytest -q
```

**Branch:** `feature/interstellar`  
**Plan:** `research/RL/implementation-plan.md`
