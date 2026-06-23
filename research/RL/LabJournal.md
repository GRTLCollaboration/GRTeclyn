# RL Campaign — Lab Journal (v1)

Reverse-chronological log for the closed-loop RadialRecipe RL stack. Newest first.

---

## Current state (2026-06-23)

### Gates passed

| Gate | What it proves | How to run | Artifact |
|------|----------------|------------|----------|
| **0A** | IVP neutrality: `rl_enabled=0` matches `rl_enabled=1` + neutral agent (rel err 1e-10 on L2_Ham / min_lapse) | `bash grteclyn-wrapper/scripts/campaigns/rl/gate0a_run.sh` | `runs/rl_gate0a/spacetime_splash_v14_eval010` |
| **0B** | ZMQ bridge smoke: dummy agent, no hang, obs_dim=14 | `GATE0B_STOP_TIME=1.0 bash …/gate0b_run.sh` | `runs/rl_gate0b/spacetime_splash_v14_eval010` |
| **1** | Tax Man T1–T4 + live `SpacetimeFtlEnv` (plot consumer, score_timeseries, clean exit) | `bash …/gate1_run.sh` | `runs/rl_gate1/spacetime_splash_v14_eval010` |

**Not yet run:** Gate 2 (kamikaze / actuation proof), PPO dry run.

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

### 1. Gate 2 — Kamikaze / actuation proof

Run non-neutral actions (random or max-amplitude) and confirm:

- Pump visibly perturbs the scalar field (not IVP no-op path)
- Reward responds to actions (not flat −21 every step)
- Governor / fences fire under stress (horizon, L2_Ham threshold)
- Sim stays stable or terminates gracefully (no NaN / segfault)

### 2. `ftl_geo` NaN triage

Headline FTL reward term uses Weyl4 extraction; historically NaN in splash campaigns. Without finite `f_geo`, PPO can only learn “don’t crash.” Fix extraction or confirm it works on boson chassis at RL grid scale.

### 3. `train_motor.py` dry run

After Gate 2 + live `ftl_geo`: short PPO (e.g. 50k steps, 1 GPU) — verify gradients, VecNormalize, return > random.

### 4. Boson-star `general_ftl` chassis (original plan item)

Generate elite initial data tuned for the RL control regime (not just QD splash elite reused as-is).

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
