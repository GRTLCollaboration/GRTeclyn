# RotatingWormholeCollapse — Debug notes

Investigation and fix for the crashing Teo weak-spin run (`weak_spin_effective`).
**Status: root cause confirmed; stable production path validated (config only, no C++ changes).**

## TL;DR

- The original crash is **not** bad Teo physics in the `.gridinit`. The file is
  constraint-satisfying at its native resolution.
- It fails when AMR refines **below** the file resolution: `ExternalGridInitialData`
  fills every level by **trilinear** interpolation, so fine levels get C0 kinks and
  CCZ4 blows up at the throat (checkerboard → NaN in `h11`).
- **Fix in practice:** use a `dx = 0.25` gridinit and evolve on a grid that does
  **not** go finer than the file. The validated setup is **unigrid**
  `N1=N2=256`, `N3=128`, `max_level=0` (`params_rotating_unigrid_dx025.txt`).
- **`max_level=1` with `N1=128` still NaNs** at t ≈ 0.17 on level 1 (past the old
  t = 0.13 crash). Do not use that AMR config until evolution/interpolation is
  improved.

**Governing rule:** finest evolved `dx` must be **≥** file `dx` (no sub-file
interpolation on refined levels).

---

## Update (4 June 2026): controls comparison + extended run

**Key finding: the unigrid weak-spin run does not carry a usable spin signal.**
By `t = 0.2` it is dominated by a *spin-independent* constraint blow-up, identical
to the spin-0 baseline.

### Three-way comparison at `t = 0.2` (raw plotfile data, full resolution)

| Run | `K` RMS | Ham L2 (t=0.005 → 0.2) | Verdict |
|-----|---------|------------------------|---------|
| `no_matter` (vacuum) | 9.8e-4 | 0.007 → **0.009** (flat) | clean & stable |
| `weak_spin` (Teo source) | 9.73e-2 | 0.005 → **1.79** | blown up |
| `a0` (spin-0 baseline) | 9.73e-2 | 0.006 → **1.79** | blown up *identically* |

Direct difference test (the criterion for claiming dynamics):

```
RMS(K_spin)            = 9.7323e-2
RMS(K_a0)              = 9.7318e-2
RMS(K_spin − K_a0)     = 8.1e-5     ← residual spin signal = 0.08% of field
RMS(K_spin − K_nomatter) = 9.7323e-2
```

The spin signal does **not** survive subtracting the `a=0` baseline: ~99.9% of the
`weak_spin` field is the same instability present without spin. The "strange"
blocky `K` frames are **not** a visualisation bug — they faithfully render a real
domain-filling constraint violation (per-MPI-box shading of an O(0.1) field). The
vacuum control renders clean because it stays O(1e-3).

### Extended restart (`Chk00040`, `stop_time = 5.0`)

Restarted cleanly from step 40 (`t = 0.2`) and continued stepping, but diverged
further and **NaN-aborted at `t ≈ 0.9`** (`nan_check = 1`): Ham L2 grew
`1.79 (t=0.2) → ~8.0 (t=0.74)`, roughly linearly. Extending this config produces a
longer movie of a spin-independent blow-up, not wormhole spin dynamics.

### Operational note — why disk fills fast

`params_rotating_unigrid_dx025.txt` uses `plot_interval = 1` and each plotfile is
**707 MB** (256×256×128 × 12 vars). At the restart rate (~0.4 s wall/step) that is
**~1.2 GB/s**. The `plot_run.sh` watcher (`--delete --keep-last 2`) cannot GC that
fast, so plotfiles pile up — one short extended run left **112 GB** of
`RotatingWormholePlt*` backlog. The extracted results (`frames/`, `small_data/`)
are tiny (~20 MB).

- For any long run, raise `plot_interval` (e.g. 10–50) or reduce `amr.plot_vars`.
- The backlog is disposable: frames + Psi4 are already extracted. Safe to delete
  `output/RotatingWormholePlt*` at any time; keep `frames/`, `small_data/`, `data/`,
  and checkpoints.

### Fixed-gauge diagnostic (`params_rotating_fixedgauge_dx025.txt`)

To separate the gauge from the source as the blow-up driver, the unigrid weak-spin
run was repeated with moving-puncture driving **disabled** (`lapse_coeff = 0`,
`shift_Gamma_coeff = 0`, `lapse_advec_coeff = 0`, `shift_advec_coeff = 0`), so
lapse/shift stay at their stationary Teo values. Same data, same frozen source.

| t | Ham (moving puncture) | Ham (fixed gauge) |
|------|-----------------------|-------------------|
| 0.005 | 5.3e-3 | 5.3e-3 |
| 0.05  | 0.386  | 0.194 |
| 0.10  | 0.853  | 0.422 |
| 0.20  | **1.786** | **0.719** |
| 0.50  | (NaN earlier) | 2.21 |

**Conclusion:** the gauge accounts for ~half the growth, but the constraint still
diverges with a frozen gauge. The Teo data + frozen source is **not a discrete
equilibrium**, so no gauge choice rescues it. This rules out the
prescribed-frozen-source design as a route to a stable evolution, and confirms the
real fix must make the source consistent / co-evolving. Output cleaned (kept
`data/constraint_norms.dat`).

### Root-cause summary (after the 4 Jun diagnostics)

Two independent problems, established by the controls + fixed-gauge tests:

- **P1 — AMR cross-level kinks.** `ExternalGridInitialData` trilinearly interpolates
  the `dx=0.25` file onto finer levels → C0 kinks → NaN. Only active with AMR.
- **P2 — frozen prescribed source (the unigrid killer).** `EffectiveTeoMatter`
  reads `teo_rho/j/S` from the state and has an **empty `add_matter_rhs`**, so the
  source is frozen at t=0 while the geometry evolves. `G_ab ≠ 8πT_ab` grows →
  constraints diverge. Independent of spin (`weak_spin ≡ a0`) and of interpolation
  (it blows up on the AMR-free unigrid). The fixed-gauge test (above) confirms no
  gauge choice fixes it: the data + frozen source is not a discrete equilibrium.

**Therefore: analytic-ID and gauge tweaks address only P1 / part of the symptom.
The decisive fix is to make the source consistent / co-evolving (attacks P2).**

### Next steps (superseded order — do not chase the prescribed-source run)

1. **Do not pursue physics from the unigrid weak-spin run as-is** — signal is at the
   0.08% noise floor and the run is constraint-violating by `t = 0.2`.
2. ~~Smooth initial data is the real fix~~ — **disproved** by the fixed-gauge test.
   Analytic ID fixes P1 (AMR smoothness) but the unigrid run still diverges from P2.
   Pursue analytic ID only as a *layer on top of* a consistent source, for AMR.
3. The correct observable remains the **difference** `weak_spin − a0` (and
   `− no_matter`), but only once the evolution is stable — today it is noise.

### Recommended path to a stable run (future work — the real fix for P2)

Goal: a stable evolution where a spin signal can rise above the `a0` baseline.
Reuse the machinery that already works (the non-rotating phantom-scalar wormhole is
stable *because* its matter co-evolves), in this order:

1. **Constraint-satisfying rotating initial data via the GRTresna elliptic solver.**
   - Driver: `grteclyn-wrapper/scripts/.../make_rotating_wormhole_id.py` (elliptic
     route — distinct from the pure-Python `make_teo_wormhole_gridinit.py`, which
     only FD-samples the analytic metric and produces the frozen source).
   - Solve the Hamiltonian **and** momentum constraints for the rotational `K_ij`
     (the current `.gridinit` route never solves the momentum constraint — that is
     why even fixed-gauge drifts).
   - Prereq: GRTresna MPI libs were rebuilt with consistent Chombo MPI
     (`grteclyn-wrapper/scripts/rebuild_grtresna_mpi.sh`) after a SIGILL from
     mismatched `-march=native` objects. Verify it builds/runs first.
2. **Evolve with a co-evolving matter model**, not the frozen tensor.
   - Use the proven `ExoticScalarField<PhantomDecayPotential>` path (the
     `wormhole_matter_model = exotic_scalar` branch in `SupportedWormholeLevel.cpp`),
     which has a non-empty `add_matter_rhs` and self-consistently supports the throat.
   - Introduce spin as genuine initial angular momentum (odd-parity `K_ij` / shift)
     from the elliptic solve, rather than `EffectiveTeoMatter`'s frozen `teo_*`.
   - Retire `EffectiveTeoMatter` for production, or keep it only for the static-ID
     regression test.
3. **Analytic C++ Teo `initData` (per-level, no interpolation)** — only now worth
   building, as the AMR-smoothness layer on top of constraint-satisfying co-evolving
   data. Mirror `SupportedWormholeInitialData` (closed-form `chi, h_ij, lapse,
   shift`); compute `K_ij, A_ij, Gamma^i` by FD over analytically-filled ghost cells
   so every level is smooth at its own `dx`.
4. **Fix the `Gamma^i` / constraint-order mismatch.** The Python ID reports
   `ham_l2 ≈ 1e-6` (interior, masked) but C++ level-0 sees `~5e-3` at t=0.005 and
   growing; the FD-computed `Gamma^i` is not consistent with GRTeclyn's BSSN
   definition / stencil. Recompute `Gamma^i` from `h_ij` with the evolution's own
   derivative order, or enforce it at init.

**Interim option (stepping stone, not physics):** make `EffectiveTeoMatter`
recompute `G_ab/8π` from the *current* evolving metric each step instead of reading
frozen `teo_*`. This keeps constraints satisfied by construction and bounds the
unigrid run, but the resulting "dynamics" are gauge/numerical, not matter response —
use only to validate the evolution pipeline, not for a spin claim.

### Retained / removed outputs (cleanup 4 June 2026)

Kept: `weak_spin_unigrid_dx025`, `weak_spin_a0_unigrid_dx025`,
`weak_spin_no_matter_unigrid_dx025` (frames + checkpoints + diagnostics), plus
`frames/` of the superseded runs. Removed: all `RotatingWormholePlt*` backlogs,
`weak_spin_ml1_dx025`, bulk of `weak_spin_effective` / `weak_spin_ml0_control`
(frames preserved), and empty stub dirs. Freed ~127 GB.

---

## Original failure (`weak_spin_effective`)

| Item | Value |
|------|--------|
| Params | `params_rotating.txt` |
| ID file | `runs/teo_wormhole/teo_weak_spin.gridinit` (`dx = 0.5`) |
| Grid | `N1=128`, `L=64` → level-0 `dx = 0.5` |
| AMR | `max_level = 5` → finest `dx ≈ 0.0156` (32× below file) |

Symptoms at t ≈ 0.13:

- NaN on **level 5**, component `h11` (then `phi`).
- Central checkerboard in `K` / `Weyl4` frames on the AMR patch only.
- `constraint_norms.dat`: monotonic L2 growth → numerical instability.

`SupportedWormholeCollapse` at the same `max_level` works because ID is **analytic
per level**, not interpolated from a coarse file.

---

## What we fixed (June 2026)

### 1. Finer Teo initial data (`dx = 0.25`)

Generated from repo root:

```bash
uv run python grteclyn-wrapper/scripts/wormhole/make_teo_wormhole_gridinit.py \
  --output runs/teo_wormhole/teo_weak_spin_dx025.gridinit \
  --nx 256 --ny 256 --nz 128 --lx 64 --ly 64 --lz 32 \
  --center 32 32 0 --spin 0.05 --check
```

Artifacts:

| File | Size | Notes |
|------|------|--------|
| `runs/teo_wormhole/teo_weak_spin_dx025.gridinit` | ~2.4 GB | `dx = 0.25` |
| `runs/teo_wormhole/teo_weak_spin_dx025.manifest.json` | | `ham_l2 ≈ 1.6e-6`, `mom_l2 ≈ 1.7e-6` |

Spin-0 baseline for junk subtraction:

```bash
uv run python grteclyn-wrapper/scripts/wormhole/make_teo_wormhole_gridinit.py \
  --output runs/teo_wormhole/teo_a0_dx025.gridinit \
  --nx 256 --ny 256 --nz 128 --lx 64 --ly 64 --lz 32 \
  --center 32 32 0 --spin 0.0 --check
```

The CLI prints progress lines (`Building…`, `Writing…`, `Evaluating…`). This path
is **pure Python** (not GRTresna).

### 2. New params files (recommended runs)

| Params file | Purpose | `max_level` | Grid / file match |
|-------------|---------|-------------|-------------------|
| `params_rotating_unigrid_dx025.txt` | **Production** weak spin | 0 | `N=256³/2`, `dx=0.25` = file |
| `params_rotating_a0_unigrid_dx025.txt` | Spin-0 control | 0 | same |
| `params_rotating_no_matter_unigrid_dx025.txt` | Vacuum control (`no_matter`) | 0 | same ID, source off |
| `params_rotating_ml0_control.txt` | Sanity: old `dx=0.5` file, no AMR | 0 | `N=128`, unigrid |
| `params_rotating_ml1_dx025.txt` | **Failed** AMR test | 1 | `N=128`, finest `dx=0.25` |
| `params_rotating.txt` | Legacy (unsafe) | 5 | original crash config |

### 3. Validation runs (stop_time = 0.2)

| Run directory | Result | Notes |
|---------------|--------|--------|
| `runs/teo_wormhole/weak_spin_ml0_control/` | OK | Confirms crash was AMR/sub-file, not bad `dx=0.5` ID |
| `runs/teo_wormhole/weak_spin_unigrid_dx025/` | OK to t=0.2; **NaN ≈ t=0.9** | Reaches t=0.2; diverges when extended (see 4 Jun update) |
| `runs/teo_wormhole/weak_spin_a0_unigrid_dx025/` | OK | Non-spin baseline — **matches weak_spin to ~0.08%** |
| `runs/teo_wormhole/weak_spin_no_matter_unigrid_dx025/` | OK | Ham L2 ≈ 0.009 at t=0.2 (vacuum) |
| `runs/teo_wormhole/weak_spin_ml1_dx025/` | **NaN** t≈0.17, level 1 | Passed t=0.13; AMR still unstable |

At t = 0.2, weak-spin unigrid constraint L2 Ham ≈ 1.79 (evolution growth, not a
static ID artifact). No-matter control stays O(10⁻³).

---

## How to run

### Build

```bash
cd Examples/RotatingWormholeCollapse
make -j 8 USE_CUDA=TRUE USE_MPI=TRUE COMP=gnu CUDA_ARCH=90
```

### Production evolution (recommended)

Use two terminals from repo root.

**Terminal 1 — live frames, small-data extraction, plotfile deletion:**

```bash
./grteclyn-wrapper/scripts/plot/plot_run.sh \
  /home/jovyan/nachevsky/test/simulation/GRTeclyn/runs/teo_wormhole/weak_spin_unigrid_dx025/output
```

`plot_run.sh` watches `RotatingWormholePlt*`, writes frames under `output/frames/`,
extracts Psi4/areal-radius to `output/small_data/`, and deletes processed HDF5
plot dirs (`--delete --keep-last 2`).

**Terminal 2 — simulation:**

```bash
cd Examples/RotatingWormholeCollapse
CUDA_VISIBLE_DEVICES=0,1 mpirun -n 2 ./main3d.gnu.MPI.CUDA.ex \
  params_rotating_unigrid_dx025.txt \
  2>&1 | tee ../../runs/teo_wormhole/weak_spin_unigrid_dx025/run.log
```

For a longer run, increase `stop_time` in the params file (or override on the
command line). Finer timestep: with `N=256`, coarse `dt ≈ 0.005` (40 steps to
t = 0.2).

**If the watcher was not started:** batch-drain plotfiles after the run:

```bash
cd grteclyn-wrapper
uv run python -m grteclyn_wrapper.visualisation.process_wave.consume_plotfiles \
  --data ../runs/teo_wormhole/weak_spin_unigrid_dx025/output \
  --out ../runs/teo_wormhole/weak_spin_unigrid_dx025/output/small_data \
  --radii 12 16 20 24 --areal-radius --embedding --embedding-rmax 5.0 \
  --frames-fields K Weyl4_Re --frames-axis z --frames-corner \
  --frames-out ../runs/teo_wormhole/weak_spin_unigrid_dx025/output/frames \
  --delete --keep-last 2 -j 8
```

### Comparison controls

```bash
# Spin-0 baseline (same grid, teo_a0_dx025.gridinit)
mpirun -n 2 ./main3d.gnu.MPI.CUDA.ex params_rotating_a0_unigrid_dx025.txt

# Vacuum relaxation (wormhole_matter_model=no_matter)
mpirun -n 2 ./main3d.gnu.MPI.CUDA.ex params_rotating_no_matter_unigrid_dx025.txt
```

Point `plot_run.sh` at each run’s `output/` directory the same way.

### Restart from checkpoint

After extending `stop_time`, resume from e.g. step 40:

```bash
cd /path/to/GRTeclyn
./grteclyn-wrapper/scripts/wormhole/rollback --step 40 \
  --data runs/teo_wormhole/weak_spin_unigrid_dx025/output

./grteclyn-wrapper/scripts/plot/plot_run.sh --not-remove \
  runs/teo_wormhole/weak_spin_unigrid_dx025/output

cd Examples/RotatingWormholeCollapse
CUDA_VISIBLE_DEVICES=0,1 mpirun -n 2 ./main3d.gnu.MPI.CUDA.ex \
  params_rotating_unigrid_dx025.txt \
  amr.restart=../../runs/teo_wormhole/weak_spin_unigrid_dx025/output/RotatingWormholeChk00040
```

### Legacy / unsafe config

`params_rotating.txt` + `teo_weak_spin.gridinit` + `max_level=5` reproduces the
original failure. Keep for regression only.

---

## Resolution and AMR constraints

With file spacing `dx_file` and level-0 spacing `dx_0 = L/N1`:

- **Safe:** `dx_file · 2^max_level ≤ dx_0` (finest level not finer than file).
- With `dx_file = 0.25` and `N1 = 128` (`dx_0 = 0.5`): `max_level = 1` has
  finest `dx = 0.25` (at file resolution in theory), but **still NaN'd** in our
  test — likely AMR prolongation/restriction + evolution, not ID check alone.
- With `dx_file = 0.25` and `N1 = 256` (`dx_0 = 0.25`): `max_level = 0` matches
  file exactly (**stable**).

File-size guide (`L=64`, `Lz=32`, 37 components):

| `dx_file` | (nx, ny, nz) | ~size |
|-----------|----------------|-------|
| 0.5 | 128, 128, 64 | 0.3 GB |
| 0.25 | 256, 256, 128 | 2.4 GB |
| 0.125 | 512, 512, 256 | ~20 GB |

---

## Root-cause reference (unchanged mechanism)

1. **ID is good at native resolution** — Python `--check` and direct
   `constraint_residuals()` give ~1e-6 L2 away from the throat.
2. **Noise is localized** to the AMR-refined throat patch in frames.
3. **Loader:** `ExternalGridInitialData.hpp` trilinear interpolation at each
   level’s `dx` in `SupportedWormholeLevel::initData()`.
4. **Tagging:** `ChiTagger` refines the throat to `max_level`.

C++ level-0 constraint norms (~0.1) can exceed Python (~1e-7) due to level-0-only
evaluation, boundary layers, and `Gamma^i` / derivative order mismatch — secondary
to the refined-level instability.

---

## Future work (not done)

> See **"Recommended path to a stable run"** above for the ordered plan. The table
> below classifies each candidate by which problem it actually addresses. Note:
> the items that only address **P1** will **not** stabilise the unigrid run on
> their own — **P2 (consistent/co-evolving source) is the prerequisite.**

| Approach | Addresses | Benefit / caveat |
|----------|-----------|------------------|
| **GRTresna elliptic rotating ID** (`make_rotating_wormhole_id.py`) + co-evolving `ExoticScalarField` matter | **P2** | The real fix: constraint-satisfying data + a source that evolves. Highest priority. |
| Consistent source: recompute `G_ab/8π` from the evolving metric in `EffectiveTeoMatter` each step | **P2** (partial) | Bounds the unigrid run, but "dynamics" are gauge/numerical — stepping stone only, not a spin claim. |
| C++ analytic Teo `initData` (like `SupportedWormholeInitialData`) | P1 | Smooth on every AMR level. Only useful *after* P2 is solved. |
| Fix FD `Gamma^i` to match GRTeclyn's BSSN stencil at init | P1/P2 | Removes the `~5e-3 → ` growing level-0 residual seen with the FD gridinit. |
| Tricubic (or higher) interpolation in `ExternalGridInitialData` | P1 | Reduces kinks without huge files. |
| `dx = 0.125` file + `max_level = 2` with `N1 = 128` | P1 | More AMR depth — only if evolution stabilises (i.e. after P2). |
| Tune `sigma`, `dt_multiplier`, `regrid_threshold` | symptom | Fixed-gauge test showed gauge ~halves growth but does not stop it; tuning alone is insufficient. |

**GRTresna:** elliptic ID is a separate route (`make_rotating_wormhole_id.py`); the
pure-Python `make_teo_wormhole_gridinit.py` does **not** use it and only FD-samples
the analytic metric (producing the frozen source). GRTresna MPI on this machine was
rebuilt with consistent Chombo MPI libs
(`grteclyn-wrapper/scripts/rebuild_grtresna_mpi.sh`) after SIGILL from mismatched
`-march=native` objects — verify it builds/runs before relying on it.

---

## Key paths

**Python / ID**

- `grteclyn-wrapper/scripts/wormhole/make_teo_wormhole_gridinit.py`
- `grteclyn-wrapper/src/grteclyn_wrapper/initial_data/teo.py`

**C++**

- `Examples/RotatingWormholeCollapse/ExternalGridInitialData.hpp` — trilinear load
- `Examples/RotatingWormholeCollapse/SupportedWormholeLevel.cpp` — evolution driver
- `Source/Matter/EffectiveTeoMatter.hpp` — prescribed `teo_*` source

**Stable run outputs**

- `runs/teo_wormhole/weak_spin_unigrid_dx025/output/` — Chk00040, `data/`, `frames/`
- `runs/teo_wormhole/weak_spin_unigrid_dx025/run.log`

**Reproduce Python constraint check**

```bash
cd grteclyn-wrapper
uv run python - <<'PY'
from grteclyn_wrapper.initial_data.teo import TeoWormholeConfig, constraint_residuals
cfg = TeoWormholeConfig(nx=256, ny=256, nz=128, Lx=64, Ly=64, Lz=32,
                        center=(32,32,0), spin=0.05, source="einstein")
r = constraint_residuals(cfg)
print(r.ham_l2, r.ham_max, r.mom_l2, r.mom_max, r.valid_fraction)
PY
```
