# RotatingWormholeCollapse — Debug notes

## Update (8 June 2026): Route B — spinning complex phantom-scalar wormhole

This is the new modelling route for the **rotating** wormhole collapse paper. It
replaces the failed single-real-scalar GRTresna run (which produced a four-lobed
pattern, dispersing matter, and a frozen embedding) with a self-consistent,
co-evolving **complex** phantom scalar that rotates smoothly.

### Why the previous rotating run failed

- The GRTresna `exotic_scalar` config solved a weak Gaussian lump (`amp=0.1`,
  `width=12`, `bh_mass=0`) with no throat-forming term, so `chi ≈ 1` everywhere
  → the embedding never changed.
- A *real* scalar carries angular momentum only through an azimuthal `cos(mφ)`
  density modulation. That modulation **is** the four-lobed pattern, and being
  non-stationary it radiates/disperses outward ("matter flying away"). The 4-leaf
  and the dispersal share this one root cause.

### What was implemented

A complex phantom scalar `Φ = f(r,θ) e^{i(m φ_az − ω t)}` stored as two real
components `(φ₁,Π₁)` and `(φ₂,Π₂)`. Its modulus `|Φ|² = φ₁²+φ₂²` is axisymmetric
(no lobes) and it carries `J_z` through the phase winding (stationary → no
dispersal). Collapse is triggered exactly as in the Ellis–Bronnikov paper by
reducing the phantom support `S_support < 1`; the inherent rotational quadrupole
seeds `ℓ=2` GW with no artificial `A_φ`.

Files added/changed:

| File | Change |
|------|--------|
| `Source/Matter/ComplexScalarFieldVars.hpp` | New: `φ₁/Π₁/φ₂/Π₂` state accessors |
| `Source/Matter/ComplexScalarFieldD1Vars.hpp` | New: first derivatives of `φ₁,φ₂` |
| `Source/Matter/ComplexScalarFieldD2Vars.hpp` | New: second derivatives of `φ₁,φ₂` |
| `Source/Matter/ComplexScalarFieldAdvecVars.hpp` | New: advection terms for `φ₁,Π₁,φ₂,Π₂` |
| `Source/Matter/ComplexExoticScalarField.hpp/.impl.hpp` | New: phantom complex matter model (summed EM tensor with `−S_support` flip; two KG RHSs) |
| `Source/Matter/Make.package` | Registers the new headers |
| `StateVariables.hpp` | Adds `c_phi2`, `c_Pi2` (even parity) |
| `SupportedWormholeInitialData.hpp` | `complex_scalar` toroidal spinning ansatz |
| `SimulationParameters.hpp` | `complex_scalar` model + `wormhole_azimuthal_m`, `wormhole_rotation_omega` |
| `SupportedWormholeLevel.cpp` | `complex_scalar` branches in `variableSetUp` / `specificEvalRHS` / constraints; new `J_z` diagnostic column |
| `PhantomDecayPotential.hpp` | `compute_potential_value(V,dV,φ)` value overload (per-component, exact for the separable quadratic potential) |
| `params_rotating_complex_equilibrium.txt` | `S_support=1` equilibrium test |
| `params_rotating_complex_collapse.txt` | `S_support=0.5` collapse run |

The initial data is **analytic per-level** (no `.gridinit`), mirroring what made
the EB collapse stable: per-cell `chi` throat, flat lapse `α=1`, toroidal
`f = φ_EB(r) (sinθ)^m`, `φ₁=f cos(mφ_az)`, `φ₂=f sin(mφ_az)`, `Π₁=ω φ₂`,
`Π₂=−ω φ₁`, with `K_ij=0` (accept the small `O(ω)` momentum defect at `t=0`, as
the EB paper accepts an `O(A_φ)` Hamiltonian defect). The optional fully
constraint-satisfying route (elliptic solve of the winding scalar via GRTresna)
is future work.

### Build (NOTE: not buildable on the current machine)

```bash
cd Examples/RotatingWormholeCollapse
make -j 8 USE_CUDA=TRUE USE_MPI=TRUE COMP=gnu CUDA_ARCH=90
```

A CPU build (no GPU) for a quick compile check:

```bash
make -j 8 USE_CUDA=FALSE USE_MPI=FALSE COMP=gnu DIM=3
```

The code in this update was written without an available build/CI on this
machine; it was reviewed for interface consistency against the existing
`ExoticScalarField` / `ScalarField*Vars` matter machinery (same `emtensor_t`,
`CCZ4Vars`/`CCZ4D1Vars`/`CCZ4D2Vars` bases, and `FourthOrderDerivatives`
`diff1_scalar`/`diff2_scalar`/`advec_scalar` calls). Compile before any
production run.

### How to test

**1. Equilibrium (validates the rotating data is self-consistent).**

```bash
# Terminal 1 — live frames (start before evolution)
GRTECLYN_FRAMES_ZOOM="48" GRTECLYN_FRAMES_CENTER="8 8 0" \
GRTECLYN_EXTRACTION_CENTER="32 32 0" \
./grteclyn-wrapper/scripts/plot/plot_run.sh \
  runs/rotating_wormhole/complex_equilibrium/output

# Terminal 2 — evolution
source grteclyn-wrapper/scripts/lib/env.sh
cd Examples/RotatingWormholeCollapse
mpirun -n 2 bash -c \
  'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; \
   exec ./main3d.gnu.MPI.CUDA.ex params_rotating_complex_equilibrium.txt'
```

Pass criteria:
- `phi`, `Pi`, `phi2`, `Pi2` frames each show a single rotating lobe pair, but
  the **modulus** `phi^2+phi2^2` (and `chi`, the embedding) is **axisymmetric** —
  no four-leaf clover.
- `data/collapse_diagnostics.dat`: `min_chi`, `max_ah_r`, and the throat stay
  ~stationary for several `M`; `J_z` is **nonzero and roughly constant** (the
  spin is held), unlike the old run where matter dispersed.
- `data/constraint_norms.dat`: `L2_Ham`, `L2_Mom` stay bounded (small `O(ω)`
  initial momentum defect that does not blow up).

**2. Collapse (the physics run for the paper).**

```bash
cd Examples/RotatingWormholeCollapse
mpirun -n 2 bash -c \
  'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; \
   exec ./main3d.gnu.MPI.CUDA.ex params_rotating_complex_collapse.txt'
```

Pass criteria:
- `min_chi` and `min_lapse` drop, `max_ah_r > 0` and `min_theta_plus ≤ 0`
  (trapped surface forms) — throat collapse.
- `Weyl4` extraction (`data/`) shows an `ℓ=2` signal; check propagation speed
  `v ≈ c` across `extraction_radii = 12 16 20 24` (same analysis as the EB
  paper). A dynamically-grown `m=2` bar mode may appear for larger `omega`.
- `J_z` is carried through the collapse rather than radiated to zero immediately.

**3. Parameter sweeps.** Vary `wormhole_rotation_omega` (e.g. 0.0, 0.05, 0.1,
0.2) and `wormhole_azimuthal_m` (1, 2). `omega=0` is the non-rotating control
(should reduce to a real-scalar-like wormhole with `J_z≈0`).

**4. Loader-fix regression.** The `.gridinit` loader now maps columns by name.
Verify with any existing `effective_teo` `.gridinit`: on load it prints how many
columns were skipped, and `teo_rho`/`teo_S*` now populate the correct
`c_teo_*` slots (previously `teo_S11` landed in `c_teo_rho`). The GRTresna
`exotic_scalar` path is unchanged (its `chi…Pi` columns are name-matched too).

### Known limitations / next steps

- Initial data carries an `O(ω)` momentum-constraint defect (`K_ij=0`). For
  publishable, fully constraint-satisfying rotating data, solve the momentum
  (and Hamiltonian) constraints elliptically for the winding scalar (extend the
  GRTresna lump ansatz to a complex/phase-winding field) and load via the
  name-mapped `.gridinit`.
- `max_level=0` (unigrid) in the provided params to avoid the documented AMR
  cross-level interpolation kinks; enable AMR only after confirming stability.

---

## Update (8 June 2026): loader alignment bug + real-scalar rotation limit

- **Gridinit loader misalignment (fixed).** The `dust_*` state variables were
  inserted into `StateVariables.hpp` *after* the Teo `.gridinit` writer
  (`GRTECLYN_STATE_VARS` in `grtresna/io.py`) was built. The C++ loader copied
  file column `c` to state slot `c` by index, so the Teo effective source loaded
  into the wrong slots (`teo_S11 → c_teo_rho`, …) and `EffectiveTeoMatter`
  sourced a scrambled stress tensor. `ExternalGridInitialData` now maps columns
  to state slots **by name** via the `component_names` header (identity fallback
  for legacy V1 files). The GRTresna `chi…Pi` block is index-0-aligned, so the
  `exotic_scalar` path was unaffected.
- **A single real scalar cannot rotate smoothly.** A real field carries angular
  momentum only through an azimuthal `cos(mφ)` density modulation — that *is* the
  four-lobed pattern, and being non-stationary it radiates/disperses outward
  ("matter flying away"). The fix (Route B) uses a **complex** phantom scalar
  whose `|Φ|²` is axisymmetric and carries `J` through phase winding.


Investigation and fix for the crashing Teo weak-spin run (`weak_spin_effective`).
**Status: Teo failure root cause confirmed; GRTresna + co-evolving scalar
production path now validated to `t = 4.0`.**

## TL;DR

- The working production solution is now **GRTresna solved ID + co-evolving
  `exotic_scalar` matter**, validated to `t = 4.0` with no NaN.
- The original crash is **not** bad Teo physics in the `.gridinit`. The file is
  constraint-satisfying at its native resolution.
- It fails when AMR refines **below** the file resolution: `ExternalGridInitialData`
  fills every level by **trilinear** interpolation, so fine levels get C0 kinks and
  CCZ4 blows up at the throat (checkerboard → NaN in `h11`).
- **Teo diagnostic fix in practice:** use a `dx = 0.25` gridinit and evolve on a
  grid that does **not** go finer than the file. This addresses the AMR
  interpolation failure only; it does not make the frozen Teo source a production
  physics run.
- **`max_level=1` with `N1=128` still NaNs** at t ≈ 0.17 on level 1 (past the old
  t = 0.13 crash). Do not use that AMR config until evolution/interpolation is
  improved.

**Governing rule:** finest evolved `dx` must be **≥** file `dx` (no sub-file
interpolation on refined levels).

---

## Update (6 June 2026): GRTresna production path validated

The production path now uses GRTresna rotating exotic-scalar initial data loaded
through `.gridinit`, then evolves with `wormhole_matter_model = exotic_scalar`.
This attacks the frozen-source problem directly: scalar fields `phi` and `Pi`
co-evolve, instead of freezing a prescribed `teo_*` tensor.

The working solution is one rotating exotic scalar lump (`omega = 0.05`,
amplitude `0.2`, width `8`, mode `m = 2`). GRTresna solves the elliptic
constraints for this lump, converts the Chombo output to `.gridinit`, and
GRTeclyn evolves that solved geometry plus the scalar fields. This is the
production route; the analytic Teo + `effective_teo` route is kept only as a
diagnostic/regression case.

### Fixes made

- `grteclyn-wrapper/src/grteclyn_wrapper/grtresna/solver.py`
  - Resolves relative `work_dir` / `gridinit_path` to absolute paths before
    launching GRTresna. Without this, Chombo `ParmParse` was launched from the
    case directory with a now-invalid relative params path and aborted before
    the solver loop.
  - Forces `CONDA_PREFIX` to the resolved `grtresna` env when using that env's
    `mpirun`, matching the wrapper README's build/run assumptions.
- `grteclyn-wrapper/src/grteclyn_wrapper/grtresna/io.py`
  - Fixed the Chombo-to-`.gridinit` converter to fill every target cell covered
    by a Chombo source cell. The old one-cell paint left most of a finer
    `.gridinit` as zeros (`chi/lapse/h_ij = 0`), causing immediate `h11` NaNs.
- `Source/Matter/CCZ4RHSWithMatter.impl.hpp`
  - Zero-initializes each RHS cell before CCZ4, matter, and dissipation `+=`
    operations. This prevents passive fields such as `teo_*` from carrying
    uninitialized RHS values in `exotic_scalar` / `no_matter` runs.
- `grteclyn-wrapper/scripts/wormhole/rollback`
  - Recognizes `RotatingWormholeChk*` and `RotatingWormholePlt*` prefixes.
- `grteclyn-wrapper/scripts/plot/plot_run.sh`
  - Resolves the run directory to an absolute path before invoking the Python
    consumer. This avoids `uv run --directory grteclyn-wrapper` interpreting a
    relative `runs/...` path under the wrapper directory.
  - Uses `--n-points 128` and renders `chi`, `chi_minus_1`, `K`, `lapse`, `phi`,
    `Pi`, `scalar_activity`, `Weyl4_Re`, `Weyl4_Im`, and `Weyl4_Mag` frames by
    default while deleting processed plotfiles.
- `grteclyn-wrapper/src/grteclyn_wrapper/visualisation/process_wave/consume_plotfiles.py`
  - Tightened the default `K` frame scale from `[-5e-2, 5e-2]` to
    `[-5e-4, 5e-4]`. The old scale hid the evolved signal, whose magnitude is
    O(`1e-4`), making `K_z` appear blank even though the data was present.

### Generated ID

Command:

```bash
uv run python grteclyn-wrapper/scripts/wormhole/make_rotating_wormhole_id.py \
  --out-dir runs/rotating_wormhole_id \
  --omegas 0.0,0.05 \
  --ranks 2 \
  --nx 64 --ny 64 --nz 32 \
  --gridinit-nx 256 --gridinit-ny 256 --gridinit-nz 128 \
  --target-center-x 32 --target-center-y 32 --target-center-z 0 \
  --iterations 30 --timeout 1800
```

Results:

| Case | GRTresna convergence | Output |
|------|----------------------|--------|
| `omega=0.05` | Ham `0.7427%`, Mom `0.0115%` | `runs/rotating_wormhole_id/rotwh_omega_p0p05_amp_0p2_w_8/initial_data.gridinit` |
| `omega=0` | Ham `0.5761%`, Mom `NaN` because zero momentum gives an undefined relative norm | `runs/rotating_wormhole_id/rotwh_omega_p0_amp_0p2_w_8/initial_data.gridinit` |

### Short evolution validation

All three short unigrid checks run to `t = 0.22` with no NaN:

| Run | Params | Ham L2 at `t=0.22` | Mom L2 at `t=0.22` | Verdict |
|-----|--------|--------------------|--------------------|---------|
| weak spin | `params_rotating_grtresna_smoke.txt` | `3.49e-4` | `1.02e-5` | clean |
| omega=0 | `params_rotating_grtresna_a0_exotic.txt` | `2.88e-4` | `8.11e-6` | clean |
| no matter | `params_rotating_grtresna_no_matter.txt` | `3.62e-4` | `3.89e-5` | clean |

Run with rank-local GPU binding to avoid multi-rank GPU mapping stalls:

```bash
mpirun -n 2 bash -c 'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; exec ./main3d.gnu.MPI.CUDA.ex params_rotating_grtresna_smoke.txt'
```

This is the current recommended path. The old Teo `effective_teo` route remains
diagnostic/regression-only.

### Longer validation with live frames

Production weak-spin was rerun to `t = 4.0` using the GRTresna `.gridinit` and
`params_rotating_grtresna_exotic.txt` physics settings. The run used the CUDA MPI
binary with rank-local GPU binding:

```bash
mpirun -n 2 bash -c 'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; exec ./main3d.gnu.MPI.CUDA.ex ../../runs/rotating_wormhole/grtresna_spin_exotic/params_long_quality.txt'
```

Validation result:

| Run | End time | Ham L2 | Mom L2 | Verdict |
|-----|----------|--------|--------|---------|
| `grtresna_spin_exotic` | `t = 4.0` | `1.27e-4` | `7.44e-6` | clean |

The run reached `STEP = 200 TIME = 4`, wrote `RotatingWormholeChk00200` and
`RotatingWormholePlt00200`, and AMReX finalized normally. The live consumer
generated 21 frames per field under
`runs/rotating_wormhole/grtresna_spin_exotic/output/frames/` for:

`chi`, `chi_minus_1`, `K`, `lapse`, `phi`, `Pi`, `scalar_activity`, `Weyl4_Re`,
`Weyl4_Im`, `Weyl4_Mag`, and `embedding`.

Plotfile deletion worked as intended: only `RotatingWormholePlt00190` and
`RotatingWormholePlt00200` remained after processing (`--keep-last 2`).

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

### Recommended path to a stable run (completed on 6 June 2026)

Goal: a stable evolution where a spin signal can rise above the `a0` baseline.
Reuse the machinery that already works (the non-rotating phantom-scalar wormhole is
stable *because* its matter co-evolves). This plan is now implemented by the
GRTresna + `exotic_scalar` route described above:

1. **Constraint-satisfying rotating initial data via the GRTresna elliptic solver.**
   - Driver: `grteclyn-wrapper/scripts/.../make_rotating_wormhole_id.py` (elliptic
     route — distinct from the pure-Python `make_teo_wormhole_gridinit.py`, which
     only FD-samples the analytic metric and produces the frozen source).
   - Solve the Hamiltonian **and** momentum constraints for the rotating scalar
     source. The accepted `omega=0.05` ID converged to Ham `0.7427%`, Mom
     `0.0115%`.
   - Prereq: GRTresna MPI libs were rebuilt with consistent Chombo MPI
     (`grteclyn-wrapper/scripts/build/rebuild_grtresna_mpi.sh`) after a SIGILL from
     mismatched `-march=native` objects. Verify it builds/runs first.
2. **Evolve with a co-evolving matter model**, not the frozen tensor.
   - Use the proven `ExoticScalarField<PhantomDecayPotential>` path (the
     `wormhole_matter_model = exotic_scalar` branch in `SupportedWormholeLevel.cpp`),
     which has a non-empty `add_matter_rhs` and self-consistently supports the throat.
   - Introduce spin as genuine initial angular momentum (odd-parity `K_ij` / shift)
     from the elliptic solve, rather than `EffectiveTeoMatter`'s frozen `teo_*`.
   - Retire `EffectiveTeoMatter` for production, or keep it only for the static-ID
     regression test.
3. **Analytic C++ Teo `initData` (per-level, no interpolation)** — still optional
   future work, and only useful as an AMR-smoothness layer after the source is
   made consistent / co-evolving.
   Mirror `SupportedWormholeInitialData` (closed-form `chi, h_ij, lapse,
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

### 2. Historical Teo params files (diagnostic, superseded)

| Params file | Purpose | `max_level` | Grid / file match |
|-------------|---------|-------------|-------------------|
| `params_rotating_unigrid_dx025.txt` | Teo weak-spin diagnostic | 0 | `N=256³/2`, `dx=0.25` = file |
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
  runs/rotating_wormhole/grtresna_spin_exotic/output
```

`plot_run.sh` watches `RotatingWormholePlt*`, writes frames under `output/frames/`,
extracts Psi4/areal-radius to `output/small_data/`, and deletes processed HDF5
plot dirs (`--delete --keep-last 2`). It renders `chi`, `chi_minus_1`, `K`,
`lapse`, `phi`, `Pi`, `scalar_activity`, `Weyl4_Re`, `Weyl4_Im`, and
`Weyl4_Mag` frames by default.

**Terminal 2 — simulation:**

```bash
cd Examples/RotatingWormholeCollapse
mpirun -n 2 bash -c 'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; exec ./main3d.gnu.MPI.CUDA.ex params_rotating_grtresna_exotic.txt' \
  2>&1 | tee ../../runs/rotating_wormhole/grtresna_spin_exotic/run.log
```

For a longer run, increase `stop_time` in the params file (or override on the
command line). The validated long check used `stop_time = 4.0`, `plot_interval = 10`,
and `checkpoint_interval = 50`.

**If the watcher was not started:** batch-drain plotfiles after the run:

```bash
cd grteclyn-wrapper
uv run python -m grteclyn_wrapper.visualisation.process_wave.consume_plotfiles \
  --data ../runs/rotating_wormhole/grtresna_spin_exotic/output \
  --out ../runs/rotating_wormhole/grtresna_spin_exotic/output/small_data \
  --radii 12 16 20 24 --areal-radius --embedding --embedding-rmax 5.0 \
  --frames-fields chi chi_minus_1 K lapse phi Pi scalar_activity Weyl4_Re Weyl4_Im Weyl4_Mag \
  --frames-axis z --frames-corner \
  --frames-out ../runs/rotating_wormhole/grtresna_spin_exotic/output/frames \
  --delete --keep-last 2 -j 8
```

### Comparison controls

```bash
# Spin-0 baseline (same GRTresna scalar-lump setup with omega=0)
mpirun -n 2 bash -c 'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; exec ./main3d.gnu.MPI.CUDA.ex params_rotating_grtresna_a0_exotic.txt'

# Vacuum relaxation control
mpirun -n 2 bash -c 'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; exec ./main3d.gnu.MPI.CUDA.ex params_rotating_grtresna_no_matter.txt'
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
