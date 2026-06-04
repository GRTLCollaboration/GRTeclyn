# RotatingWormholeCollapse

This example evolves a **rotating traversable wormhole** with the CCZ4 + matter
machinery shared with `SupportedWormholeCollapse`. The default, recommended path
uses **analytic Teo initial data** (Teo 1998, gr-qc/9803098) loaded from a
`.gridinit` file, which starts the run close to a constraint-satisfying state and
minimises the junk radiation seen with constraint-violating seeds.

## How it works

The pipeline has three stages: generate initial data, load it, evolve it.

### 1. Analytic initial data (Python)

The Teo geometry is a stationary, axisymmetric wormhole. Rather than solving an
elliptic problem, we sample the *known* metric on a uniform Cartesian grid and
do the 3+1 (ADM) decomposition in Python:

- `grteclyn_wrapper.initial_data.teo` (in `grteclyn-wrapper/src/`) is the physics
  module. Given a `TeoWormholeConfig` (grid size, throat radius `b0`, spin
  `a_hat`, source model) it computes:
  - the ADM variables — lapse `alpha`, shift `beta^i` (the frame-dragging
    `omega`), spatial metric `gamma_ij`;
  - the CCZ4 state GRTeclyn evolves — `chi`, conformal metric `h_ij`, `K`,
    `A_ij`, conformal connection `Gamma^i`;
  - the **effective stress-energy** `(rho, j_i, S_ij) = G_ab / 8pi`, computed
    numerically from the 4-metric. Teo specifies a *geometry* and infers its
    matter content, so the source is stored as data, not derived from a
    fundamental field.
- `grteclyn-wrapper/scripts/wormhole/make_teo_wormhole_gridinit.py` is a thin CLI that
  builds a config, calls the module, and writes the `.gridinit` plus a
  `.manifest.json` of diagnostics.

The throat is regularised at the puncture; the generator asserts finite values
and positive `chi`/`lapse` before writing.

### 2. Loading into GRTeclyn (C++)

`params_rotating.txt` sets `recipe_initial_data_file` to the generated
`.gridinit`. In `SupportedWormholeLevel::initData`, when that key is non-empty,
`ExternalGridInitialData` reads the file and fills every state component
(including the `teo_*` source fields) by interpolation onto the simulation grid.
If the key is empty, it falls back to the analytic `SupportedWormholeInitialData`
(the phantom-scalar wormhole).

### 3. Matter model and evolution (C++)

The `wormhole_matter_model` parameter selects how the matter source enters the
CCZ4 right-hand side and the constraints. It is wired in `variableSetUp`,
`specificEvalRHS`, and `specificPostTimeStep`:

| `wormhole_matter_model` | Class | Source of `T_ab` |
| --- | --- | --- |
| `effective_teo` (Teo runs) | `EffectiveTeoMatter` | reads the prescribed `teo_rho/j/S` fields, scaled by `teo_source_strength` |
| `no_matter` | `NoMatter` | identically zero (vacuum relaxation control) |
| `exotic_scalar` (default) | `ExoticScalarField<PhantomDecayPotential>` | evolves a phantom scalar `phi`, `Pi` |

`EffectiveTeoMatter` is a *prescribed* source: it returns the stored `teo_*`
values and has an empty `add_matter_rhs`, because the Teo source is fixed initial
data, not a dynamically evolved field. `teo_source_strength` lets you scale it
(e.g. a support-threshold scan).

Each step, `specificPostTimeStep` recomputes the Hamiltonian and momentum
constraint L2 norms with the matching matter model and writes
`data/constraint_norms.dat`, plus `data/collapse_diagnostics.dat` (min lapse/chi,
`max|K|`, an apparent-horizon `theta_+` proxy, scalar field range). Weyl4 is
extracted on the configured radii for gravitational-wave output.

### Symmetry / boundary conditions

Genuine z-axis rotation is **incompatible with octant symmetry**: the shift
`beta^x = omega y`, `beta^y = -omega x` has the wrong parity across the `x=0` and
`y=0` planes. Spinning runs therefore use full `x,y` domains with Sommerfeld
(non-reflective) outer boundaries and only an **equatorial z-reflection**:

```text
lo_boundary = 1 1 2     # Sommerfeld in x,y; reflective in z
hi_boundary = 1 1 1
```

The generated `.gridinit` `center`/`origin` must match the simulation `center`
so the throat lands on the same coordinate.

## Run variants (no extra params files)

There is a single `params_rotating.txt`. The scenario is chosen when generating
the `.gridinit` (`--spin`, `--source`); the matter behaviour is an AMReX
ParmParse command-line override. Examples:

```bash
# Default effective-source Teo run
./main3d.gnu.MPI.CUDA.ex params_rotating.txt

# Matched a=0 baseline (regenerate gridinit with --spin 0.0) for junk subtraction
./main3d.gnu.MPI.CUDA.ex params_rotating.txt \
  recipe_initial_data_file=../../runs/teo_wormhole/teo_a0.gridinit

# Vacuum-relaxation control: same data, source switched off
./main3d.gnu.MPI.CUDA.ex params_rotating.txt wormhole_matter_model=no_matter

# Support-threshold scan
./main3d.gnu.MPI.CUDA.ex params_rotating.txt teo_source_strength=0.8
```

## Generate Teo initial data

From the repo root:

```bash
# Weak-spin effective-source data (the params_rotating.txt default)
uv run python grteclyn-wrapper/scripts/wormhole/make_teo_wormhole_gridinit.py \
  --output runs/teo_wormhole/teo_weak_spin.gridinit \
  --nx 128 --ny 128 --nz 64 --lx 64 --ly 64 --lz 32 \
  --center 32 32 0 --spin 0.05

# Matched non-spinning baseline
uv run python grteclyn-wrapper/scripts/wormhole/make_teo_wormhole_gridinit.py \
  --output runs/teo_wormhole/teo_a0.gridinit \
  --nx 128 --ny 128 --nz 64 --lx 64 --ly 64 --lz 32 \
  --center 32 32 0 --spin 0.0

# Same geometry but zeroed source (pair with wormhole_matter_model=no_matter)
uv run python grteclyn-wrapper/scripts/wormhole/make_teo_wormhole_gridinit.py \
  --output runs/teo_wormhole/teo_weak_spin_zero_source.gridinit \
  --nx 128 --ny 128 --nz 64 --lx 64 --ly 64 --lz 32 \
  --center 32 32 0 --spin 0.05 --source zero
```

The ergoregion threshold is near `|a_hat| = 0.5`; keep `--spin` below that for
the first controlled runs. Each call prints and stores summary diagnostics
(min/max `chi`, `lapse`, `max|shift|`, `teo_rho` range) in the manifest.

### Validate the initial data (self-consistency)

Pass `--check` to also evaluate the full 3D ADM constraint residuals from the
generated fields before trusting the data:

```bash
uv run python grteclyn-wrapper/scripts/wormhole/make_teo_wormhole_gridinit.py \
  --output runs/teo_wormhole/teo_weak_spin.gridinit \
  --nx 128 --ny 128 --nz 64 --lx 64 --ly 64 --lz 32 \
  --center 32 32 0 --spin 0.05 --check
```

This computes

```text
Hamiltonian:  H   = R + K^2 - K_ij K^ij - 16 pi rho
Momentum:     M_i = D_j (K^j_i - delta^j_i K) - 8 pi j_i
```

excluding the boundary finite-difference layer and a small ball around the
throat, then reports L2 (RMS) and max norms (added to the `.manifest.json` under
`constraint_residuals`). Because the effective source is `T_ab = G_ab / 8pi` from
the same metric, the constraints hold *by construction* up to truncation error,
so this is a self-consistency check (it confirms the `K_ij` and Einstein-tensor
derivative paths agree away from the throat). Small residuals are expected; the
throat region is a puncture and is intentionally not resolved. A genuine
resolution-convergence test (fixed physical mask radius across resolutions) is a
separate, stronger check.

> Note: `make_rotating_wormhole_id.py` (GRTresna elliptic solve) remains in the
> tree as a secondary/legacy initial-data route for cross-validation. The Teo
> analytic path above is the primary one.

## Build

```bash
# Production multi-GPU build
make -j 8 USE_CUDA=TRUE USE_MPI=TRUE COMP=gnu CUDA_ARCH=90
```

## Run

**Multi-GPU:**

```bash
# 2 GPUs
CUDA_VISIBLE_DEVICES=0,1 mpirun -n 2 ./main3d.gnu.MPI.CUDA.ex params_rotating.txt

# 4 GPUs
CUDA_VISIBLE_DEVICES=0,1,2,3 mpirun -n 4 ./main3d.gnu.MPI.CUDA.ex params_rotating.txt

# 8 GPUs
CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7 mpirun -n 8 ./main3d.gnu.MPI.CUDA.ex params_rotating.txt
```

**Recommended GPU binding (avoids mapping warnings):**
```bash
mpirun -n 8 bash -c 'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; exec ./main3d.gnu.MPI.CUDA.ex params_rotating.txt'
```

## Restart from Checkpoint

To resume from any checkpoint (for example step **3000** or **4000**):

1. **Rollback** the data directory to the desired checkpoint (use paths matching your clone of this repo and the `output_path` set in `params_rotating.txt`):
   ```bash
   cd /path/to/GRTeclyn
   ./grteclyn-wrapper/scripts/wormhole/rollback --step 3000 --data /path/to/your/simulation/output
   ```
   (Replace `3000` with the desired checkpoint number, e.g. `4000`.)

2. **Start the watcher** (in a second terminal — use `--not-remove`):
   ```bash
   cd /path/to/GRTeclyn
   ./grteclyn-wrapper/scripts/plot/plot_run.sh --not-remove /path/to/your/simulation/output
   ```

3. **Restart the simulation** using `amr.restart` (this is the correct flag):
   ```bash
   cd /path/to/GRTeclyn/Examples/RotatingWormholeCollapse

   # 2 GPUs example
   CUDA_VISIBLE_DEVICES=0,1 mpirun -n 2 ./main3d.gnu.MPI.CUDA.ex params_rotating.txt \
     amr.restart=/path/to/your/simulation/output/RotatingWormholeChk03000
   ```

> **Tip**: Always replace both the `--step N` in rollback **and** the `ChkN` number in the restart command with the same value (e.g. 4000).

## Plotting (normal run)

```bash
# In a second terminal (from the GRTeclyn root; set RUN_DIR to your simulation output)
./grteclyn-wrapper/scripts/plot/plot_run.sh /path/to/your/simulation/output
```

## Quick MPI Setup (if needed)

Point these at **your** OpenMPI install prefix (the directory that contains `bin/mpirun` and `lib`):

```bash
export PATH=/path/to/your/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/path/to/your/openmpi/lib:${LD_LIBRARY_PATH:-}
```
