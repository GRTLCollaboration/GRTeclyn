# grteclyn-wrapper

Run isolated GRTeclyn episodes from the repo root (`GRTeclyn/`). For RadialRecipe GPU smoke tests, the wrapper can stream plotfiles into small data during the run and delete heavy HDF5 dirs afterward.

## Prerequisites

From the GRTeclyn repo root:

```bash
cd /path/to/GRTeclyn
uv sync   # Python deps including yt for plotfile extraction
```

First build (single GPU, no MPI):

```bash
BUILD=1 bash grteclyn-wrapper/scripts/run_radialrecipe_gpu_smoke.sh
```

## Single GPU run (one guessed shape)

Pick **one** initial-data source:

| Env var | Example |
|---------|---------|
| `SEED_NAME` | `ellis_bronnikov` |
| `CANDIDATE_ID` | `bubble_wall_016`, `random_000` |
| `NONSPHERICAL_ID` | `quadrupole_bubble_001`, `dipole_lopsided_000` |

```bash
# Known seed
BUILD=0 SEED_NAME=ellis_bronnikov CUDA_VISIBLE_DEVICES_OVERRIDE=0 \
  bash grteclyn-wrapper/scripts/run_radialrecipe_gpu_smoke.sh

# Spherical guesser candidate
BUILD=0 CANDIDATE_ID=bubble_wall_016 CUDA_VISIBLE_DEVICES_OVERRIDE=1 \
  bash grteclyn-wrapper/scripts/run_radialrecipe_gpu_smoke.sh

# Non-spherical guessed shape
BUILD=0 NONSPHERICAL_ID=quadrupole_bubble_001 CUDA_VISIBLE_DEVICES_OVERRIDE=2 \
  bash grteclyn-wrapper/scripts/run_radialrecipe_gpu_smoke.sh
```

Outputs go to `runs/radialrecipe_gpu_smoke/<name>_gpu_t<stop_time>_<stamp>/`.

## Batch: 7 non-spherical shapes on GPUs 0–6

```bash
BUILD=0 bash grteclyn-wrapper/scripts/run_nonspherical_gpu_batch.sh
```

Outputs: `runs/radialrecipe_nonspherical/`. One log per candidate: `<id>_<stamp>.log`.

## Plotfile consumer (default on)

With `CONSUME_PLOTFILES=1` (default), each run:

1. Starts a sidecar `consume_plotfiles` process while the GPU simulation runs.
2. Appends `small_data/shell_profiles.dat`, `small_data/areal_radius.dat`, optional PNG frames.
3. Deletes processed plotfile dirs (`--keep-last 1`).
4. Runs a **post-sim drain** for any backlog.

Useful env vars:

| Variable | Default | Meaning |
|----------|---------|---------|
| `CONSUME_PLOTFILES` | `1` | Enable streaming extraction |
| `CONSUMER_DELETE` | `1` | Delete HDF5 plot dirs after extract |
| `CONSUMER_RADII` | `4 8` | Extraction radii |
| `PLOT_INTERVAL` | `10` if consumer on, else `1` | Plotfile write cadence |
| `STOP_TIME` | `2.0` | Simulation stop time |
| `N_FULL` | `64` | Grid resolution |
| `CUDA_VISIBLE_DEVICES_OVERRIDE` | `0` | GPU index for single run |

Disable consumer (keep all plotfiles):

```bash
CONSUME_PLOTFILES=0 bash grteclyn-wrapper/scripts/run_radialrecipe_gpu_smoke.sh
```

## Post-run plots (no GW / Psi4)

```bash
bash src/scripts/plot_diagnostic_radial.sh runs/radialrecipe_nonspherical/<episode_dir>
```

Writes EPS/PNG to `src/visualisation/plots/radial/` (constraints, collapse, shell profiles).

## Manual plotfile drain (if needed)

If a run finished before extraction completed:

```bash
bash src/scripts/plot_run_radial.sh runs/radialrecipe_nonspherical/<episode_dir> --no-delete
# or one-shot batch consume (no watch):
uv run python -m src.visualisation.process_wave.consume_plotfiles \
  --data runs/radialrecipe_nonspherical/<episode_dir> \
  --out runs/radialrecipe_nonspherical/<episode_dir>/small_data \
  --radii 4 8 --no-psi4 --shell-fields chi lapse K --areal-radius \
  --delete --keep-last 1 -j 4
```
