# GRTeclyn Standalone Python Wrapper

This standalone wrapper runs isolated `SupportedWormholeCollapse` episodes from Python while keeping the C++ solver unchanged. Python owns episode setup, params generation, subprocess launch, diagnostics parsing, and scoring. `GRTeclyn` remains the numerical-relativity engine.

Set `GRTECLYN_ROOT` to your GRTeclyn checkout, or run commands from the GRTeclyn repository root:

```bash
export GRTECLYN_ROOT=/home/jovyan/nachevsky/test/simulation/GRTeclyn
```

## What It Creates

Each run gets its own episode directory under `runs/`:

```text
runs/episode_000001/
  params.txt
  metadata.json
  run.log
  score.json
  data/
    collapse_diagnostics.dat
    constraint_norms.dat
  small_data/
    psi4_mode_l2m0.dat
    areal_radius.dat
    consume_state.json
```

The wrapper rewrites only filesystem paths and requested parameter overrides. In particular, it sets:

```text
output_path
amr.check_file
amr.plot_file
```

to episode-local absolute paths so automated runs do not share output folders.

## Modules

- `config.py`: repository paths and executable name resolution.
- `episode.py`: episode directory creation and metadata writing.
- `params.py`: loads `params_2gpu.txt`, applies overrides, writes `params.txt`.
- `runner.py`: runs `check_params=1`, launches the solver, tees output to `run.log`, and can start `consume_plotfiles`.
- `metrics.py`: reads machine-readable diagnostics from an episode.
- `score.py`: computes an initial diagnostic score.
- `__main__.py`: command-line entry point.

## Basic Usage

Run commands from the GRTeclyn repository root, or set `GRTECLYN_ROOT`.

Create an episode without launching the C++ executable:

```bash
python -m grteclyn_wrapper --dry-run --set stop_time=0.1 --set N_full=64 reproduce
```

Run one reproduction episode:

```bash
python -m grteclyn_wrapper --set stop_time=5 --set N_full=64 reproduce
```

Run a random sweep over existing wormhole parameters:

```bash
python -m grteclyn_wrapper --set stop_time=2 sweep --count 10 --seed 1
```

Build a failure-atlas dataset from many cheap low-resolution runs:

```bash
python -m grteclyn_wrapper --cuda-devices 0 atlas --count 20 --seed 1
```

Use a specific executable or GPU assignment:

```bash
python -m grteclyn_wrapper \
  --executable Examples/SupportedWormholeCollapse/main3d.gnu.CUDA.ex \
  --cuda-devices 0 \
  reproduce
```

For MPI runs, set the rank count. The default executable name becomes `main3d.gnu.MPI.CUDA.ex`.

```bash
python -m grteclyn_wrapper --mpi-ranks 2 --cuda-devices 0,1 reproduce
```

## Shell Scripts

The `scripts/` directory contains ready-to-run commands for the staged flow tests:

```bash
./scripts/run_medium_flow.sh
./scripts/run_big_single_gpu.sh
./scripts/run_full_single_gpu.sh
./scripts/run_atlas_gpu.sh
./scripts/dry_run_atlas.sh
```

The scripts default to GPU `0`. Override the GPU without editing files:

```bash
CUDA_VISIBLE_DEVICES_OVERRIDE=1 ./scripts/run_medium_flow.sh
```

Atlas scripts also accept simple environment overrides:

```bash
ATLAS_COUNT=100 ATLAS_SEED=42 ./scripts/run_atlas_gpu.sh
```

## Plotfile Consumer

The wrapper can run the existing plotfile consumer as a side process:

```bash
python -m grteclyn_wrapper \
  --consume-plotfiles \
  --consumer-radii 8 16 \
  reproduce
```

This calls:

```bash
python -m src.visualisation.process_wave.consume_plotfiles
```

with episode-local `--data` and `--out` paths. Use `--consumer-delete` if old plotfiles should be removed while the run is active.

## Scoring

After a simulation exits, the wrapper reads:

- `data/collapse_diagnostics.dat`
- `data/constraint_norms.dat`

and writes `score.json`. The first score is deliberately simple:

- rewards survival toward `stop_time`;
- rewards bounded Hamiltonian and momentum constraints;
- rewards nontrivial geometry;
- penalizes horizon-proxy formation;
- penalizes severe lapse collapse.

This score is a starting point for sweeps and optimizer development, not a final physics validation metric.

## Failure Atlas

The `atlas` command is the dataset-generating search mode. It samples existing wormhole parameters, runs one isolated episode per sample, classifies the result, and writes batch-level files under `runs/atlas_<timestamp>/`.

In simple terms, atlas is a batch experiment runner:

1. Pick random wormhole parameters.
2. Create a separate episode folder.
3. Write a custom `params.txt`.
4. Run `GRTeclyn`.
5. Read diagnostics such as constraints, lapse, `chi`, and horizon proxy.
6. Compute a score.
7. Add labels such as `completed`, `horizon_formed`, or `solver_failed`.
8. Save the result to `atlas.csv` and `atlas.jsonl`.

The CSV is the quick table:

```text
parameters -> what happened -> score -> success/failure labels
```

Failures are useful. If a candidate forms a horizon and crashes with NaNs, atlas records that parameter region as bad or unstable. Later optimizers can use this failure map instead of blindly trying the same bad regions again.

```text
runs/atlas_YYYYMMDDTHHMMSSZ/
  metadata.json
  atlas.jsonl
  atlas.csv
  summary.json
  episode_000001/
  episode_000002/
```

By default, atlas runs use cheap screening settings:

```text
N_full = 32
max_level = 0
stop_time = 0.04
plot_interval = 1000
checkpoint_interval = 1000
```

Each sampled candidate varies:

- `wormhole_phi_perturbation_amplitude`
- `wormhole_support_strength`
- `wormhole_phi_perturbation_width`
- `wormhole_phi_monopole_amplitude`

You can override any sampled or fixed key with `--set`. For example:

```bash
python -m grteclyn_wrapper \
  --cuda-devices 0 \
  --set N_full=64 \
  --set stop_time=0.2 \
  atlas --count 5 --seed 10
```

Preview sampled candidates without launching `GRTeclyn`:

```bash
python -m grteclyn_wrapper --dry-run atlas --count 3 --seed 1
```

Atlas labels are intentionally simple and may be combined:

- `completed`
- `missing_diagnostics`
- `constraint_blowup`
- `lapse_collapse`
- `horizon_formed`
- `trivial_geometry`
- `solver_failed`

Use `atlas.jsonl` for complete nested records and `atlas.csv` for quick spreadsheet-style inspection.

## Build Requirement

The wrapper does not build `GRTeclyn`. Build the example first, for example:

```bash
cd Examples/SupportedWormholeCollapse
make -j 8 USE_CUDA=TRUE USE_MPI=FALSE COMP=gnu CUDA_ARCH=90
```

Then run the wrapper from the repository root.
