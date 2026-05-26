# GRTeclyn Python Wrapper

This wrapper runs isolated `SupportedWormholeCollapse` episodes from Python while keeping the C++ solver unchanged. Python owns episode setup, params generation, subprocess launch, diagnostics parsing, and scoring. `GRTeclyn` remains the numerical-relativity engine.

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

Run commands from the repository root.

Create an episode without launching the C++ executable:

```bash
python -m src.wrapper --dry-run --set stop_time=0.1 --set N_full=64 reproduce
```

Run one reproduction episode:

```bash
python -m src.wrapper --set stop_time=5 --set N_full=64 reproduce
```

Run a random sweep over existing wormhole parameters:

```bash
python -m src.wrapper --set stop_time=2 sweep --count 10 --seed 1
```

Use a specific executable or GPU assignment:

```bash
python -m src.wrapper \
  --executable Examples/SupportedWormholeCollapse/main3d.gnu.CUDA.ex \
  --cuda-devices 0 \
  reproduce
```

For MPI runs, set the rank count. The default executable name becomes `main3d.gnu.MPI.CUDA.ex`.

```bash
python -m src.wrapper --mpi-ranks 2 --cuda-devices 0,1 reproduce
```

## Plotfile Consumer

The wrapper can run the existing plotfile consumer as a side process:

```bash
python -m src.wrapper \
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

## Build Requirement

The wrapper does not build `GRTeclyn`. Build the example first, for example:

```bash
cd Examples/SupportedWormholeCollapse
make -j 8 USE_CUDA=TRUE USE_MPI=FALSE COMP=gnu CUDA_ARCH=90
```

Then run the wrapper from the repository root.
