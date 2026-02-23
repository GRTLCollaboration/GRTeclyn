# process_wave — Streaming Psi4 extraction and plotting

Extract \(\Psi_4\) from plotfiles, append to a small `.dat` file, optionally delete plotfiles. Plot later from the `.dat` without needing plotfiles.

## Quick Start (Psi4 + frames, cleans K_z etc. before run)

```bash
python -m src.visualisation.process_wave.consume_plotfiles \
  --data "/home/jovyan/nachevsky/test/simulation/data_2gpu" \
  --out  "/home/jovyan/nachevsky/test/simulation/data_2gpu/small_data" \
  --radii 12 14 \
  --n-points 12 \
  --frames-fields chi K Weyl4_Re \
  --frames-axis z \
  --frames-out "/home/jovyan/nachevsky/test/simulation/GRTeclyn/src/visualisation/visualize" \
  --watch --delete --keep-last 2 \
  --verbose
```

## Requirements

Plotfiles must contain `Weyl4_Re` and `Weyl4_Im` (in params: `amr.derive_plot_vars = Weyl4`).

## 1. consume_plotfiles — Extract and optionally delete

**One-off (no delete):**
```bash
python -m src.visualisation.process_wave.consume_plotfiles \
  --data "/home/jovyan/nachevsky/test/simulation/data_2gpu" \
  --out  "/home/jovyan/nachevsky/test/simulation/data_2gpu/small_data" \
  --radii 10 15 --n-points 64 --stable-seconds 0 --keep-last 999999
```

**Watch + delete (run while simulation runs):**
```bash
python -m src.visualisation.process_wave.consume_plotfiles \
  --data "/home/jovyan/nachevsky/test/simulation/data_2gpu" \
  --out  "/home/jovyan/nachevsky/test/simulation/data_2gpu/small_data" \
  --radii 10 14 \
  --n-points 64 \
  --verbose \
  --watch \
  --delete \
  --keep-last 2
```

**Outputs:** `.../small_data/psi4_mode_l2m0.dat`, `consume_state.json`

### Rendering field frames while consuming plotfiles

You can render 2D slice frames (like `src.visualisation.visualize`) **while** consuming plotfiles, and still delete plotfiles to avoid disk overload.

Frames are written under:

- `.../src/visualisation/visualize/<field>_<axis>/frames/frame_<axis>_####.png`

**Psi4 + frames (chi, K, Weyl4_Re), then delete plotfiles:**

```bash
python -m src.visualisation.process_wave.consume_plotfiles \
  --data "/home/jovyan/nachevsky/test/simulation/data_2gpu" \
  --out  "/home/jovyan/nachevsky/test/simulation/data_2gpu/small_data" \
  --radii 12 14 \
  --n-points 12 \
  --frames-fields chi K Weyl4_Re \
  --frames-axis z \
  --frames-out "/home/jovyan/nachevsky/test/simulation/GRTeclyn/src/visualisation/visualize" \
  --watch --delete --keep-last 2 \
  --verbose
```

**Frames only (no Psi4 `.dat`), then delete plotfiles:**

```bash
python -m src.visualisation.process_wave.consume_plotfiles \
  --no-psi4 \
  --data "/home/jovyan/nachevsky/test/simulation/data_2gpu" \
  --frames-fields chi K Weyl4_Mag \
  --frames-axis z \
  --frames-out "/home/jovyan/nachevsky/test/simulation/GRTeclyn/src/visualisation/visualize" \
  --watch --delete --keep-last 2 \
  --verbose
```

### What `--keep-last` means

`--keep-last N` means **never delete the newest N plotfiles** in `--data`.

This is a safety buffer for live runs:
- avoids deleting a plotfile that might still be getting written
- keeps a couple of recent plotfiles for debugging

Example: `--keep-last 2` keeps the newest 2 `WormholePlt*` folders even when `--delete` is enabled.

### Restart behavior (auto-reset)

If you restart a simulation in the **same** `--data` folder (and plotfile indices start again at `...00000`), the consumer will automatically:
- truncate `psi4_mode_l2m0.dat`
- reset `consume_state.json`

This prevents accidental “append to old run” when reusing the same output directory.

## 2. plot_extracted_psi4 — Plot from .dat (no plotfiles)

```bash
python -m src.visualisation.process_wave.plot_extracted_psi4 \
  "/home/jovyan/nachevsky/test/simulation/data_2gpu/small_data/psi4_mode_l2m0.dat" \
  --time-axis retarded --out "/home/jovyan/nachevsky/test/simulation/GRTeclyn/src/visualisation/process_wave"
```

**Output:** `psi4_extracted_R10_R15.png`

## Options

### consume_plotfiles arguments

- `--data PATH`: directory containing `WormholePlt*` / `plt*`
- `--out PATH`: output directory (default: `<data>/small_data`)
- `--radii R1 R2 ...`: extraction radii
- `--n-points N`: angular resolution (\(N\times N\) sphere sampling grid)
- `--center x y z`: extraction center (default `0 0 0`)
- `--stable-seconds S`: only process plotfiles whose `Header` is older than `S` seconds
- `--poll-seconds S`: polling interval when `--watch` is enabled
- `--watch`: keep running and process new plotfiles as they appear
- `--verbose`: print per-plotfile details (name, time, delete/keep, timings)
- `--psi4/--no-psi4`: enable/disable writing `psi4_mode_l2m0.dat` (default: enabled)
- `--frames-fields F1 F2 ...`: render slice frames for these fields (e.g. `chi K Weyl4_Re`)
- `--frames-axis {x,y,z}`: slice axis for frames (default `z`)
- `--frames-coord VAL`: slice coordinate along that axis
- `--frames-zoom WIDTH`: plot width (code units)
- `--frames-center x y z`: explicit plot center
- `--frames-corner`: corner mode for symmetry-reduced domains
- `--frames-out PATH`: base output dir for frames (default: `src/visualisation/visualize`)
- `--delete`: delete processed plotfile directories after successful extraction
- `--keep-last N`: keep newest N plotfiles (don’t delete them)

**plot_extracted_psi4:** `input`, `--out`, `--radii`, `--time-axis`, `--t-min`, `--t-max`, `--psd-smooth-window`, `--psd-smooth-polyorder`, `--psd-hide-raw`

## Workflow

1. Set `plot_interval = 1` and `amr.derive_plot_vars = Weyl4` in params.
2. Start consumer (Psi4 + frames, cleans K_z etc. before run):

```bash
python -m src.visualisation.process_wave.consume_plotfiles \
  --data "/home/jovyan/nachevsky/test/simulation/data_2gpu" \
  --out  "/home/jovyan/nachevsky/test/simulation/data_2gpu/small_data" \
  --radii 12 14 \
  --n-points 12 \
  --frames-fields chi K Weyl4_Re \
  --frames-axis z \
  --frames-out "/home/jovyan/nachevsky/test/simulation/GRTeclyn/src/visualisation/visualize" \
  --watch --delete --keep-last 2 \
  --verbose
```

3. Run simulation in another terminal.
4. Plot: `python -m src.visualisation.process_wave.plot_extracted_psi4 <path/to/psi4_mode_l2m0.dat>`
