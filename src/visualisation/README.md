# GRTeclyn Visualization

Scripts for visualizing GRTeclyn simulation results (plotfiles in BoxLib/AMReX format).

## Components at a glance

| Location | Role |
|----------|------|
| **`visualize/`** | 2D field slices from plotfiles; optional MP4 via `ffmpeg` or separate stitcher. |
| **`visualize/make_movies.py`** | Stitch existing PNG frame folders into MP4 (slice or embedding layouts). |
| **`make_evolution_panel/`** | Multi-panel strip figure from saved frame PNGs. |
| **`extract_wave/`** | Extract \(\Psi_4\) from plotfiles and plot time series + PSD. |
| **`process_wave/`** | Stream plotfiles to `psi4_mode_l2m0.dat`, render frames while consuming, plot from `.dat`, LIGO/strain propagation extras. |
| **`diagnostic/`** | Collapse diagnostics multi-panel plot (`collapse_diagnostics.dat`, optional `areal_radius.dat`). |
| **`constraines/`** | Constraint norms \(L_2\) of Hamiltonian and momentum (`constraint_norms.dat`). |
| **`figures/`** | Standalone publication-style schematic figures (not driven by simulation dumps). |
| **`src/scripts/`** | Shell automation: live plotfile processing (`plot_run.sh`), post-run figures (`plot_diagnostic.sh`), archive to `SimResults/` (`move_files.sh`). |

More detail for live processing and `consume_plotfiles` options: [`process_wave/README.md`](process_wave/README.md).

## Prerequisites

1. Install dependencies (from the project root):

   ```bash
   uv sync
   # or: source .venv/bin/activate && pip install -e .
   ```

2. Required packages: `yt`, `numpy`, `matplotlib`, `scipy` (from `pyproject.toml`). Optional for MPI: `uv sync --extra mpi` or `pip install mpi4py`.

3. For animations: `ffmpeg` must be installed.

4. Simulation data: plotfiles (`plt*` / `WormholePlt*`) in the default location `../data` relative to the project root. If your run writes elsewhere (e.g. `output_path = ".../data_2gpu"`), pass `--data` to the relevant tool.

---

## 1. `visualize/` — Field slices and movies

Plots slice views of GR fields (`yt`). Choose slice axis and coordinate, zoom, and symmetry/corner mode. Default output: `src/visualisation/visualize/`.

**Available fields**

| Field | Role |
|-------|------|
| `chi` | Conformal factor |
| `K` | Trace of extrinsic curvature |
| `Theta` | Z4 constraint |
| `lapse` | Lapse \(\alpha\) |
| `Ham` | Hamiltonian constraint |
| `A11`, `A12`, `A22`, `A33` | Extrinsic curvature components |
| `GW_Plus`, `GW_Cross` | GW strain proxies from \(A_{ij}\) (no Weyl4 required) |

**Basic usage**

```bash
# Last timestep only (default field: chi)
python -m src.visualisation.visualize --field K

# Animation from all plotfiles
python -m src.visualisation.visualize --field K --animate --zoom 100

# Custom data directory (e.g. multi-GPU run)
python -m src.visualisation.visualize --field K --animate --zoom 100 --data /path/to/data_2gpu
```

| Option | Default | Description |
|--------|---------|-------------|
| `--field` | `chi` | Field to plot |
| `--axis` | `z` | Normal to slice: `x`, `y`, `z` |
| `--coord` | — | Slice coordinate along that axis |
| `--animate` | off | Build MP4 after writing all frames |
| `--zoom` | full domain | Width in code units |
| `--data` | `../data` | Directory containing plotfolders |
| `--out` | `visualize/` package dir | Base output directory |
| `--center` / `--corner` | — | Full domain center vs symmetry corner mode |
| `--mpi` | off | Parallel frame generation with `mpi4py` |

**MPI**

```bash
mpirun -np 8 python -m src.visualisation.visualize --field K --animate --zoom 100 --mpi
```

**Output layout:** `<out>/<field>_<axis>/frames/frame_<axis>_NNNN.png` and, if `--animate`, `movie_<field>_<axis>.mp4`.

### 1a. `visualize/make_movies.py` — Stitch frames to MP4

Use when you already have PNGs (e.g. from `visualize` without `--animate`, or from `consume_plotfiles` frames, or `embedding/frames/`). Scans `--root` for:

- `<field>_<axis>/frames/frame_<axis>_NNNN.png` → `movie_<field>_<axis>.mp4`
- `<name>/frames/frame_NNNN.png` (e.g. embedding) → `movie_<name>.mp4`

```bash
# Default root: directory containing this script (visualize/)
python -m src.visualisation.visualize.make_movies

python -m src.visualisation.visualize.make_movies --root /path/to/visualize --framerate 10

# Only selected folders under root
python -m src.visualisation.visualize.make_movies --only K_z chi_z embedding
```

Requires `ffmpeg` on `PATH` (uses `subprocess`, not `os.system`).

---

## 2. `make_evolution_panel/` — Strip figure from frames

Builds a horizontal panel (one column per timestep) from `frame_z_NNNN.png` or `frame_NNNN.png` in `--frame_dir`. Writes `.png` and `.pdf` under `src/visualisation/plots/` when `--out` is relative.

```bash
python src/visualisation/make_evolution_panel/make_evolution_panel.py \
  --frame_dir src/visualisation/visualize/K_z/frames \
  --out evolution_K_z_panel \
  --frames 0 500 1000 1500

python src/visualisation/make_evolution_panel/make_evolution_panel.py \
  --frame_dir src/visualisation/visualize/K_z/frames \
  --out evolution_K_z_gray \
  --frames 0 20 40 60 \
  --grayscale

# Absolute --out writes exactly there
python src/visualisation/make_evolution_panel/make_evolution_panel.py \
  --frame_dir src/visualisation/visualize/K_z/frames \
  --out /tmp/my_panel \
  --frames 0 100 200
```

| Option | Description |
|--------|-------------|
| `--frame_dir` | Folder with frame PNGs |
| `--out` | Output stem (`.png`/`.pdf` added); relative → under `src/visualisation/plots/` |
| `--frames` | Indices (default `0 20 40 60`) |
| `--grayscale` | Grayscale panels |
| `--no-grid` | Disable dashed overlay grid |

---

## 3. `extract_wave/` — \(\Psi_4\) from plotfiles

Extracts the \(l=2,m=0\) mode of \(\Psi_4\) on spheres at chosen radii and plots waveform + PSD. Does **not** integrate to strain \(h(t)\) (ill-conditioned); frequency content comes from \(\Psi_4\) in the frequency domain. Supports symmetry-reduced domains via partial-sphere integration normalized to \(4\pi\).

`python -m src.visualisation.extract_wave` forwards to `plot_psi4.main`.

```bash
python -m src.visualisation.extract_wave --radii 14 30 --n-points 24

python -m src.visualisation.extract_wave --radii 14 30 --n-points 24 --time-axis retarded

python -m src.visualisation.extract_wave --data /path/to/data_2gpu --radii 14 30 --n-points 24
```

| Option | Description |
|--------|-------------|
| `--data` | Plotfile directory |
| `--radii` | Extraction radii |
| `--n-points` | Angular grid \(N \times N\) |
| `--out` | Output directory (default: `extract_wave/`) |
| `--time-axis` | `simulation` or `retarded` (\(t - R\)) |
| `--t-min`, `--t-max` | Plot range |
| `--psd-smooth-window`, `--psd-smooth-polyorder` | PSD smoothing |
| `--psd-hide-raw` | Hide raw PSD points |

**Output:** e.g. `psi4_analysis_R14.0_R30.0_n24.png`.

**Retarded time:** negative abscissa values mean \(t < R_{\rm ext}\): usually junk / near-zone, not clean wavezone GW.

---

## 4. `process_wave/` — Streaming extraction and post-processing

- **`consume_plotfiles`**: while a run produces plotfiles, extract \(\Psi_4\) to `small_data/psi4_mode_l2m0.dat`, optionally render slice frames under `visualize/`, optionally delete processed plotfiles, optionally areal radius / embedding. See [`process_wave/README.md`](process_wave/README.md).
- **`plot_extracted_psi4`**: plot from an existing `.dat` without plotfiles (waveform + PSD; optional strain, LIGO overlay, propagation speed).

**Quick start (watch + delete + frames)**

```bash
python -m src.visualisation.process_wave.consume_plotfiles \
  --data /path/to/data_2gpu \
  --out /path/to/data_2gpu/small_data \
  --radii 8 12 16 --n-points 32 \
  --frames-fields chi K Weyl4_Re Weyl4_Mag \
  --frames-axis z --frames-corner \
  --frames-out "$(pwd)/src/visualisation/visualize" \
  --watch --delete --keep-last 2 --verbose
```

**Plot from extracted `.dat`**

```bash
python -m src.visualisation.process_wave.plot_extracted_psi4 \
  /path/to/small_data/psi4_mode_l2m0.dat \
  --radii 10 14 --time-axis retarded --plot-psd \
  --out "$(pwd)/src/visualisation/process_wave"
```

**Advanced analysis** (strain PSD vs LIGO, propagation speed across radii):

```bash
python -m src.visualisation.process_wave.plot_extracted_psi4 \
  /path/to/small_data/psi4_mode_l2m0.dat \
  --plot-psd --strain --propagation-speed --mass-msun 30 --distance-mpc 10
```

Strain scaling: \(|\tilde{h}| = |\tilde{\Psi}_4| / (2\pi f)^2\); characteristic strain and detector overlays follow the implementation in `plot_extracted_psi4.py`.

**Areal radius + embedding** (when consuming plotfiles): add `--areal-radius`, `--embedding`, `--embedding-rmax`; embedding frames go under `visualize/embedding/frames/` and match the layout expected by `make_movies.py`.

---

## 5. `diagnostic/` — Collapse diagnostics

Reads `collapse_diagnostics.dat` (SmallDataIO ASCII under `<run>/data/` or path you pass). Optionally overlays `areal_radius.dat` (e.g. from `consume_plotfiles --areal-radius`) for throat radius, expansion velocity, and K-decay fit.

```bash
# Infer collapse_diagnostics.dat from run directory
python -m src.visualisation.diagnostic.diagnostic --data /path/to/data_2gpu

# Explicit paths
python -m src.visualisation.diagnostic.diagnostic \
  /path/to/data_2gpu/data/collapse_diagnostics.dat \
  --data /path/to/data_2gpu \
  --areal-radius-file /path/to/small_data/areal_radius.dat

python -m src.visualisation.diagnostic.diagnostic --data /path/to/data_2gpu --no-fit-lifetime --mass-msun 30
```

| Option | Description |
|--------|-------------|
| `input` (positional) | Optional explicit `collapse_diagnostics.dat` |
| `--data` | Run root (used to find diagnostics and default `areal_radius.dat`) |
| `--out` | Output directory |
| `--name` | Output filename (default `collapse_diagnostics_plot.eps`; also `.png`/`.pdf` as implemented) |
| `--areal-radius-file` | Override auto-detected `areal_radius.dat` |
| `--no-fit-lifetime` | Skip exponential fit to K decay |
| `--mass-msun` | Mass for physical time conversion |

---

## 6. `constraines/` — Constraint norms

Plots \(L_2\) norms of Hamiltonian and momentum constraints from `constraint_norms.dat` (typically `<run>/data/constraint_norms.dat`).

```bash
python -m src.visualisation.constraines /path/to/data/data/constraint_norms.dat

python -m src.visualisation.constraines /path/to/constraint_norms.dat -o constraints_plot.eps
```

If you pass only `-o myplot.eps` (no directory component), outputs go under `src/visualisation/constraines/`. The module also writes matching `.png` and `.pdf` next to the EPS stem.

---

## 7. `figures/` — Paper figures (schematic)

Helper scripts for publication-style visuals **not** tied to a single run’s plotfiles. Example: `plot_collapse_stages.py` builds a 2×3 panel of schematic embedding-style wormhole stages (collapse and expansion rows); saves `wormhole_collapse_stages.png` and `.pdf`.

```bash
python src/visualisation/figures/plot_collapse_stages.py
```

The script’s `__main__` block sets the output directory; edit that path if you want a different destination.

---

## 8. GW proxies without Weyl4

From extrinsic curvature on a slice (approximate \(h_+\), \(h_\times\) for propagation along \(z\)):

```bash
python -m src.visualisation.visualize --field GW_Plus --axis y --zoom 100 --animate
python -m src.visualisation.visualize --field GW_Cross --axis z --animate --zoom 40
```

---

## 9. Symmetry-reduced domains

For octant-style runs (`lo_boundary = 2 2 2`, etc.), use `--center` and/or `--corner` so the origin is readable on the slice.

```bash
python -m src.visualisation.visualize --field K --axis z --coord 0 --zoom 84 --corner --animate
```

---

## Output layout

```
src/visualisation/
├── visualize/
│   ├── <field>_<axis>/
│   │   ├── frames/
│   │   └── movie_<field>_<axis>.mp4
│   └── embedding/frames/     # from consume_plotfiles --embedding
├── plots/                    # diagnostic.sh outputs; evolution panels (default relative out)
├── extract_wave/             # default --out for plot_psi4
├── process_wave/             # plots from plot_extracted_psi4
├── diagnostic/               # optional direct output from diagnostic.py
├── constraines/              # default location for relative -o
├── figures/                  # schematic figure outputs
└── README.md

<run_dir>/small_data/
├── psi4_mode_l2m0.dat
└── areal_radius.dat
```

---

## Automation scripts

Shell helpers live under **`src/scripts/`**. Run them from the **repository root** (`GRTeclyn/`) so `python -m src.visualisation...` resolves correctly. Paths below use `./src/scripts/...`.

### `plot_run.sh` — Live processing during a simulation

Watches the run directory for new plotfiles, runs `consume_plotfiles` to extract \(\Psi_4\), optional throat **areal radius**, **3D embedding** frames, and **slice frames** (`K`, `Weyl4_Re` on the \(z\) slice with `--frames-corner`). Processed plotfile folders can be **deleted** to save disk (`--keep-last 2` retains the newest two folders).

**Data directory:** If you omit the path, the script picks the first existing directory among (parent of the repo): `data_2gpu`, then `data_supported`, then `data`. If none exist, it defaults to `data_2gpu` under that parent.

**Flags**

| Flag | Meaning |
|------|---------|
| `--not-remove` | Do **not** wipe stale plotfiles, `small_data` extraction state, run diagnostics ASCII, or previously generated images under `visualisation/`; pass `--keep-existing-frames` to the consumer. |
| `-j N` / `--jobs N` | Parallel workers for `consume_plotfiles` (default **64**). |

**Commands**

```bash
# Default data dir (auto-detect) and default 64 jobs
./src/scripts/plot_run.sh

# Explicit run/output folder
./src/scripts/plot_run.sh /path/to/data_2gpu

# Keep existing plotfiles and frames (e.g. resume processing)
./src/scripts/plot_run.sh --not-remove /path/to/data_2gpu

# Fewer parallel workers
./src/scripts/plot_run.sh -j 16 /path/to/data_2gpu

./src/scripts/plot_run.sh --not-remove -j 16
```

**What `plot_run` passes to `consume_plotfiles`** (fixed in the script): extraction radii `12 16 20 24`, `--n-points 64`, `--areal-radius`, `--embedding --embedding-rmax 5.0`, `--frames-fields K Weyl4_Re`, `--frames-axis z`, `--frames-corner`, `--frames-out` → `src/visualisation/visualize`, `--watch --delete --keep-last 2`, `--verbose`.

**Movies:** This script does **not** call `ffmpeg`. After frames exist under `src/visualisation/visualize/`, build MP4s with:

```bash
python -m src.visualisation.visualize.make_movies --root src/visualisation/visualize
```

---

### `plot_diagnostic.sh` — Post-run plots (one folder)

After the run has written `data/constraint_norms.dat`, `data/collapse_diagnostics.dat`, and `small_data/psi4_mode_l2m0.dat`, this script **removes and recreates** `src/visualisation/plots/` and writes fresh figures there.

**Arguments:** optional `RUN_DIR`, then optional extraction **radii** passed through to `plot_extracted_psi4` (if omitted, all radii in the `.dat` file are used).

**Default `RUN_DIR`:** Among `data_2gpu`, `data_supported`, and `data` (sibling of the repo), chooses the directory whose three required files are present and have the **newest combined modification time**; if that fails, falls back to an existing `data_2gpu` or `data`.

**Outputs in `src/visualisation/plots/`**

| File stem | Content |
|-----------|---------|
| `constraints_plot.*` | Hamiltonian and momentum constraint norms |
| `collapse_diagnostics_plot.*` | Collapse diagnostics (+ areal radius and K-decay fit when data exist) |
| `psi4_analysis.*` | Combined panel (`--combined`): waveforms, PSD, propagation, strain, LIGO overlay (`--strain`, 30 M⊙ at 10 Mpc) |

```bash
./src/scripts/plot_diagnostic.sh

./src/scripts/plot_diagnostic.sh /path/to/data_2gpu

./src/scripts/plot_diagnostic.sh /path/to/data_2gpu 12 16 20 24
```

---

### `move_files.sh` — Archive run + visuals to `SimResults/`

Copies visualization frames (`chi*`, `K*`, `Weyl4_Mag*`, `Weyl4_Re*`, `embedding*`), everything from `src/visualisation/plots/`, key ASCII data from the run tree, and the chosen parameter file into `src/SimResults/<auto_folder>/`. The folder name is derived from wormhole parameters read from the params file (e.g. `Run_R..._A0..._A2..._sigma...`).

**First argument — run type** (optional, default `SupportedWormholeCollapse`):

| Value | Data directory | Params file |
|-------|----------------|-------------|
| `SupportedWormholeCollapse` (default) | `<repo>/../data_supported` | `Examples/SupportedWormholeCollapse/params_2gpu.txt` |
| Any other string (e.g. `WormholeCollapse`) | `<repo>/../data_2gpu` | `Examples/WormholeCollapse/params_2gpu.txt` |

```bash
# Supported wormhole run + data_supported
./src/scripts/move_files.sh
./src/scripts/move_files.sh SupportedWormholeCollapse

# Other run type (uses data_2gpu and WormholeCollapse params)
./src/scripts/move_files.sh WormholeCollapse
```

If `src/visualisation/plots/` is missing, the script warns to run `plot_diagnostic.sh` first.

---

### Suggested workflow

```bash
# Terminal 1 — while the simulation runs
./src/scripts/plot_run.sh /path/to/your_run_output

# After the run — refresh diagnostic figures
./src/scripts/plot_diagnostic.sh /path/to/your_run_output

# Optional — stitch frame folders to MP4
python -m src.visualisation.visualize.make_movies --root src/visualisation/visualize

# Archive
./src/scripts/move_files.sh SupportedWormholeCollapse
```

---

## Custom paths (examples)

```bash
python -m src.visualisation.visualize --field K --data /path/to/plt --out /path/to/output
```

```bash
python -m src.visualisation.constraines /path/to/data/data/constraint_norms.dat
```

```bash
python -m src.visualisation.visualize.make_movies --root /path/to/visualize
```

Archive runs: see **`move_files.sh`** in [Automation scripts](#automation-scripts) above.
