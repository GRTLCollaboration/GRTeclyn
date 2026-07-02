# GRTeclyn Visualization

Scripts for visualizing GRTeclyn simulation results (plotfiles in BoxLib/AMReX format).

## Quick start — post-run figures in the run directory

From the **GRTeclyn repo root**, after a run has produced the required ASCII files:

```bash
bash grteclyn-wrapper/scripts/plot/plot_diagnostic.sh RUN_DIR [RADIUS ...]
```

Example (HQ promotion, Q-ball at box center, Psi4 radii 8/12/24):

```bash
bash grteclyn-wrapper/scripts/plot/plot_diagnostic.sh \
  runs/grtresna_promote/qball_traj_spiral_v2_t30_hq_eval000118 \
  8 12 24
```

**Output directory:** `<RUN_DIR>/plots/` (cleared and recreated each invocation).

| Figure | Input data | Module |
|--------|------------|--------|
| `constraints_plot.*` | `data/constraint_norms.dat` (`t`, `L2_Ham`, `L2_Mom`) | `constraines/` |
| `collapse_diagnostics_plot.*` | `data/collapse_diagnostics.dat` (+ optional `small_data/areal_radius.dat`) | `diagnostic/` |
| `psi4_analysis_M1000_D1.*` etc. | `small_data/psi4_mode_l2m0.dat` | `process_wave/plot_extracted_psi4.py` (`--combined --strain`) |

Frame movies: `bash grteclyn-wrapper/scripts/plot/make_movies.sh RUN_DIR --framerate 10` → `<RUN_DIR>/movies/`.

**Not the same metric:** `grtresna/Ham_and_Mom_errors.txt` reports GRTresna percent errors on the
initial-data solve; `constraints_plot` uses GPU `constraint_norms.dat` L2 norms during evolution.

See also the wrapper overview: [`grteclyn-wrapper/README.md`](../../../README.md#visualization).

---

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
| **`search/`** | QD / search campaign analytics from `trajectory.jsonl` (batch improvement, saturation). |
| **`grteclyn-wrapper/scripts/plot/`** | Shell automation: live plotfile processing (`plot_run.sh`), post-run figures (`plot_diagnostic.sh`), frame movies (`make_movies.sh`). |
| **`grteclyn-wrapper/scripts/wormhole/`** | Archive wormhole runs + visuals to `SimResults/` (`move_files.sh`). |

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

Plots slice views of GR fields (`yt`). Choose slice axis and coordinate, zoom, and symmetry/corner mode. Default output: `grteclyn-wrapper/src/grteclyn_wrapper/visualisation/visualize/`.

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
uv run python -m grteclyn_wrapper.visualisation.visualize --field K

# Animation from all plotfiles
uv run python -m grteclyn_wrapper.visualisation.visualize --field K --animate --zoom 100

# Custom data directory (e.g. multi-GPU run)
uv run python -m grteclyn_wrapper.visualisation.visualize --field K --animate --zoom 100 --data /path/to/data_2gpu
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
mpirun -np 8 uv run python -m grteclyn_wrapper.visualisation.visualize --field K --animate --zoom 100 --mpi
```

**Output layout:** `<out>/<field>_<axis>/frames/frame_<axis>_NNNN.png` and, if `--animate`, `movie_<field>_<axis>.mp4`.

### 1a. `visualize/make_movies.py` — Stitch frames to MP4

Use when you already have PNGs (e.g. from `visualize` without `--animate`, or from `consume_plotfiles` frames, or `embedding/frames/`). Scans `--root` for:

- `<field>_<axis>/frames/frame_<axis>_NNNN.png` → `movie_<field>_<axis>.mp4`
- `<name>/frames/frame_NNNN.png` (e.g. embedding) → `movie_<name>.mp4`

```bash
# Default root: directory containing this script (visualize/)
uv run python -m grteclyn_wrapper.visualisation.visualize.make_movies

uv run python -m grteclyn_wrapper.visualisation.visualize.make_movies --root /path/to/visualize --framerate 10

# Only selected folders under root
uv run python -m grteclyn_wrapper.visualisation.visualize.make_movies --only K_z chi_z embedding
```

Requires `ffmpeg` on `PATH` (uses `subprocess`, not `os.system`).

---

## 2. `make_evolution_panel/` — Strip figure from frames

Builds a horizontal panel (one column per timestep) from `frame_z_NNNN.png` or `frame_NNNN.png` in `--frame_dir`. Writes `.png` and `.pdf` under `grteclyn-wrapper/src/grteclyn_wrapper/visualisation/plots/` when `--out` is relative.

```bash
# K_z frames (default mode: k_z)
python grteclyn-wrapper/src/grteclyn_wrapper/visualisation/make_evolution_panel/make_evolution_panel.py \
  --frame_dir grteclyn-wrapper/src/grteclyn_wrapper/visualisation/visualize/K_z/frames \
  --out evolution_K_z_panel \
  --frames 0 500 1000 1500

python grteclyn-wrapper/src/grteclyn_wrapper/visualisation/make_evolution_panel/make_evolution_panel.py \
  --frame_dir grteclyn-wrapper/src/grteclyn_wrapper/visualisation/visualize/K_z/frames \
  --out evolution_K_z_gray \
  --frames 0 20 40 60 \
  --grayscale

# Embedding frames (3D embedding snapshots): use --mode embedding
# (expects frame_NNNN.png in the embedding frames folder)
python grteclyn-wrapper/src/grteclyn_wrapper/visualisation/make_evolution_panel/make_evolution_panel.py \
  --frame_dir grteclyn-wrapper/src/grteclyn_wrapper/visualisation/visualize/embedding/frames \
  --mode embedding \
  --out evolution_embedding_panel \
  --frames 2 3001 4002

# Keep the title area ("Embedding Diagram  t=...") in each panel
python grteclyn-wrapper/src/grteclyn_wrapper/visualisation/make_evolution_panel/make_evolution_panel.py \
  --frame_dir grteclyn-wrapper/src/grteclyn_wrapper/visualisation/visualize/embedding/frames \
  --mode embedding \
  --keep-title \
  --out evolution_embedding_with_title \
  --frames 2 3001 4002

# 4 panels: just provide 4 frame indices
python grteclyn-wrapper/src/grteclyn_wrapper/visualisation/make_evolution_panel/make_evolution_panel.py \
  --frame_dir grteclyn-wrapper/src/grteclyn_wrapper/visualisation/visualize/embedding/frames \
  --mode embedding \
  --keep-title \
  --out evolution_embedding_4panels \
  --frames 2 1002 3001 4002

# Absolute --out writes exactly there
python grteclyn-wrapper/src/grteclyn_wrapper/visualisation/make_evolution_panel/make_evolution_panel.py \
  --frame_dir grteclyn-wrapper/src/grteclyn_wrapper/visualisation/visualize/K_z/frames \
  --out /tmp/my_panel \
  --frames 0 100 200
```

| Option | Description |
|--------|-------------|
| `--frame_dir` | Folder with frame PNGs |
| `--mode` | Frame type/layout: `k_z` (default) or `embedding` |
| `--out` | Output stem (`.png`/`.pdf` added); relative → under `grteclyn-wrapper/src/grteclyn_wrapper/visualisation/plots/` |
| `--frames` | Indices (default `0 20 40 60`) |
| `--grayscale` | Grayscale panels |
| `--no-grid` | Disable dashed overlay grid (K_z mode only) |
| `--keep-title` | Keep title area at the top (embedding mode only) |

---

## 3. `extract_wave/` — \(\Psi_4\) from plotfiles

Extracts the \(l=2,m=0\) mode of \(\Psi_4\) on spheres at chosen radii and plots waveform + PSD. Does **not** integrate to strain \(h(t)\) (ill-conditioned); frequency content comes from \(\Psi_4\) in the frequency domain. Supports symmetry-reduced domains via partial-sphere integration normalized to \(4\pi\).

`uv run python -m grteclyn_wrapper.visualisation.extract_wave` forwards to `plot_psi4.main`.

```bash
uv run python -m grteclyn_wrapper.visualisation.extract_wave --radii 14 30 --n-points 24

uv run python -m grteclyn_wrapper.visualisation.extract_wave --radii 14 30 --n-points 24 --time-axis retarded

uv run python -m grteclyn_wrapper.visualisation.extract_wave --data /path/to/data_2gpu --radii 14 30 --n-points 24
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
uv run python -m grteclyn_wrapper.visualisation.process_wave.consume_plotfiles \
  --data /path/to/data_2gpu \
  --out /path/to/data_2gpu/small_data \
  --radii 8 12 16 --n-points 32 \
  --frames-fields chi K Weyl4_Re Weyl4_Mag \
  --frames-axis z --frames-corner \
  --frames-out "$(pwd)/grteclyn-wrapper/src/grteclyn_wrapper/visualisation/visualize" \
  --watch --delete --keep-last 2 --verbose
```

**Plot from extracted `.dat` (stacked waveform + PSD)**

```bash
uv run python -m grteclyn_wrapper.visualisation.process_wave.plot_extracted_psi4 \
  /path/to/small_data/psi4_mode_l2m0.dat \
  --radii 10 14 --time-axis retarded --plot-psd \
  --out /path/to/run/plots
```

**Article-style 6-panel GW analysis** (same layout as wormhole paper figures:
waveforms, retarded-time + QNM fit, PSD, propagation speed, spectrogram, LIGO strain):

```bash
uv run python -m grteclyn_wrapper.visualisation.process_wave.plot_extracted_psi4 \
  /path/to/small_data/psi4_mode_l2m0.dat \
  --radii 8 12 24 \
  --combined --strain \
  --mass-msun 1000 --distance-mpc 1 \
  --out /path/to/run/plots \
  --name psi4_analysis_M1000_D1.eps
```

Panels (a)–(f): simulation-time waveforms; retarded-time waveforms + damped-sinusoid QNM fit;
\(\Psi_4\) PSD; peak-tracking propagation speed (needs \(\geq 2\) radii); spectrogram at
`--spectrogram-radius` (default: innermost radius); strain vs Advanced LIGO at `--mass-msun` /
`--distance-mpc`. Use `--ligo-quantity asd` (default) or `hchar`.

**Propagation / QNM caveats:** meaningful results need enough timesteps in the `.dat` file and a
clear wave burst. Early-run snapshots (few plotfiles consumed) produce empty or misleading panels.

Strain scaling: \(|\tilde{h}| = |\tilde{\Psi}_4| / (2\pi f)^2\); characteristic strain and detector
overlays follow the implementation in `plot_extracted_psi4.py`.

**Areal radius + embedding** (when consuming plotfiles): add `--areal-radius`, `--embedding`, `--embedding-rmax`; embedding frames go under `visualize/embedding/frames/` and match the layout expected by `make_movies.py`.

---

## 5. `diagnostic/` — Collapse diagnostics

Reads `collapse_diagnostics.dat` (SmallDataIO ASCII under `<run>/data/` or path you pass). Optionally overlays `areal_radius.dat` (e.g. from `consume_plotfiles --areal-radius`) for throat radius, expansion velocity, and K-decay fit.

```bash
# Infer collapse_diagnostics.dat from run directory
uv run python -m grteclyn_wrapper.visualisation.diagnostic.diagnostic --data /path/to/data_2gpu

# Explicit paths
uv run python -m grteclyn_wrapper.visualisation.diagnostic.diagnostic \
  /path/to/data_2gpu/data/collapse_diagnostics.dat \
  --data /path/to/data_2gpu \
  --areal-radius-file /path/to/small_data/areal_radius.dat

uv run python -m grteclyn_wrapper.visualisation.diagnostic.diagnostic --data /path/to/data_2gpu --no-fit-lifetime --mass-msun 30
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

Plots \(L_2\) norms of Hamiltonian and momentum constraints from `constraint_norms.dat`
(typically `<run>/data/constraint_norms.dat`). Two log-scale panels: \(\|\mathcal{H}\|_{L^2}\) and
\(\|\mathcal{M}\|_{L^2}\) vs \(t\,[M]\) — same style as wormhole article `constraints_plot` figures.

Requires `calculate_constraint_norms=1` during GPU evolution (wrapper campaigns enable this).

```bash
uv run python -m grteclyn_wrapper.visualisation.constraines \
  /path/to/run/data/constraint_norms.dat \
  -o /path/to/run/plots/constraints_plot.eps
```

If you pass only `-o myplot.eps` (no directory component), outputs go under
`grteclyn-wrapper/src/grteclyn_wrapper/visualisation/constraines/`. The module also writes matching
`.png` and `.pdf` next to the EPS stem.

**Not plotted here:** `grtresna/Ham_and_Mom_errors.txt` (GRTresna percent errors on initial data).

---

## 7. `figures/` — Paper figures (schematic)

Helper scripts for publication-style visuals **not** tied to a single run’s plotfiles. Example: `plot_collapse_stages.py` builds a 2×3 panel of schematic embedding-style wormhole stages (collapse and expansion rows); saves `wormhole_collapse_stages.png` and `.pdf`.

```bash
python grteclyn-wrapper/src/grteclyn_wrapper/visualisation/figures/plot_collapse_stages.py
```

The script’s `__main__` block sets the output directory; edit that path if you want a different destination.

---

## 7a. `search/` — QD campaign batch progress

Marginal archive score gain per GPU batch (bar) plus cumulative best (saturation curve). Reads `trajectory.jsonl` from a finished or in-progress MAP-Elites run.

```bash
uv run python -m grteclyn_wrapper.visualisation.search \
  runs/grtresna_qd/ftl_4d/ftl_4d_v1

uv run python -m grteclyn_wrapper.visualisation.search \
  runs/grtresna_qd/ftl_4d/ftl_4d_v1 \
  --batch-size 8 --rolling 3 \
  --out /tmp/qd_saturation.png
```

Default output: `grteclyn-wrapper/src/grteclyn_wrapper/visualisation/plots/qd_batch_progress_<campaign>.png`.

---

## 8. GW proxies without Weyl4

From extrinsic curvature on a slice (approximate \(h_+\), \(h_\times\) for propagation along \(z\)):

```bash
uv run python -m grteclyn_wrapper.visualisation.visualize --field GW_Plus --axis y --zoom 100 --animate
uv run python -m grteclyn_wrapper.visualisation.visualize --field GW_Cross --axis z --animate --zoom 40
```

---

## 9. Symmetry-reduced domains

For octant-style runs (`lo_boundary = 2 2 2`, etc.), use `--center` and/or `--corner` so the origin is readable on the slice.

```bash
uv run python -m grteclyn_wrapper.visualisation.visualize --field K --axis z --coord 0 --zoom 84 --corner --animate
```

---

## Output layout

```
<run_dir>/
├── data/
│   ├── constraint_norms.dat      # constraines/
│   └── collapse_diagnostics.dat    # diagnostic/
├── small_data/
│   ├── psi4_mode_l2m0.dat          # plot_extracted_psi4 (from consume_plotfiles --psi4)
│   ├── shell_profiles.dat
│   └── areal_radius.dat
├── frames/<field>_<axis>/frames/    # slice PNGs from consume_plotfiles
├── movies/                         # make_movies.sh (MP4 stitch)
└── plots/                          # plot_diagnostic.sh (post-run bundle)

grteclyn-wrapper/src/grteclyn_wrapper/visualisation/
├── visualize/                      # plot_run.sh default frame root (wormhole workflows)
│   ├── <field>_<axis>/frames/
│   └── embedding/frames/
├── plots/                          # QD batch progress; radial diagnostics (plot_diagnostic_radial.sh)
├── extract_wave/                   # default --out for plot_psi4 (plotfile-direct, 2-panel)
├── process_wave/                   # optional default --out for plot_extracted_psi4
├── diagnostic/
├── constraines/                    # default location for relative -o on constraines CLI
├── figures/
├── search/
└── README.md
```

---

## Automation scripts

Shell helpers live under **`grteclyn-wrapper/scripts/plot/`** (and `scripts/wormhole/move_files.sh`).
Run them from the **GRTeclyn repository root** so `uv run --directory grteclyn-wrapper python -m grteclyn_wrapper.visualisation...` resolves correctly.

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
./grteclyn-wrapper/scripts/plot/plot_run.sh

# Explicit run/output folder
./grteclyn-wrapper/scripts/plot/plot_run.sh /path/to/data_2gpu

# Keep existing plotfiles and frames (e.g. resume processing)
./grteclyn-wrapper/scripts/plot/plot_run.sh --not-remove /path/to/data_2gpu

# Fewer parallel workers
./grteclyn-wrapper/scripts/plot/plot_run.sh -j 16 /path/to/data_2gpu

./grteclyn-wrapper/scripts/plot/plot_run.sh --not-remove -j 16
```

**What `plot_run` passes to `consume_plotfiles`** (fixed in the script): extraction radii `12 16 20 24`, `--n-points 64`, `--areal-radius`, `--embedding --embedding-rmax 5.0`, `--frames-fields K Weyl4_Re`, `--frames-axis z`, `--frames-corner`, `--frames-out` → `grteclyn-wrapper/src/grteclyn_wrapper/visualisation/visualize`, `--watch --delete --keep-last 2`, `--verbose`.

**Movies:** This script does **not** call `ffmpeg`. After frames exist under `grteclyn-wrapper/src/grteclyn_wrapper/visualisation/visualize/`, build MP4s with:

```bash
uv run python -m grteclyn_wrapper.visualisation.visualize.make_movies --root grteclyn-wrapper/src/grteclyn_wrapper/visualisation/visualize
```

---

### `plot_diagnostic.sh` — Post-run plots (one run folder)

After the run has written `data/constraint_norms.dat`, `data/collapse_diagnostics.dat`, and
`small_data/psi4_mode_l2m0.dat`, this script **removes and recreates** `<RUN_DIR>/plots/` and writes
fresh figures there.

**Arguments:** optional `RUN_DIR`, then optional extraction **radii** passed through to
`plot_extracted_psi4` (if omitted, all radii columns in the `.dat` file are used).

**Default `RUN_DIR`:** Among `data_2gpu`, `data_supported`, and `data` (sibling of the repo),
chooses the directory whose three required files are present and have the **newest combined
modification time**; if that fails, falls back to an existing `data_2gpu` or `data`.

**Outputs in `<RUN_DIR>/plots/`**

| File stem | Content |
|-----------|---------|
| `constraints_plot.*` | Hamiltonian and momentum L2 constraint norms |
| `collapse_diagnostics_plot.*` | Collapse diagnostics (+ areal radius and K-decay fit when data exist) |
| `psi4_analysis_M1000_D1.*` etc. | 6-panel GW analysis (`--combined --strain`) at several mass/distance configs |

Also attempts K_z and embedding evolution panels from `frames/` (skipped with a warning if frame
indices are missing — common while a run is still in progress).

**Environment:** `MASS_MSUN`, `DISTANCE_MPC` (defaults include `1000:1`, `1000:0.002`, `30:10`);
`ESD_FMAX` (PSD frequency cap); `LIGO_QUANTITY` (`asd` or `hchar`).

```bash
./grteclyn-wrapper/scripts/plot/plot_diagnostic.sh

./grteclyn-wrapper/scripts/plot/plot_diagnostic.sh runs/grtresna_promote/my_hq_run

./grteclyn-wrapper/scripts/plot/plot_diagnostic.sh runs/grtresna_promote/my_hq_run 8 12 24

# Wormhole octant (corner origin)
./grteclyn-wrapper/scripts/plot/plot_diagnostic.sh /path/to/data_supported 12 16 20 24
```

### `make_movies.sh` — Stitch episode frames to MP4

Reads `<EPISODE_DIR>/frames/<field>_<axis>/frames/*.png`, writes
`<EPISODE_DIR>/movies/movie_<field>_<axis>.mp4`. Handles gapped frame numbering (sim step indices).

```bash
./grteclyn-wrapper/scripts/plot/make_movies.sh runs/grtresna_promote/my_hq_run --framerate 10
./grteclyn-wrapper/scripts/plot/make_movies.sh runs/grtresna_promote/my_hq_run --only chi_z Weyl4_Mag_z
```

### `plot_diagnostic_radial.sh` — RadialRecipe (no Psi4)

Constraints + collapse + optional shell profiles → `visualisation/plots/radial/` (not run dir).

```bash
./grteclyn-wrapper/scripts/plot/plot_diagnostic_radial.sh runs/radialrecipe_nonspherical/<episode>
```

---

### `move_files.sh` — Archive wormhole run + visuals to `SimResults/`

Wormhole-specific helper. Copies frames from `visualisation/visualize/`, ASCII from the data
directory, and plot figures into `grteclyn-wrapper/output/SimResults/<auto_folder>/`.

**Plots:** run `plot_diagnostic.sh` with the wormhole data directory as `RUN_DIR` first — figures
land in `<DATA_DIR>/plots/`. `move_files.sh` still copies from `visualisation/plots/` (legacy
path); also copy `<DATA_DIR>/plots/*` into the SimResults folder if needed.

The folder name is derived from wormhole parameters in the params file (e.g.
`Run_R..._A0..._A2..._sigma...`).

**First argument — run type** (optional, default `SupportedWormholeCollapse`):

| Value | Data directory | Params file |
|-------|----------------|-------------|
| `SupportedWormholeCollapse` (default) | `<repo>/../data_supported` | `Examples/SupportedWormholeCollapse/params_2gpu.txt` |
| Any other string (e.g. `WormholeCollapse`) | `<repo>/../data_2gpu` | `Examples/WormholeCollapse/params_2gpu.txt` |

```bash
# Supported wormhole run + data_supported
./grteclyn-wrapper/scripts/wormhole/move_files.sh
./grteclyn-wrapper/scripts/wormhole/move_files.sh SupportedWormholeCollapse

# Other run type (uses data_2gpu and WormholeCollapse params)
./grteclyn-wrapper/scripts/wormhole/move_files.sh WormholeCollapse
```

If `visualisation/plots/` is missing, the script warns to run `plot_diagnostic.sh` first (and
copy `<DATA_DIR>/plots/` if figures were written there).

---

### Suggested workflow

**HQ promotion / search episode** (frames under `RUN_DIR/frames/`):

```bash
# While GPU runs — consumer is usually started by the campaign launcher (CONSUME_PLOTFILES=1)

# After the run (or partial refresh while consumer catches up)
./grteclyn-wrapper/scripts/plot/plot_diagnostic.sh runs/grtresna_promote/my_hq_run 8 12 24

./grteclyn-wrapper/scripts/plot/make_movies.sh runs/grtresna_promote/my_hq_run --framerate 10
```

**Wormhole standalone run** (frames under `visualisation/visualize/` via `plot_run.sh`):

```bash
# Terminal 1 — while the simulation runs
./grteclyn-wrapper/scripts/plot/plot_run.sh /path/to/data_supported

# After the run — diagnostic figures into the data directory's plots/ (if RUN_DIR=data_supported)
./grteclyn-wrapper/scripts/plot/plot_diagnostic.sh /path/to/data_supported 12 16 20 24

uv run python -m grteclyn_wrapper.visualisation.visualize.make_movies \
  --root grteclyn-wrapper/src/grteclyn_wrapper/visualisation/visualize

# Archive wormhole SimResults bundle
./grteclyn-wrapper/scripts/wormhole/move_files.sh SupportedWormholeCollapse
```

**Extraction center:** wormhole runs use `--frames-corner` and Psi4 radii from `(0,0,0)`; full-box
Q-ball / RadialRecipe HQ runs use physics `center` (e.g. `64 64 64` for L=128) — set via
`CONSUMER_RADII` / `--consumer-radii` on promotion replay, not via the plotting scripts themselves.

---

## Custom paths (examples)

```bash
uv run python -m grteclyn_wrapper.visualisation.visualize --field K --data /path/to/plt --out /path/to/output
```

```bash
uv run python -m grteclyn_wrapper.visualisation.constraines /path/to/data/data/constraint_norms.dat
```

```bash
uv run python -m grteclyn_wrapper.visualisation.visualize.make_movies --root /path/to/visualize
```

Archive runs: see **`move_files.sh`** in [Automation scripts](#automation-scripts) above.
