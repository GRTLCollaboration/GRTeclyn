# GRTeclyn Visualization

Scripts for visualizing GRTeclyn simulation results (plotfiles in BoxLib/AMReX format).

- **`visualize/`** — Field slice plots and movies. Outputs go to `src/visualisation/visualize/`.
- **`extract_wave/`** — GW extraction from Weyl scalar \(\Psi_4\) (Weyl4) and spectra. Outputs go to `src/visualisation/extract_wave/`.

## Prerequisites

1. Install dependencies (from the project root):

   ```bash
   uv sync
   # or: source .venv/bin/activate && pip install -e .
   ```

2. Required packages: `yt`, `numpy`, `matplotlib`, `scipy` (from `pyproject.toml`).
   Optional for MPI: `uv sync --extra mpi` or `pip install mpi4py`

3. For animations: `ffmpeg` must be installed.

4. Simulation data: plotfiles (`plt*` folders) in the default location `../data` relative to the project root.
   If your simulation writes somewhere else (e.g. `output_path = ".../data_2gpu"` in `params_2gpu.txt`), pass `--data` to point to that folder.

---

## Scripts

### 1. `visualize/` — Main field visualization

Plots slice views of GR fields. You can choose the slice plane (axis and coordinate). Supports multiple fields and can create MP4 animations. Outputs to `src/visualisation/visualize/` by default.

**Available fields:**

| Field      | Description                          |
|------------|--------------------------------------|
| `chi`      | Conformal factor                     |
| `K`        | Trace of extrinsic curvature         |
| `Theta`    | Z4 constraint                        |
| `lapse`    | Lapse (alpha)                        |
| `Ham`      | Hamiltonian constraint               |
| `A11`, `A12`, `A22`, `A33` | Extrinsic curvature components |
| `GW_Plus`  | GW strain proxy (+): A11 − A22       |
| `GW_Cross` | GW strain proxy (×): 2 A12            |

**Understanding the variables**

Numerical relativity (like GRTeclyn) treats the universe as a "stack" of 3D snapshots. These variables describe the geometry and the "flow of time" for those snapshots.

1. **Chi (χ) — The "dents" in space**
   - *What it is:* The Conformal Factor.
   - *Think of it as:* The "depth" of the gravitational well.
   - In vacuum: χ = 1.0. Space is flat.
   - Near a mass: χ drops below 1.0.
   - At the center of a black hole: χ → 0.
   - *Best for:* Finding where the black holes are. Dark "pits" in the plot show low χ.

2. **Lapse (α) — The speed of time**
   - *What it is:* The Lapse Function.
   - *Think of it as:* A "local clock speed."
   - Far away: Lapse ≈ 1.0. Time flows normally.
   - Near a black hole: Lapse drops toward 0.
   - *Why it matters:* Time slows down near massive objects. The "moving puncture" trick lets the lapse go to near-zero at the singularity so the simulation stays stable.
   - *Best for:* Seeing where time is "slowing down." Looks similar to Chi.

3. **K — The "breathing" of space**
   - *What it is:* The Trace of Extrinsic Curvature.
   - *Think of it as:* How much the 3D slice of space is "curving" into the 4th dimension (Time).
   - *Physics:* Represents the local expansion or contraction of space.
   - *Best for:* Seeing gravitational waves and the "wake" left by moving black holes. Chi shows the hole itself; K shows the "splashes."

4. **Theta (Θ) — The "quality check"**
   - *What it is:* The Z4 Constraint variable.
   - *Think of it as:* Mathematical "noise" or error tracking.
   - *Physics:* In a perfect universe, Θ should be exactly 0. Computers make rounding errors; Einstein's equations are sensitive. The CCZ4 formulation "pushes" errors away from the holes toward the box edges.
   - *Best for:* Checking simulation health. Large Θ (bright colors) = possible instability. Like a "check engine light."

**Summary table**

| Variable    | Normal value | BH value   | Best used for                                |
|-------------|--------------|------------|----------------------------------------------|
| Chi (χ)     | 1.0          | → 0.0      | Finding the center/position of the holes.    |
| Lapse       | 1.0          | → 0.0      | Seeing where time is "slowing down."         |
| K           | 0.0          | High +/-   | Seeing the ripples and motion of space.      |
| Theta (Θ)   | 0.0          | Small +/-  | Checking if the simulation is healthy/accurate. |

**Basic usage:**

```bash
# Single-frame plot of the last timestep (default: chi)
python -m src.visualisation.visualize --field K

# Animation for a specific field
python -m src.visualisation.visualize --field K --animate --zoom 100

# If your simulation output folder is not the default ../data (e.g. output_path=.../data_2gpu),
# pass it explicitly:
python -m src.visualisation.visualize --field K --animate --zoom 100 --data /home/jovyan/nachevsky/test/simulation/data_2gpu
```

**Options:**

| Option       | Default  | Description                              |
|--------------|----------|------------------------------------------|
| `--field`    | `chi`    | Field to plot                            |
| `--axis`     | `z`      | Axis normal to slice: `x`, `y`, or `z`   |
| `--coord`    | see below| Coordinate along axis (default: 0 for z, 64 for x/y) |
| `--animate`  | off      | Create MP4 animation from all plotfiles  |
| `--zoom`     | full domain | Plot width in code units            |
| `--data`     | `../data`| Path to plt folders                      |
| `--out`      | `visualize/` | Output directory                         |
| `--mpi`      | off      | Use MPI for parallel frame generation   |

**Zoom explanation** (for a box size L_full = 128):

| `--zoom` | What you see |
|----------|--------------|
| `128` | The entire domain—from the left wall to the right wall. |
| `64` | Half of the domain. Black holes appear twice as large. |
| `10` | A small 10×10 area. May show only one black hole (BH separation is ~12–20 units). |

**Choosing the slice plane:**

1. **Standard top-down view (default)** — Looks at the orbital plane (z=0).

   ```bash
   python -m src.visualisation.visualize --field chi --animate
   # or explicitly:
   python -m src.visualisation.visualize --field chi --axis z --animate
   ```

2. **Side view** — Cutting through the Y-axis. Shows "gravity wells" from the side. Center of the 128-unit box is at 64.

   ```bash
   python -m src.visualisation.visualize --field chi --axis y --animate
   # or with explicit coordinate:
   python -m src.visualisation.visualize --field chi --axis y --coord 64 --animate
   ```

3. **Front view** — Cutting through the X-axis.

   ```bash
   python -m src.visualisation.visualize --field chi --axis x --animate
   ```

4. **Visualizing waves (side view)** — Gravitational waves propagating upwards look best from the side.

   ```bash
   python -m src.visualisation.visualize --field GW_Plus --axis y --zoom 100 --animate
   # or use K for extrinsic curvature ripples
   python -m src.visualisation.visualize --field K --axis y --zoom 100 --animate
   ```

**MPI (parallel) usage:**

```bash
mpirun -np 8 python -m src.visualisation.visualize --field K --animate --zoom 100 --mpi
```

**Output:** `{field}_{axis}/frames/frame_XXXX.png` and `{field}_{axis}/movie_{field}_{axis}.mp4` (if `--animate`).

### 1b. Evolution panel (`make_evolution_panel.py`)

Horizontal strip figure from existing PNG frames (one column per timestep). Looks for `frame_z_####.png` or `frame_####.png` inside `--frame_dir`. Writes both `.png` and `.pdf` next to `--out` (extension is replaced).

**Default frames** if you omit `--frames`: `0 20 40 60`.

**With explicit frame indices** (any count; indices are the `####` in the filename):

```bash
# From the repository root; example: K on the z slice
python src/visualisation/make_evolution_panel.py \
  --frame_dir src/visualisation/visualize/K_z/frames \
  --out src/visualisation/evolution_K_z_panel \
  --frames 1600 1610 1620 1630

# Keep original colours (default is grayscale + white background cleanup)
python src/visualisation/make_evolution_panel.py \
  --frame_dir src/visualisation/visualize/K_z/frames \
  --out src/visualisation/evolution_K_z_panel_colour \
  --frames 0 500 1000 1500 \
  --no-grayscale
```

| Option | Description |
|--------|-------------|
| `--frame_dir` | Directory containing the frame PNGs (e.g. `visualize/K_z/frames`) |
| `--out` | Output base path (`.png` / `.pdf` are written) |
| `--frames` | Space-separated indices (default `0 20 40 60`) |
| `--no-grayscale` | Do not convert panels to grayscale |

---

### 2. Gravitational Wave Visualization (without Weyl4)

In the absence of full Weyl scalar ($\Psi_4$) extraction, we can visualize gravitational waves using proxies derived from the **Extrinsic Curvature ($A_{ij}$)**.

Since gravitational waves are perturbations in the spacetime metric, they manifest as ripples in the extrinsic curvature components. For a wave propagating along the **z-axis**, the transverse-traceless (TT) gauge strain components $h_+$ and $h_\times$ can be approximated by:

- **$h_+$ (Plus Polarization):** Proportional to $A_{11} - A_{22}$ (or $A_{xx} - A_{yy}$).
- **$h_\times$ (Cross Polarization):** Proportional to $2 A_{12}$ (or $2 A_{xy}$).

The `visualize` module provides derived fields `GW_Plus` and `GW_Cross` to plot these directly.

**Usage:**

```bash
# Visualize the Plus polarization
python -m src.visualisation.visualize --field GW_Plus --axis z --animate --zoom 40

# Visualize the Cross polarization
python -m src.visualisation.visualize --field GW_Cross --axis z --animate --zoom 40
```

### 3. Wave Extraction Tool (`extract_wave/`)

This module extracts the gravitational wave signal from the Weyl scalar \(\Psi_4\) (Weyl4) at multiple radii and plots \(\Psi_4(t)\) and its spectrum (PSD).

We intentionally **do not reconstruct strain \(h(t)\)** in this workflow, because double time-integration is ill-conditioned and tends to amplify low-frequency gauge/near-zone contamination. Instead, we infer the GW frequency content directly from the \(\Psi_4\) modes.

The extractor supports **symmetry-reduced domains** (e.g. octant symmetry with `lo_boundary = 2 2 2`) by sampling only the in-domain part of the extraction sphere and renormalizing the surface integral to \(4\pi\).

**Usage:**

```bash
# Plot Psi4 (time series + PSD). Data defaults to ../data relative to project root.
python -m src.visualisation.extract_wave --radii 14 30 --n-points 24

# Plot with retarded time (t - R):
python -m src.visualisation.extract_wave --radii 14 30 --n-points 24 --time-axis retarded

# With custom data path:
python -m src.visualisation.extract_wave --data /path/to/data --radii 14 30 --n-points 24

# Example: WormholeCollapse multi-GPU run that writes to data_2gpu/
python -m src.visualisation.extract_wave --data /home/jovyan/nachevsky/test/simulation/data_2gpu --radii 14 30 --n-points 24
```

**Options (both scripts):**
- `--data PATH`: Path to simulation data (default: `../data` from project root)
- `--out PATH`: Output directory (default: the script folder)
- `--radii R1 R2 ...`: Radii to extract from (e.g. `14 30`)
- `--n-points N`: Angular resolution (samples an \(N\\times N\) grid on the sphere)
- `--time-axis {simulation,retarded}`: Plot x-axis as simulation time \(t\) or retarded time \(t-R\)
- `--t-min`, `--t-max`: x-range to display
- `--psd-hide-raw`: Hide the raw (dotted) PSD points and show only the smoothed curve

**Output:**
- `psi4_analysis_R14.0_R30.0_n24.png` (time series + PSD)

**Note on “negative time” in plots**

If you use retarded time \(t_{\\rm ret} = t - R_{\\rm ext}\), the x-axis can be negative: this corresponds to early simulation times (\(t < R_{\\rm ext}\)) before the wave reaches that radius. Any structure there is usually **junk radiation / gauge transient** or near-zone curvature.

### 3b. Strain PSD & Propagation Speed (`process_wave/plot_extracted_psi4.py`)

The `plot_extracted_psi4.py` script can convert the \(\Psi_4\) PSD into a **strain PSD** and overlay it on the Advanced LIGO design sensitivity curve.

It can also measure the **propagation speed** of the signal across extraction radii to separate physical gravitational waves (\(v = c\)) from CCZ4 constraint-cleaning modes (\(v \neq c\)).

**Strain conversion math:** \(|\tilde{h}(f)| = |\tilde{\Psi}_4(f)| / (2\pi f)^2\), so the strain PSD is \(S_h(f) = S_{\Psi_4}(f) / (2\pi f)^4\). Physical frequency scaling uses \(f_{\text{phys}} = f_{\text{code}} / (M \cdot G M_\odot / c^3)\).

**Usage:**

```bash
# Basic waveform + PSD (existing)
python -m src.visualisation.process_wave.plot_extracted_psi4 \
  /path/to/small_data/psi4_mode_l2m0.dat --plot-psd

# Strain PSD + LIGO overlay + SNR (for 30 M_sun at 10 Mpc)
python -m src.visualisation.process_wave.plot_extracted_psi4 \
  /path/to/small_data/psi4_mode_l2m0.dat \
  --plot-psd --strain --mass-msun 30 --distance-mpc 10

# Propagation speed analysis (GW vs junk separation)
python -m src.visualisation.process_wave.plot_extracted_psi4 \
  /path/to/small_data/psi4_mode_l2m0.dat \
  --propagation-speed

# All analyses combined
python -m src.visualisation.process_wave.plot_extracted_psi4 \
  /path/to/small_data/psi4_mode_l2m0.dat \
  --plot-psd --strain --propagation-speed --mass-msun 30 --distance-mpc 10
```

**Strain options:**

| Option             | Default | Description                                       |
|--------------------|---------|---------------------------------------------------|
| `--strain`         | off     | Enable strain PSD conversion + LIGO overlay        |
| `--mass-msun`      | `30.0`  | Total mass in solar masses for physical scaling    |
| `--distance-mpc`   | `10.0`  | Luminosity distance in Mpc for SNR calculation     |
| `--propagation-speed` | off  | Measure signal speed across extraction radii       |

**Output panels** (when all flags active):

1. Time-domain waveform \(r \cdot \text{Re}(\Psi_4^{2,0})\)
2. \(\Psi_4\) PSD in code units
3. Strain PSD \(S_h(f)\) in code units
4. Characteristic strain \(h_{\text{char}}(f) = \sqrt{f \cdot S_h}\) vs Advanced LIGO \(\sqrt{f \cdot S_n}\)
5. Propagation speed panel with peak markers and \(v/c\) annotations

**Console output** includes peak frequency (code + Hz), SNR, frequency band classification, and per-radius-pair propagation speeds with GW/constraint-mode classification.

### 3c. Areal Radius & Embedding Extraction (`process_wave/consume_plotfiles.py`)

The plotfile consumer can also extract **areal radius** and render **embedding diagram frames** from each plotfile during live processing.

#### Areal radius

At each timestep the script shoots a `yt` ray along the positive x-axis (at y=0, z=0), reads the conformal factor \(\chi\), and computes the areal radius \(R_{\text{areal}}(r) = r / \sqrt{\chi(r)}\). The minimum of this profile (excluding the coordinate singularity at r=0) is the physical throat size. It is appended to `areal_radius.dat` as `(t, R_{\text{areal,min}}, r_{\text{at\,min}})`.

#### Embedding diagrams — algorithm

The embedding diagram visualises the spatial geometry around the wormhole throat as a 3D surface of revolution. The construction works as follows:

1. **Extract \(\chi(r)\) along x-axis.** A `yt` ray from the domain centre to the right boundary yields cell-centred values of \(\chi\) at the AMR cell positions.

2. **Interpolate onto a uniform grid.** The raw ray data has uneven spacing (fine near the throat, coarse far away) and step-like jumps at AMR refinement boundaries. The values are linearly interpolated onto a uniform grid of 1024+ points, then smoothed with a Gaussian filter whose width matches the coarsest AMR cell size. This eliminates derivative discontinuities at level boundaries while preserving physical structure.

3. **Compute \(R_{\text{areal}}(r)\) and \(z_{\text{embed}}(r)\).** In conformal-flat isotropic coordinates the spatial metric is \(dl^2 = \chi^{-1} dr^2\). The isometric embedding into flat 3D space requires a surface of revolution \((R(r),\, z(r))\) satisfying \(dR^2 + dz^2 = dl^2\). With \(R_{\text{areal}} = r/\sqrt{\chi}\):

   \[
   \frac{dz}{dr} = \sqrt{\frac{1}{\chi} - \left(\frac{dR_{\text{areal}}}{dr}\right)^{\!2}}
   \]

   \(z_{\text{embed}}\) is obtained by cumulative trapezoidal integration of this expression.

4. **Mirror to show both sheets.** Octant symmetry places the throat at x=0. The one-sided profile (x > 0) is reflected to negative \(z_{\text{embed}}\), producing the classic symmetric two-funnel wormhole shape with the throat at \(z=0\).

5. **Surface of revolution.** The 1D profile \((R_{\text{areal}},\, z_{\text{embed}})\) is swept through \(\phi \in [0, 2\pi]\) to produce a 3D surface, coloured by \(z_{\text{embed}}\) (viridis colourmap).

6. **Fixed axis limits.** The x, y, and z axes are locked to values derived from `--embedding-rmax` so that the coordinate box stays static across all frames. Only the wormhole geometry changes between frames, producing smooth animations without axis jumps.

**Usage:**

```bash
# Live extraction with areal radius + embedding (during simulation)
python -m src.visualisation.process_wave.consume_plotfiles \
  --data /path/to/data_2gpu \
  --radii 8 12 16 --n-points 32 \
  --areal-radius \
  --embedding --embedding-rmax 5.0 \
  --frames-fields chi K \
  --frames-corner \
  --watch --delete --keep-last 2 --verbose -j 64
```

| Option | Description |
|---|---|
| `--areal-radius` | Extract min areal radius \(R_{\text{areal}} = r/\sqrt{\chi}\) along x-axis to `areal_radius.dat` |
| `--embedding` | Render 3D embedding diagram frames (surface of revolution from \(\chi\) profile) |
| `--embedding-rmax` | Max coordinate radius for the embedding (default: full domain). Choose ~10× the throat radius for a good view of the funnel without the flat region dominating. |
| `-j`, `--jobs` | Run plotfile processing and rendering in parallel across N workers (e.g. `-j 64`). Drastically reduces rendering time. |
| `--keep-existing-frames` | Do not clear pre-existing visualization frames from output folders upon startup (useful when resuming processing). |

**Choosing `--embedding-rmax`:** If the throat radius is \(b\), a good starting value is `--embedding-rmax` \(\approx 10b\). Too large a value makes the flat asymptotic region dominate and the throat becomes invisible; too small clips the funnel before it flares out.

**Output files:**
- `small_data/areal_radius.dat` — time-series of throat areal radius
- `visualize/embedding/frames/frame_NNNN.png` — 3D embedding diagram frames (can be stitched with `make_movies.py`)

### 4. Collapse diagnostics (`diagnostic/`)

Plots streaming collapse diagnostics: min lapse, min chi, max|K|, trapped-surface proxy, null expansion. When `areal_radius.dat` is available (from `consume_plotfiles.py --areal-radius`), also plots **throat areal radius**, **expansion velocity** \(dR_{\text{areal}}/dt\), and fits the **K-decay lifetime** \(\tau\).

**Usage:**

```bash
# Standard diagnostics (auto-detects areal_radius.dat from --data)
python3 src/visualisation/diagnostic/diagnostic.py \
  --data /path/to/data_2gpu

# With explicit areal radius file
python3 src/visualisation/diagnostic/diagnostic.py \
  /path/to/collapse_diagnostics.dat \
  --data /path/to/data_2gpu \
  --areal-radius-file /path/to/areal_radius.dat

# Disable K-decay lifetime fit
python3 src/visualisation/diagnostic/diagnostic.py \
  --data /path/to/data_2gpu --no-fit-lifetime

# Custom mass for physical unit conversion
python3 src/visualisation/diagnostic/diagnostic.py \
  --data /path/to/data_2gpu --mass-msun 30
```

| Option                | Default | Description                                       |
|-----------------------|---------|---------------------------------------------------|
| `--data`              | auto    | Run output directory                               |
| `--areal-radius-file` | auto   | Path to `areal_radius.dat` (auto-detected)         |
| `--no-fit-lifetime`   | off     | Disable exponential fit to K decay                 |
| `--mass-msun`         | `30.0`  | Total mass for physical unit conversion            |

**Plot layout** (3x3 when areal radius data available):
- Row 1: min(lapse), min(chi), max|K|
- Row 2: trapped-surface radius, min(theta+), radius at min(theta+)
- Row 3: R\_areal\_min(t), expansion velocity dR/dt (with c=1 line), K-decay with exponential fit

**Console output** includes expansion velocity (mean/peak in units of c), K-decay lifetime \(\tau\) in code units and physical seconds.

### 5. Constraint norms (`constraines/`)

Some simulations write constraint norms to a `constraint_norms.dat` file (SmallDataIO format), typically under a `data/` subfolder next to your plotfiles.

The constraints plotting script accepts an explicit path, so if your simulation output folder changes you do **not** need to edit any code — just pass the new file location.

**Usage:**

```bash
# Default-ish case (simulation output in ../data relative to repo root):
python -m src.visualisation.constraines ../data/data/constraint_norms.dat

# Example: WormholeCollapse multi-GPU run output in data_2gpu/
python -m src.visualisation.constraines /home/jovyan/nachevsky/test/simulation/data_2gpu/data/constraint_norms.dat
```

### 6. Visualizing Symmetry-Reduced Domains

When a simulation uses symmetry (e.g. octant symmetry with `lo_boundary = 2 2 2`), only a portion of the physical space is simulated (e.g. \(x\\ge 0, y\\ge 0, z\\ge 0\)). The physical “center” is at the symmetry planes (the origin).

The `visualize` tool supports:
- `--center x y z`: explicit plot center
- `--corner`: place the symmetry origin at the bottom-left corner of the 2D plot (helpful for octant/quadrant runs)

**Example (octant run, domain 0→84):**

```bash
python -m src.visualisation.visualize --field K --axis z --coord 0 --zoom 84 --corner --animate
```

- With `--corner`, the symmetry origin (0,0) appears at the bottom-left of the image, so outward-propagating waves are easier to interpret visually.

---

## Output layout

```
src/visualisation/
├── visualize/                  # Field slice plots and movies
│   ├── chi_z/                  # From --field chi --axis z
│   │   ├── frames/
│   │   └── movie_chi_z.mp4
│   ├── K_z/
│   ├── lapse_z/
│   ├── Theta_z/
│   ├── embedding/              # 3D embedding diagram frames (from --embedding)
│   │   └── frames/
│   │       └── frame_NNNN.png
│   └── ...
├── extract_wave/               # GW waveform extraction
│   └── psi4_analysis_R14.0_R30.0_n24.png
├── process_wave/               # Processed waveform plots
│   ├── psi4_extracted_R*.png   # Time-domain + PSD
│   ├── psi4_strain_analysis.*  # Strain PSD + LIGO overlay + SNR
│   └── psi4_propagation_speed.*# Propagation speed analysis
├── diagnostic/                 # Collapse diagnostics
│   └── collapse_diagnostics_plot.png  # Multi-panel (up to 3x3)
└── README.md

<data_dir>/small_data/
├── psi4_mode_l2m0.dat          # Extracted Psi4 time series
└── areal_radius.dat            # Throat areal radius time series
```
---

## Automation Scripts

Two shell scripts automate the full visualization pipeline. Both auto-detect the most recent run directory.

### `plot_run.sh` — Live plotfile processing (run during simulation)

Watches for new plotfiles, extracts Ψ4 waveforms, areal radius, renders field frames and embedding diagrams. Deletes consumed plotfiles to save disk space. It automatically processes plotfiles in parallel for maximum speed.

```bash
# Default (auto-detects data_2gpu or data, uses 64 jobs):
./src/scripts/plot_run.sh

# Explicit data path:
./src/scripts/plot_run.sh /path/to/data_2gpu

# Keep existing plotfiles and generated frames (useful for restarting processing):
./src/scripts/plot_run.sh --not-remove /path/to/data_2gpu

# Specify a custom number of parallel jobs (e.g. 16 workers instead of 64):
./src/scripts/plot_run.sh -j 16 /path/to/data_2gpu

# Run with a custom number of jobs AND keep existing files (order of arguments does not matter):
./src/scripts/plot_run.sh --not-remove -j 16
./src/scripts/plot_run.sh --not-remove -j 16 /path/to/data_2gpu
```

What it does:
1. Removes stale plotfiles, old extracted small-data, and old image frames (unless `--not-remove` is provided, which preserves them).
2. Launches `consume_plotfiles.py` in watch mode with:
   - `--radii 8 12 16 --n-points 32`
   - `--areal-radius` (throat areal radius extraction)
   - `--embedding --embedding-rmax 5.0` (3D embedding diagram frames)
   - `--frames-fields chi K Weyl4_Re Weyl4_Mag`
   - `--watch --delete --keep-last 2`
   - `-j 64` (uses 64 parallel worker processes by default, adjustable via `-j` / `--jobs`)
   - `--keep-existing-frames` (if `--not-remove` was specified)

### `plot_diagnostic.sh` — Post-run diagnostics (run after simulation)

Generates all diagnostic plots from the extracted data files. All output goes into a single folder `src/visualisation/plots/` which is **wiped and recreated on each run** so results are always fresh.

```bash
# Default (auto-detects most recent run):
./src/scripts/plot_diagnostic.sh

# Explicit run directory:
./src/scripts/plot_diagnostic.sh /path/to/data_2gpu

# With specific extraction radii:
./src/scripts/plot_diagnostic.sh /path/to/data_2gpu 8 12 16
```

What it does (6 steps):
1. Constraint norms plot
2. Collapse diagnostics (+ areal radius + K-decay lifetime if `areal_radius.dat` exists)
3. Ψ4 waveform + PSD in retarded time
4. Ψ4 waveform + PSD in simulation time
5. Strain PSD + LIGO noise overlay + SNR (30 M⊙ at 10 Mpc)
6. Propagation speed analysis (GW vs constraint-mode separation)

**Output files** (all in `src/visualisation/plots/`):

| File                            | Content                          |
|---------------------------------|----------------------------------|
| `constraints_plot.*`            | Constraint norms (Ham + Mom)     |
| `collapse_diagnostics_plot.*`   | Collapse diagnostics (up to 3x3)|
| `psi4_retarded.*`               | Ψ4 waveform + PSD (retarded)    |
| `psi4_simulation.*`             | Ψ4 waveform + PSD (simulation)  |
| `psi4_strain_analysis.*`        | Strain PSD + LIGO + SNR          |
| `psi4_propagation_speed.*`      | Propagation speed analysis       |

Each plot is saved as `.png`, `.eps`, and `.pdf`.

### Typical workflow

```bash
# 1. Start simulation + live processing in one terminal:
./src/scripts/plot_run.sh /path/to/data_2gpu

# 2. After simulation completes, generate all diagnostic plots:
./src/scripts/plot_diagnostic.sh /path/to/data_2gpu

# 3. Copy results to SimResults archive:
./src/scripts/move_files.sh SupportedWormholeCollapse
```

---

## Custom paths

To override data or output paths:

```bash
python -m src.visualisation.visualize --field K --data /path/to/plt --out /path/to/output
```

```bash
python -m src.visualisation.constraines /home/jovyan/nachevsky/test/simulation/data_2gpu/data/constraint_norms.dat
```

```bash
python3 /home/jovyan/nachevsky/test/simulation/GRTeclyn/src/visualisation/visualize/make_movies.py \
  --root /home/jovyan/nachevsky/test/simulation/GRTeclyn/src/visualisation/visualize
```

Copy the results to a run-labelled folder under `src/SimResults/` (folder name is derived from `params_2gpu.txt`; optional first arg is run type):

```bash
./src/scripts/move_files.sh
./src/scripts/move_files.sh SupportedWormholeCollapse
```