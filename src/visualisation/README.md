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

### 4. Constraint norms (`constraines/`)

Some simulations write constraint norms to a `constraint_norms.dat` file (SmallDataIO format), typically under a `data/` subfolder next to your plotfiles.

The constraints plotting script accepts an explicit path, so if your simulation output folder changes you do **not** need to edit any code — just pass the new file location.

**Usage:**

```bash
# Default-ish case (simulation output in ../data relative to repo root):
python -m src.visualisation.constraines ../data/data/constraint_norms.dat

# Example: WormholeCollapse multi-GPU run output in data_2gpu/
python -m src.visualisation.constraines /home/jovyan/nachevsky/test/simulation/data_2gpu/data/constraint_norms.dat
```

### 5. Visualizing Symmetry-Reduced Domains

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
├── visualize/           # Field slice plots and movies
│   ├── chi_z/          # From --field chi --axis z
│   │   ├── frames/
│   │   └── movie_chi_z.mp4
│   ├── chi_x/
│   ├── chi_y/
│   ├── K_z/
│   ├── lapse_z/
│   ├── Theta_z/
│   └── ...
├── extract_wave/        # GW waveform extraction
│   └── psi4_analysis_R14.0_R30.0_n24.png
└── README.md
```
---

## Custom paths

To override data or output paths:

```bash
python -m src.visualisation.visualize --field K --data /path/to/plt --out /path/to/output
```

