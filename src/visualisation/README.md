# GRTeclyn Visualization

Scripts for visualizing GRTeclyn simulation results (plotfiles in BoxLib/AMReX format). All outputs are saved in this folder.

## Prerequisites

1. Activate the virtual environment (from the project root):

   ```bash
   source .venv/bin/activate
   ```

2. Required Python packages: `yt`, `numpy`, `matplotlib` (installed via `pyproject.toml`).  
   Optional for MPI: `pip install mpi4py`

3. For animations: `ffmpeg` must be installed.

4. Simulation data: plotfiles (`plt*` folders) in the default location `../data` relative to the project root, or pass `--data` to override.

---

## Scripts

### 1. `visualize.py` — Main field visualization

Plots slice views of GR fields. You can choose the slice plane (axis and coordinate). Supports multiple fields and can create MP4 animations.

**Available fields:**

| Field      | Description                          |
|------------|--------------------------------------|
| `chi`      | Conformal factor                     |
| `K`        | Trace of extrinsic curvature         |
| `Theta`    | Z4 constraint                        |
| `lapse`    | Lapse (alpha)                        |
| `Weyl4_Re` | Gravitational waves (Re Psi4)        |

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
python src/visualisation/visualize.py --field K

# Animation for a specific field
python src/visualisation/visualize.py --field K --animate --zoom 100
```

**Options:**

| Option       | Default  | Description                              |
|--------------|----------|------------------------------------------|
| `--field`    | `chi`    | Field to plot                            |
| `--axis`     | `z`      | Axis normal to slice: `x`, `y`, or `z`   |
| `--coord`    | see below| Coordinate along axis (default: 0 for z, 64 for x/y) |
| `--animate`  | off      | Create MP4 animation from all plotfiles  |
| `--zoom`     | `60`     | Plot width in code units                 |
| `--data`     | `../data`| Path to plt folders                      |
| `--out`      | this dir | Output directory                         |
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
   python src/visualisation/visualize.py --field chi --animate
   # or explicitly:
   python src/visualisation/visualize.py --field chi --axis z --animate
   ```

2. **Side view** — Cutting through the Y-axis. Shows "gravity wells" from the side. Center of the 128-unit box is at 64.

   ```bash
   python src/visualisation/visualize.py --field chi --axis y --animate
   # or with explicit coordinate:
   python src/visualisation/visualize.py --field chi --axis y --coord 64 --animate
   ```

3. **Front view** — Cutting through the X-axis.

   ```bash
   python src/visualisation/visualize.py --field chi --axis x --animate
   ```

4. **Visualizing waves (side view)** — Gravitational waves propagating upwards look best from the side.

   ```bash
   python src/visualisation/visualize.py --field Weyl4_Re --axis y --zoom 100 --animate
   # or use K if Weyl4_Re is not available
   python src/visualisation/visualize.py --field K --axis y --zoom 100 --animate
   ```

**MPI (parallel) usage:**

```bash
mpirun -np 8 python src/visualisation/visualize.py --field K --animate --zoom 100 --mpi
```

**Output:** `{field}_{axis}/frames/frame_XXXX.png` and `{field}_{axis}/movie_{field}_{axis}.mp4` (if `--animate`).

---

### 2. `visualize_sim.py` — Binary black hole chi movie

Creates a binary BH movie using the conformal factor `chi`, tuned for puncture/BH evolution.

```bash
python src/visualisation/visualize_sim.py
```

**Output:** `binary_bh/frames/frame_XXXX.png` and `binary_bh/binary_bh_movie.mp4`.

---

### 3. `visualize_waves.py` — Gravitational wave plots

Plots the `Weyl4_Re` field (real part of Weyl4) to visualize gravitational waves. Skips `plt00000` (usually no wave data at t=0).

```bash
python src/visualisation/visualize_waves.py
```

**Output:** `waves/wave_frame_XXXX.png`.

---

## Output layout

```
src/visualisation/
├── chi_z/          # From visualize.py --field chi --axis z
│   ├── frames/
│   └── movie_chi_z.mp4
├── chi_x/
├── chi_y/
├── K_z/
├── lapse/
├── Theta/
├── binary_bh/      # From visualize_sim.py
│   ├── frames/
│   └── binary_bh_movie.mp4
└── waves/          # From visualize_waves.py
    └── wave_frame_XXXX.png
```

---

## Custom paths

To override data or output paths:

```bash
python src/visualisation/visualize.py --field K --data /path/to/plt --out /path/to/output
```
