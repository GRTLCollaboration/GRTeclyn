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

Plots slice views of GR fields at the z=0 mirror plane. Supports multiple fields and can create MP4 animations.

**Available fields:**

| Field   | Description                          |
|---------|--------------------------------------|
| `chi`   | Conformal factor                     |
| `K`     | Trace of extrinsic curvature         |
| `Theta` | Z4 constraint                        |
| `lapse` | Lapse (alpha)                        |

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
| `--animate`  | off      | Create MP4 animation from all plotfiles  |
| `--zoom`     | `60`     | Plot width in code units                 |
| `--data`     | `../data`| Path to plt folders                      |
| `--out`      | this dir | Output directory                         |
| `--mpi`      | off      | Use MPI for parallel frame generation   |

**MPI (parallel) usage:**

```bash
mpirun -np 8 python src/visualisation/visualize.py --field K --animate --zoom 100 --mpi
```

**Output:** `{field}/frames/frame_XXXX.png` and `{field}/movie_{field}.mp4` (if `--animate`).

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
├── chi/           # From visualize.py --field chi
│   ├── frames/
│   └── movie_chi.mp4
├── K/
├── lapse/
├── Theta/
├── binary_bh/     # From visualize_sim.py
│   ├── frames/
│   └── binary_bh_movie.mp4
└── waves/         # From visualize_waves.py
    └── wave_frame_XXXX.png
```

---

## Custom paths

To override data or output paths:

```bash
python src/visualisation/visualize.py --field K --data /path/to/plt --out /path/to/output
```
