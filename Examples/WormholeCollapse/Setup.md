# WormholeCollapse Setup

## Prerequisites

**CUDA must be in your PATH.** The build uses `nvcc` directly (not just `CUDA_HOME`).
If you see `nvcc: not found`, it means CUDA is installed but not activated in your shell
(or CUDA is provided via environment modules).

```bash
# If your system uses modules (common on clusters), do this first:
# module avail cuda
# module load cuda

# Otherwise, add CUDA to PATH (adjust to your local install):
export PATH=/usr/local/cuda/bin:$PATH

# Common alternatives you may need:
# export PATH=/usr/local/cuda-12/bin:$PATH
# export PATH=/usr/local/cuda-12.*/bin:$PATH
# export PATH=/opt/cuda/bin:$PATH

command -v nvcc && nvcc --version
```

Add to `~/.bashrc` if needed:
```bash
echo 'export PATH=/usr/local/cuda/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

## MPI Prerequisite For Multi-GPU Runs

The MPI executable `main3d.gnu.MPI.CUDA.ex` needs both:

- `mpirun` in your `PATH`
- `libmpi.so.40` in your library path

If you see either

```bash
bash: mpirun: command not found
```

or

```bash
ldd ./main3d.gnu.MPI.CUDA.ex | grep libmpi
libmpi.so.40 => not found
```

then CUDA is active but the OpenMPI runtime is not.

### Local OpenMPI under `~/nachevsky`

If you installed OpenMPI locally at
`/home/jovyan/nachevsky/test/simulation/local/openmpi-5.0.8`, activate it in
the current shell before running the MPI executable:

```bash
export PATH=$HOME/nachevsky/test/simulation/local/openmpi-5.0.8/bin:$PATH
export LD_LIBRARY_PATH=$HOME/nachevsky/test/simulation/local/openmpi-5.0.8/lib:${LD_LIBRARY_PATH:-}
```

Verify that the shell sees the correct runtime:

```bash
which mpirun
ldd ./main3d.gnu.MPI.CUDA.ex | grep libmpi
```

Expected output should point into your local install, for example:

```bash
/home/jovyan/nachevsky/test/simulation/local/openmpi-5.0.8/bin/mpirun
libmpi.so.40 => /home/jovyan/nachevsky/test/simulation/local/openmpi-5.0.8/lib/libmpi.so.40
```

## Build

### Single-GPU (No MPI)
```bash
make -j 8 USE_CUDA=TRUE USE_MPI=FALSE COMP=gnu CUDA_ARCH=90
```

### Multi-GPU (MPI)
To run on multiple GPUs, you must enable MPI.
```bash
make -j 8 USE_CUDA=TRUE USE_MPI=TRUE COMP=gnu CUDA_ARCH=90
```

For H100 GPUs, `CUDA_ARCH=90` is correct.

## Run

### Single-GPU
```bash
# Optional: use a specific GPU
CUDA_VISIBLE_DEVICES=0 ./main3d.gnu.CUDA.ex params.txt
```

If `params.txt` aborts with a `ParmParse` error about `modes`, use
`params_2gpu.txt` instead, or add `\` line continuations to the multi-line
`modes` entry in `params.txt`.

### Multi-GPU (e.g. 2 GPUs)
Use `mpirun` to launch the MPI version of the executable.
```bash
# Use GPUs 0 and 1
CUDA_VISIBLE_DEVICES=0,1 mpirun -n 2 ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt
```

If you are using the local OpenMPI install above, activate it in the shell
first, then run the command from this directory.

For the vacuum-collapse example, the compressive kick is seeded directly in the
initial data at `t=0` through `wormhole_k_monopole_amplitude`,
`wormhole_k_quadrupole_amplitude`, and `wormhole_k_width` in the params file.
The delayed-kick parameters are still supported by the code for alternate
experiments, but they are disabled by default in the example inputs.

### Stabilization Parameters: `kappa`, `sigma`, and `shift_advec_coeff`

Three commonly tuned stabilization parameters in the wormhole runs are `kappa1`,
`sigma`, and `shift_advec_coeff`, but they act in completely different ways.

- `kappa1` is a CCZ4 constraint-damping parameter. It damps violations of the
  Z4/CCZ4 constraint variables during evolution. Increasing `kappa1` means
  stronger damping of those constraint-violating modes.
- `sigma` is the Kreiss-Oliger numerical dissipation coefficient. It adds
  high-frequency smoothing to the finite-difference RHS and mainly suppresses
  grid-scale numerical noise.
- `shift_advec_coeff` controls whether the spatial shift vector is advected by
  itself ($\beta^j \partial_j \beta^i$). In strong-field / moving-puncture
  simulations, the shift vector gets very steep. Multiplying a steep shift by
  its own derivative creates extreme numerical gradients. 

So while they all affect stability, they are not interchangeable:

- raise `kappa1` if the run looks dominated by growing CCZ4 constraint
  violations,
- raise `sigma` if the run looks contaminated by short-wavelength numerical
  noise or interpolation junk,
- set `shift_advec_coeff = 0.0` (disabled) instead of `1.0` if the run suffers
  from extreme gauge-shock NaNs near the throat/puncture. Dropping this term
  is technically a gauge "hack", but it is a standard necessity for surviving
  violent collapses in the Gamma-driver gauge.

`kappa1` is not a physical observable like mass or throat radius; it is a
formulation parameter. So using `kappa1 = 3.0` is not "unphysical" by itself.
It simply means stronger CCZ4 constraint damping than `kappa1 = 0.1`, and can
be a reasonable choice if it improves stability without spoiling convergence or
the qualitative behavior of the solution.

### Multi-GPU (e.g. 4 GPUs)
This example also works with 4 GPUs (4 MPI ranks):

```bash
CUDA_VISIBLE_DEVICES=0,1,2,3 mpirun -n 4 ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt
```

You may see a warning like:
“Multiple GPUs are visible to each MPI rank... rank-to-GPU mapping...”.
It can still run correctly, but for the most robust rank→GPU mapping, use the
binding method below.

### Multi-GPU (e.g. 8 GPUs)
On an 8×H100 node you can run with 8 GPUs (8 MPI ranks) like this:

```bash
CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7 mpirun -n 8 ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt
```

#### Quick MPI troubleshooting

Before running the MPI executable, it is worth checking:

```bash
which mpirun
ldd ./main3d.gnu.MPI.CUDA.ex | grep libmpi
```

If `which mpirun` is empty, or `libmpi.so.40` is reported as `not found`, the
OpenMPI runtime is not active in your shell yet.

Run on some random gpus 
```bash
CUDA_VISIBLE_DEVICES=4,5,6,7 mpirun -n 8 ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt
```

#### Recommended: bind MPI ranks to GPUs (avoids mapping warnings)
AMReX may warn that multiple GPUs are visible to each rank. A robust way to
ensure rank 0 uses GPU 0 and rank 1 uses GPU 1 is to make each rank see only a
single GPU:

```bash
mpirun -n 2 bash -c 'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; exec ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt'
```

For 4 GPUs (ranks 0–3 → GPUs 0–3):

```bash
mpirun -n 4 bash -c 'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; exec ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt'
```

For 8 GPUs (ranks 0–7 → GPUs 0–7):

```bash
mpirun -n 8 bash -c 'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; exec ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt'
```

#### Will it be faster with more GPUs?
It depends — not automatically.

- You will almost always get more total throughput than 1 GPU, but **8 GPUs will not be 8× faster** for this kind of AMR code.
- With many small grids (e.g. lots of `8×8×8` up to `32×32×32`), scaling is often limited by:
  - MPI halo exchanges (more ranks → more communication),
  - regridding / load-balancing overhead,
  - kernel launch overhead on small tiles.

**What you should expect**
- **2 → 4 GPUs** often gives a noticeable speedup.
- **4 → 8 GPUs** may give a smaller gain, and sometimes can be close to flat if communication dominates.

**How to know for sure (quick test)**
Run the same setup for a short time (e.g. set `stop_time = 1.0`) with `-n 2`, `-n 4`, `-n 8` and compare the printed “average evolution speed” / walltime per step. That’s the only reliable answer on your specific problem + node.

#### If it looks “stuck” after `AMReX ... initialized`
It can spend noticeable time on **CPU-side initialization** (grid generation,
distribution mapping, regridding setup) before the first evolution step prints.
During this phase GPU utilization can be low even though the run is progressing.

If you want to see which rank is printing what, you can also run:

```bash
CUDA_VISIBLE_DEVICES=0,1,2,3 mpirun --tag-output -n 4 ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt
```

#### Why GPU memory can look identical on 1 vs 2 GPUs
AMReX uses a GPU memory pool (arena). If you set
`amrex.the_arena_init_size` very large, **each MPI rank (each GPU) will
pre-allocate that amount**, so `nvidia-smi` can show similarly high memory use
per GPU even when you add GPUs.

- In `params_2gpu.txt` we set a conservative default:
  `amrex.the_arena_init_size = 20000000000` (20GB per GPU).
- If you previously used something like 70GB, that would reserve ~70GB on *each*
  GPU and can still lead to OOM because it leaves little headroom.

You can override it at runtime without editing files, e.g.:

```bash
mpirun -n 2 bash -lc 'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; exec ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt amrex.the_arena_init_size=20000000000'
```

#### If it segfaults right after reading inputs
Some MPI builds are **not CUDA-aware**. If you see a segfault during AMReX MPI
communication (e.g. in `FabArray::ParallelCopy` / `Waitall`), disable GPU-aware
MPI:

- In the params file: `amrex.use_gpu_aware_mpi = 0` (recommended default here), or
- As a one-off override on the command line:

```bash
CUDA_VISIBLE_DEVICES=0,1 mpirun -n 2 ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt amrex.use_gpu_aware_mpi=0
```
