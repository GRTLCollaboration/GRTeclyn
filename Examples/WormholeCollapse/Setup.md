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

### Multi-GPU (e.g. 2 GPUs)
Use `mpirun` to launch the MPI version of the executable.
```bash
# Use GPUs 0 and 1
CUDA_VISIBLE_DEVICES=0,1 mpirun -n 2 ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt
```

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
