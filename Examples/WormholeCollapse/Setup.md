# WormholeCollapse Setup

## Prerequisites

**CUDA must be in your PATH.** The build uses `nvcc` directly (not `CUDA_HOME`). If you see `nvcc: not found`:

```bash
export PATH=/usr/local/cuda/bin:$PATH
# Or for CUDA 12: export PATH=/usr/local/cuda-12/bin:$PATH
```

Add to `~/.bashrc` if needed:
```bash
echo 'export PATH=/usr/local/cuda/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

## Build

```bash
make -j 8 USE_CUDA=TRUE USE_MPI=FALSE COMP=gnu CUDA_ARCH=90
```

For H100 GPUs, `CUDA_ARCH=90` is correct.

## Run

```bash
# Optional: use a specific GPU
CUDA_VISIBLE_DEVICES=0 ./main3d.gnu.CUDA.ex params.txt
```
