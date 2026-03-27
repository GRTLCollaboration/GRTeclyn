# SupportedWormholeCollapse

## Build

```bash
# Single GPU
make -j 8 USE_CUDA=TRUE USE_MPI=FALSE COMP=gnu CUDA_ARCH=90

# Multi-GPU
make -j 8 USE_CUDA=TRUE USE_MPI=TRUE COMP=gnu CUDA_ARCH=90
```

## Run

**Single GPU:**
```bash
CUDA_VISIBLE_DEVICES=0 ./main3d.gnu.CUDA.ex params_2gpu.txt
```

**Multi-GPU:**

```bash
# 2 GPUs
CUDA_VISIBLE_DEVICES=0,1 mpirun -n 2 ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt

# 4 GPUs
CUDA_VISIBLE_DEVICES=0,1,2,3 mpirun -n 4 ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt

# 8 GPUs
CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7 mpirun -n 8 ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt
```

**Recommended GPU binding (avoids mapping warnings):**
```bash
mpirun -n 8 bash -c 'export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK; exec ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt'
```

## Restart from Checkpoint

To resume from a checkpoint (e.g. step 3000):

1. **Rollback** the data directory (removes later plotfiles, truncates logs, cleans frames):
   ```bash
   cd /home/jovyan/nachevsky/test/simulation/GRTeclyn
   ./src/scripts/rollback --step 3000 --data /home/jovyan/nachevsky/test/simulation/data_supported
   ```

2. **Start watcher** (in a second terminal):
   ```bash
   cd /home/jovyan/nachevsky/test/simulation/GRTeclyn
   ./src/scripts/plot_run.sh --not-remove /home/jovyan/nachevsky/test/simulation/data_supported
   ```

3. **Restart simulation** using `amr.restart` (critical):
   ```bash
   cd /home/jovyan/nachevsky/test/simulation/GRTeclyn/Examples/SupportedWormholeCollapse

   # 2 GPUs
   CUDA_VISIBLE_DEVICES=0,1 mpirun -n 2 ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt \
     amr.restart=/home/jovyan/nachevsky/test/simulation/data_supported/SupportedWormholeChk03000
   ```

   **Single GPU:**
   ```bash
   CUDA_VISIBLE_DEVICES=0 ./main3d.gnu.CUDA.ex params_2gpu.txt \
     amr.restart=/home/jovyan/nachevsky/test/simulation/data_supported/SupportedWormholeChk03000
   ```

## Plotting (normal run)

```bash
# In a second terminal
./src/scripts/plot_run.sh /home/jovyan/nachevsky/test/simulation/data_supported
```

## Quick MPI Setup (if needed)

```bash
export PATH=$HOME/nachevsky/test/simulation/local/openmpi-5.0.8/bin:$PATH
export LD_LIBRARY_PATH=$HOME/nachevsky/test/simulation/local/openmpi-5.0.8/lib:${LD_LIBRARY_PATH:-}
```

