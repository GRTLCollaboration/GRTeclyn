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

To resume from any checkpoint (for example step **3000** or **4000**):

1. **Rollback** the data directory to the desired checkpoint (use paths matching your clone of this repo and the `output_path` set in `params_2gpu.txt`):
   ```bash
   cd /path/to/GRTeclyn
   ./grteclyn-wrapper/scripts/rollback --step 3000 --data /path/to/your/simulation/output
   ```
   (Replace `3000` with the desired checkpoint number, e.g. `4000`.)

2. **Start the watcher** (in a second terminal — use `--not-remove`):
   ```bash
   cd /path/to/GRTeclyn
   ./grteclyn-wrapper/scripts/plot_run.sh --not-remove /path/to/your/simulation/output
   ```

3. **Restart the simulation** using `amr.restart` (this is the correct flag):
   ```bash
   cd /path/to/GRTeclyn/Examples/SupportedWormholeCollapse

   # 2 GPUs example
   CUDA_VISIBLE_DEVICES=0,1 mpirun -n 2 ./main3d.gnu.MPI.CUDA.ex params_2gpu.txt \
     amr.restart=/path/to/your/simulation/output/SupportedWormholeChk03000
   ```

   **Single GPU example:**
   ```bash
   CUDA_VISIBLE_DEVICES=0 ./main3d.gnu.CUDA.ex params_2gpu.txt \
     amr.restart=/path/to/your/simulation/output/SupportedWormholeChk03000
   ```

> **Tip**: Always replace both the `--step N` in rollback **and** the `ChkN` number in the restart command with the same value (e.g. 4000).

## Plotting (normal run)

```bash
# In a second terminal (from the GRTeclyn root; set RUN_DIR to your simulation output)
./grteclyn-wrapper/scripts/plot_run.sh /path/to/your/simulation/output
```

## Quick MPI Setup (if needed)

Point these at **your** OpenMPI install prefix (the directory that contains `bin/mpirun` and `lib`):

```bash
export PATH=/path/to/your/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/path/to/your/openmpi/lib:${LD_LIBRARY_PATH:-}
```

