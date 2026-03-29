Run build 
make -j 8 USE_CUDA=TRUE USE_MPI=FALSE COMP=gnu CUDA_ARCH=90

# 1. Make sure your data folder is empty if you want a fresh start (set DATA_DIR to match params.txt output_path)
rm -rf /path/to/your/simulation/data/*

# 2. Run the simulation
CUDA_VISIBLE_DEVICES=0 ./main3d.gnu.CUDA.ex params.txt