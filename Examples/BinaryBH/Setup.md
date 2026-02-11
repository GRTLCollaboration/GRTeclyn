Run build 
make -j 8 USE_CUDA=TRUE USE_MPI=FALSE COMP=gnu CUDA_ARCH=90

# 1. Make sure your data folder is empty if you want a fresh start
rm -rf /home/jovyan/nachevsky/test/simulation/data/*

# 2. Run the simulation
CUDA_VISIBLE_DEVICES=0 ./main3d.gnu.CUDA.ex params.txt