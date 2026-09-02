#!/bin/bash

#SBATCH -J a100-regression-test
#SBATCH -p ampere
#SBATCH -A SHELLARD-SL3-GPU
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:1
#SBATCH --time=00:10:00
#SBATCH --mail-type=NONE
#SBATCH --no-requeue
#SBATCH --qos=INTR

#Pass environment variables as
#sbatch --export=MODULEPATHS=${MODULEPATHS_A100},MODULES=${MODULES_A100} jobs.sb

set -euo pipefail

cd ${HOME}/${CI_PROJECT_DIR}/Tests
srun ./Tests3d.gnu.DEBUG.MPI.CUDA.ex -dt-d=yes

cd ${HOME}/${CI_PROJECT_DIR}/Examples/BinaryBH
srun ./BinaryBH3d.gnu.DEBUG.MPI.CUDA.ex ./params_test.txt amr.plot_file=${OUTPUT_DIR}/BinaryBH_

cd ${HOME}/${CI_PROJECT_DIR}/Examples/KleinGordon
srun ./KleinGordon3d.gnu.DEBUG.MPI.CUDA.ex ./params_test.txt amr.plot_file=${OUTPUT_DIR}/KleinGordon_
