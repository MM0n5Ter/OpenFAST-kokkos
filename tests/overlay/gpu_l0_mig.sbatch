#!/bin/bash
# L0 of a CUDA test tree on an A100 MIG slice (no queue): ctest in the given tree.
#   sbatch -J ofk-infra-<what> tests/overlay/gpu_l0_mig.sbatch <build tree under opt/build>
#SBATCH -p interactive_gpu
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 00:45:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/%x.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
export LD_LIBRARY_PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH
export OMP_PROC_BIND=false OMP_NUM_THREADS=4
hostname; nvidia-smi -L
cd /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/$1 && ctest --timeout 1800 2>&1 | grep -E "Test +#|tests passed|Failed"
