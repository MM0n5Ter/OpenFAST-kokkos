#!/bin/bash
# Build a fork tree in a job: sbatch -J ofk-induced-velocities-build-<x> build.sbatch <tree under opt/build> [targets...]
#SBATCH -p med
#SBATCH -c 8
#SBATCH -t 01:00:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/iv-build.%x.%j.log
set -eu
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28
export NVCC_WRAPPER_DEFAULT_COMPILER=g++
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
tree=$1; shift
if [ $# -gt 0 ]; then cmake --build $B/$tree -j 8 --target "$@"; else cmake --build $B/$tree -j 8; fi
echo "build $tree OK"; date
