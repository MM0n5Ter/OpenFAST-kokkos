#!/bin/bash
#SBATCH -J ofk-biotsavart-direct-gpu2
#SBATCH -p h100,h200
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH -t 00:30:00
#SBATCH -o gpu2.%j.log
# biotsavart-direct v1, second GPU run:
#  1. R7 across processes on CUDA (bitwise dumps of the port outputs, 4 runs per mode);
#  2. the full fork built with Kokkos-CUDA (fork-cuda): aerodyn_utest and the OLAF regression
#     tests (without the openfast 5MW case and the Python-binding cases), fast then strict,
#     plus the fork-off reference outputs for a .outb comparison.
#   sbatch --chdir=$OPT/build $REPO/tests/biotsavart-direct/gpu2.sbatch
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
SRC=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos/tests/biotsavart-direct
OUT=$B/gpu2-$SLURM_JOB_ID; mkdir -p $OUT
export OMP_PROC_BIND=spread OMP_PLACES=threads
hostname; nvidia-smi --query-gpu=name --format=csv,noheader; date
echo "=== R7 on CUDA across processes"
cd $B/biotsavart-direct-v1-cuda
for m in strict fast; do
  for exe in l0_ui_seg l0_ui_part_nograd; do bash $SRC/r7_threads.sh ./$exe $m $OUT/r7 | tail -5; done
done
R='aerodyn_utest|OLAF'; E='5MW|_py$|^py_'
for m in fast strict; do
  echo "=== fork-cuda, OFK_MODE=$m: ctest -R '$R' -E '$E'"
  ( cd $B/fork-cuda && OMP_NUM_THREADS=1 OFK_MODE=$m ctest -R "$R" -E "$E" -j 4 --timeout 1500 > $OUT/cuda-$m.ctest.log 2>&1 )
  grep -E "Test +#|tests passed" $OUT/cuda-$m.ctest.log
  grep -h "OFK_MODE" $B/fork-cuda/reg_tests/modules/aerodyn/ad_B1n2_OLAF/*.log 2>/dev/null | head -1
  mkdir -p $OUT/cuda-$m; ( cd $B/fork-cuda/reg_tests && find . -name "*.outb" -newer $OUT -print0 | xargs -0 -I{} cp --parents {} $OUT/cuda-$m/ )
done
date
