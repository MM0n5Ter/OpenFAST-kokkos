#!/bin/bash
#SBATCH -J ofk-biotsavart-direct-gpu
#SBATCH -p h100,h200
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH -t 00:20:00
#SBATCH -o gpu.%j.log
# biotsavart-direct v1 on one H100/H200 (sm_90), Kokkos-CUDA build of tests/biotsavart-direct.
# Submit from the CUDA test build tree:
#   cd $OPT/build/biotsavart-direct-v1-cuda && sbatch $REPO/tests/biotsavart-direct/gpu.sbatch
# Chain claim 2 (strict GPU vs Fortran OPENMP=OFF, which is strict CPU by claim 1), claim 3 on
# the GPU (strict vs fast, same backend), R7 on the GPU (run to run, in-process and across
# processes), and timings against the Fortran -fopenmp build on the job's 8 cores.
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28
SRC=$(dirname "$(readlink -f "$0")")
[ -f "$SRC/r7_threads.sh" ] || SRC=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos/tests/biotsavart-direct
export OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=threads
hostname; nvidia-smi --query-gpu=name,driver_version --format=csv,noheader; date

echo "=== chain claim 2: Fortran (OPENMP=OFF) vs strict on CUDA; err / sum|contributions|"
for exe in l0_ui_seg l0_ui_part_nograd; do
  for seed in 1 7; do
    OFK_MODE=strict OFK_L0_SEED=$seed OFK_L0_RTOL=1e-14 ./$exe | tail -7
  done
done
echo "=== chain claim 3 on CUDA: strict vs fast (shape auto, then forced skinny / square for ui_seg)"
for seed in 1 7; do
  OFK_L0_STRICTFAST=1 OFK_L0_SEED=$seed OFK_L0_RTOL=1e-12 ./l0_ui_part_nograd | tail -5
  for shape in auto skinny square; do
    echo "--- ui_seg OFK_BSD_SHAPE=$shape seed=$seed"
    OFK_BSD_SHAPE=$shape OFK_L0_STRICTFAST=1 OFK_L0_SEED=$seed OFK_L0_RTOL=1e-12 ./l0_ui_seg | grep -v "^FAIL" | tail -7
    OFK_BSD_SHAPE=$shape OFK_L0_STRICTFAST=1 OFK_L0_SEED=$seed OFK_L0_RTOL=1e-12 ./l0_ui_seg | grep -c "^FAIL" | sed 's/^/    failing cases at 1e-12: /'
  done
done
echo "=== R7 on CUDA: port vs port in-process (rtol 0), both modes"
for m in strict fast; do
  for exe in l0_ui_seg l0_ui_part_nograd; do
    OFK_MODE=$m OFK_L0_PORTPORT=1 OFK_L0_RTOL=0 ./$exe | tail -1
  done
done
echo "=== R7 on CUDA: across processes (dumps, 4 runs; the thread count only affects the host)"
for m in strict fast; do
  for exe in l0_ui_seg l0_ui_part_nograd; do
    bash "$SRC/r7_threads.sh" ./$exe $m r7gpu | tail -1
  done
done
echo "=== ctest"
ctest --output-on-failure -E "api_fvw" 2>&1 | tail -30
ctest -R api_fvw --output-on-failure 2>&1 | tail -5
echo "=== timings: routine,nCP,nSrc,reg,ref_s(Fortran -fopenmp, 8 cores),strict_s,fast_s,strict_speedup,fast_speedup,relmax_strict,relmax_fast"
for r in 0 2 3; do ./perf_direct seg 10000 10000 $r 3; done
for r in 0 1 2; do ./perf_direct part 10000 10000 $r 3; done
./perf_direct seg 150 300000 3 3
OFK_BSD_PROFILE=1 ./perf_direct seg 150 300000 3 1 2>&1 | grep "fast\|^seg" | tail -3
./perf_direct seg 100000 100000 2 1
./perf_direct part 100000 100000 1 1
OFK_BSD_PROFILE=1 ./perf_direct part 100000 100000 1 1 2>&1 | grep "\[ofk_bsd\]" | tail -2
date
