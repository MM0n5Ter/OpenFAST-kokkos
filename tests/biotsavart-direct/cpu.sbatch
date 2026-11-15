#!/bin/bash
#SBATCH -J ofk-biotsavart-direct-cpu
#SBATCH -p med
#SBATCH -c 8
#SBATCH -t 00:40:00
#SBATCH -o cpu.%j.log
# biotsavart-direct v1 on 8 dedicated cores, Kokkos-OpenMP build of tests/biotsavart-direct.
# Submit from the OpenMP test build tree:
#   cd $OPT/build/biotsavart-direct-v1-omp && sbatch $REPO/tests/biotsavart-direct/cpu.sbatch
# Chain claims 1 and 3 (seeds 1 and 7), the Fortran's own FMA envelope, R7 across runs and
# 1/4/8 threads, and P2: Fortran built with -fopenmp vs strict and fast, 8 threads.
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28
SRC=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos/tests/biotsavart-direct
export OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=threads
hostname; lscpu | grep "Model name"; date

for seed in 1 7; do
  echo "=== seed $seed: chain claim 1, Fortran (OPENMP=OFF) vs strict, rtol 0 (bitwise)"
  for exe in l0_ui_seg l0_ui_part_nograd; do OFK_MODE=strict OFK_L0_SEED=$seed OFK_L0_RTOL=0 ./$exe | grep -v "^\[ofk\]" | tail -7; done
  echo "=== seed $seed: chain claim 3, strict vs fast (err / sum|c|), declared R1 bound 1e-12"
  for exe in l0_ui_seg l0_ui_part_nograd; do OFK_L0_STRICTFAST=1 OFK_L0_SEED=$seed OFK_L0_RTOL=1e-12 ./$exe | grep -v "^FAIL\|^\[ofk\]" | tail -7; done
  echo "=== seed $seed: envelope, Fortran built with FMA (-march=x86-64-v3) vs Fortran without (= strict)"
  for exe in l0_ui_seg_ffma l0_ui_part_nograd_ffma; do OFK_MODE=strict OFK_L0_SEED=$seed OFK_L0_RTOL=1 ./$exe | grep -v "^\[ofk\]" | tail -7; done
  echo "=== seed $seed: informational, Fortran built with -fopenmp vs strict"
  for exe in l0_ui_seg_fomp l0_ui_part_nograd_fomp; do OFK_MODE=strict OFK_L0_SEED=$seed OFK_L0_RTOL=1 ./$exe | grep -v "^\[ofk\]" | tail -1; done
done
echo "=== R7: bitwise across processes and OMP_NUM_THREADS=1/4/8"
for m in strict fast; do
  for exe in l0_ui_seg l0_ui_part_nograd; do bash $SRC/r7_threads.sh ./$exe $m r7cpu; done
done
echo "=== ctest"
ctest 2>&1 | tail -12

echo "=== P2: routine,nCP,nSrc,reg,ref_s(Fortran -fopenmp),strict_s,fast_s,strict_speedup,fast_speedup,relmax_strict,relmax_fast"
for r in 0 1 2 3 4; do ./perf_direct seg 10000 10000 $r 5 | grep -v "^\[ofk\]"; done
for r in 0 1 2; do ./perf_direct part 10000 10000 $r 5 | grep -v "^\[ofk\]"; done
./perf_direct seg 150 300000 3 5 | grep -v "^\[ofk\]"
for r in 0 1 2; do ./perf_direct part 100000 100000 $r 1 | grep -v "^\[ofk\]"; done
for r in 0 1 2 3 4; do ./perf_direct seg 100000 100000 $r 1 | grep -v "^\[ofk\]"; done
date
