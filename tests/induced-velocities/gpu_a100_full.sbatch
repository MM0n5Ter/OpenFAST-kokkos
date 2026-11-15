#!/bin/bash
# GPU plan of units 2 and 3 on a full A100 (partition a100), built from the overlay (notes/overlay.md):
#  1. L0 on CUDA: unit 2 (ov-tests-iv-cuda-a100), unit 3 (u3-cuda-a100);
#  2. BAR cases with ov-cuda-a100 (the current code: Kokkos induced velocities + Kokkos treecode),
#     fast and strict; profile (P3 split) and Kokkos-Tools memory events (P6); R7 across processes;
#     R8 against base1 = Fortran 1 thread (job 1860700) and the FMA-built pin envelope (job 1860993);
#  3. the same-node CPU baseline: the plain pin with Fortran OpenMP (openfast-pin-fomp), 8 threads,
#     so P1 compares like with like on this node's EPYC CPU.
#   sbatch tests/induced-velocities/gpu_a100_full.sbatch
#SBATCH -J ofk-induced-velocities-a100full
#SBATCH -p a100
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=96000M
#SBATCH -t 02:30:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/a100full.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
# The cluster openblas is a zen4 (AVX-512) build; A100 nodes are Zen3. Use the portable build.
export LD_LIBRARY_PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
T=$R/tests/induced-velocities; C=$B/iv-cases
OUT=$B/a100full-$SLURM_JOB_ID; mkdir -p $OUT/cases
PY=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin/python; L=$R/deps/openfast/reg_tests/lib
export PYTHONDONTWRITEBYTECODE=1 OMP_PROC_BIND=spread OMP_PLACES=threads OMP_NUM_THREADS=8
hostname; nvidia-smi -L; lscpu | grep "Model name"; date

echo "=== 1. L0 on CUDA"
( cd $B/ov-tests-iv-cuda-a100 && ctest --timeout 900 > $OUT/l0.iv.ctest.log 2>&1 ); grep -E "tests passed|Failed" $OUT/l0.iv.ctest.log
( cd $B/u3-cuda-a100 && ctest -E ffma --timeout 1800 > $OUT/l0.tc.ctest.log 2>&1 ); grep -E "tests passed|Failed" $OUT/l0.tc.ctest.log

echo "=== 2. BAR cases on the GPU (ov-cuda-a100)"
D=$B/ov-cuda-a100/modules/aerodyn/aerodyn_driver; MT=$B/ov-tests-iv-cuda-a100/libofk_memtool.so
bash $T/run_case.sh $D $C/bar-r8-direct $OUT/cases/bar-r8-direct.fast.mem OFK_MODE=fast OFK_IV_PROFILE=1 KOKKOS_TOOLS_LIBS=$MT OFK_MEMTOOL_OUT=$OUT/memtool.direct.txt
tail -1 $OUT/memtool.direct.txt
bash $T/run_case.sh $D $C/bar-r8-tree $OUT/cases/bar-r8-tree.fast.mem OFK_MODE=fast OFK_IV_PROFILE=1 KOKKOS_TOOLS_LIBS=$MT OFK_MEMTOOL_OUT=$OUT/memtool.tree.txt
tail -1 $OUT/memtool.tree.txt
for c in bar-r8-direct bar-r8-tree; do
  for m in strict fast; do bash $T/run_case.sh $D $C/$c $OUT/cases/$c.gpu-$m OFK_MODE=$m; done
  bash $T/run_case.sh $D $C/$c $OUT/cases/$c.gpu-fast-b OFK_MODE=fast
  echo "--- R7 $c (GPU fast, two processes): $(cmp -s <(tail -n +8 $OUT/cases/$c.gpu-fast/ad_driver.out) <(tail -n +8 $OUT/cases/$c.gpu-fast-b/ad_driver.out) && echo identical || echo DIFFER)"
  echo "--- R8 $c on GPU: base1 = Fortran 1 thread (job 1860700), envelope = Fortran built with FMA (job 1860993)"
  $PY $T/r8.py $L $B/iv-r8-1860700/$c.fomp1/ad_driver.out $B/iv-r8ffma-1860993/$c.ffma/ad_driver.out \
      gpu-strict=$OUT/cases/$c.gpu-strict/ad_driver.out gpu-fast=$OUT/cases/$c.gpu-fast/ad_driver.out 2>&1 | grep -E "max drift|envelope max|over"
done
bash $T/run_case.sh $D $C/bar-perf-direct $OUT/cases/bar-perf-direct.gpu-fast        OFK_MODE=fast
bash $T/run_case.sh $D $C/bar-perf-direct $OUT/cases/bar-perf-direct.gpu-fast.prof   OFK_MODE=fast OFK_IV_PROFILE=1
bash $T/run_case.sh $D $C/bar-perf-direct $OUT/cases/bar-perf-direct.gpu-strict      OFK_MODE=strict
bash $T/run_case.sh $D $C/bar-perf-part   $OUT/cases/bar-perf-part.gpu-fast.prof     OFK_MODE=fast OFK_IV_PROFILE=1

echo "=== 3. same-node CPU baseline: plain pin, Fortran OpenMP, 8 threads"
F=$B/openfast-pin-fomp/modules/aerodyn/aerodyn_driver
for c in bar-perf-direct bar-r8-tree bar-r8-direct; do bash $T/run_case.sh $F $C/$c $OUT/cases/$c.fortran-omp8; done

echo "=== wall times"
for d in $OUT/cases/*/; do printf "%-40s %s s\n" "$(basename $d)" "$(cat $d/time.txt 2>/dev/null)"; done
date
