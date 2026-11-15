#!/bin/bash
# GPU plan of unit induced-velocities (run by gpu.sbatch; output dir $B/iv-gpu-<jobid>):
#  1. L0 on CUDA (claims 2 and 3, R7 in process and across processes, A7, P6 counters), both
#     source paths (packed, lattice);
#  2. fork-cuda: aerodyn_utest + OLAF r-tests (no openfast 5MW, no Python), fast then strict;
#  3. scratch BAR cases (make_cases.sh): direct wake (the GPU carries the wake), tree wake;
#     fast and strict; OFK_IV_PROFILE (P3 split) and Kokkos-Tools memory events (P6).
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28
# Cluster openblas is a zen4 (AVX-512) build; use the portable one so this also runs on Zen3 GPU nodes.
export LD_LIBRARY_PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
T=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos/tests/induced-velocities
OUT=$B/iv-gpu-${SLURM_JOB_ID:-manual}; mkdir -p $OUT
export OMP_PROC_BIND=spread OMP_PLACES=threads OMP_NUM_THREADS=8
nvidia-smi --query-gpu=name --format=csv,noheader; date
echo "=== 1. L0 on CUDA"
( cd $B/iv-cuda && ctest --timeout 900 > $OUT/l0.ctest.log 2>&1 )
grep -E "Test +#|tests passed" $OUT/l0.ctest.log
grep -E "^l0_iv: .*passed|^  (Basic|Part|A7) |P6 \(L0\)|R7 |FAIL" $B/iv-cuda/Testing/Temporary/LastTest.log | head -80
D=$B/fork-cuda/modules/aerodyn/aerodyn_driver
C=$B/iv-cases
MT=$B/iv-omp/libofk_memtool.so
echo "=== 3. scratch BAR cases on the GPU"
bash $T/run_case.sh $D $C/bar-r8-direct  $OUT/cases/bar-r8-direct.fast.prof   OFK_MODE=fast OFK_IV_PROFILE=1 KOKKOS_TOOLS_LIBS=$MT OFK_MEMTOOL_OUT=$OUT/memtool.bar-r8-direct.fast.txt
tail -1 $OUT/memtool.bar-r8-direct.fast.txt
bash $T/run_case.sh $D $C/bar-r8-direct  $OUT/cases/bar-r8-direct.fast        OFK_MODE=fast
bash $T/run_case.sh $D $C/bar-r8-direct  $OUT/cases/bar-r8-direct.strict      OFK_MODE=strict
bash $T/run_case.sh $D $C/bar-r8-tree    $OUT/cases/bar-r8-tree.fast          OFK_MODE=fast
bash $T/run_case.sh $D $C/bar-r8-tree    $OUT/cases/bar-r8-tree.strict        OFK_MODE=strict
bash $T/run_case.sh $D $C/bar-perf-direct $OUT/cases/bar-perf-direct.fast     OFK_MODE=fast
bash $T/run_case.sh $D $C/bar-perf-direct $OUT/cases/bar-perf-direct.fast.prof OFK_MODE=fast OFK_IV_PROFILE=1
bash $T/run_case.sh $D $C/bar-perf-direct $OUT/cases/bar-perf-direct.fast.packed.prof OFK_MODE=fast OFK_IV_PROFILE=1 OFK_IV_SOURCES=packed
bash $T/run_case.sh $D $C/bar-perf-part  $OUT/cases/bar-perf-part.fast.prof    OFK_MODE=fast OFK_IV_PROFILE=1
bash $T/run_case.sh $D $C/bar-perf-direct $OUT/cases/bar-perf-direct.strict   OFK_MODE=strict
R='aerodyn_utest|OLAF'; E='5MW|_py$|^py_'
for m in fast strict; do
  echo "=== 2. fork-cuda, OFK_MODE=$m: ctest -R '$R' -E '$E'"
  touch $OUT/stamp.$m
  ( cd $B/fork-cuda && OMP_NUM_THREADS=1 OFK_MODE=$m ctest -R "$R" -E "$E" -j 4 --timeout 1500 > $OUT/cuda-$m.ctest.log 2>&1 )
  grep -E "Test +#|tests passed" $OUT/cuda-$m.ctest.log
  mkdir -p $OUT/cuda-$m; ( cd $B/fork-cuda/reg_tests && find . -name "*.outb" -newer $OUT/stamp.$m -print0 | xargs -0 -I{} cp --parents {} $OUT/cuda-$m/ )
done
date
