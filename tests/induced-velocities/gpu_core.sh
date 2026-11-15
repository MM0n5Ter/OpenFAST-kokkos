#!/bin/bash
# Short GPU plan (backfill-sized, <30 min): the core evidence of gpu_steps.sh only.
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
( cd $B/iv-cuda && ctest -E "kkoff|ffma|refref" --timeout 600 > $OUT/l0.ctest.log 2>&1 )
grep -E "Test +#|tests passed" $OUT/l0.ctest.log
grep -E "^l0_iv: .*passed|^  (Basic|Part|A7|pack) |P6 \(L0\)|R7 |FAIL" $B/iv-cuda/Testing/Temporary/LastTest.log | head -80
D=$B/fork-cuda/modules/aerodyn/aerodyn_driver; C=$B/iv-cases; MT=$B/iv-omp/libofk_memtool.so
bash $T/run_case.sh $D $C/bar-r8-direct  $OUT/cases/bar-r8-direct.fast.prof OFK_MODE=fast OFK_IV_PROFILE=1 KOKKOS_TOOLS_LIBS=$MT OFK_MEMTOOL_OUT=$OUT/memtool.bar-r8-direct.fast.txt
tail -1 $OUT/memtool.bar-r8-direct.fast.txt
bash $T/run_case.sh $D $C/bar-r8-direct  $OUT/cases/bar-r8-direct.fast   OFK_MODE=fast
bash $T/run_case.sh $D $C/bar-r8-direct  $OUT/cases/bar-r8-direct.strict OFK_MODE=strict
bash $T/run_case.sh $D $C/bar-r8-tree    $OUT/cases/bar-r8-tree.fast     OFK_MODE=fast
bash $T/run_case.sh $D $C/bar-perf-direct $OUT/cases/bar-perf-direct.fast.prof OFK_MODE=fast OFK_IV_PROFILE=1
bash $T/run_case.sh $D $C/bar-perf-direct $OUT/cases/bar-perf-direct.fast OFK_MODE=fast
date
