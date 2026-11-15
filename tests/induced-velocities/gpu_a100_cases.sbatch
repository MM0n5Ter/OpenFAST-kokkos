#!/bin/bash
# GPU correctness for units 1 and 2 on an A100 MIG slice (interactive_gpu: no queue, 4 h limit).
# Correctness only (chain claims 2/3 on GPU, R7 on GPU, P6 on GPU, R8 on GPU). Timings need a full GPU (h100 jobs).
#SBATCH -J ofk-induced-velocities-a100cases
#SBATCH -p interactive_gpu
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 01:30:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/a100cases.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
T=$R/tests/induced-velocities; C=$B/iv-cases
OUT=$B/a100cases-$SLURM_JOB_ID; mkdir -p $OUT
PY=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin/python; L=$R/deps/openfast/reg_tests/lib
export OMP_PROC_BIND=false OMP_NUM_THREADS=4
hostname; nvidia-smi -L; date
# The cluster openblas is a zen4 (AVX-512) build; A100 nodes are Zen3. Use the portable build (Haswell baseline, dynamic arch).
export LD_LIBRARY_PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH
echo "=== unit 2 full cases on GPU (fork-cuda-a100 aerodyn_driver)"
D=$B/fork-cuda-a100/modules/aerodyn/aerodyn_driver; MT=$B/iv-cuda-a100/libofk_memtool.so
bash $T/run_case.sh $D $C/bar-r8-direct $OUT/bar-r8-direct.gpu-fast-mem OFK_MODE=fast OFK_IV_PROFILE=1 KOKKOS_TOOLS_LIBS=$MT OFK_MEMTOOL_OUT=$OUT/memtool.direct.txt
tail -1 $OUT/memtool.direct.txt
for c in bar-r8-direct bar-r8-tree; do
  for m in strict fast; do bash $T/run_case.sh $D $C/$c $OUT/$c.gpu-$m OFK_MODE=$m; done
  bash $T/run_case.sh $D $C/$c $OUT/$c.gpu-fast-b OFK_MODE=fast
  echo "--- R7 (GPU fast, two processes): $(cmp -s <(tail -n +8 $OUT/$c.gpu-fast/ad_driver.out) <(tail -n +8 $OUT/$c.gpu-fast-b/ad_driver.out) && echo identical || echo DIFFER)"
  echo "--- R8 $c on GPU: base1 = Fortran 1 thread (job 1860700), envelope = Fortran built with FMA (job 1860993)"
  $PY $T/r8.py $L $B/iv-r8-1860700/$c.fomp1/ad_driver.out $B/iv-r8ffma-1860993/$c.ffma/ad_driver.out \
      gpu-strict=$OUT/$c.gpu-strict/ad_driver.out gpu-fast=$OUT/$c.gpu-fast/ad_driver.out 2>&1 | grep -E "max drift|envelope max"
done
date
