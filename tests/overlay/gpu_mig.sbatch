#!/bin/bash
# "No change in behaviour" evidence for the overlay layout (notes/overlay.md), GPU part: L0 of units 1
# and 2 on the overlay-built CUDA test trees (ov-tests-bsd-cuda-a100, ov-tests-iv-cuda-a100) on an
# A100 MIG slice. Compare with job 1862373 (a100mig.1862373.log, fork-built trees).
#   sbatch -J ofk-infra-gpu tests/overlay/gpu_mig.sbatch
#SBATCH -p interactive_gpu
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 01:30:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/ov-gpu.%x.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
export LD_LIBRARY_PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
OUT=$B/ov-gpu-$SLURM_JOB_ID; mkdir -p $OUT
export OMP_PROC_BIND=false OMP_NUM_THREADS=4
hostname; nvidia-smi -L; date
echo "=== unit 1 L0 on GPU: the test set of job 1862373 (seed 1)"
cd $B/ov-tests-bsd-cuda-a100
ctest -E "kkoff|fomp|ffma|refref|api_fvw" --timeout 900 > $OUT/bsd.old-set.ctest.log 2>&1
grep -E "Test +#|tests passed" $OUT/bsd.old-set.ctest.log
grep -hE "^=== |^  reg |passed, .* failed" Testing/Temporary/LastTest.log | head -60
cp Testing/Temporary/LastTest.log $OUT/bsd.old-set.LastTest.log
for seed in 1 7; do
  echo "=== unit 1 full ctest on GPU, OFK_L0_SEED=$seed (claim 3 bound = max(1e-12, 1.1 x Fortran FMA spread))"
  OFK_L0_SEED=$seed ctest --timeout 900 > $OUT/bsd.seed$seed.ctest.log 2>&1
  grep -E "Test +#|tests passed" $OUT/bsd.seed$seed.ctest.log
  grep -hE "^  reg .* bound=|^l0: per-branch" Testing/Temporary/LastTest.log
  cp Testing/Temporary/LastTest.log $OUT/bsd.seed$seed.LastTest.log
done
echo "=== unit 2 L0 on GPU: the test set of job 1862373"
cd $B/ov-tests-iv-cuda-a100
ctest -E "kkoff|ffma|refref" --timeout 1200 > $OUT/iv.ctest.log 2>&1
grep -E "Test +#|tests passed" $OUT/iv.ctest.log
grep -hE "^l0_iv: .*passed|^  (Basic|Part|TreePart|A7|pack) |P6 \(L0|R7 |FAIL" Testing/Temporary/LastTest.log | head -80
cp Testing/Temporary/LastTest.log $OUT/iv.LastTest.log
date
