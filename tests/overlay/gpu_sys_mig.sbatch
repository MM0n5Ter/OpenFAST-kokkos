#!/bin/bash
# Overlay layout, GPU system level (notes/overlay.md §4.3): aerodyn_driver r-tests on the overlay
# CUDA build ov-cuda-a100 (A100 MIG), strict and fast; .outb channel data vs the CPU OFF build.
#   sbatch -J ofk-infra-gpusys tests/overlay/gpu_sys_mig.sbatch <dir with off/ .outb files, e.g. opt/build/ov-sys-<job>/off>
#SBATCH -p interactive_gpu
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 01:30:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/ov-gpu.%x.%j.log
set -u
REF=$1
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
export LD_LIBRARY_PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OMP_PROC_BIND=false
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
PY=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin/python; L=$R/deps/openfast/reg_tests/lib
OUT=$B/ov-gpusys-$SLURM_JOB_ID; mkdir -p $OUT
hostname; nvidia-smi -L | grep MIG; date
RX='aerodyn_utest|ad_.*OLAF'; EX='py_'
for m in strict fast; do
  touch $OUT/stamp.$m
  ( cd $B/ov-cuda-a100 && OFK_MODE=$m ctest -R "$RX" -E "$EX" -j 2 --timeout 2400 > $OUT/cuda-$m.ctest.log 2>&1 )
  echo "=== ov-cuda-a100 OFK_MODE=$m"; grep -E "Test +#|tests passed" $OUT/cuda-$m.ctest.log
  mkdir -p $OUT/cuda-$m; ( cd $B/ov-cuda-a100/reg_tests && find . -name "*.outb" -newer $OUT/stamp.$m -print0 | xargs -0 -I{} cp --parents {} $OUT/cuda-$m/ )
done
echo "=== .outb channel data vs the CPU OFF build ($REF)"
$PY $R/tests/biotsavart-direct/cmp_outb.py $L $REF $OUT/cuda-strict $OUT/cuda-fast
date
