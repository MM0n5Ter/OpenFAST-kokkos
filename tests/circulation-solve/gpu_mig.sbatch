#!/bin/bash
# Unit circulation-solve on a CUDA build (A100 MIG 2g.20gb, interactive_gpu): the kernel runs in the Kokkos
# HOST space (OpenMP) of the CUDA build; this checks that the CUDA build and link work and give the same answers.
#   sbatch -J ofk-circulation-solve-mig tests/circulation-solve/gpu_mig.sbatch
# Trees: u4a-cuda-a100 (L0 test tree, nvcc_wrapper), ${OFK_KC:-ov4-cuda-a100} (aerodyn_driver).
#SBATCH -p interactive_gpu
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 01:30:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u4a-jobs/%x.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
export LD_LIBRARY_PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
TI=$R/tests/induced-velocities; L=$R/deps/openfast/reg_tests/lib
PY=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin/python
OUT=$B/u4a-mig-$SLURM_JOB_ID; mkdir -p $OUT
cp $R/tests/circulation-solve/iters_cmp.py $TI/run_case.sh $TI/r8.py $OUT/; H=$OUT
export OMP_PROC_BIND=false OMP_NUM_THREADS=4
hostname; nvidia-smi -L; lscpu | grep "Model name"; date
echo "=== L0 on the CUDA build (host space): LL by the Fortran segment tree (isolates this unit), seeds 1 and 7"
for seed in 1 7; do
  OFK_L0_SEED=$seed OFK_L0_LLVEL=4 $B/u4a-cuda-a100/l0_cs > $OUT/l0.llvel4.s$seed.log 2>&1; echo "  seed $seed rc=$?"
  grep -E "variant|pin with FMA|knife-edge ConvCrit|A7|P6|l0_cs:" $OUT/l0.llvel4.s$seed.log | cut -c1-170 | sed 's/^/    /'
done
echo "=== L0 on the CUDA build, LL by unit 2's GPU kernel (strict on the GPU: claim 2, a few ulp; informational)"
OFK_L0_LLVEL=1 OFK_L0_INFO=1 $B/u4a-cuda-a100/l0_cs > $OUT/l0.llvel1.log 2>&1
grep -E "variant|l0_cs:" $OUT/l0.llvel1.log | cut -c1-170 | sed 's/^/    /'
echo "=== bar-r8-tree on ov4-cuda-a100: F0 (fast), strict; iteration counts vs the pin (job: u4a-iters); R8"
c=bar-r8-tree; KC=$B/${OFK_KC:-ov4-cuda-a100}/modules/aerodyn/aerodyn_driver
bash $H/run_case.sh $KC $B/iv-cases/$c $OUT/$c.F0 OFK_MODE=fast OFK_CS_FMA=0 OFK_CS_PROFILE=1 OFK_CS_ITERS=$OUT/$c.F0.iters
bash $H/run_case.sh $KC $B/iv-cases/$c $OUT/$c.F0b OFK_MODE=fast OFK_CS_FMA=0 OFK_CS_ITERS=$OUT/$c.F0b.iters
bash $H/run_case.sh $KC $B/iv-cases/$c $OUT/$c.strict OFK_MODE=strict OFK_CS_ITERS=$OUT/$c.strict.iters
bash $H/run_case.sh $KC $B/iv-cases/$c $OUT/$c.F0-csfortran OFK_MODE=fast OFK_CS=fortran
grep -h "ofk_cs\]" $OUT/$c.F0/stdout.log | head -3
same() { cmp -s <(tail -n +8 $1) <(tail -n +8 $2) && echo IDENTICAL || echo DIFFER; }
echo "--- R7 (GPU, two processes): $(same $OUT/$c.F0/ad_driver.out $OUT/$c.F0b/ad_driver.out)"
echo "--- this unit alone on the CUDA build: F0 vs OFK_CS=fortran (same GPU kernels elsewhere): $(same $OUT/$c.F0/ad_driver.out $OUT/$c.F0-csfortran/ad_driver.out)"
PIN=$(ls -d $B/u4a-iters-*/ 2>/dev/null | tail -1)
[ -n "$PIN" ] && $PY $H/iters_cmp.py $PIN/$c.pin.iters F0-gpu=$OUT/$c.F0.iters strict-gpu=$OUT/$c.strict.iters
$PY $H/r8.py $L $B/iv-r8-1860700/$c.fomp1/ad_driver.out $B/iv-r8ffma-1860993/$c.ffma/ad_driver.out \
    gpu-F0=$OUT/$c.F0/ad_driver.out gpu-strict=$OUT/$c.strict/ad_driver.out
echo "=== ad_Sphere_OLAF (no lifting-line control point: nCP_tot = 0) on the CUDA build, fast and strict"
RT=$R/deps/openfast/reg_tests/r-test/modules/aerodyn
for m in fast strict; do
  d=$OUT/sphere.$m; mkdir -p $d; for x in $RT/*; do ln -sfn $x $d/; done; rm $d/ad_Sphere_OLAF; cp -r $RT/ad_Sphere_OLAF $d/
  ( cd $d/ad_Sphere_OLAF && OFK_MODE=$m $KC ad_driver.dvr > stdout.log 2>&1; echo "  $m rc=$?"; tail -2 stdout.log )
done
date
