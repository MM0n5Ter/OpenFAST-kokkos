#!/bin/bash
# Check of the default change (kDefaultFma = 0) on the CUDA build, A100 MIG (interactive_gpu): the log line shows
# fma=0; the default gives bitwise the OFK_CS_FMA=0 run and the OFK_CS=fortran run (this unit alone), bar-r8-tree.
#   sbatch -J ofk-circulation-solve-migdefault tests/circulation-solve/gpu_default.sbatch
#SBATCH -p interactive_gpu
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH -t 00:45:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u4a-jobs/%x.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
export LD_LIBRARY_PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build; R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
OUT=$B/u4a-migdefault-$SLURM_JOB_ID; mkdir -p $OUT; cp $R/tests/induced-velocities/run_case.sh $OUT/
export OMP_PROC_BIND=false OMP_NUM_THREADS=4
hostname; nvidia-smi -L; date
c=bar-r8-tree; KC=$B/ov4b-cuda-a100/modules/aerodyn/aerodyn_driver
bash $OUT/run_case.sh $KC $B/iv-cases/$c $OUT/$c.default OFK_MODE=fast
bash $OUT/run_case.sh $KC $B/iv-cases/$c $OUT/$c.F0 OFK_MODE=fast OFK_CS_FMA=0
bash $OUT/run_case.sh $KC $B/iv-cases/$c $OUT/$c.csfortran OFK_MODE=fast OFK_CS=fortran
grep -h "ofk_cs\] circulation" $OUT/$c.default/stdout.log
same() { cmp -s <(tail -n +8 $1) <(tail -n +8 $2) && echo IDENTICAL || echo DIFFER; }
echo "--- $c on the CUDA build: default = OFK_CS_FMA=0: $(same $OUT/$c.default/ad_driver.out $OUT/$c.F0/ad_driver.out); default = OFK_CS=fortran: $(same $OUT/$c.default/ad_driver.out $OUT/$c.csfortran/ad_driver.out)"
date
