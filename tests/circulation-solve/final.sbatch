#!/bin/bash
# Final-source check of unit circulation-solve: ov4b-fomp (nCP_tot = 0 fix, default OFK_CS_FMA=1) must give
# bitwise the outputs the FMA experiment measured with ov4-fomp (job u4a-iters): default = F1, OFK_CS_FMA=0 = F0.
#   sbatch -J ofk-circulation-solve-final final.sbatch <iters job output dir>
#SBATCH -p med
#SBATCH -c 8
#SBATCH -t 00:30:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u4a-jobs/%x.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 openblas/0.3.28 >/dev/null 2>&1
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build; I=$1
OUT=$B/u4a-final-$SLURM_JOB_ID; mkdir -p $OUT; cp /weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos/tests/induced-velocities/run_case.sh $OUT/
KF=$B/ov4b-fomp/modules/aerodyn/aerodyn_driver
export OMP_PROC_BIND=spread OMP_PLACES=threads OMP_NUM_THREADS=8
same() { cmp -s <(tail -n +8 $1) <(tail -n +8 $2) && echo IDENTICAL || echo DIFFER; }
for c in bar-r8-tree bar-r8-direct; do
  bash $OUT/run_case.sh $KF $B/iv-cases/$c $OUT/$c.default OFK_MODE=fast
  bash $OUT/run_case.sh $KF $B/iv-cases/$c $OUT/$c.F0 OFK_MODE=fast OFK_CS_FMA=0
  bash $OUT/run_case.sh $KF $B/iv-cases/$c $OUT/$c.strict OFK_MODE=strict
  grep -h "ofk_cs\] circulation" $OUT/$c.default/stdout.log
  echo "--- $c: default = experiment F1: $(same $OUT/$c.default/ad_driver.out $I/$c.kP1-8/ad_driver.out); F0 = experiment F0: $(same $OUT/$c.F0/ad_driver.out $I/$c.kP0-8/ad_driver.out); strict = base1: $(same $OUT/$c.strict/ad_driver.out $B/iv-r8-1860700/$c.fomp1/ad_driver.out)"
done
date
