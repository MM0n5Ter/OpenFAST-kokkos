#!/bin/bash
# Check of the default change (kDefaultFma = 0, decisions.md #15) on CPU:
#  1. BAR whole runs on ov4b-fomp (8 threads): the log line shows fma=0; the default (no OFK_CS_FMA) gives bitwise
#     the FMA experiment's F0 outputs, OFK_CS_FMA=1 its F1 outputs (opt-in still works), strict = base1;
#     per-solve iteration counts of the default = the pin's (via the experiment's traced counts);
#  2. R6 on ov4b-omp (1 thread): ctest aerodyn_utest|OLAF|5MW_MRSemi, default fast mode and strict.
#   sbatch -J ofk-circulation-solve-default tests/circulation-solve/default.sbatch
#SBATCH -p med
#SBATCH -c 8
#SBATCH -t 01:30:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u4a-jobs/%x.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 openblas/0.3.28 >/dev/null 2>&1
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build; R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
I=$B/u4a-iters-1925585                       # the FMA experiment's outputs and traced pin counts
OUT=$B/u4a-default-$SLURM_JOB_ID; mkdir -p $OUT; cp $R/tests/induced-velocities/run_case.sh $R/tests/circulation-solve/iters_cmp.py $OUT/
KF=$B/ov4b-fomp/modules/aerodyn/aerodyn_driver
export OMP_PROC_BIND=spread OMP_PLACES=threads OMP_NUM_THREADS=8
hostname; date
same() { cmp -s <(tail -n +8 $1) <(tail -n +8 $2) && echo IDENTICAL || echo DIFFER; }
echo "=== 1. BAR whole runs, ov4b-fomp, 8 threads"
for c in bar-r8-tree bar-r8-direct; do
  bash $OUT/run_case.sh $KF $B/iv-cases/$c $OUT/$c.default OFK_MODE=fast OFK_CS_ITERS=$OUT/$c.default.iters
  bash $OUT/run_case.sh $KF $B/iv-cases/$c $OUT/$c.F1 OFK_MODE=fast OFK_CS_FMA=1
  bash $OUT/run_case.sh $KF $B/iv-cases/$c $OUT/$c.strict OFK_MODE=strict
  grep -h "ofk_cs\] circulation" $OUT/$c.default/stdout.log
  echo "--- $c: default = experiment F0: $(same $OUT/$c.default/ad_driver.out $I/$c.kP0-8/ad_driver.out); OFK_CS_FMA=1 = experiment F1: $(same $OUT/$c.F1/ad_driver.out $I/$c.kP1-8/ad_driver.out); strict = base1: $(same $OUT/$c.strict/ad_driver.out $B/iv-r8-1860700/$c.fomp1/ad_driver.out)"
  echo "--- $c: default iteration counts vs the pin's (traced):"; python3 $OUT/iters_cmp.py $I/$c.pin.iters default=$OUT/$c.default.iters 2>&1 | sed 's/^/    /'
done
echo "=== 2. R6 on ov4b-omp, 1 thread: default fast mode and strict"
KO=$B/ov4b-omp
for v in default strict; do
  case $v in default) e="OFK_MODE=fast";; strict) e="OFK_MODE=strict";; esac
  ( cd $KO && env OMP_NUM_THREADS=1 $e ctest -R 'aerodyn_utest|OLAF|5MW_MRSemi_DLL_WSt_WavesIrr' --timeout 3000 > $OUT/ctest.$v.log 2>&1 )
  echo "  $v: $(grep -E 'tests passed' $OUT/ctest.$v.log)"; grep -E "\*\*\*Failed|Not Run|Timeout" $OUT/ctest.$v.log | sed 's/^/    /'
done
date
