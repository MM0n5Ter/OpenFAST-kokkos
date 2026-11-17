#!/bin/bash
# Unit-4 analysis: exact call counts (gprof, -pg build u4-fomp-pg, 1 thread) on the BAR cases.
#SBATCH -p med
#SBATCH -c 2
#SBATCH -t 01:30:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u4-prof/%x.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 openblas/0.3.28 >/dev/null 2>&1
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
TI=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos/tests/induced-velocities
X=$B/u4-fomp-pg/modules/aerodyn/aerodyn_driver; OUT=$B/u4-prof/gprof-$SLURM_JOB_ID; mkdir -p $OUT
hostname; date
for c in bar-r8-tree bar-r8-direct bar-perf-tree; do
  d=$B/iv-cases/$c; [ $c = bar-perf-tree ] && d=$B/u3-cases/$c
  bash $TI/run_case.sh $X $d $OUT/$c OMP_NUM_THREADS=1 OFK_MODE=fast
  ( cd $OUT/$c && gprof -b $X gmon.out > gprof.txt 2>&1 && gprof -b -p $X gmon.out > flat.txt 2>&1 )
  grep -E "ui_quad_n1|circulationfrompolardata|wings_computecirculationpolardata|liftinglineinducedvelocities|wakeinducedvelocities|ifw_flowfield_getvelacc|ad_calcwind |twrinflarray|setinputsforfvw|fvw_updatestates|fvw_calcoutput|propagatewake|fvw_euler1|fvw_calccontstatederiv|fvw_copycontstate|ui_seg_11" $OUT/$c/flat.txt
done
date
