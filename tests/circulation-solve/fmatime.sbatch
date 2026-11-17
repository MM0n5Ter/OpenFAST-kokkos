#!/bin/bash
# FMA experiment of unit circulation-solve, the speed half of the decision rule: F0 (no FMA) vs F1 (FMA),
# same binary (ov4-fomp, OFK_CS_FMA=0|1), interleaved A/B so drifts of the node hit both alike.
#   bar-r8-tree: 5 interleaved pairs at 8 threads and at 1 thread; bar-r8-direct: 3 pairs at 8 threads.
#   OFK_CS_PROFILE gives the solve kernel time (cache build + replays) per run; wall time from run_case.sh.
#   perf_cs: replay/build time per call at nCP = 87 (BAR) and 150, 10 repetitions, 1 and 8 threads.
#   sbatch -J ofk-circulation-solve-fmatime fmatime.sbatch
#SBATCH -p med
#SBATCH -c 8
#SBATCH -t 03:00:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u4a-jobs/%x.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 openblas/0.3.28 >/dev/null 2>&1
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
OUT=$B/u4a-fmatime-$SLURM_JOB_ID; mkdir -p $OUT; cp $R/tests/induced-velocities/run_case.sh $OUT/
KF=$B/ov4-fomp/modules/aerodyn/aerodyn_driver; C=$B/iv-cases
export OMP_PROC_BIND=spread OMP_PLACES=threads
hostname; lscpu | grep "Model name"; date
pair() {  # case threads rep
  for v in 0 1; do
    bash $OUT/run_case.sh $KF $C/$1 $OUT/$1.F$v-$2.$3 OMP_NUM_THREADS=$2 OFK_MODE=fast OFK_CS_FMA=$v OFK_CS_PROFILE=1 > /dev/null
    echo "  run,$1,F$v,$2,$3,$(cat $OUT/$1.F$v-$2.$3/time.txt),$(grep -ho 'per solve=[0-9.e+-]*' $OUT/$1.F$v-$2.$3/stdout.log | head -1 | cut -d= -f2),$(grep -ho 'build=[0-9.]* replay=[0-9.]*' $OUT/$1.F$v-$2.$3/stdout.log | head -1)"
  done
}
for rep in 1 2 3 4 5; do pair bar-r8-tree 8 $rep; done
for rep in 1 2 3 4 5; do pair bar-r8-tree 1 $rep; done
for rep in 1 2 3; do pair bar-r8-direct 8 $rep; done
echo "=== P6 (full run): Kokkos memory events, bar-r8-tree, 8 threads, F0 (unit 2 memtool)"
bash $OUT/run_case.sh $KF $C/bar-r8-tree $OUT/bar-r8-tree.mem OMP_NUM_THREADS=8 OFK_MODE=fast KOKKOS_TOOLS_LIBS=$B/ov-tests-iv-omp/libofk_memtool.so OFK_MEMTOOL_OUT=$OUT/memtool.txt > /dev/null
tail -1 $OUT/memtool.txt; echo "  ofk_cs allocations (label, bytes, ofk_iv evaluations before):"; grep "^alloc.*ofk_cs" $OUT/memtool.txt | sed 's/^/    /'
echo "=== perf_cs (u4a-omp, same kernels), nCP 87 and 150, 10 reps (best), 8 and 1 threads"
for t in 8 1; do OMP_NUM_THREADS=$t $B/u4a-jobs/perf_cs.v2 10 3 87 150 | grep csv; done
date
