#!/bin/bash
# Unit-4 analysis: wall-clock stack profiles (main thread, LD_PRELOAD libwsamp.so) of the BAR OLAF cases.
#   sbatch -J ofk-unit4-prof prof.sbatch [set]
# set=all (default): ov-fomp (P2 config: Kokkos-OpenMP + Fortran OpenMP), openfast-pin-fomp (baseline),
#   ov-omp (Kokkos-OpenMP, Fortran OPENMP=OFF = the host side of every CUDA build) -- 8 threads, sequential.
# Sampler: gcc -O2 -fPIC -shared -o $P/libwsamp.so scripts/wsamp/wsamp.c -lrt ; report: scripts/wsamp/wsamp_report.py <prefix> <hz> [--buckets] [--csv label]
#SBATCH -p med
#SBATCH -c 8
#SBATCH -t 02:30:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u4-prof/%x.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 openblas/0.3.28 >/dev/null 2>&1
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
TI=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos/tests/induced-velocities
P=$B/u4-prof; OUT=$P/run-$SLURM_JOB_ID; mkdir -p $OUT
export OMP_PROC_BIND=spread OMP_PLACES=threads
hostname; lscpu | grep "Model name"; date
C=$B/iv-cases; C3=$B/u3-cases
declare -A EXE=([kf]=$B/ov-fomp/modules/aerodyn/aerodyn_driver [pin]=$B/openfast-pin-fomp/modules/aerodyn/aerodyn_driver [ko]=$B/ov-omp/modules/aerodyn/aerodyn_driver)
cdir() { [ $1 = bar-perf-tree ] && echo $C3/$1 || echo $C/$1; }
run() { # build case threads extra...
  local b=$1 c=$2 t=$3; shift 3; local o=$OUT/$c.$b$t
  bash $TI/run_case.sh ${EXE[$b]} $(cdir $c) $o OMP_NUM_THREADS=$t OFK_MODE=fast "$@" \
       LD_PRELOAD=$P/libwsamp.so WSAMP_OUT=$o/prof WSAMP_HZ=500 OFK_IV_PROFILE=1
  grep -ao "Total Real Time.*" $o/stdout.log; grep -ao "\[ofk_iv\].*" $o/stdout.log | cut -c1-420
}
for c in bar-r8-tree bar-r8-direct bar-perf-tree bar-perf-part bar-perf-direct; do
  run kf $c 8; run ko $c 8; run pin $c 8
done
echo "=== sampler overhead check: ov-fomp bar-r8-tree unsampled (job 1887061: 40.9/41.0 s)"
bash $TI/run_case.sh ${EXE[kf]} $C/bar-r8-tree $OUT/bar-r8-tree.kf8-nosamp OMP_NUM_THREADS=8 OFK_MODE=fast
for d in $OUT/*/; do [ -f $d/prof.stacks ] && python3 $P/wsamp_report.py $d/prof 500 --buckets > $d/report.txt 2>&1; done
date
