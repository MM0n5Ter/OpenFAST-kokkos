#!/bin/bash
# CPU jobs of unit treecode (med partition, 8 cores; never two concurrent spread-bound runs). One of:
#   sbatch -J ofk-treecode-sweep cpu.sbatch sweep   P1/P5: perf_tc N=1e3..1e6, Fortran+OpenMP vs Kokkos-OpenMP fast/strict
#   sbatch -J ofk-treecode-l0    cpu.sbatch l0      L0 of units 1, 2 (with the tree ported) and 3 on the CPU, R7 threads
#   sbatch -J ofk-treecode-sys   cpu.sbatch sys     R8/R7 on bar-r8-tree (full-precision text) with ov-fomp, P2 (ov-fomp vs
#                                                  pin+OpenMP, 8 threads, bar-r8-tree and bar-perf-tree, OFK_TC=fortran A/B), P6.
#                                                  OFK_KF=<aerodyn_driver> overrides the Kokkos build (default ov-fomp).
#SBATCH -p med
#SBATCH -c 8
#SBATCH -t 06:00:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u3-cpu.%x.%j.log
set -u
what=$1
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
T=$R/tests/treecode; TI=$R/tests/induced-velocities; L=$R/deps/openfast/reg_tests/lib; export PYTHONDONTWRITEBYTECODE=1
PY=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin/python
export PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin:$PATH
OUT=$B/u3-$what-$SLURM_JOB_ID; mkdir -p $OUT
export OMP_PROC_BIND=spread OMP_PLACES=threads
hostname; lscpu | grep "Model name"; date
C=$B/iv-cases; C3=$B/u3-cases
FOMP=$B/openfast-pin-fomp/modules/aerodyn/aerodyn_driver   # the pin, Fortran OPENMP=ON
KF=${OFK_KF:-$B/ov-fomp/modules/aerodyn/aerodyn_driver}           # overlay, Kokkos-OpenMP + Fortran OPENMP=ON (like for like)
case $what in
sweep)
  cp $B/u3-omp/perf_tc $OUT/
  for n in 1000 3000 10000 30000 100000 300000 1000000; do
    OMP_NUM_THREADS=8 $OUT/perf_tc $n 3 1 1.5 fortran | grep csv | tee -a $OUT/sweep.csv
    for m in fast strict; do
      OMP_NUM_THREADS=8 OFK_MODE=$m $OUT/perf_tc $n 3 1 1.5 port | grep csv | tee -a $OUT/sweep.csv
      OMP_NUM_THREADS=8 OFK_MODE=$m OFK_IV_PROFILE=1 $OUT/perf_tc $n 3 1 1.5 port 2>&1 | grep "report" | sed "s/^/  $m n=$n /"
    done
  done
  echo "=== accuracy check of the sweep inputs: port fast vs Fortran tree (same tree, rounding only)"
  for n in 10000 100000; do OMP_NUM_THREADS=8 OFK_MODE=fast $OUT/perf_tc $n 1 1 1.5 both | grep "max |"; done
  ;;
l0)
  export OMP_NUM_THREADS=8
  ( cd $B/u3-omp && ctest --timeout 3000 > $OUT/tc.ctest.log 2>&1 ); echo "=== unit 3 l0_tc"; grep -E "Test +#|tests passed" $OUT/tc.ctest.log
  grep -hE "^l0_tc|claim|R2:|P6|FAIL|R7" $B/u3-omp/Testing/Temporary/LastTest.log | head -40
  ( cd $B/u3-iv-omp && ctest --timeout 3000 > $OUT/iv.ctest.log 2>&1 ); echo "=== unit 2 l0_iv (tree ported)"; grep -E "Test +#|tests passed" $OUT/iv.ctest.log
  grep -hE "^l0_iv: .*passed|^  TreePart|P6 \(L0" $B/u3-iv-omp/Testing/Temporary/LastTest.log | head -40
  ( cd $B/u3-bsd-omp && ctest --timeout 3000 > $OUT/bsd.ctest.log 2>&1 ); echo "=== unit 1 l0 (unchanged code, rebuilt from the src/ overlay)"; grep -E "Test +#|tests passed" $OUT/bsd.ctest.log
  ;;
sys)
  # R6 and claim 1 at system level are done by job 1863342 (notes/overlay.md 4.3-4.4): ov-off vs ov-omp strict
  # (Kokkos tree) identical from line 8. This block adds R8/R7 on the P2 build, P2 and P6.
  echo "=== R8 + R7 (system): bar-r8-tree on ov-fomp, full-precision text output"
  c=bar-r8-tree
  bash $TI/run_case.sh $KF $C/$c $OUT/$c.kf-fast1 OMP_NUM_THREADS=1 OMP_PROC_BIND=false OFK_MODE=fast &
  bash $TI/run_case.sh $KF $C/$c $OUT/$c.kf-strict1 OMP_NUM_THREADS=1 OMP_PROC_BIND=false OFK_MODE=strict &
  wait
  bash $TI/run_case.sh $KF $C/$c $OUT/$c.kf-fast8 OMP_NUM_THREADS=8 OFK_MODE=fast
  same() { cmp -s <(tail -n +8 $1) <(tail -n +8 $2) && echo IDENTICAL || echo DIFFER; }
  echo "--- R7 (system): fast 1 vs 8 threads, from line 8: $(same $OUT/$c.kf-fast1/ad_driver.out $OUT/$c.kf-fast8/ad_driver.out)"
  echo "--- strict on ov-fomp (Fortran OPENMP=ON elsewhere) vs base1 (pin OPENMP=ON, 1 thread): $(same $OUT/$c.kf-strict1/ad_driver.out $B/iv-r8-1860700/$c.fomp1/ad_driver.out)"
  echo "--- R8: base1 = pin OPENMP=ON 1 thread (job 1860700); envelope = pin built with FMA (job 1860993)"
  $PY $TI/r8.py $L $B/iv-r8-1860700/$c.fomp1/ad_driver.out $B/iv-r8ffma-1860993/$c.ffma/ad_driver.out \
      kf-strict1=$OUT/$c.kf-strict1/ad_driver.out kf-fast1=$OUT/$c.kf-fast1/ad_driver.out kf-fast8=$OUT/$c.kf-fast8/ad_driver.out
  echo "=== P2 (8 threads, sequential, 2 repetitions interleaved): pin+OpenMP vs ov-fomp"
  for c in bar-r8-tree bar-perf-tree; do
    d=$C; [ $c = bar-perf-tree ] && d=$C3
    for rep in 1 2; do
      bash $TI/run_case.sh $FOMP $d/$c $OUT/$c.p2-fomp8.$rep OMP_NUM_THREADS=8
      bash $TI/run_case.sh $KF   $d/$c $OUT/$c.p2-kf-fast8.$rep OMP_NUM_THREADS=8 OFK_MODE=fast
      bash $TI/run_case.sh $KF   $d/$c $OUT/$c.p2-kf-tcfortran8.$rep OMP_NUM_THREADS=8 OFK_MODE=fast OFK_TC=fortran
      bash $TI/run_case.sh $KF   $d/$c $OUT/$c.p2-kf-strict8.$rep OMP_NUM_THREADS=8 OFK_MODE=strict
    done
    bash $TI/run_case.sh $KF   $d/$c $OUT/$c.p2-kf-fast8.prof OMP_NUM_THREADS=8 OFK_MODE=fast OFK_IV_PROFILE=1
    for r in p2-fomp8 p2-kf-fast8 p2-kf-tcfortran8 p2-kf-strict8; do
      echo "    csv,$c,$r,$(cat $OUT/$c.$r.1/time.txt),$(cat $OUT/$c.$r.2/time.txt)"
    done
    grep -h "ofk_iv\]" $OUT/$c.p2-kf-fast8.prof/stdout.log | tail -3
  done
  echo "=== P6 (CPU, full run): Kokkos memory events"
  MT=$B/ov-tests-iv-omp/libofk_memtool.so
  bash $TI/run_case.sh $KF $C/bar-r8-tree $OUT/bar-r8-tree.mem OMP_NUM_THREADS=8 OFK_MODE=fast KOKKOS_TOOLS_LIBS=$MT OFK_MEMTOOL_OUT=$OUT/memtool.tree.txt
  tail -1 $OUT/memtool.tree.txt
  ;;
esac
date
