#!/bin/bash
# "No change in behaviour" evidence for the overlay layout (notes/overlay.md), CPU part (med, 8 cores).
#   sbatch -J ofk-infra-l0  tests/overlay/cpu.sbatch l0    L0 of units 1 and 2 on the overlay-built test trees
#                                                          (ov-tests-bsd-omp, ov-tests-iv-omp): the per-branch
#                                                          lines of unit 1's job 1859729, then full ctest
#   sbatch -J ofk-infra-sys tests/overlay/cpu.sbatch sys   ctest 'aerodyn_utest|OLAF' on ov-off, ov-omp strict, ov-omp
#                                                          fast (+ fast with OFK_TC=fortran); .outb channels vs OFF;
#                                                          BAR full-precision cases on ov-omp vs the fork-omp outputs
#SBATCH -p med
#SBATCH -c 8
#SBATCH -t 03:00:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/ov-cpu.%x.%j.log
set -u
what=$1
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
PY=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin/python
export PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin:$PATH
L=$R/deps/openfast/reg_tests/lib
export PYTHONDONTWRITEBYTECODE=1   # the r-test Python lives in deps/openfast: leave no __pycache__ there
OUT=$B/ov-$what-$SLURM_JOB_ID; mkdir -p $OUT
hostname; lscpu | grep "Model name"; date
case $what in
l0)
  export OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=threads
  cd $B/ov-tests-bsd-omp
  for seed in 1 7; do   # the same commands and filters as unit 1's cpu.sbatch (job 1859729)
    echo "=== seed $seed: chain claim 1, Fortran (OPENMP=OFF) vs strict, rtol 0 (bitwise)"
    for exe in l0_ui_seg l0_ui_part_nograd; do OFK_MODE=strict OFK_L0_SEED=$seed OFK_L0_RTOL=0 ./$exe | grep -v "^\[ofk\]" | tail -7; done
    echo "=== seed $seed: chain claim 3, strict vs fast (err / sum|c|), declared R1 bound 1e-12"
    for exe in l0_ui_seg l0_ui_part_nograd; do OFK_L0_STRICTFAST=1 OFK_L0_SEED=$seed OFK_L0_RTOL=1e-12 ./$exe | grep -v "^FAIL\|^\[ofk\]" | tail -7; done
    echo "=== seed $seed: envelope, Fortran built with FMA (-march=x86-64-v3) vs Fortran without (= strict)"
    for exe in l0_ui_seg_ffma l0_ui_part_nograd_ffma; do OFK_MODE=strict OFK_L0_SEED=$seed OFK_L0_RTOL=1 ./$exe | grep -v "^\[ofk\]" | tail -7; done
    echo "=== seed $seed: informational, Fortran built with -fopenmp vs strict"
    for exe in l0_ui_seg_fomp l0_ui_part_nograd_fomp; do OFK_MODE=strict OFK_L0_SEED=$seed OFK_L0_RTOL=1 ./$exe | grep -v "^\[ofk\]" | tail -1; done
  done
  for seed in 1 7; do
    echo "=== unit 1 ctest, OFK_L0_SEED=$seed (claim 3 bound = max(1e-12, 1.1 x Fortran FMA spread) per branch)"
    if [ $seed = 1 ]; then sel=(); else sel=(-R 'ffma|strictfast|_fast$|_strict$'); fi
    OFK_L0_SEED=$seed ctest "${sel[@]}" --timeout 900 > $OUT/bsd.seed$seed.ctest.log 2>&1
    grep -E "Test +#|tests passed" $OUT/bsd.seed$seed.ctest.log
    grep -hE "^  reg .* bound=|^l0: per-branch" Testing/Temporary/LastTest.log
    cp Testing/Temporary/LastTest.log $OUT/bsd.seed$seed.LastTest.log
  done
  echo "=== unit 2 ctest (ov-tests-iv-omp)"
  cd $B/ov-tests-iv-omp
  ctest --timeout 900 > $OUT/iv.ctest.log 2>&1
  grep -E "Test +#|tests passed" $OUT/iv.ctest.log
  cp Testing/Temporary/LastTest.log $OUT/iv.LastTest.log
  grep -hE "^l0_iv: .*passed|^  (Basic|Part|TreePart|A7|pack) |P6 \(L0|R7 |FAIL" $OUT/iv.LastTest.log | head -80
  ;;
sys)
  export OMP_NUM_THREADS=1 OMP_PROC_BIND=false
  RX='aerodyn_utest|OLAF'
  run() {  # <tree> <label> [env...]
    local tree=$1 label=$2; shift 2
    touch $OUT/stamp.$label
    ( cd $B/$tree && env "$@" ctest -R "$RX" -j 8 --timeout 3000 > $OUT/$label.ctest.log 2>&1 )
    echo "=== $label ($tree; $*)"; grep -E "Test +#|tests passed" $OUT/$label.ctest.log
    mkdir -p $OUT/$label; ( cd $B/$tree/reg_tests && find . -name "*.outb" -newer $OUT/stamp.$label -print0 | xargs -0 -I{} cp --parents {} $OUT/$label/ )
  }
  run ov-off off                OFK_MODE=fast
  run ov-omp omp-strict         OFK_MODE=strict
  run ov-omp omp-fast           OFK_MODE=fast
  run ov-omp omp-fast-tcfortran OFK_MODE=fast OFK_TC=fortran
  echo "=== .outb channel data vs ov-off"
  $PY $R/tests/biotsavart-direct/cmp_outb.py $L $OUT/off $OUT/omp-strict $OUT/omp-fast $OUT/omp-fast-tcfortran
  echo "=== .outb channel data: ov-* vs the fork builds (job 1860600: iv-r6-1860600/{off,omp-strict,omp-fast})"
  for l in off omp-strict omp-fast; do
    echo "--- fork $l vs ov $l (and ov omp-fast-tcfortran)"
    $PY $R/tests/biotsavart-direct/cmp_outb.py $L $B/iv-r6-1860600/$l $OUT/$l $OUT/omp-fast-tcfortran
  done
  echo "=== BAR full-precision cases on ov-omp (data lines vs fork-omp, job 1860700: iv-r8-1860700)"
  T=$R/tests/induced-velocities; C=$B/iv-cases; KK=$B/ov-omp/modules/aerodyn/aerodyn_driver; OFF=$B/ov-off/modules/aerodyn/aerodyn_driver
  for c in bar-r8-direct bar-r8-tree; do
    bash $T/run_case.sh $KK  $C/$c $OUT/$c.kk-strict OMP_NUM_THREADS=2 OFK_MODE=strict &
    bash $T/run_case.sh $KK  $C/$c $OUT/$c.kk-fast   OMP_NUM_THREADS=2 OFK_MODE=fast &
    bash $T/run_case.sh $OFF $C/$c $OUT/$c.off       OMP_NUM_THREADS=1 &
    bash $T/run_case.sh $KK  $C/$c $OUT/$c.kk-fast-tcfortran OMP_NUM_THREADS=2 OFK_MODE=fast OFK_TC=fortran &
    wait
    for m in kk-strict kk-fast off; do
      a=$B/iv-r8-1860700/$c.$m/ad_driver.out; b=$OUT/$c.$m/ad_driver.out
      echo "--- $c $m: ov vs fork data lines (from line 8): $(cmp -s <(tail -n +8 $a) <(tail -n +8 $b) && echo IDENTICAL || echo DIFFER) ($(tail -n +8 $b | wc -l) lines)"
    done
    a=$B/iv-r8-1860700/$c.kk-fast/ad_driver.out; b=$OUT/$c.kk-fast-tcfortran/ad_driver.out
    echo "--- $c kk-fast with OFK_TC=fortran: ov vs fork kk-fast: $(cmp -s <(tail -n +8 $a) <(tail -n +8 $b) && echo IDENTICAL || echo DIFFER)"
    echo "--- $c ov kk-strict vs ov off: $(cmp -s <(tail -n +8 $OUT/$c.off/ad_driver.out) <(tail -n +8 $OUT/$c.kk-strict/ad_driver.out) && echo IDENTICAL || echo DIFFER)"
  done
  ;;
esac
date
