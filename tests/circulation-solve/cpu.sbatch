#!/bin/bash
# CPU jobs of unit circulation-solve (med partition, 8 cores). One of:
#   sbatch -J ofk-circulation-solve-iters  cpu.sbatch iters   FMA experiment on bar-r8-tree / bar-r8-direct: per-solve
#        iteration counts of the pin, the pin with FMA (itertrace, unmodified binaries) and the port variants
#        (OFK_CS_ITERS); claim 1 at system level (strict, F0 isolated); R7; R8 of F0 and F1 (r8.py, FMA envelope).
#   sbatch -J ofk-circulation-solve-timing cpu.sbatch timing  P2 + FMA-experiment speed: pin+OpenMP vs port F0 vs F1,
#        1 and 8 threads, best of 3, sequential; OFK_CS_PROFILE per-solve times; wsamp solve buckets (1 run each).
#   sbatch -J ofk-circulation-solve-rtest  cpu.sbatch rtest   R6 (ctest aerodyn_utest|OLAF|5MW_MRSemi on ov4-omp,
#        strict / F0 / F1), r-test OLAF iteration counts (pin, pin+FMA, F0, F1), P4 wall times (ov-off vs ov4-omp).
# Trees: OFK_KF (default ov4-fomp aerodyn_driver), OFK_KO (ov4-omp). Every output dir: opt/build/u4a-<what>-<job>.
#SBATCH -p med
#SBATCH -c 8
#SBATCH -t 08:00:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u4a-jobs/%x.%j.log
set -u
what=$1
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
T=$R/tests/circulation-solve; TI=$R/tests/induced-velocities; L=$R/deps/openfast/reg_tests/lib
PY=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin/python
export PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin:$PATH PYTHONDONTWRITEBYTECODE=1
OUT=$B/u4a-$what-$SLURM_JOB_ID; mkdir -p $OUT
# frozen copies of the helpers this job calls (never edit a script a running job uses)
cp $T/iters_cmp.py $TI/run_case.sh $TI/r8.py $R/scripts/wsamp/wsamp_report.py $OUT/; H=$OUT
hostname; lscpu | grep "Model name"; date
C=$B/iv-cases
FOMP=$B/openfast-pin-fomp/modules/aerodyn/aerodyn_driver   # the pin, Fortran OPENMP=ON (base1 at 1 thread)
FFMA=$B/openfast-pin-ffma/modules/aerodyn/aerodyn_driver   # the pin built with -march=x86-64-v3 (FMA), OPENMP=OFF
KF=${OFK_KF:-$B/ov4-fomp/modules/aerodyn/aerodyn_driver}   # port, Kokkos-OpenMP + Fortran OPENMP=ON (P2 config)
KO=${OFK_KO:-$B/ov4-omp}                                    # port, Kokkos-OpenMP, Fortran OPENMP=OFF (as ov-off)
IT=$B/u4a-jobs/itertrace
same() { cmp -s <(tail -n +8 $1) <(tail -n +8 $2) && echo IDENTICAL || echo DIFFER; }
# breakpoint addresses of an executable: entry of Wings_ComputeCirculationPolarData, the in-loop call of
# CirculationFromPolarData (the last of its two call sites in that routine)
addrs() {
  local s c
  s=$(nm $1 | awk '/ __fvw_wings_MOD_wings_computecirculationpolardata$/{print $1}')
  c=$(objdump -d --no-show-raw-insn $1 | awk '/<__fvw_wings_MOD_wings_computecirculationpolardata>:/{f=1} f&&/^$/{exit} f' \
      | grep "call.*circulationfrompolardata" | tail -1 | awk '{print $1}' | tr -d :)
  echo "$s $c"
}
case $what in
iters)
  export OMP_PROC_BIND=false
  for c in bar-r8-tree bar-r8-direct; do
    echo "=== $c: single-thread runs (7 concurrent, unbound)"
    jobs_=(
      "pin|$IT $OUT/$c.pin.iters $(addrs $FOMP) -- $FOMP|OMP_NUM_THREADS=1"
      "pinfma|$IT $OUT/$c.pinfma.iters $(addrs $FFMA) -- $FFMA|OMP_NUM_THREADS=1"
      "kS|$KF|OMP_NUM_THREADS=1 OFK_MODE=strict OFK_CS_ITERS=$OUT/$c.kS.iters"
      "kI0|$KF|OMP_NUM_THREADS=1 OFK_MODE=strict OFK_CS_MODE=fast OFK_CS_FMA=0 OFK_CS_ITERS=$OUT/$c.kI0.iters"
      "kI1|$KF|OMP_NUM_THREADS=1 OFK_MODE=strict OFK_CS_MODE=fast OFK_CS_FMA=1 OFK_CS_ITERS=$OUT/$c.kI1.iters"
      "kP0|$KF|OMP_NUM_THREADS=1 OFK_MODE=fast OFK_CS_FMA=0 OFK_CS_ITERS=$OUT/$c.kP0.iters"
      "kP1|$KF|OMP_NUM_THREADS=1 OFK_MODE=fast OFK_CS_FMA=1 OFK_CS_ITERS=$OUT/$c.kP1.iters"
    )
    for j in "${jobs_[@]}"; do
      IFS='|' read -r lab exe envs <<< "$j"
      # run_case.sh <exe> <case> <out> VAR=...: the exe may be "itertrace ... -- driver" (word-split on purpose)
      ( cp -r $C/$c $OUT/$c.$lab; ln -sfn $(readlink -f $C/BAR_Baseline) $OUT/BAR_Baseline; cd $OUT/$c.$lab
        t0=$(date +%s.%N); env $envs $exe ad_driver.dvr > stdout.log 2>&1; rc=$?; t1=$(date +%s.%N)
        echo "  $c.$lab rc=$rc wall=$(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.1f", b-a}') s" ) &
    done
    wait
    echo "=== $c: 8-thread and second-process runs (R7), sequential"
    for lab in kP0-8 kP1-8 kP0-b; do
      case $lab in kP0-8) e="OMP_NUM_THREADS=8 OFK_CS_FMA=0";; kP1-8) e="OMP_NUM_THREADS=8 OFK_CS_FMA=1";; kP0-b) e="OMP_NUM_THREADS=1 OFK_CS_FMA=0";; esac
      bash $H/run_case.sh $KF $C/$c $OUT/$c.$lab $e OFK_MODE=fast OFK_CS_ITERS=$OUT/$c.$lab.iters >/dev/null
    done
    b1=$B/iv-r8-1860700/$c.fomp1/ad_driver.out; bf=$B/iv-r8ffma-1860993/$c.ffma/ad_driver.out
    echo "--- traced pin = base1 (job 1860700, untraced): $(same $OUT/$c.pin/ad_driver.out $b1)"
    echo "--- traced pin+FMA = envelope run (job 1860993, untraced): $(same $OUT/$c.pinfma/ad_driver.out $bf)"
    echo "--- claim 1 (system): port strict = base1: $(same $OUT/$c.kS/ad_driver.out $b1)"
    echo "--- claim 1 (system): port F0 isolated (others strict) = base1: $(same $OUT/$c.kI0/ad_driver.out $b1)"
    echo "--- F1 isolated (others strict) = base1: $(same $OUT/$c.kI1/ad_driver.out $b1)"
    echo "--- R7: F0 1 vs 8 threads: $(same $OUT/$c.kP0/ad_driver.out $OUT/$c.kP0-8/ad_driver.out); F0 two processes: $(same $OUT/$c.kP0/ad_driver.out $OUT/$c.kP0-b/ad_driver.out); F1 1 vs 8 threads: $(same $OUT/$c.kP1/ad_driver.out $OUT/$c.kP1-8/ad_driver.out)"
    echo "--- per-solve iteration counts (reference: the pin = base1)"
    $PY $H/iters_cmp.py $OUT/$c.pin.iters pinFMA=$OUT/$c.pinfma.iters strict=$OUT/$c.kS.iters F0iso=$OUT/$c.kI0.iters \
        F1iso=$OUT/$c.kI1.iters F0=$OUT/$c.kP0.iters F1=$OUT/$c.kP1.iters F0-8=$OUT/$c.kP0-8.iters F1-8=$OUT/$c.kP1-8.iters
    echo "--- R8: base1 = pin OPENMP=ON 1 thread (job 1860700); envelope = pin built with FMA (job 1860993)"
    $PY $H/r8.py $L $b1 $bf strict=$OUT/$c.kS/ad_driver.out F0iso=$OUT/$c.kI0/ad_driver.out F1iso=$OUT/$c.kI1/ad_driver.out \
        F0=$OUT/$c.kP0/ad_driver.out F1=$OUT/$c.kP1/ad_driver.out F0-8=$OUT/$c.kP0-8/ad_driver.out F1-8=$OUT/$c.kP1-8/ad_driver.out
  done
  echo "=== itertrace cross-check on a port binary (strict, bar-r8-tree): traced count = OFK_CS_ITERS"
  ( c=bar-r8-tree; cp -r $C/$c $OUT/$c.kS-traced; cd $OUT/$c.kS-traced
    env OMP_NUM_THREADS=1 OFK_MODE=strict OFK_CS_ITERS=$OUT/$c.kS-traced.cs.iters $IT $OUT/$c.kS-traced.iters $(addrs $KF) -- $KF ad_driver.dvr > stdout.log 2>&1
    $PY $H/iters_cmp.py $OUT/$c.kS-traced.iters ofk_cs_counter=$OUT/$c.kS-traced.cs.iters pin=$OUT/$c.pin.iters
    echo "  traced port output = untraced: $(same ad_driver.out $OUT/$c.kS/ad_driver.out)" )
  ;;
timing)
  export OMP_PROC_BIND=spread OMP_PLACES=threads
  for c in bar-r8-tree bar-r8-direct; do
    for t in 8 1; do
      for rep in 1 2 3; do
        bash $H/run_case.sh $FOMP $C/$c $OUT/$c.pin$t.$rep OMP_NUM_THREADS=$t
        bash $H/run_case.sh $KF $C/$c $OUT/$c.F0-$t.$rep OMP_NUM_THREADS=$t OFK_MODE=fast OFK_CS_FMA=0 OFK_CS_PROFILE=1
        bash $H/run_case.sh $KF $C/$c $OUT/$c.F1-$t.$rep OMP_NUM_THREADS=$t OFK_MODE=fast OFK_CS_FMA=1 OFK_CS_PROFILE=1
      done
      for v in pin F0- F1-; do
        echo "  csv,$c,$v$t,$(cat $OUT/$c.$v$t.1/time.txt),$(cat $OUT/$c.$v$t.2/time.txt),$(cat $OUT/$c.$v$t.3/time.txt)"
      done
      for v in F0- F1-; do
        for rep in 1 2 3; do echo "  cs,$c,$v$t.$rep,$(grep -ho 'per solve=[0-9.e+-]*' $OUT/$c.$v$t.$rep/stdout.log | head -1),$(grep -ho 'build=[0-9.]* replay=[0-9.]*' $OUT/$c.$v$t.$rep/stdout.log | head -1)"; done
      done
    done
    echo "=== $c: wsamp solve buckets (one sampled run each, 8 and 1 threads)"
    for t in 8 1; do
      for v in pin F0 F1; do
        case $v in pin) exe=$FOMP; e="";; F0) exe=$KF; e="OFK_MODE=fast OFK_CS_FMA=0";; F1) exe=$KF; e="OFK_MODE=fast OFK_CS_FMA=1";; esac
        o=$OUT/$c.$v$t.prof
        bash $H/run_case.sh $exe $C/$c $o OMP_NUM_THREADS=$t $e LD_PRELOAD=$B/u4-prof/libwsamp.so WSAMP_OUT=$o/prof WSAMP_HZ=500 >/dev/null
        $PY $H/wsamp_report.py $o/prof 500 --buckets --csv $c.$v$t | grep "^csv"
      done
    done
  done
  ;;
rtest)
  export OMP_PROC_BIND=false OMP_NUM_THREADS=1
  echo "=== R6: ctest aerodyn_utest|OLAF|5MW_MRSemi on ov4-omp (OMP_NUM_THREADS=1): strict, fast F0, fast F1"
  for v in strict F0 F1; do
    case $v in strict) e="OFK_MODE=strict";; F0) e="OFK_MODE=fast OFK_CS_FMA=0";; F1) e="OFK_MODE=fast OFK_CS_FMA=1";; esac
    ( cd $KO && env $e ctest -R 'aerodyn_utest|OLAF|5MW_MRSemi_DLL_WSt_WavesIrr' --timeout 3000 > $OUT/ctest.$v.log 2>&1 )
    echo "  $v: $(grep -E 'tests passed' $OUT/ctest.$v.log)"; grep -E "\*\*\*Failed|Not Run|Timeout" $OUT/ctest.$v.log | sed 's/^/    /'
  done
  echo "=== r-test OLAF cases: iteration counts and P4 wall times (ov-off = the pin path, OPENMP=OFF, vs ov4-omp)"
  RT=$R/deps/openfast/reg_tests/r-test
  # cases as ctest prepared them in the ov4-omp tree (the R6 run above; glue cases need its DISCON builds)
  mkrun() {  # <glue|ad> <case> <dest>: a copy of the case with its siblings linked
    local kind=$1 c=$2 d=$3 src; [ $kind = ad ] && src=$KO/reg_tests/modules/aerodyn || src=$KO/reg_tests/glue-codes/openfast
    mkdir -p $d; for s in $src/*; do [ "$(basename $s)" = $c ] || ln -sfn $s $d/; done
    rm -rf $d/$c; cp -r $src/$c $d/$c
  }
  run1() {  # <kind> <case> <label> <command words...> (the input file is appended) -> wall seconds
    local kind=$1 c=$2 lab=$3; shift 3; local d=$OUT/rt/$lab
    mkrun $kind $c $d; cd $d/$c
    local in=ad_driver.dvr; [ $kind = glue ] && in=$c.fst
    local t0=$(date +%s.%N); "$@" $in > stdout.log 2>&1; local rc=$?; local t1=$(date +%s.%N)
    cd - >/dev/null; awk -v a=$t0 -v b=$t1 -v rc=$rc 'BEGIN{printf "%.2f%s", b-a, (rc?"(rc!)":"")}'
  }
  OFFA=$B/ov-off/modules/aerodyn/aerodyn_driver; OFFG=$B/ov-off/glue-codes/openfast/openfast
  KOA=$KO/modules/aerodyn/aerodyn_driver; KOG=$KO/glue-codes/openfast/openfast
  for cc in ad:ad_BAR_OLAF ad:ad_BAR_OLAF_RegFunctionPart2 ad:ad_B1n2_OLAF ad:ad_EllipticalWingInf_OLAF ad:ad_HelicalWakeInf_OLAF \
            ad:ad_Kite_OLAF ad:ad_QuadRotor_OLAF ad:ad_VerticalAxis_OLAF ad:ad_Sphere_OLAF \
            glue:HelicalWake_OLAF glue:EllipticalWing_OLAF glue:5MW_Land_QuadRotor_DLL_WSt_OLAF; do
    kind=${cc%%:*}; c=${cc#*:}
    if [ $kind = ad ]; then OE=$OFFA; KE=$KOA; else OE=$OFFG; KE=$KOG; fi
    a=$(run1 $kind $c pin $IT $OUT/$c.pin.iters $(addrs $OE) -- $OE)
    f="-"; [ $kind = ad ] && f=$(run1 $kind $c pinfma $IT $OUT/$c.pinfma.iters $(addrs $FFMA) -- $FFMA)
    s=$(run1 $kind $c F0iso env OFK_MODE=strict OFK_CS_MODE=fast OFK_CS_FMA=0 OFK_CS_ITERS=$OUT/$c.F0iso.iters $KE)
    g=$(run1 $kind $c F1iso env OFK_MODE=strict OFK_CS_MODE=fast OFK_CS_FMA=1 OFK_CS_ITERS=$OUT/$c.F1iso.iters $KE)
    h=$(run1 $kind $c F1 env OFK_MODE=fast OFK_CS_FMA=1 OFK_CS_ITERS=$OUT/$c.F1.iters $KE)
    echo "    .outb channel data, F0iso and F1iso vs the traced pin:"; $PY $R/tests/biotsavart-direct/cmp_outb.py $L $OUT/rt/pin/$c $OUT/rt/F0iso/$c $OUT/rt/F1iso/$c 2>&1 | tail -4 | sed 's/^/      /'
    echo "--- $c (traced pin ${a}s, pin+FMA ${f}s)"
    if [ -s $OUT/$c.pin.iters ]; then
      $PY $H/iters_cmp.py $OUT/$c.pin.iters pinFMA=$OUT/$c.pinfma.iters F0iso=$OUT/$c.F0iso.iters F1iso=$OUT/$c.F1iso.iters F1=$OUT/$c.F1.iters | sed 's/^/    /'
    else echo "    no circulation solve (prescribed circulation or none)"; fi
    # P4: wall time, best of 3, untraced: pin path vs port (fast F0, default threshold), 1 and 8 threads
    for t in 1 8; do
      bo=1e9; bk=1e9
      for rep in 1 2 3; do
        x=$(run1 $kind $c off$t env OMP_NUM_THREADS=$t $OE); bo=$(awk -v a=$bo -v b=${x%%(*} 'BEGIN{print (b<a)?b:a}')
        y=$(run1 $kind $c kk$t env OMP_NUM_THREADS=$t OFK_MODE=fast $KE); bk=$(awk -v a=$bk -v b=${y%%(*} 'BEGIN{print (b<a)?b:a}')
      done
      echo "    P4 $c threads=$t: pin path $bo s, port fast $bk s, ratio $(awk -v a=$bk -v b=$bo 'BEGIN{printf "%.3f", a/b}')"
    done
  done
  ;;
esac
date
