#!/bin/bash
# CPU jobs of unit induced-velocities (med partition, 8 cores). One of:
#   sbatch -J ofk-induced-velocities-r6   cpu.sbatch r6    R6 + L0 ctest: fork-off / fork-omp fast / strict, .outb vs OFF
#   sbatch -J ofk-induced-velocities-r8   cpu.sbatch r8    R8: Fortran (pin, OPENMP=ON) OMP 1 vs 8, fork-off, fork-omp fast/strict
#   sbatch -J ofk-induced-velocities-perf cpu.sbatch perf  P1/P2/P3: Fortran+OpenMP vs fork-omp fast/strict, profile split
#   sbatch -J ofk-induced-velocities-r8ffma cpu.sbatch r8ffma <r8 output dir>
#                                                    R8 with a second envelope: the pin built with FMA (-march=x86-64-v3)
#   sbatch -J ofk-induced-velocities-perf2 cpu.sbatch perf2  the fork with Kokkos-OpenMP AND Fortran OpenMP
#                                                    (fork-omp-fomp) vs the pin with OpenMP; P6 memory events
#SBATCH -p med
#SBATCH -c 8
#SBATCH -t 04:00:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/iv-cpu.%x.%j.log
set -u
what=$1
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
T=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos/tests/induced-velocities
L=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos/deps/openfast/reg_tests/lib
PY=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin/python
export PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin:$PATH
OUT=$B/iv-$what-$SLURM_JOB_ID; mkdir -p $OUT
export OMP_PROC_BIND=spread OMP_PLACES=threads
hostname; lscpu | grep "Model name"; date
C=$B/iv-cases
FOMP=$B/openfast-pin-fomp/modules/aerodyn/aerodyn_driver   # the pin, Fortran OPENMP=ON
OFF=$B/fork-off/modules/aerodyn/aerodyn_driver             # fork, OPENFAST_KOKKOS=OFF (= pin, OPENMP=OFF)
KK=$B/fork-omp/modules/aerodyn/aerodyn_driver              # fork, Kokkos-OpenMP
case $what in
r6)
  ( cd $B/iv-omp && OMP_NUM_THREADS=8 ctest --timeout 900 > $OUT/l0.ctest.log 2>&1 ); grep -E "Test +#|tests passed" $OUT/l0.ctest.log
  R='aerodyn_utest|OLAF'
  for cfg in "fork-off off fast" "fork-omp omp-fast fast" "fork-omp omp-strict strict"; do
    set -- $cfg
    touch $OUT/stamp.$2
    ( cd $B/$1 && OMP_NUM_THREADS=1 OMP_PROC_BIND=false OFK_MODE=$3 ctest -R "$R" -j 8 --timeout 3000 > $OUT/$2.ctest.log 2>&1 )
    echo "=== $2"; grep -E "Test +#|tests passed" $OUT/$2.ctest.log
    mkdir -p $OUT/$2; ( cd $B/$1/reg_tests && find . -name "*.outb" -newer $OUT/stamp.$2 -print0 | xargs -0 -I{} cp --parents {} $OUT/$2/ )
  done
  echo "=== .outb channel data vs OFF"; $PY /weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos/tests/biotsavart-direct/cmp_outb.py $L $OUT/off $OUT/omp-strict $OUT/omp-fast
  ;;
r8)
  # concurrent runs must not pin to the same cores: OMP_PROC_BIND=false for the 4 parallel runs
  for c in bar-r8-tree bar-r8-direct; do
    bash $T/run_case.sh $FOMP $C/$c $OUT/$c.fomp1 OMP_NUM_THREADS=1 OMP_PROC_BIND=false &
    bash $T/run_case.sh $OFF  $C/$c $OUT/$c.off   OMP_NUM_THREADS=1 OMP_PROC_BIND=false &
    bash $T/run_case.sh $KK   $C/$c $OUT/$c.kk-strict OMP_NUM_THREADS=2 OMP_PROC_BIND=false OFK_MODE=strict &
    bash $T/run_case.sh $KK   $C/$c $OUT/$c.kk-fast   OMP_NUM_THREADS=2 OMP_PROC_BIND=false OFK_MODE=fast &
    wait
    bash $T/run_case.sh $FOMP $C/$c $OUT/$c.fomp8 OMP_NUM_THREADS=8
    bash $T/run_case.sh $KK   $C/$c $OUT/$c.kk-fast8 OMP_NUM_THREADS=8 OFK_MODE=fast
    echo "=== R8 $c (CPU): base1 = Fortran (OPENMP=ON) 1 thread, baseN = 8 threads"
    $PY $T/r8.py $L $OUT/$c.fomp1/ad_driver.out $OUT/$c.fomp8/ad_driver.out off=$OUT/$c.off/ad_driver.out \
        kk-strict=$OUT/$c.kk-strict/ad_driver.out kk-fast2=$OUT/$c.kk-fast/ad_driver.out kk-fast8=$OUT/$c.kk-fast8/ad_driver.out
  done
  ;;
r8ffma)
  R8=$2   # output dir of the r8 job (base1, runs)
  FFMA=$B/openfast-pin-ffma/modules/aerodyn/aerodyn_driver
  for c in bar-r8-tree bar-r8-direct; do
    bash $T/run_case.sh $FFMA $C/$c $OUT/$c.ffma OMP_NUM_THREADS=1 OMP_PROC_BIND=false &
  done
  wait
  for c in bar-r8-tree bar-r8-direct; do
    echo "=== R8 $c: base1 = Fortran 1 thread, envelope = Fortran built with FMA (same pin, -march=x86-64-v3)"
    $PY $T/r8.py $L $R8/$c.fomp1/ad_driver.out $OUT/$c.ffma/ad_driver.out \
        kk-strict=$R8/$c.kk-strict/ad_driver.out kk-fast8=$R8/$c.kk-fast8/ad_driver.out
  done
  ;;
perf2)
  KF=$B/fork-omp-fomp/modules/aerodyn/aerodyn_driver   # Kokkos-OpenMP + Fortran OPENMP=ON
  MT=$B/iv-omp/libofk_memtool.so
  bash $T/run_case.sh $FOMP $C/bar-r8-tree $OUT/bar-r8-tree.fomp8 OMP_NUM_THREADS=8
  bash $T/run_case.sh $KF   $C/bar-r8-tree $OUT/bar-r8-tree.kf-fast8 OMP_NUM_THREADS=8 OFK_MODE=fast
  bash $T/run_case.sh $KF   $C/bar-r8-tree $OUT/bar-r8-tree.kf-strict8 OMP_NUM_THREADS=8 OFK_MODE=strict
  bash $T/run_case.sh $KF   $C/bar-perf-direct $OUT/bar-perf-direct.kf-fast.prof OMP_NUM_THREADS=8 OFK_MODE=fast OFK_IV_PROFILE=1
  bash $T/run_case.sh $KF   $C/bar-perf-direct $OUT/bar-perf-direct.kf-fast OMP_NUM_THREADS=8 OFK_MODE=fast
  for s in packed lattice; do
    bash $T/run_case.sh $KF $C/bar-r8-direct $OUT/bar-r8-direct.kf-fast.mem.$s OMP_NUM_THREADS=8 OFK_MODE=fast OFK_IV_SOURCES=$s \
         OFK_IV_PROFILE=1 KOKKOS_TOOLS_LIBS=$MT OFK_MEMTOOL_OUT=$OUT/memtool.bar-r8-direct.$s.txt
    tail -1 $OUT/memtool.bar-r8-direct.$s.txt
  done
  ;;
perf)
  for c in bar-perf-direct bar-perf-part; do
    bash $T/run_case.sh $FOMP $C/$c $OUT/$c.fomp8 OMP_NUM_THREADS=8
    bash $T/run_case.sh $KK   $C/$c $OUT/$c.kk-fast OMP_NUM_THREADS=8 OFK_MODE=fast
    bash $T/run_case.sh $KK   $C/$c $OUT/$c.kk-fast.prof OMP_NUM_THREADS=8 OFK_MODE=fast OFK_IV_PROFILE=1
    bash $T/run_case.sh $KK   $C/$c $OUT/$c.kk-strict OMP_NUM_THREADS=8 OFK_MODE=strict
    bash $T/run_case.sh $KK   $C/$c $OUT/$c.kk-fast.lattice.prof OMP_NUM_THREADS=8 OFK_MODE=fast OFK_IV_PROFILE=1 OFK_IV_SOURCES=lattice
  done
  ;;
esac
date
