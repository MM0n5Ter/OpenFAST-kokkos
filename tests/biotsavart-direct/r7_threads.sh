#!/bin/bash
# R7 (criteria.md): the port's outputs are bitwise identical across runs (separate processes)
# and, for both modes, across thread counts 1, 4 and 8 on Kokkos-OpenMP (the order is fixed
# per target: strict = Fortran order, fast = 8 lanes + fixed tree, independent of threads).
#   r7_threads.sh <l0 executable> <strict|fast> <scratch dir>
# Uses OFK_L0_DUMP (l0_util.f90): every port output array of the L0 run, in order.
set -u
exe=$1; mode=$2; dir=$3
mkdir -p "$dir"
name=$(basename "$exe")_$mode
ref=""
rc=0
for run in t1a t4 t8 t8b; do
  case $run in t1a) nt=1;; t4) nt=4;; t8|t8b) nt=8;; esac
  f="$dir/$name.$run.bin"
  OMP_NUM_THREADS=$nt OMP_PROC_BIND=spread OMP_PLACES=threads OFK_MODE=$mode OFK_L0_RTOL=1 \
    OFK_L0_DUMP="$f" "$exe" > "$dir/$name.$run.log" 2>&1 || { echo "run $run failed"; tail -3 "$dir/$name.$run.log"; rc=1; }
  [ -s "$f" ] || { echo "run $run: empty dump $f (OFK_L0_DUMP not honoured)"; rc=1; }
  sum=$(md5sum < "$f" | cut -c1-32)
  echo "$name OMP_NUM_THREADS=$nt ($run): $(stat -c %s "$f") bytes md5 $sum"
  if [ -z "$ref" ]; then ref=$sum; elif [ "$sum" != "$ref" ]; then echo "MISMATCH vs t1a"; rc=1; fi
done
[ $rc -eq 0 ] && echo "R7 $name: PASS (bitwise identical across runs and 1/4/8 threads)"
exit $rc
