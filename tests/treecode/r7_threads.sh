#!/bin/bash
# R7 (criteria.md) for unit treecode: l0_tc dumps the strict AND the fast outputs of every case
# (OFK_L0_DUMP); they must be bitwise identical across processes and 1/4/8 OpenMP threads (host
# backends: each target's order is fixed; the tree does not depend on threads). On a device
# backend the thread count only changes the host build, which is serial: same check.
#   r7_threads.sh <l0_tc executable> <scratch dir>
set -u
exe=$1; dir=$2
mkdir -p "$dir"
ref=""; rc=0
for run in t1a t4 t8 t8b; do
  case $run in t1a) nt=1;; t4) nt=4;; t8|t8b) nt=8;; esac
  f="$dir/l0_tc.$run.bin"
  OMP_NUM_THREADS=$nt OMP_PROC_BIND=spread OMP_PLACES=threads OFK_L0_RTOL=1 OFK_L0_INFO=1 OFK_L0_BIG=0 \
    OFK_L0_DUMP="$f" "$exe" > "$dir/l0_tc.$run.log" 2>&1 || { echo "run $run failed"; tail -3 "$dir/l0_tc.$run.log"; rc=1; }
  [ -s "$f" ] || { echo "run $run: empty dump $f"; rc=1; continue; }
  sum=$(md5sum < "$f" | cut -c1-32)
  echo "l0_tc OMP_NUM_THREADS=$nt ($run): $(stat -c %s "$f") bytes md5 $sum"
  if [ -z "$ref" ]; then ref=$sum; elif [ "$sum" != "$ref" ]; then echo "MISMATCH vs t1a"; rc=1; fi
done
[ $rc -eq 0 ] && echo "R7 l0_tc (strict and fast): PASS (bitwise identical across runs and 1/4/8 threads)"
exit $rc
