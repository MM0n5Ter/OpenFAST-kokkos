#!/bin/bash
# R7 (criteria.md) for unit circulation-solve: every port output of the L0 run (strict, F0 and F1 variants,
# serial and parallel_for paths; OFK_L0_DUMP) is bitwise identical across processes and across
# OMP_NUM_THREADS 1, 4 and 8. Rejects empty dumps.
#   r7_threads.sh <l0_cs executable> <label> <scratch dir>
set -u
exe=$1; label=$2; dir=$3
mkdir -p "$dir"
ref=""; rc=0
for run in t1 t4 t8 t8b; do
  case $run in t1) nt=1;; t4) nt=4;; t8|t8b) nt=8;; esac
  f="$dir/l0_cs_$label.$run.bin"
  OMP_NUM_THREADS=$nt OMP_PROC_BIND=spread OMP_PLACES=threads OFK_L0_CASES=60 OFK_L0_INFO=1 OFK_L0_DUMP="$f" \
    "$exe" > "$dir/l0_cs_$label.$run.log" 2>&1 || { echo "run $run failed"; rc=1; }
  [ -s "$f" ] || { echo "run $run: empty dump $f"; rc=1; continue; }
  sum=$(md5sum < "$f" | cut -c1-32)
  echo "l0_cs OMP_NUM_THREADS=$nt ($run): $(stat -c %s "$f") bytes md5 $sum"
  if [ -z "$ref" ]; then ref=$sum; elif [ "$sum" != "$ref" ]; then echo "MISMATCH vs t1"; rc=1; fi
done
[ $rc -eq 0 ] && echo "R7 l0_cs: PASS (bitwise identical across runs and 1/4/8 threads)"
exit $rc
