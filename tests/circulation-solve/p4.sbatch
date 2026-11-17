#!/bin/bash
# P4 of unit circulation-solve: perf_cs (build and replay time per call, serial loop vs Kokkos host
# parallel_for, strict / F0 / F1) at 1, 2, 4, 8 threads, nCP from 3 to 300. med partition, 8 cores.
#   sbatch -J ofk-circulation-solve-p4 p4.sbatch <perf_cs binary>
#SBATCH -p med
#SBATCH -c 8
#SBATCH -t 00:30:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u4a-jobs/%x.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 openblas/0.3.28 >/dev/null 2>&1
hostname; lscpu | grep "Model name"; date
export OMP_PROC_BIND=spread OMP_PLACES=threads
for t in 1 2 4 8; do
  OMP_NUM_THREADS=$t $1 5 3 3 6 9 12 15 18 24 30 36 48 60 87 120 150 300 | grep csv
done
date
