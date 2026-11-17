#!/bin/bash
# Unit-4 analysis: F10 and IEEE microbenchmarks (scripts/unit4-mb), one core, med.
#SBATCH -p med
#SBATCH -c 2
#SBATCH -t 00:20:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u4-prof/%x.%j.log
set -u
source /etc/profile.d/zz-dsai_lmod.sh; module load gcc/12.3.0 >/dev/null 2>&1
hostname; lscpu | grep "Model name"
D=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u4-prof/mb-$SLURM_JOB_ID; bash /weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos/scripts/unit4-mb/build.sh $D; cd $D
for s in 29 50 100; do for r in 0 1 2 3 4; do ./f10mb $s $r; done; done
./ieeemb; ./ieeemb
