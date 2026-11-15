#!/bin/bash
#SBATCH -J ofk-biotsavart-direct-r6
#SBATCH -p med
#SBATCH -c 8
#SBATCH -t 01:30:00
#SBATCH -o r6.%j.log
# R6 (criteria.md): aerodyn_utest and the OLAF regression tests on the fork built with
# OPENFAST_KOKKOS=OFF (the baseline: identical to the pin) and ON with Kokkos-OpenMP, in fast
# (default) and strict mode, OMP_NUM_THREADS=1 (testing.md). Also: the strict run's .outb files
# must be bitwise identical to the OFF build's (chain claim 1 at simulation level).
#   sbatch --chdir=$OPT/build $REPO/tests/biotsavart-direct/r6.sbatch
set -u
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
OUT=$B/r6-$SLURM_JOB_ID
mkdir -p $OUT
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false
hostname; date
R='aerodyn_utest|OLAF'
run() {  # <build> <label> [env...]
  local tree=$1 label=$2; shift 2
  echo "=== $label: ctest -R '$R' in $tree ($*)"
  ( cd $B/$tree && env "$@" ctest -R "$R" -j 8 --timeout 3000 > $OUT/$label.ctest.log 2>&1 )
  grep -E "tests passed|Passed|Failed|\*\*\*|Not Run|Timeout" $OUT/$label.ctest.log | sed 's/^ *//' | grep -v "^Start" | tail -25
  mkdir -p $OUT/$label
  ( cd $B/$tree/reg_tests && find . -name "*.outb" -newer $OUT/stamp.$label -print0 | xargs -0 -I{} cp --parents {} $OUT/$label/ )
}
for l in off omp-fast omp-strict; do touch $OUT/stamp.$l; done
touch $OUT/stamp.off;        run fork-off off        OFK_MODE=fast
touch $OUT/stamp.omp-fast;   run fork-omp omp-fast   OFK_MODE=fast
touch $OUT/stamp.omp-strict; run fork-omp omp-strict OFK_MODE=strict
echo "=== .outb bitwise comparison against the OFF build"
for l in omp-strict omp-fast; do
  same=0; diff=0
  while IFS= read -r f; do
    if cmp -s "$OUT/off/$f" "$OUT/$l/$f"; then same=$((same+1)); else diff=$((diff+1)); echo "  $l differs: $f"; fi
  done < <(cd $OUT/off && find . -name "*.outb")
  echo "$l: $same .outb files bitwise identical to OFF, $diff different"
done
date
