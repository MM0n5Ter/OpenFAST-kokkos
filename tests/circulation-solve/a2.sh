#!/bin/bash
# A2 (criteria.md) for the overlay: every overlaid Fortran file, preprocessed WITHOUT OPENFAST_KOKKOS
# (double precision, as OpenFAST builds it), is the pin's file preprocessed the same way, line for line
# (blank lines and cpp line markers ignored: the #ifdef blocks leave blank lines behind).
#   tests/circulation-solve/a2.sh
set -u
R=$(cd "$(dirname "$0")/../.." && pwd)
module load gcc/12.3.0 >/dev/null 2>&1 || true
rc=0
for f in $(cd $R/src && find modules -name '*.f90' | sort); do
  pp() { gfortran -E -cpp -DOPENFAST_DOUBLE_PRECISION "$1" 2>/dev/null | grep -v '^#' | grep -v '^[[:space:]]*$'; }
  if cmp -s <(pp $R/src/$f) <(pp $R/deps/openfast/$f); then
    echo "A2 $f: identical to the pin without OPENFAST_KOKKOS ($(diff <(cat $R/deps/openfast/$f) <(cat $R/src/$f) | grep -c '^>') lines added under #ifdef)"
  else
    echo "A2 $f: DIFFERS from the pin without OPENFAST_KOKKOS"; rc=1
  fi
done
exit $rc
