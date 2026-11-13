#!/usr/bin/env bash
# Lines per OpenFAST module directory, split into hand-written Fortran,
# registry-generated Fortran (*_Types.f90), C/C++, and the number of
# registry files (generate_f90_types calls in its CMakeLists.txt). Covers modules/*, glue-codes/* and openfast_io.
#
# Usage: scripts/loc-by-module.sh [OPENFAST_ROOT]   (default: deps/openfast)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${1:-$HERE/../deps/openfast}"
cd "$ROOT"
lines() { find "$1" -type f \( "${@:2}" \) -print0 | xargs -0 cat 2>/dev/null | wc -l; }
printf '%9s %9s %8s %4s  %s\n' hand_f90 types_f90 c_cxx reg dir
for d in modules/*/ glue-codes/*/; do
  d="${d%/}"
  all=$(lines "$d" -iname '*.f90' -o -iname '*.f')
  gen=$(lines "$d" -name '*_Types.f90')
  c=$(lines "$d" -name '*.c' -o -name '*.cpp' -o -name '*.cu' -o -name '*.h' -o -name '*.hpp')
  reg=$(cat "$d"/CMakeLists.txt 2>/dev/null | grep -c 'generate_f90_types *(' || true)
  printf '%9d %9d %8d %4d  %s\n' $((all-gen)) "$gen" "$c" "$reg" "$d"
done | sort -n | awk '{print} {h+=$1; g+=$2; c+=$3; r+=$4}
  END {printf "%9d %9d %8d %4d  %s\n", h, g, c, r, "TOTAL"}'
