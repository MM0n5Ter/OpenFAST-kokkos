#!/usr/bin/env bash
# Where is a routine called from? Counts the Fortran files outside the owning
# module directory that call each routine, and the module directories they are in.
#
# Usage: scripts/callers.sh [-o OWNER] ROUTINE...
#   -o OWNER   module directory under modules/ whose own files are excluded (default aerodyn)
#   OPENFAST_ROOT environment variable overrides deps/openfast
# Example: scripts/callers.sh -o inflowwind IfW_FlowField_GetVelAcc
#          scripts/callers.sh -o nwtc-library Transfer_Line2_to_Line2
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${OPENFAST_ROOT:-$HERE/../deps/openfast}"
OWNER=aerodyn
if [[ "${1:-}" == -o ]]; then OWNER="$2"; shift 2; fi
[[ $# -gt 0 ]] || { sed -n '2,10p' "$0"; exit 1; }
cd "$ROOT"
for r in "$@"; do
  hits=$(grep -rliE --include='*.f90' --include='*.F90' --include='*.f' \
           "(call +$r\\b|\\b$r *\\()" -- modules glue-codes unit_tests \
         | grep -v "^modules/$OWNER/" | grep -v '_Types\.f90$' || true)
  nf=$(printf '%s' "$hits" | grep -c . || true)
  dirs=$(printf '%s\n' "$hits" | cut -d/ -f1-2 | { grep . || true; } | sort -u | tr '\n' ' ')
  nd=$(printf '%s' "$dirs" | wc -w)
  printf '%-32s files=%-4d dirs=%-3d %s\n' "$r" "$nf" "$nd" "$dirs"
done
