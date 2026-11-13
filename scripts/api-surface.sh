#!/usr/bin/env bash
# External API surface in the source tree:
#  1. C-callable symbols: every bind(C, name='...') per file (C bindings, FAST_Library, ...)
#  2. Framework entry points: <Mod>_Init/_UpdateStates/_CalcOutput/_End per module
#  3. Controller DLL interface (ServoDyn BladedInterface)
#
# Usage: scripts/api-surface.sh [-v]   (-v also lists the C symbol names)
#   OPENFAST_ROOT environment variable overrides deps/openfast
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${OPENFAST_ROOT:-$HERE/../deps/openfast}"
V="${1:-}"
cd "$ROOT"
echo "== 1. bind(C, name=...) symbols per file (count  file)"
grep -rIoiE --include='*.f90' --include='*.F90' "bind *\( *c *, *name *= *['\"][A-Za-z0-9_]+['\"]" modules glue-codes \
  | grep -v '_Types\.f90:' | cut -d: -f1 | sort | uniq -c | sort -rn
if [[ "$V" == -v ]]; then
  echo "-- names"
  grep -rIoiE --include='*.f90' --include='*.F90' "bind *\( *c *, *name *= *['\"][A-Za-z0-9_]+['\"]" modules glue-codes \
    | grep -v '_Types\.f90:' | sed -E "s/.*name *= *['\"]([A-Za-z0-9_]+).*/\1/I" | sort -u | column -c 120
fi
echo "== 2. framework entry points (_Init/_UpdateStates/_CalcOutput/_End) per module"
grep -rIhoiE --include='*.f90' "^ *subroutine +[A-Za-z0-9]+_(Init|UpdateStates|CalcOutput|End) *\(" modules \
  | sed -E 's/^ *subroutine +//I; s/ *\($//' | sed -E 's/_(Init|UpdateStates|CalcOutput|End)$//I' \
  | sort | uniq -c | awk '$1>=3{printf "%s(%d) ", $2, $1} END{print ""}'
echo "== 3. controller DLL interface"
grep -nIiE "subroutine +(BladedInterface_|CallBladedDLL|DLL_controller_call)" modules/servodyn/src/BladedInterface.f90 || true
