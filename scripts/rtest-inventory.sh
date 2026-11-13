#!/usr/bin/env bash
# Regression-test inventory from reg_tests/CTestList.cmake: number of cases per
# registration function (of_regression, ad_regression, hd_regression, ...),
# then the cases whose name or labels match a pattern (default: OLAF).
#
# Usage: scripts/rtest-inventory.sh [PATTERN]
#   OPENFAST_ROOT environment variable overrides deps/openfast
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${OPENFAST_ROOT:-$HERE/../deps/openfast}"
PAT="${1:-olaf}"
L="$ROOT/reg_tests/CTestList.cmake"
cases() { grep -E '^[[:space:]]*[a-z_]+_regression[a-z_]*\("' "$L"; }
echo "== cases per registration function"
cases | sed -E 's/^[[:space:]]*([a-z_]+)\(.*/\1/' | sort | uniq -c | sort -rn
echo "   total: $(cases | wc -l)"
echo "== cases matching /$PAT/ (case-insensitive)"
cases | grep -i -- "$PAT" | sed -E 's/^[[:space:]]*([a-z_]+)\("([^"]+)".*/\1 \2/' | sort
echo "   total: $(cases | grep -ic -- "$PAT" || true)"
