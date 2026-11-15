#!/bin/bash
# Run one aerodyn_driver case in a private copy:  run_case.sh <driver exe> <case dir> <out dir> [VAR=value ...]
# Writes <out dir>/{ad_driver.outb, stdout.log, time.txt (wall seconds)}; the case's ../BAR_Baseline
# is linked next to the copy. OFK_IV_PROFILE / KOKKOS_TOOLS_LIBS etc. come in as VAR=value.
set -u
exe=$1; case=$2; out=$3; shift 3
rm -rf "$out"; mkdir -p "$(dirname "$out")"
cp -r "$case" "$out"
[ -e "$(dirname "$out")/BAR_Baseline" ] || ln -sfn "$(readlink -f "$(dirname "$case")/BAR_Baseline")" "$(dirname "$out")/BAR_Baseline"
cd "$out"
t0=$(date +%s.%N)
env "$@" "$exe" ad_driver.dvr > stdout.log 2>&1
rc=$?
t1=$(date +%s.%N)
awk -v a="$t0" -v b="$t1" "BEGIN{printf \"%.3f\\n\", b-a}" > time.txt
echo "$(basename "$out"): rc=$rc wall=$(cat time.txt) s  ($*)"
grep -h "^\[ofk_iv\]\|^\[ofk\]" stdout.log | head -3
exit $rc
