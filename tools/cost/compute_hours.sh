#!/bin/bash
# CPU core-hours and GPU-hours used by this project's Slurm jobs, for the cost ledger.
#
# Convention: every job for this project is named "ofk-<unit>-<what>" (sbatch/srun -J ofk-...).
# Usage: tools/cost/compute_hours.sh [SINCE=2026-11-01] [NAME_PREFIX=ofk-]
# Prints one line per job and a total. Pending/cancelled-before-start jobs count 0.
set -euo pipefail
SINCE=${1:-2026-11-01}
PREFIX=${2:-ofk-}
sacct -u "$USER" -S "$SINCE" -X -n -P --format=JobID,JobName,Partition,State,ElapsedRaw,AllocCPUS,AllocTRES \
| awk -F'|' -v p="$PREFIX" '
  index($2, p) == 1 {
    gpus = 0
    if (match($7, /gres\/gpu[^=]*=[0-9]+/)) { s = substr($7, RSTART, RLENGTH); sub(/.*=/, "", s); gpus = s + 0 }
    cpu_h = $5 * $6 / 3600; gpu_h = $5 * gpus / 3600
    printf "%-12s %-24s %-10s %-10s cpu_core_h=%.3f gpu_h=%.3f\n", $1, $2, $3, $4, cpu_h, gpu_h
    tc += cpu_h; tg += gpu_h
  }
  END { printf "TOTAL since '"$SINCE"': cpu_core_h=%.3f gpu_h=%.3f\n", tc, tg }'
