#!/bin/bash
# Unit induced-velocities: GPU job. Submitted early (the h100 queue waits hours); it runs the plan
# in gpu_steps.sh, which is filled in while the job is pending. If the plan is not marked ready
# when the job starts, it polls for the flag (at most OFK_WAIT_MIN minutes) before running.
#SBATCH -p h100,h200
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH -t 02:00:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/iv-gpu.%j.log
set -u
T=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos/tests/induced-velocities
W=${OFK_WAIT_MIN:-60}
for i in $(seq 1 $W); do [ -f $T/.gpu_ready ] && break; sleep 60; done
echo "== start $(date) on $(hostname) waited=${i}min"; nvidia-smi -L
bash $T/${OFK_GPU_PLAN:-gpu_steps.sh}
echo "== end $(date)"
