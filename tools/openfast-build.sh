#!/bin/bash
# Configure (and build) OpenFAST from the frozen pin deps/openfast with the OpenFAST-kokkos overlay
# (src/cmake/ofk_overlay.cmake; notes/overlay.md). Login node; at most -j 12.
#
#   tools/openfast-build.sh <config> <builddir> [extra cmake -D options...]
#
#   config     OPENFAST_KOKKOS  Kokkos install              Fortran OPENMP  notes
#   pin        (no overlay)     -                           OFF             the plain pin: the baseline
#   off        OFF              -                           OFF             overlay, upstream bodies (A5)
#   omp        ON               kokkos-5.3.0-omp            OFF             as the old fork-omp
#   fomp       ON               kokkos-5.3.0-omp            ON              as fork-omp-fomp (P2 config)
#   cuda-h100  ON               kokkos-5.3.0-cuda-h100      OFF             sm_90 (h100/h200 nodes)
#   cuda-a100  ON               kokkos-5.3.0-cuda-a100      OFF             sm_80 (interactive_gpu MIG)
#
# Every config: Release, DOUBLE_PRECISION=ON, BUILD_TESTING=ON, OpenBLAS, g++/gfortran 12.3,
# Python_EXECUTABLE=$OPT/venv-rtest/bin/python (numpy, bokeh, pyOpenFAST for the r-tests).
# CUDA keeps CMAKE_CXX_COMPILER=g++ (ofk_kernels calls kokkos_compilation(GLOBAL)).
# Environment: OFK_JOBS (default 12); OFK_BUILD=0 configures only; OFK_TARGETS="t1 t2" builds
# only those targets (default: all). Logs: <builddir>.cfg.log and <builddir>.build.log.
# GPU runs of a CUDA build need the portable OpenBLAS first on LD_LIBRARY_PATH (EPYC nodes):
#   export LD_LIBRARY_PATH=$OPT/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH
set -euo pipefail
if [ $# -lt 2 ]; then sed -n '2,22p' "$0"; exit 2; fi
cfg=$1; bld=$2; shift 2
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
OPT=${OFK_OPT:-/weka/scratch/jhu/zli426/zzhan381/Transformation/opt}
source /etc/profile.d/zz-dsai_lmod.sh >/dev/null 2>&1 || true
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
export NVCC_WRAPPER_DEFAULT_COMPILER=g++

# ctest's r-test scripts need numpy and bokeh (and pyOpenFAST for the *_py cases): OpenFAST records
# the Python found at configure time, so point it at the venv (OFK_PYTHON overrides).
PY=${OFK_PYTHON:-$OPT/venv-rtest/bin/python}
common=(-DCMAKE_BUILD_TYPE=Release -DDOUBLE_PRECISION=ON -DBUILD_TESTING=ON -DBLA_VENDOR=OpenBLAS
        -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_CXX_COMPILER=g++ -DPython_EXECUTABLE=$PY)
overlay=(-DCMAKE_PROJECT_INCLUDE=$REPO/src/cmake/ofk_overlay.cmake)
kk() { echo -DOPENFAST_KOKKOS=ON -DKokkos_ROOT=$OPT/kokkos-5.3.0-$1; }
case $cfg in
  pin)       args=("${common[@]}" -DOPENMP=OFF) ;;
  off)       args=("${common[@]}" "${overlay[@]}" -DOPENMP=OFF -DOPENFAST_KOKKOS=OFF) ;;
  omp)       args=("${common[@]}" "${overlay[@]}" -DOPENMP=OFF $(kk omp)) ;;
  fomp)      args=("${common[@]}" "${overlay[@]}" -DOPENMP=ON  $(kk omp)) ;;
  cuda-h100) args=("${common[@]}" "${overlay[@]}" -DOPENMP=OFF $(kk cuda-h100)) ;;
  cuda-a100) args=("${common[@]}" "${overlay[@]}" -DOPENMP=OFF $(kk cuda-a100)) ;;
  *) echo "unknown config '$cfg' (pin|off|omp|fomp|cuda-h100|cuda-a100)" >&2; exit 2 ;;
esac

echo "cmake -S $REPO/deps/openfast -B $bld ${args[*]} $*"
cmake -S "$REPO/deps/openfast" -B "$bld" "${args[@]}" "$@" > "$bld.cfg.log" 2>&1 \
  || { tail -30 "$bld.cfg.log"; echo "configure FAILED: $bld.cfg.log" >&2; exit 1; }
grep -E '^-- (ofk overlay|OPENFAST_KOKKOS|ofk_kernels)' "$bld.cfg.log" || true
[ "${OFK_BUILD:-1}" = 0 ] && exit 0
targets=()
for t in ${OFK_TARGETS:-}; do targets+=(--target "$t"); done
cmake --build "$bld" -j "${OFK_JOBS:-12}" "${targets[@]}" > "$bld.build.log" 2>&1 \
  || { tail -30 "$bld.build.log"; echo "build FAILED: $bld.build.log" >&2; exit 1; }
echo "build OK: $bld ($cfg)"
