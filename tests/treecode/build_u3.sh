#!/bin/bash
# Build trees of unit treecode (outside the repo, opt/build/u3-* and opt/build/ov3-*). Login node.
# Since 2026-11-14 the port comes from the overlay (src/, notes/overlay.md); the full OpenFAST
# builds are tools/openfast-build.sh builds of deps/openfast, named ov3-* (the old fork3-* trees
# were built from fork/ and are left as they are).
#   [U3=u3d] build_u3.sh <what>...   what: omp | bsd-omp | cuda-a100 | cuda-h100 | iv-omp | iv-cuda-a100 |
#                                 iv-cuda-h100 | ov3-off | ov3-fomp | ov3-cuda-a100 | ov3-cuda-h100
# U3 (default u3) prefixes the test trees: never rebuild a tree a queued job uses (unit 3 used u3d-*
# for the device build, because u3-cuda-* were held by queued jobs).
set -u
source /etc/profile.d/zz-dsai_lmod.sh >/dev/null 2>&1; module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
O=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt; R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
export NVCC_WRAPPER_DEFAULT_COMPILER=g++
U3=${U3:-u3}
common=(-DCMAKE_BUILD_TYPE=Release -DBLA_VENDOR=OpenBLAS -DCMAKE_Fortran_COMPILER=gfortran
        -DOFK_FAST_FMA=ON "-DOFK_FAST_HOST_FLAGS=-march=x86-64-v3" -DOFK_FAST_VECTOR_LIBM=OFF)
pin=(-DOPENFAST_PIN_BUILD=$O/build/openfast-pin -DOPENFAST_SOURCE=$R/deps/openfast)
run() {  # name src extra-args...
  local name=$1 src=$2; shift 2
  cmake -S $src -B $O/build/$name "${common[@]}" "$@" > $O/build/$name.cfg.log 2>&1 &&
    cmake --build $O/build/$name -j 12 > $O/build/$name.build.log 2>&1
  echo "$name rc=$?"
}
for w in "$@"; do case $w in
  omp)          run $U3-omp $R/tests/treecode -DKokkos_ROOT=$O/kokkos-5.3.0-omp "${pin[@]}" ;;
  cuda-a100)    run $U3-cuda-a100 $R/tests/treecode -DKokkos_ROOT=$O/kokkos-5.3.0-cuda-a100 -DCMAKE_CXX_COMPILER=$O/kokkos-5.3.0-cuda-a100/bin/nvcc_wrapper "${pin[@]}" -DOFK_L0_RTOL=1e-13 ;;
  cuda-h100)    run $U3-cuda-h100 $R/tests/treecode -DKokkos_ROOT=$O/kokkos-5.3.0-cuda-h100 -DCMAKE_CXX_COMPILER=$O/kokkos-5.3.0-cuda-h100/bin/nvcc_wrapper "${pin[@]}" -DOFK_L0_RTOL=1e-13 ;;
  bsd-omp)      run $U3-bsd-omp $R/tests/biotsavart-direct -DKokkos_ROOT=$O/kokkos-5.3.0-omp "${pin[@]}" -DOFK_DOUBLE_PRECISION=ON ;;
  iv-omp)       run $U3-iv-omp $R/tests/induced-velocities -DKokkos_ROOT=$O/kokkos-5.3.0-omp "${pin[@]}" -DOFK_L0_FAST_RTOL=1e-12 ;;
  iv-cuda-a100) run $U3-iv-cuda-a100 $R/tests/induced-velocities -DKokkos_ROOT=$O/kokkos-5.3.0-cuda-a100 -DCMAKE_CXX_COMPILER=$O/kokkos-5.3.0-cuda-a100/bin/nvcc_wrapper "${pin[@]}" -DOFK_L0_FAST_RTOL=1e-12 -DOFK_L0_RTOL=1e-14 ;;
  iv-cuda-h100) run $U3-iv-cuda-h100 $R/tests/induced-velocities -DKokkos_ROOT=$O/kokkos-5.3.0-cuda-h100 -DCMAKE_CXX_COMPILER=$O/kokkos-5.3.0-cuda-h100/bin/nvcc_wrapper "${pin[@]}" -DOFK_L0_FAST_RTOL=1e-12 -DOFK_L0_RTOL=1e-14 ;;
  ov3-off|ov3-fomp|ov3-cuda-a100|ov3-cuda-h100)   # full OpenFAST: the pin + the overlay
    OFK_JOBS=12 $R/tools/openfast-build.sh ${w#ov3-} $O/build/$w; echo "$w rc=$?" ;;
esac; done
