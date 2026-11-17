# ofk_kernels, the target (see CMakeLists.txt for the overview). Included by CMakeLists.txt (test
# trees, stand-alone) and by src/cmake/ofk_overlay.cmake (inside OpenFAST, at the end of its
# top-level directory). Every path is relative to this file, not to CMAKE_CURRENT_SOURCE_DIR.
include_guard(GLOBAL)
set(_ofk_dir ${CMAKE_CURRENT_LIST_DIR})

if(NOT TARGET Kokkos::kokkos)
  find_package(Kokkos 5.3...<6 REQUIRED)  # decision: Kokkos 5.3.0, C++20
endif()
kokkos_check(DEVICES CUDA RETURN_VALUE OFK_KOKKOS_CUDA)

# CUDA: OpenFAST has C++ of its own (e.g. modules/awae/src/vtk.cpp with fast_float.h) that nvcc
# does not compile, so CMAKE_CXX_COMPILER=nvcc_wrapper for the whole build fails. Keep g++ as the
# C++ compiler and let Kokkos's launcher swap in nvcc_wrapper only for compile and link commands
# that carry -DKOKKOS_DEPENDENCE, i.e. ofk_kernels and the links of anything that uses it.
# (Fortran commands start with gfortran and are never redirected.)
if(OFK_KOKKOS_CUDA AND NOT CMAKE_CXX_COMPILER MATCHES "nvcc_wrapper$")
  kokkos_compilation(GLOBAL)
  message(STATUS "ofk_kernels: CUDA via kokkos_compilation(GLOBAL) launcher, C++ compiler ${CMAKE_CXX_COMPILER}")
endif()

# OFK_FAST_HOST_FLAGS, OFK_FAST_FMA, OFK_FAST_VECTOR_LIBM: defined once, in src/cmake/ofk_options.cmake.
include(${_ofk_dir}/../../cmake/ofk_options.cmake)
separate_arguments(_ofk_fast_host_flags UNIX_COMMAND "${OFK_FAST_HOST_FLAGS}")
add_library(ofk_kernels STATIC
  ${_ofk_dir}/src/ofk_runtime.cpp
  ${_ofk_dir}/src/biotsavart_direct.cpp
  ${_ofk_dir}/src/biotsavart_direct_strict.cpp
  ${_ofk_dir}/src/biotsavart_direct_fast.cpp
  ${_ofk_dir}/src/induced_velocities.cpp
  ${_ofk_dir}/src/induced_velocities_pack.cpp
  ${_ofk_dir}/src/treecode_strict.cpp
  ${_ofk_dir}/src/treecode_fast.cpp
  ${_ofk_dir}/src/circulation.cpp
  ${_ofk_dir}/src/circulation_strict.cpp
  ${_ofk_dir}/src/circulation_fast.cpp
  ${_ofk_dir}/src/circulation_fastfma.cpp)
target_include_directories(ofk_kernels PUBLIC $<BUILD_INTERFACE:${_ofk_dir}/src>)
target_link_libraries(ofk_kernels PUBLIC Kokkos::kokkos)
set_target_properties(ofk_kernels PROPERTIES
  POSITION_INDEPENDENT_CODE ON   # also linked into shared libraries (C bindings)
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF)
# -fno-math-errno (gfortran's default) lets sqrt inline; -fno-trapping-math matches the pin's
# per-file flag for FVW_BiotSavart.f90 (aerodyn/CMakeLists.txt). Neither changes a value.
# Contraction is OFF by default for every TU; only *_fast.cpp turns it on (below).
target_compile_options(ofk_kernels PRIVATE
  $<$<COMPILE_LANGUAGE:CXX>:-fno-math-errno -fno-trapping-math -ffp-contract=off>)

# The FMA switch is per translation unit (architecture.md §3). The flags are appended after
# the target's, so -ffp-contract=fast overrides the default above for the fast TUs.
# Every TU with device arithmetic must be in one of these two lists: nvcc contracts by default
# (--fmad=true), so a TU in neither list would silently get FMA in device code. The pack TU of
# induced-velocities is class B in both modes, hence "strict".
set(OFK_STRICT_SOURCES ${_ofk_dir}/src/biotsavart_direct_strict.cpp ${_ofk_dir}/src/induced_velocities_pack.cpp ${_ofk_dir}/src/treecode_strict.cpp)
set(OFK_FAST_SOURCES ${_ofk_dir}/src/biotsavart_direct_fast.cpp ${_ofk_dir}/src/treecode_fast.cpp)
set(_ofk_strict_opts -ffp-contract=off)
if(OFK_FAST_FMA)
  set(_ofk_fast_opts -ffp-contract=fast ${_ofk_fast_host_flags})
else()
  set(_ofk_fast_opts -ffp-contract=off ${_ofk_fast_host_flags})
endif()
if(OFK_FAST_VECTOR_LIBM)
  list(APPEND _ofk_fast_opts -DOFK_VECTOR_LIBM=1)
endif()
if(OFK_KOKKOS_CUDA)
  list(APPEND _ofk_strict_opts --fmad=false)   # device code: no FMA contraction
  if(OFK_FAST_FMA)
    list(APPEND _ofk_fast_opts --fmad=true)    # nvcc's default, stated explicitly
  else()
    list(APPEND _ofk_fast_opts --fmad=false)
  endif()
endif()
set_source_files_properties(${OFK_STRICT_SOURCES} PROPERTIES COMPILE_OPTIONS "${_ofk_strict_opts}")
set_source_files_properties(${OFK_FAST_SOURCES} PROPERTIES COMPILE_OPTIONS "${_ofk_fast_opts}")
# Unit circulation-solve (host execution space only): one body, three TUs, for the FMA experiment
# of its fast mode (tests/circulation-solve/unit.md), independent of OFK_FAST_FMA:
#   circulation_strict.cpp   baseline ISA, no contraction             (strict)
#   circulation_fast.cpp     OFK_FAST_HOST_FLAGS, no contraction      (fast, OFK_CS_FMA=0: F0)
#   circulation_fastfma.cpp  OFK_FAST_HOST_FLAGS, contraction         (fast, OFK_CS_FMA=1: F1)
# (--fmad only matters for device code, which these TUs do not have; stated for the rule above.)
set(_ofk_cs_f0_opts -ffp-contract=off ${_ofk_fast_host_flags})
set(_ofk_cs_f1_opts -ffp-contract=fast ${_ofk_fast_host_flags})
if(OFK_KOKKOS_CUDA)
  list(APPEND _ofk_cs_f0_opts --fmad=false)
  list(APPEND _ofk_cs_f1_opts --fmad=true)
endif()
set_source_files_properties(${_ofk_dir}/src/circulation_strict.cpp PROPERTIES COMPILE_OPTIONS "${_ofk_strict_opts}")
set_source_files_properties(${_ofk_dir}/src/circulation_fast.cpp PROPERTIES COMPILE_OPTIONS "${_ofk_cs_f0_opts}")
set_source_files_properties(${_ofk_dir}/src/circulation_fastfma.cpp PROPERTIES COMPILE_OPTIONS "${_ofk_cs_f1_opts}")
message(STATUS "ofk_kernels: Kokkos ${Kokkos_VERSION} (CUDA=${OFK_KOKKOS_CUDA}); "
               "strict TUs: ${_ofk_strict_opts}; fast TUs: ${_ofk_fast_opts}")

# Inside OpenFAST (added by src/cmake/ofk_overlay.cmake): aerodynlib (and later other modules)
# link ofk_kernels, and OpenFAST exports its libraries (install(EXPORT OpenFASTLibraries)), so
# ofk_kernels must be in that export set, with a $<BUILD_INTERFACE:...> include directory.
# Consumers of the installed OpenFAST then need find_package(Kokkos) for Kokkos::kokkos.
if(CMAKE_PROJECT_NAME STREQUAL "OpenFAST")
  install(TARGETS ofk_kernels
    EXPORT "${CMAKE_PROJECT_NAME}Libraries"
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
endif()
