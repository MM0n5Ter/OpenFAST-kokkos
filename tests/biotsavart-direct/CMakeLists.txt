# Unit biotsavart-direct: oracle (pin's FVW_BiotSavart.f90 renamed *_ref), port (the overlay's
# src/modules/aerodyn/src/FVW_BiotSavart.f90 shim + src/modules/kokkos), L0 differential tests for the equivalence
# chain (architecture.md §4), API-compile check of the pin's callers, timing driver.
#
#   cmake -S tests/biotsavart-direct -B <build> -DKokkos_ROOT=<kokkos install> \
#         -DOPENFAST_PIN_BUILD=<pin build tree with nwtclibs (+ aerodynlib for the API check)> \
#         [-DCMAKE_CXX_COMPILER=<kokkos install>/bin/nvcc_wrapper -DOFK_L0_RTOL=1e-12]
#   cmake --build <build> -j && ctest --test-dir <build> --output-on-failure
# See unit.md in this directory for the exact commands used.
cmake_minimum_required(VERSION 3.22)
project(ofk_test_biotsavart_direct Fortran CXX)
enable_testing()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)   # the pin's default (-O3)
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

get_filename_component(OFK_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../.. ABSOLUTE)
set(OPENFAST_SOURCE ${OFK_ROOT}/deps/openfast CACHE PATH "The frozen baseline (pin)")
set(OPENFAST_PIN_BUILD "" CACHE PATH "Build tree of the pin (nwtclibs; aerodynlib for the API check)")
option(OFK_DOUBLE_PRECISION "ReKi = R8Ki; must match the pin build" ON)
set(OFK_L0_RTOL 0 CACHE STRING "Chain claim 1 (Fortran vs strict): 0 (class B) on host backends; claim 2 bound on CUDA")
set(OFK_L0_FAST_RTOL 1e-12 CACHE STRING "Chain claim 3 (strict vs fast): declared class R1 bound, relative to sum|contributions|")
set(OFK_L0_FMA_FACTOR 1.1 CACHE STRING "Chain claim 3 with FMA: bound = max(OFK_L0_FAST_RTOL, factor x the pin Fortran's FMA spread), per branch")
if(NOT EXISTS ${OPENFAST_PIN_BUILD}/ftnmods/nwtc_library.mod)
  message(FATAL_ERROR "OPENFAST_PIN_BUILD=${OPENFAST_PIN_BUILD} has no ftnmods/nwtc_library.mod")
endif()

find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_package(OpenMP REQUIRED COMPONENTS Fortran)

# One ofk_kernels library for every unit (one Kokkos runtime per process): src/modules/kokkos.
# The port's Fortran: src/<path> where the repo overlays the file, else the pin's (ofk_sources.cmake).
include(${OFK_ROOT}/src/cmake/ofk_sources.cmake)
ofk_add_kernels()

# ---- The pin's Fortran settings (deps/openfast/cmake/OpenfastFortranOptions.cmake, set_fast_gfortran;
# per-file -fno-trapping-math from deps/openfast/modules/aerodyn/CMakeLists.txt:92-96).
set(PIN_FFLAGS -ffree-line-length-none -cpp -fstack-reuse=none)
set(PIN_FDEFS HAS_FORTRAN2008_FEATURES)
if(OFK_DOUBLE_PRECISION)
  list(APPEND PIN_FFLAGS -fdefault-real-8 -fdefault-double-8)
  list(APPEND PIN_FDEFS OPENFAST_DOUBLE_PRECISION)
endif()
set(BS_FFLAGS -fno-trapping-math)   # FVW_BiotSavart.f90 / FVW_VortexTools.f90 only

# The pin's modules, minus the OLAF modules that the port and the API check rebuild here, so no
# stale pin fvw_biotsavart.mod can be picked up by mistake (FVW_Types stays: shared, not renamed).
set(PINMODS ${CMAKE_BINARY_DIR}/mod/pin)
file(MAKE_DIRECTORY ${PINMODS})
file(GLOB _pinmods ${OPENFAST_PIN_BUILD}/ftnmods/*.mod)
foreach(m IN LISTS _pinmods)
  get_filename_component(n ${m} NAME)
  if(NOT n MATCHES "^fvw(_biotsavart|_vortextools|_subs|_wings|_io|_tests)?\\.mod$")
    file(COPY ${m} DESTINATION ${PINMODS})
  endif()
endforeach()
set(PIN_INCS ${PINMODS} ${OPENFAST_SOURCE}/modules/nwtc-library/include)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
set(PIN_NWTC_LIBS ${OPENFAST_PIN_BUILD}/modules/nwtc-library/libnwtclibs.a
                  ${OPENFAST_PIN_BUILD}/modules/version/libversioninfolib.a ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

function(ofk_fortran_target tgt moddir)
  set_target_properties(${tgt} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/${moddir})
  target_include_directories(${tgt} PUBLIC ${CMAKE_BINARY_DIR}/mod/${moddir})
  target_include_directories(${tgt} PRIVATE ${PIN_INCS})
  target_compile_options(${tgt} PRIVATE ${PIN_FFLAGS})
  target_compile_definitions(${tgt} PRIVATE ${PIN_FDEFS})
endfunction()

# ---- 1. Oracle: the pin's file with its module renamed FVW_BiotSavart -> FVW_BiotSavart_ref.
set(REF_SRC ${OPENFAST_SOURCE}/modules/aerodyn/src/FVW_BiotSavart.f90)
set(REF_RENAMED ${CMAKE_BINARY_DIR}/oracle/FVW_BiotSavart_ref.f90)
add_custom_command(
  OUTPUT ${REF_RENAMED}
  COMMAND Python3::Interpreter ${OFK_ROOT}/examples/fortran-kokkos/rename_modules.py
          ${REF_SRC} ${REF_RENAMED} FVW_BiotSavart
  DEPENDS ${REF_SRC} ${OFK_ROOT}/examples/fortran-kokkos/rename_modules.py)
# (a) as the pin builds it by default (OPENMP=OFF): the oracle for class B
add_library(oracle_biotsavart_direct STATIC ${REF_RENAMED})
ofk_fortran_target(oracle_biotsavart_direct oracle)
target_compile_options(oracle_biotsavart_direct PRIVATE ${BS_FFLAGS})
target_link_libraries(oracle_biotsavart_direct PUBLIC ${PIN_NWTC_LIBS})
# (b) with -fopenmp (pin OPENMP=ON): the timing baseline. Its !$OMP SIMD reductions reorder sums.
add_library(oracle_biotsavart_direct_fomp STATIC ${REF_RENAMED})
ofk_fortran_target(oracle_biotsavart_direct_fomp oracle_fomp)
target_compile_options(oracle_biotsavart_direct_fomp PRIVATE ${BS_FFLAGS})
target_link_libraries(oracle_biotsavart_direct_fomp PUBLIC OpenMP::OpenMP_Fortran ${PIN_NWTC_LIBS})

# (c) with FMA (-march=x86-64-v3; gfortran contracts by default): how much the Fortran itself
#     moves when FMA is allowed. Informational envelope for chain claim 3 (fast mode has FMA).
add_library(oracle_biotsavart_direct_ffma STATIC ${REF_RENAMED})
ofk_fortran_target(oracle_biotsavart_direct_ffma oracle_ffma)
target_compile_options(oracle_biotsavart_direct_ffma PRIVATE ${BS_FFLAGS} -march=x86-64-v3)
target_link_libraries(oracle_biotsavart_direct_ffma PUBLIC ${PIN_NWTC_LIBS})

# ---- 2. Port: same module name and interface; bodies forward to ofk_kernels.
ofk_port_sources(PORT_SRC modules/aerodyn/src/FVW_BiotSavart.f90)
add_library(port_biotsavart_direct STATIC ${PORT_SRC})
ofk_fortran_target(port_biotsavart_direct port)
target_compile_options(port_biotsavart_direct PRIVATE ${BS_FFLAGS})
target_compile_definitions(port_biotsavart_direct PRIVATE OPENFAST_KOKKOS)
target_link_libraries(port_biotsavart_direct PUBLIC ofk_kernels ${PIN_NWTC_LIBS})
# A5: the same file with OPENFAST_KOKKOS undefined must be the upstream code path.
add_library(port_biotsavart_direct_off STATIC ${PORT_SRC})
ofk_fortran_target(port_biotsavart_direct_off port_off)
target_compile_options(port_biotsavart_direct_off PRIVATE ${BS_FFLAGS})
target_link_libraries(port_biotsavart_direct_off PUBLIC ${PIN_NWTC_LIBS})

# ---- 3. L0 differential tests (l0_util.f90 documents the environment)
# Every L0 executable links ofk_kernels, for the ofk_set_mode test hook.
set(L0_BIND "OMP_PROC_BIND=spread;OMP_PLACES=threads")
function(ofk_add_l0 name src oracle port)
  add_executable(${name} l0_util.f90 ${src})
  ofk_fortran_target(${name} ${name})
  target_link_libraries(${name} PRIVATE ${oracle} ${port} ofk_kernels)
  set_target_properties(${name} PROPERTIES LINKER_LANGUAGE CXX)
endfunction()
function(ofk_l0_test name exe env)
  add_test(NAME ${name} COMMAND ${exe})
  set_tests_properties(${name} PROPERTIES ENVIRONMENT "${L0_BIND};${env}")
endfunction()
foreach(r ui_seg ui_part_nograd)
  ofk_add_l0(l0_${r}        l0_${r}.f90 oracle_biotsavart_direct      port_biotsavart_direct)
  ofk_add_l0(l0_${r}_kkoff  l0_${r}.f90 oracle_biotsavart_direct      port_biotsavart_direct_off)
  ofk_add_l0(l0_${r}_fomp   l0_${r}.f90 oracle_biotsavart_direct_fomp port_biotsavart_direct)
  ofk_add_l0(l0_${r}_ffma   l0_${r}.f90 oracle_biotsavart_direct_ffma port_biotsavart_direct)
  # chain claim 1 (CPU) / claim 2 (GPU): Fortran (OPENMP=OFF) vs strict
  ofk_l0_test(l0_${r}_strict     l0_${r}       "OFK_MODE=strict;OFK_L0_RTOL=${OFK_L0_RTOL}")
  # chain claim 3: strict vs fast, same backend. Bound per branch (user decision 2026-11-14, FMA
  # kept): max(OFK_L0_FAST_RTOL, OFK_L0_FMA_FACTOR x the pin Fortran's own FMA-vs-no-FMA spread on
  # the same inputs and seed), measured first by l0_${r}_ffma (fixture fma_${r}).
  set(_fma "OFK_L0_ENVELOPE=${CMAKE_BINARY_DIR}/envelope/l0_${r}.ffma.txt;OFK_L0_ENVELOPE_FACTOR=${OFK_L0_FMA_FACTOR}")
  ofk_l0_test(l0_${r}_strictfast l0_${r}       "OFK_L0_STRICTFAST=1;OFK_L0_RTOL=${OFK_L0_FAST_RTOL};${_fma}")
  # Fortran vs fast (claim 1 + claim 3 composed), same bound
  ofk_l0_test(l0_${r}_fast       l0_${r}       "OFK_MODE=fast;OFK_L0_RTOL=${OFK_L0_FAST_RTOL};${_fma}")
  set_tests_properties(l0_${r}_strictfast l0_${r}_fast PROPERTIES FIXTURES_REQUIRED fma_${r})
  # R7 in-process: the port twice, each mode, bitwise
  ofk_l0_test(l0_${r}_portport_strict l0_${r}  "OFK_L0_PORTPORT=1;OFK_MODE=strict;OFK_L0_RTOL=0")
  ofk_l0_test(l0_${r}_portport_fast   l0_${r}  "OFK_L0_PORTPORT=1;OFK_MODE=fast;OFK_L0_RTOL=0")
  # oracle against itself: must be exactly 0
  ofk_l0_test(l0_${r}_refref     l0_${r}       "OFK_L0_REFREF=1;OFK_L0_RTOL=0")
  # A5: the shim file with OPENFAST_KOKKOS undefined is the upstream code path
  ofk_l0_test(l0_${r}_kkoff      l0_${r}_kkoff "OFK_L0_RTOL=0")
  # informational: Fortran built with -fopenmp (its !$OMP SIMD reductions reorder sums) vs strict
  ofk_l0_test(l0_${r}_fomp       l0_${r}_fomp  "OFK_MODE=strict;OFK_L0_RTOL=1e-12")
  # Fortran compiled with FMA vs strict (= Fortran without FMA, claim 1): the Fortran's own FMA
  # spread per branch, written for the claim-3 bound above (fixture setup; it passes unless it crashes)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/envelope)
  ofk_l0_test(l0_${r}_ffma       l0_${r}_ffma  "OFK_MODE=strict;OFK_L0_RTOL=1;OFK_L0_ENVELOPE_OUT=${CMAKE_BINARY_DIR}/envelope/l0_${r}.ffma.txt")
  set_tests_properties(l0_${r}_ffma PROPERTIES FIXTURES_SETUP fma_${r})
  # R7 across processes and thread counts (1, 4, 8): bitwise identical dumps of the port output
  foreach(m strict fast)
    add_test(NAME l0_${r}_threads_${m}
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/r7_threads.sh $<TARGET_FILE:l0_${r}> ${m} ${CMAKE_BINARY_DIR}/r7)
  endforeach()
endforeach()

# ---- 4. Timing driver (Fortran with -fopenmp vs the port in strict and fast mode)
add_executable(perf_direct perf_direct.f90)
ofk_fortran_target(perf_direct perf)
target_link_libraries(perf_direct PRIVATE oracle_biotsavart_direct_fomp port_biotsavart_direct ofk_kernels)
set_target_properties(perf_direct PROPERTIES LINKER_LANGUAGE CXX)

# ---- 5. API: the pin's callers, compiled UNCHANGED against the port's module, then the pin's
# own OLAF unit tests (FVW_RunTests) run on the port.
set(PIN_AD ${OPENFAST_PIN_BUILD}/modules/aerodyn/libaerodynlib.a)
if(EXISTS ${PIN_AD})
  set(AD_SRC ${OPENFAST_SOURCE}/modules/aerodyn/src)
  add_library(api_callers OBJECT
    ${AD_SRC}/FVW_VortexTools.f90 ${AD_SRC}/FVW_Subs.f90 ${AD_SRC}/FVW_Wings.f90
    ${AD_SRC}/FVW_IO.f90 ${AD_SRC}/FVW_Tests.f90 ${AD_SRC}/FVW.f90)
  ofk_fortran_target(api_callers api)
  target_include_directories(api_callers PRIVATE ${CMAKE_BINARY_DIR}/mod/port)
  set_source_files_properties(${AD_SRC}/FVW_VortexTools.f90 PROPERTIES COMPILE_OPTIONS "${BS_FFLAGS}")
  add_dependencies(api_callers port_biotsavart_direct)
  add_executable(api_fvw_runtests api_fvw_runtests.f90 $<TARGET_OBJECTS:api_callers>)
  ofk_fortran_target(api_fvw_runtests api_exe)
  target_include_directories(api_fvw_runtests PRIVATE ${CMAKE_BINARY_DIR}/mod/api)
  set(PIN_AD_LIBS ${PIN_AD} ${OPENFAST_PIN_BUILD}/modules/aerodyn/libbasicaerolib.a
      ${OPENFAST_PIN_BUILD}/modules/aerodyn/libaeroacousticslib.a
      ${OPENFAST_PIN_BUILD}/modules/seastate/libseastlib.a
      ${OPENFAST_PIN_BUILD}/modules/inflowwind/libifwlib.a)
  target_link_libraries(api_fvw_runtests PRIVATE port_biotsavart_direct
    -Wl,--start-group ${PIN_AD_LIBS} ${PIN_NWTC_LIBS} -Wl,--end-group)
  set_target_properties(api_fvw_runtests PROPERTIES LINKER_LANGUAGE CXX)
  add_test(NAME api_fvw_runtests COMMAND api_fvw_runtests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
