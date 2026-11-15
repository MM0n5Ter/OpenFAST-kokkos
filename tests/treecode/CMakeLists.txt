# Unit treecode: oracle (the pin's FVW_BiotSavart / FVW_VortexTools with their modules renamed
# *_ref), port (src/modules/kokkos through the induced-velocities handle's C API, the calls the
# FVW_Subs shims make), L0 (claims 1-3, R2, R7, P6) and the P1/P5 scaling driver.
#
#   cmake -S tests/treecode -B <build> -DKokkos_ROOT=<kokkos install> \
#         -DOPENFAST_PIN_BUILD=<pin build tree (nwtclibs)> -DBLA_VENDOR=OpenBLAS \
#         [-DCMAKE_CXX_COMPILER=<kokkos install>/bin/nvcc_wrapper -DOFK_L0_RTOL=1e-13]
#   cmake --build <build> -j && ctest --test-dir <build> --output-on-failure
# The system paths (InducedVelocitiesAll_Init/Calc, LiftingLineInducedVelocities with
# VelocityMethod 2) are covered by tests/induced-velocities (l0_iv), built from the same src/ overlay.
cmake_minimum_required(VERSION 3.22)
project(ofk_test_treecode Fortran CXX)
enable_testing()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

get_filename_component(OFK_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../.. ABSOLUTE)
set(OPENFAST_SOURCE ${OFK_ROOT}/deps/openfast CACHE PATH "The frozen baseline (pin)")
set(OPENFAST_PIN_BUILD "" CACHE PATH "Build tree of the pin")
set(OFK_L0_RTOL 0 CACHE STRING "Chain claim 1 (Fortran vs strict): 0 on host backends; claim 2 bound on CUDA")
if(NOT EXISTS ${OPENFAST_PIN_BUILD}/ftnmods/nwtc_library.mod)
  message(FATAL_ERROR "OPENFAST_PIN_BUILD=${OPENFAST_PIN_BUILD} has no ftnmods/nwtc_library.mod")
endif()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_package(OpenMP REQUIRED COMPONENTS Fortran)

include(${OFK_ROOT}/src/cmake/ofk_sources.cmake)   # ofk_add_kernels(): src/modules/kokkos (notes/overlay.md)
ofk_add_kernels()

# The pin's Fortran settings (as tests/biotsavart-direct, tests/induced-velocities)
set(PIN_FFLAGS -ffree-line-length-none -cpp -fstack-reuse=none -fdefault-real-8 -fdefault-double-8)
set(PIN_FDEFS HAS_FORTRAN2008_FEATURES OPENFAST_DOUBLE_PRECISION)
set(BS_FFLAGS -fno-trapping-math)   # FVW_BiotSavart.f90 / FVW_VortexTools.f90 (pin aerodyn/CMakeLists.txt:92-96)

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
set(PIN_LIBS
  -Wl,--start-group
  ${OPENFAST_PIN_BUILD}/modules/nwtc-library/libnwtclibs.a ${OPENFAST_PIN_BUILD}/modules/version/libversioninfolib.a
  -Wl,--end-group ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

function(ofk_fortran_target tgt moddir)
  set_target_properties(${tgt} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/${moddir})
  target_include_directories(${tgt} PUBLIC ${CMAKE_BINARY_DIR}/mod/${moddir})
  target_include_directories(${tgt} PRIVATE ${PIN_INCS})
  target_compile_options(${tgt} PRIVATE ${PIN_FFLAGS})
  target_compile_definitions(${tgt} PRIVATE ${PIN_FDEFS})
endfunction()

# ---- 1. Oracle: the pin's two files, modules renamed X -> X_ref
set(AD_PIN ${OPENFAST_SOURCE}/modules/aerodyn/src)
set(RENAMED_MODS FVW_BiotSavart FVW_VortexTools)
set(REF_SRCS)
foreach(f FVW_BiotSavart FVW_VortexTools)
  set(out ${CMAKE_BINARY_DIR}/oracle/${f}_ref.f90)
  add_custom_command(OUTPUT ${out}
    COMMAND Python3::Interpreter ${OFK_ROOT}/examples/fortran-kokkos/rename_modules.py ${AD_PIN}/${f}.f90 ${out} ${RENAMED_MODS}
    DEPENDS ${AD_PIN}/${f}.f90 ${OFK_ROOT}/examples/fortran-kokkos/rename_modules.py)
  list(APPEND REF_SRCS ${out})
endforeach()
set_source_files_properties(${REF_SRCS} PROPERTIES COMPILE_OPTIONS "${BS_FFLAGS}")
# (a) the pin's default build (OPENMP=OFF): the class-B oracle
add_library(oracle_tc STATIC ${REF_SRCS})
ofk_fortran_target(oracle_tc oracle)
# (b) the pin with OpenMP (OPENMP=ON adds -fopenmp): the P1/P5 baseline "Fortran+OpenMP"
add_library(oracle_tc_omp STATIC ${REF_SRCS})
ofk_fortran_target(oracle_tc_omp oracle_omp)
target_compile_options(oracle_tc_omp PRIVATE ${OpenMP_Fortran_FLAGS})
# (c) with FMA (-march=x86-64-v3): the Fortran tree's own FMA spread (informational)
add_library(oracle_tc_ffma STATIC ${REF_SRCS})
ofk_fortran_target(oracle_tc_ffma oracle_ffma)
target_compile_options(oracle_tc_ffma PRIVATE -march=x86-64-v3)

# ---- 2. L0 and the scaling driver
function(ofk_tc_exe name src oracle)
  add_executable(${name} tc_api.f90 ${src})
  ofk_fortran_target(${name} ${name})
  add_dependencies(${name} ${oracle})
  target_link_libraries(${name} PRIVATE ${oracle} ofk_kernels ${PIN_LIBS})
  set_target_properties(${name} PROPERTIES LINKER_LANGUAGE CXX)
endfunction()
ofk_tc_exe(l0_tc l0_tc.f90 oracle_tc)
ofk_tc_exe(l0_tc_ffma l0_tc.f90 oracle_tc_ffma)
ofk_tc_exe(perf_tc perf_tc.f90 oracle_tc_omp)
target_link_libraries(perf_tc PRIVATE OpenMP::OpenMP_Fortran)

set(L0_BIND "OMP_PROC_BIND=spread;OMP_PLACES=threads")
add_test(NAME l0_tc COMMAND l0_tc)                       # claims 1 (CPU) / 2 (GPU), R2, claim 3 info, P6
set_tests_properties(l0_tc PROPERTIES ENVIRONMENT "${L0_BIND};OFK_L0_RTOL=${OFK_L0_RTOL}")
add_test(NAME l0_tc_seed7 COMMAND l0_tc)
set_tests_properties(l0_tc_seed7 PROPERTIES ENVIRONMENT "${L0_BIND};OFK_L0_RTOL=${OFK_L0_RTOL};OFK_L0_SEED=7;OFK_L0_CASES=80;OFK_L0_BIG=0")
add_test(NAME l0_tc_ffma COMMAND l0_tc_ffma)             # informational: Fortran-with-FMA vs strict
set_tests_properties(l0_tc_ffma PROPERTIES ENVIRONMENT "${L0_BIND};OFK_L0_RTOL=1;OFK_L0_INFO=1;OFK_L0_BIG=0")
add_test(NAME l0_tc_r7 COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/r7_threads.sh $<TARGET_FILE:l0_tc> ${CMAKE_BINARY_DIR}/r7)
