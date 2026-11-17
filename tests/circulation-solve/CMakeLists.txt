# Unit circulation-solve (OLAF unit 4a): oracle = the pin's FVW_BiotSavart / FVW_VortexTools / FVW_Subs /
# FVW_Wings with their modules renamed *_ref (OPENMP=OFF flags: the class-B oracle), a second copy renamed
# *_rfm built with FMA (-march=x86-64-v3: the Fortran's own FMA envelope, FMA experiment), and the port
# (the overlay's FVW_Wings + FVW_Subs + FVW_BiotSavart with OPENFAST_KOKKOS, + src/modules/kokkos).
# AirfoilInfo, FVW_Types and NWTC_Library are shared, from the pin build.
#
#   cmake -S tests/circulation-solve -B <build> -DKokkos_ROOT=<kokkos install> \
#         -DOPENFAST_PIN_BUILD=<pin build tree> -DBLA_VENDOR=OpenBLAS [-DCMAKE_CXX_COMPILER=<install>/bin/nvcc_wrapper]
#   cmake --build <build> -j 12 && ctest --test-dir <build> --output-on-failure
cmake_minimum_required(VERSION 3.22)
project(ofk_test_circulation_solve Fortran C CXX)
enable_testing()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

get_filename_component(OFK_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../.. ABSOLUTE)
set(OPENFAST_SOURCE ${OFK_ROOT}/deps/openfast CACHE PATH "The frozen baseline (pin)")
set(OPENFAST_PIN_BUILD "" CACHE PATH "Build tree of the pin")
if(NOT EXISTS ${OPENFAST_PIN_BUILD}/ftnmods/nwtc_library.mod)
  message(FATAL_ERROR "OPENFAST_PIN_BUILD=${OPENFAST_PIN_BUILD} has no ftnmods/nwtc_library.mod")
endif()
find_package(Python3 REQUIRED COMPONENTS Interpreter)

include(${OFK_ROOT}/src/cmake/ofk_sources.cmake)
ofk_add_kernels()

set(PIN_FFLAGS -ffree-line-length-none -cpp -fstack-reuse=none -fdefault-real-8 -fdefault-double-8)
set(PIN_FDEFS HAS_FORTRAN2008_FEATURES OPENFAST_DOUBLE_PRECISION)
set(BS_FFLAGS -fno-trapping-math)   # FVW_BiotSavart.f90 / FVW_VortexTools.f90 (pin aerodyn/CMakeLists.txt)

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
  ${OPENFAST_PIN_BUILD}/modules/aerodyn/libaerodynlib.a ${OPENFAST_PIN_BUILD}/modules/aerodyn/libbasicaerolib.a
  ${OPENFAST_PIN_BUILD}/modules/aerodyn/libaeroacousticslib.a ${OPENFAST_PIN_BUILD}/modules/seastate/libseastlib.a
  ${OPENFAST_PIN_BUILD}/modules/inflowwind/libifwlib.a ${OPENFAST_PIN_BUILD}/modules/nwtc-library/libnwtclibs.a
  ${OPENFAST_PIN_BUILD}/modules/version/libversioninfolib.a
  -Wl,--end-group ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

function(ofk_fortran_target tgt moddir)
  set_target_properties(${tgt} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/${moddir})
  target_include_directories(${tgt} PUBLIC ${CMAKE_BINARY_DIR}/mod/${moddir})
  target_include_directories(${tgt} PRIVATE ${PIN_INCS})
  target_compile_options(${tgt} PRIVATE ${PIN_FFLAGS})
  target_compile_definitions(${tgt} PRIVATE ${PIN_FDEFS})
endfunction()

# ---- 1. Oracles: the pin's four files, modules renamed (suffix _ref; the FMA copy _rfm)
set(AD_PIN ${OPENFAST_SOURCE}/modules/aerodyn/src)
set(FILES FVW_BiotSavart FVW_VortexTools FVW_Subs FVW_Wings)
foreach(sfx ref rfm)
  set(SRCS_${sfx})
  foreach(f IN LISTS FILES)
    set(out ${CMAKE_BINARY_DIR}/oracle_${sfx}/${f}_${sfx}.f90)
    set(cmd COMMAND Python3::Interpreter ${OFK_ROOT}/examples/fortran-kokkos/rename_modules.py ${AD_PIN}/${f}.f90 ${out} ${FILES})
    if(sfx STREQUAL rfm)   # _ref -> _rfm on the renamed module names only
      list(APPEND cmd COMMAND sed -i -E "s/\\b(FVW_BiotSavart|FVW_VortexTools|FVW_Subs|FVW_Wings)_ref\\b/\\1_rfm/Ig" ${out})
    endif()
    add_custom_command(OUTPUT ${out} ${cmd} DEPENDS ${AD_PIN}/${f}.f90 ${OFK_ROOT}/examples/fortran-kokkos/rename_modules.py VERBATIM)
    list(APPEND SRCS_${sfx} ${out})
  endforeach()
  list(GET SRCS_${sfx} 0 _bs)
  list(GET SRCS_${sfx} 1 _vt)
  set_source_files_properties(${_bs} ${_vt} PROPERTIES COMPILE_OPTIONS "${BS_FFLAGS}")
endforeach()
add_library(oracle_cs STATIC ${SRCS_ref})       # the pin's default build (OPENMP=OFF): class-B oracle
ofk_fortran_target(oracle_cs oracle)
add_library(oracle_cs_ffma STATIC ${SRCS_rfm})  # the pin built with FMA: the natural envelope
ofk_fortran_target(oracle_cs_ffma oracle_ffma)
target_compile_options(oracle_cs_ffma PRIVATE -march=x86-64-v3)

# ---- 2. Port: the overlay's files (FVW_VortexTools is the pin's), with OPENFAST_KOKKOS
ofk_port_sources(PORT_SRCS modules/aerodyn/src/FVW_BiotSavart.f90 modules/aerodyn/src/FVW_VortexTools.f90
                           modules/aerodyn/src/FVW_Subs.f90 modules/aerodyn/src/FVW_Wings.f90)
list(GET PORT_SRCS 0 _port_bs)
list(GET PORT_SRCS 1 _port_vt)
set_source_files_properties(${_port_bs} ${_port_vt} PROPERTIES COMPILE_OPTIONS "${BS_FFLAGS}")
add_library(port_cs STATIC ${PORT_SRCS})
ofk_fortran_target(port_cs port)
target_compile_definitions(port_cs PRIVATE OPENFAST_KOKKOS)
target_link_libraries(port_cs PUBLIC ofk_kernels)

# ---- 3. L0: AirfoilInfo calls are counted (ld --wrap) to read every side's iteration count
add_executable(l0_cs l0_cs.f90 afi_count.c)
ofk_fortran_target(l0_cs l0_cs)
add_dependencies(l0_cs oracle_cs oracle_cs_ffma port_cs)
target_link_libraries(l0_cs PRIVATE oracle_cs oracle_cs_ffma port_cs ofk_kernels ${PIN_LIBS}
                      -Wl,--wrap=__airfoilinfo_MOD_afi_computeairfoilcoefs)
target_compile_definitions(l0_cs PRIVATE "OFK_AFDIR=\"${OPENFAST_SOURCE}/reg_tests/r-test/modules/aerodyn/BAR_Baseline/Airfoils\"")
set_target_properties(l0_cs PROPERTIES LINKER_LANGUAGE CXX)

# ---- 4. P4: serial vs parallel_for threshold and the replay/build cost per size (perf_cs)
add_executable(perf_cs perf_cs.cpp)
target_link_libraries(perf_cs PRIVATE ofk_kernels)

set(L0_BIND "OMP_PROC_BIND=spread;OMP_PLACES=threads")
function(ofk_l0_test name env)
  add_test(NAME ${name} COMMAND l0_cs)
  set_tests_properties(${name} PROPERTIES ENVIRONMENT "${L0_BIND};${env}")
endfunction()
ofk_l0_test(l0_cs        "OFK_L0_SEED=1")
ofk_l0_test(l0_cs_seed7  "OFK_L0_SEED=7")
ofk_l0_test(l0_cs_1thread "OFK_L0_SEED=1;OMP_NUM_THREADS=1")
ofk_l0_test(l0_cs_llvel4  "OFK_L0_SEED=3;OFK_L0_LLVEL=4")   # LL induction by the Fortran segment tree (CUDA isolation)
add_test(NAME l0_cs_threads COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/r7_threads.sh $<TARGET_FILE:l0_cs> all ${CMAKE_BINARY_DIR}/r7)
