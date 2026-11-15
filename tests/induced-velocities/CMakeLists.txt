# Unit induced-velocities: oracle (the pin's FVW_BiotSavart / FVW_VortexTools / FVW_Subs with
# their modules renamed *_ref), port (the three files from the overlay src/ where overlaid, else the
# pin: FVW_Subs with the OPENFAST_KOKKOS shims, FVW_BiotSavart with unit 1's, FVW_VortexTools
# unchanged; + src/modules/kokkos), L0 differential tests for the
# equivalence chain (architecture.md §4), A7 sequence and P6 counters, a Kokkos-Tools memory
# event logger (memtool) for P6 in the full OpenFAST build (overlay).
#
#   cmake -S tests/induced-velocities -B <build> -DKokkos_ROOT=<kokkos install> \
#         -DOPENFAST_PIN_BUILD=<pin build tree (nwtclibs, aerodynlib)> -DBLA_VENDOR=OpenBLAS \
#         [-DCMAKE_CXX_COMPILER=<kokkos install>/bin/nvcc_wrapper -DOFK_L0_RTOL=1e-12]
#   cmake --build <build> -j && ctest --test-dir <build> --output-on-failure
# FVW_Types is shared (the port changes no registry type), taken from the pin build.
cmake_minimum_required(VERSION 3.22)
project(ofk_test_induced_velocities Fortran CXX)
enable_testing()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

get_filename_component(OFK_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../.. ABSOLUTE)
set(OPENFAST_SOURCE ${OFK_ROOT}/deps/openfast CACHE PATH "The frozen baseline (pin)")
set(OPENFAST_PIN_BUILD "" CACHE PATH "Build tree of the pin")
set(OFK_L0_RTOL 0 CACHE STRING "Chain claim 1 (Fortran vs strict): 0 on host backends; claim 2 bound on CUDA")
set(OFK_L0_FAST_RTOL 1e-12 CACHE STRING "Chain claim 3 (strict vs fast), relative to sum|contributions|")
if(NOT EXISTS ${OPENFAST_PIN_BUILD}/ftnmods/nwtc_library.mod)
  message(FATAL_ERROR "OPENFAST_PIN_BUILD=${OPENFAST_PIN_BUILD} has no ftnmods/nwtc_library.mod")
endif()
find_package(Python3 REQUIRED COMPONENTS Interpreter)

include(${OFK_ROOT}/src/cmake/ofk_sources.cmake)   # ofk_port_sources(), ofk_add_kernels()
ofk_add_kernels()

# The pin's Fortran settings (as tests/biotsavart-direct)
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

# ---- 1. Oracle: the pin's three files, modules renamed X -> X_ref
set(AD_PIN ${OPENFAST_SOURCE}/modules/aerodyn/src)
set(RENAMED_MODS FVW_BiotSavart FVW_VortexTools FVW_Subs)
set(REF_SRCS)
foreach(f FVW_BiotSavart FVW_VortexTools FVW_Subs)
  set(out ${CMAKE_BINARY_DIR}/oracle/${f}_ref.f90)
  add_custom_command(OUTPUT ${out}
    COMMAND Python3::Interpreter ${OFK_ROOT}/examples/fortran-kokkos/rename_modules.py ${AD_PIN}/${f}.f90 ${out} ${RENAMED_MODS}
    DEPENDS ${AD_PIN}/${f}.f90 ${OFK_ROOT}/examples/fortran-kokkos/rename_modules.py)
  list(APPEND REF_SRCS ${out})
endforeach()
set_source_files_properties(${CMAKE_BINARY_DIR}/oracle/FVW_BiotSavart_ref.f90 ${CMAKE_BINARY_DIR}/oracle/FVW_VortexTools_ref.f90
  PROPERTIES COMPILE_OPTIONS "${BS_FFLAGS}")
# (a) the pin's default build (OPENMP=OFF): the class-B oracle
add_library(oracle_iv STATIC ${REF_SRCS})
ofk_fortran_target(oracle_iv oracle)
# (b) with FMA (-march=x86-64-v3): the Fortran's own FMA spread (informational, FMA decision)
add_library(oracle_iv_ffma STATIC ${REF_SRCS})
ofk_fortran_target(oracle_iv_ffma oracle_ffma)
target_compile_options(oracle_iv_ffma PRIVATE -march=x86-64-v3)

# ---- 2. Port: the three files with the shims (src/ overlay, else the pin), and (A5) without OPENFAST_KOKKOS
ofk_port_sources(PORT_SRCS modules/aerodyn/src/FVW_BiotSavart.f90 modules/aerodyn/src/FVW_VortexTools.f90
                           modules/aerodyn/src/FVW_Subs.f90)
list(GET PORT_SRCS 0 _port_bs)
list(GET PORT_SRCS 1 _port_vt)
set_source_files_properties(${_port_bs} ${_port_vt} PROPERTIES COMPILE_OPTIONS "${BS_FFLAGS}")
add_library(port_iv STATIC ${PORT_SRCS})
ofk_fortran_target(port_iv port)
target_compile_definitions(port_iv PRIVATE OPENFAST_KOKKOS)
target_link_libraries(port_iv PUBLIC ofk_kernels)
add_library(port_iv_off STATIC ${PORT_SRCS})
ofk_fortran_target(port_iv_off port_off)
target_link_libraries(port_iv_off PUBLIC ofk_kernels)   # ofk_set_mode / counters hooks of the test

# ---- 3. L0
set(L0_BIND "OMP_PROC_BIND=spread;OMP_PLACES=threads")
function(ofk_add_l0 name oracle port)
  add_executable(${name} l0_iv_util.f90 l0_iv.f90)
  ofk_fortran_target(${name} ${name})
  add_dependencies(${name} ${oracle} ${port})
  target_link_libraries(${name} PRIVATE ${oracle} ${port} ofk_kernels ${PIN_LIBS})
  set_target_properties(${name} PROPERTIES LINKER_LANGUAGE CXX)
endfunction()
ofk_add_l0(l0_iv        oracle_iv      port_iv)
ofk_add_l0(l0_iv_kkoff  oracle_iv      port_iv_off)
ofk_add_l0(l0_iv_ffma   oracle_iv_ffma port_iv)
function(ofk_l0_test name exe env)
  add_test(NAME ${name} COMMAND ${exe})
  set_tests_properties(${name} PROPERTIES ENVIRONMENT "${L0_BIND};${env}")
endfunction()
foreach(src packed lattice)
  # chain claim 1 (CPU) / 2 (GPU): Fortran (OPENMP=OFF) vs strict
  ofk_l0_test(l0_iv_strict_${src}      l0_iv "OFK_IV_SOURCES=${src};OFK_MODE=strict;OFK_L0_RTOL=${OFK_L0_RTOL}")
  # chain claim 3: strict vs fast, same backend
  ofk_l0_test(l0_iv_strictfast_${src}  l0_iv "OFK_IV_SOURCES=${src};OFK_L0_STRICTFAST=1;OFK_L0_RTOL=${OFK_L0_FAST_RTOL}")
  # R7 in process: the port twice (and A7: long-lived vs fresh handles), bitwise
  ofk_l0_test(l0_iv_portport_strict_${src} l0_iv "OFK_IV_SOURCES=${src};OFK_L0_PORTPORT=1;OFK_MODE=strict;OFK_L0_RTOL=0")
  ofk_l0_test(l0_iv_portport_fast_${src}   l0_iv "OFK_IV_SOURCES=${src};OFK_L0_PORTPORT=1;OFK_MODE=fast;OFK_L0_RTOL=0")
  # informational: Fortran vs fast (claims 1+3)
  ofk_l0_test(l0_iv_fast_${src}        l0_iv "OFK_IV_SOURCES=${src};OFK_MODE=fast;OFK_L0_RTOL=1e-10")
endforeach()
ofk_l0_test(l0_iv_refref  l0_iv       "OFK_L0_REFREF=1;OFK_L0_RTOL=0")
ofk_l0_test(l0_iv_kkoff   l0_iv_kkoff "OFK_L0_RTOL=0")
ofk_l0_test(l0_iv_ffma    l0_iv_ffma  "OFK_MODE=strict;OFK_L0_RTOL=1;OFK_L0_INFO=1")
foreach(m strict fast)
  add_test(NAME l0_iv_threads_${m}
           COMMAND bash ${OFK_ROOT}/tests/biotsavart-direct/r7_threads.sh $<TARGET_FILE:l0_iv> ${m} ${CMAKE_BINARY_DIR}/r7)
endforeach()

# ---- 4. Kokkos-Tools memory-event logger (P6 in the full OpenFAST build): KOKKOS_TOOLS_LIBS=<this .so>
add_library(ofk_memtool SHARED memtool.cpp)
set_target_properties(ofk_memtool PROPERTIES CXX_STANDARD 17)
