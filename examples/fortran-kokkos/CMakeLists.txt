# Toy for the Fortran <-> C++/Kokkos pattern used for every OpenFAST unit.
# See README.md here and notes/build-kokkos-fortran.md.
#
#   cmake -S . -B build -DKokkos_ROOT=<kokkos install> [-DCMAKE_CXX_COMPILER=<install>/bin/nvcc_wrapper]
#   cmake --build build -j && ctest --test-dir build --output-on-failure
cmake_minimum_required(VERSION 3.22)            # Kokkos 5 needs >= 3.22
project(ofk_fortran_kokkos_toy Fortran CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)         # objects also go into a shared library
find_package(Kokkos 5 REQUIRED)                 # -DKokkos_ROOT=<install prefix>
find_package(Python3 REQUIRED COMPONENTS Interpreter)  # oracle renamer, C-API test

# Shared kinds (stands in for NWTC_Library): one copy, never renamed.
add_library(toy_kinds STATIC ref/toy_kinds.f90)
set_target_properties(toy_kinds PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/kinds)
target_include_directories(toy_kinds PUBLIC ${CMAKE_BINARY_DIR}/mod/kinds)

# 1. C++/Kokkos side. In OpenFAST: ONE library for all ported units, so the
#    process has one Kokkos runtime.
add_library(ofk_kernels STATIC src/ofk_runtime.cpp src/toy_biotsavart_kokkos.cpp)
target_link_libraries(ofk_kernels PUBLIC Kokkos::kokkos)

# 2. The ported unit: same module name and interface as the reference, body is a shim.
add_library(toy_port STATIC src/toy_biotsavart.f90)
set_target_properties(toy_port PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/port)
target_include_directories(toy_port PUBLIC ${CMAKE_BINARY_DIR}/mod/port)
target_link_libraries(toy_port PUBLIC toy_kinds ofk_kernels)

# 3. The oracle: the frozen reference with its module renamed (_ref), so both can
#    be linked into one test. Only the unit's own module names are rewritten.
set(REF_RENAMED ${CMAKE_BINARY_DIR}/oracle/toy_biotsavart_ref.f90)
add_custom_command(
  OUTPUT ${REF_RENAMED}
  COMMAND Python3::Interpreter ${CMAKE_SOURCE_DIR}/rename_modules.py
          ${CMAKE_SOURCE_DIR}/ref/toy_biotsavart.f90 ${REF_RENAMED} toy_biotsavart
  DEPENDS ${CMAKE_SOURCE_DIR}/ref/toy_biotsavart.f90 ${CMAKE_SOURCE_DIR}/rename_modules.py)
add_library(toy_oracle STATIC ${REF_RENAMED})
set_target_properties(toy_oracle PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/oracle)
target_include_directories(toy_oracle PUBLIC ${CMAKE_BINARY_DIR}/mod/oracle)
target_link_libraries(toy_oracle PUBLIC toy_kinds)

# L0: differential test. CMake links it with the C++ compiler, because a static
# dependency is C++ (CXX has a higher linker preference than Fortran) and adds
# the gfortran runtime itself.
add_executable(l0_ui_part tests/l0_ui_part.f90)
target_link_libraries(l0_ui_part PRIVATE toy_oracle toy_port)
add_test(NAME l0_ui_part COMMAND l0_ui_part)
set_tests_properties(l0_ui_part PROPERTIES ENVIRONMENT "OMP_PROC_BIND=spread;OMP_PLACES=threads")

# API gate on a shared library loaded from Python (like libaerodyn_inflow_c_binding.so).
add_library(toy_c_binding SHARED src/toy_c_binding.f90)
target_link_libraries(toy_c_binding PRIVATE toy_port)
add_test(NAME capi_smoke COMMAND Python3::Interpreter ${CMAKE_SOURCE_DIR}/tests/capi_smoke.py $<TARGET_FILE:toy_c_binding>)
set_tests_properties(capi_smoke PROPERTIES ENVIRONMENT "OMP_PROC_BIND=spread;OMP_PLACES=threads")
