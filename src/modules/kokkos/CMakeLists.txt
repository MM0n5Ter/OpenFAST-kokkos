# ofk_kernels: the C++/Kokkos side of every ported OpenFAST unit, in ONE library, so a
# process has one Kokkos runtime (OpenFAST-kokkos notes/architecture.md §1, §5).
#
# Built inside OpenFAST by the overlay (src/cmake/ofk_overlay.cmake, -DOPENFAST_KOKKOS=ON), or
# stand-alone / from a test tree with add_subdirectory():
#   cmake -S modules/kokkos -B <build> -DKokkos_ROOT=<install> [-DCMAKE_CXX_COMPILER=<install>/bin/nvcc_wrapper]
#
# Layout: src/<unit>/ per ported unit (same name as tests/<unit>/), src/common/ for what every unit
# shares. The include root is src/, so a cross-unit include names its directory
# ("biotsavart-direct/biotsavart_direct_body.hpp"); same-directory includes stay bare.
#   src/common/              ofk_runtime.{hpp,cpp}: Kokkos lifetime, OFK_MODE; ofk_vector_libm.hpp
#   src/biotsavart-direct/   FVW_BiotSavart::ui_seg, ::ui_part_nograd. The direct Biot-Savart kernels,
#                            also used by induced-velocities and treecode.
#                            biotsavart_direct.cpp (entry points), _strict.cpp / _fast.cpp (one body, two flag sets)
#   src/induced-velocities/  FVW_Subs::InducedVelocitiesAll_*, ::LiftingLineInducedVelocities: per-instance
#                            handles, lattice mirror, sources, evaluation; _pack.cpp = device packing (class B)
#   src/treecode/            FVW_VortexTools::grow_tree_part, ::ui_tree_part, through the induced-velocities
#                            handle; treecode_{strict,fast}.cpp, treecode_dbuild.hpp (device tree build)
#   src/circulation-solve/   FVW_Wings::Wings_ComputeCirculationPolarData, the Vvar block (host space):
#                            circulation.cpp (handles, C API), circulation_{strict,fast,fastfma}.cpp
# The target itself is defined in ofk_kernels.cmake, so that the overlay can include() it from a
# deferred call at the end of OpenFAST's top-level directory, where add_subdirectory() is not
# allowed ("Subdirectories may not be created during deferred execution").
cmake_minimum_required(VERSION 3.22)  # Kokkos 5 needs >= 3.22
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(ofk_kernels CXX)
endif()
include(${CMAKE_CURRENT_LIST_DIR}/ofk_kernels.cmake)

