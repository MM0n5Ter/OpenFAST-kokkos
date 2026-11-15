# ofk_kernels: the C++/Kokkos side of every ported OpenFAST unit, in ONE library, so a
# process has one Kokkos runtime (OpenFAST-kokkos notes/architecture.md §1, §5).
#
# Built inside OpenFAST by the overlay (src/cmake/ofk_overlay.cmake, -DOPENFAST_KOKKOS=ON), or
# stand-alone / from a test tree with add_subdirectory():
#   cmake -S modules/kokkos -B <build> -DKokkos_ROOT=<install> [-DCMAKE_CXX_COMPILER=<install>/bin/nvcc_wrapper]
#
# Units in this library (each: <unit>_body.hpp = the physics, compiled twice):
#   biotsavart-direct  FVW_BiotSavart::ui_seg, ::ui_part_nograd
#     src/biotsavart_direct.cpp         entry points, views, OFK_MODE switch
#     src/biotsavart_direct_strict.cpp  strict mode: Fortran order, FMA off
#     src/biotsavart_direct_fast.cpp    fast mode: blocked order, FMA on
#   induced-velocities FVW_Subs::InducedVelocitiesAll_*, ::LiftingLineInducedVelocities (direct methods)
#     src/induced_velocities.cpp        per-instance handles, lattice mirror, sources, evaluation
#     src/induced_velocities_pack.cpp   device packing lattice -> segments -> particles (F1, class B)
#   treecode           FVW_VortexTools::grow_tree_part, ::ui_tree_part (VelocityMethod 2), through
#                      the induced-velocities handle (induced_velocities.cpp)
#     src/treecode_strict.cpp           strict mode: the Fortran tree bit for bit, FMA off
#     src/treecode_fast.cpp             fast mode: same tree, FMA on, vectorized leaves (class R2)
# The target itself is defined in ofk_kernels.cmake, so that the overlay can include() it from a
# deferred call at the end of OpenFAST's top-level directory, where add_subdirectory() is not
# allowed ("Subdirectories may not be created during deferred execution").
cmake_minimum_required(VERSION 3.22)  # Kokkos 5 needs >= 3.22
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(ofk_kernels CXX)
endif()
include(${CMAKE_CURRENT_LIST_DIR}/ofk_kernels.cmake)

