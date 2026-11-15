# OpenFAST-kokkos overlay: builds the UNMODIFIED OpenFAST at deps/openfast with our changes.
#
#   cmake -S deps/openfast -B <build> -DCMAKE_PROJECT_INCLUDE=<repo>/src/cmake/ofk_overlay.cmake \
#         [-DOPENFAST_KOKKOS=ON -DKokkos_ROOT=<kokkos install>] <the usual OpenFAST options>
#   (or: tools/openfast-build.sh <config> <builddir>)
#
# What it does (notes/overlay.md):
#  1. Source overlay, ALWAYS (also with OPENFAST_KOKKOS=OFF): every file src/modules/<mod>/src/<f>
#     replaces deps/openfast/modules/<mod>/src/<f> in every target that compiles it. The overlay
#     files are full copies of the upstream files with the ported bodies under
#     `#ifdef OPENFAST_KOKKOS`; their `#else` branches are the upstream bodies, so an OFF build
#     is the pin's build (gate A5). Per-file source properties set by OpenFAST on the upstream file
#     (e.g. -fno-trapping-math on FVW_BiotSavart.f90) are copied to the overlay file.
#     src/modules/kokkos is not an overlay (it has no upstream counterpart): it is ofk_kernels.
#  2. OPENFAST_KOKKOS=ON: find_package(Kokkos 5.3...<6), include(src/modules/kokkos/ofk_kernels.cmake)
#     (the ofk_kernels target, in OpenFAST's top-level directory); every target that received an overlay file gets
#     OPENFAST_KOKKOS (PRIVATE) and links ofk_kernels (plain signature, as OpenFAST links them).
#     ofk_kernels installs itself into OpenFAST's export set (OpenFASTLibraries), because the
#     exported aerodynlib links it.
#
# All of 1 and 2 runs in a call deferred to the end of OpenFAST's top-level CMakeLists.txt, when
# every target exists and OpenFAST's global settings (add_definitions(-DOPENFAST_DOUBLE_PRECISION),
# OpenMP flags and link_libraries, LAPACK link options) are final, so ofk_kernels sees the same
# settings as it did when the fork added it from its top-level CMakeLists.txt.
# Needs CMake >= 3.19 (cmake_language(DEFER), set_source_files_properties(TARGET_DIRECTORY)).

# CMAKE_PROJECT_INCLUDE runs after every project() call: act once, for OpenFAST's top level only
# (not for nested projects or for the aerodyn_trapping fixture, which is a separate cmake run).
if(NOT PROJECT_NAME STREQUAL "OpenFAST" OR NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  return()
endif()
include_guard(GLOBAL)
if(CMAKE_VERSION VERSION_LESS 3.19)
  message(FATAL_ERROR "ofk_overlay.cmake needs CMake >= 3.19 (have ${CMAKE_VERSION})")
endif()

get_filename_component(OFK_REPO "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(OFK_SRC "${OFK_REPO}/src")

option(OPENFAST_KOKKOS "Run ported kernels through C++/Kokkos (OpenFAST-kokkos src/modules/kokkos, Kokkos 5.3)" OFF)
include("${OFK_SRC}/cmake/ofk_options.cmake")   # OFK_FAST_FMA, OFK_FAST_HOST_FLAGS, OFK_FAST_VECTOR_LIBM

# Source properties copied from the upstream file to its overlay copy (OpenFAST sets
# COMPILE_FLAGS on FVW_BiotSavart.f90/FVW_VortexTools.f90; the rest are for future overlays).
set(_OFK_SOURCE_PROPERTIES
  COMPILE_FLAGS COMPILE_OPTIONS COMPILE_DEFINITIONS INCLUDE_DIRECTORIES
  Fortran_FORMAT Fortran_PREPROCESS OBJECT_DEPENDS SKIP_PRECOMPILE_HEADERS)

function(_ofk_collect_targets dir out)
  get_property(_t DIRECTORY "${dir}" PROPERTY BUILDSYSTEM_TARGETS)
  get_property(_subs DIRECTORY "${dir}" PROPERTY SUBDIRECTORIES)
  foreach(_s IN LISTS _subs)
    _ofk_collect_targets("${_s}" _st)
    list(APPEND _t ${_st})
  endforeach()
  set(${out} ${_t} PARENT_SCOPE)
endfunction()

function(_ofk_overlay_apply)
  # ---- 1. the overlay files: src/modules/<mod>/src/<file>, except the ofk_kernels library
  file(GLOB _files RELATIVE "${OFK_SRC}" LIST_DIRECTORIES false "${OFK_SRC}/modules/*/src/*")
  list(FILTER _files EXCLUDE REGEX "^modules/kokkos/")
  list(SORT _files)
  set(_up_paths)
  foreach(_rel IN LISTS _files)
    if(NOT EXISTS "${CMAKE_SOURCE_DIR}/${_rel}")
      message(FATAL_ERROR "ofk overlay: src/${_rel} has no upstream file ${CMAKE_SOURCE_DIR}/${_rel}. "
                          "An overlay replaces upstream files; a new file needs its own target_sources() here.")
    endif()
    list(APPEND _up_paths "${CMAKE_SOURCE_DIR}/${_rel}")
  endforeach()

  # ---- 2. replace them in every target of the build that compiles them
  _ofk_collect_targets("${CMAKE_SOURCE_DIR}" _targets)
  set(_patched_targets)
  set(_used)
  foreach(_tgt IN LISTS _targets)
    get_target_property(_type ${_tgt} TYPE)
    if(_type STREQUAL "UTILITY" OR _type STREQUAL "GLOBAL_TARGET")
      continue()
    endif()
    get_target_property(_srcs ${_tgt} SOURCES)
    if(NOT _srcs)
      continue()
    endif()
    get_target_property(_tdir ${_tgt} SOURCE_DIR)
    set(_new)
    set(_changed FALSE)
    foreach(_s IN LISTS _srcs)
      set(_entry "${_s}")
      if(NOT _s MATCHES "^\\$<")
        get_filename_component(_abs "${_s}" ABSOLUTE BASE_DIR "${_tdir}")
        list(FIND _up_paths "${_abs}" _i)
        if(_i GREATER -1)
          list(GET _files ${_i} _rel)
          set(_entry "${OFK_SRC}/${_rel}")
          foreach(_p IN LISTS _OFK_SOURCE_PROPERTIES)
            get_source_file_property(_v "${_abs}" TARGET_DIRECTORY ${_tgt} ${_p})
            if(NOT _v STREQUAL "NOTFOUND")
              set_source_files_properties("${_entry}" TARGET_DIRECTORY ${_tgt} PROPERTIES ${_p} "${_v}")
              message(STATUS "ofk overlay:   ${_tgt}: ${_p} \"${_v}\" carried over to src/${_rel}")
            endif()
          endforeach()
          set(_changed TRUE)
          list(APPEND _used "${_rel}")
          message(STATUS "ofk overlay: ${_tgt}: ${_rel} <- src/${_rel}")
        endif()
      endif()
      list(APPEND _new "${_entry}")
    endforeach()
    if(_changed)
      set_property(TARGET ${_tgt} PROPERTY SOURCES "${_new}")
      list(APPEND _patched_targets ${_tgt})
    endif()
  endforeach()
  foreach(_rel IN LISTS _files)
    if(NOT _rel IN_LIST _used)
      message(FATAL_ERROR "ofk overlay: src/${_rel} replaces nothing: no target compiles ${CMAKE_SOURCE_DIR}/${_rel}")
    endif()
  endforeach()

  # ---- 3. OPENFAST_KOKKOS: ofk_kernels, and the macro + link for the targets that got shims
  if(OPENFAST_KOKKOS)
    find_package(Kokkos 5.3...<6 REQUIRED)
    message(STATUS "OPENFAST_KOKKOS: Kokkos ${Kokkos_VERSION} from ${Kokkos_DIR}")
    # ofk_kernels installs itself into EXPORT OpenFASTLibraries (src/modules/kokkos/ofk_kernels.cmake,
    # CMAKE_PROJECT_NAME = OpenFAST), as the exported aerodynlib now links it. A deferred call may not
    # add_subdirectory(), so the target is include()d into OpenFAST's top-level directory; it then
    # sees exactly the directory settings a subdirectory added here would inherit.
    include("${OFK_SRC}/modules/kokkos/ofk_kernels.cmake")
    foreach(_tgt IN LISTS _patched_targets)
      # ofk_kernels is a static C++ library: CMake then links every executable or shared library
      # that uses ${_tgt} with the C++ linker and adds the Fortran runtime itself.
      target_compile_definitions(${_tgt} PRIVATE OPENFAST_KOKKOS)
      target_link_libraries(${_tgt} ofk_kernels)
      message(STATUS "ofk overlay: ${_tgt}: OPENFAST_KOKKOS, links ofk_kernels")
    endforeach()
  else()
    message(STATUS "ofk overlay: OPENFAST_KOKKOS=OFF (overlay sources compile their #else = upstream bodies)")
  endif()
endfunction()

cmake_language(DEFER DIRECTORY "${CMAKE_SOURCE_DIR}" CALL _ofk_overlay_apply)
