# For test trees (tests/<unit>/CMakeLists.txt): where the PORT's copy of an OpenFAST file lives.
# The same rule as the OpenFAST overlay (ofk_overlay.cmake): src/<rel> if the repo overlays it,
# otherwise the frozen pin's <OPENFAST_SOURCE>/<rel>.
#
#   include(<repo>/src/cmake/ofk_sources.cmake)
#   ofk_port_sources(PORT_SRCS modules/aerodyn/src/FVW_BiotSavart.f90 modules/aerodyn/src/FVW_Subs.f90)
#   ofk_add_kernels()     # add_subdirectory(<repo>/src/modules/kokkos) as ${CMAKE_BINARY_DIR}/ofk_kernels
include_guard(GLOBAL)
get_filename_component(OFK_SRC_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

function(ofk_port_sources out)
  if(NOT OPENFAST_SOURCE)
    message(FATAL_ERROR "ofk_port_sources: set OPENFAST_SOURCE (the pin) first")
  endif()
  set(_r)
  foreach(_rel IN LISTS ARGN)
    if(EXISTS "${OFK_SRC_DIR}/${_rel}")
      list(APPEND _r "${OFK_SRC_DIR}/${_rel}")
    elseif(EXISTS "${OPENFAST_SOURCE}/${_rel}")
      list(APPEND _r "${OPENFAST_SOURCE}/${_rel}")
    else()
      message(FATAL_ERROR "ofk_port_sources: ${_rel} is neither in ${OFK_SRC_DIR} nor in ${OPENFAST_SOURCE}")
    endif()
  endforeach()
  set(${out} ${_r} PARENT_SCOPE)
endfunction()

# One ofk_kernels library for every unit (one Kokkos runtime per process).
macro(ofk_add_kernels)
  add_subdirectory("${OFK_SRC_DIR}/modules/kokkos" "${CMAKE_BINARY_DIR}/ofk_kernels")
endmacro()
