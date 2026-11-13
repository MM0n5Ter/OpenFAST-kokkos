!> Shared kinds. Stands in for NWTC_Library (ReKi, IntKi): it is NOT renamed
!! in the oracle build, so the reference and the port see the same kinds.
module toy_kinds
   use iso_c_binding, only: c_double, c_int
   implicit none
   integer, parameter :: ReKi  = c_double   ! OpenFAST: ReKi = R8Ki when DOUBLE_PRECISION=ON
   integer, parameter :: IntKi = c_int      ! OpenFAST: IntKi = B4Ki
end module toy_kinds
