!> A bind(C) entry point in a SHARED library, like AeroDyn_Inflow_C_Binding.f90
!! (libaerodyn_inflow_c_binding.so, loaded by pyOpenFAST through ctypes).
!! It calls the unit through its Fortran interface and does not know it was ported.
module toy_c_binding
   use iso_c_binding, only: c_int, c_double
   use toy_biotsavart, only: ui_part
   implicit none
contains
subroutine Toy_C_Induced(nCPs, CPs, nPart, Part, Alpha, RegParam, Uind) bind(C, name='Toy_C_Induced')
   integer(c_int), value         :: nCPs, nPart
   real(c_double), intent(in)    :: CPs(3,nCPs), Part(3,nPart), Alpha(3,nPart), RegParam(nPart)
   real(c_double), intent(inout) :: Uind(3,nCPs)
   call ui_part(nCPs, CPs, nPart, Part, Alpha, RegParam, Uind)
end subroutine Toy_C_Induced
end module toy_c_binding
