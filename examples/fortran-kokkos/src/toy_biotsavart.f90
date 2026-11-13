!> PORTED UNIT. Same module name and same public interface as
!! ref/toy_biotsavart.f90, so every caller compiles unchanged (API gate).
!! The body is now a thin shim that forwards to C++/Kokkos through bind(C).
module toy_biotsavart
   use iso_c_binding, only: c_int, c_double
   use toy_kinds,     only: ReKi, IntKi
   implicit none
   private
   public :: ui_part

   interface
      !> Implemented in toy_biotsavart_kokkos.cpp. Explicit-shape (3,*) dummies:
      !! a contiguous actual is passed by address; a non-contiguous one is copied in/out
      !! by the compiler. real(ReKi) must be real(c_double): if ReKi changes (single
      !! precision build), this interface stops compiling, which is what we want.
      subroutine ofk_toy_ui_part(nCPs, CPs, nPart, Part, Alpha, RegParam, Uind) &
            bind(C, name='ofk_toy_ui_part')
         import :: c_int, c_double
         integer(c_int), value         :: nCPs, nPart
         real(c_double), intent(in)    :: CPs(3,*), Part(3,*), Alpha(3,*), RegParam(*)
         real(c_double), intent(inout) :: Uind(3,*)
      end subroutine ofk_toy_ui_part
   end interface

contains

subroutine ui_part(nCPs, CPs, nPart, Part, Alpha, RegParam, Uind)
   integer(IntKi),             intent(in)    :: nCPs, nPart
   real(ReKi), dimension(:,:), intent(in)    :: CPs
   real(ReKi), dimension(:,:), intent(in)    :: Part
   real(ReKi), dimension(:,:), intent(in)    :: Alpha
   real(ReKi), dimension(:),   intent(in)    :: RegParam
   real(ReKi), dimension(:,:), intent(inout) :: Uind
   if (nCPs <= 0 .or. nPart <= 0) return
   call ofk_toy_ui_part(int(nCPs, c_int), CPs, int(nPart, c_int), Part, Alpha, RegParam, Uind)
end subroutine ui_part

end module toy_biotsavart
