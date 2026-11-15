!> Fortran view of the unit-treecode C API (fork/modules/kokkos/src/induced_velocities.hpp) for the
!! L0 and performance drivers: the same calls the FVW_Subs shims make (ofk_tc_grow / ofk_tc_ui),
!! with an explicit key (any object's address; the shims use m%Sgmt).
module tc_api
   use NWTC_Library, only: ReKi, IntKi
   use, intrinsic :: iso_c_binding, only: c_int, c_long_long
   implicit none
   integer(c_int), parameter :: OFK_RK = int(storage_size(1.0_ReKi)/8, c_int)
   interface
      subroutine ofk_iv_sources_packed(key, rk, method, reg, nSeg, ldPts, nPts, Pts, ldC, Connct, Gam, Eps, &
            regPart, nPart, ldP, PartP, ldA, PartAlpha, PartReg) bind(C, name='ofk_iv_sources_packed')
         import :: c_int
         type(*),        intent(in) :: key
         integer(c_int), value      :: rk, method, reg, nSeg, ldPts, nPts, ldC, regPart, nPart, ldP, ldA
         type(*),        intent(in) :: Pts(*), Gam(*), Eps(*), PartP(*), PartAlpha(*), PartReg(*)
         integer(c_int), intent(in) :: Connct(*)
      end subroutine ofk_iv_sources_packed
      subroutine ofk_iv_tree_build(key, rk, degenerate) bind(C, name='ofk_iv_tree_build')
         import :: c_int
         type(*),        intent(in) :: key
         integer(c_int), value      :: rk
         type(*),        intent(in) :: degenerate
      end subroutine ofk_iv_tree_build
      subroutine ofk_iv_eval_tree(key, rk, nCPs, capCP, ldCP, CPs, ldU, Uind, minnorm, bf, dd) bind(C, name='ofk_iv_eval_tree')
         import :: c_int
         type(*),        intent(in)    :: key
         integer(c_int), value         :: rk, nCPs, capCP, ldCP, ldU
         type(*),        intent(in)    :: CPs(*)
         type(*),        intent(inout) :: Uind(*)
         type(*),        intent(in)    :: minnorm, bf, dd
      end subroutine ofk_iv_eval_tree
      subroutine ofk_iv_reserve(key, rk, latTot, capSeg, capPart, nBlocks) bind(C, name='ofk_iv_reserve')
         import :: c_int, c_long_long
         type(*),              intent(in) :: key
         integer(c_int),       value      :: rk, nBlocks
         integer(c_long_long), value      :: latTot, capSeg, capPart
      end subroutine ofk_iv_reserve
      subroutine ofk_iv_free(key) bind(C, name='ofk_iv_free')
         type(*), intent(in) :: key
      end subroutine ofk_iv_free
      subroutine ofk_iv_counters(key, out) bind(C, name='ofk_iv_counters')
         import :: c_long_long
         type(*), intent(in) :: key
         integer(c_long_long), intent(out) :: out(9)
      end subroutine ofk_iv_counters
      subroutine ofk_iv_report() bind(C, name='ofk_iv_report')
      end subroutine ofk_iv_report
      subroutine ofk_set_mode(m) bind(C, name='ofk_set_mode')
         import :: c_int
         integer(c_int), value :: m
      end subroutine ofk_set_mode
      integer(c_int) function ofk_get_mode() bind(C, name='ofk_get_mode')
         import :: c_int
      end function ofk_get_mode
   end interface
contains
   !> The shims' ofk_tc_grow, with an explicit key: particle source set + flat tree.
   subroutine tc_grow(key, n, P, A, R, reg, capPart)
      integer,    intent(in), target :: key
      integer,    intent(in) :: n, reg, capPart
      real(ReKi), intent(in) :: P(:,:), A(:,:), R(:)
      real(ReKi) :: dR(3), degenerate
      integer(c_int) :: dI(4)
      dR = 0.0_ReKi; dI = 1
      degenerate = real(1.0e-5, ReKi)   ! as the shim (FVW_VortexTools.f90:814 literal)
      call ofk_iv_reserve(key, OFK_RK, 0_c_long_long, 0_c_long_long, int(capPart, c_long_long), 1_c_int)
      call ofk_iv_sources_packed(key, OFK_RK, 2_c_int, 0_c_int, 0_c_int, 3_c_int, 1_c_int, dR, 4_c_int, dI, dR, dR, &
         int(reg, c_int), int(n, c_int), int(size(P,1), c_int), P, int(size(A,1), c_int), A, R)
      call ofk_iv_tree_build(key, OFK_RK, degenerate)
   end subroutine tc_grow

   !> The shims' ofk_tc_ui: accumulate the tree velocity at CPs(:,1:n) into U.
   subroutine tc_ui(key, n, CPs, U, bf, dd, minnorm)
      integer,    intent(in), target :: key
      integer,    intent(in) :: n
      real(ReKi), intent(in) :: CPs(:,:), bf, dd, minnorm
      real(ReKi), intent(inout) :: U(:,:)
      if (n < 1) return
      call ofk_iv_eval_tree(key, OFK_RK, int(n, c_int), int(size(CPs,2), c_int), int(size(CPs,1), c_int), CPs, &
         int(size(U,1), c_int), U, minnorm, bf, dd)
   end subroutine tc_ui
end module tc_api
