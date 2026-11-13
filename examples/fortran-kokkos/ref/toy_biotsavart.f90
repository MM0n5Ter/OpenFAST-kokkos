!> FROZEN REFERENCE. Stands in for an unmodified file in deps/openfast
!! (a cut-down OLAF ui_part_nograd with the exponential core, FVW_BiotSavart.f90).
!! Never edit this file; the L0 test compiles it under a renamed module as the oracle.
module toy_biotsavart
   use toy_kinds, only: ReKi, IntKi
   implicit none
   private
   public :: ui_part
   real(ReKi), parameter :: fourpi_inv = 0.25_ReKi / acos(-1.0_ReKi)
   real(ReKi), parameter :: MINNORM2   = 1.0e-8_ReKi
contains

!> Induced velocity of nPart vortex particles at nCPs control points.
!! Uind is ACCUMULATED into (intent inout, "side effects" as in OLAF), and only
!! columns 1:nCPs are touched although the arrays may be larger ("3 x nCPs++").
subroutine ui_part(nCPs, CPs, nPart, Part, Alpha, RegParam, Uind)
   integer(IntKi),             intent(in)    :: nCPs, nPart
   real(ReKi), dimension(:,:), intent(in)    :: CPs       !< 3 x nCPs++
   real(ReKi), dimension(:,:), intent(in)    :: Part      !< 3 x nPart++
   real(ReKi), dimension(:,:), intent(in)    :: Alpha     !< 3 x nPart++, particle intensity
   real(ReKi), dimension(:),   intent(in)    :: RegParam  !< nPart++, core radius
   real(ReKi), dimension(:,:), intent(inout) :: Uind      !< 3 x nCPs++
   integer    :: icp, ip
   real(ReKi) :: dx, dy, dz, r2, r, q, sp, ux, uy, uz
   do icp = 1, nCPs
      ux = 0.0_ReKi; uy = 0.0_ReKi; uz = 0.0_ReKi
      do ip = 1, nPart
         dx = CPs(1,icp) - Part(1,ip)
         dy = CPs(2,icp) - Part(2,ip)
         dz = CPs(3,icp) - Part(3,ip)
         r2 = dx*dx + dy*dy + dz*dz
         if (r2 < MINNORM2) cycle                       ! singularity guard
         r  = sqrt(r2)
         q  = r / RegParam(ip)
         sp = fourpi_inv / (r2*r) * (1.0_ReKi - exp(-q*q*q))
         ux = ux + sp * (Alpha(2,ip)*dz - Alpha(3,ip)*dy)
         uy = uy + sp * (Alpha(3,ip)*dx - Alpha(1,ip)*dz)
         uz = uz + sp * (Alpha(1,ip)*dy - Alpha(2,ip)*dx)
      end do
      Uind(1,icp) = Uind(1,icp) + ux
      Uind(2,icp) = Uind(2,icp) + uy
      Uind(3,icp) = Uind(3,icp) + uz
   end do
end subroutine ui_part

end module toy_biotsavart
