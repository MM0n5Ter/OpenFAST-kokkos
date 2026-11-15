!> API gate evidence: the pin's own OLAF unit tests (FVW_Tests.f90, compiled unchanged from
!! deps/openfast against the PORT's FVW_BiotSavart module) run on the Kokkos-backed kernels.
!! FVW_RunTests covers Test_BiotSavart_Sgmt/Part/PartTree, Test_SegmentsToPart and
!! FVW_Test_WakeInducedVelocities, among others.
program api_fvw_runtests
   use NWTC_Library
   use FVW_Tests, only: FVW_RunTests
   implicit none
   integer(IntKi)       :: ErrStat
   character(ErrMsgLen) :: ErrMsg
   call NWTC_Init()
   call FVW_RunTests(ErrStat, ErrMsg)
   write(*,'(a,i0)') 'FVW_RunTests ErrStat = ', ErrStat
   if (len_trim(ErrMsg) > 0) write(*,'(a)') trim(ErrMsg)
   if (ErrStat >= AbortErrLev) error stop 1
end program api_fvw_runtests
