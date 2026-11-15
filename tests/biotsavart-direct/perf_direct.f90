!> Timing driver for P2 (not the full P gate): the pin's Fortran (FVW_BiotSavart_ref, built with
!! -fopenmp in the perf target, i.e. as OPENMP=ON builds it) vs the port in strict and in fast
!! mode, same inputs, same process.
!!   perf_direct <seg|part> <nTargets> <nSources> <RegFunction> <reps>
!! Prints one CSV row:
!!   routine,nCP,nSrc,reg,ref_s,strict_s,fast_s,strict_speedup,fast_speedup,relmax_strict,relmax_fast
!! (min over reps; times are wall clock per call; speedup = ref_s/port_s; relmax_* =
!! max|port - ref| / max|ref|, informational: the -fopenmp Fortran is itself reordered).
!! The first call of each mode is a warm-up (Kokkos initialization, first allocation).
!! OFK_BSD_PROFILE=1 prints the port's alloc+H2D / kernel / D2H split on stderr.
program perf_direct
   use NWTC_Library,       only: ReKi, IntKi
   use FVW_BiotSavart_ref, only: seg_ref => ui_seg, part_ref => ui_part_nograd
   use FVW_BiotSavart,     only: seg_port => ui_seg, part_port => ui_part_nograd
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none
   interface
      subroutine ofk_set_mode(m) bind(C, name='ofk_set_mode')
         import :: c_int
         integer(c_int), value :: m
      end subroutine ofk_set_mode
   end interface
   character(len=16) :: arg, routine
   integer :: nCP, nSrc, reg, reps, i, it
   real(ReKi), allocatable :: CPs(:,:), Pts(:,:), Gam(:), Eps(:), U1(:,:), U2(:,:), U3(:,:), Alpha(:,:)
   integer(IntKi), allocatable :: Connct(:,:)
   real(ReKi) :: rel_s, rel_f, sc
   real(8) :: tref, ts, tf
   integer(8) :: rate

   call get_command_argument(1, routine)
   call get_command_argument(2, arg); read(arg,*) nCP
   call get_command_argument(3, arg); read(arg,*) nSrc
   call get_command_argument(4, arg); read(arg,*) reg
   call get_command_argument(5, arg); read(arg,*) reps
   call random_seed()
   allocate(CPs(3,nCP), U1(3,nCP), U2(3,nCP), U3(3,nCP), Eps(nSrc), Gam(nSrc))
   call random_number(CPs); CPs = 10*CPs - 5
   call random_number(Eps); Eps = 0.05_ReKi + 0.5_ReKi*Eps
   call random_number(Gam); Gam = 2*Gam - 1
   if (routine == 'seg') then
      allocate(Pts(3, nSrc+1), Connct(4, nSrc))
      call random_number(Pts); Pts = 10*Pts - 5
      do i = 1, nSrc            ! a random polyline: neighbour points, like a wake filament
         Connct(:, i) = [i, i+1, 1, 1]
      end do
   else
      allocate(Pts(3, nSrc), Alpha(3, nSrc))
      call random_number(Pts); Pts = 10*Pts - 5
      call random_number(Alpha); Alpha = 2*Alpha - 1
   end if
   call system_clock(count_rate=rate)

   call ofk_set_mode(0_c_int); U2 = 0; call port(U2)   ! warm-up
   call ofk_set_mode(1_c_int); U3 = 0; call port(U3)
   tref = huge(1.0d0); ts = huge(1.0d0); tf = huge(1.0d0)
   do it = 1, reps
      U1 = 0; tref = min(tref, timed_ref())
      call ofk_set_mode(0_c_int); U2 = 0; ts = min(ts, timed_port(U2))
      call ofk_set_mode(1_c_int); U3 = 0; tf = min(tf, timed_port(U3))
   end do
   sc = max(maxval(abs(U1)), tiny(1.0_ReKi))
   rel_s = maxval(abs(U2 - U1)) / sc
   rel_f = maxval(abs(U3 - U1)) / sc
   write(*,'(a,3(",",i0),3(",",es11.4),2(",",f8.3),2(",",es9.2))') trim(routine), nCP, nSrc, reg, &
        tref, ts, tf, tref/ts, tref/tf, rel_s, rel_f

contains
   real(8) function timed_ref()
      integer(8) :: c0, c1
      call system_clock(c0)
      if (routine == 'seg') then
         call seg_ref(1, nCP, CPs, 1, nSrc, Pts, Connct, Gam, reg, Eps, U1)
      else
         call part_ref(nCP, CPs, nSrc, Pts, Alpha, reg, Eps, U1)
      end if
      call system_clock(c1)
      timed_ref = real(c1 - c0, 8)/rate
   end function timed_ref
   real(8) function timed_port(U)
      real(ReKi), intent(inout) :: U(:,:)
      integer(8) :: c0, c1
      call system_clock(c0)
      call port(U)
      call system_clock(c1)
      timed_port = real(c1 - c0, 8)/rate
   end function timed_port
   subroutine port(U)
      real(ReKi), intent(inout) :: U(:,:)
      if (routine == 'seg') then
         call seg_port(1, nCP, CPs, 1, nSrc, Pts, Connct, Gam, reg, Eps, U)
      else
         call part_port(nCP, CPs, nSrc, Pts, Alpha, reg, Eps, U)
      end if
   end subroutine port
end program perf_direct
