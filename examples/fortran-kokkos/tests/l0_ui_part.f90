!> L0 differential test: frozen Fortran (module renamed to toy_biotsavart_ref at
!! build time) vs the port (module toy_biotsavart, backed by Kokkos), in one process,
!! on randomized inputs plus edge cases. Exit code 0 = pass.
!!
!! Environment: OFK_L0_SEED (default 1), OFK_L0_CASES (default 50),
!!              OFK_L0_RTOL (default 1e-12; use ~1e-10 for GPU backends).
program l0_ui_part
   use toy_kinds,          only: ReKi, IntKi
   use toy_biotsavart_ref, only: ui_ref  => ui_part
   use toy_biotsavart,     only: ui_port => ui_part
   implicit none
   integer, parameter :: EXTRA = 3                    ! "nCPs++": columns the routine must not touch
   real(ReKi), parameter :: SENTINEL = -12345.0_ReKi
   real(ReKi) :: rtol, worst
   integer :: seed, ncases, icase, nfail, nCPs, nPart, kind
   character(len=64) :: buf

   seed = 1; ncases = 50; rtol = 1.0e-12_ReKi
   call get_environment_variable('OFK_L0_SEED',  buf); if (len_trim(buf) > 0) read(buf,*) seed
   call get_environment_variable('OFK_L0_CASES', buf); if (len_trim(buf) > 0) read(buf,*) ncases
   call get_environment_variable('OFK_L0_RTOL',  buf); if (len_trim(buf) > 0) read(buf,*) rtol
   call seed_rng(seed)

   nfail = 0; worst = 0
   ! Fixed edge cases first, then random sizes.
   call run_case(0,   10, 0, 'nCPs = 0')
   call run_case(10,   0, 0, 'nPart = 0')
   call run_case(1,    1, 0, 'one-to-one')
   call run_case(40,  40, 1, 'CPs on top of particles (singularity guard)')
   call run_case(40,  60, 2, 'tiny and huge core radii')
   call run_case(40,  60, 3, 'non-contiguous array sections')
   do icase = 1, ncases
      nCPs = rand_int(1, 300); nPart = rand_int(1, 500); kind = mod(icase, 4)
      call run_case(nCPs, nPart, kind, 'random')
   end do

   write(*,'(a,i0,a,i0,a,es9.2,a,es9.2,a)') 'l0_ui_part: ', ncases + 6 - nfail, ' passed, ', nfail, &
        ' failed; worst rel err ', worst, ' (rtol ', rtol, ')'
   if (nfail > 0) error stop 1

contains

   subroutine run_case(nCPs, nPart, kind, label)
      integer, intent(in) :: nCPs, nPart, kind
      character(*), intent(in) :: label
      real(ReKi), allocatable :: CPs(:,:), Part(:,:), Alpha(:,:), Reg(:), U0(:,:), Uref(:,:), Uport(:,:)
      real(ReKi) :: err, scale
      integer :: i
      allocate(CPs(3, nCPs+EXTRA), Part(3, nPart+EXTRA), Alpha(3, nPart+EXTRA), Reg(nPart+EXTRA))
      allocate(U0(3, nCPs+EXTRA))
      call random_number(CPs);  CPs  = 10*CPs  - 5
      call random_number(Part); Part = 10*Part - 5
      call random_number(Alpha); Alpha = 2*Alpha - 1
      call random_number(Reg);  Reg = 0.05_ReKi + Reg
      call random_number(U0);   U0 = U0 - 0.5_ReKi       ! nonzero: the routine accumulates
      U0(:, nCPs+1:) = SENTINEL
      select case (kind)
      case (1)
         do i = 1, min(nCPs, nPart), 2
            CPs(:, i) = Part(:, i)
         end do
      case (2)
         if (nPart > 1) then
            Reg(1:nPart:2) = 1.0e-6_ReKi
            Reg(2:nPart:2) = 1.0e+3_ReKi
         end if
      end select
      allocate(Uref(3, size(U0,2)), Uport(3, size(U0,2)))
      Uref = U0; Uport = U0

      if (kind == 3 .and. nCPs > 1) then
         ! Pass every other column: the port's explicit-shape interface gets a copy.
         call ui_ref (int((nCPs+1)/2, IntKi), CPs(:, 1:nCPs:2), int(nPart, IntKi), Part, Alpha, Reg, Uref(:, 1:nCPs:2))
         call ui_port(int((nCPs+1)/2, IntKi), CPs(:, 1:nCPs:2), int(nPart, IntKi), Part, Alpha, Reg, Uport(:, 1:nCPs:2))
      else
         call ui_ref (int(nCPs, IntKi), CPs, int(nPart, IntKi), Part, Alpha, Reg, Uref)
         call ui_port(int(nCPs, IntKi), CPs, int(nPart, IntKi), Part, Alpha, Reg, Uport)
      end if

      scale = max(maxval(abs(Uref)), tiny(1.0_ReKi))
      err = maxval(abs(Uport - Uref)) / scale
      worst = max(worst, err)
      if (err > rtol .or. any(Uport(:, nCPs+1:) /= SENTINEL)) then
         nfail = nfail + 1
         write(*,'(a,a,a,i0,a,i0,a,es10.3)') 'FAIL [', label, '] nCPs=', nCPs, ' nPart=', nPart, ' rel err=', err
         if (any(Uport(:, nCPs+1:) /= SENTINEL)) write(*,'(a)') '     wrote past column nCPs'
      end if
   end subroutine run_case

   subroutine seed_rng(s)
      integer, intent(in) :: s
      integer :: n
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n)); put = s + 37*[(n - 1, n = 1, size(put))]
      call random_seed(put=put)
   end subroutine seed_rng

   integer function rand_int(lo, hi)
      integer, intent(in) :: lo, hi
      real :: r
      call random_number(r)
      rand_int = lo + int(r * (hi - lo + 1))
      rand_int = min(rand_int, hi)
   end function rand_int

end program l0_ui_part
