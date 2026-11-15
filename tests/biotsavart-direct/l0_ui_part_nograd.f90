!> L0 differential test for FVW_BiotSavart::ui_part_nograd: the pin's Fortran (module renamed
!! FVW_BiotSavart_ref) vs the port (module FVW_BiotSavart, Kokkos-backed). Exit code 0 = pass.
!! Audit items covered (tests/biotsavart-direct/unit.md): RegFunction 0..2, MINNORM guard
!! straddled at BOTH candidate values (1e-4 as double and as single-precision literal), Exp
!! r^3 = 8 rc^3 cut straddled, Compact r^2 = rc^2 and RegParam = 0 clamps, accumulation,
!! "n++" tails (NaN poisoned), non-contiguous sections, padded leading dimension, NaN points.
program l0_ui_part_nograd
   use NWTC_Library,       only: ReKi, IntKi
   use FVW_BiotSavart_ref, only: ui_ref => ui_part_nograd, MINNORM, PART_REG_C2
   use FVW_BiotSavart,     only: ui_port => ui_part_nograd
   use l0_util
   implicit none
   character(len=12), parameter :: names(0:2) = [character(len=12) :: 'None', 'Exp', 'Compact']
   integer :: reg, icase

   call l0_init()
   write(*,'(a,es24.17,a,es24.17)') 'MINNORM as compiled = ', MINNORM, '   MINNORM**2 = ', MINNORM*MINNORM
   do reg = 0, 2
      call case_random(reg, 10, 0, 'nPart = 0')
      call case_random(reg, 0, 10, 'nCPs = 0')
      call case_random(reg, 1, 1, 'one-to-one')
      call case_random(reg, 60, 80, 'CPs on particles (r2 = 0)', kind=1)
      call case_random(reg, 60, 80, 'tiny, huge and zero core radii', kind=2)
      call case_random(reg, 61, 80, 'non-contiguous CPs/UIout sections', kind=3)
      call case_random(reg, 60, 80, 'padded (4,n) arrays', kind=4)
      call case_random(reg, 40, 60, 'NaN particle', kind=5)
      call case_straddle(reg)
      call case_random(reg, 700, 900, 'medium square')
   end do
   do icase = 1, ncases
      reg = mod(icase, 3)
      call case_random(reg, rand_int(1, 300), rand_int(1, 500), 'random', kind=rand_int(0, 4))
   end do
   call summary('l0_ui_part_nograd', names)
   if (nfail > 0) error stop 1

contains

   !> Call A and B (see l0_util) on the same inputs and compare; stride=2 passes every other
   !! column of every array. The scale sum_j|c_j| comes from the oracle called per particle.
   subroutine run(label, reg, nCPs, CPs, nPart, Part, Alpha, Eps, U0, stride)
      character(*), intent(in) :: label
      integer,      intent(in) :: reg, nCPs, nPart, stride
      real(ReKi),   intent(in) :: CPs(:,:), Part(:,:), Alpha(:,:), Eps(:), U0(:,:)
      real(ReKi), allocatable :: Uref(:,:), Uport(:,:), S(:,:), Ut(:,:)
      logical, allocatable :: active(:)
      integer :: j, ip
      allocate(Uref, source=U0); allocate(Uport, source=U0); allocate(active(size(U0,2)))
      allocate(S(size(U0,1), size(U0,2)), Ut(size(U0,1), size(U0,2)))
      call invoke_part(roleA, Uref, 1, nPart, reg, nCPs, CPs, Part, Alpha, Eps, stride)
      call invoke_part(roleB, Uport, 1, nPart, reg, nCPs, CPs, Part, Alpha, Eps, stride)
      S = 0.0_ReKi
      do ip = 1, nPart
         Ut = 0.0_ReKi
         call invoke_part(ROLE_REF, Ut, ip, ip, reg, nCPs, CPs, Part, Alpha, Eps, stride)
         S = S + abs(Ut)
      end do
      if (stride == 1) then
         do j = 1, size(U0,2); active(j) = (j <= nCPs); end do
      else
         do j = 1, size(U0,2); active(j) = (mod(j,2) == 1 .and. (j+1)/2 <= nCPs); end do
      end if
      call compare(label, reg, U0, Uref, Uport, active, S)
   end subroutine run

   !> Particles p1..p2 (p1 = 1 for the full call: the arrays are passed whole, with tails).
   subroutine invoke_part(role, U, p1, p2, reg, nCPs, CPs, Part, Alpha, Eps, stride)
         integer,    intent(in)    :: role, p1, p2, reg, nCPs, stride
         real(ReKi), intent(inout) :: U(:,:)
         real(ReKi), intent(in)    :: CPs(:,:), Part(:,:), Alpha(:,:), Eps(:)
         integer :: q, n
         n = p2 - p1 + 1
         q = stride*(p1 - 1) + 1     ! first column of particle p1
         if (use_ref(role)) then
            if (stride == 1) then
               call ui_ref (nCPs, CPs, n, Part(:,q:), Alpha(:,q:), reg, Eps(q:), U)
            else
               call ui_ref (nCPs, CPs(:,1::2), n, Part(:,q::2), Alpha(:,q::2), reg, Eps(q::2), U(:,1::2))
            end if
         else
            if (stride == 1) then
               call ui_port(nCPs, CPs, n, Part(:,q:), Alpha(:,q:), reg, Eps(q:), U)
            else
               call ui_port(nCPs, CPs(:,1::2), n, Part(:,q::2), Alpha(:,q::2), reg, Eps(q::2), U(:,1::2))
            end if
         end if
   end subroutine invoke_part

   subroutine case_random(reg, nCPs, nPart, label, kind)
      integer, intent(in) :: reg, nCPs, nPart
      character(*), intent(in) :: label
      integer, intent(in), optional :: kind
      integer, parameter :: EXTRA = 3
      integer :: k, stride, rows, i
      real(ReKi), allocatable :: CPs(:,:), Part(:,:), Alpha(:,:), Eps(:), U0(:,:)
      k = 0; if (present(kind)) k = kind
      stride = 1; if (k == 3) stride = 2
      rows = 3; if (k == 4) rows = 4
      allocate(CPs(rows, stride*(nCPs+EXTRA)), U0(rows, stride*(nCPs+EXTRA)))
      allocate(Part(rows, stride*(nPart+EXTRA)), Alpha(rows, stride*(nPart+EXTRA)), Eps(stride*(nPart+EXTRA)))
      call random_number(CPs);   CPs = 10*CPs - 5
      call random_number(Part);  Part = 10*Part - 5
      call random_number(Alpha); Alpha = 2*Alpha - 1
      call random_number(Eps);   Eps = 0.05_ReKi + Eps
      call random_number(U0);    U0 = U0 - 0.5_ReKi
      ! "n++" tails: NaN inputs, sentinel outputs
      CPs(:, stride*nCPs+1:) = qnan(); Part(:, stride*nPart+1:) = qnan()
      Alpha(:, stride*nPart+1:) = qnan(); Eps(stride*nPart+1:) = qnan()
      U0(:, stride*nCPs+1:) = SENTINEL
      if (stride == 2) U0(:, 2::2) = SENTINEL
      if (rows == 4) then
         CPs(4,:) = qnan(); Part(4,:) = qnan(); Alpha(4,:) = qnan(); U0(4,:) = SENTINEL
      end if
      select case (k)
      case (1)
         do i = 1, min(nCPs, nPart), 2
            CPs(1:3, i) = Part(1:3, i)
         end do
      case (2)
         if (nPart > 2) then
            Eps(1:nPart:3) = 1.0e-6_ReKi
            Eps(2:nPart:3) = 1.0e+3_ReKi
            Eps(3:nPart:3) = 0.0_ReKi
         end if
      case (5)
         if (nPart > 0) Part(1:3, 1) = qnan()
      end select
      call run(label, reg, nCPs, CPs, nPart, Part, Alpha, Eps, U0, stride)
   end subroutine case_random

   !> Thresholds: r^2 around MINNORM**2 for MINNORM = 1d-4 AND MINNORM = real(1e-4,4)
   !! (so the case is meaningful whatever -fdefault-real-8 did to the literal), the Exp cut
   !! r^3 = 8 rc^3 and the Compact edge r^2 = (1.6 rc)^2.
   subroutine case_straddle(reg)
      integer, intent(in) :: reg
      integer, parameter :: NS = 41
      real(ReKi) :: CPs(3, 4*NS), Part(3, 1), Alpha(3, 1), Eps(1), U0(3, 4*NS)
      real(ReKi) :: f, rc, mA, mB, d
      integer :: i
      Part(:,1) = 0.0_ReKi; Alpha(:,1) = [0.3_ReKi, -0.7_ReKi, 0.5_ReKi]
      rc = 0.2_ReKi; Eps(1) = rc
      mA = 1.0e-4_ReKi                       ! double literal
      mB = real(1.0e-4_4, ReKi)              ! single-precision literal widened (9.9999997473787516e-05)
      do i = 1, NS
         f = (real(i-1, ReKi)/real(NS-1, ReKi) - 0.5_ReKi) * 4.0e-8_ReKi   ! relative offsets -2e-8 .. 2e-8
         CPs(:, i)      = [mA*(1.0_ReKi + f), 0.0_ReKi, 0.0_ReKi]
         CPs(:, NS+i)   = [0.0_ReKi, mB*(1.0_ReKi + f), 0.0_ReKi]
         d = 2.0_ReKi*rc*(1.0_ReKi + f*1.0e-6_ReKi)                        ! r = 2 rc: r^3 vs 8 rc^3
         CPs(:, 2*NS+i) = [d/sqrt(3.0_ReKi), d/sqrt(3.0_ReKi), d/sqrt(3.0_ReKi)]
         d = PART_REG_C2*rc*(1.0_ReKi + f*1.0e-6_ReKi)                     ! Compact edge
         CPs(:, 3*NS+i) = [0.0_ReKi, 0.0_ReKi, d]
      end do
      ! the gap between the two MINNORM**2 values, exactly in the middle
      CPs(:, NS/2) = [sqrt(0.5_ReKi*(mA*mA + mB*mB)), 0.0_ReKi, 0.0_ReKi]
      U0 = 0.125_ReKi
      call run('straddle thresholds (MINNORM x2, 2rc, 1.6rc)', reg, 4*NS, CPs, 1, Part, Alpha, Eps, U0, 1)
   end subroutine case_straddle

end program l0_ui_part_nograd
