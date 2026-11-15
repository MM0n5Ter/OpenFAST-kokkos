!> L0 differential test for FVW_BiotSavart::ui_seg: the pin's Fortran (module renamed
!! FVW_BiotSavart_ref at build time) vs the port (module FVW_BiotSavart, Kokkos-backed),
!! in one process, through their real Fortran interfaces. Exit code 0 = pass.
!! Audit items covered (tests/biotsavart-direct/unit.md): every RegFunction 0..4, each guard
!! (denominator / norm2_orth / norm2_r0 vs PRECISION_UI, straddled), Rankine r_bar2<1,
!! LambOseen exp_value<MIN_EXP_VALUE, NaN guard style (Vatistas cycle vs nested if),
!! accumulation into Uind_out, untouched columns outside iCPStart:iCPEnd, "n++" tails of
!! every input, sub-ranges iCPStart>1 / iSegStart>1, 1-based SegConnct with depth/span rows
!! and mirror offset, non-contiguous sections, padded leading dimension, skinny shape.
program l0_ui_seg
   use NWTC_Library,       only: ReKi, IntKi
   use FVW_BiotSavart_ref, only: ui_ref => ui_seg, PRECISION_UI
   use FVW_BiotSavart,     only: ui_port => ui_seg
   use l0_util
   implicit none
   character(len=12), parameter :: names(0:4) = [character(len=12) :: 'None', 'Rankine', 'LambOseen', 'Vatistas', 'Offset']
   integer :: reg, icase, nCP, nSeg
   logical :: skinny_only
   character(len=16) :: buf

   call l0_init()
   call get_environment_variable('OFK_L0_SKINNY_ONLY', buf); skinny_only = (buf(1:1) == '1')

   if (.not. skinny_only) then
      do reg = 0, 4
         call case_random(reg, 10, 0, 'nCP = 0 (iCPEnd < iCPStart)', nCPoverride=0)
         call case_random(reg, 10, 0, 'nSeg = 0 (iSegEnd < iSegStart)', nSegoverride=0)
         call case_random(reg, 1, 1, 'one-to-one')
         call case_random(reg, 60, 80, 'CPs on segment points (guards)', kind=1)
         call case_random(reg, 60, 80, 'tiny and huge core radii', kind=2)
         call case_random(reg, 60, 80, 'zero-length segments (repeated index)', kind=3)
         call case_random(reg, 60, 80, 'sub-range iCPStart>1, iSegStart>1', kind=4)
         call case_random(reg, 61, 80, 'non-contiguous CPs/Uind_out sections', kind=5)
         call case_random(reg, 60, 80, 'padded SegPoints(4,:) and Uind_out(4,:)', kind=6)
         call case_random(reg, 60, 80, 'non-contiguous SegPoints(1:3,:) of (4,n)', kind=7)
         call case_random(reg, 40, 60, 'NaN point in one segment (guard style)', kind=8)
         call case_random(reg, 40, 60, 'negative zero in Uind_out', kind=9)
         call case_straddle(reg)
         call case_lattice(reg, .false.)
         call case_lattice(reg, .true.)
         call case_random(reg, 700, 900, 'medium square')
      end do
   end if
   do reg = 0, 4   ! skinny shape (lifting line: few targets, many sources)
      call case_random(reg, 150, 20000, 'skinny 150 x 20000')
      call case_random(reg, 7, 5000, 'skinny 7 x 5000, sub-range', kind=4)
   end do
   if (.not. skinny_only) then
      do icase = 1, ncases
         reg  = mod(icase, 5)
         nCP  = rand_int(1, 300); nSeg = rand_int(1, 500)
         call case_random(reg, nCP, nSeg, 'random', kind=rand_int(0, 7))
      end do
   end if

   call summary('l0_ui_seg', names)
   if (nfail > 0) error stop 1

contains

   !> Call A and B (see l0_util: oracle, port, strict, fast) on the same inputs and compare.
   !! stride=2 passes every other column of CPs and Uind_out (non-contiguous sections);
   !! ptsRows=3 with a (4,n) Pts passes the non-contiguous section Pts(1:3,:).
   !! The scale sum_j|c_j| comes from the oracle called once per segment.
   subroutine run(label, reg, iCPs, iCPe, CPs, iSs, iSe, Pts, Connct, Gam, Eps, U0, stride, ptsRows)
      character(*), intent(in) :: label
      integer,      intent(in) :: reg, iCPs, iCPe, iSs, iSe, stride, ptsRows
      real(ReKi),   intent(in) :: CPs(:,:), Pts(:,:), Gam(:), Eps(:), U0(:,:)
      integer(IntKi), intent(in) :: Connct(:,:)
      real(ReKi), allocatable :: Uref(:,:), Uport(:,:), S(:,:), Ut(:,:)
      logical, allocatable :: active(:)
      integer :: j, is
      allocate(Uref, source=U0); allocate(Uport, source=U0)
      allocate(S(size(U0,1), size(U0,2)), Ut(size(U0,1), size(U0,2)))
      allocate(active(size(U0,2)))
      call invoke_seg(roleA, Uref, iSs, iSe, reg, iCPs, iCPe, CPs, Pts, Connct, Gam, Eps, stride, ptsRows)
      call invoke_seg(roleB, Uport, iSs, iSe, reg, iCPs, iCPe, CPs, Pts, Connct, Gam, Eps, stride, ptsRows)
      S = 0.0_ReKi
      do is = iSs, iSe
         Ut = 0.0_ReKi
         call invoke_seg(ROLE_REF, Ut, is, is, reg, iCPs, iCPe, CPs, Pts, Connct, Gam, Eps, stride, ptsRows)
         S = S + abs(Ut)
      end do
      if (stride == 1) then
         do j = 1, size(U0,2); active(j) = (j >= iCPs .and. j <= iCPe); end do
      else
         do j = 1, size(U0,2); active(j) = (mod(j,2) == 1 .and. (j+1)/2 >= iCPs .and. (j+1)/2 <= iCPe); end do
      end if
      call compare(label, reg, U0, Uref, Uport, active, S)
   end subroutine run

   subroutine invoke_seg(role, U, s1, s2, reg, iCPs, iCPe, CPs, Pts, Connct, Gam, Eps, stride, ptsRows)
         integer,    intent(in)    :: role, s1, s2, reg, iCPs, iCPe, stride, ptsRows
         real(ReKi), intent(inout) :: U(:,:)
         real(ReKi), intent(in)    :: CPs(:,:), Pts(:,:), Gam(:), Eps(:)
         integer(IntKi), intent(in) :: Connct(:,:)
         if (use_ref(role)) then
            if (stride == 2) then
               call ui_ref (iCPs, iCPe, CPs(:,1::2), s1, s2, Pts, Connct, Gam, reg, Eps, U(:,1::2))
            else if (ptsRows == 3 .and. size(Pts,1) == 4) then
               call ui_ref (iCPs, iCPe, CPs, s1, s2, Pts(1:3,:), Connct, Gam, reg, Eps, U)
            else
               call ui_ref (iCPs, iCPe, CPs, s1, s2, Pts, Connct, Gam, reg, Eps, U)
            end if
         else
            if (stride == 2) then
               call ui_port(iCPs, iCPe, CPs(:,1::2), s1, s2, Pts, Connct, Gam, reg, Eps, U(:,1::2))
            else if (ptsRows == 3 .and. size(Pts,1) == 4) then
               call ui_port(iCPs, iCPe, CPs, s1, s2, Pts(1:3,:), Connct, Gam, reg, Eps, U)
            else
               call ui_port(iCPs, iCPe, CPs, s1, s2, Pts, Connct, Gam, reg, Eps, U)
            end if
         end if
   end subroutine invoke_seg

   !> Random cloud of nSeg segments over nSeg+1 points, nCP targets. All inputs have "n++"
   !! tails (EXTRA columns) filled with NaN or poison, so reading past the range is detected.
   subroutine case_random(reg, nCPin, nSegIn, label, kind, nCPoverride, nSegoverride)
      integer, intent(in) :: reg, nCPin, nSegIn
      character(*), intent(in) :: label
      integer, intent(in), optional :: kind, nCPoverride, nSegoverride
      integer, parameter :: EXTRA = 3
      integer :: k, nCP, nSeg, nPts, i, iCPs, iCPe, iSs, iSe, stride, ptsRows, nColU, rowsP, rowsU, nHead
      real(ReKi), allocatable :: CPs(:,:), Pts(:,:), Gam(:), Eps(:), U0(:,:)
      integer(IntKi), allocatable :: Connct(:,:)
      real(ReKi) :: r
      k = 0; if (present(kind)) k = kind
      nCP = max(nCPin, 1); nSeg = max(nSegIn, 1)
      nPts = nSeg + 1
      nHead = 0; if (k == 4) nHead = 5                    ! columns before iCPStart / iSegStart
      rowsP = 3; if (k == 6 .or. k == 7) rowsP = 4
      rowsU = 3; if (k == 6) rowsU = 4
      stride = 1; if (k == 5) stride = 2
      ptsRows = rowsP; if (k == 7) ptsRows = 3

      allocate(CPs(3, stride*(nHead+nCP+EXTRA)), Pts(rowsP, nPts+EXTRA), Gam(nHead+nSeg+EXTRA), Eps(nHead+nSeg+EXTRA))
      allocate(Connct(4, nHead+nSeg+EXTRA))
      nColU = size(CPs,2)
      allocate(U0(rowsU, nColU))
      call random_number(CPs); CPs = 10*CPs - 5
      call random_number(Pts); Pts = 10*Pts - 5
      Pts(:, nPts+1:) = qnan()                            ! unreferenced tail
      if (rowsP == 4) Pts(4,:) = qnan()                   ! row 4 is never read
      call random_number(Gam); Gam = 2*Gam - 1
      call random_number(Eps); Eps = 0.05_ReKi + Eps
      do i = 1, size(Connct,2)
         call random_number(r); Connct(1,i) = 1 + int(r*nPts)
         call random_number(r); Connct(2,i) = 1 + int(r*nPts)
         if (Connct(2,i) == Connct(1,i)) Connct(2,i) = 1 + mod(Connct(1,i), nPts)
         Connct(3,i) = -7; Connct(4,i) = 99999            ! depth/span rows: never used as indices
      end do
      Connct(1:2, :) = min(Connct(1:2, :), nPts)
      ! poison the segments outside the selected range: NaN circulation
      Gam(1:nHead) = qnan(); Gam(nHead+nSeg+1:) = qnan()
      Eps(1:nHead) = qnan(); Eps(nHead+nSeg+1:) = qnan()
      call random_number(U0); U0 = U0 - 0.5_ReKi          ! nonzero: the routine accumulates
      iCPs = 1 + nHead; iCPe = nHead + nCP
      iSs  = 1 + nHead; iSe  = nHead + nSeg
      if (present(nCPoverride)) iCPe = iCPs + nCPoverride - 1
      if (present(nSegoverride)) iSe = iSs + nSegoverride - 1
      ! targets outside the range are NaN, outputs outside the range are sentinels
      if (stride == 1) then
         CPs(:, 1:iCPs-1) = qnan(); CPs(:, iCPe+1:) = qnan()
         U0(:, 1:iCPs-1) = SENTINEL; U0(:, iCPe+1:) = SENTINEL
      else
         U0(:, 2::2) = SENTINEL; U0(:, 2*iCPe+1:) = SENTINEL
      end if
      if (rowsU == 4) U0(4,:) = SENTINEL

      select case (k)
      case (1)   ! CPs exactly on segment end points and on segment midpoints
         do i = iCPs, iCPe, 2
            CPs(:, i) = Pts(1:3, Connct(1 + mod(i,2), iSs + mod(i, nSeg)))
         end do
         do i = iCPs+1, iCPe, 4
            CPs(:, i) = 0.5_ReKi*(Pts(1:3, Connct(1, iSs + mod(i, nSeg))) + Pts(1:3, Connct(2, iSs + mod(i, nSeg))))
         end do
      case (2)   ! core radii from 1e-6 to 1e3 (Rankine both sides of r_bar2=1, LambOseen exp cut)
         Eps(iSs:iSe:2) = 1.0e-6_ReKi
         Eps(iSs+1:iSe:2) = 1.0e+3_ReKi
         Eps(iSs+2:iSe:5) = 0.0_ReKi
      case (3)   ! zero-length segments
         Connct(2, iSs:iSe:3) = Connct(1, iSs:iSe:3)
      case (8)   ! one NaN point used by one segment: Vatistas (`<=`/cycle) propagates it, the others skip it
         Pts(:, Connct(1, iSs)) = qnan()
      case (9)   ! -0.0 in Uind_out: Fortran adds +0.0 when every segment is skipped
         U0(1:3, iCPs:iCPe) = -0.0_ReKi
         Gam(iSs:iSe) = 0.0_ReKi
      end select
      if (k == 7) then
         call run(label, reg, iCPs, iCPe, CPs, iSs, iSe, Pts, Connct, Gam, Eps, U0, stride, 3)
      else
         call run(label, reg, iCPs, iCPe, CPs, iSs, iSe, Pts, Connct, Gam, Eps, U0, stride, rowsP)
      end if
   end subroutine case_random

   !> Inputs straddling every threshold: PRECISION_UI on norm2_orth, denominator and norm2_r0;
   !! Rankine r_bar2 = 1; LambOseen exp_value = MIN_EXP_VALUE.
   subroutine case_straddle(reg)
      integer, intent(in) :: reg
      integer, parameter :: NS = 41
      real(ReKi) :: CPs(3, 4*NS), Pts(3, 2+2*NS), Gam(1+NS), Eps(1+NS), U0(3, 4*NS)
      integer(IntKi) :: Connct(4, 1+NS)
      real(ReKi) :: f, d, rc
      integer :: i
      Pts(:,1) = [0.0_ReKi, 0.0_ReKi, 0.0_ReKi]; Pts(:,2) = [1.0_ReKi, 0.0_ReKi, 0.0_ReKi]
      Connct(:,1) = [1, 2, 1, 1]; Gam(1) = 1.0_ReKi; Eps(1) = 1.0e-3_ReKi
      do i = 1, NS
         f = 0.5_ReKi + real(i-1, ReKi)/real(NS-1, ReKi)               ! 0.5 .. 1.5
         d = sqrt(PRECISION_UI) * f
         CPs(:, i)        = [0.5_ReKi, d, 0.0_ReKi]                     ! norm2_orth = d^2 ~ PRECISION_UI
         CPs(:, NS+i)     = [0.0_ReKi, 0.0_ReKi, d]                     ! denominator ~ d^2 near point 1
         rc = 1.0e-3_ReKi
         if (i <= NS/2) then
            CPs(:, 2*NS+i) = [0.5_ReKi, rc*(1.0_ReKi + (f-1.0_ReKi)*1.0e-12_ReKi), 0.0_ReKi]  ! Rankine r_bar2 ~ 1
         else
            CPs(:, 2*NS+i) = [0.5_ReKi, rc*sqrt(10.0_ReKi/1.25643_ReKi)*(1.0_ReKi + (f-1.0_ReKi)*1.0e-12_ReKi), 0.0_ReKi]
         end if
         CPs(:, 3*NS+i)   = [0.3_ReKi, 0.7_ReKi, -0.2_ReKi]
         ! short segments of length d ~ sqrt(PRECISION_UI): norm2_r0 straddles PRECISION_UI
         Pts(:, 2+2*i-1) = [2.0_ReKi, 0.0_ReKi, 0.0_ReKi]
         Pts(:, 2+2*i)   = [2.0_ReKi + d, 0.0_ReKi, 0.0_ReKi]
         Connct(:, 1+i)  = [2+2*i-1, 2+2*i, 1, 1+i]
         Gam(1+i) = 1.0_ReKi; Eps(1+i) = 1.0e-3_ReKi
      end do
      U0 = 0.25_ReKi
      call run('straddle thresholds (PRECISION_UI, r_bar2=1, exp=-10)', reg, 1, 4*NS, CPs, 1, 1+NS, Pts, Connct, Gam, Eps, &
               U0, 1, 3)
   end subroutine case_straddle

   !> A realistic wake: one helical vortex lattice (nSpan x nDepth panels) packed into segments
   !! like LatticeToSegments (spanwise + streamwise, Gamma differences, depth/span rows), with
   !! optional ground mirror (indices + nSegP, -Gamma). Targets = all lattice points (wake on
   !! wake, so every target sits on segment end points) + lifting-line midpoints.
   subroutine case_lattice(reg, mirror)
      integer, intent(in) :: reg
      logical, intent(in) :: mirror
      integer, parameter :: nS = 12, nD = 40
      real(ReKi), allocatable :: Pts(:,:), CPs(:,:), Gam(:), Eps(:), U0(:,:), GamPanel(:,:)
      integer(IntKi), allocatable :: Connct(:,:)
      integer :: iS, iD, ip, iseg, nP, nSegL, nM, nCPl
      real(ReKi) :: psi, rr
      nP = (nS+1)*(nD+1)
      nM = 1; if (mirror) nM = 2
      nSegL = nS*(nD+1) + (nS+1)*nD
      allocate(Pts(3, nM*nP), Connct(4, nM*nSegL), Gam(nM*nSegL), Eps(nM*nSegL), GamPanel(nS, nD))
      do iD = 0, nD
         do iS = 0, nS
            ip = 1 + iS + (nS+1)*iD
            psi = 0.15_ReKi*iD
            rr = 1.0_ReKi + 3.0_ReKi*iS/real(nS, ReKi) + 0.002_ReKi*iD*iS
            Pts(:, ip) = [0.3_ReKi*iD, rr*cos(psi), rr*sin(psi) + 4.0_ReKi]
         end do
      end do
      do iD = 1, nD
         do iS = 1, nS
            GamPanel(iS, iD) = sin(3.1_ReKi*iS/nS) * (1.0_ReKi + 0.01_ReKi*iD)
         end do
      end do
      iseg = 0
      do iD = 0, nD          ! spanwise segments (shed vorticity)
         do iS = 1, nS
            iseg = iseg + 1
            Connct(:, iseg) = [1 + (iS-1) + (nS+1)*iD, 1 + iS + (nS+1)*iD, iD+1, iS]
            if (iD == 0) then
               Gam(iseg) = GamPanel(iS, 1)
            else if (iD == nD) then
               Gam(iseg) = -GamPanel(iS, nD)
            else
               Gam(iseg) = GamPanel(iS, iD+1) - GamPanel(iS, iD)
            end if
         end do
      end do
      do iD = 1, nD          ! streamwise segments (trailed vorticity)
         do iS = 0, nS
            iseg = iseg + 1
            Connct(:, iseg) = [1 + iS + (nS+1)*(iD-1), 1 + iS + (nS+1)*iD, iD, iS+1]
            if (iS == 0) then
               Gam(iseg) = -GamPanel(1, iD)
            else if (iS == nS) then
               Gam(iseg) = GamPanel(nS, iD)
            else
               Gam(iseg) = GamPanel(iS, iD) - GamPanel(iS+1, iD)
            end if
         end do
      end do
      Eps(1:nSegL) = 0.05_ReKi + 0.002_ReKi*Connct(3, 1:nSegL)   ! core grows with depth
      if (mirror) then       ! ground mirror: z -> -z, indices + nSegP, opposite circulation
         Pts(:, nP+1:2*nP) = Pts(:, 1:nP); Pts(3, nP+1:2*nP) = -Pts(3, 1:nP)
         Connct(:, nSegL+1:2*nSegL) = Connct(:, 1:nSegL)
         Connct(1:2, nSegL+1:2*nSegL) = Connct(1:2, 1:nSegL) + nP
         Gam(nSegL+1:2*nSegL) = -Gam(1:nSegL)
         Eps(nSegL+1:2*nSegL) = Eps(1:nSegL)
      end if
      nCPl = nP + nS + 2
      allocate(CPs(3, nCPl), U0(3, nCPl))
      CPs(:, 1:nP) = Pts(:, 1:nP)
      do iS = 1, nS
         CPs(:, nP+iS) = 0.5_ReKi*(Pts(:, iS) + Pts(:, iS+1))
      end do
      CPs(:, nCPl-1:nCPl) = qnan()                          ! "nCPs++" tail
      call random_number(U0); U0(:, nCPl-1:nCPl) = SENTINEL
      if (mirror) then
         call run('lattice + ground mirror (wake on wake)', reg, 1, nP+nS, CPs, 1, 2*nSegL, Pts, Connct, Gam, Eps, U0, 1, 3)
      else
         call run('lattice (wake on wake)', reg, 1, nP+nS, CPs, 1, nSegL, Pts, Connct, Gam, Eps, U0, 1, 3)
      end if
   end subroutine case_lattice

end program l0_ui_seg
