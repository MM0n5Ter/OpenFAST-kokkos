!> L0 for unit treecode: the pin's grow_tree_part / ui_tree_part (FVW_VortexTools_ref, built with
!! the pin's flags and OPENMP=OFF) and the pin's direct sum ui_part_nograd (FVW_BiotSavart_ref)
!! against the port (C API of the induced-velocities handle, the same calls as the FVW_Subs shims),
!! in one process. Every case runs the port in BOTH modes:
!!   claim 1/2 : strict vs the Fortran tree, max_k |S-F| / sum|c| per target  <= OFK_L0_RTOL
!!               (0 = bitwise, including -0, NaN tails, untouched columns; CPU)
!!   R2        : E(fast) <= max(1.1 * E(Fortran tree), OFK_L0_R2FLOOR), E(X) = max|X - D| / max|D|,
!!               D = the direct sum, X and D from zero (criteria.md R2)
!!   claim 3   : fast vs strict, max |K-S| / sum|c| (informational: the tree is class R2, a
!!               rounding-level change of a node center may flip one acceptance test)
!! Cases (unit.md audit list): RegFunctionPart 0..2; BranchFactor 1.0 1.5 2.0 3.0; DistanceDirect
!! 0, small, 2.2 x size (the LL value); clouds: uniform box, helical wake filaments, regular grid
!! (ties with cell centers: `>` vs `>=`), duplicated points, a clump of identical points, clumps
!! straddling the degenerate-extent threshold 1.0e-5 by 1e-9 relative, zero-vorticity particles (dead nodes), all
!! zero, tiny vorticity straddling EqualRealNos' tolerance, large eps (radius <= maxRegParam leaf
!! rule), eps = 0; nPart = 0, 1, 2, 3, 31, 32, 33 (TREE_LEAF_NCRIT), 65, 500, 3000, 20000; targets
!! on particles, inside, far; nonzero initial Uind with -0, NaN tail columns, leading dimension 4.
!! P6: a fixed-size grow/eval sequence after warm-up allocates nothing (ofk_iv_counters).
!! Environment: OFK_L0_SEED (1), OFK_L0_CASES (40 random), OFK_L0_RTOL (0), OFK_L0_R2FLOOR (1e-12),
!!   OFK_L0_DUMP=<file> (strict and fast outputs, R7 across processes/threads), OFK_L0_INFO=1,
!!   OFK_L0_BIG=0 (skip the 20000-particle cases).
program l0_tc
   use NWTC_Library, only: ReKi, IntKi, ErrMsgLen
   use FVW_VortexTools_ref, only: T_Tree, grow_tree_part, ui_tree_part, cut_tree
   use FVW_BiotSavart_ref, only: ui_part_nograd, MINNORM
   use tc_api
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
   use, intrinsic :: iso_c_binding, only: c_int, c_long_long
   implicit none
   real(ReKi), parameter :: Pi = 3.14159265358979323846_ReKi  ! (NWTC_Library%Pi is set by NWTC_Init only)
   integer, target :: key = 0
   integer :: seed = 1, ncases = 40, nrun = 0, nfail = 0, dumpunit = 0, nbig = 1
   real(ReKi) :: rtol = 0.0_ReKi, r2floor = 1.0e-12_ReKi
   logical :: info = .false., dumping = .false.
   real(ReKi) :: worst_c1 = 0, worst_c3 = 0, worst_r2ratio = 0, worst_EF = 0, worst_EK = 0
   integer :: nbit_c1 = 0, nbit_c3 = 0, nflip = 0, nsing = 0
   real(ReKi), parameter :: BFS(4) = [1.0_ReKi, 1.5_ReKi, 2.0_ReKi, 3.0_ReKi]
   integer :: k, kind, n, reg, ib, id, icase
   integer, parameter :: NS(13) = [0, 1, 2, 3, 31, 32, 33, 65, 200, 500, 3000, 20000, 20000]
   character(len=1024) :: buf

   call get_environment_variable('OFK_L0_SEED',  buf); if (len_trim(buf) > 0) read(buf,*) seed
   call get_environment_variable('OFK_L0_CASES', buf); if (len_trim(buf) > 0) read(buf,*) ncases
   call get_environment_variable('OFK_L0_RTOL',  buf); if (len_trim(buf) > 0) read(buf,*) rtol
   call get_environment_variable('OFK_L0_R2FLOOR', buf); if (len_trim(buf) > 0) read(buf,*) r2floor
   call get_environment_variable('OFK_L0_INFO', buf); info = (buf(1:1) == '1')
   call get_environment_variable('OFK_L0_BIG', buf); if (len_trim(buf) > 0) read(buf,*) nbig
   call get_environment_variable('OFK_L0_DUMP', buf)
   if (len_trim(buf) > 0) then
      open(newunit=dumpunit, file=trim(buf), access='stream', form='unformatted', status='replace', action='write')
      dumping = .true.
   end if
   call seed_rng(seed)

   ! ---- fixed cases: every cloud kind x size x regularization; BF/DD cycled
   k = 0
   do kind = 1, 11
      do n = 1, size(NS)
         ! the 20000-particle clouds take the parallel subtree build (kSplitMinN = 8192)
         if (NS(n) >= 20000 .and. (nbig == 0 .or. .not. any(kind == [1, 2, 3, 4, 8, 10]))) cycle
         if (kind >= 6 .and. NS(n) > 3000 .and. NS(n) < 20000) cycle
         do reg = 0, 2
            k = k + 1
            ib = 1 + mod(k, 4); id = mod(k / 4, 3)
            call one_case(kind, NS(n), reg, BFS(ib), id)
         end do
      end do
   end do
   ! every BF x DD on the wake cloud, all regularizations
   do ib = 1, 4
      do id = 0, 2
         do reg = 0, 2
            call one_case(2, 3000, reg, BFS(ib), id)
         end do
      end do
   end do
   ! ---- random cases
   do icase = 1, ncases
      call one_case(rand_int(1, 11), rand_int(0, 2500), rand_int(0, 2), BFS(rand_int(1, 4)), rand_int(0, 2))
   end do
   call p6_sequence()

   write(*,'(a,i0,a,i0,a,es9.2,a,es9.2)') 'l0_tc: seed=', seed, ' random cases=', ncases, ' rtol=', rtol, &
        ' r2floor=', r2floor
   write(*,'(a,es10.3,a,i0)') '  claim 1/2 strict vs Fortran tree: worst err=', worst_c1, '  cases with bit diffs=', nbit_c1
   write(*,'(a,es10.3,a,es10.3,a,f8.4,a,i0)') '  R2: worst E(Fortran tree)=', worst_EF, ' worst E(fast)=', worst_EK, &
        ' worst E(fast)/max(1.1E(F),floor)=', worst_r2ratio, '  Fortran-singular cases (excluded, port singular too)=', nsing
   write(*,'(a,es10.3,a,i0,a,i0)') '  claim 3 fast vs strict (info): worst=', worst_c3, '  cases with bit diffs=', nbit_c3, &
        '  cases with > 1e-12 (acceptance flips)=', nflip
   write(*,'(a,i0,a,i0,a)') 'l0_tc: ', nrun - nfail, ' passed, ', nfail, ' failed'
   if (dumping) close(dumpunit)
   call ofk_iv_free(key)
   if (nfail > 0 .and. .not. info) error stop 1

contains

   real(ReKi) function rnd()
      call random_number(rnd)
   end function rnd
   integer function rand_int(lo, hi)
      integer, intent(in) :: lo, hi
      rand_int = min(lo + int(rnd() * (hi - lo + 1)), hi)
   end function rand_int
   subroutine seed_rng(s)
      integer, intent(in) :: s
      integer :: m, j
      integer, allocatable :: put(:)
      call random_seed(size=m)
      allocate(put(m)); put = s + 37*[(j - 1, j = 1, m)]
      call random_seed(put=put)
   end subroutine seed_rng
   real(ReKi) function qnan()
      qnan = ieee_value(1.0_ReKi, ieee_quiet_nan)
   end function qnan
   logical function same_bits(a, b)
      real(ReKi), intent(in) :: a, b
      if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
         same_bits = ieee_is_nan(a) .and. ieee_is_nan(b)
      else
         same_bits = (a == b) .and. (sign(1.0_ReKi, a) == sign(1.0_ReKi, b))
      end if
   end function same_bits

   !> Particle clouds. kind: 1 box, 2 helical wake, 3 regular grid, 4 duplicated points, 5 identical
   !! clump, 6 clump extent 0.9e-5, 7 clump extent 1.1e-5, 8 zero vorticity, 9 tiny vorticity,
   !! 10 large eps, 11 eps = 0.
   subroutine cloud(kind, n, P, A, R, L)
      integer,    intent(in)  :: kind, n
      real(ReKi), intent(out) :: P(:,:), A(:,:), R(:), L
      integer :: i, b, nb, m
      real(ReKi) :: s, th, t(3), ext
      L = 10.0_ReKi
      do i = 1, n
         P(:,i) = L*[rnd(), rnd(), rnd()]
         A(:,i) = [2*rnd()-1, 2*rnd()-1, 2*rnd()-1]
         R(i)   = 0.05_ReKi + 0.5_ReKi*rnd()
      end do
      select case (kind)
      case (2)   ! three helical filaments, alpha along the tangent (a rotor wake)
         nb = 3; L = 20.0_ReKi
         do i = 1, n
            b = mod(i-1, nb); m = (i-1)/nb
            s = real(m, ReKi) / real(max(n/nb, 1), ReKi)
            th = 2*Pi*4*s + b*2*Pi/nb
            P(:,i) = [L*s, 5*cos(th), 5*sin(th)] + 1.0e-3_ReKi*[rnd(), rnd(), rnd()]
            t = [L/(2*Pi*4), -5*sin(th), 5*cos(th)]
            A(:,i) = (0.5_ReKi + rnd()) * t / sqrt(sum(t**2)) * (L / max(n/nb, 1))
            R(i) = 0.1_ReKi + 0.3_ReKi*s
         end do
      case (3)   ! regular grid: coordinates on cell boundaries and centers
         m = max(1, nint(real(n, ReKi)**(1.0_ReKi/3.0_ReKi)))
         do i = 1, n
            P(:,i) = 0.5_ReKi*[real(mod(i-1, m), ReKi), real(mod((i-1)/m, m), ReKi), real((i-1)/(m*m), ReKi)]
         end do
      case (4)   ! every point twice
         do i = 2, n, 2
            P(:,i) = P(:,i-1)
         end do
      case (5)   ! a clump of identical points (40) among the rest
         do i = 1, min(n, 40)
            P(:,i) = [3.0_ReKi, 4.0_ReKi, 5.0_ReKi]
         end do
      case (6, 7)   ! a clump of 40 with extent just below / above the threshold 1.0e-5, by 1e-9
         ! relative: tight enough to tell the double 1.0e-5 from float(1.0e-5) = 1.0e-5 - 2.5e-14
         ! (mutation check: a hard-coded 1.0e-5f survived the earlier 0.9e-5 / 1.1e-5 clumps)
         ext = merge(1.0e-5_ReKi*(1 - 1.0e-9_ReKi), 1.0e-5_ReKi*(1 + 1.0e-9_ReKi), kind == 6)
         do i = 1, min(n, 40)
            P(:,i) = [3.0_ReKi, 4.0_ReKi, 5.0_ReKi] + ext*[real(mod(i,2), ReKi), 0.5_ReKi*rnd(), 0.0_ReKi]
         end do
      case (8)   ! zero vorticity: 30 % of the particles, and every particle in the lower half box
         do i = 1, n
            if (rnd() < 0.3_ReKi .or. P(3,i) < 0.5_ReKi*L) A(:,i) = 0.0_ReKi
         end do
      case (9)   ! tiny vorticity around EqualRealNos' tolerance (1.1e-14 on |sum of norms|)
         do i = 1, n
            A(:,i) = A(:,i) * 1.0e-16_ReKi * (1 + 3*rnd())
         end do
      case (10)  ! eps larger than the cells: radius <= maxRegParam stops the subdivision early
         R(1:n) = 2.0_ReKi + 3.0_ReKi*R(1:n)
      case (11)
         R(1:n) = 0.0_ReKi
      end select
   end subroutine cloud

   subroutine one_case(kind, n, reg, bf, idd)
      integer,    intent(in) :: kind, n, reg, idd
      real(ReKi), intent(in) :: bf
      real(ReKi), allocatable :: P(:,:), A(:,:), R(:), CP(:,:), CP4(:,:), U0(:,:), UF(:,:), US(:,:), UK(:,:)
      real(ReKi), allocatable :: D(:,:), F0(:,:), K0(:,:), S(:,:), Ut(:,:)
      real(ReKi) :: L, dd, err, EF, EK, dmax, c3, smax, bound
      integer :: m, i, j, nb1, nb3, ErrStat
      character(ErrMsgLen) :: ErrMsg
      character(len=160) :: lab
      type(T_Tree) :: Tree
      logical :: bad
      allocate(P(3, max(n,1)), A(3, max(n,1)), R(max(n,1)))
      call cloud(kind, n, P, A, R, L)
      select case (idd)
      case (0); dd = 0.0_ReKi
      case (1); dd = 0.3_ReKi*L
      case default; dd = 2.2_ReKi*L
      end select
      m = min(max(n/4, 7), 400)
      allocate(CP(3, m+2), CP4(4, m+2), U0(3, m+2))
      do i = 1, m+2
         CP(:,i) = L*[1.4_ReKi*rnd()-0.2_ReKi, 1.4_ReKi*rnd()-0.2_ReKi, 1.4_ReKi*rnd()-0.2_ReKi]
         U0(:,i) = [rnd(), -rnd(), 0.0_ReKi]
      end do
      U0(3,1) = -0.0_ReKi
      if (n > 0) then
         CP(:,2) = P(:,1)                ! on a particle
         CP(:,3) = P(:,n)
      end if
      CP(:,4) = 30*L*[1.0_ReKi, 1.0_ReKi, 1.0_ReKi]   ! far from everything
      U0(:, m+1:) = qnan()               ! tail: untouched
      CP4(1:3,:) = CP; CP4(4,:) = qnan()
      write(lab,'(a,i0,a,i0,a,i0,a,f4.2,a,es8.1)') 'kind=', kind, ' n=', n, ' reg=', reg, ' bf=', bf, ' dd=', dd
      allocate(UF, source=U0); allocate(US, source=U0); allocate(UK, source=U0)
      ! --- Fortran tree (oracle) from U0
      call grow_tree_part(Tree, n, P(:,1:n), A(:,1:n), reg, R(1:n), 0)
      call ui_tree_part(Tree, m, CP, bf, dd, UF, ErrStat, ErrMsg)
      call cut_tree(Tree)
      ! --- port, strict, from U0 (CP with leading dimension 4 on odd cases)
      call ofk_set_mode(0_c_int)
      call tc_grow(key, n, P(:,1:n), A(:,1:n), R(1:n), reg, max(n,1))
      if (mod(n, 2) == 1) then
         call tc_ui(key, m, CP4, US, bf, dd, MINNORM)
      else
         call tc_ui(key, m, CP, US, bf, dd, MINNORM)
      end if
      ! --- port, fast, from U0
      call ofk_set_mode(1_c_int)
      call tc_grow(key, n, P(:,1:n), A(:,1:n), R(1:n), reg, max(n,1))
      call tc_ui(key, m, CP, UK, bf, dd, MINNORM)
      if (dumping) write(dumpunit) US, UK
      ! --- scale: sum over particles of |contribution|, and the direct sum D (from zero)
      allocate(S(3, m), D(3, m), Ut(3, m), F0(3, m), K0(3, m))
      S = 0; D = 0
      call ui_part_nograd(m, CP, n, P(:,1:n), A(:,1:n), reg, R(1:n), D)
      do j = 1, n
         Ut = 0
         call ui_part_nograd(m, CP, 1, P(:,j:j), A(:,j:j), reg, R(j:j), Ut)
         S = S + abs(Ut)
      end do
      ! --- claim 1/2: strict vs Fortran tree
      nrun = nrun + 1
      bad = .false.; err = 0; nb1 = 0
      do i = 1, m+2
         smax = maxval(S(:, min(i, m)))
         do j = 1, 3
            if (.not. same_bits(UF(j,i), US(j,i))) nb1 = nb1 + 1
            if (i > m) then
               if (.not. same_bits(UF(j,i), US(j,i))) bad = .true.
            else if (ieee_is_nan(UF(j,i)) .neqv. ieee_is_nan(US(j,i))) then
               bad = .true.
            else if (.not. ieee_is_nan(UF(j,i))) then
               err = max(err, abs(UF(j,i) - US(j,i)) / max(smax, tiny(1.0_ReKi)))
            end if
         end do
      end do
      if (nb1 > 0) nbit_c1 = nbit_c1 + 1
      worst_c1 = max(worst_c1, err)
      if (bad .or. err > rtol .or. (rtol == 0 .and. nb1 > 0)) then
         nfail = nfail + 1
         write(*,'(a,a,a,es10.3,a,i0,a,l1)') 'FAIL [claim1 strict vs Fortran tree | ', trim(lab), '] err=', err, &
              ' bitdiffs=', nb1, ' tail/NaN mismatch=', bad
      end if
      ! --- claim 3 (info): fast vs strict
      c3 = 0; nb3 = 0
      do i = 1, m
         smax = maxval(S(:, i))
         do j = 1, 3
            if (.not. same_bits(US(j,i), UK(j,i))) nb3 = nb3 + 1
            c3 = max(c3, abs(US(j,i) - UK(j,i)) / max(smax, tiny(1.0_ReKi)))
         end do
      end do
      do i = m+1, m+2
         do j = 1, 3
            if (.not. same_bits(U0(j,i), UK(j,i))) then
               nfail = nfail + 1; write(*,'(a,a)') 'FAIL [fast wrote a tail column | ', trim(lab)//']'
            end if
         end do
      end do
      if (nb3 > 0) nbit_c3 = nbit_c3 + 1
      if (c3 > huge(1.0_ReKi) .or. ieee_is_nan(c3)) c3 = 0   ! singular cases (see R2 below)
      if (c3 > 1.0e-12_ReKi) nflip = nflip + 1
      worst_c3 = max(worst_c3, c3)
      ! --- R2: from zero, against the direct sum
      nrun = nrun + 1
      F0 = 0; K0 = 0
      call grow_tree_part(Tree, n, P(:,1:n), A(:,1:n), reg, R(1:n), 0)
      call ui_tree_part(Tree, m, CP, bf, dd, F0, ErrStat, ErrMsg)
      call cut_tree(Tree)
      call tc_ui(key, m, CP, K0, bf, dd, MINNORM)   ! the fast tree is still in the handle (re-evaluation)
      dmax = max(maxval(abs(D)), tiny(1.0_ReKi))
      EF = maxval(abs(F0 - D)) / dmax
      EK = maxval(abs(K0 - D)) / dmax
      if (n == 0) then
         EF = 0; EK = maxval(abs(K0))
      end if
      if (EF > 1.0e-2_ReKi .and. info) write(*,'(a,a,a,es10.3,a,es10.3)') 'INFO large Fortran tree error | ', trim(lab), &
           ' E(Fortran tree)=', EF, ' E(fast)=', EK
      bound = max(1.1_ReKi*EF, r2floor)
      if (.not. (EF <= 1.0_ReKi)) then
         ! The Fortran tree is singular here (error larger than the field: a multipole evaluated at
         ! r ~ 1e-15 from a zero-radius node, RegFunctionPart=0 and DistanceDirect=0, audit T4).
         ! R2 is meaningless; require the port to be singular too, and count the case.
         nsing = nsing + 1
         if (EK <= 1.0_ReKi) then
            nfail = nfail + 1
            write(*,'(a,a,a,es10.3,a,es10.3)') 'FAIL [R2 singular case not reproduced | ', trim(lab), '] E(fast)=', EK, &
                 ' E(Fortran tree)=', EF
         end if
         return
      end if
      worst_EF = max(worst_EF, EF); worst_EK = max(worst_EK, EK)
      worst_r2ratio = max(worst_r2ratio, EK / bound)
      if (.not. (EK <= bound)) then
         nfail = nfail + 1
         write(*,'(a,a,a,es10.3,a,es10.3)') 'FAIL [R2 | ', trim(lab), '] E(fast)=', EK, ' E(Fortran tree)=', EF
      end if
   end subroutine one_case

   !> P6 at L0: after a warm-up, a fixed-size grow/eval sequence (both modes) allocates nothing.
   subroutine p6_sequence()
      integer, parameter :: n = 3000, m = 500
      real(ReKi) :: P(3,n), A(3,n), R(n), CP(3,m), U(3,m), L
      integer(c_long_long) :: c0(9), c1(9)
      integer :: it
      call cloud(2, n, P, A, R, L)
      do it = 1, m
         CP(:,it) = L*[rnd(), rnd(), rnd()]
      end do
      do it = 1, 24
         if (it == 4) call ofk_iv_counters(key, c0)
         call ofk_set_mode(int(mod(it, 2), c_int))
         A = A * (1.0_ReKi + 1.0e-3_ReKi*mod(it, 3))
         U = 0
         call tc_grow(key, n, P, A, R, mod(it, 3), n)
         call tc_ui(key, m, CP, U, 1.5_ReKi, 0.0_ReKi, MINNORM)
      end do
      call ofk_iv_counters(key, c1)
      nrun = nrun + 1
      write(*,'(a,i0,a,i0,a,i0)') 'P6 (L0): device allocations during 20 steady grow+eval: ', c1(1)-c0(1), &
           ' (instance total ', c1(1), '); evaluations ', c1(4)-c0(4)
      if (c1(1) /= c0(1)) then
         nfail = nfail + 1; write(*,'(a)') 'FAIL [P6] device allocation after warm-up'
      end if
      call ofk_set_mode(-1_c_int)
   end subroutine p6_sequence
end program l0_tc
