!> Helpers for l0_iv.f90 (unit induced-velocities): environment, RNG, case configuration,
!! comparisons and bookkeeping. Same environment and roles as tests/biotsavart-direct/l0_util.f90:
!!   OFK_L0_SEED (1), OFK_L0_CASES (60), OFK_L0_RTOL (0 = bitwise, incl. -0 and NaN pattern),
!!   OFK_L0_REFREF=1    A = oracle, B = oracle         (must be 0)
!!   OFK_L0_PORTPORT=1  A = port,   B = port           (R7 in process)
!!   OFK_L0_STRICTFAST=1 A = port strict, B = port fast (chain claim 3)
!!   default            A = oracle (pin, OPENMP=OFF flags), B = port in OFK_MODE (claim 1 / 2)
!!   OFK_L0_INFO=1      informational run: print, always exit 0
!!   OFK_L0_DUMP=<file> every B array compared, appended as unformatted stream (R7 across processes)
!! Error metric (architecture.md §4): per target column j, max_k |A-B| / max_k S(k,j), where
!! S = sum over sources of |contribution| from the oracle called once per source; entries with
!! S < 0 (tree methods: the same Fortran on both sides) and the packing outputs must be bitwise.
!! Columns outside the active range (sentinel tails, skipped strided columns) must be bitwise.
module l0_iv_util
   use NWTC_Library, only: ReKi, IntKi
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
   use, intrinsic :: iso_c_binding, only: c_int, c_long_long
   use FVW_Types, only: FVW_MiscVarType, T_Sgmt
   implicit none
   interface
      subroutine ofk_set_mode(m) bind(C, name='ofk_set_mode')
         import :: c_int
         integer(c_int), value :: m
      end subroutine ofk_set_mode
      subroutine ofk_iv_counters_c(key, out) bind(C, name='ofk_iv_counters')
         import :: c_long_long
         type(*), intent(in) :: key
         integer(c_long_long), intent(out) :: out(9)
      end subroutine ofk_iv_counters_c
      subroutine ofk_iv_free(key) bind(C, name='ofk_iv_free')
         type(*), intent(in) :: key
      end subroutine ofk_iv_free
   end interface
   type cfg_t
      integer :: nWings = 2, nSpan(3) = [3, 5, 2], nNWMax = 5, nNW = 4, nFWMax = 3, nFW = 2
      logical :: FWShed = .false., mirror = .false.
      integer :: reg = 3, regPart = 1, ppsW = 1, ppsLL = 1, velW = 1, velLL = 1
      integer :: nNWFree = 5, nFWFree = 3, iDS = 2
   end type cfg_t
   integer, parameter :: ROLE_REF = 0, ROLE_PORT = 1, ROLE_STRICT = 2, ROLE_FAST = 3
   integer, parameter :: CAT_PACK = 5, CAT_A7 = 6, NCAT = 6
   character(len=10), parameter :: catname(NCAT) = [character(len=10) :: 'Basic', 'TreePart', 'Part', &
                                                    'TreeSeg', 'pack', 'A7']
   integer    :: roleA = ROLE_REF, roleB = ROLE_PORT
   logical    :: refref = .false., portport = .false., strictfast = .false., dumping = .false., info = .false.
   integer    :: dumpunit = 0, seed = 1, ncases = 60, nrun = 0, nfail = 0
   real(ReKi) :: rtol = 0.0_ReKi, worst_all = 0.0_ReKi
   real(ReKi) :: worst(NCAT) = 0.0_ReKi
   integer    :: ncmp(NCAT) = 0, nbitcases(NCAT) = 0
contains
   subroutine l0_init()
      character(len=1024) :: buf
      call get_environment_variable('OFK_L0_SEED',  buf); if (len_trim(buf) > 0) read(buf,*) seed
      call get_environment_variable('OFK_L0_CASES', buf); if (len_trim(buf) > 0) read(buf,*) ncases
      call get_environment_variable('OFK_L0_RTOL',  buf); if (len_trim(buf) > 0) read(buf,*) rtol
      call get_environment_variable('OFK_L0_REFREF', buf); refref = (buf(1:1) == '1')
      call get_environment_variable('OFK_L0_PORTPORT', buf); portport = (buf(1:1) == '1')
      call get_environment_variable('OFK_L0_STRICTFAST', buf); strictfast = (buf(1:1) == '1')
      call get_environment_variable('OFK_L0_INFO', buf); info = (buf(1:1) == '1')   ! report only, exit 0
      if (refref) then
         roleA = ROLE_REF; roleB = ROLE_REF
      else if (portport) then
         roleA = ROLE_PORT; roleB = ROLE_PORT
      else if (strictfast) then
         roleA = ROLE_STRICT; roleB = ROLE_FAST
      end if
      call get_environment_variable('OFK_L0_DUMP', buf)
      if (len_trim(buf) > 0) then
         open(newunit=dumpunit, file=trim(buf), access='stream', form='unformatted', status='replace', action='write')
         dumping = .true.
      end if
      call seed_rng(seed)
   end subroutine l0_init

   logical function use_ref(role)
      integer, intent(in) :: role
      use_ref = (role == ROLE_REF)
      if (role == ROLE_STRICT) call ofk_set_mode(0_c_int)
      if (role == ROLE_FAST)   call ofk_set_mode(1_c_int)
   end function use_ref

   subroutine counters(Sgmt, out)
      type(T_Sgmt), intent(in) :: Sgmt
      integer(8), intent(out) :: out(9)
      integer(c_long_long) :: o(9)
      call ofk_iv_counters_c(Sgmt, o)
      out = int(o, 8)
   end subroutine counters

   real(ReKi) function qnan()
      qnan = ieee_value(1.0_ReKi, ieee_quiet_nan)
   end function qnan

   real(ReKi) function rnd()
      call random_number(rnd)
   end function rnd

   integer function rand_int(lo, hi)
      integer, intent(in) :: lo, hi
      rand_int = min(lo + int(rnd() * (hi - lo + 1)), hi)
   end function rand_int

   subroutine seed_rng(s)
      integer, intent(in) :: s
      integer :: n, k
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n)); put = s + 37*[(k - 1, k = 1, n)]
      call random_seed(put=put)
   end subroutine seed_rng

   logical function same_bits(a, b)
      real(ReKi), intent(in) :: a, b
      if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
         same_bits = ieee_is_nan(a) .and. ieee_is_nan(b)
      else
         same_bits = (a == b) .and. (sign(1.0_ReKi, a) == sign(1.0_ReKi, b))
      end if
   end function same_bits

   !> Compare A (reference role) and B. Columns 1..nAct (every `stride`-th) are "targets": rows
   !! 1..3 use the metric with scale S (S<0: bitwise); all other entries must be bitwise.
   !! Without S: max-normalized |A-B|/max|A| over the targets. exact=.true.: all bitwise.
   subroutine cmp2(label, cat, A, B, S, nAct, stride, exact)
      character(*), intent(in) :: label
      integer,      intent(in) :: cat
      real(ReKi),   intent(in) :: A(:,:), B(:,:)
      real(ReKi),   intent(in), optional :: S(:,:)
      integer,      intent(in), optional :: nAct, stride
      logical,      intent(in), optional :: exact
      integer :: i, j, na, st, nbit
      real(ReKi) :: err, dmax, smax, amax, dall
      logical :: tgt, bad, ex
      na = size(A,2); if (present(nAct)) na = nAct
      st = 1; if (present(stride)) st = stride
      ex = .false.; if (present(exact)) ex = exact
      err = 0; nbit = 0; bad = .false.; amax = 0; dall = 0
      if (any(shape(A) /= shape(B))) then
         nfail = nfail + 1; write(*,'(a,a)') 'FAIL [shape] ', label; return
      end if
      do j = 1, size(A,2)
         tgt = (j <= na) .and. (mod(j-1, st) == 0)
         dmax = 0; smax = 0
         do i = 1, size(A,1)
            if (.not. same_bits(A(i,j), B(i,j))) nbit = nbit + 1
            if (tgt .and. i <= 3 .and. .not. ex) then
               if (ieee_is_nan(A(i,j)) .neqv. ieee_is_nan(B(i,j))) then
                  bad = .true.
               else if (.not. ieee_is_nan(A(i,j))) then
                  dmax = max(dmax, abs(A(i,j) - B(i,j)))
                  amax = max(amax, abs(A(i,j)))
                  if (present(S)) then
                     if (S(min(i,size(S,1)), j) < 0) then
                        if (.not. same_bits(A(i,j), B(i,j))) bad = .true.
                     else
                        smax = max(smax, S(min(i,size(S,1)), j))
                     end if
                  end if
               end if
            else if (.not. same_bits(A(i,j), B(i,j))) then
               bad = .true.     ! non-target entries, packing outputs, exact comparisons
            end if
         end do
         if (present(S)) then
            if (dmax > 0) err = max(err, dmax / max(smax, tiny(1.0_ReKi)))
         else
            dall = max(dall, dmax)
         end if
      end do
      if (.not. present(S)) err = dall / max(amax, tiny(1.0_ReKi))
      if (dumping) write(dumpunit) B
      nrun = nrun + 1; ncmp(cat) = ncmp(cat) + 1
      if (nbit > 0) nbitcases(cat) = nbitcases(cat) + 1
      worst(cat) = max(worst(cat), err); worst_all = max(worst_all, err)
      if (bad .or. err > rtol .or. (rtol == 0 .and. nbit > 0)) then
         nfail = nfail + 1
         write(*,'(a,a,a,es10.3,a,i0,a,l1)') 'FAIL [', label, '] err=', err, ' bitdiffs=', nbit, ' exact-mismatch=', bad
      end if
   end subroutine cmp2

   subroutine cmp3(label, cat, A, B, S)
      character(*), intent(in) :: label
      integer,      intent(in) :: cat
      real(ReKi),   intent(in) :: A(:,:,:), B(:,:,:), S(:,:,:)
      call cmp2(label, cat, reshape(A, [size(A,1), size(A,2)*size(A,3)]), reshape(B, [size(B,1), size(B,2)*size(B,3)]), &
                reshape(S, [size(S,1), size(S,2)*size(S,3)]))
   end subroutine cmp3

   !> The Fortran packing outputs (Sgmt, Part) are API-visible misc: bitwise in every mode.
   subroutine cmp_sgmt(label, mA, mB)
      character(*), intent(in) :: label
      type(FVW_MiscVarType), intent(in) :: mA, mB
      call cmp2(label//' Points', CAT_PACK, mA%Sgmt%Points, mB%Sgmt%Points, exact=.true.)
      call cmp2(label//' Connct', CAT_PACK, real(mA%Sgmt%Connct, ReKi), real(mB%Sgmt%Connct, ReKi), exact=.true.)
      call cmp2(label//' Gamma', CAT_PACK, reshape(mA%Sgmt%Gamma, [1, size(mA%Sgmt%Gamma)]), &
                reshape(mB%Sgmt%Gamma, [1, size(mB%Sgmt%Gamma)]), exact=.true.)
      call cmp2(label//' Epsilon', CAT_PACK, reshape(mA%Sgmt%Epsilon, [1, size(mA%Sgmt%Epsilon)]), &
                reshape(mB%Sgmt%Epsilon, [1, size(mB%Sgmt%Epsilon)]), exact=.true.)
      if (mA%Sgmt%nAct /= mB%Sgmt%nAct .or. mA%Sgmt%nActP /= mB%Sgmt%nActP) then
         nfail = nfail + 1; write(*,'(a,a)') 'FAIL [nAct] ', label
      end if
      if (allocated(mA%Part%P)) then
         call cmp2(label//' Part%P', CAT_PACK, mA%Part%P, mB%Part%P, exact=.true.)
         call cmp2(label//' Part%Alpha', CAT_PACK, mA%Part%Alpha, mB%Part%Alpha, exact=.true.)
         call cmp2(label//' Part%RegParam', CAT_PACK, reshape(mA%Part%RegParam, [1, size(mA%Part%RegParam)]), &
                   reshape(mB%Part%RegParam, [1, size(mB%Part%RegParam)]), exact=.true.)
         if (mA%Part%nAct /= mB%Part%nAct) then
            nfail = nfail + 1; write(*,'(a,a)') 'FAIL [Part nAct] ', label
         end if
      end if
   end subroutine cmp_sgmt

   subroutine summary(prog)
      character(*), intent(in) :: prog
      integer :: k
      write(*,'(a,a,i0,a,i0,a,l1,a,l1,a,l1)') prog, ': seed=', seed, ' random cases=', ncases, ' refref=', refref, &
           ' portport=', portport, ' strictfast=', strictfast
      do k = 1, NCAT
         write(*,'(a,a10,a,i6,a,es10.3,a,i6)') '  ', catname(k), ' comparisons=', ncmp(k), ' worst err=', worst(k), &
              ' with bit diffs=', nbitcases(k)
      end do
      write(*,'(a,a,i0,a,i0,a,es9.2,a,es9.2,a)') prog, ': ', nrun - nfail, ' passed, ', nfail, ' failed; worst err ', &
           worst_all, ' (rtol ', rtol, ')'
      if (dumping) close(dumpunit)
   end subroutine summary
end module l0_iv_util
