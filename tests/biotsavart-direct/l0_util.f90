!> Shared helpers for the biotsavart-direct L0 tests: settings from the environment, RNG,
!! comparison of port vs reference (relative error, bitwise, NaN pattern, untouched columns),
!! and per-branch bookkeeping.
!!
!! Environment: OFK_L0_SEED (default 1), OFK_L0_CASES (random cases, default 200),
!!              OFK_L0_RTOL (default 0 = bitwise, which also requires equal signs of zero;
!!                           otherwise the bound on the error metric below),
!!              OFK_L0_METRIC = sum (default) | max: the gated error is
!!                 sum: max over targets of max_k|A_k - B_k| / max_k sum_j|c_jk|, the error
!!                      relative to the sum of the absolute contributions per target
!!                      (architecture.md §4; cancellation does not cause false failures);
!!                 max: max|A - B| / max|A| over the call (v0's metric). Both are printed.
!!              Which two runs are compared (A = "reference", B = "port"):
!!                 default            A = Fortran oracle,     B = port (mode from OFK_MODE)
!!                 OFK_L0_REFREF=1    A = Fortran oracle,     B = Fortran oracle (must be 0)
!!                 OFK_L0_PORTPORT=1  A = port,               B = port (determinism, R7)
!!                 OFK_L0_STRICTFAST=1 A = port, OFK_MODE=strict; B = port, OFK_MODE=fast
!!                                      (chain claim 3)
!!              OFK_L0_DUMP=<file>: append every B output array (unformatted stream), so runs
!!                 with different thread counts or processes can be compared with cmp (R7).
!!              FMA envelope (chain claim 3, user decision 2026-11-14: fast mode keeps FMA):
!!                 OFK_L0_ENVELOPE_OUT=<file>: after the run, write the worst gated error per
!!                    branch (with seed and case count) to <file>. Used by the *_ffma run, whose
!!                    error is the pin Fortran's own FMA-vs-no-FMA spread on these inputs.
!!                 OFK_L0_ENVELOPE=<file>: read such a file (same seed and case count, else
!!                    abort) and gate branch r at max(OFK_L0_RTOL, OFK_L0_ENVELOPE_FACTOR x
!!                    spread(r)); the factor defaults to 1.1.
module l0_util
   use NWTC_Library, only: ReKi, IntKi
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none
   interface
      subroutine ofk_set_mode(m) bind(C, name='ofk_set_mode')   ! ofk_runtime.hpp (test hook)
         import :: c_int
         integer(c_int), value :: m
      end subroutine ofk_set_mode
   end interface
   integer, parameter :: ROLE_REF = 0, ROLE_PORT = 1, ROLE_STRICT = 2, ROLE_FAST = 3
   integer    :: roleA = ROLE_REF, roleB = ROLE_PORT
   logical    :: strictfast = .false.
   logical    :: metric_sum = .true.
   integer    :: dumpunit = 0
   logical    :: dumping = .false.         ! newunit= returns NEGATIVE unit numbers
   real(ReKi) :: worst_max(0:4) = 0.0_ReKi            ! worst max-normalized error per branch
   integer, parameter :: NREG = 5                    ! branches 0..NREG-1
   real(ReKi), parameter :: SENTINEL = -12345.0_ReKi
   real(ReKi) :: rtol = 0.0_ReKi
   real(ReKi) :: rtol_reg(0:NREG-1) = 0.0_ReKi       ! per-branch bound: rtol, or the FMA envelope
   logical    :: have_envelope = .false.
   character(len=1024) :: envelope_out = ''
   integer    :: seed = 1, ncases = 200
   logical    :: refref = .false.
   logical    :: portport = .false.
   integer    :: nrun = 0, nfail = 0
   real(ReKi) :: worst_all = 0.0_ReKi
   real(ReKi) :: worst(0:NREG-1) = 0.0_ReKi          ! worst relative error per branch
   integer    :: ncase_reg(0:NREG-1) = 0              ! cases per branch
   integer    :: nbitdiff_reg(0:NREG-1) = 0           ! cases with at least one bit difference
   integer    :: nanCases(0:NREG-1) = 0               ! cases where the reference had NaN outputs
contains

   subroutine l0_init()
      character(len=1024) :: buf   ! was 64: truncated OFK_L0_DUMP paths (unit 2 finding)
      call get_environment_variable('OFK_L0_SEED',  buf); if (len_trim(buf) > 0) read(buf,*) seed
      call get_environment_variable('OFK_L0_CASES', buf); if (len_trim(buf) > 0) read(buf,*) ncases
      call get_environment_variable('OFK_L0_RTOL',  buf); if (len_trim(buf) > 0) read(buf,*) rtol
      call get_environment_variable('OFK_L0_REFREF', buf); refref = (len_trim(buf) > 0 .and. buf(1:1) == '1')
      call get_environment_variable('OFK_L0_PORTPORT', buf); portport = (len_trim(buf) > 0 .and. buf(1:1) == '1')
      call get_environment_variable('OFK_L0_STRICTFAST', buf); strictfast = (len_trim(buf) > 0 .and. buf(1:1) == '1')
      call get_environment_variable('OFK_L0_METRIC', buf); metric_sum = .not. (buf(1:3) == 'max')
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
      rtol_reg = rtol
      call get_environment_variable('OFK_L0_ENVELOPE_OUT', envelope_out)
      call get_environment_variable('OFK_L0_ENVELOPE', buf)
      if (len_trim(buf) > 0) call read_envelope(trim(buf))
      call seed_rng(seed)
   end subroutine l0_init

   !> Per-branch bound from the Fortran's own FMA spread (written by the *_ffma run).
   subroutine read_envelope(fname)
      character(*), intent(in) :: fname
      character(len=64) :: buf
      real(ReKi) :: factor, spread
      integer :: u, ios, eseed, ecases, nreg_file, r, rr
      factor = 1.1_ReKi
      call get_environment_variable('OFK_L0_ENVELOPE_FACTOR', buf); if (len_trim(buf) > 0) read(buf,*) factor
      open(newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(a,a)') 'l0: cannot read the FMA envelope ', fname; error stop 2
      end if
      read(u,*) eseed, ecases, nreg_file
      if (eseed /= seed .or. ecases /= ncases) then
         write(*,'(a,i0,a,i0,a,i0,a,i0)') 'l0: FMA envelope is for seed ', eseed, ' cases ', ecases, &
              '; this run has seed ', seed, ' cases ', ncases
         error stop 2
      end if
      do rr = 1, nreg_file
         read(u,*) r, spread
         rtol_reg(r) = max(rtol, factor * spread)
      end do
      close(u)
      have_envelope = .true.
      write(*,'(a,es9.2,a,f5.2,a,a)') 'l0: per-branch bound = max(OFK_L0_RTOL=', rtol, ', ', factor, &
           ' x Fortran FMA spread) from ', trim(fname)
   end subroutine read_envelope

   real(ReKi) function qnan()
      qnan = ieee_value(1.0_ReKi, ieee_quiet_nan)
   end function qnan

   !> Before a call in role `role`: select the port's mode for the strict/fast roles.
   !! Returns .true. if the Fortran oracle must be called, .false. for the port.
   logical function use_ref(role)
      integer, intent(in) :: role
      use_ref = (role == ROLE_REF)
      if (role == ROLE_STRICT) call ofk_set_mode(0_c_int)
      if (role == ROLE_FAST)   call ofk_set_mode(1_c_int)
   end function use_ref

   !> Compare one call. `active(j)` marks columns the routine may change; every other column,
   !! and rows 4.. of active columns, must equal U0 bit for bit (sentinels, "n++" tails).
   !! `S` (same shape as U0) holds sum_j |contribution_j| per component of active columns.
   subroutine compare(label, reg, U0, Uref, Uport, active, S)
      character(*), intent(in) :: label
      integer,      intent(in) :: reg
      real(ReKi),   intent(in) :: U0(:,:), Uref(:,:), Uport(:,:), S(:,:)
      logical,      intent(in) :: active(:)
      real(ReKi) :: scale, err, errsum, dmax, smax, gated
      integer :: i, j, nbit
      logical :: nanmismatch, touched, refnan
      scale = 0; err = 0; errsum = 0; nbit = 0; nanmismatch = .false.; touched = .false.; refnan = .false.
      do j = 1, size(U0,2)
         dmax = 0; smax = 0
         do i = 1, size(U0,1)
            if (active(j) .and. i <= 3) then
               if (ieee_is_nan(Uref(i,j))) refnan = .true.
               if (ieee_is_nan(Uref(i,j)) .neqv. ieee_is_nan(Uport(i,j))) then
                  nanmismatch = .true.
               else if (.not. ieee_is_nan(Uref(i,j))) then
                  scale = max(scale, abs(Uref(i,j)))
                  err   = max(err, abs(Uport(i,j) - Uref(i,j)))
                  dmax  = max(dmax, abs(Uport(i,j) - Uref(i,j)))
                  if (.not. ieee_is_nan(S(i,j))) smax = max(smax, S(i,j))
                  if (.not. same_bits(Uref(i,j), Uport(i,j))) nbit = nbit + 1
               end if
            else
               if (.not. same_bits(U0(i,j), Uport(i,j)) .or. .not. same_bits(U0(i,j), Uref(i,j))) touched = .true.
            end if
         end do
         if (dmax > 0) errsum = max(errsum, dmax / max(smax, tiny(1.0_ReKi)))
      end do
      err = err / max(scale, tiny(1.0_ReKi))
      if (dumping) write(dumpunit) Uport
      nrun = nrun + 1
      ncase_reg(reg) = ncase_reg(reg) + 1
      if (nbit > 0) nbitdiff_reg(reg) = nbitdiff_reg(reg) + 1
      if (refnan) nanCases(reg) = nanCases(reg) + 1
      gated = errsum; if (.not. metric_sum) gated = err
      worst(reg) = max(worst(reg), gated)
      worst_max(reg) = max(worst_max(reg), err)
      worst_all  = max(worst_all, gated)
      if (gated > rtol_reg(reg) .or. (rtol_reg(reg) == 0 .and. nbit > 0) .or. nanmismatch .or. touched) then
         nfail = nfail + 1
         write(*,'(a,a,a,i0,a,es10.3,a,es10.3,a,i0,2(a,l1))') 'FAIL [', label, '] reg=', reg, ' err/sum|c|=', errsum, &
              ' err/max|ref|=', err, ' bitdiffs=', nbit, ' nan-mismatch=', nanmismatch, ' wrote-outside-range=', touched
      end if
   end subroutine compare

   !> Equal value and equal sign (catches -0 vs +0); NaN handled by the caller.
   logical function same_bits(a, b)
      real(ReKi), intent(in) :: a, b
      if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
         same_bits = ieee_is_nan(a) .and. ieee_is_nan(b)
      else
         same_bits = (a == b) .and. (sign(1.0_ReKi, a) == sign(1.0_ReKi, b))
      end if
   end function same_bits

   subroutine summary(prog, names)
      character(*), intent(in) :: prog
      character(*), intent(in) :: names(0:)
      integer :: r
      character(len=8) :: mname
      integer :: u
      mname = 'sum|c|'; if (.not. metric_sum) mname = 'max|ref|'
      write(*,'(a,a,i0,a,i0,a,i0,a,l1,a,l1,a,l1,a,a)') prog, ': seed=', seed, ' cases=', nrun, ' (random ', ncases, &
           ') refref=', refref, ' portport=', portport, ' strictfast=', strictfast, ' metric=err/', trim(mname)
      do r = 0, size(names) - 1
         write(*,'(a,a12,a,i5,a,es10.3,a,es10.3,a,i5,a,i4)') '  reg ', names(r), ' cases=', ncase_reg(r), &
              ' worst err=', worst(r), ' (err/max|ref| ', worst_max(r), ') cases with bit diffs=', nbitdiff_reg(r), &
              ' NaN-output cases=', nanCases(r)
         if (have_envelope) write(*,'(a,a12,a,es10.3,a,f6.3)') '  reg ', names(r), ' bound=', rtol_reg(r), &
              ' worst/bound=', worst(r) / max(rtol_reg(r), tiny(1.0_ReKi))
      end do
      if (len_trim(envelope_out) > 0) then
         open(newunit=u, file=trim(envelope_out), status='replace', action='write')
         write(u,'(i0,1x,i0,1x,i0)') seed, ncases, size(names)
         do r = 0, size(names) - 1
            write(u,'(i0,1x,es24.16e3)') r, worst(r)
         end do
         close(u)
         write(*,'(a,a)') 'l0: wrote the per-branch FMA spread to ', trim(envelope_out)
      end if
      write(*,'(a,a,i0,a,i0,a,es9.2,a,es9.2,a)') prog, ': ', nrun - nfail, ' passed, ', nfail, &
           ' failed; worst err ', worst_all, ' (rtol ', rtol, ')'
      if (dumping) close(dumpunit)
   end subroutine summary

   subroutine seed_rng(s)
      integer, intent(in) :: s
      integer :: n, k
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n)); put = s + 37*[(k - 1, k = 1, n)]
      call random_seed(put=put)
   end subroutine seed_rng

   integer function rand_int(lo, hi)
      integer, intent(in) :: lo, hi
      real(ReKi) :: r
      call random_number(r)
      rand_int = min(lo + int(r * (hi - lo + 1)), hi)
   end function rand_int

end module l0_util
