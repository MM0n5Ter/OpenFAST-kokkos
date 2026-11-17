!> L0 differential test of unit circulation-solve (OLAF unit 4a, F10-B): the pin's
!! Wings_ComputeCirculationPolarData (modules renamed *_ref at build time; OPENMP=OFF flags) against the
!! port's (the overlay's FVW_Wings with the OPENFAST_KOKKOS shim), in one process, through the real Fortran
!! interface and registry types. A second pin copy built with FMA (*_rfm) gives the Fortran's own FMA
!! envelope for the FMA experiment.
!!
!! Each solve's iteration count is read without touching either side: AirfoilInfo::AFI_ComputeAirfoilCoefs
!! is wrapped (ld --wrap, afi_count.c), and CirculationFromPolarData calls it once per control point.
!!
!! Per case (unit.md audit list): RegFunction 0-4 (C2, C5), 1-3 wings with unequal nSpan (C8), iNWStart 1-3,
!! InductionAtCP T/F, FirstCall and GammaScale = 0 / 0.5 / 1 (C10), degenerate geometry (C3: CP on a ring
!! corner and on a bound segment, zero-length span segment, collapsed ring, Gamma = +-0 rings, NaN in r_NW),
!! Eps_NW rows 2-3 poisoned (C7), knife-edge ConvCrit (C4: the two adjacent doubles where the pin's
!! iteration count changes), A7 sequences (C12: re-evaluation, misc copies, changed Eps / r_NW / CP), P6
!! (C11: no allocation after the first solve). Compared bitwise (NaN = NaN): z%W%Gamma_LL and m%W%{alpha_LL,
!! Vreln_LL, Vind_CP, Vtot_CP, Vind_LL}, ErrStat, and the iteration count.
!!
!! Variants of the port (ofk_cs_set_variant), every one against the pin:
!!   strict (threshold / serial / parallel_for), fast F0 (serial / parallel_for): must be bitwise, same
!!   iteration counts (claim 1); fast F1 (FMA): informational, as is the pin-with-FMA envelope.
!! The rest of the process (unit 2's LL induction inside the solve) runs in strict mode (ofk_set_mode(0)),
!! so the comparison isolates this unit.
!!
!! Environment: OFK_L0_LLVEL (1): p%VelocityMethod(2) for the solve's LL induction (4 = the Fortran segment
!! tree on both sides: isolates this unit on CUDA, where unit 2's LL kernel is a few ulp from the Fortran),
!! OFK_L0_SEED (1), OFK_L0_CASES (200), OFK_L0_DUMP=<file> (every port output, stream, for R7
!! across processes and thread counts), OFK_L0_INFO=1 (always exit 0).
module l0_cs_afi
   use, intrinsic :: iso_c_binding, only: c_long_long
   implicit none
   integer(c_long_long), bind(C, name='ofk_afi_calls') :: afi_calls   ! defined in afi_count.c
end module l0_cs_afi

program l0_cs
   use l0_cs_afi
   use NWTC_Library
   use FVW_Types
   use AirfoilInfo,       only: AFI_Init
   use AirfoilInfo_Types, only: AFI_ParameterType, AFI_InitInputType
   use FVW_Subs_ref,      only: FVW_InitStates, FVW_InitMiscVarsPostParam
   use FVW_Wings_ref,     only: r_Solve => Wings_ComputeCirculationPolarData
   use FVW_Wings_rfm,     only: f_Solve => Wings_ComputeCirculationPolarData
   use FVW_Wings,         only: p_Solve => Wings_ComputeCirculationPolarData
   use, intrinsic :: iso_c_binding, only: c_int, c_long_long
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
   implicit none

   interface
      subroutine ofk_set_mode(m) bind(C, name='ofk_set_mode')
         import :: c_int
         integer(c_int), value :: m
      end subroutine ofk_set_mode
      subroutine ofk_cs_set_variant(mode, fma, par) bind(C, name='ofk_cs_set_variant')
         import :: c_int
         integer(c_int), value :: mode, fma, par
      end subroutine ofk_cs_set_variant
      subroutine ofk_cs_counters(key, out) bind(C, name='ofk_cs_counters')
         import :: c_long_long
         type(*), intent(in) :: key
         integer(c_long_long), intent(out) :: out(8)
      end subroutine ofk_cs_counters
      subroutine ofk_cs_free(key) bind(C, name='ofk_cs_free')
         type(*), intent(in) :: key
      end subroutine ofk_cs_free
   end interface

   ! sides
   integer, parameter :: S_PIN = 0, S_PINFMA = 1, S_PORT = 2
   ! port variants: name, mode, fma, par
   integer, parameter :: NV = 6
   character(len=16), parameter :: vname(NV) = [character(len=16) :: 'strict', 'strict-serial', 'strict-parfor', &
                                                'F0-serial', 'F0-parfor', 'F1']
   integer, parameter :: vmode(NV) = [0, 0, 0, 1, 1, 1], vfma(NV) = [-1, -1, -1, 0, 0, 1], vpar(NV) = [-1, 0, 1, 0, 1, -1]
   logical, parameter :: vexact(NV) = [.true., .true., .true., .true., .true., .false.]

   type cfg_t
      integer    :: nWings = 3, nSpan(3) = [8, 6, 7], iNWStart = 2, nNWMax = 6, nNW = 6, nFWMax = 2, nFW = 2
      integer    :: reg = 3, maxIter = 30, degen = 0
      real(ReKi) :: conv = 0.001_ReKi, relax = 0.1_ReKi, U = 10.0_ReKi, alpha = 0.10_ReKi, chord = 1.0_ReKi
      real(ReKi) :: epsScale = 1.0_ReKi, gscale = 1.0_ReKi, gpert = 0.2_ReKi
      logical    :: firstCall = .false., atCP = .true., onLL = .false.
   end type cfg_t

   type res_t
      integer :: errStat = 0, iters = -1
      type(FVW_ConstraintStateType) :: z
      type(FVW_MiscVarType)         :: m
   end type res_t

   type(AFI_ParameterType), allocatable :: AFInfo(:)
   integer    :: llvel = 1
   integer    :: seed = 1, ncases = 200, nfail = 0, ncase = 0, nsolve = 0, dumpunit = 0
   logical    :: info = .false., dumping = .false.
   integer(8) :: rngs = 88172645463325252_8
   ! per-variant bookkeeping
   integer    :: vcmp(NV) = 0, vbit(NV) = 0, vitd(NV) = 0, vitdmax(NV) = 0
   real(ReKi) :: vrel(NV) = 0.0_ReKi
   integer    :: ecmp = 0, eitd = 0, eitdmax = 0, eitdhist(-30:30) = 0, vitdhist(-30:30) = 0
   real(ReKi) :: erel = 0.0_ReKi
   integer    :: ithist(0:31) = 0
   integer    :: c3n = 0, c3fail = 0            ! claim-3-style check of F1, per case
   real(ReKi) :: c3worst = 0.0_ReKi              ! worst rel(F1) / max(1e-12, 1.1 rel(pin+FMA))
   integer    :: nrefref = 0                ! cases where the pin is not reproducible (core outputs gated)
   integer    :: kn = 0, knF1 = 0, knFma = 0   ! knife-edge solves; with F1 / pin+FMA count /= pin
   ! Recorded F1 claim-3 failures (tests/circulation-solve/unit.md, "Recorded F1 failures"): the random cases of
   ! the 400-case runs of seeds 6 and 8 where F1 exceeded the per-case bound. A random case is fully determined by
   ! the generator state before rand_cfg(); replayed every run (after the random cases, generator restored).
   integer, parameter :: NREC = 2
   integer(8), parameter :: rec_rng(NREC) = [-3469895770392615260_8, -629506827529011285_8]
   character(len=24), parameter :: rec_name(NREC) = [character(len=24) :: 'seed 6, 400 cases', 'seed 8, 400 cases']
   integer(8) :: case_rng = 0                 ! generator state before the current random case's rand_cfg()
   integer    :: irec = 0                     ! > 0 while replaying recorded case irec
   real(ReKi) :: rec_ratio(NREC) = -1.0_ReKi, rec_relF1(NREC) = 0.0_ReKi, rec_relF(NREC) = 0.0_ReKi
   integer    :: rec_dit(NREC) = 0
   integer    :: i

   call init()
   call fixed_cases()
   do i = 1, ncases
      case_rng = rngs
      call one_case(rand_cfg(), 'random')
   end do
   call recorded_cases()
   call knife_cases()
   call a7_sequence()
   call p6_sequence()
   call summary()
   if (dumping) close(dumpunit)
   if (nfail > 0 .and. .not. info) error stop 1

contains

   ! ===================================================================================== setup
   subroutine init()
      character(len=256) :: s
      integer :: st, k
      type(AFI_InitInputType) :: ai
      integer(IntKi) :: ErrStat
      character(ErrMsgLen) :: ErrMsg
      call NWTC_Init(EchoLibVer=.false.)
      call get_environment_variable('OFK_L0_SEED', s, status=st);  if (st == 0) read(s, *) seed
      call get_environment_variable('OFK_L0_CASES', s, status=st); if (st == 0) read(s, *) ncases
      call get_environment_variable('OFK_L0_LLVEL', s, status=st); if (st == 0) read(s, *) llvel
      call get_environment_variable('OFK_L0_INFO', s, status=st);  info = (st == 0 .and. trim(s) == '1')
      call get_environment_variable('OFK_L0_DUMP', s, status=st)
      if (st == 0 .and. len_trim(s) > 0) then
         dumping = .true.
         open(newunit=dumpunit, file=trim(s), access='stream', form='unformatted', status='replace')
      end if
      rngs = rngs + 7919_8 * int(seed, 8)
      do k = 1, 10; rngs = xorshift(rngs); end do
      ! Two BAR airfoils (r-test), as AeroDyn reads them for ad_BAR_OLAF
      allocate(AFInfo(2))
      ai%AFTabMod = 1; ai%InCol_Alfa = 1; ai%InCol_Cl = 2; ai%InCol_Cd = 3; ai%InCol_Cm = 4; ai%InCol_Cpmin = 0
      ai%UAMod = 0
      ai%FileName = OFK_AFDIR//'/BAR0_AeroDyn_Polar_15.dat'
      call AFI_Init(ai, AFInfo(1), ErrStat, ErrMsg)
      if (ErrStat >= AbortErrLev) then; print *, 'AFI_Init: ', trim(ErrMsg); error stop 3; end if
      ai%FileName = OFK_AFDIR//'/BAR0_AeroDyn_Polar_25.dat'
      call AFI_Init(ai, AFInfo(2), ErrStat, ErrMsg)
      if (ErrStat >= AbortErrLev) then; print *, 'AFI_Init: ', trim(ErrMsg); error stop 3; end if
      call ofk_set_mode(0_c_int)   ! the rest of the process (unit 2 LL induction) in strict mode
      print '(a,i0,a,i0,a,i0)', 'l0_cs: seed=', seed, ' random cases=', ncases, ' LL VelocityMethod=', llvel
   end subroutine init

   integer(8) function xorshift(x) result(y)
      integer(8), intent(in) :: x
      y = ieor(x, ishft(x, 13)); y = ieor(y, ishft(y, -7)); y = ieor(y, ishft(y, 17))
   end function xorshift

   real(ReKi) function rnd()   ! [0,1)
      rngs = xorshift(rngs)
      rnd = real(ishft(rngs, -11), ReKi) / 9007199254740992.0_ReKi
   end function rnd

   integer function irnd(lo, hi)
      integer, intent(in) :: lo, hi
      irnd = lo + min(hi - lo, int(rnd() * real(hi - lo + 1, ReKi)))
   end function irnd

   real(ReKi) function qnan()
      qnan = ieee_value(1.0_ReKi, ieee_quiet_nan)
   end function qnan

   type(cfg_t) function rand_cfg() result(c)
      integer :: iW
      c%nWings = irnd(1, 3)
      do iW = 1, 3; c%nSpan(iW) = irnd(1, 12); end do
      c%iNWStart = irnd(1, 3)
      c%nNWMax = c%iNWStart + irnd(1, 6)
      c%nNW = irnd(c%iNWStart - 1, c%nNWMax)
      c%nFWMax = irnd(0, 3); c%nFW = irnd(0, c%nFWMax)
      c%reg = irnd(0, 4)
      c%conv = 10.0_ReKi**(-1.5_ReKi - 2.0_ReKi*rnd())
      c%relax = 0.1_ReKi + 0.5_ReKi*rnd()
      c%U = 4.0_ReKi + 12.0_ReKi*rnd()
      c%alpha = -0.05_ReKi + 0.25_ReKi*rnd()
      c%chord = 0.5_ReKi + 1.5_ReKi*rnd()
      c%epsScale = 10.0_ReKi**(-1.5_ReKi + 2.0_ReKi*rnd())
      c%gpert = 0.5_ReKi*rnd()
      c%firstCall = rnd() < 0.15_ReKi
      c%atCP = rnd() < 0.7_ReKi
      c%onLL = rnd() < 0.5_ReKi
      select case (irnd(0, 5))
      case (0); c%gscale = 0.0_ReKi
      case (1); c%gscale = 0.5_ReKi
      case default; c%gscale = 1.0_ReKi
      end select
      c%degen = 0
      if (rnd() < 0.3_ReKi) c%degen = irnd(1, 7)
   end function rand_cfg

   !> A rotor-like lifting line: blade iW along e_r = (0, cos t, sin t), chord along x, wind along x with an
   !! angle of attack in the blade's normal; the near wake convects along x. Realistic enough that the pin's
   !! iteration converges in a few to 30 iterations.
   subroutine setup(c, p, x, m, zprev)
      type(cfg_t),                   intent(in)    :: c
      type(FVW_ParameterType),       intent(inout) :: p
      type(FVW_ContinuousStateType), intent(inout) :: x
      type(FVW_MiscVarType),         intent(inout) :: m
      type(FVW_ConstraintStateType), intent(inout) :: zprev
      integer :: iW, j, k, nS, ErrStat
      character(ErrMsgLen) :: ErrMsg
      real(ReKi) :: th, er(3), Tc(3), Nv(3), r0, R, ds, rr, xi, dx, G0
      p%nWings = c%nWings
      allocate(p%W(p%nWings))
      do iW = 1, p%nWings
         nS = c%nSpan(iW); p%W(iW)%nSpan = nS
         allocate(p%W(iW)%s_LL(nS+1), p%W(iW)%s_CP(nS), p%W(iW)%chord_CP(nS), p%W(iW)%AFindx(nS,1))
         do j = 1, nS+1; p%W(iW)%s_LL(j) = real(j-1, ReKi); end do
         do j = 1, nS
            p%W(iW)%s_CP(j) = real(j, ReKi) - 0.5_ReKi
            p%W(iW)%chord_CP(j) = c%chord*(1.2_ReKi - 0.6_ReKi*real(j, ReKi)/real(nS+1, ReKi))*(0.95_ReKi + 0.1_ReKi*rnd())
            p%W(iW)%AFindx(j,1) = 1 + mod(j + iW, 2)
         end do
      end do
      p%nNWMax = c%nNWMax; p%nFWMax = c%nFWMax; p%nNWFree = c%nNWMax; p%nFWFree = c%nFWMax
      p%FWShedVorticity = .false.; p%ShearModel = 0
      p%RegFunction = c%reg; p%RegFunctionPart = 0
      p%VelocityMethod(1) = 1; p%VelocityMethod(2) = llvel; p%PartPerSegment = 1; p%TreeBranchFactor = 1.5_ReKi
      p%iNWStart = c%iNWStart
      p%CircSolvMaxIter = c%maxIter; p%CircSolvConvCrit = c%conv; p%CircSolvRelaxation = c%relax
      p%InductionAtCP = c%atCP
      p%KinVisc = 1.464e-5_ReKi
      call FVW_InitStates(x, p, ErrStat, ErrMsg)
      allocate(m%W(p%nWings))
      allocate(zprev%W(p%nWings))
      do iW = 1, p%nWings
         nS = p%W(iW)%nSpan
         th = 2.0_ReKi*Pi*real(iW-1, ReKi)/real(p%nWings, ReKi) + 0.3_ReKi*rnd()
         er = [0.0_ReKi, cos(th), sin(th)]
         Tc = [1.0_ReKi, 0.0_ReKi, 0.0_ReKi]
         Nv = [0.0_ReKi, -sin(th), cos(th)]
         r0 = 0.5_ReKi + 0.5_ReKi*rnd(); R = r0 + 3.0_ReKi + 6.0_ReKi*rnd(); ds = (R - r0)/real(max(nS, 1), ReKi)
         dx = 0.3_ReKi*c%chord*(0.8_ReKi + 0.4_ReKi*rnd())
         allocate(m%W(iW)%CP(3,nS), m%W(iW)%Tang(3,nS), m%W(iW)%Norm(3,nS), m%W(iW)%Orth(3,nS), m%W(iW)%dl(3,nS))
         allocate(m%W(iW)%Area(nS), m%W(iW)%alpha_LL(nS), m%W(iW)%Vreln_LL(nS))
         allocate(m%W(iW)%Vind_CP(3,nS), m%W(iW)%Vtot_CP(3,nS), m%W(iW)%Vstr_CP(3,nS), m%W(iW)%Vwnd_CP(3,nS))
         allocate(m%W(iW)%r_LL(3,nS+1,2), m%W(iW)%Vind_LL(3,nS+1))
         allocate(m%W(iW)%Vind_NW(3,nS+1,p%nNWMax+1), m%W(iW)%Vind_FW(3,2,p%nFWMax+1))
         m%W(iW)%alpha_LL = -999999.0_ReKi; m%W(iW)%Vreln_LL = -999999.0_ReKi
         m%W(iW)%Vind_CP = -999999.0_ReKi; m%W(iW)%Vtot_CP = -999999.0_ReKi
         m%W(iW)%Vind_LL = -999999.0_ReKi; m%W(iW)%Vind_NW = -1.0_ReKi; m%W(iW)%Vind_FW = -1.0_ReKi
         do j = 1, nS+1
            rr = r0 + real(j-1, ReKi)*ds*(1.0_ReKi + 0.05_ReKi*(rnd() - 0.5_ReKi))
            m%W(iW)%r_LL(:,j,1) = rr*er
            m%W(iW)%r_LL(:,j,2) = rr*er + 0.75_ReKi*c%chord*Tc
            do k = 1, size(x%W(iW)%r_NW, 3)
               x%W(iW)%r_NW(:,j,k) = rr*er + (0.25_ReKi*c%chord + real(k-1, ReKi)*dx)*Tc &
                                    + 0.02_ReKi*c%chord*real(k-1, ReKi)*[rnd()-0.5_ReKi, rnd()-0.5_ReKi, rnd()-0.5_ReKi]
            end do
         end do
         if (c%onLL) x%W(iW)%r_NW(:,:,1) = m%W(iW)%r_LL(:,:,1)   ! CP on its own bound segment (C3)
         do j = 1, nS
            m%W(iW)%CP(:,j) = 0.5_ReKi*(m%W(iW)%r_LL(:,j,1) + m%W(iW)%r_LL(:,j+1,1))
            m%W(iW)%Tang(:,j) = Tc; m%W(iW)%Norm(:,j) = Nv
            m%W(iW)%Orth(:,j) = er
            m%W(iW)%dl(:,j) = m%W(iW)%r_LL(:,j+1,1) - m%W(iW)%r_LL(:,j,1)
            m%W(iW)%Area(j) = p%W(iW)%chord_CP(j)*sqrt(sum(m%W(iW)%dl(:,j)**2))
            m%W(iW)%Vwnd_CP(:,j) = c%U*(cos(c%alpha)*Tc + sin(c%alpha)*Nv) + 0.3_ReKi*[rnd(), rnd(), rnd()]
            m%W(iW)%Vstr_CP(:,j) = 0.2_ReKi*[rnd()-0.5_ReKi, rnd()-0.5_ReKi, rnd()-0.5_ReKi]
         end do
         ! a bound-circulation-like distribution: guesses and the wake beyond iNWStart (LL #b)
         G0 = Pi*c%chord*c%U*sin(c%alpha + 0.05_ReKi)
         allocate(zprev%W(iW)%Gamma_LL(nS))
         do j = 1, nS
            xi = (real(j, ReKi) - 0.5_ReKi)/real(nS, ReKi)
            zprev%W(iW)%Gamma_LL(j) = G0*sqrt(max(0.0_ReKi, 1.0_ReKi - (2.0_ReKi*xi - 1.0_ReKi)**2) + 0.05_ReKi) &
                                      *(1.0_ReKi + c%gpert*(rnd() - 0.5_ReKi))
            do k = 1, size(x%W(iW)%Gamma_NW, 2)
               x%W(iW)%Gamma_NW(j,k) = zprev%W(iW)%Gamma_LL(j)*(0.9_ReKi + 0.2_ReKi*rnd())
               x%W(iW)%Eps_NW(:,j,k) = c%epsScale*0.25_ReKi*p%W(iW)%chord_CP(j)*(0.9_ReKi + 0.2_ReKi*[rnd(), rnd(), rnd()])
            end do
         end do
         do k = 1, size(x%W(iW)%r_FW, 3)
            do j = 1, size(x%W(iW)%r_FW, 2)
               x%W(iW)%r_FW(:,j,k) = (r0 + real(j-1, ReKi)*(R - r0))*er + (0.25_ReKi*c%chord &
                                     + real(c%nNWMax + k, ReKi)*dx)*Tc
            end do
         end do
         x%W(iW)%Gamma_FW = 0.5_ReKi*G0; x%W(iW)%Eps_FW = c%epsScale*0.3_ReKi*c%chord
      end do
      call FVW_InitMiscVarsPostParam(p, m, ErrStat, ErrMsg)
      m%nNW = c%nNW; m%nFW = c%nFW
      m%FirstCall = c%firstCall
      call degenerate(c, p, x, m, zprev)
   end subroutine setup

   !> C3 / C7 edge cases, applied after the regular setup.
   subroutine degenerate(c, p, x, m, zprev)
      type(cfg_t),                   intent(in)    :: c
      type(FVW_ParameterType),       intent(in)    :: p
      type(FVW_ContinuousStateType), intent(inout) :: x
      type(FVW_MiscVarType),         intent(inout) :: m
      type(FVW_ConstraintStateType), intent(inout) :: zprev
      integer :: iW, s, d, nS
      iW = irnd(1, p%nWings); nS = p%W(iW)%nSpan
      s = irnd(1, nS); d = irnd(1, p%iNWStart)
      select case (c%degen)
      case (1)   ! a control point exactly on a ring corner (norm_a = 0: denominator <= PRECISION_UI)
         m%W(iW)%CP(:, s) = x%W(iW)%r_NW(:, s, d)
      case (2)   ! a zero-length spanwise segment (norm2_r0 = 0)
         x%W(iW)%r_NW(:, s+1, d) = x%W(iW)%r_NW(:, s, d)
      case (3)   ! a collapsed ring: all four corners coincide
         x%W(iW)%r_NW(:, s+1, d) = x%W(iW)%r_NW(:, s, d)
         x%W(iW)%r_NW(:, s, d+1) = x%W(iW)%r_NW(:, s, d)
         x%W(iW)%r_NW(:, s+1, d+1) = x%W(iW)%r_NW(:, s, d)
      case (4)   ! Gamma = +0 and -0 rings (guess path)
         zprev%W(iW)%Gamma_LL(s) = 0.0_ReKi
         zprev%W(iW)%Gamma_LL(nS + 1 - s) = -0.0_ReKi
      case (5)   ! NaN in a ring corner (the guards do not return on NaN: it propagates)
         x%W(iW)%r_NW(irnd(1,3), s, d+1) = qnan()
      case (6)   ! Eps_NW rows 2 and 3 are not read (C7)
         x%W(iW)%Eps_NW(2:3, :, :) = qnan()
      case (7)   ! a control point on a bound segment line (norm2_orth ~ 0), every wing
         do iW = 1, p%nWings
            do s = 1, p%W(iW)%nSpan
               m%W(iW)%CP(:, s) = 0.5_ReKi*(x%W(iW)%r_NW(:, s, 1) + x%W(iW)%r_NW(:, s+1, 1))
            end do
         end do
      end select
   end subroutine degenerate

   subroutine teardown(p, x, m, z)
      type(FVW_ParameterType),       intent(inout) :: p
      type(FVW_ContinuousStateType), intent(inout) :: x
      type(FVW_MiscVarType),         intent(inout) :: m
      type(FVW_ConstraintStateType), intent(inout) :: z
      integer :: ErrStat
      character(ErrMsgLen) :: ErrMsg
      call FVW_DestroyParam(p, ErrStat, ErrMsg)
      call FVW_DestroyContState(x, ErrStat, ErrMsg)
      call FVW_DestroyMisc(m, ErrStat, ErrMsg)
      call FVW_DestroyConstrState(z, ErrStat, ErrMsg)
   end subroutine teardown

   ! ===================================================================================== solves
   !> One solve by one side on a copy of m0 (or in place on mio when given): result in r.
   subroutine solve(side, c, p, x, m0, zprev, r, mio)
      integer,                       intent(in)    :: side
      type(cfg_t),                   intent(in)    :: c
      type(FVW_ParameterType),       intent(in)    :: p
      type(FVW_ContinuousStateType), intent(in)    :: x
      type(FVW_MiscVarType),         intent(in)    :: m0
      type(FVW_ConstraintStateType), intent(in)    :: zprev
      type(res_t),                   intent(inout) :: r
      type(FVW_MiscVarType), optional, intent(inout) :: mio
      integer :: iW, ErrStat, nCP
      integer(c_long_long) :: a0
      character(ErrMsgLen) :: ErrMsg
      call FVW_DestroyConstrState(r%z, ErrStat, ErrMsg)
      allocate(r%z%W(p%nWings))
      nCP = 0
      do iW = 1, p%nWings
         allocate(r%z%W(iW)%Gamma_LL(p%W(iW)%nSpan)); r%z%W(iW)%Gamma_LL = -777.0_ReKi
         nCP = nCP + p%W(iW)%nSpan
      end do
      if (present(mio)) then
         a0 = afi_calls
         call run_side(side, c, p, x, mio, zprev, r)
         call FVW_DestroyMisc(r%m, ErrStat, ErrMsg)
         call FVW_CopyMisc(mio, r%m, MESH_NEWCOPY, ErrStat, ErrMsg)
      else
         call FVW_DestroyMisc(r%m, ErrStat, ErrMsg)
         call FVW_CopyMisc(m0, r%m, MESH_NEWCOPY, ErrStat, ErrMsg)
         a0 = afi_calls
         call run_side(side, c, p, x, r%m, zprev, r)
      end if
      r%iters = -1
      if (nCP > 0) r%iters = int(afi_calls - a0) / nCP - merge(1, 0, c%firstCall)
      nsolve = nsolve + 1
   end subroutine solve

   subroutine run_side(side, c, p, x, m, zprev, r)
      integer,                       intent(in)    :: side
      type(cfg_t),                   intent(in)    :: c
      type(FVW_ParameterType),       intent(in)    :: p
      type(FVW_ContinuousStateType), intent(in)    :: x
      type(FVW_MiscVarType),         intent(inout) :: m
      type(FVW_ConstraintStateType), intent(in)    :: zprev
      type(res_t),                   intent(inout) :: r
      character(ErrMsgLen) :: ErrMsg
      select case (side)
      case (S_PIN);    call r_Solve(r%z, zprev, p, x, m, AFInfo, c%gscale, r%errStat, ErrMsg, 1)
      case (S_PINFMA); call f_Solve(r%z, zprev, p, x, m, AFInfo, c%gscale, r%errStat, ErrMsg, 1)
      case default;    call p_Solve(r%z, zprev, p, x, m, AFInfo, c%gscale, r%errStat, ErrMsg, 1)
      end select
   end subroutine run_side

   logical function same(a, b)   ! bitwise, any NaN equals any NaN (payloads are not compared)
      real(ReKi), intent(in) :: a, b
      if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
         same = ieee_is_nan(a) .and. ieee_is_nan(b)
      else
         same = transfer(a, 1_8) == transfer(b, 1_8)
      end if
   end function same

   !> Compare two results: number of differing values, worst |dGamma| / max |Gamma_ref|.
   subroutine compare(p, A, B, ndiff, rel, coreOnly)
      type(FVW_ParameterType), intent(in)  :: p
      type(res_t),             intent(in)  :: A, B
      integer,                 intent(out) :: ndiff
      real(ReKi),              intent(out) :: rel
      logical, optional,       intent(in)  :: coreOnly   ! Gamma_LL, alpha_LL, Vreln_LL, ErrStat only
      integer :: iW, j, k
      logical :: full
      real(ReKi) :: gmax, dmax
      full = .true.
      if (present(coreOnly)) full = .not. coreOnly
      ndiff = 0; gmax = 0.0_ReKi; dmax = 0.0_ReKi
      if (A%errStat /= B%errStat) ndiff = ndiff + 1
      do iW = 1, p%nWings
         do j = 1, p%W(iW)%nSpan
            if (.not. same(A%z%W(iW)%Gamma_LL(j), B%z%W(iW)%Gamma_LL(j))) ndiff = ndiff + 1
            if (.not. same(A%m%W(iW)%alpha_LL(j), B%m%W(iW)%alpha_LL(j))) ndiff = ndiff + 1
            if (.not. same(A%m%W(iW)%Vreln_LL(j), B%m%W(iW)%Vreln_LL(j))) ndiff = ndiff + 1
            do k = 1, 3
               if (.not. full) exit
               if (.not. same(A%m%W(iW)%Vind_CP(k,j), B%m%W(iW)%Vind_CP(k,j))) ndiff = ndiff + 1
               if (.not. same(A%m%W(iW)%Vtot_CP(k,j), B%m%W(iW)%Vtot_CP(k,j))) ndiff = ndiff + 1
            end do
            if (.not. ieee_is_nan(A%z%W(iW)%Gamma_LL(j))) gmax = max(gmax, abs(A%z%W(iW)%Gamma_LL(j)))
            if (.not. (ieee_is_nan(A%z%W(iW)%Gamma_LL(j)) .or. ieee_is_nan(B%z%W(iW)%Gamma_LL(j)))) &
               dmax = max(dmax, abs(A%z%W(iW)%Gamma_LL(j) - B%z%W(iW)%Gamma_LL(j)))
         end do
         do j = 1, p%W(iW)%nSpan + 1
            if (.not. full) exit
            do k = 1, 3
               if (.not. same(A%m%W(iW)%Vind_LL(k,j), B%m%W(iW)%Vind_LL(k,j))) ndiff = ndiff + 1
            end do
         end do
      end do
      rel = 0.0_ReKi
      if (gmax > 0.0_ReKi) rel = dmax/gmax
   end subroutine compare

   subroutine dump(p, r)
      type(FVW_ParameterType), intent(in) :: p
      type(res_t),             intent(in) :: r
      integer :: iW
      if (.not. dumping) return
      write(dumpunit) r%errStat, r%iters
      do iW = 1, p%nWings
         write(dumpunit) r%z%W(iW)%Gamma_LL, r%m%W(iW)%alpha_LL, r%m%W(iW)%Vreln_LL, r%m%W(iW)%Vind_LL
      end do
   end subroutine dump

   !> One configuration: the pin, the pin with FMA, and every port variant.
   subroutine one_case(c, label)
      type(cfg_t),  intent(in) :: c
      character(*), intent(in) :: label
      type(FVW_ParameterType)       :: p
      type(FVW_ContinuousStateType) :: x
      type(FVW_MiscVarType)         :: m0
      type(FVW_ConstraintStateType) :: zprev
      type(res_t) :: R0, RF, RV, R0b
      logical :: core
      real(ReKi) :: relF, bnd
      integer :: v, nd, dit
      real(ReKi) :: rel
      call setup(c, p, x, m0, zprev)
      ncase = ncase + 1
      call solve(S_PIN, c, p, x, m0, zprev, R0)
      ! ref vs ref: the pin twice on the same inputs must agree bitwise, or the case is not a valid oracle
      call solve(S_PIN, c, p, x, m0, zprev, R0b)
      call compare(p, R0, R0b, nd, rel)
      call teardown_res(R0b)
      core = .false.
      if (nd > 0 .or. R0b%iters /= R0%iters) then
         nrefref = nrefref + 1
         core = .true.
         print '(a,a,a,i0,a)', 'NOTE pin not reproducible (ref vs ref) in case "', label, '": ', nd, &
               ' values differ between two pin solves (uninitialised reads in the pin); gated on Gamma_LL, alpha_LL,'// &
               ' Vreln_LL, ErrStat and the iteration count only'
      end if
      ithist(max(0, min(31, R0%iters))) = ithist(max(0, min(31, R0%iters))) + 1
      call solve(S_PINFMA, c, p, x, m0, zprev, RF)
      call compare(p, R0, RF, nd, rel)
      ecmp = ecmp + 1; dit = RF%iters - R0%iters
      if (dit /= 0) then
         eitd = eitd + 1; eitdmax = max(eitdmax, abs(dit)); eitdhist(max(-30,min(30,dit))) = eitdhist(max(-30,min(30,dit))) + 1
      end if
      erel = max(erel, rel)
      relF = rel
      do v = 1, NV
         call ofk_cs_set_variant(int(vmode(v),c_int), int(vfma(v),c_int), int(vpar(v),c_int))
         call solve(S_PORT, c, p, x, m0, zprev, RV)
         call compare(p, R0, RV, nd, rel, core)
         if (.not. core) call dump(p, RV)
         vcmp(v) = vcmp(v) + 1
         dit = RV%iters - R0%iters
         if (dit /= 0) then
            vitd(v) = vitd(v) + 1; vitdmax(v) = max(vitdmax(v), abs(dit))
            if (.not. vexact(v)) vitdhist(max(-30,min(30,dit))) = vitdhist(max(-30,min(30,dit))) + 1
         end if
         vrel(v) = max(vrel(v), rel)
         if (.not. vexact(v) .and. .not. core .and. irec > 0) then   ! recorded case: reported on its own
            bnd = max(1.0e-12_ReKi, 1.1_ReKi*relF)
            rec_ratio(irec) = rel/bnd; rec_relF1(irec) = rel; rec_relF(irec) = relF; rec_dit(irec) = dit
         else if (.not. vexact(v) .and. .not. core) then   ! claim 3 with FMA (architecture.md §4): <= max(1e-12, 1.1 x pin's FMA spread), same case
            bnd = max(1.0e-12_ReKi, 1.1_ReKi*relF)
            c3n = c3n + 1
            c3worst = max(c3worst, rel/bnd)
            if (rel > bnd .or. dit /= 0) then
               c3fail = c3fail + 1
               print '(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,l1,a,es10.3,a,es10.3,a,es10.3,a,i0)', 'NOTE F1 outside the claim-3 bound (', &
                  label, ') seed=', seed, ' reg=', c%reg, ' degen=', c%degen, ' nWings=', c%nWings, ' iNWStart=', c%iNWStart, &
                  ' onLL=', c%onLL, ' epsScale=', c%epsScale, ': F1 ', rel, ' pin+FMA ', relF, &
                  ' replay state (rngs before rand_cfg; random cases only)=', case_rng
            end if
         end if
         if (nd > 0) vbit(v) = vbit(v) + 1
         if (vexact(v) .and. (nd > 0 .or. dit /= 0)) then
            nfail = nfail + 1
            if (nfail <= 20) print '(a,a,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,es10.3)', 'FAIL ', trim(vname(v)), ' (', label, &
               ') reg=', c%reg, ' degen=', c%degen, ' iNWStart=', c%iNWStart, ' differing values=', nd, &
               ' iters pin/port=', R0%iters*1000 + RV%iters, ' rel=', rel
         end if
         call ofk_cs_free(RV%m%Sgmt)
      end do
      call ofk_cs_set_variant(-1_c_int, -1_c_int, -1_c_int)
      call teardown(p, x, m0, zprev)
      call teardown_res(R0); call teardown_res(RF); call teardown_res(RV)
   end subroutine one_case

   subroutine teardown_res(r)
      type(res_t), intent(inout) :: r
      integer :: ErrStat
      character(ErrMsgLen) :: ErrMsg
      call FVW_DestroyMisc(r%m, ErrStat, ErrMsg)
      call FVW_DestroyConstrState(r%z, ErrStat, ErrMsg)
   end subroutine teardown_res

   !> Replay the recorded F1 claim-3 failures (kept as reproducible inputs; F1 is opt-in, OFK_CS_FMA=1). The exact
   !! variants are gated on them like on any case; F1's ratio is reported, not gated.
   subroutine recorded_cases()
      integer(8) :: saved
      integer :: k
      saved = rngs
      do k = 1, NREC
         if (rec_rng(k) == 0_8) cycle
         rngs = rec_rng(k); case_rng = rec_rng(k); irec = k
         call one_case(rand_cfg(), 'recorded F1 failure, '//trim(rec_name(k)))
      end do
      irec = 0; rngs = saved
   end subroutine recorded_cases

   subroutine fixed_cases()
      type(cfg_t) :: c
      integer :: reg, d, nw, ins
      do reg = 0, 4                ! every regularization (C2, C5), BAR-like default size
         c = cfg_t(); c%reg = reg; call one_case(c, 'branches')
         c%onLL = .true.; call one_case(c, 'branches, CP on bound segment')
         c%atCP = .false.; call one_case(c, 'branches, InductionAtCP=F')
      end do
      do d = 1, 7                  ! degenerate geometry (C3, C7), two regularizations
         do reg = 1, 3, 2
            c = cfg_t(); c%degen = d; c%reg = reg; call one_case(c, 'degenerate')
         end do
      end do
      do nw = 1, 3                 ! wings with unequal nSpan, including nSpan = 1 (C8)
         c = cfg_t(); c%nWings = nw; c%nSpan(1:3) = [1, 9, 4]; call one_case(c, 'wings')
         c%nSpan(1:3) = [13, 2, 7]; call one_case(c, 'wings')
      end do
      do ins = 1, 3                ! iNWStart
         c = cfg_t(); c%iNWStart = ins; c%nNWMax = ins + 2; c%nNW = c%nNWMax; call one_case(c, 'iNWStart')
         c%nNW = ins - 1; call one_case(c, 'iNWStart, short NW (empty LL #b)')
      end do
      c = cfg_t(); c%nWings = 2; c%nSpan(1:3) = [0, 5, 3]; call one_case(c, 'a wing without control points')
      c = cfg_t(); c%nWings = 1; c%nSpan(1:3) = [0, 0, 0]; call one_case(c, 'no control point (nCP_tot = 0, ad_Sphere_OLAF)')
      c = cfg_t(); c%firstCall = .true.; call one_case(c, 'FirstCall guess')           ! C10
      c = cfg_t(); c%gscale = 0.0_ReKi;  call one_case(c, 'GammaScale = 0')
      c = cfg_t(); c%gscale = 0.5_ReKi;  call one_case(c, 'GammaScale = 0.5')
      c = cfg_t(); c%maxIter = 3; c%conv = 1.0e-12_ReKi; call one_case(c, 'iteration cap reached')
      c = cfg_t(); c%nSpan(1:3) = 29; call one_case(c, 'BAR size (87 targets)')
      c%nSpan(1:3) = 50; call one_case(c, '150 targets')
   end subroutine fixed_cases

   ! ===================================================================================== C4
   !> Knife edge of the convergence test: the two adjacent doubles of ConvCrit at which the pin's iteration
   !! count changes (bisection; the iterate sequence does not depend on ConvCrit until it stops). At both,
   !! the exact variants must reproduce the pin's count; F1 and the pin with FMA are recorded.
   subroutine knife_cases()
      type(cfg_t) :: c
      integer :: k, reg
      do k = 1, 8
         c = cfg_t(); reg = mod(k, 5); c%reg = reg; c%relax = 0.1_ReKi + 0.05_ReKi*real(k, ReKi)
         c%nSpan(1:3) = [8 + k, 6, 7 - mod(k, 3)]
         call knife(c)
      end do
   end subroutine knife_cases

   subroutine knife(c0)
      type(cfg_t), intent(in) :: c0
      type(cfg_t) :: c
      type(FVW_ParameterType)       :: p
      type(FVW_ContinuousStateType) :: x
      type(FVW_MiscVarType)         :: m0
      type(FVW_ConstraintStateType) :: zprev
      type(res_t) :: R0
      real(ReKi) :: lo, hi, mid
      integer :: itlo, ithi, itmid, n, side, k0
      c = c0
      call setup(c, p, x, m0, zprev)
      ! the boundary between k0 = the count at ConvCrit = c%conv and more iterations (monotone in ConvCrit)
      lo = 1.0e-14_ReKi; hi = c%conv
      p%CircSolvConvCrit = hi; call solve(S_PIN, c, p, x, m0, zprev, R0); k0 = R0%iters; ithi = k0
      p%CircSolvConvCrit = lo; call solve(S_PIN, c, p, x, m0, zprev, R0); itlo = R0%iters
      n = 0
      do while (nearest(lo, 1.0_ReKi) < hi .and. n < 200 .and. itlo > k0)
         mid = 0.5_ReKi*(lo + hi)
         if (mid <= lo .or. mid >= hi) exit
         p%CircSolvConvCrit = mid; call solve(S_PIN, c, p, x, m0, zprev, R0); itmid = R0%iters
         if (itmid > k0) then; lo = mid; itlo = itmid; else; hi = mid; ithi = itmid; end if
         n = n + 1
      end do
      do side = 1, 2
         c%conv = merge(lo, hi, side == 1)
         p%CircSolvConvCrit = c%conv
         call knife_one(c, p, x, m0, zprev)
      end do
      print '(a,i0,a,es24.16,a,i0,a,i0,a,i0)', 'knife-edge reg=', c%reg, ': ConvCrit ', lo, ' -> iters ', itlo, &
            ' | next double -> ', ithi, ' (bisection steps ', n
      call teardown(p, x, m0, zprev); call teardown_res(R0)
   end subroutine knife

   subroutine knife_one(c, p, x, m0, zprev)
      type(cfg_t),                   intent(in) :: c
      type(FVW_ParameterType),       intent(in) :: p
      type(FVW_ContinuousStateType), intent(in) :: x
      type(FVW_MiscVarType),         intent(in) :: m0
      type(FVW_ConstraintStateType), intent(in) :: zprev
      type(res_t) :: R0, RV, RF
      integer :: v, nd, dit
      real(ReKi) :: rel
      call solve(S_PIN, c, p, x, m0, zprev, R0)
      call solve(S_PINFMA, c, p, x, m0, zprev, RF)
      call compare(p, R0, RF, nd, rel)
      kn = kn + 1
      if (RF%iters /= R0%iters) knFma = knFma + 1
      print '(a,i0,a,i0,a,es10.3)', '   knife: pin iters ', R0%iters, ', pin+FMA iters ', RF%iters, ' rel dGamma ', rel
      do v = 1, NV
         call ofk_cs_set_variant(int(vmode(v),c_int), int(vfma(v),c_int), int(vpar(v),c_int))
         call solve(S_PORT, c, p, x, m0, zprev, RV)
         call compare(p, R0, RV, nd, rel)
         dit = RV%iters - R0%iters
         if (vexact(v) .and. (nd > 0 .or. dit /= 0)) then
            nfail = nfail + 1
            print '(a,a,a,i0,a,i0)', 'FAIL knife ', trim(vname(v)), ': differing values ', nd, ' iteration delta ', dit
         end if
         if (.not. vexact(v)) then
            if (dit /= 0) knF1 = knF1 + 1
            print '(a,a,a,i0,a,es10.3)', '   knife: ', trim(vname(v)), ' iters ', RV%iters, ' rel dGamma ', rel
         end if
         call ofk_cs_free(RV%m%Sgmt)
      end do
      call ofk_cs_set_variant(-1_c_int, -1_c_int, -1_c_int)
      call teardown_res(R0); call teardown_res(RV); call teardown_res(RF)
   end subroutine knife_one

   ! ===================================================================================== C12
   !> A7: one long-lived misc mA re-used across changes of x (Eps only, r_NW only), of m%W%CP, a misc copy
   !! used in between, a re-evaluation at the first state; every result against a fresh pin solve.
   subroutine a7_sequence()
      type(cfg_t) :: c
      type(FVW_ParameterType)       :: p
      type(FVW_ContinuousStateType) :: x1, x2
      type(FVW_MiscVarType)         :: m0, mA, mC, mt
      type(FVW_ConstraintStateType) :: zprev
      type(res_t) :: R0, RA, RA1
      integer :: ErrStat, nd, step, mode
      real(ReKi) :: rel
      character(ErrMsgLen) :: ErrMsg
      integer :: na7
      na7 = 0
      do mode = 0, 1
         c = cfg_t(); c%reg = 2 + mode
         call setup(c, p, x1, m0, zprev)
         call FVW_CopyContState(x1, x2, MESH_NEWCOPY, ErrStat, ErrMsg)
         x2%W(1)%r_NW(2, :, 2) = x2%W(1)%r_NW(2, :, 2) + 0.01_ReKi
         call FVW_CopyMisc(m0, mA, MESH_NEWCOPY, ErrStat, ErrMsg)
         call ofk_cs_set_variant(int(mode,c_int), 0_c_int, -1_c_int)
         do step = 1, 7
            call FVW_DestroyMisc(mt, ErrStat, ErrMsg)
            select case (step)
            case (1, 3)   ! x1 (3: re-evaluation after a different state on another misc)
               call FVW_CopyMisc(mA, mt, MESH_NEWCOPY, ErrStat, ErrMsg)
               call solve(S_PIN, c, p, x1, mt, zprev, R0); call solve(S_PORT, c, p, x1, m0, zprev, RA, mA)
               if (step == 1) call copy_res(RA, RA1)
            case (2)      ! a misc copy (FVW_CopyMisc: new Sgmt address, new handle) on x2
               call FVW_CopyMisc(mA, mC, MESH_NEWCOPY, ErrStat, ErrMsg)
               call FVW_CopyMisc(mA, mt, MESH_NEWCOPY, ErrStat, ErrMsg)
               call solve(S_PIN, c, p, x2, mt, zprev, R0); call solve(S_PORT, c, p, x2, m0, zprev, RA, mC)
            case (4)      ! Eps_NW changed only (cache must not be reused)
               x1%W(1)%Eps_NW(1, :, 1) = 2.0_ReKi*x1%W(1)%Eps_NW(1, :, 1)
               call FVW_CopyMisc(mA, mt, MESH_NEWCOPY, ErrStat, ErrMsg)
               call solve(S_PIN, c, p, x1, mt, zprev, R0); call solve(S_PORT, c, p, x1, m0, zprev, RA, mA)
            case (5)      ! r_NW changed only
               x1%W(p%nWings)%r_NW(3, 1, 2) = x1%W(p%nWings)%r_NW(3, 1, 2) - 0.02_ReKi
               call FVW_CopyMisc(mA, mt, MESH_NEWCOPY, ErrStat, ErrMsg)
               call solve(S_PIN, c, p, x1, mt, zprev, R0); call solve(S_PORT, c, p, x1, m0, zprev, RA, mA)
            case (6)      ! a control point moved (misc input)
               mA%W(1)%CP(1, 1) = mA%W(1)%CP(1, 1) + 0.01_ReKi
               call FVW_CopyMisc(mA, mt, MESH_NEWCOPY, ErrStat, ErrMsg)
               call solve(S_PIN, c, p, x1, mt, zprev, R0); call solve(S_PORT, c, p, x1, m0, zprev, RA, mA)
            case (7)      ! the copy again, then destroyed without a free (orphaned handle)
               call FVW_CopyMisc(mC, mt, MESH_NEWCOPY, ErrStat, ErrMsg)
               call solve(S_PIN, c, p, x2, mt, zprev, R0); call solve(S_PORT, c, p, x2, m0, zprev, RA, mC)
               call FVW_DestroyMisc(mC, ErrStat, ErrMsg)
            end select
            call compare(p, R0, RA, nd, rel)
            na7 = na7 + 1
            if (nd > 0 .or. R0%iters /= RA%iters) then
               nfail = nfail + 1
               print '(a,i0,a,i0,a,i0)', 'FAIL A7 mode=', mode, ' step ', step, ': differing values ', nd
            end if
            if (step == 3) then
               call compare(p, RA1, RA, nd, rel)
               if (nd > 0) then; nfail = nfail + 1; print '(a,i0)', 'FAIL A7 re-evaluation differs from the first: ', nd; end if
            end if
         end do
         call ofk_cs_free(mA%Sgmt)
         call FVW_DestroyMisc(mA, ErrStat, ErrMsg); call FVW_DestroyMisc(mt, ErrStat, ErrMsg)
         call FVW_DestroyContState(x2, ErrStat, ErrMsg)
         call teardown(p, x1, m0, zprev)
      end do
      call ofk_cs_set_variant(-1_c_int, -1_c_int, -1_c_int)
      call teardown_res(R0); call teardown_res(RA); call teardown_res(RA1)
      print '(a,i0,a)', 'A7 (re-evaluation, misc copies, changed Eps / r_NW / CP, orphaned handle): ', na7, &
            ' solves compared with fresh pin solves'
   end subroutine a7_sequence

   subroutine copy_res(a, b)
      type(res_t), intent(in)    :: a
      type(res_t), intent(inout) :: b
      integer :: ErrStat
      character(ErrMsgLen) :: ErrMsg
      call teardown_res(b)
      b%errStat = a%errStat; b%iters = a%iters
      call FVW_CopyMisc(a%m, b%m, MESH_NEWCOPY, ErrStat, ErrMsg)
      call FVW_CopyConstrState(a%z, b%z, MESH_NEWCOPY, ErrStat, ErrMsg)
   end subroutine copy_res

   ! ===================================================================================== C11
   !> P6: 20 solves on one misc (BAR size): Kokkos allocations only in the first.
   subroutine p6_sequence()
      type(cfg_t) :: c
      type(FVW_ParameterType)       :: p
      type(FVW_ContinuousStateType) :: x
      type(FVW_MiscVarType)         :: m0, mA
      type(FVW_ConstraintStateType) :: zprev
      type(res_t) :: RA
      integer(c_long_long) :: cnt(8), a1
      integer :: k, ErrStat, it
      character(ErrMsgLen) :: ErrMsg
      c = cfg_t(); c%nSpan(1:3) = 29
      call setup(c, p, x, m0, zprev)
      call FVW_CopyMisc(m0, mA, MESH_NEWCOPY, ErrStat, ErrMsg)
      a1 = -1; it = 0
      do k = 1, 20
         x%W(1)%r_NW(1, 1, 2) = x%W(1)%r_NW(1, 1, 2) + 1.0e-3_ReKi
         call solve(S_PORT, c, p, x, m0, zprev, RA, mA)
         it = it + RA%iters
         call ofk_cs_counters(mA%Sgmt, cnt)
         if (k == 1) a1 = cnt(1)
      end do
      print '(a,i0,a,i0,a,f6.2,a,i0,a,i0,a,i0)', 'P6 (L0): allocations after solve 1: ', a1, ', after solve 20: ', &
            cnt(1), ' (', real(cnt(2))/1e6, ' MB); counter iterations ', cnt(4), ' = AFI-derived ', it, &
            '; last solve parallel_for: ', cnt(8)
      if (cnt(1) /= a1 .or. cnt(4) /= it) then; nfail = nfail + 1; print '(a)', 'FAIL P6 / iteration counter'; end if
      call ofk_cs_free(mA%Sgmt)
      call FVW_DestroyMisc(mA, ErrStat, ErrMsg)
      call teardown(p, x, m0, zprev); call teardown_res(RA)
   end subroutine p6_sequence

   ! ===================================================================================== report
   subroutine summary()
      integer :: v, k, ntot
      character(len=32) :: tag
      print '(a,i0,a,i0,a,i0)', 'cases: ', ncase, ', solves: ', nsolve, ', pin not reproducible (ref vs ref; gated on Gamma, alpha, Vreln, ErrStat, iterations): ', nrefref
      write(*, '(a)', advance='no') 'pin iteration counts (cases):'
      do k = 0, 31
         if (ithist(k) > 0) write(*, '(1x,i0,a,i0)', advance='no') k, ':', ithist(k)
      end do
      print *
      do v = 1, NV
         if (vexact(v)) then; tag = ' [must be 0: claim 1]'; else; tag = ' [info: FMA]'; end if
         print '(a,a16,a,i0,a,i0,a,i0,a,i0,a,es10.3,a)', 'variant ', vname(v), ': cases ', vcmp(v), &
               ', with any differing value ', vbit(v), ', with a different iteration count ', vitd(v), &
               ' (max |delta| ', vitdmax(v), '), max |dGamma|/max|Gamma| ', vrel(v), trim(tag)
      end do
      print '(a,i0,a,i0,a,i0,a,es10.3)', 'pin with FMA vs pin (envelope): cases ', ecmp, ', with a different iteration count ', &
            eitd, ' (max |delta| ', eitdmax, '), max |dGamma|/max|Gamma| ', erel
      print '(a,i0,a,i0,a,f8.3,a)', 'F1 claim-3 check per case (|dGamma| <= max(1e-12, 1.1 x pin+FMA spread) and same iterations): ', &
            c3n - c3fail, '/', c3n, ' within, worst ratio ', c3worst, ' [info: FMA]'
      do k = 1, NREC
         print '(a,a,a,f8.3,a,es10.3,a,es10.3,a,i0,a)', 'recorded F1 claim-3 failure (', trim(rec_name(k)), '): ratio ', &
               rec_ratio(k), ' (F1 ', rec_relF1(k), ', pin+FMA ', rec_relF(k), ', iteration delta ', rec_dit(k), &
               ') [info: F1 is opt-in; > 1 means the failure still reproduces]'
      end do
      ntot = 0
      write(*, '(a)', advance='no') 'iteration delta histogram, F1 - pin:'
      do k = -30, 30
         if (vitdhist(k) > 0) write(*, '(1x,i0,a,i0)', advance='no') k, ':', vitdhist(k)
      end do
      print *
      write(*, '(a)', advance='no') 'iteration delta histogram, pin+FMA - pin:'
      do k = -30, 30
         if (eitdhist(k) > 0) write(*, '(1x,i0,a,i0)', advance='no') k, ':', eitdhist(k)
      end do
      print *
      print '(a,i0,a,i0,a,i0,a)', 'knife-edge ConvCrit (C4): ', kn, ' solves; iteration count differs from the pin: F1 ', &
            knF1, ', pin+FMA ', knFma, ' (strict and F0: 0, checked above)'
      if (nfail == 0) then
         print '(a)', 'l0_cs: PASS (strict and F0 bitwise with the pin, same iteration counts; A7; P6)'
      else
         print '(a,i0,a)', 'l0_cs: FAIL (', nfail, ' failures)'
      end if
   end subroutine summary

end program l0_cs
