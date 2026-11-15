!> L0 differential test for unit induced-velocities: the pin's FVW_Subs/FVW_VortexTools/
!! FVW_BiotSavart (modules renamed *_ref at build time) vs the fork's (FVW_Subs with the
!! OPENFAST_KOKKOS shims, Kokkos-backed), in one process, through the real Fortran interfaces
!! and the real registry types (FVW_Types is shared: the port does not change it).
!!
!! Entry points: WakeInducedVelocities (-> InducedVelocitiesAll_Init/Calc/End, PackConvectingPoints,
!! PackPanelsToSegments, SegmentsToPartWrap), InducedVelocitiesAll_Init + several Calc + End
!! (the GRID/SRC pattern: one source set, many target sets), LiftingLineInducedVelocities.
!! Covered (unit.md audit list): every VelocityMethod 1..4 for the wake (iVel=1) and the lifting
!! line (iVel=2, normally hard-wired to 1), RegFunction 0..4, RegFunctionPart 0..2,
!! PartPerSegment 1..3, FWShedVorticity T/F, LastNWShed T/F (nNW<nNWMax or =), ground mirror,
!! iDepthStart 1/2/3, nNW < iDepthStart (empty NW), nFW = 0, nWings 1..3 with different nSpan,
!! InductionAtCP T/F, sentinel tails (m%Uind(:,nCPs+1:)=1000, Vind -9999/2222, NaN-poisoned
!! tails of user arrays), accumulation into a nonzero Uind, non-contiguous CPs/Uind sections.
!! A7: re-evaluation at the same x, a changed x, rollback to the first x, a misc copy
!! (FVW_CopyMisc) used alternately with the original, a copy destroyed without FVW_End.
!! P6: no device allocation after the first calls of a fixed-size sequence (ofk_iv_counters).
!!
!! Environment: as tests/biotsavart-direct/l0_util.f90 (OFK_L0_SEED, OFK_L0_CASES, OFK_L0_RTOL,
!! OFK_L0_REFREF, OFK_L0_PORTPORT, OFK_L0_STRICTFAST, OFK_L0_DUMP), plus
!!   OFK_IV_SOURCES=packed|lattice (C++ side) to select the device source path.
program l0_iv
   use NWTC_Library
   use FVW_Types
   use FVW_Subs,      only: p_Init => InducedVelocitiesAll_Init, p_Calc => InducedVelocitiesAll_Calc, &
                            p_End => InducedVelocitiesAll_End, p_Wake => WakeInducedVelocities, &
                            p_LL => LiftingLineInducedVelocities, p_Panl => T_Panl
   use FVW_VortexTools, only: p_Tree => T_Tree
   use FVW_Subs_ref,  only: r_Init => InducedVelocitiesAll_Init, r_Calc => InducedVelocitiesAll_Calc, &
                            r_End => InducedVelocitiesAll_End, r_Wake => WakeInducedVelocities, &
                            r_LL => LiftingLineInducedVelocities, r_Panl => T_Panl, &
                            FVW_InitStates, FVW_InitMiscVarsPostParam, FWnSpan, CountCPs
   use FVW_VortexTools_ref, only: r_Tree => T_Tree, VecToLattice
   use FVW_BiotSavart_ref, only: ui_seg_ref => ui_seg, ui_part_ref => ui_part_nograd
   use l0_iv_util
   implicit none
   integer :: icase
   call l0_init()
   ! Fixed cases: every method x branch x pack flag combination
   call fixed_cases()
   do icase = 1, ncases
      call one_case(rand_cfg(), 'random')
   end do
   call a7_sequence(1)
   call a7_sequence(2)   ! unit treecode: the same hazards with the particle tree as wake method
   call summary('l0_iv')
   if (nfail > 0 .and. .not. info) error stop 1

contains

   ! ---------------------------------------------------------------------------------------
   type(cfg_t) function rand_cfg() result(c)
      c%nWings = rand_int(1, 3)
      c%nSpan(1:3) = [rand_int(1, 7), rand_int(1, 7), rand_int(1, 7)]
      c%nNWMax = rand_int(1, 10); c%nNW = rand_int(0, c%nNWMax)
      c%nFWMax = rand_int(0, 6);  c%nFW = rand_int(0, c%nFWMax)
      c%FWShed = rand_int(0, 1) == 1; c%mirror = rand_int(0, 1) == 1
      c%reg = rand_int(0, 4); c%regPart = rand_int(0, 2)
      c%ppsW = rand_int(1, 3); c%ppsLL = rand_int(1, 3)
      c%velW = rand_int(1, 4); c%velLL = rand_int(1, 4)
      c%nNWFree = rand_int(0, c%nNWMax); c%nFWFree = rand_int(0, c%nFWMax)
      c%iDS = rand_int(1, 3)
   end function rand_cfg

   subroutine fixed_cases()
      type(cfg_t) :: c
      integer :: vel, reg, fws, mir, full, ids
      do vel = 1, 4
         do fws = 0, 1
            do mir = 0, 1
               do full = 0, 1
                  c = cfg_t()
                  c%nWings = 3; c%nSpan(1:3) = [4, 2, 5]
                  c%nNWMax = 6; c%nNW = merge(6, 3, full == 1)   ! LastNWShed: nNW<nNWMax or not
                  c%nFWMax = 4; c%nFW = merge(4, 2, full == 1)
                  c%FWShed = fws == 1; c%mirror = mir == 1
                  c%velW = vel; c%velLL = vel; c%reg = mod(vel + fws + mir, 5); c%regPart = mod(vel + mir, 3)
                  c%ppsW = 1 + mod(vel, 3); c%ppsLL = 1 + mod(vel + 1, 3)
                  c%nNWFree = c%nNWMax; c%nFWFree = c%nFWMax; c%iDS = 2
                  call one_case(c, 'flags')
               end do
            end do
         end do
      end do
      do reg = 0, 4   ! every segment branch, both direct methods
         do vel = 1, 3, 2
            c = cfg_t(); c%velW = vel; c%velLL = vel; c%reg = reg; c%regPart = mod(reg, 3)
            call one_case(c, 'branches')
         end do
      end do
      do ids = 1, 3   ! iDepthStart, including nNW < iDepthStart (no NW segments)
         c = cfg_t(); c%iDS = ids; c%nNW = 1; c%velLL = 1
         call one_case(c, 'iDepthStart, short NW')
         c%nFW = 0; c%nFWMax = 0
         call one_case(c, 'no far wake')
      end do
      c = cfg_t(); c%nNW = 0; c%nFW = 0; call one_case(c, 'empty wake')
      c = cfg_t(); c%nWings = 1; c%nSpan(1) = 1; call one_case(c, 'one panel span')
      c = cfg_t(); c%nWings = 3; c%nSpan(1:3) = 30; c%nNWMax = 60; c%nNW = 60; c%nFWMax = 20; c%nFW = 20
      c%velW = 1; c%velLL = 1; c%reg = 3; call one_case(c, 'skinny LL (5000+ segments)')
      c%velW = 3; c%velLL = 3; c%regPart = 2; call one_case(c, 'skinny LL particles')
   end subroutine fixed_cases

   ! ---------------------------------------------------------------------------------------
   !> Build p, x and a misc for a configuration (the same random lattice for both sides).
   subroutine setup(c, p, x, m)
      type(cfg_t),                   intent(in)    :: c
      type(FVW_ParameterType),       intent(inout) :: p
      type(FVW_ContinuousStateType), intent(inout) :: x
      type(FVW_MiscVarType),         intent(inout) :: m
      integer :: iW, j, k, nS, ErrStat
      character(ErrMsgLen) :: ErrMsg
      real(ReKi) :: r(3)
      p%nWings = c%nWings
      allocate(p%W(p%nWings))
      do iW = 1, p%nWings
         p%W(iW)%nSpan = c%nSpan(iW)
         nS = c%nSpan(iW)
         allocate(p%W(iW)%s_LL(nS+1), p%W(iW)%s_CP(nS))
         do j = 1, nS+1; p%W(iW)%s_LL(j) = real(j-1, ReKi)*1.3_ReKi + 0.1_ReKi*rnd(); end do
         do j = 1, nS;   p%W(iW)%s_CP(j) = 0.5_ReKi*(p%W(iW)%s_LL(j) + p%W(iW)%s_LL(j+1)); end do
      end do
      p%nNWMax = c%nNWMax; p%nFWMax = c%nFWMax; p%nNWFree = c%nNWFree; p%nFWFree = c%nFWFree
      p%FWShedVorticity = c%FWShed
      p%ShearModel = merge(1, 0, c%mirror)
      p%RegFunction = c%reg; p%RegFunctionPart = c%regPart
      p%VelocityMethod(1) = c%velW; p%VelocityMethod(2) = c%velLL
      p%PartPerSegment(1) = c%ppsW; p%PartPerSegment(2) = c%ppsLL
      p%TreeBranchFactor = 1.5_ReKi
      p%iNWStart = 2
      call FVW_InitStates(x, p, ErrStat, ErrMsg)
      do iW = 1, p%nWings
         do k = 1, size(x%W(iW)%r_NW,3)
            do j = 1, size(x%W(iW)%r_NW,2)
               r = [real(k,ReKi)*0.8_ReKi, real(j,ReKi)*1.1_ReKi + 3.0_ReKi*iW, 2.0_ReKi + 0.3_ReKi*iW]
               x%W(iW)%r_NW(:,j,k) = r + 0.25_ReKi*[rnd(), rnd(), rnd()]
            end do
         end do
         do k = 1, size(x%W(iW)%r_FW,3)
            do j = 1, size(x%W(iW)%r_FW,2)
               r = [real(c%nNWMax + k,ReKi)*0.8_ReKi, real(j,ReKi)*4.0_ReKi + 3.0_ReKi*iW, 2.5_ReKi]
               x%W(iW)%r_FW(:,j,k) = r + 0.25_ReKi*[rnd(), rnd(), rnd()]
            end do
         end do
         x%W(iW)%Gamma_NW = 0.0_ReKi; x%W(iW)%Gamma_FW = 0.0_ReKi
         do k = 1, size(x%W(iW)%Gamma_NW,2); do j = 1, size(x%W(iW)%Gamma_NW,1)
            x%W(iW)%Gamma_NW(j,k) = 2.0_ReKi*rnd() - 1.0_ReKi
            x%W(iW)%Eps_NW(:,j,k) = 0.02_ReKi + 0.5_ReKi*[rnd(), rnd(), rnd()]
         end do; end do
         do k = 1, size(x%W(iW)%Gamma_FW,2); do j = 1, size(x%W(iW)%Gamma_FW,1)
            x%W(iW)%Gamma_FW(j,k) = 2.0_ReKi*rnd() - 1.0_ReKi
            x%W(iW)%Eps_FW(:,j,k) = 0.02_ReKi + 0.5_ReKi*[rnd(), rnd(), rnd()]
         end do; end do
      end do
      call init_misc(c, p, x, m)
   end subroutine setup

   subroutine init_misc(c, p, x, m)
      type(cfg_t),                   intent(in)    :: c
      type(FVW_ParameterType),       intent(in)    :: p
      type(FVW_ContinuousStateType), intent(in)    :: x
      type(FVW_MiscVarType),         intent(inout) :: m
      integer :: iW, nS, ErrStat
      character(ErrMsgLen) :: ErrMsg
      allocate(m%W(p%nWings))
      do iW = 1, p%nWings
         nS = p%W(iW)%nSpan
         allocate(m%W(iW)%Vind_NW(3, nS+1, p%nNWMax+1), m%W(iW)%Vind_FW(3, FWnSpan+1, p%nFWMax+1))
         allocate(m%W(iW)%CP(3, nS), m%W(iW)%r_LL(3, nS+1, 2), m%W(iW)%Vind_CP(3, nS), m%W(iW)%Vind_LL(3, nS+1))
         m%W(iW)%Vind_NW = -1.0_ReKi; m%W(iW)%Vind_FW = -1.0_ReKi
         m%W(iW)%r_LL(:,:,1) = x%W(iW)%r_NW(:,:,1)
         m%W(iW)%r_LL(:,:,2) = x%W(iW)%r_NW(:,:,1) - 0.1_ReKi
         m%W(iW)%CP(:,:) = 0.5_ReKi*(x%W(iW)%r_NW(:,1:nS,1) + x%W(iW)%r_NW(:,2:nS+1,1)) + 0.05_ReKi
      end do
      call FVW_InitMiscVarsPostParam(p, m, ErrStat, ErrMsg)
      m%nNW = c%nNW; m%nFW = c%nFW
   end subroutine init_misc

   subroutine teardown(p, x, m)
      type(FVW_ParameterType),       intent(inout) :: p
      type(FVW_ContinuousStateType), intent(inout) :: x
      type(FVW_MiscVarType),         intent(inout) :: m
      integer :: ErrStat
      character(ErrMsgLen) :: ErrMsg
      call FVW_DestroyParam(p, ErrStat, ErrMsg)
      call FVW_DestroyContState(x, ErrStat, ErrMsg)
      call FVW_DestroyMisc(m, ErrStat, ErrMsg)
   end subroutine teardown

   ! ---------------------------------------------------------------------------------------
   !> One configuration: wake, Init + 2 Calcs + End, lifting line (both InductionAtCP).
   subroutine one_case(c, label)
      type(cfg_t),  intent(in) :: c
      character(*), intent(in) :: label
      type(FVW_ParameterType)       :: p
      type(FVW_ContinuousStateType) :: x
      type(FVW_MiscVarType)         :: mA, mB
      character(len=160) :: lab
      integer :: k
      call setup(c, p, x, mA)
      call FVW_CopyMisc(mA, mB, MESH_NEWCOPY, k, lab)
      write(lab,'(a,a,i0,a,i0,a,i0,a,i0,a,l1,a,l1,a,i0,a,i0)') trim(label), ' velW=', c%velW, ' velLL=', c%velLL, &
           ' reg=', c%reg, '/', c%regPart, ' FWShed=', c%FWShed, ' mirror=', c%mirror, ' nNW=', c%nNW, ' nFW=', c%nFW
      call check_wake(trim(lab), c, p, x, mA, mB)
      call check_init_calc(trim(lab), c, p, x, mA, mB)
      call check_ll(trim(lab), c, p, x, mA, mB, .true.,  c%iDS)
      call check_ll(trim(lab), c, p, x, mA, mB, .false., 1)
      call teardown(p, x, mA)
      call teardown(p, x, mB)
   end subroutine one_case

   ! --- the A side and the B side of a comparison (l0_iv_util: oracle | port, modes) ---------
   subroutine wake_role(role, p, x, m)
      integer, intent(in) :: role
      type(FVW_ParameterType),       intent(in)    :: p
      type(FVW_ContinuousStateType), intent(in)    :: x
      type(FVW_MiscVarType),         intent(inout) :: m
      integer :: ErrStat
      character(ErrMsgLen) :: ErrMsg
      if (use_ref(role)) then
         call r_Wake(p, x, m, ErrStat, ErrMsg)
      else
         call p_Wake(p, x, m, ErrStat, ErrMsg)
      end if
   end subroutine wake_role

   subroutine ll_role(role, p, x, m, atCP, iDS)
      integer, intent(in) :: role, iDS
      logical, intent(in) :: atCP
      type(FVW_ParameterType),       intent(in)    :: p
      type(FVW_ContinuousStateType), intent(in)    :: x
      type(FVW_MiscVarType),         intent(inout) :: m
      integer :: ErrStat
      character(ErrMsgLen) :: ErrMsg
      if (use_ref(role)) then
         call r_LL(p, x, atCP, iDS, m, ErrStat, ErrMsg)
      else
         call p_LL(p, x, atCP, iDS, m, ErrStat, ErrMsg)
      end if
   end subroutine ll_role

   !> Init, then Calc on two target sets (CPs1 contiguous, CPs2 a non-contiguous section), End.
   subroutine initcalc_role(role, p, x, m, CPs1, U1, CPs2, U2, n1, n2)
      integer, intent(in) :: role, n1, n2
      type(FVW_ParameterType),       intent(in), target :: p
      type(FVW_ContinuousStateType), intent(in)    :: x
      type(FVW_MiscVarType),         intent(inout), target :: m
      real(ReKi), intent(in)    :: CPs1(:,:), CPs2(:,:)
      real(ReKi), intent(inout) :: U1(:,:), U2(:,:)
      type(r_Tree) :: tr; type(r_Panl) :: pr
      type(p_Tree) :: tp; type(p_Panl) :: pp
      integer :: ErrStat
      character(ErrMsgLen) :: ErrMsg
      if (use_ref(role)) then
         call r_Init(p, x, m, m%Sgmt, m%Part, tr, pr, ErrStat, ErrMsg, allocPart=.false.)
         call r_Calc(CPs1, n1, U1, p, m%Sgmt, m%Part, tr, pr, ErrStat, ErrMsg)
         call r_Calc(CPs2(:,1:2*n2:2), n2, U2(:,1:2*n2:2), p, m%Sgmt, m%Part, tr, pr, ErrStat, ErrMsg)
         call r_End(p, tr, m%Part, pr, ErrStat, ErrMsg, deallocPart=.false.)
      else
         call p_Init(p, x, m, m%Sgmt, m%Part, tp, pp, ErrStat, ErrMsg, allocPart=.false.)
         call p_Calc(CPs1, n1, U1, p, m%Sgmt, m%Part, tp, pp, ErrStat, ErrMsg)
         call p_Calc(CPs2(:,1:2*n2:2), n2, U2(:,1:2*n2:2), p, m%Sgmt, m%Part, tp, pp, ErrStat, ErrMsg)
         call p_End(p, tp, m%Part, pp, ErrStat, ErrMsg, deallocPart=.false.)
      end if
   end subroutine initcalc_role

   ! --- sum_j |c_j| per target, from the oracle called once per source ---------------------
   !> S(k,i) = sum over the sources of |contribution_k| at CPs(:,i) of the source set the misc
   !! holds after the call (method 1: segments, 2 and 3: particles). Negative = exact comparison
   !! (method 4, the segment tree: both sides run the same Fortran). Method 2 (particle tree) is
   !! ported by unit treecode (tests/treecode): bitwise in strict (rtol 0), and strict vs fast at the
   !! strictfast bound (the same tree; a rounding-level change of a node center could flip one
   !! acceptance test, which tests/treecode's R2 gate covers; none seen in L0).
   subroutine scale_sum(method, m, CPs, n, S)
      integer,               intent(in)  :: method, n
      type(FVW_MiscVarType), intent(in)  :: m
      real(ReKi),            intent(in)  :: CPs(:,:)
      real(ReKi),            intent(out) :: S(:,:)
      real(ReKi), allocatable :: Ut(:,:)
      integer :: j
      S = 0.0_ReKi
      if (n <= 0) return
      if (method /= 1 .and. method /= 2 .and. method /= 3) then
         S = -1.0_ReKi; return
      end if
      allocate(Ut(3, n))
      if (method == 1) then
         do j = 1, m%Sgmt%nAct
            Ut = 0.0_ReKi
            call ui_seg_ref(1, n, CPs, j, j, m%Sgmt%Points, m%Sgmt%Connct, m%Sgmt%Gamma, m%Sgmt%RegFunction, m%Sgmt%Epsilon, Ut)
            S(1:3,1:n) = S(1:3,1:n) + abs(Ut)
         end do
      else
         do j = 1, m%Part%nAct
            Ut = 0.0_ReKi
            call ui_part_ref(n, CPs, 1, m%Part%P(:,j:j), m%Part%Alpha(:,j:j), m%Part%RegFunction, m%Part%RegParam(j:j), Ut)
            S(1:3,1:n) = S(1:3,1:n) + abs(Ut)
         end do
      end if
   end subroutine scale_sum

   ! ---------------------------------------------------------------------------------------
   subroutine check_wake(label, c, p, x, mA, mB)
      character(*), intent(in) :: label
      type(cfg_t),  intent(in) :: c
      type(FVW_ParameterType),       intent(in)    :: p
      type(FVW_ContinuousStateType), intent(in)    :: x
      type(FVW_MiscVarType),         intent(inout) :: mA, mB
      real(ReKi), allocatable :: S(:,:), SL(:,:,:)
      integer :: nCPs, iW, iHead, nNWEff, nFWEff
      call wake_role(roleA, p, x, mA)
      call wake_role(roleB, p, x, mB)
      nNWEff = min(mA%nNW, p%nNWFree); nFWEff = min(mA%nFW, p%nFWFree)
      nCPs = CountCPs(p, nNWEff, nFWEff)
      allocate(S(3, size(mA%Uind,2))); S = 0.0_ReKi
      call scale_sum(p%VelocityMethod(1), mA, mA%CPs, nCPs, S)
      call cmp2('wake Uind | '//label, p%VelocityMethod(1), mA%Uind, mB%Uind, S, nCPs)
      call cmp2('wake CPs | '//label, CAT_PACK, mA%CPs, mB%CPs)
      call cmp_sgmt('wake Sgmt/Part | '//label, mA, mB)
      ! Vind lattices: the unpacked Uind (scale unpacked the same way)
      iHead = 1
      do iW = 1, p%nWings
         allocate(SL, mold=mA%W(iW)%Vind_NW); SL = 0.0_ReKi
         if (S(1,1) < 0) SL = -1.0_ReKi
         if (nNWEff+1 >= 1 .and. S(1,1) >= 0) call VecToLattice(S, 1, SL(:,:,1:nNWEff+1), iHead)
         call cmp3('wake Vind_NW | '//label, p%VelocityMethod(1), mA%W(iW)%Vind_NW, mB%W(iW)%Vind_NW, SL)
         deallocate(SL)
      end do
      do iW = 1, p%nWings
         allocate(SL, mold=mA%W(iW)%Vind_FW); SL = 0.0_ReKi
         if (S(1,1) < 0) SL = -1.0_ReKi
         if (nFWEff > 0 .and. S(1,1) >= 0) call VecToLattice(S, 1, SL(:,1:FWnSpan+1,1:nFWEff+1), iHead)
         call cmp3('wake Vind_FW | '//label, p%VelocityMethod(1), mA%W(iW)%Vind_FW, mB%W(iW)%Vind_FW, SL)
         deallocate(SL)
      end do
   end subroutine check_wake

   subroutine check_init_calc(label, c, p, x, mA, mB)
      character(*), intent(in) :: label
      type(cfg_t),  intent(in) :: c
      type(FVW_ParameterType),       intent(in)    :: p
      type(FVW_ContinuousStateType), intent(in)    :: x
      type(FVW_MiscVarType),         intent(inout) :: mA, mB
      real(ReKi), allocatable :: CPs1(:,:), CPs2(:,:), U0a(:,:), U0b(:,:), U1A(:,:), U1B(:,:), U2A(:,:), U2B(:,:)
      real(ReKi), allocatable :: S1(:,:), S2(:,:), C2c(:,:)
      integer :: n1, n2, i
      n1 = rand_int(0, 40); n2 = rand_int(1, 25)
      allocate(CPs1(3, n1+3), CPs2(3, 2*n2+1), U0a(3, n1+3), U0b(3, 2*n2+1))
      do i = 1, size(CPs1,2); CPs1(:,i) = [8.0_ReKi*rnd(), 12.0_ReKi*rnd(), 4.0_ReKi*rnd()-1.0_ReKi]; end do
      do i = 1, size(CPs2,2); CPs2(:,i) = [8.0_ReKi*rnd(), 12.0_ReKi*rnd(), 4.0_ReKi*rnd()-1.0_ReKi]; end do
      if (n1 > 0) CPs1(:,1) = x%W(1)%r_NW(:,1,1)   ! a target on a lattice point (guards)
      do i = 1, size(U0a,2); U0a(:,i) = [rnd(), -rnd(), 0.0_ReKi]; end do
      U0a(3,1) = -0.0_ReKi
      do i = 1, size(U0b,2); U0b(:,i) = [rnd(), rnd(), rnd()]; end do
      U0a(:, n1+1:) = qnan(); U0b(:, 2:2*n2:2) = qnan()   ! tails and skipped columns: untouched
      allocate(U1A, source=U0a); allocate(U1B, source=U0a); allocate(U2A, source=U0b); allocate(U2B, source=U0b)
      call initcalc_role(roleA, p, x, mA, CPs1, U1A, CPs2, U2A, n1, n2)
      call initcalc_role(roleB, p, x, mB, CPs1, U1B, CPs2, U2B, n1, n2)
      allocate(S1(3, size(U0a,2)), S2(3, size(U0b,2)), C2c(3, n2))
      S1 = 0.0_ReKi; S2 = 0.0_ReKi
      call scale_sum(p%VelocityMethod(1), mA, CPs1, n1, S1)
      C2c = CPs2(:,1:2*n2:2)
      call scale_sum(p%VelocityMethod(1), mA, C2c, n2, S2(:,1:n2))
      ! scale of the strided section, spread back to its columns
      S2(:, 1:2*n2:2) = S2(:, 1:n2)
      call cmp2('calc #1 | '//label, p%VelocityMethod(1), U1A, U1B, S1, n1)
      call cmp2('calc #2 strided | '//label, p%VelocityMethod(1), U2A, U2B, S2, 2*n2, stride=2)
      call cmp_sgmt('init Sgmt/Part | '//label, mA, mB)
   end subroutine check_init_calc

   subroutine check_ll(label, c, p, x, mA, mB, atCP, iDS)
      character(*), intent(in) :: label
      type(cfg_t),  intent(in) :: c
      logical,      intent(in) :: atCP
      integer,      intent(in) :: iDS
      type(FVW_ParameterType),       intent(in)    :: p
      type(FVW_ContinuousStateType), intent(in)    :: x
      type(FVW_MiscVarType),         intent(inout) :: mA, mB
      real(ReKi), allocatable :: CPs(:,:), S(:,:)
      integer :: iW, n, h, nS
      character(len=24) :: tag
      write(tag,'(a,l1,a,i0)') 'LL atCP=', atCP, ' iDS=', iDS
      call ll_role(roleA, p, x, mA, atCP, iDS)
      call ll_role(roleB, p, x, mB, atCP, iDS)
      ! the direct targets of the LL call, and their scale
      n = 0
      do iW = 1, p%nWings
         n = n + p%W(iW)%nSpan + merge(0, 1, atCP)
      end do
      allocate(CPs(3, n), S(3, n)); h = 0
      do iW = 1, p%nWings
         nS = p%W(iW)%nSpan
         if (atCP) then
            CPs(:, h+1:h+nS) = mA%W(iW)%CP(:, 1:nS); h = h + nS
         else
            CPs(:, h+1:h+nS+1) = mA%W(iW)%r_LL(:, 1:nS+1, 1); h = h + nS + 1
         end if
      end do
      S = 0.0_ReKi
      if (mA%Sgmt%nActP > 0) call scale_sum(p%VelocityMethod(2), mA, CPs, n, S)
      h = 0
      do iW = 1, p%nWings
         nS = p%W(iW)%nSpan
         if (atCP) then
            call cmp2(trim(tag)//' Vind_CP | '//label, p%VelocityMethod(2), mA%W(iW)%Vind_CP, mB%W(iW)%Vind_CP, S(:,h+1:h+nS))
            ! Vind_LL: linear inter/extrapolation of Vind_CP; scale: 3 x the wing's max (extrapolation)
            call cmp2(trim(tag)//' Vind_LL | '//label, p%VelocityMethod(2), mA%W(iW)%Vind_LL, mB%W(iW)%Vind_LL, &
                      spread_scale(S(:,h+1:h+nS), nS+1, 3.0_ReKi))
            h = h + nS
         else
            call cmp2(trim(tag)//' Vind_LL | '//label, p%VelocityMethod(2), mA%W(iW)%Vind_LL, mB%W(iW)%Vind_LL, S(:,h+1:h+nS+1))
            call cmp2(trim(tag)//' Vind_CP | '//label, p%VelocityMethod(2), mA%W(iW)%Vind_CP, mB%W(iW)%Vind_CP, &
                      spread_scale(S(:,h+1:h+nS+1), nS, 1.0_ReKi))
            h = h + nS + 1
         end if
      end do
      call cmp_sgmt(trim(tag)//' Sgmt/Part | '//label, mA, mB)
   end subroutine check_ll

   function spread_scale(S, n, f) result(R)
      real(ReKi), intent(in) :: S(:,:), f
      integer,    intent(in) :: n
      real(ReKi) :: R(3, n)
      if (size(S,2) == 0) then
         R = 0.0_ReKi
      else if (any(S < 0)) then
         R = -1.0_ReKi
      else
         R = f*maxval(S)
      end if
   end function spread_scale

   ! ---------------------------------------------------------------------------------------
   !> A7 + P6: one long-lived instance (mS) through a time-loop-like sequence of LL and wake
   !! calls on changing states; after every call its outputs must be bitwise those of a FRESH
   !! instance (a new FVW_CopyMisc of the initial misc, hence a new handle) on the same input,
   !! in the current mode: no stale device state. In the default role (oracle vs port) each call
   !! is also compared with the oracle at OFK_L0_RTOL.
   !!   x1; x1 again (re-evaluation); x2 = x1 with the first two depth rows changed (the
   !!   Map_LL_NW pattern); x3 = everything moved (propagation); x1 again (rollback); a copy
   !!   of the misc (FVW_CopyMisc) used alternately with the original, then destroyed without
   !!   FVW_End (orphaned handle). Then 20 steady steps: no device allocation (P6).
   subroutine a7_sequence(velW)
      integer, intent(in) :: velW
      type(cfg_t) :: c
      type(FVW_ParameterType)       :: p
      type(FVW_ContinuousStateType) :: x1, x2, x3
      type(FVW_MiscVarType)         :: m0, mS, mS2
      integer :: ErrStat, k, iW
      character(ErrMsgLen) :: ErrMsg
      integer(8) :: c0(9), c1(9)
      if (refref) return
      c = cfg_t(); c%velW = velW; c%velLL = 1; c%reg = 3; c%nNW = 6; c%nNWMax = 6; c%nFW = 3; c%nFWMax = 3
      call setup(c, p, x1, m0)
      call FVW_CopyMisc(m0, mS, MESH_NEWCOPY, ErrStat, ErrMsg)
      call FVW_CopyContState(x1, x2, MESH_NEWCOPY, ErrStat, ErrMsg)
      call FVW_CopyContState(x1, x3, MESH_NEWCOPY, ErrStat, ErrMsg)
      do iW = 1, p%nWings
         x2%W(iW)%Gamma_NW(:,1:2) = x2%W(iW)%Gamma_NW(:,1:2) + 0.125_ReKi
         x2%W(iW)%r_NW(3,:,1:2)   = x2%W(iW)%r_NW(3,:,1:2) + 0.0625_ReKi
         x3%W(iW)%r_NW = x3%W(iW)%r_NW + 0.03125_ReKi
         x3%W(iW)%r_FW = x3%W(iW)%r_FW - 0.03125_ReKi
      end do
      call a7_call('A7 LL x1', 1, p, x1, m0, mS)
      call a7_call('A7 LL x1 re-evaluation', 1, p, x1, m0, mS)
      call a7_call('A7 wake x1', 0, p, x1, m0, mS)
      call a7_call('A7 LL x2 (near rows changed)', 1, p, x2, m0, mS)
      call a7_call('A7 LL x3 (all moved)', 1, p, x3, m0, mS)
      call a7_call('A7 LL x1 (rollback)', 1, p, x1, m0, mS)
      call a7_call('A7 wake x1 (rollback)', 0, p, x1, m0, mS)
      call FVW_CopyMisc(mS, mS2, MESH_NEWCOPY, ErrStat, ErrMsg)
      call a7_call('A7 copy on x2', 1, p, x2, m0, mS2)
      call a7_call('A7 original on x1', 1, p, x1, m0, mS)
      call a7_call('A7 copy on x3', 0, p, x3, m0, mS2)
      call a7_call('A7 original on x2', 0, p, x2, m0, mS)
      call FVW_DestroyMisc(mS2, ErrStat, ErrMsg)   ! no FVW_End: its handle is orphaned
      call counters(mS%Sgmt, c0)
      do k = 1, 20
         select case (mod(k, 3))
         case (0); call a7_call('A7 steady', 2, p, x1, m0, mS)
         case (1); call a7_call('A7 steady', 2, p, x2, m0, mS)
         case default; call a7_call('A7 steady', 2, p, x3, m0, mS)
         end select
      end do
      call counters(mS%Sgmt, c1)
      write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,f9.4,a,f9.4,a,i0)') 'P6 (L0, wake VelocityMethod ', velW, &
           '): device allocations during 20 steady steps: ', &
           c1(1)-c0(1), ' (instance total ', c1(1), '); source builds ', c1(3)-c0(3), ', evaluations ', c1(4)-c0(4), &
           '; lattice MB compared ', real(c1(7)-c0(7))/1e6, ', uploaded ', real(c1(8)-c0(8))/1e6, '; live handles ', c1(9)
      if (c1(1) /= c0(1)) then
         nfail = nfail + 1; write(*,'(a)') 'FAIL [P6] device allocation after warm-up'
      end if
      call ofk_iv_free(mS%Sgmt)
      call teardown(p, x1, m0); call FVW_DestroyMisc(mS, ErrStat, ErrMsg)
      call FVW_DestroyContState(x2, ErrStat, ErrMsg); call FVW_DestroyContState(x3, ErrStat, ErrMsg)
   end subroutine a7_sequence

   !> kind 0: wake; 1: LL (iDepthStart = iNWStart+1, then 1, as in one step); 2: wake then LL x2.
   subroutine a7_call(label, kind, p, x, m0, mS)
      character(*), intent(in) :: label
      integer,      intent(in) :: kind
      type(FVW_ParameterType),       intent(in)    :: p
      type(FVW_ContinuousStateType), intent(in)    :: x
      type(FVW_MiscVarType),         intent(in)    :: m0
      type(FVW_MiscVarType),         intent(inout) :: mS
      type(FVW_MiscVarType) :: mF, mR
      integer :: ErrStat, iW, ids, step
      character(ErrMsgLen) :: ErrMsg
      do step = 1, merge(3, 2, kind == 2)
         if (kind == 0 .or. (kind == 2 .and. step == 1)) then
            ids = 0
         else
            ids = merge(p%iNWStart+1, 1, step == 1 .or. (kind == 2 .and. step == 2))
         end if
         call FVW_CopyMisc(m0, mF, MESH_NEWCOPY, ErrStat, ErrMsg)
         call FVW_CopyMisc(m0, mR, MESH_NEWCOPY, ErrStat, ErrMsg)
         if (ids == 0) then
            call wake_role(roleB, p, x, mS)
            call wake_role(roleB, p, x, mF)
            call cmp2(label//' wake Uind: long-lived vs fresh', CAT_A7, mF%Uind, mS%Uind, exact=.true.)
            if (roleA == ROLE_REF .and. roleB == ROLE_PORT) then
               call wake_role(ROLE_REF, p, x, mR)
               call cmp2(label//' wake Uind: vs oracle', CAT_A7, mR%Uind, mS%Uind)
            end if
         else
            call ll_role(roleB, p, x, mS, .true., ids)
            call ll_role(roleB, p, x, mF, .true., ids)
            do iW = 1, p%nWings
               call cmp2(label//' LL Vind_CP: long-lived vs fresh', CAT_A7, mF%W(iW)%Vind_CP, mS%W(iW)%Vind_CP, exact=.true.)
            end do
            if (roleA == ROLE_REF .and. roleB == ROLE_PORT) then
               call ll_role(ROLE_REF, p, x, mR, .true., ids)
               do iW = 1, p%nWings
                  call cmp2(label//' LL Vind_CP: vs oracle', CAT_A7, mR%W(iW)%Vind_CP, mS%W(iW)%Vind_CP)
               end do
            end if
         end if
         call ofk_iv_free(mF%Sgmt)
         call FVW_DestroyMisc(mF, ErrStat, ErrMsg)
         call FVW_DestroyMisc(mR, ErrStat, ErrMsg)
      end do
   end subroutine a7_call

end program l0_iv
