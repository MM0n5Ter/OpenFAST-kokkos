!> P1/P5 driver for unit treecode: tree build + evaluation on a wake-like particle cloud (three
!! helical filaments, alpha along the tangent, eps growing with age), N particles, M = N/2 targets
!! near the particles (the wake markers of an OLAF wake call: N_part ~ N_seg ~ 2 N_markers).
!!   perf_tc <N> [reps=3] [reg=1] [bf=1.5] [which=both|fortran|port]
!! Fortran = the pin's grow_tree_part / ui_tree_part compiled with -fopenmp (OPENMP=ON build).
!! Port   = the C API the shims call, in OFK_MODE; OFK_IV_PROFILE=1 adds the per-phase split
!!          (host pack, build, tree H2D, CP H2D, kernel, D2H, host add) printed at the end.
!! Prints one CSV row per implementation: best of reps (after one warm-up).
program perf_tc
   use NWTC_Library, only: ReKi, ErrMsgLen
   use FVW_VortexTools_ref, only: T_Tree, grow_tree_part, ui_tree_part, cut_tree
   use FVW_BiotSavart_ref, only: MINNORM
   use tc_api
   use omp_lib
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none
   real(ReKi), parameter :: Pi = 3.14159265358979323846_ReKi  ! (NWTC_Library%Pi is set by NWTC_Init only)
   integer, target :: key = 0
   integer :: n, m, reps = 3, reg = 1, i, it, b, nb, mm, ErrStat
   real(ReKi) :: bf = 1.5_ReKi, s, th, t(3), L
   real(ReKi), allocatable :: P(:,:), A(:,:), R(:), CP(:,:), UF(:,:), UK(:,:)
   real(8) :: t0, t1, t2, bF_b, bF_e, bK_b, bK_e
   character(len=64) :: arg, which = 'both'
   character(ErrMsgLen) :: ErrMsg
   type(T_Tree) :: Tree
   call get_command_argument(1, arg); read(arg,*) n
   call get_command_argument(2, arg); if (len_trim(arg) > 0) read(arg,*) reps
   call get_command_argument(3, arg); if (len_trim(arg) > 0) read(arg,*) reg
   call get_command_argument(4, arg); if (len_trim(arg) > 0) read(arg,*) bf
   call get_command_argument(5, arg); if (len_trim(arg) > 0) which = arg
   m = max(n/2, 1)
   allocate(P(3,n), A(3,n), R(n), CP(3,m), UF(3,m), UK(3,m))
   ! deterministic cloud (no RNG): three helices, a 4-turn wake 20 rotor radii... scaled to N
   nb = 3; L = 20.0_ReKi
   do i = 1, n
      b = mod(i-1, nb); mm = (i-1)/nb
      s = real(mm, ReKi) / real(max(n/nb, 1), ReKi)
      th = 2*Pi*4*s + b*2*Pi/nb
      P(:,i) = [L*s, 5*cos(th), 5*sin(th)] + 1.0e-3_ReKi*[sin(3.7_ReKi*i), cos(1.3_ReKi*i), sin(0.7_ReKi*i)]
      t = [L/(2*Pi*4), -5*sin(th), 5*cos(th)]
      A(:,i) = (1.0_ReKi + 0.5_ReKi*sin(0.37_ReKi*i)) * t / sqrt(sum(t**2)) * (L / max(n/nb, 1))
      R(i) = 0.05_ReKi + 0.3_ReKi*s
   end do
   do i = 1, m
      CP(:,i) = P(:, 2*i-1) + 0.5_ReKi*(P(:, min(2*i, n)) - P(:, 2*i-1)) + 1.0e-2_ReKi*[1.0_ReKi, -1.0_ReKi, 0.5_ReKi]
   end do
   bF_b = huge(1.0d0); bF_e = huge(1.0d0); bK_b = huge(1.0d0); bK_e = huge(1.0d0)
   if (which /= 'port') then
      do it = 0, reps
         UF = 0
         t0 = omp_get_wtime()
         call grow_tree_part(Tree, n, P, A, reg, R, 0)
         t1 = omp_get_wtime()
         call ui_tree_part(Tree, m, CP, bf, 0.0_ReKi, UF, ErrStat, ErrMsg)
         t2 = omp_get_wtime()
         call cut_tree(Tree)
         if (it > 0) then
            bF_b = min(bF_b, t1-t0); bF_e = min(bF_e, t2-t1)
         end if
      end do
      write(*,'(a,i0,a,i0,a,i0,a,f4.2,a,i0,a,3(es11.4,a))') 'csv,', n, ',', m, ',', reg, ',', bf, ',fortran-omp,', &
           omp_get_max_threads(), ',', bF_b, ',', bF_e, ',', bF_b + bF_e, ''
   end if
   if (which /= 'fortran') then
      do it = 0, reps
         UK = 0
         t0 = omp_get_wtime()
         call tc_grow(key, n, P, A, R, reg, n)
         t1 = omp_get_wtime()
         call tc_ui(key, m, CP, UK, bf, 0.0_ReKi, MINNORM)
         t2 = omp_get_wtime()
         if (it > 0) then
            bK_b = min(bK_b, t1-t0); bK_e = min(bK_e, t2-t1)
         end if
      end do
      write(*,'(a,i0,a,i0,a,i0,a,f4.2,a,i0,a,i0,a,3(es11.4,a))') 'csv,', n, ',', m, ',', reg, ',', bf, ',kokkos-mode', &
           ofk_get_mode(), ',', omp_get_max_threads(), ',', bK_b, ',', bK_e, ',', bK_b + bK_e, ''
      if (which == 'both') write(*,'(a,es10.3)') 'max |port - Fortran| / max|Fortran| = ', maxval(abs(UK-UF))/maxval(abs(UF))
      call ofk_iv_report()
      call ofk_iv_free(key)
   end if
end program perf_tc
