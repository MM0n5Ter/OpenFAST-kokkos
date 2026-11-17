program f10mb
   use bs
   implicit none
   integer :: nB, nS, nD, nCP, nIt, it, i, j, d, k, q, nbad, reg
   real(ReKi), allocatable :: R(:,:,:), CP(:,:), Eps(:,:), G(:), V1(:,:), V2(:,:), V3(:,:), A(:,:,:)
   real(ReKi), allocatable :: cK(:,:,:,:), cS(:,:,:,:), cD(:,:,:,:), cC(:,:,:,:,:)
   logical, allocatable :: cOk(:,:,:,:)
   real(ReKi) :: P(3,4), U(3), Kv, t0, t1, ta, tbg, tbi, tcg, tci, z, rr
   integer(8) :: c0, c1, cr
   character(32) :: arg
   integer, parameter :: sa(4) = [1,3,2,4], sb(4) = [2,4,3,1]   ! segment order of ui_quad_n1: 1-2, 3-4, 2-3, 4-1
   nB = 3; nS = 29; nD = 2; nIt = 20; reg = 3
   if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read(arg,*) nS; end if
   if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read(arg,*) reg; end if
   nCP = nB*nS
   allocate(R(3, nS+1, nB*(nD+1)), CP(3,nCP), Eps(nS, nB*nD), G(nCP), V1(3,nCP), V2(3,nCP), V3(3,nCP), A(3,nCP,nCP*nD*4))
   allocate(cK(4,nD,nCP,nCP), cS(4,nD,nCP,nCP), cD(4,nD,nCP,nCP), cC(3,4,nD,nCP,nCP), cOk(4,nD,nCP,nCP))
   ! geometry: 3 blades, span 3..100 m, wake panels 0.6 m apart downstream; CPs mid-span, 0.25 chord ahead
   do k = 1, nB; do j = 1, nS+1; do d = 1, nD+1
      rr = 3 + 97*real(j-1,ReKi)/nS; z = 2.0944_ReKi*(k-1)
      R(:, j, (k-1)*(nD+1)+d) = [0.6_ReKi*(d-1), rr*cos(z) + 0.01_ReKi*(d-1), rr*sin(z)]
   end do; end do; end do
   do k = 1, nB; do j = 1, nS
      i = (k-1)*nS + j; CP(:, i) = 0.5_ReKi*(R(:, j, (k-1)*(nD+1)+1) + R(:, j+1, (k-1)*(nD+1)+1)) - [0.3_ReKi, 0._ReKi, 0._ReKi]
      G(i) = 10 + sin(real(i,ReKi))
      do d = 1, nD; Eps(j, (k-1)*nD+d) = 0.25_ReKi*(1 + 0.1_ReKi*d); end do
   end do; end do
   ! (a) pin, nIt iterations: FVW_Wings.f90:421-469 loop shape, ui_quad_n1 per ring per target wing
   call system_clock(c0, cr)
   do it = 1, nIt
      V1 = 0
      do k = 1, nB; do j = 1, nS; do d = 1, nD
         call ring(k, j, d)
         do q = 1, nB
            call ui_quad_n1(CP(1:3,(q-1)*nS+1:q*nS), nS, P(:,1), P(:,2), P(:,3), P(:,4), G((k-1)*nS+j), reg, Eps(j,(k-1)*nD+d), V1(1:3,(q-1)*nS+1:q*nS))
         end do
      end do; end do; end do
   end do
   call system_clock(c1); ta = real(c1-c0,ReKi)/cr/nIt
   ! (b) F10-B: cache once
   call system_clock(c0)
   do i = 1, nCP; do k = 1, nB; do j = 1, nS; do d = 1, nD
      call ring(k, j, d)
      do q = 1, 4
         call seg_geom(CP(:,i)-P(:,sa(q)), CP(:,i)-P(:,sb(q)), reg, Eps(j,(k-1)*nD+d), cOk(q,d,(k-1)*nS+j,i), cK(q,d,(k-1)*nS+j,i), &
                       cS(q,d,(k-1)*nS+j,i), cD(q,d,(k-1)*nS+j,i), cC(:,q,d,(k-1)*nS+j,i))
      end do
   end do; end do; end do; end do
   call system_clock(c1); tbg = real(c1-c0,ReKi)/cr
   call system_clock(c0)
   do it = 1, nIt   ! replay, target-outer (same per-target order as the pin: ring, depth, segment)
      do i = 1, nCP
         V2(:,i) = 0
         do j = 1, nCP; do d = 1, nD; do q = 1, 4
            if (cOk(q,d,j,i)) then
               Kv = G(j)*fourpi_inv*cK(q,d,j,i)*cS(q,d,j,i)/cD(q,d,j,i)
               V2(:,i) = V2(:,i) + Kv*cC(:,q,d,j,i)
            else
               V2(:,i) = V2(:,i) + 0.0_ReKi
            end if
         end do; end do; end do
      end do
   end do
   call system_clock(c1); tbi = real(c1-c0,ReKi)/cr/nIt
   nbad = count(V1 /= V2)
   ! (c) F10-R1: matrix with Gamma=1, then mat-vec
   call system_clock(c0)
   A = 0
   do k = 1, nB; do j = 1, nS; do d = 1, nD
      call ring(k, j, d)
      do i = 1, nCP; do q = 1, 4
         call ui_seg_11(CP(:,i)-P(:,sa(q)), CP(:,i)-P(:,sb(q)), 1.0_ReKi, reg, Eps(j,(k-1)*nD+d), U)
         A(:, i, (((k-1)*nS+j)-1)*nD*4 + (d-1)*4 + q) = U
      end do; end do
   end do; end do; end do
   call system_clock(c1); tcg = real(c1-c0,ReKi)/cr
   call system_clock(c0)
   do it = 1, nIt
      V3 = 0
      do j = 1, nCP; do q = 1, nD*4
         V3(:, :) = V3(:, :) + G(j)*A(:, :, (j-1)*nD*4+q)
      end do; end do
   end do
   call system_clock(c1); tci = real(c1-c0,ReKi)/cr/nIt
   print '(a,i4,a,i9,a,i0)', 'nCP=', nCP, '  pairs/iter=', nCP*nCP*nD*4, '  RegFunction=', reg
   print '(a,es10.3,a)', '(a) pin ui_seg_11 per iteration      ', ta, ' s'
   print '(a,es10.3,a,es10.3,a,i0)', '(b) F10-B  cache build ', tbg, ' s, replay per iteration ', tbi, ' s; bit-differing components vs pin: ', nbad
   print '(a,es10.3,a,es10.3,a,es10.3)', '(c) F10-R1 matrix build', tcg, ' s, mat-vec per iteration ', tci, ' s; max rel diff vs pin ', &
         maxval(abs(V3-V1))/maxval(abs(V1))
contains
   subroutine ring(k, j, d)
      integer, intent(in) :: k, j, d
      P(:,1) = R(:, j, (k-1)*(nD+1)+d); P(:,4) = R(:, j, (k-1)*(nD+1)+d+1)
      P(:,2) = R(:, j+1, (k-1)*(nD+1)+d); P(:,3) = R(:, j+1, (k-1)*(nD+1)+d+1)
   end subroutine
end program
