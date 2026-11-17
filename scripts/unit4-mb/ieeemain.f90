program ieeemain
   use uf
   implicit none
   integer, parameter :: n = 18623, nrep = 300
   real(ReKi), allocatable :: P(:,:), V1(:,:), V2(:,:)
   type(op_t) :: op
   integer :: i, r
   integer(8) :: c0, c1, cr
   real(ReKi) :: t1, t2
   allocate(P(3,n), V1(3,n), V2(3,n))
   do i = 1, n; P(:,i) = [real(i,ReKi)*0.01_ReKi, sin(real(i,ReKi))*60, 140 + cos(real(i,ReKi))*60]; end do
   op = op_t(8.0_ReKi, 0.6_ReKi, 0._ReKi, 1._ReKi, 0._ReKi, 0._ReKi, 0._ReKi, 1._ReKi, 0._ReKi, 0._ReKi)
   call system_clock(c0, cr)
   do r = 1, nrep; do i = 1, n; if (P(3,i) > 0) V1(:,i) = get_ieee(140._ReKi, 240._ReKi, op, P(:,i)); end do; end do
   call system_clock(c1); t1 = real(c1-c0,ReKi)/cr
   call system_clock(c0)
   do r = 1, nrep; do i = 1, n; if (P(3,i) > 0) V2(:,i) = get_plain(140._ReKi, 240._ReKi, op, P(:,i)); end do; end do
   call system_clock(c1); t2 = real(c1-c0,ReKi)/cr
   print '(a,i0,a,f8.3,a,f8.3,a,f6.1,a,i0)', 'points=', n*nrep, '  with ieee use: ', t1, ' s   without: ', t2, &
        ' s   ratio ', t1/t2, '   bit-differing components: ', count(V1 /= V2)
end program
