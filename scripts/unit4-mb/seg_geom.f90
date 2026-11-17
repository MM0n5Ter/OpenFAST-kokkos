! F10-B geometry half of ui_seg_11 (all statements before line 124 'Kv = SegGamma*...'); ok=.false. <=> an early return
subroutine seg_geom(DeltaPa, DeltaPb, RegFunction, RegParam1, ok, Kv, s, den, cp)
   real(ReKi), intent(in) :: DeltaPa(3), DeltaPb(3), RegParam1
   integer, intent(in) :: RegFunction
   logical, intent(out) :: ok; real(ReKi), intent(out) :: Kv, s, den, cp(3)
   real(ReKi) :: r_bar2, norm_a, norm_b, norm2_r0, norm2_orth, xa, ya, za, xb, yb, zb, exp_value
   ok = .false.; Kv = 0; s = 0; den = 1; cp = 0
   xa=DeltaPa(1); ya=DeltaPa(2); za=DeltaPa(3); xb=DeltaPb(1); yb=DeltaPb(2); zb=DeltaPb(3)
   norm_a = sqrt(xa*xa + ya*ya + za*za); norm_b = sqrt(xb*xb + yb*yb + zb*zb)
   den = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb)
   if (den <= PRECISION_UI) return
   cp(1) = ya*zb-za*yb; cp(2) = za*xb-xa*zb; cp(3) = xa*yb-ya*xb
   norm2_orth = cp(1)**2 + cp(2)**2 + cp(3)**2
   if (norm2_orth <= PRECISION_UI) return
   norm2_r0 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) + (za-zb)*(za-zb)
   if (norm2_r0 <= PRECISION_UI) return
   norm2_orth = norm2_orth/norm2_r0
   select case (RegFunction)
   case (idRegNone); Kv = 1.0_ReKi
   case (idRegRankine); r_bar2 = norm2_orth/ RegParam1**2; if (r_bar2<1) then; Kv=r_bar2; else; Kv=1.0_ReKi; end if
   case (idRegLambOseen); r_bar2 = norm2_orth/ RegParam1**2; exp_value = -1.25643_ReKi*r_bar2
      if (exp_value<MIN_EXP_VALUE) then; Kv = 1.0_ReKi; else; Kv = 1.0_ReKi-exp(exp_value); endif
   case (idRegVatistas); r_bar2 = norm2_orth/ RegParam1**2; Kv = r_bar2/sqrt(1.0_ReKi+r_bar2**2)
   case (idRegOffset); Kv = 1.0_ReKi; den = den + RegParam1**2*norm2_r0
   end select
   s = norm_a+norm_b; ok = .true.
end subroutine
end module
