! Unit-4 analysis microbenchmark: cost of 'use, intrinsic :: ieee_arithmetic' in UniformField_GetVel
! (IfW_FlowField.f90:385-410 @12fa3fc). Two copies of the body, identical except for the use statement.
module uf
   implicit none
   integer, parameter :: ReKi = 8
   type op_t
      real(ReKi) :: VelH, ShrV, ShrH, CosAngleH, SinAngleH, LinShrV, VelGust, CosAngleV, SinAngleV, VelV
   end type
contains
pure function get_ieee(RefHeight, RefLength, op, Position) result(Velocity)
   use, intrinsic :: ieee_arithmetic
   real(ReKi), intent(in) :: RefHeight, RefLength, Position(3)
   type(op_t), intent(in) :: op
   real(ReKi) :: Velocity(3), V1, V1_rotate, VZ_rotate
   V1 = op%VelH*((Position(3)/RefHeight)**op%ShrV &
                 + (op%ShrH*(Position(2)*op%CosAngleH + Position(1)*op%SinAngleH) &
                    + op%LinShrV*(Position(3) - RefHeight))/RefLength) + op%VelGust
   V1_rotate = op%CosAngleV*V1 - op%SinAngleV*op%VelV
   VZ_rotate = op%SinAngleV*V1 + op%CosAngleV*op%VelV
   Velocity = [V1_rotate*op%CosAngleH, -V1_rotate*op%SinAngleH, VZ_rotate]
end function
pure function get_plain(RefHeight, RefLength, op, Position) result(Velocity)
   real(ReKi), intent(in) :: RefHeight, RefLength, Position(3)
   type(op_t), intent(in) :: op
   real(ReKi) :: Velocity(3), V1, V1_rotate, VZ_rotate
   V1 = op%VelH*((Position(3)/RefHeight)**op%ShrV &
                 + (op%ShrH*(Position(2)*op%CosAngleH + Position(1)*op%SinAngleH) &
                    + op%LinShrV*(Position(3) - RefHeight))/RefLength) + op%VelGust
   V1_rotate = op%CosAngleV*V1 - op%SinAngleV*op%VelV
   VZ_rotate = op%SinAngleV*V1 + op%CosAngleV*op%VelV
   Velocity = [V1_rotate*op%CosAngleH, -V1_rotate*op%SinAngleH, VZ_rotate]
end function
end module
