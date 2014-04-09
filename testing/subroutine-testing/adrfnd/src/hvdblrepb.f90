real*8 function dblerepb( xorb1, xorb2, borb1, borb2, &
  max2e, moints2)
!==========================================================
! See dblerepa.f90
!==========================================================
  use detci5
!----------------------------------------------------------
  implicit none
  integer, intent(in) :: xorb1, xorb2, borb1, borb2, max2e

  real*8,dimension(max2e) :: moints2
!----------------------------------------------------------
!

  dblerepb = moints2(index2e2(borb1,xorb1,borb2,xorb2)) - &
                   moints2(index2e2(borb1,xorb2,borb2,xorb1))

end function dblerepb
