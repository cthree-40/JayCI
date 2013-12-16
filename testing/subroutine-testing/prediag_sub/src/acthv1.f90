subroutine acthv1( vector1, dgls, actstrlen, vector2 )
!==========================================================
! This subroutine computes the diagonal contribution to
!  Hv=u
!==========================================================
! Input:
!  vector1 = vector
!  dgls = diagonal matrix elements
!  actstrlen = length of vector1
! Output:
!  vector2 = diagonal contributions
!----------------------------------------------------------
  implicit none
  integer, intent(in) :: actstrlen
  real*8, dimension(actstrlen) :: vector1, vector2, dgls

  integer :: i
!----------------------------------------------------------
!

  do i=1, actstrlen
    vector2(i) = vector2(i) + dgls(i)*vector1(i)
  end do

  return

end subroutine
