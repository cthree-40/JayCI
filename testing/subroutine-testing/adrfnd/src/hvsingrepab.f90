real*8 function singrepab( xorba, xorbb, borba, borbb, &
  moints2, max2e)
!==========================================================
! This function computes the value of a single excitation
!  in the alpha and beta string
!==========================================================
  use detci5
!----------------------------------------------------------
! Inputs:
!  xorba = excitation of alpha string
!  xorbb = excitation of beta string
!  borba = excitation from orbital borba in alpha string
!  borbb = excitation from orbital borbb in beta string
!  moints1 = 1-e integrals
!  moints2 = 2-e integrals
!  max1e   = dimension of moints1
!  max2e   = dimension of moints2
!----------------------------------------------------------
  implicit none
  integer,intent(in) :: xorba, xorbb, borba, borbb, max2e

  real*8,dimension(max2e) :: moints2
!----------------------------------------------------------
!

  singrepab = moints2(index2e2(borba,xorba,borbb,xorbb))

end function singrepab
