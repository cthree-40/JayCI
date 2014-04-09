real*8 function dblerepa( xorb1, xorb2, borb1, borb2, &
  max2e, moints2 )
!==========================================================
! This function computes the matrix element between 
!  determinants differing by two excitations in the alpha
!  string.
!==========================================================
  use detci5
!----------------------------------------------------------
! Inputs:
! xorb1   = electron 1 is excited into this orbital
! xorb2   = electron 2 is excited into this orbital
! borb1   = electron 1 is excited from this orbital
! borb2   = electron 2 is excited from this orbital
! moints2 = 2-e integrals
!----------------------------------------------------------
  implicit none
  integer,intent(in) :: borb1, borb2, xorb1, xorb2, max2e

  real*8, dimension(max2e) :: moints2
!----------------------------------------------------------
!

  dblerepa = moints2(index2e2(borb1,xorb1,borb2,xorb2)) - &
                       moints2(index2e2(borb1,xorb2,borb2,xorb1))

end function dblerepa
