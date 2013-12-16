! Subroutine to select lowest energy E_i = < K_i | H | K_i >
!====================================================================
!Input:
! vector1   = < K | H | K >               real*8  array  1-d
! length    = legnth of diagonals()       integer scalar
! choose    = number of diagonals to pick integer scalar
!Output:
! initguess = singular vectors            real*8  array  2-d
!====================================================================
subroutine lowdiags( vector1, length, choose, initguess )
 
  implicit none

! ...input real*8 arrays...
  real*8, dimension(length) :: vector1

! ...input integer scalars...
  integer, intent(in) :: length, choose

! ...output real*8 arrays...
  real*8, dimension(length, choose) :: initguess

! ...loop integer scalars...
  integer :: i, j, k, l

! ...loop real*8 arrays...
  real*8, dimension(length) :: vector2
!--------------------------------------------------------------------

! Order the < K | H | K >
  vector2 = vector1
  call cannonreal( vector2, length, l )

  print *, "Lowest diagonals are..."
  do i=1, choose
    print *, vector2(i)
  end do
  initguess = 0d0
! Find lowest K's
  LPA: do i=1, choose
    LPB: do j=1, length
      if ( vector2( i ) .eq. vector1(j) ) then
        initguess(j, i) = 1d0
        cycle LPA
      else
        cycle LPB
      end if
    end do LPB
  end do LPA

  return

end subroutine 
