subroutine gramschmidt( matrix, vector1, length, dimen, &
  vectorf )
!==========================================================
! This subroutine executes the gram schmidt algorithm
!==========================================================
! Input:
!  matrix = collection of orthonormal vectors
!  vector1 = vector to be orthogonalized
!  length = length of vector1, leading dim. of matrix
!  dimen = dimension of matrix
! Output:
!  vectorf = orthogonalized vector
!----------------------------------------------------------
  implicit none
  integer,intent(in) :: length, dimen
  real*8,dimension(length,dimen) :: matrix
  real*8,dimension(length) :: vector1

  real*8,dimension(length) :: vector2, vector3
  real*8,dimension(length,dimen) :: hldvecs

  real*8,dimension(length) :: vector4, vectorf

  integer :: i, j
  real*8 :: ddot
!----------------------------------------------------------
  
  do i=1, dimen
    do j=1, length
      vector2(j) = matrix(j,i)
    end do
    do j=1, length
      hldvecs(j,i) = ddot(length, vector1, 1, vector2, 1)*vector2(j)
    end do
  end do
  
  vectorf = vector1
  do i=1, dimen
    vectorf = vectorf - hldvecs(:,i)
  end do

  vectorf = vectorf / ddot( length, vectorf, 1, vectorf, 1)
  return
end subroutine
