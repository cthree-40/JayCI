subroutine orthog(vector1, vecdim, spdim,krmax, space, vector5)
!==========================================================
! Subroutine to orthonormalize a real vector a space of 
!  real vectors
! V' = V - M.M^T.V
! V'' = V'/ || V' ||
!==========================================================
! Input:
!  vecdim = dimension of vector
!  vector1 = unorthogonalized vector
!  spdim = dimension of space
!  space = orthonormal collection of vectors
! Output:
!  vector5 = orthonomralized vector (output)
!----------------------------------------------------------
  implicit none
  integer :: vecdim, spdim, krmax
  real*8,dimension(vecdim) :: vector1
  real*8,dimension(vecdim,krmax) :: space

  real*8,dimension(vecdim) :: vector5

  real*8,parameter :: beta=0d0, alpha=1d0
  real*8,dimension(spdim) :: vector2
  real*8,dimension(vecdim) :: vector3

  integer :: i

  real*8 :: ddot
!----------------------------------------------------------
!

! Perform M^T.V
  call dgemv( 'T', vecdim, spdim, alpha, space, vecdim, &
    vector1, 1, beta, vector2, 1 )

! Perform M(M^T.V)
  call dgemv( 'N', vecdim, spdim-1, alpha, space, vecdim, &
    vector2, 1, beta, vector3, 1 )

! Subtract vector3 from vector1
  do i=1, vecdim
    vector1(i) = vector1(i) - vector3(i)
  end do

! Normalize vector1, call this vector vector5
  do i=1, vecdim
    vector5(i) = vector1(i) / ddot( vecdim, vector1, 1, &
      vector1, 1 )
  end do

! Test Orthogonality
#if ORTHOG_TEST > 0
  do i=1, spdim
    print *, ddot( vecdim, space(:,i), 1, vector5, 1 )
  end do
#endif

  return
end subroutine
  
