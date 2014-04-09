subroutine orthonrm( a, n, k, lda )
!==========================================================
! Subroutine orthonormalizes the subspace b(lda,k) of a
!  matrix a(lda,n) 
!==========================================================
! Input:
!  a = matrix to be orthonormalized
!  n = dimension of a
!  k = dimension of subspace
!  lda = leading dimension of a
! Output:
!  v = output matrix
!----------------------------------------------------------
  implicit none
  integer, intent(in) :: n, k, lda
  real*8, dimension(lda,n) :: a

  integer :: i,j
  real*8, dimension(lda) :: vec1
  real*8, dimension(lda) :: vec2
  real*8 :: ddot
!----------------------------------------------------------
!
  do i=1, k
    a(:,i) = a(:,i) / ddot( lda, a(:,i), 1, a(:,i), 1 )
    do j=1, i-1
      a(:,j) = a(:,j) - ddot( lda, a(:,j), 1, a(:,i), 1)*a(:,i)
    end do
  end do
  return

end subroutine
