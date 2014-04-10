subroutine spacetrunc( space, length, krmin, krmax, root, &
  eigvecs )
!==========================================================
! Subroutine truncates the kyrlov space
!==========================================================
! Input:
!  space = basis space
!  length = dimension of vectors in space(length, krmax)
!  krmin = minimum dimension of Krylov space
!  krmax = maximum dimension of Krylov space
!  root = root currently solving for
!----------------------------------------------------------
  implicit none
  integer, intent(in) :: length, krmin, krmax, root
  real*8, dimension(length, krmax) :: space  
  real*8, dimension(krmax,krmax) :: eigvecs
   
  integer :: i, j, l
  real*8 :: ddot

  character*1 ::  transa, transb
  integer :: m, n, k
  real*8 :: alpha
  real*8,dimension(length,krmax) ::  c
  integer :: lda, ldb, ldc
  real*8 :: beta

  real*8,dimension(length) :: vec1, vec2, vec3, vec4
  real*8 :: ovrlp, norm1

!----------------------------------------------------------
!
! Given space(length,krmax), truncate to space(length,krmin)
! B_ij e_jk = A_ik
!

  transa = 'n'
  transb = 'n'
  m = length
  n = krmax
  k = krmax
  alpha = 1d0
  lda = m
  ldb = k
  beta = 0d0
  ldc = m

  call dgemm( transa, transb, m, n, k, alpha, space, lda, eigvecs, &
    ldb, beta, c, ldc ) 

! Zero out space
  space = 0d0


! make c = space
  do i=1, krmin
    do j=1, length
      space(j,i) = c(j,i)
    end do
  end do
  
  return

end subroutine
