! Subroutine to compute residual
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chemistry, The Johns Hopkins University
!
! Last edit: 10-23-2013 - rewritten from Det-CI-4.0.0
!====================================================================
!Input:
! hvectors   = action of H on the basis vectors  real*8  array   2-d
! basisdim   = current dimension of space        integer scalar
! ldabasis   = LDA of hvectors and basisvecs     integer scalar
! basisvecs  = basis vectors                     real*8  array   2-d
! curroot    = current root                      integer scalar
! kreigvals  = eigenvalues of krylov space       real*8  array   1-d
! kreigvecs  = eigenvectors of krylov space      real*8  array   2-d
!Output:
! residual   = ( H - p)x                         real*8  array   1-d
!====================================================================
subroutine residgen( hvectors, basisdim, ldabasis, basisvecs, curroot,&
  kreigvals, kreigvecs, residual )

  implicit none

! ...input integer scalars...
  integer, intent(in) :: basisdim, ldabasis, curroot

! ...input real*8 arrays...
  real*8, dimension(ldabasis, basisdim) :: hvectors, basisvecs
  real*8, dimension(basisdim, basisdim) :: kreigvecs
  real*8, dimension(basisdim) :: kreigvals

! ...output real*8 array...
  real*8, dimension(ldabasis) :: residual

! ...loop integer scalars...
  integer :: i, j, k

!--------------------------------------------------------------------

! Zero out residual vector
  residual = 0d0

! Loop over dimension of space
  do i=1, basisdim
    do j=1, ldabasis
      residual(j) = residual(j) + kreigvecs(i, curroot)*( hvectors(j,i) &
                    - kreigvals(curroot)*basisvecs(j,i))
    end do
  end do

  return
 

end subroutine
