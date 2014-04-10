! Module containing various orthogonalization routines for Davidson Alg.
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chemistry, The Johns Hopkins University
!
! Last edit: 10-23-2013
!====================================================================
module orthogroutines
  implicit none
contains
!--------------------------------------------------------------------
! Orthogonalizes real*8 matrix
!====================================================================
! Input:
!  matrix   = matrix to be orthogonalized                 real*8  array  2-d
!  lda      = leading dimension of matrix                 integer scalar
!  matdim   = dimension of matrix                         integer scalar
!  arrdim   = dimension of array passed ( >matdim )       integer scalar
! Ouput:
!  matrix   = orthogonalized
!====================================================================
  subroutine modgramschmidt( matrix, matdim, arrdim, lda )

    implicit none

!   ...input integer scalars...
    integer, intent(in) :: matdim, arrdim, lda

!   ...input/output real*8 array...
    real*8, dimension(lda, matdim), intent(inout) :: matrix

!   ...loop integer scalars...
    integer :: i, j, k

!   ...loop real*8 scalars...
    real*8 :: overlap, ddot, norm
!   ...loop real*8 arrays...
    real*8, dimension( lda, matdim ) :: scratch
    real*8, dimension( lda ) :: vec_scratch
!--------------------------------------------------------------------

!  Loop over vectors
    do i=1, matdim
      norm =  sqrt( ddot( lda, matrix(1,i), 1, matrix(1,i), 1 ) )
      do j=1, lda
        vec_scratch(j) = matrix(j,i) / norm
      end do
      do j= 1, i-1
        call orthogvector( matrix, j, matdim, lda, vec_scratch )
      end do
      do j=1, lda
        matrix(j,i) = vec_scratch(j)
      end do
    end do


    return

  end subroutine


!--------------------------------------------------------------------
! Orthogonalizes a vector to a space
!==================================================================== 
!Input:
! ( See above)
! vector = vector to orthonormalize          real*8 array 1-d
!Output:
! vector = now orthogonormalize              real*8 array 1-d
!====================================================================
  subroutine orthogvector( matrix, matdim, arrdim, lda, vector )

    implicit none

!   ...input integer scalars...
    integer, intent(in) :: matdim, arrdim, lda

!   ...input real*8 arrays...
    real*8, dimension(lda, arrdim) :: matrix

!   ...input/output real*8 arrays...
    real*8, dimension(lda) :: vector

!   ...loop integer scalars...
    integer :: i, j, k

!   ...loop real*8 scalars...
    real*8 :: norm, overlap, ddot

!--------------------------------------------------------------------
 
! Loop over vectors
    do i=1, matdim
      overlap = ddot( lda, matrix(1,i), 1, vector, 1 )
      do j=1, lda
        vector(j) = vector(j) - overlap*matrix(j,i)
      end do
    end do

! Normalize
    norm = sqrt( ddot( lda, vector, 1, vector, 1 ) )
    do i=1, lda
      vector(i) = vector(i) / norm
    end do

    return
  end subroutine
end module
