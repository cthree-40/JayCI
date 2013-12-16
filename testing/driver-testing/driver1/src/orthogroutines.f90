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
    real*8, dimension(lda, arrdim) :: matrix

!   ...loop integer scalars...
    integer :: i, j, k

!   ...loop real*8 scalars...
    real*8 :: overlap, ddot, norm

!--------------------------------------------------------------------

! Loop over vectors
    do i=1, matdim
! Loop over orthogonalizing vectors
      do j=1, i-1
        overlap = ddot( lda, matrix(1,i), 1, matrix(1,j), 1 )
        do k=1, lda
          matrix(k,i) = matrix(k,i) - overlap*matrix(k,j)
        end do
      end do
    end do
! Normalize
    do i=1, matdim
      norm = sqrt( ddot( lda, matrix(1,i), 1, matrix(1,i), 1))
      do j=1, lda
        matrix(j,i) = matrix(j,i) / norm
      end do
    end do

#ifdef GSDEBUG
    print *, " Debugging the Gram-Schmidt procedure.  You chose the entire reorthogonalization."
    print *, " Here are some overlaps:  "
    do i=1, matdim
      do j=i, matdim
         print *, ddot( lda, matrix(1,i), 1, matrix(1,j), 1)
      end do
    end do
    print *, " Debugging the Gram-Schmidt has finished. I hope you like what you saw. "
#endif
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

#ifdef GSDEBUG
    print *, " Debugging the Gram-Schmidt procedure. You chose the vector orthonormalization."
    print *, " Here are some overlaps: "
    do i=1, matdim
      print *, ddot( lda, matrix(1,i), 1, vector, 1)
    end do
    print *, ddot( lda, vector, 1, vector, 1 )
    print *, " Debugging the Gram-Schmidt has finished. I hope you like what you saw. "
#endif
    return
  end subroutine
end module
