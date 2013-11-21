! Module containing prediagonalization routines
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chem., The Johns Hopkins Unversity
!
! Last edit: 11-20-13
!====================================================================
module prediag

contains
!====================================================================
!====================================================================
!> geninituvec
!
! Subroutine to generate initial guess unit vectors
!--------------------------------------------------------------------
  subroutine geninituvec( numvec, lenvec, initvecs )
! Input:
!  numvec = number of desired vectors
!  lenvec = length of desired vectors
! Output:
!  initvecs = numvec unit vectors
!--------------------------------------------------------------------
    implicit none
! ...input integer scalars...
    integer, intent(in) :: numvec, lenvec
! ...output real*8 arrays...
    real*8, dimension(lenvec,numvec), intent(out) :: initvecs
! ...loop integer scalars...
    integer :: i,j
!--------------------------------------------------------------------
    initvecs = 0d0
    do i=1, numvec
      initvecs(i,i) = 1d0
    end do
    return
  end subroutine geninituvec
!====================================================================
!====================================================================
!> lowdiagprecond
!
! Subroutine to generate intial guess unit vectors by grabbing the 
!  lowest energy diagonals, and diagonalizing the v.Hv matrix
!--------------------------------------------------------------------
  subroutine lowdiagprecond( diagonals, cidim,
      
end module
