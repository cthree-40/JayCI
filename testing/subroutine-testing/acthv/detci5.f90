module detci5
!==========================================================
! This module contains routines for locating the 1 and 2
! electron integrals
!
! Subroutines:
!  ind2val
!  index2e2
!----------------------------------------------------------
! Christopher L. Malbon
! Yarkony Group
! Dept. Chemistry, The Johns Hopkins University
!
! Last edit: 
! 09-27-13 - removed SUMINDS and INDEX2E
!==========================================================
  implicit none

contains
!----------------------------------------------------------
! ind2val(i,j) 
!
! computes index of element (i, j) packed in lower triangle
!  array
!----------------------------------------------------------
  integer function ind2val(i, j)
!----------------------------
! Input:
!  i = row number
!  j = column number
!----------------------------
    implicit none
    integer, intent(in) :: i, j
!----------------------------
!
    if ( i .ge. j ) then
      ind2val = (i - 1)*i/2 + j
    else
      ind2val = (j - 1)*j/2 + i
    end if

  end function ind2val

!----------------------------------------------------------
! index2e2( i, j, k, l)
!
! computes index of element (i, j, k, l) in twice packed
!  lower triangle array
!----------------------------------------------------------
  integer function index2e2(i, j, k, l)
!----------------------------
    implicit none
    integer, intent(in) :: i, j, k, l
!----------------------------
!
    index2e2 = ind2val(ind2val(i,j),ind2val(k, l))

  end function

!----------------------------------------------------------
!==========================================================
end module
