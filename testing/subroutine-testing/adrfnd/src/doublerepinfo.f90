! Subroutine to generate parity and index of new string after double
!  excitation
!====================================================================
!Input:
! string1       = original string
! str1len       = length of original string
! neworb1       = new orbital index1
! indxswp1      = index of old orbital index1
! neworb2       = new orbital index2
! indxswp2      = index of old orbital index2
! orbitals      = number of MO's in system
!Output:
! eps2         = parity
! indx         = index of new string
!==================================================================== 
subroutine doublerepinfo( string1, str1len, neworb1, indxswp1, neworb2, &
  indxswp2, orbitals, eps1, indx)

  use detci2, only: adrfnd

  implicit none
  
! ...input integer scalars...
  integer, intent(in) :: str1len, neworb1, indxswp1, orbitals, neworb2, &
                         indxswp2 

! ...input integer arrays...
  integer, dimension( str1len ) :: string1

! ...OUTPUT integer scalars...
  integer :: eps1, indx

! ...loop integer arrays...
  integer, dimension( str1len ) :: newstring

! ...loop integer scalars...
  integer :: i, j, k
!--------------------------------------------------------------------

! Swap orbitals
  newstring = string1
  newstring(indxswp1) = neworb1
  newstring(indxswp2) = neworb2

! Get eps1 & return new string to cannonical ordering
  call cannon( newstring, str1len, eps1 )

! Find new string index
  call adrfnd( newstring, str1len, orbitals, indx )

! Return
  return

end subroutine
