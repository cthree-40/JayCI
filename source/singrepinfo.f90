! Subroutine to generate parity and index of new string after single
!  excitation
!====================================================================
!Input:
! string1      = original string
! str1len      = length of original string
! neworb       = new orbital index
! indxwp       = index of old orbital index
! orbitals     = number of MO's in system
!Output:
! eps1         = parity
! indx         = index of new string
!==================================================================== 
subroutine singrepinfo( string1, str1len, neworb, indxswp, orbitals, &
  eps1, indx)

  use detci2, only: adrfnd

  implicit none
  
! ...input integer scalars...
  integer, intent(in) :: str1len, neworb, indxswp, orbitals

! ...input integer arrays...
  integer, dimension( str1len ), intent(in) :: string1

! ...OUTPUT integer scalars...
  integer, intent(out) :: eps1, indx

! ...loop integer arrays...
  integer, dimension( str1len ) :: newstring

! ...loop integer scalars...
  integer :: i, j, k
!--------------------------------------------------------------------

! Swap orbitas
  newstring = string1
  newstring(indxswp) = neworb

! Get eps1 & return new string to cannonical ordering
  call cannon( newstring, str1len, eps1 )

! Find new string index
! This if statement is for stringdiffs() in construct.f90
  if ( orbitals .eq. 0 ) then
          return
  end if
  call adrfnd( newstring, str1len, orbitals, indx )

! Return
  return 

end subroutine
