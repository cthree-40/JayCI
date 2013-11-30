subroutine possex1( string, orbitals, electrons, pexitslen, &
  pexits )
!==========================================================
! This subroutine generates the list of possible excitations
!  of string(:)
!==========================================================
!----------------------------
! Input:
!  string = determinant string (alpha or beta)
!  orbitals = number of orbitals
!  electrons = number of electrons (alpha or beta)
!  pexitslen = length of pexits(:)
! Output:
!  pexits(:)
!----------------------------
  implicit none
  integer,intent(in) :: orbitals, electrons, pexitslen
  integer,dimension(electrons) :: string
  integer,dimension(pexitslen) :: pexits
  integer,dimension(orbitals) :: totorbs

  integer :: i, j, point
!----------------------------
!
! Write out string of orbitals
  do i=1, orbitals
    totorbs(i) = i
  end do

  point=1

! Loop through string() and totorbs(). If totorbs(i) = string(j),
!  make totorbs(i) = 0.
  do i=1, electrons
    do j = point, orbitals
      if ( string(i) == totorbs(j) ) then
        totorbs(j) = 0
        point = 1
        exit
      else
        cycle
      end if
    end do
  end do

! Construct pexits
  do i=1, pexitslen
    do j=1, orbitals
      if ( totorbs(j) .ne. 0 ) then
        pexits(i) = totorbs(j)
        totorbs(j) = 0
        exit
      else
        cycle
      end if
    end do
  end do

  return

end subroutine
