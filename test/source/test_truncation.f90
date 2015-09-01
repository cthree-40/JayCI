program test_truncation
  ! test_truncation
  ! ---------------
  ! Purpose: Test validity of wavefunction.
  !
  !-----------------------------------------------------------------------------
  implicit none

  ! .. LOCAL scalars ..
  integer :: ndets, orbitals
  integer :: ios, i, j, k
  integer :: odocc, oactv, ndocc, nactv
  
  ! .. LOCAL arrays ..
  integer, dimension(:,:), allocatable :: detocc

  print *, "==================================================================="
  
  print *, "Enter (orbitals - nfrozen): "
  read  *, orbitals
  print *, "Enter DOCC: "
  read  *, ndocc
  print *, "Enter ACTIVE: "
  read  *, nactv
  print *, " Orbitals = ", orbitals
  print *, " ndocc    = ", ndocc
  print *, " nactv    = ", nactv
  print *, ""
  print *, "Reading det.list"
  ndets = 0
  open(file="det.list",unit=9,action="read",status="old")
  do
          read(9,10,iostat=ios)
          
          if (ios .eq. 0) then
                  ndets = ndets + 1
          else if (ios .lt. 0) then
                  exit
          else
                  stop "***ERROR***"
          end if

  end do
    
  allocate(detocc(orbitals,ndets))
  rewind(9)
  do i = 1, ndets
          read(9,10) detocc(1:orbitals,i)
          write(*,"(23(i3))") detocc(1:orbitals,i)
  end do

  ! loop through occupation strings noting occupations of each section
  write(*,"(80('-'))")
  do i = 1, ndets

          odocc = 0
          do j = 1, ndocc
                  if (detocc(j,i) .gt. 0) then
                          if (detocc(j,i) .lt. 3) then
                                  odocc = odocc + 1
                          else
                                  odocc = odocc + 2
                          end if
                  end if
          end do

          oactv = 0
          do j = ndocc + 1, nactv
                  if (detocc(j,i) .gt. 0) then
                          if (detocc(j,i) .lt. 3) then
                                  oactv = oactv + 1
                          else
                                  oactv = oactv + 2
                          end if
                  end if
          end do

          write(*,11) i, odocc, oactv
          
  end do
10 format(57x,200i4)
11 format(1x,'Determinant = ',i15,5x, 'DOCC occ =',i5, 5x &
            'ACTV occ =',i5)
end program test_truncation
