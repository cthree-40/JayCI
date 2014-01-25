subroutine acthv5( vector1, moints1, moints2, max1e, max2e,&
  actstr, actstrlen, aelec, belec, totels, orbitals, adets,&
  bdets, nfrzn, ndocc, ncas, vector2 )
!==========================================================
! Double excitations in beta strings
!==========================================================
	use detci2
!----------------------------------------------------------
! Input: 
!  vector1 = input vector
!  moints1 = 1-e integrals
!  moints2 = 2-e integrals
!  max1e = length of moints1(:)
!  max2e = length of moints2(:)
!  actstr = active determinants
!  actstrlen = length of active determinants
!  aelec = alpha electrons
!  belec = beta electrons
!  totels = alpha + beta electrons
!  orbitals = number of orbitals
!  adets = alpha determinants( non truncated )
!  bdets = beta determinants( non truncated )
!  nfrzn = number of frozen orbitals
!  ndocc = number of docc orbitals
!  ncas = number of cas orbitals
! Output:
!  vector2
!----------------------------------------------------------
	implicit none
  integer, intent(in) :: max1e, max2e, actstrlen, aelec, &
    belec, totels, orbitals, adets, bdets, nfrzn, ndocc, ncas
  integer, dimension(actstrlen) :: actstr
  real*8, dimension(actstrlen) :: vector1, vector2
  real*8, dimension(max1e) :: moints1
  real*8, dimension(max2e) :: moints2

  integer, dimension(:,:), allocatable :: alphmat, betamat
  integer, dimension(belec) :: qstr1, nqstr
  integer, dimension(orbitals-belec) :: qexits
  integer :: p, q, i, j, k, l, m, n, eps1, testk, testq, dindx

  real*8 :: dblerepb
!----------------------------------------------------------
!

! Allocate and construct the alpha and beta string matrices
  allocate(alphmat(adets, aelec))
  allocate(betamat(bdets, belec))

  call strfnd(adets, aelec, orbitals, adets, alphmat)
  call strfnd(bdets, belec, orbitals, bdets, betamat)

! Loop through active determinants
  do i=1, actstrlen

! Obtain deterimant alpha and beta string indices
    call k2indc( actstr(i), belec, orbitals, p, q ) 

! Set pstr1
    do j=1, belec
      qstr1(j) = betamat(q,j)
    end do

! Obtain set of possible excitations
    call possex1( qstr1, orbitals, belec, orbitals-belec, qexits )

! Loop over electrons j
    do j=1, belec

! Loop over electrons k > j
      do k=j+1, belec

! Loop over excitations l
        do l=1, orbitals-belec

! Loop over excitations m > l
          do m=l+1, orbitals-belec
            nqstr = qstr1
            nqstr(j) = qexits(l)
            nqstr(k) = qexits(m)
            call cannon( nqstr, belec, eps1 )

! Check if < (npstr), q | is in expansion
						call adrfnd( nqstr, belec, orbitals, testq )
						testk = 0
						dindx = indxk( p, testq, totels, orbitals )
						do n=1, actstrlen
						  if ( dindx == actstr(n) ) then
						    testk=1
						    exit
! A determinant match has been found, so leave loop and 
!  compute contribution
						  else
						    cycle

! Check next determinant
						  end if
						end do
						if ( testk == 0 ) then
						  cycle
						else
						  vector2(i) = vector2(i) + eps1*dblerepb( &
                qexits(l), qexits(m), qstr1(j), qstr1(k), &
                max2e, moints2 )*vector1(n)
						  testk=0
						end if
				  end do
				end do
		  end do
		end do
  end do

  return

end subroutine 
