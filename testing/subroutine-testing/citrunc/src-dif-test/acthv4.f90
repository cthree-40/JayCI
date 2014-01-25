subroutine acthv4( vector1, moints1, moints2, max1e, &
  max2e, actstr, actstrlen, aelec, belec, totels, &
  orbitals, adets, bdets, nfrzn, ndocc, ncas, vector2 )
!==========================================================
! Subroutine to compute the contribution of double 
!  excitations in alpha strings
!==========================================================
  use detci2
!----------------------------------------------------------
! Input:
!  vector1 = input vector
!  moints1 = 1-e integrals
!  moints2 = 2-e integrals
!  max1e = length of 1-e integrals
!  max2e = length of 2-e integrals
!  actstr = string of active determinants
!  actstrlen = length of actstr(:)
!  aelec = alpha electrons
!  belec = beta electrons
!  totels = alpha + beta electrons
!  orbitals = number of orbitals
!  adets = alpha determinants( non-truncated )
!  bdets = beta determinants( non-truncated )
!  nfrzn = number of frozen electrons
!  ndocc = number of docc electrons
!  ncas  = number of cas electrons
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
  integer, dimension(aelec) :: pstr1, npstr
  integer, dimension(orbitals-aelec) :: pexits
  integer :: p, q, i, j, k, l, m, n, eps1, testk, testp, dindx, v

  real*8 :: dblerepa
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
    do j=1, aelec
      pstr1(j) = alphmat(p,j)
    end do

! Obtain set of possible excitations
    call possex1( pstr1, orbitals, aelec, orbitals-aelec, pexits )

! Loop over electrons j
    do j=1, aelec

! Loop over electrons k > j
      do k=j+1, aelec

! Loop over excitations l
        do l=1, orbitals-aelec

! Loop over excitations m > l
          do m=l+1, orbitals-aelec
            npstr = pstr1
            npstr(j) = pexits(l)
            npstr(k) = pexits(m)
            call cannon( npstr, aelec, eps1 )

! Check if < (npstr), q | is in expansion
						call adrfnd( npstr, aelec, orbitals, testp )
						testk = 0
						dindx = indxk( testp, q, totels, orbitals )
						do n=1, actstrlen
						  if ( dindx == actstr(n) ) then
						    testk=1
								v=n
						    exit

! A determinant match has been found, so leave loop and compute contribution
						  else
						    cycle

! Check next determinant
						  end if
						end do
						if ( testk == 0 ) then
						  cycle
						else
						  vector2(i) = vector2(i) + eps1*dblerepa( &
                pexits(l), pexits(m), pstr1(j), pstr1(k), &
                max2e, moints2 )*vector1(v)
						  testk=0
						end if
				  end do
				end do
		  end do
		end do
  end do
  return
end subroutine
  
  
