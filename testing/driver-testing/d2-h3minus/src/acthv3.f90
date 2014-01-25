subroutine acthv3( vector1, moints1, moints2, max1e, max2e,&
  actstr, actstrlen, aelec, belec, totels, orbitals, adets,&
  bdets, nfrzn, ndocc, ncas, vector2 )
!==========================================================
! This subroutine computes the contribution of single 
!  excitations in beta strings to Hv=u
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
!  aelec = number of alpha electrons
!  belec = number of beta electrons
!  totels = aelec + belec
!  orbitals = number of orbitals
!  adets = alpha determinants( non-truncated )
!  bdets = beta determinants(  non-truncated )
!  nfrzn = number of frozen orbitals
!  ndocc = number of docc orbitals
!  ncas = number of cas orbitals
! Output:
!  vector2
!----------------------------------------------------------
  implicit none
  integer, intent(in) :: max1e, max2e, actstrlen, aelec, &
    belec, totels, orbitals, adets, bdets, nfrzn, ndocc, ncas
  real*8, dimension(actstrlen) :: vector2, vector1
  integer, dimension(actstrlen) :: actstr
  
  integer, dimension(:,:), allocatable :: alphmat, betamat
  integer :: p, q, i, j, k, testq, testk, dindx, l, eps1, v
  integer, dimension(belec) :: qstr1, nqstr
  integer, dimension(orbitals-belec) :: qexits
  real*8, dimension(max1e) :: moints1
  real*8, dimension(max2e) :: moints2
  
  real*8 :: singrepb
!----------------------------------------------------------
!

! Allocate and construct the alpha and beta string matrices
  allocate(alphmat(adets, aelec))
  allocate(betamat(bdets, belec))

  call strfnd(adets, aelec, orbitals, adets, alphmat)
  call strfnd(bdets, belec, orbitals, bdets, betamat)
 
  do i=1, actstrlen
    call k2indc(actstr(i), belec, orbitals, p, q)
    do j=1, belec
      qstr1(j) = betamat(q,j)
    end do
    call possex1( qstr1, orbitals, belec, orbitals-belec, qexits)

! Loop through the electrons in pstr1()
		do j=1, belec

! Loop though possible excitations
		  do k=1, orbitals-belec
		    nqstr = qstr1
		    nqstr(j) = qexits(k)
		    call cannon( nqstr, belec, eps1)

! Check if < (pstr1), q | is a valid determinant in expansion
				call adrfnd(nqstr, belec, orbitals, testq)
				testk = 0
				dindx = indxk( p, testq, totels, orbitals )
				do l=1, actstrlen
				  if (dindx == actstr(l)) then
				    testk=1
						v=l
				    exit
				  else
				    cycle
				  end if
				end do
				if ( testk == 0 ) then
				  cycle
				else
				  vector2(i) = vector2(i) + eps1 * singrepb( p, q, &
            totels, aelec, belec, orbitals, alphmat, &
            betamat, adets, bdets, max1e, max2e, moints1, &
            moints2, qexits(k), qstr1(j), actstr, &
            actstrlen ) * vector1(v)
				  testk=0
				end if
		  end do
		end do
	end do
  return
end subroutine
