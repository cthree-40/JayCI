subroutine acthv2( vector1, moints1, moints2, max1e, &
  max2e, actstr, actstrlen, aelec, belec, totels, &
  orbitals, adets, bdets, nfrzn, ndocc, ncas, vector2 )
!==========================================================
! This subroutine computes the contribution of single
!  excitations in alpha strings to Hv=u
!==========================================================
  use detci2
!----------------------------------------------------------
! Input:
!  vector1 = initial vector
!  moints1 = one electron integrals
!  moints2 = two electron integrals
!  max1e = length of moints1(:)
!  max2e = length of moints2(:)
!  actstr = string of active determinants
!  actstrlen = length of actstr
!  aelec = number of alpha electrons
!  belec = number of beta electrons
!  totels = aelec + belec
!  orbitals = number of orbitals
!  adets = alpha determinants( non truncated )
!  bdets = beta determinants( non truncated )
!  nfrzn = number of frozen orbitals
!  ndocc = number of doubly occupied orbitals
!  ncas = number of cas orbitals
! Output:
!  vector2
!----------------------------------------------------------
  implicit none
  integer, intent(in) :: max1e, max2e, actstrlen, aelec, &
    belec, totels, orbitals, adets, bdets, nfrzn, ndocc, ncas
  real*8, dimension(max1e) :: moints1
  real*8, dimension(max2e) :: moints2
  real*8, dimension(actstrlen) :: vector2, vector1
  integer, dimension(actstrlen) :: actstr
  
  integer, dimension(:,:), allocatable :: alphmat, betamat
  integer :: p, q, i, j, k, testp, testk, dindx, l, eps1, v
  integer, dimension(aelec) :: pstr1, npstr
  integer, dimension(orbitals-aelec) :: pexits
  
  real*8 :: singrepa

  integer :: debug1, debug2
!----------------------------------------------------------
!

! Allocate and construct the alpha and beta string matrices
  allocate(alphmat(adets, aelec))
  allocate(betamat(bdets, belec))

  call strfnd(adets, aelec, orbitals, adets, alphmat)
  call strfnd(bdets, belec, orbitals, bdets, betamat)
! 
  do i=1, actstrlen
    call k2indc(actstr(i), belec, orbitals, p, q)

    do j=1, aelec
      pstr1(j) = alphmat(p,j)
    end do

    call possex1( pstr1, orbitals, aelec, orbitals-aelec, pexits)

! Loop through the electrons in pstr1()
	  do j=1, aelec

! Loop though possible excitations
		  do k=1, orbitals-aelec
		    npstr = pstr1
		    npstr(j) = pexits(k)
		    call cannon( npstr, aelec, eps1)

! Check if < (pstr1), q | is a valid determinant in expansion
				call adrfnd(npstr, aelec, orbitals, testp)
				testk = 0
				dindx = indxk( testp, q, totels, orbitals )

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
#if ACTION=1
          print *, " Single excitation contribution...determinant ", actstr(v)
          print *, " The alpha string in this determinant looks like:    "
          call k2indc( actstr(v), belec, orbitals, debug1, debug2 )
          print *, alphmat( debug1, :) 
#endif
				  vector2(i) = vector2(i) + eps1 * singrepa( p, q, &
            totels, aelec, belec, orbitals, alphmat, betamat, &
            adets, bdets, max1e, max2e, moints1, moints2, &
            pexits(k), pstr1(j), actstr, actstrlen ) * vector1(v)
				  testk=0
				end if
		  end do
		end do
	end do
  return
end subroutine
    
