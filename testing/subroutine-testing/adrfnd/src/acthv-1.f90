subroutine acthv(vector1, moints1, moints2, max1e, max2e, &
  actstr, actstrlen, aelec, belec, totels, orbitals, &
  adets, bdets, dgls, nfrzn, ndocc, ncas, vector2)
!==========================================================
	use detci2
!----------------------------------------------------------
!Input:
! vector1  = input vector
! moints1  = 1-e integrals
! moints2  = 2-e integrals
! max1e    = length of moints1()
! max2e    = length of moints2()
! actstr   = list of active determinants
! actstrlen= length of actstr()
! aelec    = alpha electrons
! belec    = beta electrons
! totels   = total electrons
! orbitals = number of MOs
! adets    = alpha determinants
! bdets    = beta determinants
! dgls     = < K | H | K >
!Output:
! vector2  = output vector; Hv = vector2
!
! vector2(K) = SUM_L < K | H | L > vector2(L)
!----------------------------------------------------------
	implicit none
	integer,intent(in) :: max1e, max2e, actstrlen, aelec, &
    belec, totels, orbitals, adets, bdets, nfrzn, ndocc, ncas
	real*8,dimension(max1e) :: moints1
	real*8,dimension(max2e) :: moints2
	integer,dimension(actstrlen) :: actstr
	real*8,dimension(actstrlen) :: vector1,vector2,dgls
	
  integer :: i
!----------------------------------------------------------
!

! Diagonal contribution
	call acthv1( vector1, dgls, actstrlen, vector2 )

! Single replacements in alpha strings
	call acthv2( vector1, moints1, moints2, max1e, max2e, &
    actstr, actstrlen, aelec, belec, totels, orbitals, &
    adets, bdets, nfrzn, ndocc, ncas, vector2)
  
! Single replacements in beta strings
	call acthv3( vector1, moints1, moints2, max1e, max2e, &
    actstr, actstrlen, aelec, belec, totels, orbitals, &
    adets, bdets, nfrzn, ndocc, ncas, vector2)

! Double replacements in alpha strings
	call acthv4( vector1, moints1, moints2, max1e, max2e, &
    actstr, actstrlen, aelec, belec, totels, orbitals, &
    adets, bdets, nfrzn, ndocc, ncas, vector2)

! Double replacements in beta strings
	call acthv5( vector1, moints1, moints2, max1e, max2e, &
    actstr, actstrlen, aelec, belec, totels, orbitals, &
    adets, bdets, nfrzn, ndocc, ncas, vector2 )

! Single replacements in alpha and beta strings
 	call acthv6( vector1, moints1, moints2, max1e, max2e, &
    actstr, actstrlen, aelec, belec, totels, orbitals, &
    adets, bdets, nfrzn, ndocc, ncas, vector2 )

  return

end subroutine acthv

