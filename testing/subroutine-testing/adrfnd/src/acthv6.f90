subroutine acthv6( vector1, moints1, moints2, max1e, max2e,&
  actstr, actstrlen, aelec, belec, totels, orbitals, adets,&
  bdets, nfrzn, ndocc, ncas, vector2 )
!==========================================================
! Subroutine to compute the contribution of single 
!  replacements in alpha and beta strings to Hv=u
!==========================================================
	use detci2
!----------------------------------------------------------
! Input:
!  vector1 = input vector
!  moints1 = 1-e integrals
!  moints2 = 2-e integrals
!  max1e = length of moints1
!  max2e = length of moints2
!  actstr = list of active determinants
!  actstrlen = length of actstr(:)
!  aelec = alpha electrons
!  belec = beta electrons
!  totels = alpha + beta electrons
!  orbitals = number of orbitals
!  adets = alpha determinants( non-truncated )
!  bdets = beta determinants( non-truncated )
!  nfrzn = number of frozen orbitals
!  ndocc = number of docc orbitals
!  ncas = number of cas orbitals
! Output:
!  vector2 
!----------------------------------------------------------
	implicit none
	integer, intent(in) :: max1e, max2e, actstrlen, aelec, &
    belec, totels, orbitals, adets, bdets, nfrzn, ndocc, ncas
  real*8, dimension(actstrlen) :: vector1, vector2
  integer, dimension(actstrlen) :: actstr
  real*8, dimension(max1e) :: moints1
  real*8, dimension(max2e) :: moints2

  integer, dimension(:,:), allocatable :: alphmat, betamat
  integer, dimension(aelec) :: pstr1, npstr
  integer, dimension(belec) :: qstr1, nqstr
  integer, dimension(orbitals - aelec) :: pexits
  integer, dimension(orbitals - belec) :: qexits

  real*8 :: singrepab

  integer :: p, q, testp, testq, testk, i, j, k, l, &
    m, n, dindx, eps1, eps2, v
!----------------------------------------------------------
!
  allocate( alphmat(adets, aelec))
  allocate( betamat(bdets, belec))

  call strfnd(adets, aelec, orbitals, adets, alphmat)
  call strfnd(bdets, belec, orbitals, bdets, betamat)

  do i=1, actstrlen
    call k2indc( actstr(i), belec, orbitals, p, q )

    do j=1, aelec
      pstr1(j) = alphmat(p, j)
    end do

    call possex1( pstr1, orbitals, aelec, orbitals-aelec, pexits )
    do j=1, belec
      qstr1(j) = betamat(q, j)
    end do

    call possex1( qstr1, orbitals, belec, orbitals-belec, qexits )
    do j=1, aelec
      do k=1, belec
        do l=1, orbitals-aelec
          do m=1, orbitals-belec
            do n=1, aelec
              npstr(n) = pstr1(n)
            end do
            do n=1, belec
              nqstr(n) = qstr1(n)
            end do
            npstr(j) = pexits(l)
            nqstr(k) = qexits(m)

            call cannon( npstr, aelec, eps1 )
            call cannon( nqstr, belec, eps2 )
            call adrfnd( nqstr, belec, orbitals, testq )
            call adrfnd( npstr, aelec, orbitals, testp )
            testk = 0
            dindx = indxk( testp, testq, totels, orbitals )
            do n=1, actstrlen
              if ( dindx == actstr(n) ) then
                testk=1
                v=n
								exit
              else
                cycle
              end if
            end do
            if ( testk == 0 ) then
              cycle
            else
              vector2(i) = vector2(i) + eps1 * eps2 * &
                singrepab( pexits(l), qexits(m), pstr1(j),&
                qstr1(k), moints2, max2e ) * vector1(v)
            end if
          end do
        end do
      end do
    end do
  end do

  return

end subroutine
            
