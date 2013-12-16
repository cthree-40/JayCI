real*8 function singrepa( p1,q1, elecs, aelec, belec,  &
  orbtls, alphmat, betamat, adets, bdets, max1e, max2e,&
  moints1, moints2, xorb, borb, nactstr, live )
!==========================================================
! This function computes the single alpha replacement term
!
! < K | H | L > = < p,q | H | r(pjv*),q 
!               = term1 + term2 + term3 
!               = singrepa()
!==========================================================
	use detci5
!----------------------------------------------------------
! Input:
!  detindx1 = index  of bra determinant
!  elecs    = number of electrons (alpha + beta)
!  aelec    = number of alpha electrons
!  belec    = number of beta electrons
!  orbtls   = number of orbitals
!  alphmat  = matrix of alpha strings
!  betamat  = matrix of beta  strings
!  adets    = number of alpha determinants
!  bdets    = number of beta determinants
!  max1e    = dimension of moints1(:)
!  max2e    = dimension of moints2(:)
!  moints1  = array  of 1-e integrals
!  moints2  = array  of 2-e integrals
!  xorb     = orbital exciting INTO
!  borb     = orbital exciting 
!  nactstr  = string of active determinants
!  live     = number of active determinants
! Variables:
!  p1       = index of alpha string for bra determinant
!  q1       = index of beta string for bra determinant
!  bound    = bound alpha string
!  singx    = single excitation string
!  term1    = < v(pj) | h | v*(p) >
!  term2    = E_{i .ne. j} [v(pi)v(pi)|v(pj)v*(p)] 
!                        - [v(pi)v(pj)|v(pi)v*(p)]
!  term3    = E_n [v(qn)v(qn)|v(pj)v*(p)]
!---------------------------------------------------------- 
!
  implicit none
  integer,intent(in) :: p1,q1, elecs, aelec, belec, orbtls, &
    adets, bdets, max1e, max2e, xorb, borb, live
  integer,dimension(adets,aelec),intent(in) :: alphmat 
  integer,dimension(bdets,belec),intent(in) :: betamat
  real*8,dimension(max1e) :: moints1
  real*8,dimension(max2e) :: moints2
  integer,dimension(live) :: nactstr

  integer :: a, b, c, d, n, i
  integer,dimension(:),allocatable :: bound, singx
  real*8 :: term1, term2, term3
!----------------------------------------------------------
!

! Evaluation of 1 electron term: term1
  term1 = moints1(ind2val(xorb,borb))

! 2-e term 1 (Second term)
  term2 = 0.0
  do i = 1, aelec
 
    if ( alphmat(p1,i) .ne. borb ) then
      term2 = term2 + moints2(index2e2(alphmat(p1,i),alphmat(p1,i), &
         borb, xorb)) - moints2( index2e2(alphmat(p1,i),borb,alphmat(p1,i),xorb))
    end if
 
  end do

! 2-e term 2 (Third term)
  term3 = 0.0
  do n = 1, belec
    term3 = term3 + moints2(index2e2(betamat(q1,n),betamat(q1,n), &
      borb,xorb))
  end do

! Add terms
  singrepa = term1 + term2 + term3

end function singrepa
