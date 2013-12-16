real*8 function diagonal( detindx, elecs, aelec, belec, &
  orbtls, alphmat, betamat, adets, bdets, max1e, max2e, &
  moints1, moints2)
!==========================================================
! Integer function to compute the diagonal matrix elments 
!
!  < K | H | K > = < p,q | H | p, q>
!                = diag1 + diag2 + diag3 + diag4 + diag5 
!                 = diagonal()
!
!==========================================================
  use detci1
  use detci2
  use detci5
!----------------------------------------------------------
! Inputs:
!  detindx = index of determinant
!  elecs   = total number of electrons
!  aelec   = alpha electrons
!  belec   = beta electrons
!  orbtls  = number of orbitals
!  alphmat = alpha determinant array
!  betamat = beta determinant array
!  adets   = number of alpha determinant strings
!  bdets   = number of beta determinant strings
!  max1e   = dimension of moints1
!  max2e   = dimension of moints2
!  moints1 = 1-e integrals (unique)
!  moints2 = 2-e integrals (unique)
! Variables:
!  p       = index of the alpha determinant string
!  q       = index of the beta determinant string
!  diag1   = E_j <v(pj)|h|v(pj)>
!  diag2   = E_n <v(qn)|h|v(qn)>
!  diag3   = E_i E_{j<i} {[v(pi)v(pi)|v(pj)v(pj)] -
!                            [v(pi)v(pj)|v(pi)v(pj)]}
!  diag4   = E_m E_{n<m} {[v(qm)v(qm)|v(qn)v(qn)] -
!                            [v(qm)v(qn)|v(qm)v(qn)]}
!  diag5   = E_j E_n [v(pj)v(pj)|v(qn)v(qn)]
!----------------------------------------------------------
  implicit none
  integer,intent(in) :: detindx, elecs, aelec, belec, orbtls, &
    adets, bdets, max1e, max2e
  integer,dimension(adets,aelec) :: alphmat
  integer,dimension(bdets,belec) :: betamat
  real*8,dimension(max1e) :: moints1
  real*8,dimension(max2e) :: moints2
  integer :: p, q !index2e, index2e2, ind2val
  real*8 :: diag1, diag2, diag3, diag4, diag5
  
  integer :: i, j, k, l, m !dummy variables
!----------------------------------------------------------
!

#ifdef HV_DIAG
  print *, "We are debugging diagonal()."
  print *, "Determinant index"
  print *, detindx
#endif

! Call k2indc to find the indices of the alpha and beta strings
  call k2indc(detindx, belec, orbtls, p, q)

#ifdef HV_DIAG
  write(*,fmt=10),  &
    "Indices of alpha and beta strings comprising det. indx", p, q
10 format(1x, A, I10, I10)
#endif


! Compute diag1
  diag1 = 0.0
  do i = 1, aelec
    
    j = ind2val(alphmat(p,i),alphmat(p,i))
    diag1 = diag1 + moints1(j)
  
  end do

! Compute diag2
  diag2 = 0.0
  do i = 1, belec
    j = ind2val(betamat(q,i),betamat(q,i))
    diag2 = diag2 + moints1(j)
  end do

! Compute diag3
  diag3 = 0.0
  do i = 2, aelec
    do j = 1, i-1
      k = index2e2(alphmat(p,i),alphmat(p,i),alphmat(p,j),alphmat(p,j))
      l = index2e2(alphmat(p,i),alphmat(p,j),alphmat(p,i),alphmat(p,j))
      diag3 = diag3 + moints2(k) - moints2(l)
    end do
  end do  

! Compute diag4
  diag4 = 0.0
  do i = 2, belec
    do j = 1, i-1
      k = index2e2(betamat(q,i),betamat(q,i),betamat(q,j),betamat(q,j))
      l = index2e2(betamat(q,i),betamat(q,j),betamat(q,i),betamat(q,j))
      diag4 = diag4 + moints2(k) - moints2(l)
    end do
  end do

! Compute diag5
  diag5 = 0.0
  do i = 1, aelec
    do j = 1, belec
      k = alphmat(p,i)
      l = betamat(q,j)
      if ( l > k ) then
        m = index2e2(l,l,k,k)
        diag5 = diag5 + moints2(m)
      else
 	      m = index2e2(k,k,l,l)
	      diag5 = diag5 + moints2(m)
      end if
    end do
  end do

! Add 5 diagonal terms together
  diagonal = diag1 + diag2 + diag3 + diag4 + diag5
 
end function diagonal
