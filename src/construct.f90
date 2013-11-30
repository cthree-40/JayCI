subroutine hv_construct( spdim, actstrlen,  moints1, &
  moints2, max1e, max2e, actstr, aelec, belec, totels,  &
  orbitals, adets, bdets, dgls, nfrzn, ndocc, ncas )
!==========================================================
! Explicitly construct hamiltonian
!==========================================================
!
  implicit none
  integer, intent(in) :: spdim, actstrlen, max1e, max2e,&
    aelec, belec, totels, orbitals, adets, bdets, &
    nfrzn, ndocc, ncas
  integer, dimension(actstrlen) :: actstr
  real*8, dimension(max1e) :: moints1
  real*8, dimension(max2e) :: moints2
  real*8, dimension(actstrlen) :: dgls


  real*8, dimension(actstrlen,actstrlen) :: unitmat, hammat

  integer :: i, j, k


  character*1 :: uplo, jobz, rnge
  real*8 :: abstol, ddot, dlamch
  integer :: lwork, liwork, il, iu, info, m, vu, vl
  integer,dimension(:),allocatable :: isuppz
  real*8, dimension(:),allocatable :: work
  real*8, dimension(:),allocatable :: iwork
  real*8, dimension(actstrlen, actstrlen) :: eigvec
  real*8, dimension(actstrlen) :: eigval
!----------------------------------------------------------
  do i=1, actstrlen
    do j=1, actstrlen
      if ( i .eq. j ) then
        unitmat(i,j) = 1d+0
      else
        unitmat(i,j) = 0d+0
      end if
    end do
  end do

  hammat = 0d0

  do i=1, actstrlen
    call acthv( unitmat(1,i), moints1, moints2, max1e, max2e, &
      actstr, actstrlen, aelec, belec, totels, orbitals, adets,&
      bdets, dgls, nfrzn, ndocc, ncas, hammat(1,i) )
  end do

  jobz = 'v'
  rnge = 'a'
  uplo = 'u'
  abstol = dlamch( 'safe minimum' )
  liwork = 12*actstrlen
  lwork = 27*actstrlen
  allocate(isuppz(2*actstrlen))
  allocate(work(lwork))
  allocate(iwork(liwork))

  call dsyevr( jobz, rnge, uplo, actstrlen, hammat, actstrlen, vl, vu, il, &
    iu, abstol, m, eigval, eigvec, actstrlen, isuppz, work, lwork, iwork, liwork, &
    info )
  print *, eigval(1)

  print *, "Finished"
  return
end subroutine hv_construct
!====================================================================
!====================================================================
!> exp_construct
!
! Subroutine to explicitly construct H by finding value of each matrix
!  element.
!--------------------------------------------------------------------
subroutine exp_construct( moints1, moints1len, moints2, moints2len, 
  cidim, aelec, belec, orbitals, pdets, qdets, pdetslen, qdetslen )
  implicit none

!--------------------------------------------------------------------
! Construct hamiltonian
  do i=1, cidim
    do j=1, cidim
      hamiltonian(j,i) = ham_element( j, i, moints1, moints1len, moints2,
                         moints2len, cidim, &
                         aelec, belec, orbitals, pdets, pdetslen, qdets, &
                         qdetslen, ....  )
    end do
  end do
  return
end subroutine exp_construct
!====================================================================
!====================================================================
!> ham_element
!
! real*8 function to compute hamiltonian element(i,j)
!--------------------------------------------------------------------
real*8 function ham_element( ind1, ind2, moints1, moints1len, moints2, &
  moints2len, cidimension, aelec, belec, orbitals, pdets, pdetslen,    &
  qdets, qdetslen, .... )
  use detci2
  implicit none
  integer, intent(in) :: ind1, ind2, moints1len, moints2len, cidimension, &
                         aelec, belec, orbitals, pdetslen, qdetslen

  integer :: p1, q1, p2, q2
  integer, dimension( aelec ) :: pstring1, pstring2
  integer, dimension( belec ) :: qstring1, qstring2

!--------------------------------------------------------------------
! Find determinant string indices for ind1 and ind2
  call k2indc( ind1, belec, orbitals, p1, q1 )
  call k2indc( ind2, belec, orbitals, p2, q2 )
! Find respective strings for p1, q1, p2, q2
  call genorbstring( p1, aelec, orbitals, adets, pstring1 )
  call genorbstring( p2, aelec, orbitals, adets, pstring2 )
  call genorbstring( q1, belec, orbitals, bdets, qstring1 )
  call genorbstring( q2, belec, orbitals, bdets, qstring2 )
! Test differences in strings. If > 2 orbitals, element is 0
