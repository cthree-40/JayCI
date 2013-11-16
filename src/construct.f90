subroutine construct( spdim, actstrlen,  moints1, &
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
end subroutine
