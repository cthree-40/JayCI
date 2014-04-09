!Subroutine to compute inital guess vectors for davison algorithm by
! constructing & diagonalizing the MxM upper block of the Hamiltonian
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chemistry, The Johns Hopkins University
!
! Last edit: 10-24-2013
!====================================================================
!Input:
! m           = size of space to diagonalize            integer scalar
! actdetslen  = length of actdets array (No. of dets)   integer scalar
! moints1     = 1-e integrals                           real*8  array
! moints2     = 2-e integrals                           real*8  array
! max1e       = length of moints1                       integer scalar
! max2e       = length of moints2                       integer scalar
! actdets     = array of 'active' determinants          integer array
! aelec       = number of alpha electrons               integer scalar
! belec       = number of beta electrons                integer scalar
! totels      = aelec + belec                           integer scalar
! orbitals    = number of MO's                          integer scalar
! adets       = number of alpha determinants            integer scalar
! bdets       = number of beta determinants             integer scalar
! ndocc       = number of doubly occupied orbitals      integer scalar
! nfrzn       = number of frozen orbitals               integer scalar
! ncas        = number of CAS orbitals                  integer scalar
! xlevel      = excitations level                       integer scalar
! krmin       = minimum dimension of krylov subspace    integer scalar
!               Columns of _initguess_
! dgls        = < k | H | k >                           real*8  array
!Output:
! initguess   = initial guess eigenvectors of H         real*8  array
!====================================================================
subroutine subdiag1( m, actdetslen, moints1, moints2, max1e, max2e,   &
  actdets, aelec, belec, totels, orbitals, adets, bdets, ndocc, nfrzn,&
  ncas, xlevel, krmin, dgls, initguess )

  implicit none

! ...input integer scalars...
  integer, intent(in) :: m, actdetslen, max1e, max2e, aelec, belec,    &
                         totels, orbitals, adets, bdets, ndocc, nfrzn, &
                         ncas, xlevel, krmin

! ...input integer arrays...
  integer, dimension(actdetslen) :: actdets

! ...input real*8 arrays...
  real*8, dimension(max1e) :: moints1
  real*8, dimension(max2e) :: moints2
  real*8, dimension(actdetslen) :: dgls
  
! ...output real*8 arrays...
  real*8, dimension(actdetslen, krmin) :: initguess

! ...loop integer scalars...
  integer :: i, j, k, l

! ...loop real*8 scalars...
  real*8 :: ddot

! ...loop real*8 arrays...
  real*8, dimension(m, m) :: subblockham
  real*8, dimension(actdetslen, m) :: unitguess, hvectors1, &
                                      vectors1

! ...DSYEVR real*8 scalars...
  real*8 :: abstol, dlamch, vl, vu

! ...DSYEVR characters...
  character*1 :: uplo, jobz, rnge

! ...DSYEVR integer scalars...
  integer :: lwork, liwork, il, iu, info, eigfound
  integer, parameter :: lwmax=50000

! ...DSYEVR integer arrays...
  integer, dimension(20) :: isuppz

! ...DSYEVR real*8 arrays...
  real*8, dimension(lwmax) :: work
  real*8, dimension(lwmax) :: iwork
  real*8, dimension( m ) :: eigvalues
  real*8, dimension(m, m) :: eigvectors

! ...DGEMM characters...
  character*1 :: transa, transb

! ...DGEMM integer scalars...
  integer :: p, q, r, lda, ldb, ldc

! ...DGEMM real*8 scalars...
  real*8 :: alpha, beta

! ...initial guess gen. real*8 arrays...
  real*8, dimension(actdetslen) :: hartreefock, honhartfock, &
                                   vector2

! ...initial guess gen. integer arrays...
  integer, dimension(m) :: hfelements

! ...initial guess gen. integer scalars...
  integer :: indice

!--------------------------------------------------------------------

! Generate 'unit' matrix 
! Select m greatest elements of H on HF
   hartreefock = 0d0
   hartreefock(1) = 1d0
   call acthv( hartreefock, moints1, moints2, max1e, max2e, actdets,&
               actdetslen, aelec, belec, totels, orbitals, adets,   &
               bdets, nfrzn, ndocc, ncas, dgls, honhartfock )

! Construct vector guesses
   unitguess=0d0
   call lowdiags( dgls, actdetslen, m, unitguess ) 
! Perform Hv, generating hvectors1
  do i=1, m
    call acthv( unitguess(1,i), moints1, moints2, max1e, max2e,      &
                actdets, actdetslen, aelec, belec, totels, orbitals, &
                adets, bdets, nfrzn, ndocc, ncas, dgls, hvectors1(1,i))
  end do

! Construct upper triangle of v.Hv matrix
  do i=1, m
    do j=1, m
      subblockham(j,i) = ddot( actdetslen, unitguess(1,j), 1, &
                               hvectors1(1,i), 1 )
    end do
  end do
  print *, " Printing subblock Hamiltonian, this should have overlaps..."
  do i=1, m
      print *, subblockham(i,:)
  end do

! Set DSYEVR variables
  jobz = 'v'  ! Return eigenvectos
  rnge = 'a'  ! Return full range of eigenvectors and eigenvalues
  uplo = 'l'  ! upper triangle is stored in matrix A
  abstol = dlamch( 'Safe minimum' )
  lwork  = lwmax
  liwork = lwmax

! Diagonalize subblockham
  call dsyevr( jobz, rnge, uplo, m, subblockham, m, vl, vu,     &
               il, iu, abstol, eigfound, eigvalues, eigvectors, &
               m, isuppz, work, lwork, iwork, liwork, info )
  print *, " Number of eigenvalues found:  ", eigfound
  do i=1, m
    print *, " Eigenvalue   ", i, " is:  ", eigvalues(i)
  end do

  print *, " The first eigenvector is..."
  do i=1, m
    print *, eigvectors(i,1)
  end do 
  if ( info .ne. 0 ) then
    print *, " DSYEVR did not converge!!!!! ( initdiag1 ) "
  end if

! Generate eigenvector guesses...
  transa = 'n'
  transb = 'n' 
  p = actdetslen
  q = m
  r = m
  alpha = 1d0
  beta = 0d0
  lda = p
  ldb = r
  ldc = p

! Call DGEMM
  print *, " Calling DGEMM... "
  call dgemm( transa, transb, p, q, r, alpha, unitguess, lda, &
              eigvectors, ldb, beta, vectors1, ldc )
  print *, " Finished DGEMM..."

! Generate initial guess initguess( actdetslen, krmin)
  do i=1, krmin
    print *, "Writing ", i, "th vector...."
    do j=1, actdetslen
      initguess(j,i) = vectors1(j,i)
    end do
  end do

! Return initial guess
  return

end subroutine
