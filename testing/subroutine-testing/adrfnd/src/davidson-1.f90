! Davidson algorithm subroutine
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chemistry, The Johns Hopkins University
!
! Last edit: 10-22-13, This is rewritten from 4.0.0
!====================================================================
!Input:
!  initialspace = initial vectors guess            real*8  array  2-d
!  initdim      = initial dimension of guess       integer scalar
!  actdetslen   = number of active determinants    integer scalar
!  moints1      = 1-e integrals                    real*8  array  1-d
!  moints2      = 2-e integrals                    real*8  array  1-d
!  max1e        = length of 1-e integrals          integer scalar
!  max2e        = length of 2-e integrals          integer scalar
!  actdets      = active determinants              integer array  1-d
!  restol       = residual tolerance               real*8  scalar
!  aelec        = alpha electrons                  integer scalar
!  belec        = beta electrons                   integer scalar
!  totels       = aelec + belec                    integer scalar
!  orbitals     = number of MO's                   integer scalar
!  adets        = number of alpha determinants     integer scalar
!  bdets        = number of beta determinants      integer scalar
!  roots        = number of roots desired          integer scalar
!  krmin        = minimum dimension of kry. space  integer scalar
!  krmax        = maximum dimension of kry. space  integer scalar
!
!Output:
!  energies     = eigenvalues                      real*8  array  1-d
!  eigenvectors = eigenvectors of H                real*8  array  2-d
!====================================================================
#ifdef NOPREDIAG
subroutine davidson( initialspace, initdim, actdetslen, moints1, moints2,  &
  max1e, max2e, actdets, restol, aelec, belec, totels, orbitals, adets,    &
  bdets, ndocc, nfrzn, ncas, xlevel, roots, krmin, krmax, diagonals, miter,&
  energies, eigenvectors )
#endif

#ifdef PREDIAG
subroutine davidson( initdim, actdetslen, moints1, moints2, max1e, max2e, &
  actdets, restol, aelec, belec, totels, orbitals, adets, bdets, ndocc,   &
  nfrzn, ncas, xlevel, roots, krmin, krmax, diagonals, miter,             &
  energies, eigenvectors )
#endif

#ifdef ORTHVEC
  use orthogroutines, only: orthogvector
#endif
#ifdef MODGRAM
  use orthogroutines, only: modgramschmidt
#endif

  implicit none

! ...input integer scalars...
  integer, intent(in) :: initdim, actdetslen, max1e, max2e, aelec, belec, &
                         totels, orbitals, adets, bdets, roots, krmin,    &
                         krmax, nfrzn, ndocc, ncas, xlevel, miter

! ...input real*8 scalars...
  real*8, intent(in) :: restol

! ...input integer arrays...
  integer, dimension(actdetslen) :: actdets
  
! ...input real*8 arrays...
#ifdef NOPREDIAG
  real*8, dimension(actdetslen, initdim) :: initialspace
#endif

  real*8, dimension(actdetslen) :: diagonals
  real*8, dimension(max1e) :: moints1
  real*8, dimension(max2e) :: moints2
  
! ...output real*8 arrays...
  real*8, dimension(roots) :: energies
  real*8, dimension(actdetslen, roots) :: eigenvectors

! ...loop control integer scalars...
  integer :: i, j, l, k, m

! ...loop integer scalars...
  integer :: spdim, curroot
  integer :: oiter

! ...loop real*8 arrays...
  real*8, dimension(actdetslen) :: residual 

#ifdef ORTHVEC
  real*8, dimension(actdetslen) :: newvector
#endif

#ifdef PREDIAG
  real*8, dimension(actdetslen, initdim) :: initialspace
#endif

  real*8, dimension(:), allocatable :: ovlptest
  real*8, dimension(actdetslen, krmax) :: basisvecs, hvectors
  real*8, dimension(krmax) :: kreigvals
  real*8, dimension(krmax,krmax) :: kreigvecs, subhammat

!  real*8, dimension(krmax) :: kreigvals
!  real*8, dimension(krmax, krmax) :: kreigvecs, subhammat

! ...loop real*8 scalars...
  real*8 :: residnorm, ddot
!  real*8 :: ovlptest
#ifdef ORTHVEC
! ...loop real*8 parameters...
  real*8, parameter :: maxovrlp = 1.0D-16
#endif

! ...dsyevr real*8 scalars...
  real*8 :: abstol, dlamch, vl, vu
  
! ...dsyevr characters...
  character*1 :: uplo, jobz, rnge

! ...dsyevr integer scalars...
  integer :: lwork, liwork, il, iu, info, eigfnd
  integer, parameter :: lwmax=100000

! ...dsyevr integer arrays...
  integer, dimension(2*krmax) :: isuppz

! ...dsyevr real*8 arrays...
  real*8, dimension(lwmax) :: work
  real*8, dimension(lwmax) :: iwork


!--------------------------------------------------------------------
!

#ifdef DACTHV
  call dacthv( moints1, moints2, max1e, max2e, actdets, actdetslen,    &
               aelec, belec, totels, orbitals, adets, bdets, diagonals,&
               nfrzn, ndocc, ncas )
  return
#endif 

! Zero out basis space and h-on-basis space
  hvectors = 0d0
  basisvecs = 0d0

! No loop over states, yet - 10-22-2013
  curroot = 1

! Set basis space to initial guess vectors
  spdim = initdim
#ifdef NOPREDIAG
  print *, " You did not select a prediagonalization routine."
  do i=1, spdim
    do j=1, actdetslen
      basisvecs(j,i) = initialspace(j,i)
    end do
  end do
#endif

#ifdef PREDIAG
! Calling prediagonalization routine
  print *, " Calling prediagonalization routine..."
  call subdiag1( initdim, actdetslen, moints1, moints2, max1e, max2e, &
                 actdets, aelec, belec, totels, orbitals, adets, bdets,  &
                 ndocc, nfrzn, ncas, xlevel, krmin, diagonals, basisvecs )
#endif 
  do i=spdim+1,krmax
    do j=1, actdetslen
      basisvecs(j,i) = 0d0
    end do
  end do
! Debug flags for pretty debug print out
#if DEBUG > 0
  print *, "*********************************************************"
  print *, "************ Debugging **********************************"
#endif

  print *, "*********************************************************"
  print *, "**************** DAVIDSON ALGORITHM *********************"
  print *, "*********************************************************"
 
! ...MAIN LOOP...
#if DEBUG > 0
  print *, int( miter/krmax )
#endif
  MAINLP: do i=1, int( miter/krmax )
! Call david1. This subroutine performs H.v for all v in basisvecs
#if DEBUG > 0
    print *, " Calling david1 to perform H.v... "
#endif
    call initdavid( basisvecs, actdetslen, spdim, moints1, moints2,   &
                    max1e, max2e, actdets, aelec, belec, totels,      &
                    orbitals, adets, bdets, diagonals, nfrzn, ndocc,  &
                    ncas, krmax, hvectors )

! Construct the upper triangle of the subspace Hamiltonian v.Hv
! Zero subhammat
    subhammat = 0d0
#if DEBUG > 0
    print *, " Generating upper triangle of subspace Hamiltonian... "
#endif
    do j=1, spdim
      do k=1, spdim
        subhammat( j, k ) = ddot( actdetslen, basisvecs(1,j), 1,  &
                                  hvectors(1, k), 1 )
      end do
    end do

! Diagonalize this matrix using the dsyevr subroutine. This will return
!  both eigenvalues and eigenvectors of v.Hv
! DSYEVR variables must first be set
    jobz = 'v'  ! Return eigenvectors
    rnge = 'a'  ! Return full range of eigenvectors and eigenvalues
    uplo = 'l'  ! upper triangle is stored in matrix A ( subhammat )
    abstol = dlamch( 'Safe minimum' )
    lwork  = lwmax
    liwork = lwmax
! Zero out krylov arrays
    kreigvals = 0d0
    kreigvecs = 0d0
#if DEBUG > 0
    print *, " Diagonalizing subspace Hamiltonian... "
#endif
#if DSYEV = 1
    call dsyevr( jobz, rnge, uplo, spdim, subhammat, krmax, vl, vu,  &
                 il, iu, abstol, eigfnd, kreigvals, kreigvecs, krmax,&
                 isuppz, work, lwork, iwork, liwork, info )
    print *, info
    print *, " Eigenvalues... "
    do j=1, spdim
      print *, kreigvals(j)
    end do
    print *, " Lowest eigenvalue:   ", kreigvals(1)
#else
    print *, " YOU DID NOT SELECT A DIAGONALIZATION SUBROUTINE. "
    print *, " PLEASE DO SO WHEN BUILDING THIS PROGRAM. -DDSYEV="
    stop
#endif    

! Zero out subhammat, because it has been destroyed by dsyev$
    subhammat = 0d0

! ...MAIN SUB-LOOP...
    SUBLP: do j=1, krmax + 1
! Form the residual vector, calling the subroutine residgen()
! Zero out residual
      residual = 0d0
#if DEBUG > 0
      print *, " Constructing residual vector... "
#endif
      do k=1, spdim

        print *, " Coefficient for eigenvector ", k, " is ", kreigvecs(k, curroot)

        do l=1, actdetslen
          residual(l) = residual(l) + kreigvecs(k, curroot)* &
                        ( hvectors(l,k) - kreigvals(curroot)*&
                          basisvecs(l,k))
        end do
      end do
#if DEBUG > 1
      open(unit=10, file='residual.vectors', status='unknown', position='append')
      write(unit=10,fmt=120) "********************************"
120 format(1x, A)
      do k=1, actdetslen
        write(unit=10,fmt=110) residual(k)
      end do
      close(unit=10)
110 format(1x, F15.10)
#endif
! Compute the norm of this residual.
      residnorm = sqrt(ddot( actdetslen, residual, 1, residual, 1 ))
      write(*,fmt=100) "Norm of residual vector:     ", residnorm
100 format(1x, A, F15.10)

! Test norm of residual. If it is lower than rtol, exit subroutine
      if ( residnorm < restol ) then
        exit MAINLP
      end if

! Generate new vector
#if DEBUG > 0
      print *, " Constructing new vector... "
#endif

#ifdef ORTHVEC

! Zero out newvector
      newvector = 0d0
      do k=1, actdetslen
        newvector(k) = - residual(k) / ( diagonals(k) - kreigvals(curroot) + 1D-8)
      end do

      open(unit=11, file='newvector.vectors', status='unknown',position='append')
      write(unit=11,fmt=120) "**********************************"
      do k=1, actdetslen
        write(unit=11,fmt=110) newvector(k)
      end do
      close(unit=11)
    
      print *, " Calling orthogvector()"
      call orthogvector( basisvecs, spdim, krmax, actdetslen, newvector )

! Test orthogonality.
      if (allocated(ovlptest)) deallocate(ovlptest)
      allocate(ovlptest(spdim))
      do k=1, spdim
        ovlptest(k) = ddot( actdetslen, basisvecs(1,k), 1, newvector, 1 )
      end do
      oiter = 1
      do while ( maxval(ovlptest) > maxovrlp .and. oiter .le. 10 )
         call orthogvector( basisvecs, spdim, krmax, actdetslen, newvector )
         do l=1, spdim
           ovlptest(j) = ddot( actdetslen, basisvecs(1,l), 1, newvector, 1)
         end do
         oiter = oiter + 1
      end do
      if (allocated(ovlptest)) deallocate(ovlptest)

      print *, " basisvecs(:,spdim+1) BEFORE addition of newvector..."
      do k=1, actdetslen
        print *, basisvecs(k,spdim+1)
      end do
! Add new vector to space
      do k=1, actdetslen
        basisvecs(k, spdim+1) = 0d0
        basisvecs(k, spdim+1) = basisvecs(k,spdim+1) + newvector(k)
      end do
      print *, " basisvecs(:,spdim+1) AFTER addition of newvector..."
      do k=1, actdetslen
        print *, basisvecs(k,spdim+1)
      end do
! spdim ==> spdim +1
      spdim = spdim + 1
      print *, " spdim is now ", spdim
#endif

#ifdef MODGRAM
! Add new vector to basis space
      do k=1, actdetslen
        basisvecs(k, spdim+1) = - residual(k) / (diagonals(k) - kreigvals(curroot) + 1d-8 )
      end do

      open(unit=11, file='newvector.vectors', status='unknown',position='append')
      write(unit=11,fmt=120) "***********************************"
      do k=1, actdetslen
        write(unit=11, fmt=110) basisvecs(k,spdim+1)
      end do
      close(unit=11)


! spdim ==> spdim +1
      spdim = spdim + 1
      print *, " spdim is now ", spdim
! Orthogonalize space. Call modgramschmidt()
      print *, "Calling modgramschmidt()"
      call modgramschmidt( basisvecs, spdim, krmax, actdetslen )

#endif

! Perform Hv for new vector
      call acthv( basisvecs(1,spdim), moints1, moints2, max1e, max2e, actdets,&
                  actdetslen, aelec, belec, totels, orbitals, adets, bdets,   &
                  nfrzn, ndocc, ncas, diagonals, hvectors(1,spdim ))


! Construct subspace Hamiltonian
! Zero out subhammat
      subhammat = 0d0
        
      do k=1, spdim
        do l=1, spdim
          subhammat(l,k) = ddot( actdetslen, basisvecs(1,l), 1, hvectors(1,k), 1 )
        end do
      end do

! Diagonalize subhammat
! First, set DSYEVR variables
      jobz = 'v'  ! Return eigenvectors
      rnge = 'a'  ! Return full range of eigenvectors and eigenvalues
      uplo = 'l'  ! upper triangle is stored in matrix A ( subhammat )
      abstol = dlamch( 'Safe minimum' )
      lwork  = lwmax
      liwork = lwmax
! Zero out krylov arrays
      kreigvals = 0d0
      kreigvecs = 0d0
#if DEBUG > 0
      print *, " Diagonalizing subspace Hamiltonian... "
#endif
#if DSYEV = 1
      call dsyevr( jobz, rnge, uplo, spdim, subhammat, krmax, vl, vu, &
                   il, iu, abstol, eigfnd, kreigvals, kreigvecs, krmax,&
                   isuppz, work, lwork, iwork, liwork, info ) 
      print *, " Eigenvalues... "
      do k=1,spdim
        print *, kreigvals(k)
      end do
      print *, " Lowest eigenvalue:  ", kreigvals(1)
#else
      print *, " YOU DID NOT SELECT A DIAGONALIZATION SUBROUTINE. "
      print *, " PLEASE DO SO WHEN BUILDING THIS PROGRAM. -DDSYEV="
      stop
#endif

! Zero out subhammat, because it has been destroyed by dsyev$
      subhammat = 0d0

! Check size of krylov space
      if ( spdim .eq. krmax ) then
#if DEBUG > 0
        print *, " Truncating the space... "
#endif
        call spacetrunc( basisvecs, actdetslen, krmin, krmax, curroot, &
          kreigvecs )
        spdim = krmin
! Zero out hspace
        hvectors = 0d0
        exit SUBLP
      end if
    end do SUBLP
  end do MAINLP

  return
end subroutine

        
