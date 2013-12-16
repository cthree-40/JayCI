! Module containing prediagonalization routines
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chem., The Johns Hopkins Unversity
!
! Last edit: 11-20-13
!====================================================================
module prediag

contains
!====================================================================
!====================================================================
!> geninituvec
!
! Subroutine to generate initial guess unit vectors
!--------------------------------------------------------------------
  subroutine geninituvec( numvec, lenvec, initvecs )
! Input:
!  numvec = number of desired vectors
!  lenvec = length of desired vectors
! Output:
!  initvecs = numvec unit vectors
!--------------------------------------------------------------------
    implicit none
! ...input integer scalars...
    integer, intent(in) :: numvec, lenvec
! ...output real*8 arrays...
    real*8, dimension(lenvec,numvec), intent(out) :: initvecs
! ...loop integer scalars...
    integer :: i,j
!--------------------------------------------------------------------
    initvecs = 0d0
    do i=1, numvec
      initvecs(i,i) = 1d0
    end do
    return
  end subroutine geninituvec
!====================================================================
!====================================================================
!> lowdiagprecond
!
! Subroutine to generate intial guess unit vectors by grabbing the 
!  lowest energy diagonals, and diagonalizing the v.Hv matrix
!--------------------------------------------------------------------
  subroutine lowdiagprecond( diagonals, cidim, moints1, moints2, moints1len,      &
    moints2len, pstring, pstep, plocate, qstring, qstep, qlocate, pdets,          &
    qdets, pdetslen, qdetslen, adets, bdets, aelec, belec, orbitals, initguessdim,&
    intdiagsguess, outvectors )
!Input:
! diagonals          = <K|H|K>                                 real*8  array  1-d
! cidim              = dimension of CI space                   integer scalar
! moints1            = 1-e integtrals                          real*8  array  1-d
! moints2            = 2-e integerals                          real*8  array
! moints1len         = length of moints1                       integer scalar
! moints2len         = length of moints2                       integer scalar
! pstring            = det string pairs ( p leading )          integer array  2-d
! pstep              = pstring step array                      integer array  1-d
! plocate            = location of p in pstring                integer array  1-d
! qstring            = det string pairs ( q leading )          integer array  2-d
! qstep              = qstring step array                      integer array  1-d
! qlocate            = location of q in qstring                integer array  1-d
! pdets              = p-strings                               integer array  1-d
! qdets              = q-strings                               integer array  1-d
! pdetslen           = p-strings length                        integer scalar
! qdetslen           = q-strings length                        integer scalar
! adets              = # of untruncated alpha strings          integer scalar
! bdets              = # of untruncated beta  strings          integer scalar
! orbitals           = # of MO's                               integer scalar
! intdiagsguess      = # of diagonals to generate guess with   integer scalar
! initguessdim       = # of guess vectors to return            integer scalar
!Output:
! outvectors         = output guess vectors
!--------------------------------------------------------------------
    implicit none
! ...input integer scalars...
    integer, intent(in) :: cidim, moints1len, moints2len, pdetslen, qdetslen, adets, &
                           bdets, orbitals, intdiagsguess, initguessdim
! ...input integer arrays...
    integer, dimension( cidim, 2 ), intent(in)           :: pstring, qstring
    integer, dimension( pdetslen ), intent(in)           :: pdets, plocate
    integer, dimension( qdetslen ), intent(in)           :: qdets, qlocate
! ...input real*8 arrays...
    real*8, dimension( moints1len ), intent(in)          :: moints1
    real*8, dimension( moints2len ), intent(in)          :: moints2
    real*8, dimension(   cidim    ), intent(in)          :: diagonals
! ...output real*8 arrays...
    real*8, dimension( cidim, initguessdim), intent(out) :: outvectors
! ...loop integer scalars...
    integer :: i, j
! ...real*8 scalars...
    real*8 :: ddot
! ...real*8 arrays...
    real*8, dimension( cidim, initguessdim )             :: intguess, hvectors
    real*8, dimension( initguessdim, initguessdim )      :: subham
! ...DSYEVR VARIABLES...
    character(len=1)   :: uplo, jobz, rnge
    real*8             :: abstol, dlamch, vl, vu
    integer            :: lwork, liwork, il, iu, info, eigfound
    integer, dimension(:),   allocatable :: isuppz 
    real*8,  dimension(:,:), allocatable :: eigvectors
    real*8,  dimension(:),   allocatable :: eigvalues, work, iwork
! ...DGEMM VARIABLES...
    character(len=1)   :: transa, transb
    integer            :: lda, ldb, ldc, p, q, r
    real*8             :: alpha, beta
!--------------------------------------------------------------------
! Get lowest diagonals
    intguess = 0d0
    call lowdiags( diagonals, cidim, intdiagsguess, intguess )
! Perform Hv, generating hvectors
    hvectors = 0d0
    do i=1, intdiagsguess
      call acthv( intguess(1,i), moints1, moints2, moints1len, moints2len, &
                  pstring, pstep, plocate, qstring, qstep, qlocate, cidim,  &
                  pdets, qdets, pdetslen, qdetslen, adets, bdets, aelec,    &
                  belec, orbitals, diagonals, hvectors(1,i) )
    end do
! Construct v.Hv
    do i=1, intdiagsguess
      do j=1, intdiagsguess
        subham(j,i) = ddot( cidim, intguess(1,j), 1, hvectors(1,i), 1 ) 
      end do
    end do
! Set DSYEVR variabls
    jobz  = 'v'    ! Return eigenvalues
    rnge  = 'a'    ! Return full range of eigenvectors and eigenvalues
    uplo  = 'l'    ! Lower triangle is stored in matrix A
    abstol= dlamch( 'Safe minimum' )
    lwork = 26*intdiagsguess+10
    liwork= 11*intdiagsguess
    allocate( isuppz(intdiagsguess) )
    allocate( eigvectors(intdiagsguess, intdiagsguess) )
    allocate(  work( lwork) )
    allocate( iwork(liwork) )
    allocate( eigvalues(intdiagsguess) )
    call dsyevr( jobz, rnge, uplo, intdiagsguess, subham, intdiagsguess, vl, &
                vu, il, iu, abstol, eigfound, eigvalues, eigvectors, intdiagsguess,&
                isuppz, work, lwork, iwork, liwork, info )
! Generate eigenvector guess...
    transa = 'n'
    transb = 'n'
    p = cidim
    q = intdiagsguess
    r = q
    alpha = 1d0
    beta  = 0d0
    lda   = p
    ldb   = initguessdim
    ldc   = p
    call dgemm( transa, transb, p, q, r, alpha, intguess, lda, eigvectors, ldb, beta, &
                outvectors, ldc )
! Return
    return
  end subroutine
!====================================================================
!====================================================================
!> prediagsubblock
!
! Subroutine that explicitly constructs matrix elements Hij and diagonalizes
!  the subblock.
!--------------------------------------------------------------------
  subroutine prediagsubblock( diagonals, cidim, moints1, moints2, moints1len,  &
    moints2len, pdets, pdetslen, qdets, qdetslen, tdets, subblockdim, initgdim,&
    outvectors )
!Input:
! diagonals        = <K|H|K>
! cidim            = dimension of CI space
! moints1          = 1-e integrals
! moints2          = 2-e integrals
! moints1len       = length of 1-e integrals
! moints2len       = length of 2-e integrals
! pdets            = alpha determinants (truncated)
! pdetslen         = length of pdets
! qdets            = beta determinants  (truncated)
! qdetslen         = length of qdets
! tdets            = determinants
! initgdim         = number of vectors to return
! subblockdim      = dimension of H subblock
!Output:
! outvectors      = output vectors
    implicit none
! ...input integer scalars...
    integer, intent(in) :: cidim, moints1len, moints2len, pdetslen, qdetslen, &
                           subblockdim, initgdim
! ...input integer arrays...
    integer, dimension( pdetslen ), intent(in)   :: pdets
    integer, dimension( qdetslen ), intent(in)   :: qdets
    integer, dimension( cidim),     intent(in)   :: tdets
! ...input real*8 arrays...
    real*8,  dimension( moints1len ), intent(in) :: moints1
    real*8,  dimension( moints2len ), intent(in) :: moints2
! ...output real*8 arrays...
    real*8, dimension( cidim, initgdim ), intent(inout) :: outvectors
! ...loop integer scalars...
    integer :: i, j
! ...real*8 arrays...
    real*8, dimension( subblockdim, subblockdim ) :: hamblock
    real*8, dimension( cidim, inigdim )           :: unitvecs
! ...DSYEVR variables...
    character(len=1) :: jobz, uplo, rnge
    integer          :: lwork, liwork, il, iu, info, eigfound
    real*8           :: abstol, dlamch, vl, vu
    integer, dimension(:),   allocatable :: isuppz 
    real*8,  dimension(:,:), allocatable :: eigvectors
    real*8,  dimension(:),   allocatable :: eigvalues, work, iwork
! ...DGEMM variables...
! ...DGEMM VARIABLES...
    character(len=1)   :: transa, transb
    integer            :: lda, ldb, ldc, p, q, r
    real*8             :: alpha, beta
!--------------------------------------------------------------------
! Construct subblockdim x subblockdim hamiltonian
    do i=1, subblockdim
      do j=1, subblockdim
        hamblock(j,i) = ham_element( tdets(j), tdets(i), moints1, moints1len, &
                                     moints2, moints2len, aelec, belec, orbitals )
      end do
    end do
! Diagonalize hamblock using DSYEVR
! Set DSYEVR variables
    jobz  = 'v'    ! Return eigenvalues
    rnge  = 'a'    ! Return full range of eigenvectors and eigenvalues
    uplo  = 'l'    ! Lower triangle is stored in matrix A
    abstol= dlamch( 'Safe minimum' )
    lwork = 26*subblockdim+10
    liwork= 11*subblockdim
    allocate( isuppz(subblockdim) )
    allocate( eigvectors(subblockdim, subblockdim) )
    allocate(  work( lwork) )
    allocate( iwork(liwork) )
    allocate( eigvalues(subblockdim) )
    call dsyevr( jobz, rnge, uplo, subblockdim, hamblock, subblockdim, vl, &
                vu, il, iu, abstol, eigfound, eigvalues, eigvectors, subblockdim,&
                isuppz, work, lwork, iwork, liwork, info )
! Generate unit vectors
    unitvec = 0d0
    do i=1, initgdim
      unitvecs(i,i) = 1d0
    end do
! Set DGEMM variables
    transa = 'n'
    transb = 'n'
    p = cidim
    q = subblockdim
    r = q
    alpha = 1d0
    beta  = 0d0
    lda   = p
    ldb   = initgdim
    ldc   = p
    call dgemm( transa, transb, p, q, r, alpha, unitvecs, lda, eigvectors, ldb, beta, &
                outvectors, ldc )
! Return
    return
  end subroutine
!====================================================================
!====================================================================
end module 
