! Module containing prediagonalization routines
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chem., The Johns Hopkins Unversity
!
! Last edit: 11-20-13
!====================================================================
module prediag
  implicit none
contains
!====================================================================
!====================================================================
!> geninituvec
!
! Subroutine to generate initial guess unit vectors
!--------------------------------------------------------------------
  subroutine geninituvec( numvec, lenvec, initvecs )
! Input:
! numvec = number of desired vectors
! lenvec = length of desired vectors
! Output:
! initvecs = numvec unit vectors
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
! lowest energy diagonals, and diagonalizing the v.Hv matrix
!--------------------------------------------------------------------
  subroutine lowdiagprecond( diagonals, diaglen, moints1, moints1len, &
    moints2, moints2len, pstring, pstep, plocate, qstring, &
    qstep, qlocate, pdets, qdets, pdetslen, qdetslen, adets,&
    bdets, aelec, belec, orbitals, num_diags, init_dim, outvectors )
    implicit none
    integer, intent(in) :: diaglen, moints1len, moints2len , pdetslen, &
                           qdetslen, adets, bdets, aelec, belec, orbitals,&
                           num_diags, init_dim
    integer, dimension( pdetslen, 2 ), intent(in) :: pstring
    integer, dimension( qdetslen, 2 ), intent(in) :: qstring
    integer, dimension( pdetslen ), intent(in) :: pstep, plocate, pdets
    integer, dimension( qdetslen ), intent(in) :: qstep, qlocate, qdets
    real*8, dimension( moints1len ), intent(in) :: moints1
    real*8, dimension( moints2len ), intent(in) :: moints2
    real*8, dimension( diaglen ), intent(in) :: diagonals
    real*8, dimension( diaglen, init_dim ), intent(inout) :: outvectors
    real*8, dimension( diaglen, num_diags) :: vectors1, hvectors
    integer :: i,j, cidim
    real*8 :: ddot
    real*8, dimension( num_diags, num_diags ) :: sub_hammat, eig_vectors
  !--------------------------------------------------------------------
    cidim = diaglen
    call lowdiags( diagonals, diaglen, num_diags, vectors1 )
    ! Perform Hv on hvectors
    do i=1, num_diags
      call acthv( vectors1(1,i), moints1, moints2, moints1len, moints2len, &
                  pstring, pstep, plocate, qstring, qstep, qlocate, &
                  cidim, pdets, qdets, pdetslen, qdetslen, adets, &
                  bdets, aelec, belec, orbitals, diagonals, hvectors(1,i) )
    end do
    ! Construct vHv
    do i=1, num_diags
      do j=1, num_diags
        sub_hammat(j,i) = ddot( cidim, vectors1(1,j), 1, hvectors(1,i), 1 )
      end do
    end do
    ! Diagonalize this matrix. Return init_dim vectors
    call diag_hamsub( sub_hammat, num_diags, init_dim, eig_vectors )
    return
  end subroutine
!======================================================================
!======================================================================
!>diag_hamsub
!
! Subroutine to diagonalize subblock of hamiltonian using DSYEVR
!----------------------------------------------------------------------
  subroutine diag_hamsub( matrix, mat_dim, a, eig_vectors )
    implicit none
    integer, intent(in) :: mat_dim, a
    real*8, dimension( mat_dim, mat_dim ), intent(in) :: matrix
    real*8, dimension( mat_dim, a ), intent(in) :: eig_vectors
    character(len=1) :: jobz, rnge, uplo
    real*8 :: abstol, dlamch, vl, vu
    integer :: lwork, liwork, il, iu, info, eigfound
    integer, dimension(:), allocatable :: isuppz
    real*8, dimension(a) :: eig_values
    real*8, dimension(:), allocatable :: work, iwork
!----------------------------------------------------------------------
    ! Set dsyever variables
    jobz = 'v'
    rnge = 'i'
    il = 1
    iu = a
    uplo = 'l'
    abstol = dlamch( 'Safe minimum' )
    lwork = 26*mat_dim+10
    liwork = 11*mat_dim
    allocate( isuppz( mat_dim ) )
    allocate( work( lwork ) )
    allocate( iwork(liwork) )
    call dsyevr( jobz, rnge, uplo, mat_dim, matrix, mat_dim, vl, &
                 vu, il, iu, abstol, eigfound, eig_values, eig_vectors, &
                 mat_dim, isuppz, work, lwork, iwork, liwork, info )
    return
  end subroutine
!======================================================================
!======================================================================
    
!====================================================================
!====================================================================
!> prediagsubblock
!
! Subroutine that explicitly constructs matrix elements Hij and diagonalizes
! the subblock.
!--------------------------------------------------------------------
  subroutine prediagsubblock( cidim, moints1, moints2, moints1len, &
    moints2len, pdets, pdetslen, qdets, qdetslen, tdets, subblockdim, initgdim,&
    aelec, belec, orbitals, outvectors )
!Input:
! diagonals = <K|H|K>
! cidim = dimension of CI space
! moints1 = 1-e integrals
! moints2 = 2-e integrals
! moints1len = length of 1-e integrals
! moints2len = length of 2-e integrals
! pdets = alpha determinants (truncated)
! pdetslen = length of pdets
! qdets = beta determinants (truncated)
! qdetslen = length of qdets
! tdets = determinants
! initgdim = number of vectors to return
! subblockdim = dimension of H subblock
!Output:
! outvectors = output vectors
    implicit none
! ...input integer scalars...
    integer, intent(in) :: cidim, moints1len, moints2len, pdetslen, qdetslen, &
                           subblockdim, initgdim, aelec, belec, orbitals
! ...input integer arrays...
    integer, dimension( pdetslen ), intent(in) :: pdets
    integer, dimension( qdetslen ), intent(in) :: qdets
    integer, dimension( cidim), intent(in) :: tdets
! ...input real*8 arrays...
    real*8, dimension( moints1len ), intent(in) :: moints1
    real*8, dimension( moints2len ), intent(in) :: moints2
! ...output real*8 arrays...
    real*8, dimension( cidim, initgdim ), intent(inout) :: outvectors
! ...loop integer scalars...
    integer :: i, j
! ...real*8 arrays...
    real*8, dimension( subblockdim, subblockdim ) :: hamblock
    real*8, dimension( cidim, initgdim ) :: unitvecs
! ...DSYEVR variables...
    character(len=1) :: jobz, uplo, rnge
    integer :: lwork, liwork, il, iu, info, eigfound
    real*8 :: abstol, dlamch, vl, vu
    integer, dimension(:), allocatable :: isuppz
    real*8, dimension(:,:), allocatable :: eigvectors
    real*8, dimension(:), allocatable :: eigvalues, work, iwork
! ...DGEMM variables...
! ...DGEMM VARIABLES...
    character(len=1) :: transa, transb
    integer :: lda, ldb, ldc, p, q, r
    real*8 :: alpha, beta
! ...real*8 scalars...
    real*8 :: ham_element
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
    jobz = 'v' ! Return eigenvalues
    rnge = 'a' ! Return full range of eigenvectors and eigenvalues
    uplo = 'l' ! Lower triangle is stored in matrix A
    abstol= dlamch( 'Safe minimum' )
    lwork = 26*subblockdim+10
    liwork= 11*subblockdim
    allocate( isuppz(subblockdim) )
    allocate( eigvectors(subblockdim, subblockdim) )
    allocate( work( lwork) )
    allocate( iwork(liwork) )
    allocate( eigvalues(subblockdim) )
    call dsyevr( jobz, rnge, uplo, subblockdim, hamblock, subblockdim, vl, &
                vu, il, iu, abstol, eigfound, eigvalues, eigvectors, subblockdim,&
                isuppz, work, lwork, iwork, liwork, info )
! Generate unit vectors
    unitvecs = 0d0
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
    beta = 0d0
    lda = p
    ldb = initgdim
    ldc = p
    call dgemm( transa, transb, p, q, r, alpha, unitvecs, lda, eigvectors, ldb, beta, &
                outvectors, ldc )
! Return
    return
end subroutine
!====================================================================
!====================================================================
end module 
