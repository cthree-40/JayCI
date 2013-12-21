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
  subroutine lowdiagprecond( diagonals, diaglen, moints1, moints1len, &
    moints2, moints2len, pstring, pstep, plocate, pxreflist, qstring, &
    qstep, qlocate, qxreflist, pdets, qdets, pdetslen, qdetslen, adets,&
    bdets, aelec, belec, orbitals, num_diags, init_dim, outvectors )
    implicit none
    integer, intent(in) :: diaglen, moints1len, moints2len , pdetslen, &
                           qdetslen, adets, bdets, aelec, belec, orbitals,&
                           num_diags, init_dim
    integer, dimension( pdetslen, 2 ), intent(in) :: pstring
    integer, dimension( qdetslen, 2 ), intent(in) :: qstring
    integer, dimension( pdetslen ),    intent(in) :: pstep, plocate, pdets
    integer, dimension( qdetslen ),    intent(in) :: qstep, qlocate, qdets
    integer, dimension( diaglen ),     intent(in) :: qxreflist, pxreflist
    real*8,  dimension( diaglen ),     intent(in) :: diagonals
    real*8,  dimension( moints1len ),  intent(in) :: moints1
    real*8,  dimension( moints2len ),  intent(in) :: moints2
    real*8,  dimension( diaglen, init_dim ), intent(inout) :: outvectors
    real*8,  dimension( diaglen, num_diags) :: vectors1, hvectors
    integer :: i, j, cidim
    real*8  :: ddot
    real*8, dimension( num_diags, num_diags ) :: sub_hammat
    real*8, dimension( num_diags, num_diags ) :: eig_vectors
  !--------------------------------------------------------------------
    cidim = diaglen
    vectors1 = 0d0
    call lowdiags( diagonals, diaglen, num_diags, vectors1 )
    ! Perform Hv on hvectors
    do i=1, num_diags
      call acthv( vectors1(1,i), moints1, moints2, moints1len, moints2len, &
                  pstring, pstep, plocate, pxreflist, qstring, qstep, qlocate, &
                  qxreflist, cidim, pdets, qdets, pdetslen, qdetslen, adets, &
                  bdets, aelec, belec, orbitals, diagonals, hvectors(1,i) )
    end do
    ! Construct vHv
    do i=1, num_diags
      do j=1, num_diags
        sub_hammat(j,i) = ddot( cidim, vectors1(1,j), 1, hvectors(1,i), 1 )
      end do
    end do
    ! Diagonalize this matrix. Return init_dim eigenvectors
    call diag_hamsub( sub_hammat, num_diags, num_diags, eig_vectors )
    print *, "Printing eigvectors..."
    print *, "--"
    do i=1, init_dim
      print *, "--"
      do j=1, num_diags
        print *, eig_vectors(j,i)
      end do
    end do 
    print *, "-==============-"
    print *, "Calling gen_outvecs"
    ! Generate outvectors
    call gen_outvecs( eig_vectors, num_diags, vectors1, diaglen, num_diags, &
                      outvectors, diaglen, init_dim )
    ! Orthonormalize outvectors
    call orthonorm_vecs( outvectors, cidim, init_dim )
    return
  end subroutine
!======================================================================
!======================================================================
!>orthonorm_vecs
!
! Subroutine to orthonormalize set of vectors
!----------------------------------------------------------------------
  subroutine orthonorm_vecs( matrix, ldm, a )
    implicit none
    integer, intent(in) :: ldm, a
    real*8,  dimension( ldm, a ), intent(inout) :: matrix
    integer :: i, j, k
    real*8  :: overlap, ddot, norm
  !--------------------------------------------------------------------
    ! Loop over columns
    do i=1, a
      ! Loop over orthogonalizing columns
      do j=1, i-1
        overlap = ddot( ldm, matrix(1,i), 1, matrix(1,j), 1 )
        do k=1, ldm
          matrix(k,i) = matrix(k,i) - overlap*matrix(k,j)
        end do
      end do
    end do
    ! Normalize
    do i=1, a
      norm = sqrt( ddot( ldm, matrix(1,i), 1, matrix(1,i), 1 ))
      do j=1, ldm
        matrix(j,i) = matrix(j,i) / norm
      end do
    end do
    return
  end subroutine
!======================================================================
!======================================================================
!======================================================================
!======================================================================
!>gen_outvecs
!
! Subroutine to generate vector guess
!----------------------------------------------------------------------
  subroutine gen_outvecs( eig_vectors, dim_eigvec, vectors1, vectors1len, &
    vectors1dim, vectors2, vectors2len, vectors2dim )
    implicit none
    integer, intent(in) :: dim_eigvec, vectors1len, vectors1dim, vectors2len, &
                           vectors2dim
    real*8,  dimension( dim_eigvec, dim_eigvec   ), intent(in)    :: eig_vectors
    real*8,  dimension( vectors1len, vectors1dim ), intent(in)    :: vectors1
    real*8,  dimension( vectors2len, vectors2dim ), intent(out)   :: vectors2
    character(len=1) :: transa, transb
    integer :: lda, ldb, ldc, p, q, r, i, j
    real*8  :: alpha, beta
    real*8,  dimension( vectors2len, dim_eigvec ) :: vectors_int
  !--------------------------------------------------------------------
    ! Set DGEMM variables
    ! aAB + bC = C
    !  A := p x q
    !  B := q x r
    !  C := p x r
    transa = 'n'     ! Do not use the transpose of A
    transb = 'n'     ! Do not use the transpose of B
    p = vectors1len
    q = dim_eigvec
    r = dim_eigvec
    alpha = 1d0
    beta  = 0d0
    lda   = p
    ldb   = q
    ldc   = p
    print *, "Calling DGEMM()"
    call dgemm( transa, transb, p, q, r, alpha, vectors1, lda, eig_vectors, ldb, beta, &
                vectors_int, ldc )
    do i=1, vectors2dim
      do j=1, vectors2len
        vectors2(j,i) = vectors_int(j,i)
      end do
    end do
    return
  end subroutine
!----------------------------------------------------------------------
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
    real*8, dimension( mat_dim, a ),      intent(out) :: eig_vectors
    character(len=1) :: jobz, rnge, uplo
    real*8 :: abstol, dlamch, vl, vu
    integer :: lwork, liwork, il, iu, info, eigfound
    integer, dimension( mat_dim ) :: isuppz
    real*8,  dimension(a) :: eig_values
    real*8,  dimension(26*mat_dim+10) :: work
    real*8,  dimension(11*mat_dim)    :: iwork
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
    print *, "Calling dsyevr()"
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
!  the subblock.
!--------------------------------------------------------------------
  subroutine prediagsubblock( cidim, moints1, moints2, moints1len,  &
    moints2len, pdets, pdetslen, qdets, qdetslen, tdets, subblockdim, initgdim,&
    outvectors, aelec, belec, orbitals )
    implicit none
    integer, intent(in) :: moints1len, moints2len, cidim, pdetslen, qdetslen, &
                           aelec, belec, orbitals, subblockdim, initgdim
    integer, dimension( pdetslen ),   intent(in) :: pdets
    integer, dimension( qdetslen ),   intent(in) :: qdets
    integer, dimension( cidim ),      intent(in) :: tdets
    real*8,  dimension( moints1len ), intent(in) :: moints1
    real*8,  dimension( moints2len ), intent(in) :: moints2
    real*8,  dimension( cidim, initgdim ), intent(out) :: outvectors
    real*8,  dimension( subblockdim, subblockdim) :: sub_hamil
    real*8 :: ham_element
    real*8,  dimension( subblockdim, initgdim ) :: eig_vectors
    real*8,  dimension( cidim, subblockdim )    :: vectors1
    integer :: i, j
  !--------------------------------------------------------------------
    ! Construct subblockdim x subblockdim hamiltonian
    do i=1, subblockdim
      do j=1, subblockdim
        sub_hamil(j,i) = ham_element( tdets(j), tdets(i), moints1, moints1len, &
                                      moints2, moints2len, aelec, belec, orbitals )
      end do
    end do
    print *, "Calling diag_hamsub"
    ! Diagonalize this matrix
    call diag_hamsub( sub_hamil, subblockdim, initgdim, eig_vectors )
    ! Generate guess
    vectors1 = 0d0
    do i=1, subblockdim
      vectors1(i,i) = 1d0
    end do
    print *, "Calling gen_outvecs"
    call gen_outvecs( eig_vectors, subblockdim, vectors1, cidim, subblockdim, &
                      outvectors, cidim, initgdim )
    ! Orthonormalize outvectors
    call orthonorm_vecs( outvectors, cidim, initgdim )
    return
  end subroutine
!====================================================================
!====================================================================
end module 
