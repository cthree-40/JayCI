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
    real*8, dimension( diaglen, init_dim ), intent(out) :: outvectors
    real*8, dimension(:,:), allocatable :: vectors1, hvectors
    integer :: i,j, cidim
    real*8, external :: ddot
    real*8, dimension(:,:), allocatable :: sub_hammat
    real*8, dimension(:,:), allocatable :: eig_vectors
  !--------------------------------------------------------------------
  ! Allocate arrays
    allocate( vectors1( diaglen, num_diags ) )
    allocate( hvectors( diaglen, num_diags ) )
    allocate( sub_hammat( num_diags, num_diags ) )
    allocate( eig_vectors( num_diags, init_dim ) )

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

    ! Generate vectors
    call ritz_vec( eig_vectors, cidim, num_diags, init_dim, vectors1, outvectors )    
  ! Deallocate arrays
    deallocate( vectors1, hvectors, sub_hammat, eig_vectors )
    
    return
  end subroutine
!======================================================================
!======================================================================
!>ritz_vec
!
! Subroutine to generate the ritz vector. (E)_nxm (alpha)_mxo = (VEC)_nxo
!----------------------------------------------------------------------
  subroutine ritz_vec( alpha, n, m, o, Evec, VEC )
    implicit none
    integer, intent(in) :: n, m, o
    real*8,  dimension( n, m ), intent(in) :: Evec 
    real*8,  dimension( m, o ), intent(in) :: alpha
    real*8,  dimension( n, o ), intent(out):: VEC
    character(len=1) :: transa, transb
    integer :: lda, ldb, ldc, p, q, r
    real*8 :: a, b
!----------------------------------------------------------------------
    transa = 'n'   ! Use transpose of A
    transb = 'n'   ! Use transpose of B
    p = n
    q = o
    r = m
    a = 1d0
    b = 0d0
    lda = n
    ldb = m
    ldc = n
    call dgemm( transa, transb, p, q, r, a, Evec, lda, alpha, ldb, b, &
                VEC, ldc )
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
    real*8, dimension( mat_dim, mat_dim ),intent(inout) :: matrix
    real*8, dimension( mat_dim, a ),        intent(out) :: eig_vectors
    character*1 :: jobz, rnge, uplo
    real*8 :: abstol, vl, vu
    real*8, external :: dlamch
    integer :: lwork, liwork, il, iu, info, eigfound
    integer, dimension(:),  allocatable :: isuppz
    real*8, dimension(:),   allocatable :: eig_values
    real*8, dimension(:,:), allocatable :: dys_eigvec
    real*8, dimension(:),   allocatable :: work
    real*8, dimension(:),   allocatable :: iwork
    integer :: i, j
    integer :: ispzstat, eigstat, dysstat, workstat, iworkstat, &
               allocstat
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate( isuppz( 2*mat_dim ),stat=allocstat )
    print *, "isuppz=", allocstat
    print *, size(isuppz)
    allocate( eig_values( mat_dim ),stat=allocstat )
    print *, "eig_values=", allocstat
    allocate( dys_eigvec( mat_dim, mat_dim ), stat=allocstat )
    print *, "dys_eigvec=", allocstat
    allocate( work( 27*mat_dim ), stat=allocstat )
    print *, "work=", allocstat
    allocate(iwork( 12*mat_dim ),stat=allocstat )
    print *, "iwork=", allocstat
    ! Set dsyever variables
    jobz = 'v'
    rnge = 'a'
    uplo = 'u'
    abstol = dlamch( 'Safe minimum' )
    lwork = 26*mat_dim
    liwork = 11*mat_dim
    open( unit=97, file='diag-debug', status='new' )
    do i=1, mat_dim
      do j=1, mat_dim
        write(unit=97,fmt=10) j, i, matrix(j,i)
      end do
    end do
10 format(1x,I4,I4,F10.6)
    close(unit=97)
    print *, " Calling dsyevr..."
    call dsyevr( jobz, rnge, uplo, mat_dim, matrix, mat_dim, vl, &
                 vu, il, iu, abstol, eigfound, eig_values, dys_eigvec, &
                 mat_dim, isuppz, work, lwork, iwork, liwork, info )
    print *, " Finished dsyevr..."
    print *, "info=",info
    print *, "eigfound=",eigfound
    print *, "eigVALUE=",eig_values(1)
    do i=1, mat_dim
      print *, dys_eigvec(i,1)
    end do
    do i=1, a
      do j=1, mat_dim
        eig_vectors(j,i) = dys_eigvec(j,i)
      end do
    end do 
  ! Deallocate
    deallocate( isuppz, stat=ispzstat )
    print *, " isuppz=", ispzstat
    deallocate( eig_values, stat=eigstat)
    print *, " eig_values=", eigstat
    deallocate( dys_eigvec, stat=dysstat)
    print *, " dys_eigvec=", dysstat
    deallocate( work, stat=workstat )
    print *, " work=", workstat
    deallocate( iwork, stat=iworkstat)
    print *, "iwork=", iworkstat

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
    use construct, only: ham_element
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
! tdets = determinant list ( truncated )
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
    real*8, dimension( cidim, initgdim ), intent(out) :: outvectors
! ...loop integer scalars...
    integer :: i, j
! ...real*8 arrays...
    real*8, dimension(:,:), allocatable    :: hamblock
    real*8, dimension(:,:), allocatable    :: unitvecs
    real*8, dimension(:,:), allocatable    :: vectors2
! ...real*8 scalars...
!    real*8 :: ham_element
!--------------------------------------------------------------------
! Allocate arrays
    allocate( hamblock( subblockdim, subblockdim ) )
    allocate( unitvecs( cidim, subblockdim ) )
    allocate( vectors2( cidim, subblockdim ) )

! Construct subblockdim x subblockdim hamiltonian
    do i=1, subblockdim
      do j=1, subblockdim
        hamblock(j,i) = ham_element( tdets(j), tdets(i), moints1, moints1len, &
                                     moints2, moints2len, aelec, belec, orbitals )
        print *, j, i, hamblock(j,i)
      end do
    end do
   
    print *, "Calling diag_hamsub..."
! Diagonalize hamblock using DSYEVR
    call diag_hamsub( hamblock, subblockdim, initgdim, vectors2 ) 
    print *, "Finished diag_hamsub..."
! Generate unit vectors
    do i=1, subblockdim
      do j=1, cidim
        unitvecs(j,i) = 0d0
      end do
      unitvecs(i,i) = 1d0
    end do
    print *, "Calling ritz_vec..."
! Set DGEMM variables
    call ritz_vec( vectors2, cidim, subblockdim, initgdim, unitvecs, outvectors )

! Deallocate arrays
    deallocate( hamblock, unitvecs, vectors2 )

    return
  end subroutine
!====================================================================
!====================================================================
!>ovrlp_subspace
!
! Choose initial basis vectors based upon greatest overlap with HF
!----------------------------------------------------------------------
  subroutine ovrlp_subspace( m, initgdim, moints1, moints1len, moints2, &
    moints2len, pstring, pstep, plocate, qstring, qstep, qlocate, cidim,&
    pdets, qdets, pdetslen, qdetslen, adets, bdets, aelec, belec,      &
    orbitals, diagonals, out_vectors )

    implicit none
  ! ...integer input scalars...
    integer, intent(in) :: m, initgdim, moints1len, moints2len, cidim, &
                           pdetslen, qdetslen, adets, bdets, aelec, belec, &
                           orbitals
  ! ...real*8 input arrays...
    real*8,  dimension( moints1len ), intent(in) :: moints1  !1-e integrals
    real*8,  dimension( moints2len ), intent(in) :: moints2  !2-e integrals
  ! ...integer input arrays...
    integer, dimension( cidim, 2 ),   intent(in) :: pstring, qstring  ! Det. string pairs
    integer, dimension( pdetslen ),   intent(in) :: pdets, pstep, plocate
    integer, dimension( qdetslen ),   intent(in) :: qdets, qstep, qlocate
  ! ...real*8 input arrays... 
    real*8,  dimension( cidim ),      intent(in) :: diagonals
  ! ...real*8 OUTPUT arrays...
    real*8,  dimension( cidim, initgdim ), intent(out) :: out_vectors
  
  ! ...loop control scalars...
    integer :: i, j, k, l
  ! ...real*8 scalars...
    real*8, external :: ddot
  ! ...real*8 arrays...
    real*8, dimension(:,:),  allocatable     :: sub_block, eigvectors
    real*8, dimension(:,:),  allocatable     :: unitguess, hv_vectors, bs_vectors
    real*8, dimension(:),    allocatable     :: hartreefock, hv_hartree
  !--------------------------------------------------------------------
  ! Allocate arrays
    allocate(  sub_block( m, m ) )
    allocate( eigvectors( m, m ) )
    allocate(  unitguess( cidim, m ))
    allocate( hv_vectors( cidim, m) )
    allocate( bs_vectors( cidim, m) )
    allocate( hartreefock( cidim ) )
    allocate(  hv_hartree( cidim ) )

  ! Select m greatest elements of H on HF
    hartreefock    = 0d0
    hartreefock(1) = 1d0
    call acthv( hartreefock, moints1, moints2, moints1len, moints2len, pstring, &
                pstep, plocate, qstring, qstep, qlocate, cidim, pdets, qdets,   &
                pdetslen, qdetslen, adets, bdets, aelec, belec, orbitals,       &
                diagonals, hv_hartree )
  
  ! Construct vector guesses
    unitguess = 0d0
    call lowdiags( diagonals, cidim, m, unitguess )
  
  ! Perform Hv, generating hv_vectors
    do i=1, m
      call acthv( unitguess(1,i), moints1, moints2, moints1len, moints2len, &
                  pstring, pstep, plocate, qstring, qstep, qlocate, cidim,  &
                  pdets, qdets, pdetslen, qdetslen, adets, bdets, aelec, belec, &
                  orbitals, diagonals, hv_vectors(1,i) )
    end do

  ! Construct upper triangle of v.Hv matrix
    do i=1, m
      do j=1, m
        sub_block(j,i) = ddot( cidim, unitguess(1,i), 1, hv_vectors(1,i), 1 )
      end do
    end do

  ! Diagonalize this matrix
    call diag_hamsub( sub_block, m, m, eigvectors )
  
  ! Generate eigenvector guesses
    call ritz_vec( eigvectors, cidim, m, m, unitguess, bs_vectors )

  ! Generate output vectors
    do i=1, initgdim
      do j=1, cidim
        out_vectors(j,i) = bs_vectors(j,i)
      end do
    end do

  ! Deallocate arrays
    deallocate( sub_block, eigvectors, unitguess, hv_vectors, bs_vectors, &
                hartreefock, hv_hartree )

    return
  end subroutine
!======================================================================
!======================================================================
end module 
