module david_util
  implicit none
contains
  !====================================================================
  !====================================================================
  !>David_Initial
  !
  ! Subroutine to perform Hv on inital basis space
  !--------------------------------------------------------------------
  subroutine David_Initial( vectors1, num_vecs, dim_vecs, moints1, moints1len, &
    moints2, moints2len, pstring, pstep, plocate, qstring, qstep, qlocate,   &
    xreflist, pdets, qdets, pdetslen, qdetslen, adets, bdets, aelec, belec,  &
    orbitals, diagonals, nfrozen, ndocc, nactive, vectors2 )
    use construct, only: exp_construct
    implicit none
    integer, intent(in) :: num_vecs, dim_vecs, moints1len, moints2len, pdetslen, &
                           qdetslen, adets, bdets, aelec, belec, orbitals, &
                          nfrozen, ndocc, nactive
    integer, dimension(pdetslen),    intent(in) :: pdets, plocate, pstep
    integer, dimension(qdetslen),    intent(in) :: qdets, qlocate, qstep
    integer, dimension(dim_vecs),  intent(in) :: pstring, qstring
    integer, dimension(dim_vecs),    intent(in) :: xreflist
    real*8,  dimension(dim_vecs),    intent(in) :: diagonals
    real*8,  dimension(moints1len),  intent(in) :: moints1
    real*8,  dimension(moints2len),  intent(in) :: moints2
    real*8,  dimension(dim_vecs, num_vecs), intent(inout) :: vectors2
    real*8,  dimension(dim_vecs, num_vecs), intent(in)    :: vectors1
    integer :: i
  ! Loop over initial vectors
    do i=1, num_vecs
            print *, "    Vector 1"
            CALL acthv (vectors1(1,i),moints1,moints2,moints1len,moints2len,&
                  pstring,pstep,plocate,qstring,qstep,qlocate,xreflist,dim_vecs,&
                  pdets,pdetslen,qdets,qdetslen,adets,bdets,aelec,belec,&
                  orbitals,diagonals,nfrozen,ndocc,nactive,vectors2(1,i))
    end do
    RETURN
  end subroutine
  !====================================================================
  !====================================================================
  !>build_subham
  !
  ! Subroutine to build subspace hamiltonian
  !--------------------------------------------------------------------
  subroutine build_subham( bs_vectors, hv_vectors, dim_ci, dim_bs, matrix )
    implicit none
    integer, intent(in) :: dim_ci, dim_bs
    real*8, dimension( dim_ci, dim_bs ), intent(in)  :: bs_vectors, hv_vectors
    real*8, dimension( dim_bs, dim_bs ), intent(out) :: matrix
    integer :: i, j
    real*8  :: ddot
  !--------------------------------------------------------------------
  ! Loop over dimension of basis
    matrix = 0d0
    do i=1, dim_bs
      do j=1, dim_bs
        matrix(j,i) = ddot( dim_ci, bs_vectors(1,j), 1, hv_vectors(1,i), 1 )
      end do
    end do
    return
  end subroutine
  !====================================================================
  !====================================================================
  !>diag_dsyevr
  !
  ! Subroutine to diagonalize subspace hamiltonian
  !--------------------------------------------------------------------
  subroutine diag_dsyevr( matrix, dim_mat, eig_vec, eig_val )
    implicit none
    integer, intent(in) :: dim_mat
    real*8, dimension( dim_mat, dim_mat ), intent(inout) :: matrix
    real*8, dimension( dim_mat, dim_mat ), intent(out)   :: eig_vec
    real*8, dimension( dim_mat ),          intent(out)   :: eig_val
  !...DSYEVR variables...
    character(len=1) :: uplo, jobz, rnge
    real*8           :: abstol, dlamch, vl, vu
    integer          :: lwork, liwork, il, iu, info, eigfound
    integer, dimension(:), allocatable     :: isuppz
    real*8,  dimension(:), allocatable     :: work
    real*8,  dimension(:), allocatable     :: iwork 
  !--------------------------------------------------------------------
    jobz = 'v'            ! Return eigenvectors & eigenvalues
    rnge = 'a'            ! Return full range of eigenvectors & eigenvalues
    uplo = 'l'            ! Lower  triangle is stored in matrix A
    abstol = dlamch( 'Safe minimum' ) 
    lwork  = 30*dim_mat + 10
    liwork = 15*dim_mat
    info = 0
  ! Allocate arrays
    allocate( isuppz( 2*dim_mat ) )
    allocate( work( lwork ) )
    allocate( iwork( liwork) )
    call dsyevr( jobz, rnge, uplo, dim_mat, matrix, dim_mat, vl, vu, il, &
                 iu, abstol, eigfound, eig_val, eig_vec, dim_mat, isuppz,&
                 work, lwork, iwork, liwork, info )
    if ( info .ne. 0 ) then
       print *, "Info = ", info
       stop " ***Diagonalization of v.Hv failed! "
    end if
    ! Deallocate arrays
    deallocate( isuppz, work, iwork )
    return
  end subroutine diag_dsyevr
  !====================================================================
  !====================================================================
  !>gen_ResidVec
  !
  ! Subroutine to generate residual
  !--------------------------------------------------------------------
  subroutine gen_ResidVec( bs_vec, hv_vec, length, sp_dim, kry_vec, kry_val,&
    root, residual )
    implicit none
    integer, intent(in) :: length, sp_dim, root
    real*8,  dimension( length, sp_dim ), intent(in) :: bs_vec, hv_vec
    real*8,  dimension( sp_dim, sp_dim ), intent(in) :: kry_vec
    real*8,  dimension( sp_dim ),         intent(in) :: kry_val
    real*8,  dimension( length ),      intent(out) :: residual
    integer :: i, j
  !--------------------------------------------------------------------
  ! Construct residual
    residual = 0d0
    do i=1, length
      do j=1, sp_dim
        residual(i) = residual(i) + kry_vec(j,root)*&
                      ( hv_vec(i,j) - kry_val(root)*bs_vec(i,j) )
      end do
    end do
    print *, "Returning from gen_ResidVec"
    return
  end subroutine gen_ResidVec
  !====================================================================
  !====================================================================
  !>gen_newVector
  !
  ! Subroutine to generate new vector 
  !--------------------------------------------------------------------
  subroutine gen_newVector(residual, diagonals, currnt_eigval,&
    cidim, new_vector )
    implicit none
    integer, intent(in) :: cidim
    real*8,  intent(in) :: currnt_eigval
    real*8,  dimension( cidim ), intent(in)    :: residual, diagonals
    real*8,  dimension( cidim ), intent(out) :: new_vector
    integer :: i
  !--------------------------------------------------------------------
    do i=1, cidim
      new_vector(i) = - residual(i) / ( diagonals(i) - currnt_eigval )
    end do
    return
  end subroutine gen_newVector
  !====================================================================
  !====================================================================
  !>orthog_newvec
  !
  ! Subroutine to orthogonalize the new vector to the space
  !--------------------------------------------------------------------
  subroutine orthog_newvec( bs_vectors, new_vector, length, space_dim )
    use orthogroutines, only: orthogvector
    implicit none
    integer, intent(in)    :: length
    integer, intent(in)    :: space_dim
    real*8, dimension( length),                  intent(inout) :: new_vector
    real*8, dimension( length, space_dim ),      intent(in) :: bs_vectors
    integer :: i, orth_iter
    real*8, dimension(:), allocatable :: ovrlp_test
    real*8  :: ddot
    real*8, parameter :: max_ovrlp = 1.0D-15
    real*8  :: norm
  !--------------------------------------------------------------------
  ! Allocate arrays
    allocate( ovrlp_test( space_dim ) )
  ! Orthogonalize
    call orthogvector( bs_vectors, space_dim, space_dim, length, new_vector )
  ! Test orthogonality
    do i=1, space_dim
      ovrlp_test(i) = ddot( length, bs_vectors(1,i), 1, new_vector, 1 )
    end do
    orth_iter = 1
    do while ( maxval( ovrlp_test ) > max_ovrlp .and. orth_iter .le. 10 )
      call orthogvector( bs_vectors, space_dim, space_dim, length, new_vector )
      do i=1, space_dim
        ovrlp_test(i) = ddot( length, bs_vectors(1,i), 1, new_vector, 1 )
      end do
      orth_iter = orth_iter + 1
    end do
  ! Deallocate arrays
    deallocate( ovrlp_test )
    return
  end subroutine orthog_newvec
  !====================================================================
  !====================================================================
  !>append_newvec
  !
  ! Subroutine to append vector to space
  !--------------------------------------------------------------------
  subroutine append_newvec( basis, space_dim, length, new_vector )
    implicit none
    integer, intent(inout) :: space_dim
    integer, intent(in)    :: length
    real*8,  dimension(length), intent(in)                :: new_vector
    real*8,  dimension(length,space_dim+1), intent(inout) :: basis
    integer :: i
  !--------------------------------------------------------------------
    do i=1, length
      basis(i,space_dim+1) = new_vector(i)
    end do
    space_dim = space_dim+1
    return
  end subroutine append_newvec
  !====================================================================
  !====================================================================
  !>space_trunc
  !
  ! Subroutine to truncate subspace
  !--------------------------------------------------------------------
  subroutine space_trunc( basis_vec, len_basis, kry_min, kry_max, root, &
    kry_vecs )
    use orthogroutines
    implicit none
    integer, intent(in) :: len_basis, kry_min, kry_max, root
    real*8, dimension( len_basis, kry_max ), intent(inout) :: basis_vec
    real*8, dimension( kry_max,   kry_max ), intent(in)    :: kry_vecs
    integer :: i, j
    real*8  :: ddot
  !...DGEMM variables...
    character(len=1) :: transa, transb
    integer          :: m, n, k
    integer          :: lda, ldb, ldc
    real*8           :: alpha, beta
    real*8, dimension(:,:), allocatable :: c
  !--------------------------------------------------------------------
  ! Allocate arrays
    allocate( c( len_basis, kry_max ) )
  ! B_ij e_jk = A_ik
  ! Set DGEMM variables
    transa = 'n'
    transb = 'n'
    m = len_basis
    n = kry_max
    k = kry_max
    alpha = 1d0
    beta  = 0d0
    lda = m
    ldb = k
    ldc = m
  
  ! Call DGEMM
    call dgemm( transa, transb, m, n, k, alpha, basis_vec, lda, kry_vecs, &
                ldb, beta, c, ldc )
  

  ! Zero out basis space
    basis_vec = 0d0

  ! Make c = basis_vec
    do i = 1, kry_min
      do j = 1, len_basis
        basis_vec(j,i) = c(j,i)
      end do
    end do
  ! 
  ! Deallocate c
    deallocate( c )
    return
  end subroutine space_trunc
  !====================================================================
  !====================================================================
end module

