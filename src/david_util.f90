module david_util
  implicit none
contains
  !====================================================================
  !====================================================================
  !>dvd_initial
  !
  ! Subroutine to perform Hv on inital basis space
  !--------------------------------------------------------------------
  subroutine dvd_initial( vectors1, num_vecs, dim_vecs, moints1, moints1len, &
    moints2, moints2len, pstring, pstep, plocate, qstring, qstep, qlocate,   &
    pdets, qdets, pdetslen, qdetslen, adets, bdets, aelec, belec, orbitals,  &
    diagonals, vectors2 )
    implicit none
    integer, intent(in) :: num_vecs, dim_vecs, moints1len, moints2len, pdetslen, &
                           qdetslen, adets, bdets, aelec, belec, orbitals
    integer, dimension(pdetslen),  intent(in) :: pdets, plocate, pstep
    integer, dimension(qdetslen),  intent(in) :: qdets, qlocate, qstep
    integer, dimension(dim_vecs),  intent(in) :: pstring, qstring
    real*8,  dimension(dim_vecs),  intent(in) :: diagonals
    real*8,  dimension(moints1len),intent(in) :: moints1
    real*8,  dimension(moints2len),intent(in) :: moints2
    real*8,  dimension(dim_vecs, num_vecs), intent(inout) :: vectors2
    real*8,  dimension(dim_vecs, num_vecs), intent(in)    :: vectors1
    integer :: i
  !--------------------------------------------------------------------
  ! Loop over initial vectors
    do i=1, num_vecs
      call acthv( vectors1(1,i), moints1, moints2, moints1len, moints2len,    &
                  pstring, pstep, plocate, qstring, qstep, qlocate, dim_vecs, &
                  pdets, qdets, pdetslen, qdetslen, adets, bdets, aelec, belec,&
                  orbitals, diagonals, vectors2(1,i) )
    end do
    return
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
  !>diagonalize
  !
  ! Subroutine to diagonalize subspace hamiltonian
  !--------------------------------------------------------------------
  subroutine diagonalize( matrix, dim_mat, eig_vec, eig_val )
    implicit none
    integer, intent(in) :: dim_mat
    real*8, dimension( dim_mat, dim_mat ), intent(inout) :: matrix
    real*8, dimension( dim_mat, dim_mat ), intent(out)   :: eig_vec
    real*8, dimension( dim_mat ),          intent(out)   :: eig_val
  !...DSYEVR variables...
    character(len=1) :: uplo, jobz, rnge
    real*8           :: abstol, dlamch, vl, vu
    integer          :: lwork, liwork, il, iu, info, eigfound
    integer, dimension( dim_mat )          :: isuppz
    real*8,  dimension( 26*dim_mat + 10  ) :: work
    real*8,  dimension( 11*dim_mat )       :: iwork 
  !--------------------------------------------------------------------
    jobz = 'v'            ! Return eigenvectors & eigenvalues
    rnge = 'a'            ! Return full range of eigenvectors & eigenvalues
    uplo = 'l'            ! Lower  triangle is stored in matrix A
    abstol = dlamch( 'Safe minimum' ) 
    lwork  = 26*dim_mat + 10
    liwork = 11*dim_mat
    call dsyevr( jobz, rnge, uplo, dim_mat, matrix, dim_mat, vl, vu, il, &
                 iu, abstol, eigfound, eig_val, eig_vec, dim_mat, isuppz,&
                 work, lwork, iwork, liwork, info )
    if ( info .ne. 0 ) stop " ***Diagonalization of v.Hv failed! "
    return
  end subroutine
  !====================================================================
  !====================================================================
  !>gen_residual
  !
  ! Subroutine to generate residual
  !--------------------------------------------------------------------
  subroutine gen_residual( bs_vec, hv_vec, length, sp_dim, kry_vec, kry_val,&
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
    do i=1, sp_dim
      do j=1, length
        residual(j) = kry_vec(i,root)*&
                      ( hv_vec(j,i) - kry_val(root)*bs_vec(j,i) )
      end do
    end do
    return
  end subroutine
  !====================================================================
  !====================================================================
  !>gen_newvector
  !
  ! Subroutine to generate new vector 
  !--------------------------------------------------------------------
  subroutine gen_newvector(residual, diagonals, currnt_eigval,&
    cidim, new_vector )
    implicit none
    integer, intent(in) :: cidim
    real*8,  intent(in) :: currnt_eigval
    real*8,  dimension( cidim ), intent(in)    :: residual, diagonals
    real*8,  dimension( cidim ), intent(out) :: new_vector
    integer :: i
  !--------------------------------------------------------------------
    do i=1, cidim
      new_vector(i) = - residual(i) / ( diagonals(i) - currnt_eigval + 1D-6 )
    end do
    return
  end subroutine
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
    real*8, dimension( space_dim ) :: ovrlp_test
    real*8  :: ddot
    real*8, parameter :: max_ovrlp = 1.0D-15
    real*8  :: norm
  !--------------------------------------------------------------------
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
    return
  end subroutine
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
    real*8,  dimension(length)                            :: new_vector
    real*8,  dimension(length,space_dim+1), intent(inout) :: basis
    integer :: i
  !--------------------------------------------------------------------
    do i=1, length
      basis(i,space_dim+1) = new_vector(i)
    end do
    space_dim = space_dim+1
    return
  end subroutine
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
    real*8, dimension( len_basis, kry_max ) :: c
  !--------------------------------------------------------------------
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
    call modgramschmidt( basis_vec, kry_min, kry_min, len_basis )

    return
  end subroutine
  !====================================================================
  !====================================================================
end module

