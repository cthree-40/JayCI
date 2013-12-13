! Davidson algorithm subroutines
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chemistry, The Johns Hopkins University
!
! Last edit: 12-12-13
!====================================================================
!====================================================================
!>davidson
!
! Main subroutine
!--------------------------------------------------------------------
subroutine davidson( initgdim, init_vectors, diagonals, moints1, moints1len,    & 
  moints2, moints2len, cidim, pstring, pstep, plocate, qstring, qstep, qlocate, &
  pdets, pdetslen, qdets, qdetslen, adets, bdets, aelec, belec, orbitals, krmin,&
  krmax, rtol, roots, max_iter, eigvalues, eigvectors )
!Input:
! initgdim     = number of initial vector guesses                      integer scalar
! init_vectors = initial vector guess                                  real*8  array  2-d
! diagonals    = <K|H|K>                                               real*8  array  1-d
! moints1      = 1-e integrals                                         real*8  array  1-d
! moints1len   = length of moints1                                     integer scalar
! moints2      = 2-e integrals                                         real*8  array  1-d
! moints2len   = length of moints2                                     integer scalar
! cidim        = dimension of CI                                       integer scalar
! pstring      = array of pstrings and corresponding qstrings          integer array  2-d
! pstep        = number of qstrings corresponding to ith pstring       integer array  1-d
! plocate      = location of ith p string in pstring                   integer array  1-d
! qstring      
! qstep
! qlocate
! pdets        = active p strings                                      integer array  1-d
! qdets        = active q strings                                      integer array  1-d
! pdetslen     = number of active p strings                            integer scalar
! qdetslen     = number of active q strings                            integer scalar
! adets        = untruncated number of p strings                       integer scalar
! bdets        = untruncated number of q strings                       integer scalar
! aelec        = alpha electrons                                       integer scalar
! belec        = beta electrons                                        integer scalar
! orbitals     = number of orbitals                                    integer scalar
! krmin        = lowest dim of Krylov space                            integer scalar
! krmax        = greatest dimension of Krylov space                    integer scalar
! rtol         = convergence tolerance for CI                          real*8  scalar
! roots        = number of roots to find                               integer scalar
! max_iter     = maximum iterations
!Output:
! eigvalues    = eigenvalues                                           real*8  array  1-d
! eigvectors   = eigenvectors                                          real*8  array  2-d
!--------------------------------------------------------------------
  implicit none
! ...input integer scalars...
  integer, intent(in) :: initgdim, moints1len, moints2len, cidim, pdetslen, qdetslen, &
                         adets, bdets, aelec, belec, orbitals, krmin, krmax, roots, max_iter
! ...input real*8 scalars...
  real*8,  intent(in) :: rtol
! ...input real*8 arrays...
  real*8,  dimension( cidim ),           intent(in) :: diagonals
  real*8,  dimension( cidim, initgdim ), intent(in) :: init_vectors
  real*8,  dimension( moints1len ),      intent(in) :: moints1
  real*8,  dimension( moints2len ),      intent(in) :: moints2
! ...input integer arrays...
  integer, dimension( cidim, 2 ),        intent(in) :: pstring, qstring
  integer, dimension( pdetslen ),        intent(in) :: pstep, plocate, pdets
  integer, dimension( qdetslen ),        intent(in) :: qstep, qlocate, qdets
! ...OUTPUT real*8 arrays...
  real*8,  dimension( roots ),           intent(inout) :: eigenvalues
  real*8,  dimension( cidim, roots),     intent(inout) :: eigenvectors
! ...loop integer scalars...
  integer :: i, j, k, l, m
! ...davidson integer scalars...
  integer :: space_dim, currnt_root
! ...davidson real*8  scalars...
  real*8            :: resid_norm
  real*8, parameter :: max_ovrlp = 1.0D-16
! ...davidson real*8 arrays...
  real*8, dimension( cidim, krmax )   :: bs_vectors, hv_vectors
  real*8, dimension( cidim )          :: residual, new_vector
  real*8, dimension(:,:), allocatable :: sub_ham
  real*8, dimension(:,:), allocatable :: kry_eigvec
  real*8, dimension(:),   allocatable :: kry_eigval
!--------------------------------------------------------------------
! Make currnt_root = 1, for debugging and testing - CLM 12-13-13
  currnt_root = 1

! Zero out bs_vectors and hv_vectors
  bs_vectors = 0d0
  hv_vectors = 0d0
! Equate bs_vectors to the inital vector guess
  do i=1, krmin
    do j=1, cidim
      bs_vectors(j,i) = init_vectors(j,i)
    end do
  end do
!
! Main Loop of Davidson algorithm
  main_loop: do i=1, int( max_iter / (krmax-krmin) )
  ! Perform Hv on initial space to generate vectors v_i = Hb_i
 
    call dvd_initial( bs_vectors, krmin, cidim, moints1, moints1len, moints2, &
                      moints2len, pstring, pstep, plocate, qstring, qstep,    &
                      qlocate, pdets, qdets, pdetslen, qdetslen, adets, bdets,&
                      aelec, belec, orbitals, diagonals, hv_vectors )
 
  ! Construct subspace hamiltonian (v_i . H v_j ) = H_ij
    if allocated( sub_ham ) deallocate(sub_ham)
    allocate( sub_ham( krmin, krmin ) )
 
    call build_subham( bs_vectors, hv_vectors, cidim, krmin, sub_ham )
 
  ! Diagonalize this matrix with DSYEVR
    if allocated( kry_eigvec ) deallocate( kry_eigvec )
    if allocated( kry_eigval ) deallocate( kry_eigval )
    allocate( kry_eigvec( krmin, krmin ) )
    allocate( kry_eigval( krmin ) )
 
    call diagonalize( sub_ham, krmin, kry_eigvec, kry_eigval )
 
  ! Deallocate sub_ham <-- it has been destroyed by dsyevr()
    deallocate(sub_ham)
    
  ! Set space dimension, space_dim = krmin
    space_dim = krmin

  ! Enter sub_loop
    sub_loop: do j=1, krmax+1
   
    ! Form residual vector
      call gen_residual( bs_vectors, hv_vectors, cidim, space_dim, kry_eigvec, &
                         kry_eigval, currnt_root, residual )  

    ! Compute norm of residual
      resid_norm = sqrt( ddot( cidim, residual, 1, residual, 1 )

    ! Test if ||r|| < rtol
      if ( resid_norm  < rtol ) then
        exit main_loop
      end if

    ! Generate new vector
      call gen_newvector( residual, diagonals, kry_eigval(currnt_root), space_dim,&
                          cidim, new_vector )
     
    ! Orthogonalize new vector to basis space 










    
end subroutine
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
  real*8,  dimension(dim_vecs, num_vecs), intent(inout) :: num_vecs
  integer :: i
!--------------------------------------------------------------------
! Loop over initial vectors
  do i=1, num_vecs
    call acthv( vectors1(1,i), moints1, moints2, moints1len, moints2len,    &
                pstring, pstep, plocate, qstring, qstep, qlocate, dim_vecs, &
                pdets, qdets, pdetslen, qdetslen, adets, bets, aelec, belec,&
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
  real*8, dimension( dim_ci, dim_bs ), intent(in)    :: bs_vectors, hv_vectors
  real*8, dimension( dim_bs, dim_bs ), intent(inout) :: matrix
  integer :: i, j
  real*8  :: ddot
!--------------------------------------------------------------------
! Loop over dimension of basis
  matrix = 0d0
  do i=1, dim_bs
    do j=1, dim_bs
      matrix(j,i) = ddot( dim_ci, bs_vectors, 1, hv_vectors, 1 )
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
  real*8, dimension( dim_mat, dim_mat ), intent(in)    :: matrix
  real*8, dimension( dim_mat, dim_mat ), intent(inout) :: eig_vec
  real*8, dimension( dim_mat ),          intent(inout) :: eig_val
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
  if ( info .ne. 0 ) stop " Diagonalization of v.Hv failed. Aborting... info=", info
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
  real*8,  dimension( length ),      intent(inout) :: residual
  integer :: i, j
!--------------------------------------------------------------------
! Construct residual
  residual = 0d0
  do i=1, sp_dim
    do j=1, length
      residual(j) = residual(j) + kry_vec(i,root)*&
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
subroutine gen_newvector(residual, diagonals, currnt_eigval, space_dim,&
  cidim, new_vector )
  implicit none
  integer, intent(in) :: space_dim, currnt_root, cidim
  real*8,  intent(in) :: currnt_eigval
  real*8,  dimension( cidim ), intent(in)    :: residual, diagonals
  real*8,  dimension( cidim ), intent(inout) :: new_vector
  integer :: i
!--------------------------------------------------------------------
  new_vector = 0d0
  do i=1, cidim
    new_vector(i) = - residual(i) / ( diagonals(i) - currnt_eigval + 1D-9 )
  end do
  return
end subroutine
