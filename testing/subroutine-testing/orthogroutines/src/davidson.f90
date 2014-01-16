! Davidson algorithm subroutines
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chemistry, The Johns Hopkins University
!
! Last edit: 12-12-13
!====================================================================
!====================================================================
!======================================================================
!======================================================================
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
  use david_util
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
  real*8,  dimension( roots ),           intent(inout) :: eigvalues
  real*8,  dimension( cidim, roots),     intent(inout) :: eigvectors
! ...loop integer scalars...
  integer :: i, j, k, l, m
! ...davidson integer scalars...
  integer :: space_dim, currnt_root
! ...davidson real*8  scalars...
  real*8            :: resid_norm, ddot
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
    if (allocated( sub_ham )) deallocate(sub_ham)
    allocate( sub_ham( krmin, krmin ) )
 
    call build_subham( bs_vectors, hv_vectors, cidim, krmin, sub_ham )
    print *, "  "
    print *, " ========================================== "
    print *, "  Printing subspace Hamiltonian by column..."
    do k=1, krmin
      print *, " --------------------- "
      do j=1, krmin
        print *, sub_ham(j,k)
      end do
    end do
 
  ! Diagonalize this matrix with DSYEVR
    if (allocated( kry_eigvec )) deallocate( kry_eigvec )
    if (allocated( kry_eigval )) deallocate( kry_eigval )
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
      resid_norm = sqrt( ddot( cidim, residual, 1, residual, 1 ) )
      print *, " ================== "
      print *, " Eigenvalue:     ", kry_eigval(currnt_root)
      print *, " Residual norm:  ", resid_norm
      print *, "  "

    ! Test if ||r|| < rtol
      if ( resid_norm  < rtol ) then
        exit main_loop
      end if

    ! Generate new vector
      call gen_newvector( residual, diagonals, kry_eigval(currnt_root),&
                          cidim, new_vector )
      print *, " ============================================ "
      print *, "  Printing new vector... "
      do k=1, cidim
        print *, new_vector(k)
      end do
      print *, " "
     
    ! Orthogonalize new vector to basis space 
      call orthog_newvec( bs_vectors, new_vector, cidim, space_dim )
    ! Print orthogonalization information
      print *, " ============================================ "
      print *, "  Printing orthogonalization information... "
      do k=1, space_dim
        print *, " Overlap of vector ", k, " with new vector: "
        print *, ddot( cidim, bs_vectors(1,k), 1, new_vector, 1 )
      end do
    ! Add new vector to basis. *Note: This subroutine makes space_dim = space_dim+1
      call append_newvec( bs_vectors, space_dim, cidim, new_vector )
      print *, " "
      print *, " Dimension of space is now  ", space_dim
    ! Perform Hv on new vector
      call acthv( bs_vectors(1,space_dim), moints1, moints2, moints1len, moints2len,    &
                  pstring, pstep, plocate, qstring, qstep, qlocate, cidim, &
                  pdets, qdets, pdetslen, qdetslen, adets, bdets, aelec, belec,&
                  orbitals, diagonals, hv_vectors(1,space_dim) )

    ! Build vHv subspace Hamiltonian
      if (allocated( sub_ham )) deallocate(sub_ham)
      allocate( sub_ham( space_dim, space_dim ) )
      
      call build_subham( bs_vectors, hv_vectors, cidim, space_dim, sub_ham )
      
      print *, "  "
      print *, " ============================================== "
      print *, "  Printing subspace hamiltonian... "
      do k=1, space_dim
        print *, " ------------------------ "
        do l=1, space_dim
          print *, sub_ham(l,k)
        end do
      end do
   

    ! Diagonalize this matrix with DSYEVR
      if (allocated( kry_eigvec )) deallocate(kry_eigvec)
      if (allocated( kry_eigval )) deallocate(kry_eigval)
      allocate( kry_eigvec( space_dim, space_dim ) )
      allocate( kry_eigval( space_dim ) )

      call diagonalize( sub_ham, space_dim, kry_eigvec, kry_eigval )

    ! Deallocate sub_ham <-- it has been destroyed by dsyevr()
      deallocate( sub_ham )

    ! Test size of subspace
      if ( space_dim .eq. krmax ) then

        call space_trunc( bs_vectors, cidim, krmin, krmax, currnt_root, kry_eigvec )
        space_dim = krmin

        ! Zero out hv_vectors
        hv_vectors = 0d0

        ! Exit sub loop
        exit sub_loop

      end if
    end do sub_loop
  end do main_loop
  return
   
end subroutine
!======================================================================
!======================================================================

  







