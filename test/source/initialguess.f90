module initialguess
  !=============================================================================
  ! MODULE: initialguess
  ! --------------------
  !
  ! Purpose: perform inital guess routines for jayci.x
  !
  ! Contains:
  !  diagrefblock
  !
  ! By Christopher L Malbon
  ! Dept. of Chemistry, The Johns Hopkins University
  !-----------------------------------------------------------------------------
  implicit none
contains
  !*
  !*
  subroutine diagrefblock(cidim, refdim, moints1, m1len, moints2, m2len, &
    detlist, aelec, belec, orbs, ivdim, init_vecs, tot_frz)
    !===========================================================================
    ! prediagrefblock
    ! ---------------
    ! Purpose: diagonalize of reference block of size (refdim x refdim)
    !
    ! Input:
    !  cidim   = ci expansion size
    !  refdim  = dimension of subblock
    !  moints1 = 1-e integrals
    !  moints2 = 2-e integrals
    !  m1len   = len(moints1)
    !  m2len   = len(moints2)
    !  detlist = list of determinants in expansion
    !  aelec   = alpha electrons
    !  belec   = beta  electrons
    !  orbs    = orbitals
    !  ivdim   = number of initial vectors
    !
    ! Output:
    !  init_vecs = initial vectors for ci
    !---------------------------------------------------------------------------
    use construct,  only: ham_element
    use david_util, only: diag_dsyevr
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: cidim, refdim, ivdim
    integer, intent(in) :: m1len, m2len
    integer, intent(in) :: aelec, belec, orbs
    real*8,  intent(in) :: tot_frz
    integer, dimension(cidim), intent(in) :: detlist
    real*8,  dimension(m1len), intent(in) :: moints1
    real*8,  dimension(m2len), intent(in) :: moints2

    ! .. OUTPUT arguments ..
    real*8, dimension(cidim, ivdim), intent(out) :: init_vecs

    ! .. LOCAL arrays ..
    real*8, dimension(:,:), allocatable :: hblock
    real*8, dimension(:,:), allocatable :: evecs, tmpvs
    real*8, dimension(:),   allocatable :: evals

    ! .. LOCAL scalars ..
    integer :: i, j
    integer :: ierr

    ! allocate arrays
    allocate(hblock(refdim, refdim))
    allocate(evecs(refdim, refdim))
    allocate(evals(refdim))
    allocate(tmpvs(cidim, refdim))

    ! construct refdim x refdim block of H
    do i = 1, refdim
            ! diagonal elements
            hblock(i,i) = ham_element(detlist(i), detlist(i), moints1, &
              m1len, moints2, m2len, aelec, belec, orbs)
            do j = (i + 1), refdim
                    ! off-diagonal elements
                    hblock(j,i) = ham_element(detlist(j), detlist(i), moints1, &
                      m1len, moints2, m2len, aelec, belec, orbs)
                    hblock(i,j) = hblock(j,i)
            end do
    end do
    
    ! diagonalize hblock using dsyevr
    call diag_dsyevr(hblock, refdim, evecs, evals)

    write (*,10) evals(1) + tot_frz
10  format(1x,/,'Lowest eigenvalue in reference space =', f15.8,/)
    
    ! construct temporary vectors
    tmpvs = 0d0
    do i = 1, refdim
            tmpvs(1:refdim, i) = evecs(1:refdim, i)
    end do

    ! generate new vectors
    call genvecs_dgemm(tmpvs, cidim, refdim, ivdim, init_vecs)

    ! deallocate arrays
    deallocate(hblock, tmpvs, evecs, evals)
    
    return

  contains
    !*
    subroutine genvecs_dgemm(in_vecs, ivlen, ivnum, ovnum, out_vecs)
      !=========================================================================
      ! genvecs_dgemm
      ! -------------
      ! Purpose: generate vectors using dgemm
      !
      ! Input:
      !  in_vecs = input vectors
      !  ivlen   = leading dimension of in_vecs
      !  ivnum   = number of input  vectors
      !  ovnum   = number of output vectors
      !
      ! Output:
      !  out_vecs = output vectors
      !-------------------------------------------------------------------------
      implicit none

      ! .. INPUT arguments ..
      integer, intent(in) :: ivlen, ivnum, ovnum
      real*8,  dimension(ivlen, ivnum), intent(in) :: in_vecs

      ! .. OUTPUT arguments ..
      real*8, dimension(ivlen, ovnum), intent(out) :: out_vecs

      ! .. LOCAL scalars ..
      character*1 :: transa = 'n', transb = 'n'
      real*8      :: alpha = 1d0,  beta = 1d0
      integer     :: a_rows, b_cols, a_cols
      integer     :: i
      
      ! .. LOCAL arrays..
      real*8, dimension(ivnum, ivnum) :: unitmat
      real*8, dimension(ivlen, ivnum) :: tmp_vec

      ! set unitmat
      unitmat = 0d0
      do i = 1, ivnum
              unitmat(i,i) = 1d0
      end do

      ! set dgemm variables and call dgemm.
      !  Perform: c = alpha * A * B
      !           c = tmp_vec, A = in_vecs, B = unitmat
      a_rows = ivlen
      b_cols = ivnum
      a_cols = ivnum

      call dgemm(transa, transb, a_rows, b_cols, a_cols, alpha, in_vecs, &
        ivlen, unitmat, ivnum, beta, tmp_vec, ivlen)

      do i = 1, ovnum
              out_vecs(1:ivlen, i) = tmp_vec(1:ivlen, i)
      end do

      return
    end subroutine genvecs_dgemm
    !*
  end subroutine diagrefblock
  !*
  !*
end module initialguess
