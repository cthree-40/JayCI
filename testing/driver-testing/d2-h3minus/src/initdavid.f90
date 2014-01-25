! Subroutine 1 of Davidson Algorithm
! This subroutine is responsible for computing Hv of the initial
!  vector space, or truncated space
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chemistry, The Johns Hopkins University
!
! Last edit: 10-22-2013 - rewritten from Det-CI-4.0.0
!====================================================================
!Input:
! basisvecs   = basis vectors                   real*8  array  2-d
! basislda    = LDA of basisvecs                integer scalar
! basisdim    = number of basis vectors         integer scalar
! moints1     = 1-e integrals                   real*8  array  1-d
! moints2     = 2-e integrals                   real*8  array  1-d
! max1e       = length of moints1               integer scalar
! max2e       = length of moints2               integer scalar
! actdets     = active determinant list         integer array  1-d
! aelec       = alpha electrons                 integer scalar 
! belec       = beta electrons                  integer scalar
! totels      = aelec + belec                   integer scalar
! orbitals    = number of MO's                  integer scalar
! adets       = alpha determinants              integer scalar
! bdets       = beta determinants               integer scalar
! diagonals   = < K | H | K >                   real*8  array  1-d
! nfrzn       = number of frozen orbitals       integer scalar
! ndocc       = number of docc orbitals         integer scalar
! ncas        = number of active orbitals       integer scalar
! krmax       = maximum dimension of Kry space  integer scalar
!Output:
! hvectors    = vectors from Hv on basis        real*8  array   2-d
!===================================================================
subroutine initdavid( basisvecs, basislda, basisdim, moints1, moints2,  &
  max1e, max2e, actdets, aelec, belec, totels, orbitals, adets, bdets,  &
  diagonals, nfrzn, ndocc, ncas, krmax, hvectors )

  implicit none

! ...input integer scalars...
  integer, intent(in) :: basislda, basisdim, max1e, max2e, aelec, &
                         belec, totels, orbitals, adets, bdets,   &
                         nfrzn, ndocc, ncas, krmax

! ...input integer arrays...
  integer, dimension(basislda) :: actdets

! ...input real*8 arrays...
  real*8, dimension(max1e) :: moints1
  real*8, dimension(max2e) :: moints2
  real*8, dimension(basislda, basisdim) :: basisvecs
  real*8, dimension(basislda) :: diagonals

! ...output real*8 arrays...
  real*8, dimension(basislda, basisdim) :: hvectors

! ...loop integer scalars...
  integer :: i, j, k
!--------------------------------------------------------------------

! Zero out hvectors
  hvectors = 0d0

! Loop over basisvecs
  do i=1, basisdim
    call acthv( basisvecs(1,i), moints1, moints2, max1e, max2e,    &
                actdets, basislda, aelec, belec, totels, orbitals, &
                adets, bdets, nfrzn, ndocc, ncas, diagonals, hvectors(1,i))
  end do

  return

end subroutine




