! Subroutine to compute diagonal contribution <K|H|K>
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. Chem., The Johns Hopkins University
!
! Last edit: 11-19-13
!====================================================================
!Input:
! pstring   = determinant pairings by alpha strings       integer array
! pstep     = alpha string step array                     integer array
! qstring   = determinant pairings by beta strings        integer array
! qstep     = beta string step array                      integer array
! cidim     = dimension of space                          integer scalar
! totels    = total electrons                             integer scalar
! pdets     = truncated alpha strings
! qdets     = truncated beta strings
! orbitals  = number of MO's                              integer scalar
! aelec     = alpha electrons                             integer scalar
! belec     = beta electrons                              integer scalar
! adets     = alpha determinants ( non-truncated )        integer scalar
! bdets     = beta determinants ( non-truncated )         integer scalar
! moints1   = 1-e integrals                               real*8  array
! moints1len= length of moints1()                         integer scalar
! moints2   = 2-e integrals                               real*8  array
! moints2len= length of moints2()                         integer scalar
!Output:
! diagonals = < K | H | K >                               real*8  array
!====================================================================
subroutine diagonal( pstring, pstep, psteplen, qstring, qstep, qsteplen,   & 
  pdets, pdetslen, qdets, qdetslen, detlist, cidim, totels, orbitals, aelec, belec, &
  adets, bdets, moints1, moints1len, moints2, moints2len, plocate, qlocate, diagonals )
  use construct, only: ham_element_diag
  implicit none
! ...input integer scalars...
  integer,intent(in) :: psteplen, qsteplen, cidim, totels, orbitals, aelec, &
                        belec, adets, bdets, moints1len, moints2len, pdetslen,&
                        qdetslen
! ...input integer arrays...
  integer, dimension(cidim,2),  intent(in) :: pstring, qstring
  integer, dimension(cidim),    intent(in) :: detlist
  integer, dimension(psteplen), intent(in) :: pstep, plocate
  integer, dimension(qsteplen), intent(in) :: qstep, qlocate
  integer, dimension(pdetslen), intent(in) :: pdets
  integer, dimension(qdetslen), intent(in) :: qdets
! ...input real*8 arrays...
  real*8, dimension(moints1len), intent(in) :: moints1
  real*8, dimension(moints2len), intent(in) :: moints2
! ...output real*8 arrays...
  real*8, dimension(cidim), intent(inout) :: diagonals
! ...integer scalars...
  integer :: i
! ...integer arrays...
!  integer, dimension(cidim ) :: alphadiagindex, betadiagindex
!--------------------------------------------------------------------
  diagonals = 0d0
  do i=1, cidim
    diagonals(i) = ham_element_diag( detlist(i), moints1, moints1len, moints2, moints2len, &
                                     aelec, belec, orbitals )
  end do
  
  return
end subroutine
!====================================================================
!====================================================================
          
    
