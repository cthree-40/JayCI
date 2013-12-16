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
  adets, bdets, moints1, moints1len, moints2, moints2len, pxreflist, qxreflist, &
  xreflistlen, plocate, qlocate, diagonals )

  implicit none
! ...input integer scalars...
  integer,intent(in) :: psteplen, qsteplen, cidim, totels, orbitals, aelec, &
                        belec, adets, bdets, moints1len, moints2len, pdetslen,&
                        qdetslen, xreflistlen
! ...input integer arrays...
  integer, dimension(cidim,2),  intent(in) :: pstring, qstring
  integer, dimension(cidim),    intent(in) :: detlist
  integer, dimension(psteplen), intent(in) :: pstep, plocate
  integer, dimension(qsteplen), intent(in) :: qstep, qlocate
  integer, dimension(pdetslen), intent(in) :: pdets
  integer, dimension(qdetslen), intent(in) :: qdets
  integer, dimension(xreflistlen),intent(in)::pxreflist,qxreflist
! ...input real*8 arrays...
  integer, dimension(moints1len), intent(in) :: moints1
  integer, dimension(moints2len), intent(in) :: moints2
! ...output real*8 arrays...
  real*8, dimension(cidim), intent(inout) :: diagonals
! ...integer arrays...
!  integer, dimension(cidim ) :: alphadiagindex, betadiagindex
!--------------------------------------------------------------------
  diagonals = 0d0
! Alpha string contribution
  call alpha_diag( pdets, pdetslen, pstring, cidim, pstep, plocate, &
                   moints1, moints1len, &
                   moints2, moints2len, adets, bdets, orbitals, aelec,   &
                   belec, diagonals )

!====================================================================
module diag_util 
  implicit none
contains
!====================================================================
!====================================================================
!>alpha_diag
!
! Subroutine to compute diagonal contribution of alpha strings
!--------------------------------------------------------------------
  subroutine alpha_diag( pdets, pdetslen, pstring, cidim, pstep, plocate, &
    moints1, moints1len, moints2, moints2len, adets, bdets, orbitals,     &
    aelec, belec, diagonals )
    use detci2
    use detci5
    implicit none
    integer, intent(in) :: pdetslen, cidim, moints1len, moints2len, adets,&
                           bdets, orbitals, aelec, belec
    integer, dimension( pdetslen ), intent(in) :: pdets, plocate
    integer, dimension( cidim, 2 ), intent(in) :: pstring
    real*8,  dimension(moints1len), intent(in) :: moints1
    real*8,  dimension(moints2len), intent(in) :: moints2
    real*8,  dimension(cidim),   intent(inout) :: diagonals
    integer :: i, j
    integer, dimension( aelec ) :: astring1
    real*8  :: int1e1, int2e1, int2e2
  !------------------------------------------------------------------
  ! Loop over p strings in expansion
    loop_pstrings: do i=1, pdetslen
      ! Generate orbital string
      call genorbstring( pdets(i), aelec, orbitals, adets, astring1 )
      ! 1-e contribution
      int1e1 = 0d0
      do j=1, aelec
        int1e1 = int1e1 + moints1( ind2val( astring1(j), astring1(j) ) )
      end do
      ! 2-e contribution
      int2e1 = 0d0
      do j=1, aelec
        do k=j, aelec
          int2e1 = int2e1 + moints2( index2e2( astring1(j), astring1(j),     &
                                     astring1(k), astring1(k) ) ) - moints2( &
                                     index2e2( astring1(j), astring1(k),     &
                                     astring1(j), astring1(k) ) )
        end do
      end do
      ! Loop over beta strings corresponding to p's in expansion
      do 
