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
  adets, bdets, moints1, moints1len, moints2, moints2len, &
  plocate, qlocate, diagonals )

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
  integer, dimension(moints1len), intent(in) :: moints1
  integer, dimension(moints2len), intent(in) :: moints2
! ...output real*8 arrays...
  real*8, dimension(cidim), intent(inout) :: diagonals
! ...integer arrays...
!  integer, dimension(cidim ) :: alphadiagindex, betadiagindex
!--------------------------------------------------------------------
  diagonals = 0d0
! Alpha string contribution
!  call formdiagindexes( pstring, cidim, belec, orbitals, alphadiagindex )
!  call contribdiag( 'a', detlist, cidim, pstring, cidim, qstring, cidim, &
!                     pstep, psteplen, qstep, qsteplen, & 
!                     aelec, belec, orbitals, adets, bdets, pdets, &
!                     pdetslen, qdets, qdetslen, moints1, moints1len,     &
!                     moints2, moints2len, pxreflist,qxreflist, xreflistlen, &
!                     plocate, qlocate, diagonals )
! Beta  string contribution
!  call contribdiag( 'b', detlist, cidim, pstring, cidim, qstring, cidim, &
!                     pstep, psteplen, qstep, qsteplen, & 
!                     aelec, belec, orbitals, adets, bdets, pdets, &
!                     pdetslen, qdets, qdetslen, moints1, moints1len,     &
!                     moints2, moints2len, pxreflist,qxreflist, xreflistlen, &
!                     plocate, qlocate, diagonals )
  do i=1, cidim
    diagonals(i) = ham_element_diag( i, moints1, moints1len, moints2, moints2len, &
                   aelec, belec, orbitals )
  end do
  return
end subroutine
!====================================================================
!====================================================================
! form index of determinants
!--------------------------------------------------------------------
subroutine formdiagindexes( string, length, belec, orbitals, indexstring )
!Input:
! string     = p/q string pairs                           integer array  2-d
! length     = dimension of CI                            integer scalar
!Output:
! indextring = determinant of corresponding string pairs  integer array  1-d
!--------------------------------------------------------------------
  use detci2
  implicit none
! ...input integer scalars...
  integer, intent(in) :: length, belec, orbitals
! ...input integer arrays...
  integer, dimension(length,2),intent(in) :: string
! ...output integet arrays...
  integer, dimension(length),  intent(out) :: indexstring
! ...loop integer scalars...
  integer :: i
!--------------------------------------------------------------------
! Loop over string pairs
  do i=1, length
    indexstring(i) = indxk( string(i,1), string(i,2), belec, orbitals )
  end do
  return
end subroutine
!====================================================================
!====================================================================
! Subroutine to compute diagonal contribution of alpha/beta strings
!--------------------------------------------------------------------
subroutine contribdiag( spin, detlist, detlistlen, pstring, pstringlen,  &
  qstring, qstringlen, pstep, psteplen, qstep, qsteplen,  aelec, belec,  &
  orbitals, adets, bdets, pdets, pdetslen,    &
  qdets, qdetslen, moints1, moints1len, moints2, moints2len, pxreflist,   &
  qxreflist, xreflistlen, plocate, qlocate, diagonals )
!Input:
! See above.
!--------------------------------------------------------------------
  use detci2
  use detci5
  implicit none
! ...input integer scalars...
  integer, intent(in) :: detlistlen, pstringlen, qstringlen, psteplen, &
                         qsteplen, moints1len, moints2len,pdetslen, qdetslen,&
                         adets, bdets, aelec, belec, orbitals, xreflistlen
! ...input integer arrays...
  integer, dimension(detlistlen),  intent(in) :: detlist
  integer, dimension(pstringlen,2),intent(in) :: pstring
  integer, dimension(qstringlen,2),intent(in) :: qstring
  integer, dimension(psteplen), intent(in)    :: pstep, plocate
  integer, dimension(qsteplen), intent(in)    :: qstep, qlocate
  integer, dimension(pdetslen), intent(in)    :: pdets
  integer, dimension(qdetslen), intent(in)    :: qdets
  integer, dimension(xreflistlen),intent(in)  :: pxreflist, qxreflist
! ...input real*8 arrays...
  real*8, dimension(moints1len), intent(in)   :: moints1
  real*8, dimension(moints2len), intent(in)   :: moints2
! ...input character strings...
  character(len=1), intent(in) :: spin
! ...output real*8 arrays...
  real*8, dimension(detlistlen), intent(inout) :: diagonals
! ...loop integer scalars...
  integer :: i, j, k, l
! ...integer arrays...
  integer, dimension(aelec) :: astring1
  integer, dimension(belec) :: bstring1
! ...real*8 scalars...
  real*8 :: int1e1, int2e1, int2e2
! ...integer scalars...
  integer :: vecindx1
!--------------------------------------------------------------------

  if ( spin .eq. 'a' ) then
! Loop over p strings in expansion
    do i=1, pdetslen
! Generate orbital string
      call genorbstring( pdets(i), aelec, orbitals, adets, astring1 )
! 1-e contribution
      int1e1=0d0
      do j=1, aelec
        int1e1 = int1e1 + moints1( ind2val( astring1(j), astring1(j)) )
      end do
! 2-e contribution
      int2e1=0d0
      do j=1, aelec
        do k=1, j
          int2e1 = int2e1 + moints2(index2e2(astring1(j),astring1(j),    &
                            astring1(k),astring1(k))) - moints2(index2e2(&
                            astring1(j),astring1(k),astring1(j),astring1(k)))
        end do
      end do
! Loop through beta strings corresponding to p in expasion
      do j=1, pstep(i)
! Generate orbital string
        call genorbstring( pstring(plocate(i)+j, 2), belec, orbitals, bdets, bstring1 )
! 2-e contribution from alpha and beta strings
        int2e2 = 0d0
        do k=1, aelec
          do l=1, belec
            int2e2 = int2e2 + moints2( index2e2( astring1(k), astring1(k), &
                              bstring1(l), bstring1(l)))
          end do
        end do

! Add sums to diagonal element <p,q| H |p,q>
        diagonals(pxreflist(plocate(i)+j)) = diagonals(pxreflist(plocate(i)+j)) + int1e1 + int2e1 + int2e2
      end do
    end do
  else
! Loop over q strings in expansion
    do i=1, qdetslen
! Generate orbital string
      call genorbstring( qdets(i), belec, orbitals, bdets, bstring1 )
! 1-e contribution
      int1e1=0d0
      do j=1, belec
        int1e1 = int1e1 + moints1( ind2val( bstring1(j), bstring1(j)) )
      end do
! 2-e contribution
      int2e1=0d0
      do j=1, belec
        do k=1, j
          int2e1 = int2e1 + moints2(index2e2(bstring1(j),bstring1(j),    &
                            bstring1(k),bstring1(k))) - moints2(index2e2(&
                            bstring1(j),bstring1(k),bstring1(j),bstring1(j)))
        end do
      end do
! Loop through alpha strings corresponding to q in expansion
      do j=1, qstep(i)
! Generate orbital string
        call genorbstring( qstring(qlocate(i)+j,2), belec, orbitals, bdets, bstring1 )
! 2-e contribution from alpha and beta strings
        int2e2 = 0d0
!        do k=1, belec
!          do l=1, aelec
!            int2e2 = int2e2 + moints2( index2e2( bstring1(k), bstring1(k), &
!                              astring1(l), astring1(l)) )
!          end do
!        end do
! Add sums to diagonal element <p,q|H|p,q>
        diagonals(qxreflist(qlocate(i)+j)) = diagonals(qxreflist(qlocate(i)+j)) + int1e1 + int2e1 + int2e2
      end do
    end do
  end if
  return
end subroutine
!====================================================================
!====================================================================
          
    
