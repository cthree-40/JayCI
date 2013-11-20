! Subroutine to compute diagonal contribution <K|H|K>
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. Chem., The Johns Hopkins University
!
! Last edit: 11-13-13
!====================================================================
!Input:
! actdets   = active determinants                 integer array   1-d
! actdetslen= length of actdets()                 integer scalar
! totels    = total number of electrons           integer scalar
! aelec     = alpha electrons                     integer scalar
! belec     = beta electrons                      integer scalar
! orbitals  = number of MO's                      integer scalar
! adets     = alpha determinants                  integer scalar
! bdets     = beta determinants                   integer scalar
! max1e     = number of 1-e integrals             integer scalar
! max2e     = number of 2-e integerals            integer scalar
! moints1   = 1-e integrals                       real*8  array   1-d
! moints2   = 2-e integrals                       real*8  array   1-d
!Output:
! dglsout   = diagonals
!
!====================================================================
subroutine diagonal( actdets, actdetslen, totels, aelec, belec, orbitals, &
  adets, bdets, max1e, max2e, moints1, moints2, dglsout )

  use detci5
  use detci2

  implicit none


! ...input integer scalars...
  integer, intent(in) :: actdetslen, totels, aelec, belec, orbitals, adets,&
                         bdets, max1e, max2e

! ...input real*8 arrays...
  real*8, dimension(max1e) :: moints1
  real*8, dimension(max2e) :: moints2

! ...input integer arrays...
  integer, dimension(actdetslen) :: actdets

! ...OUTPUT real*8 arrays...
  real*8, dimension(actdetslen) :: dglsout

! ...loop integer scalars...
  integer :: i, j, k, l
  integer :: vecindx1

! ...loop integer arrays...
  integer, dimension(aelec) :: astring1, astring2
  integer, dimension(belec) :: bstring1, bstring2

! ...loop real*8 scalars...
  real*8 :: int1e1, int2e1, int2e2

!--------------------------------------------------------------------

  dglsout = 0d0

! Loop over the alpha strings
  do i=1, adets

! Generate orbital string
    if ( i .eq. 1 ) then
      do j=1, aelec
        astring1(j) = j
      end do
    else
      call nstrfnd( astring2, aelec, orbitals, adets, astring1 )
    end if

! 1-e alpha contribution
    int1e1=0d0
    do j=1, aelec
      int1e1 = int1e1 + moints1( ind2val( astring1(j), astring1(j) ))
    end do

! 2-e alpha contribution
    int2e1=0d0
    do j=1, aelec
      do k=1, j
        int2e1 = int2e1 + moints2( index2e2( astring1(j), astring1(j), &
                                   astring1(k), astring1(k) ) )  -     &
                          moints2( index2e2( astring1(j), astring1(k), &
                                   astring1(j), astring1(k) ) )
      end do
    end do

! Loop through beta strings
    do j=1, bdets
      if ( j .eq. 1 ) then
        do k=1, belec
          bstring1(k) = k
        end do
      else
        call nstrfnd( bstring2, belec, orbitals, bdets, bstring1 )
      end if

! 2-e contribution from alpha and beta strings
      int2e2 = 0d0
      do k=1, aelec
        do l=1, belec
          int2e2 = int2e2 + moints2( index2e2( astring1(k), astring1(k), &
                                     bstring1(l), bstring1(l) ) )
        end do
      end do

! Add sums to diagonal elment <p,q| H |p,q>
      vecindx1 = indxk( i, j, belec, orbitals )

      dglsout(vecindx1) = dglsout(vecindx1) + int1e1 + int2e1 + int2e2

! Close loop over beta strings
      bstring2 = bstring1
      bstring1 = 0
    end do

! Close loop over alpha strings
    astring2 = astring1
    astring1 = 0
  end do


! Loop over beta determinants
  do i=1, bdets
    if ( i .eq. 1 ) then
      do j=1, belec
        bstring1(j) = j
      end do
    else
      call nstrfnd( bstring2, belec, orbitals, bdets, bstring1 )
    end if

! 1-e beta contribution
    int1e1 = 0d0
    do j=1, belec
      int1e1 = int1e1 + moints1( ind2val( bstring1(j), bstring1(j) ))
    end do

! 2-e beta contribution
    int2e1 = 0d0
    do j=1, belec
      do k=1, j
        int2e1 = int2e1 + moints2( index2e2( bstring1(j), bstring1(j), &
                                   bstring1(k), bstring1(k) ) ) -      &
                          moints2( index2e2( bstring1(j), bstring1(k), &
                                   bstring1(j), bstring1(k) ) )
      end do
    end do 

! Loop through alpha strings w/o generating new orbital index sets
    do j=1, adets
      vecindx1 = indxk( j, i, belec, orbitals )
      dglsout(vecindx1) = dglsout(vecindx1) + int1e1 + int2e1
    end do

! Close loop over beta strings
    bstring2 = bstring1
    bstring1 = 0
  end do

  return

end subroutine









