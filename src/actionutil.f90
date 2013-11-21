! Module containing subroutines for computing the action of H on 
!  some vector
!====================================================================
module actionutil
contains
!====================================================================
!====================================================================
!> acthv_diag
!
! Diagonal contribution
!--------------------------------------------------------------------
  subroutine acthv_diag( vector1, dgls, length, vector2 )
    implicit none
    integer, intent(in) :: length
    real*8, dimension(length), intent(in)  :: vector1, dgls
    real*8, dimension(length), intent(inout) :: vector2
    integer :: i
    !----------------------------------------------------------------
    do i=1, length
      vector2(i) = vector2(i) + dgls(i)*vector1(i)
    end do
    return
  end subroutine acthv_diag
!====================================================================
!====================================================================
!> acthv_alpha
!
! Single & double excitations in alpha strings and single excitations
!  in alpha and beta strings
!--------------------------------------------------------------------
  subroutine acthv_alpha( vector1, moints1, moints2, moints1len,    &
    moints2len, pstring, pstep, plocate, pxreflist, pdets, qstring, &
    qstep, qlocate, qxreflist, qdets, cidim, pdetstrunc, qdetstrunc,&
    adets, bdets, aelec, belec, orbitals, vector2 )
    use detci2
    use detci5
    implicit none
    integer, intent(in) :: moints1len, moints2len, cidim, pdetstrunc,&
                           qdetstrunc, adets, bdets, aelec, belec,   &
                           orbitals
    
    integer, dimension(cidim, 2), intent(in)  :: pstring, qstring
    integer, dimension(pdetstrunc), intent(in):: pstep, plocate, pdets
    integer, dimension(qdetstrunc), intent(in):: qstep, qlocate, qdets
    integer, dimension(cidim), intent(in)     :: pxreflist, qxreflist
    
    real*8, dimension(moints1len), intent(in) :: moints1
    real*8, dimension(moints2len), intent(in) :: moints2
    real*8, dimension(cidim), intent(in)      :: vector1
  
    real*8, dimension(cidim), intent(inout)   :: vector2

    integer :: i, j, k, l, m, n
   
    integer, dimension(aelec) :: astring
    integer, dimension(belec) :: bstring
    integer, dimension(orbitals-aelec) :: pexits1
    
    integer, dimension(:,:), allocatable :: srepinfo

    integer :: singexlen, eps1, xindx1, vecindx1, vecinx2

    real*8 :: int1e1, int2e1
    !----------------------------------------------------------------
    ! Loop over alpha strings in expansion
    do i=1, pdetstrunc
    ! Generate orbital index string
      call genorbstring( pdets(i), aelec, orbitals, adets, astring )

    ! Generate list of excitations
      call possex1( astring, orbitals, aelec, (orbitals-aelec), pexits1 )

    ! Allocate info array
      singexlen = aelec*(orbitals-aelec)
      if ( allocated(srepfinfo) ) deallocate(srepinfo)
      allocate( srepinfo( singexlen, 4 ))

    ! Loop over single excitations
      loopelec: do j=1, aelec
        loopexite: do k=1, orbitals-aelec
          call singrepinfo( astring1, aelec, pexits1(k), j, orbitals, &
                            eps1, xindx1 )      
    ! Test if xindx1 is in pdets
          do l=1, pdetstrunc
            if ( xindx1 .eq. pdets(l) ) then
              srepinfo((j-1)*(orbitals-aelec)+k,1) = astring(j) ! Orbital exciting from
              srepinfo((j-1)*(orbitals-aelec)+k,2) = pexits1(k) ! Orbital exciting to
              srepinfo((j-1)*(orbitals-aelec)+k,3) = eps1       ! Parity
              srepinfo((j-1)*(orbitals-aelec)+k,4) = xindx1     ! Index of new string

              int1e1 = moints1(ind2val(astring(j),pexits(k)))
              
              int2e1 = 0d0
              do m=1, aelec
                if ( m .ne. j ) then
                  int2e1 = int2e1 + moints2( index2e2( astring(m), astring(m), astring(j),   &
                                    pexits1(k))) - moints2( index2e2( astring(m), astring(j),&
                                    astring(m), pexits1(k)))
                end if
              end do
            
    ! Loop over corresponding q strings. This info is in pstrings()
    ! These q strings correspond to both p and r(pjv*)
              
    

end module
