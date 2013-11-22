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

    integer :: i, j, k, l, m, n, o, p, q, r, s
   
    integer, dimension(aelec) :: astring, astring2
    integer, dimension(belec) :: bstring
    integer, dimension(orbitals-aelec) :: pexits1
    
    integer, dimension(:,:), allocatable :: srepinfo

    integer :: singexlen, eps1, eps2, xindx1, xindx2, vecindx1, vecindx2

    real*8 :: int1e1, int2e1, int2e2, int3elweps
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
    ! We are evaluating <p,q|H|p*,q>; therefore, |p,q> and |p*,q> must be in
    ! expansion, so q must correspond to both p and p*
    ! We use the arrays plocate and pstep to find our p's and search for our
    ! corresponding q's.
              
              do m=1, pstep(l)
                do n=1, pstep(i)
                ! Test if the corresponding q of pdets(i) corresponds to
                ! pdets(l)
                  if ( pstring(plocate(l)+m,2) .eq. pstring(plocate(i)+n,2) ) then
                  ! Generate orbital index string
                    call genorbstring(pstring(plocate(i)+n,2), belec, orbitals, bdets, bstring )
                    int2e2=0d0
                    do o=1, belec
                      int2e2 = int2e2 + moints2( index2e2( bstring(o), bstring(o), astring(j), &
                                        pexits1(k) ) )
                    end do
                  ! Find index of determinant C[r(pjv*,q)]
                    vecindx1 = indxk( xindx1, pstring(plocate(l)+m,2), belec, orbitals )
                  ! Find index of determinant V[(p,q)]
                    vecindx2 = indxk( pdets(i), pstring(plocate(l)+m,2), belec, orbitals )
                  ! Add the stored integrals and multiply by parity
                    vector2(vecindx2) = vector2(vecindx2) + eps1*( int1e1 + int2e1 + int2e2 )*&
                                         vector1(vecindx1)
                  end if
                end do
              end do

    ! Loop over additional replacements in alpha strings. ( m > j; n > k )
              do m=j+1, aelec      
                do n=k+1, orbitals-aelec
                  call doublerepinfo( astring, aelec, pexits1(k), j, pexits1(n), m, &
                                      orbitals, eps2, xindx2 )
                  int3elweps = eps2*(moints2(index2e2(astring(j),pexits1(k),astring(m),&
                                     pexits1(n) )) - moints2( index2e2(astring(j),pexits1(n),&
                                     astring(m),pexits1(k) ))

                  ! Test if xindx2 is in expansion
                  do o=1, pdetstrunc
                    if ( xindx2 .eq. pdets(o) ) then
                    ! Loop over q string without generating orbital index set
                    ! Again, q must correspond to both p and p**
                      do p=1, pstep(o)
                        do q=1, pstep(i)
                          if ( pstring(plocate(o)+p,2) .eq. pstring(plocate(i)+q,2) ) then
                          ! Do NOT generate orbital index strings. All that is needed is the index.
                            vecindx1 = indxk( xindx2, pstring(plocate(o)+p), belec, orbitals )
                            vecindx2 = indxk( pdets(i), pstring(plocate(o)+p), belec, orbitals )
                            vector2(vecindx2) = vector2(vecindx2) + int3elweps*vector1(vecindx1)
                          end if ! If q corresponds to both p and p**
                        end do ! Testing loop
                      end do ! Testing loop
                    end if ! p** is in expansion
                  end do ! Testing loop
                end do ! Loop over excitations
              end do ! Loop over electrons
            end if ! If p* is in expansion
          end do ! Testing loop
        end do loopexite
      end do loopelec
      
      ! We now will loop over single excitations in both alpha and beta strings.
      ! Loop over beta strings in expansion associated with alpha string p.
      ! This information is in pstring.  We are looping over the second column.
    






end module
