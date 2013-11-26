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
    integer, dimension(orbitals-aelec) :: pexits1, qexits1
    
    !integer, dimension(:,:), allocatable :: srepinfo
    integer, dimension(:), allocatable :: srepinfo

    integer :: singexlen, eps1, eps2, eps3, xindx1, xindx2, vecindx1, vecindx2

    real*8 :: int1e1, int2e1, int2e2, int3e2, int3elweps
#ifdef HVACTION
    integer :: openstat
#endif
    !----------------------------------------------------------------
    ! Loop over alpha strings in expansion
    do i=1, pdetstrunc
    ! Generate orbital index string
      call genorbstring( pdets(i), aelec, orbitals, adets, astring )

    ! Generate list of excitations
      call possex1( astring, orbitals, aelec, (orbitals-aelec), pexits1 )

    ! Allocate info array
      singexlen = aelec*(orbitals-aelec)
      !if ( allocated(srepfinfo) ) deallocate(srepinfo)
      !allocate( srepinfo( singexlen, 4 ))
      !srepinfo = 0
      if ( allocated(srepinfo) ) deallocate(srepinfo)
      allocate( srepinfo(2))
    ! Loop over single excitations
      loopelec: do j=1, aelec
        loopexite: do k=1, orbitals-aelec
          call singrepinfo( astring, aelec, pexits1(k), j, orbitals, &
                            eps1, xindx1 )      
    ! Test if xindx1 is in pdets
          do l=1, pdetstrunc
            if ( xindx1 .eq. pdets(l) ) then
             ! srepinfo((j-1)*(orbitals-aelec)+k,1) = astring(j) ! Orbital exciting from
             ! srepinfo((j-1)*(orbitals-aelec)+k,2) = pexits1(k) ! Orbital exciting to
              srepinfo(1) = eps1       ! Parity
              srepinfo(2) = xindx1     ! Index of new string

              int1e1 = moints1(ind2val(astring(j),pexits1(k)))
              
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
#ifdef HVACTION
                    open( unit=50,file='alpha-cont.dets',status='unknown',position='append', iostat=openstat)
                    if ( openstat .ne. 0 ) stop "COULD NOT OPEN FILE"
                    write(unit=50,fmt=10) "From ", vecindx1," to ", vecindx2
                    close(unit=50)
10 format( 1x, A, I8, A, I8 )
#endif
                  ! Add the stored integrals and multiply by parity
                    vector2(vecindx2) = vector2(vecindx2) + eps1*( int1e1 + int2e1 + int2e2 )*&
                                         vector1(vecindx1)
                  end if
                end do
              end do
#ifdef HVACTION
              open( unit=51, file='alpha-single.cont',status='unknown', position='append', iostat=openstat )
              if ( openstat .ne. 0 ) stop "COULD NOT OPEN alpha-single.cont"
              write( unit=51 , fmt=10) "Looping over det ", pdets(i), " det num. ", i
              do m=1, cidim
                write( unit=51, fmt=11) m, vector2(m)
              end do
              close( unit=51 )
11 format( 1x, I9, F10.6 )
#endif
    ! Loop over additional replacements in alpha strings. ( m > j; n > k )
              do m=j+1, aelec      
                do n=k+1, orbitals-aelec
                  call doublerepinfo( astring, aelec, pexits1(k), j, pexits1(n), m, &
                                      orbitals, eps2, xindx2 )
                  int3elweps = eps2*(moints2(index2e2(astring(j),pexits1(k),astring(m),&
                                     pexits1(n) )) - moints2( index2e2(astring(j),pexits1(n),&
                                     astring(m),pexits1(k) )))

                  ! Test if xindx2 is in expansion
                  do o=1, pdetstrunc
                    if ( xindx2 .eq. pdets(o) ) then
                    ! Loop over q string without generating orbital index set
                    ! Again, q must correspond to both p and p**
                      do p=1, pstep(o)
                        do q=1, pstep(i)
                          if ( pstring(plocate(o)+p,2) .eq. pstring(plocate(i)+q,2) ) then
                          ! Do NOT generate orbital index strings. All that is needed is the index.
                            vecindx1 = indxk( xindx2, pstring(plocate(o)+p,2), belec, orbitals )
                            vecindx2 = indxk( pdets(i), pstring(plocate(o)+p,2), belec, orbitals )
                            vector2(vecindx2) = vector2(vecindx2) + int3elweps*vector1(vecindx1)
                          end if ! If q corresponds to both p and p**
                        end do ! Testing loop
                      end do ! Testing loop
                    end if ! p** is in expansion
                  end do ! Testing loop
                end do ! Loop over excitations
              end do ! Loop over electrons
      
      ! We now will loop over single excitations in both alpha and beta strings.
      ! Loop over beta strings in expansion associated with alpha string p.
      ! This information is in pstring.  We are looping over the second column.
              do m=1, pstep(i)
                ! Generate an orbital string set for q
                call genorbstring( pstring(plocate(i)+m,2), belec, orbitals, bdets, bstring )
                ! Generate a list of possible excitations
                call possex1( bstring, orbitals, belec, ( orbitals-belec ), qexits1 )

                ! Loop over single excitations in beta strings 
                do n=1, belec
                  do o=1, orbitals -belec
                    call singrepinfo( bstring, belec, qexits1(o), n, orbitals, eps3, xindx1 )

                    ! Loop over q indices that correspond to the single excitation r(pjv*) in p strings
                    do p=1, pstep(l)  ! l is the indice of p* in our p string list
                      if ( xindx1 .eq. pstring(plocate(l)+p,2) ) then
                        int3e2 = moints2( index2e2( astring(j), pexits1(k), bstring(n), qexits1(o) ))
                        vecindx1 = indxk(srepinfo(2), xindx1, belec, orbitals )
                        vecindx2 = indxk( pdets(i),pstring(plocate(i)+m,2), belec, orbitals )
                        vector2(vecindx2) = vector2(vecindx2) + srepinfo(1)*eps3*int3e2
                      end if
                    end do
                  end do ! Excitations
                end do ! Electrons
              end do ! q strings
            end if ! Test
          end do ! Test loop
        end do loopexite
      end do loopelec
    end do ! alpha strings in expansion
    return
  end subroutine
!====================================================================
!====================================================================
!>acthv_beta
!
! Subroutine to evaluate the contribution to Hv of single and double
!  excitations in beta strings.  This subroutine follows the above
!  subroutine, with beta and alpha strings swapped.
!--------------------------------------------------------------------
  subroutine acthv_beta( vector1, moints1, moints2, moints1len, moints2len,&
    pstring, pstep, plocate, pxreflist, pdets, qstring, qstep, qlocate,    &
    qxreflist, qdets, cidim, pdetstrunc, qdetstrunc, adets, bdets, aelec,  &
    belec, orbitals, vector2 )
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
   
    integer, dimension(aelec) :: astring
    integer, dimension(belec) :: bstring, bstring2
    integer, dimension(orbitals-aelec) :: pexits1, qexits1
    
    !integer, dimension(:,:), allocatable :: srepinfo
    integer, dimension(:), allocatable :: srepinfo

    integer :: singexlen, eps1, eps2, eps3, xindx1, xindx2, vecindx1, vecindx2

    real*8 :: int1e1, int2e1, int2e2, int3elweps 
    !----------------------------------------------------------------
    do i=1, qdetstrunc
      ! Generate orbital index string
      call genorbstring ( qdets(i), belec, orbitals, bdets, bstring )

      ! Genererate list of excitations
      call possex1( bstring, orbitals, belec, ( orbitals-belec ), pexits1 )

      ! Allocate info array
      if ( allocated(srepinfo) ) deallocate(srepinfo)
      allocate( srepinfo(2) )
      ! Loop over single excitations
      loopelec: do j=1, belec
        loopexite: do k=1, orbitals - belec
          call singrepinfo( bstring, belec, pexits1(k), j, orbitals, &
                            eps1, xindx1 )
          ! Test if xindx1 is in qdets
          do l=1, qdetstrunc
            if ( xindx1 .eq. qdets(l) ) then
              srepinfo(1)=eps1
              srepinfo(2)=xindx1

              int1e1 = moints1(ind2val(bstring(j),pexits1(k)))

              int2e1 = 0d0
              do m=1, belec
                if ( m .ne. j ) then
                  int2e1 = int2e1 + moints2( index2e2( bstring(m), bstring(m), bstring(j),   &                                     pexits1(k))) - moints2( index2e2( bstring(m), bstring(j),&
                                      bstring(m), pexits1(k)))
                end if
              end do

                ! Loop over corresponding p strings
              do m=1, qstep(l)
                do n=1, qstep(i)
                ! Test if corresponding p of qdets(i) corresponds to qdets(l)
                  if ( qstring(qlocate(l)+m,2) .eq. qstring(qlocate(i)+n,2)) then
                    !Generate orbital index string
                    call genorbstring(qstring(qlocate(i)+n,2),aelec, orbitals, adets, astring )
                    int2e2 = 0d0
                    do o=1, aelec
                      int2e2 = int2e2 + moints2( index2e2( astring(o), astring(o), bstring(j), &
                                          pexits1(k) ) )
                    end do
                    ! Find index of determinant C[r(p,qjv*)]
                    vecindx1 = indxk( qstring(qlocate(l)+m,2),xindx1, belec, orbitals )
                    ! Find index of determinant V[(p,q)]
                    vecindx2 = indxk( qstring(qlocate(l)+m,2),qdets(i), belec, orbitals )
                    ! Add the stored integerals and multiply by parity
                    vector2(vecindx2) = vector2(vecindx2) + eps1*(int1e1+int2e1+int2e2)*&
                                        vector1(vecindx1)
                  end if
                end do
              end do

              ! Loop over additional replacements in beta strings ( m > j; n> k )
              do m=j+1, belec
                do n=k+1, orbitals-belec
                  call doublerepinfo( bstring, belec, pexits1(k), j, pexits1(n), m, &
                                      orbitals, eps2, xindx2 )
                  int3elweps = eps2*( moints2(index2e2(bstring(j),pexits1(k),bstring(m),                                        pexits1(n) ) ) - moints2( index2e2( bstring(j),   &
                                        pexits1(n), bstring(m), pexits1(k) )))
                  ! Test if xindx2 is in expansion
                  do o=1, qdetstrunc
                    if ( xindx2 .eq. qdets(o) ) then
                    ! Loop over p string w/o generating orbital index set
                    ! p must correspond to both q and q**
                      do p=1, qstep(o)
                        do q=1, qstep(i)
                          if ( qstring(qlocate(o)+p,2) .eq. qstring(qlocate(i)+q,2) ) then
                            vecindx1 = indxk( qstring(qlocate(o)+p,2), xindx2, belec, orbitals )
                            vecindx2 = indxk( qstring(qlocate(o)+p,2), qdets(i),belec, orbitals)
                            vector2(vecindx2) = vector2(vecindx2) + int3elweps*vector1(vecindx1)
                          end if
                        end do
                      end do
                    end if
                  end do
                end do
              end do
            end if
          end do
        end do loopexite
      end do loopelec
    end do
    return
  end subroutine
!====================================================================
!====================================================================
end module
