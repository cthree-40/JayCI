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
    moints2len, pstring, pstep, plocate, pdets, qstring, &
    qstep, qlocate, qdets, cidim, pdetstrunc, qdetstrunc,&
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
    character(len=11) :: debugfl
    integer :: debugunit
    integer :: openstat
#endif
#ifdef HVACTION
9  format(1x,A)
10 format(1x,A,I8)
11 format(1x,I3,I3)
12 format(1x,I2)
13 format(1x, A, I8, A, I8 )
14 format(1x, A, F10.6)
#endif
    !----------------------------------------------------------------
    ! Loop over alpha strings in expansion
    do i=1, pdetstrunc
!___________________________________________________
#ifdef HVACTION
    if ( pdets(i) .eq. 27 ) then
      debugunit=5
      debugfl='debug.acthv'
      open(unit=debugunit,file=debugfl, status='unknown', iostat=openstat )
      if ( openstat .ne. 0 ) stop "**** COULD NOT OPEN DEBUGGING FILE *****"
      write( unit=debugunit, fmt=10 ) "Looping over alpha string ", pdets(i)
      write( unit=debugunit, fmt=14 ) "VECTOR2(1) : " , vector2(1)
    end if
#endif
!___________________________________________________
    ! Generate orbital index string
      call genorbstring( pdets(i), aelec, orbitals, adets, astring )
!___________________________________________________
#ifdef HVACTION
      if ( pdets(i) .eq. 27 ) then
        write( unit=debugunit, fmt=9) "This string looks like : "
        write( unit=debugunit, fmt=11)  astring(1)
        write( unit=debugunit, fmt=9) " "
      end if
#endif 
    ! Generate list of excitations
      call possex1( astring, orbitals, aelec, (orbitals-aelec), pexits1 )
!___________________________________________________
#ifdef HVACTION
      if ( pdets(i) .eq. 27 ) then
        write( unit=debugunit, fmt=9) "Possible excitations for this string are: "
        do j=1, orbitals-aelec
          write( unit=debugunit, fmt=10) "    ", pexits1(j)
        end do
        write( unit=debugunit, fmt=9) " "
      end if
#endif
!____________________________________________________      
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
!__________________________________________________
#ifdef HVACTION
          if ( pdets(i) .eq. 27 ) then
            write( unit=debugunit,fmt=10 ) " Exciting electron: ", j
            write( unit=debugunit,fmt=10 ) " Into orbital:      ", pexits1(k)
          end if
#endif
!__________________________________________________
    ! Test if xindx1 is in pdets
          test_p: do l=1, pdetstrunc
            if ( xindx1 .eq. pdets(l) ) then
!__________________________________________________
#ifdef HVACTION
              if ( pdets(i) .eq. 27 ) then
                write( unit=debugunit, fmt=9 ) "Entering contribution Loop..."
                write( unit=debugunit, fmt=10) "Index of new string is ( found from singrepinfo ) : ", pdets(l)
                write( unit=debugunit, fmt=10) " xindx1 : ", xindx1
              end if
#endif
!___________________________________________________
             ! srepinfo((j-1)*(orbitals-aelec)+k,1) = astring(j) ! Orbital exciting from
             ! srepinfo((j-1)*(orbitals-aelec)+k,2) = pexits1(k) ! Orbital exciting to
              srepinfo(1) = eps1       ! Parity
              srepinfo(2) = xindx1     ! Index of new string

             call eval_singlex1( astring, pexits1, aelec, moints1, moints2, moints1len, &
                                 moints2len, orbitals, j, k, int1e1, int2e1 )             
    ! Loop over corresponding q strings. This info is in pstrings()
    ! These q strings correspond to both p and r(pjv*)
    ! We are evaluating <p,q|H|p*,q>; therefore, |p,q> and |p*,q> must be in
    ! expansion, so q must correspond to both p and p*
    ! We use the arrays plocate and pstep to find our p's and search for our
    ! corresponding q's.
              
              do m=1, pstep(l)
                test_p2: do n=1, pstep(i)
                ! Test if the corresponding q of pdets(i) corresponds to
                ! pdets(l)
                  if ( pstring(plocate(l)+m,2) .eq. pstring(plocate(i)+n,2) ) then
!_______________________________________________
#ifdef HVACTION
                    if ( pdets(i) .eq. 27 ) then
                      write( unit=debugunit, fmt=10 ) " Q string is ", pstring(plocate(l)+m,2)
                      write( unit=debugunit, fmt=9) " "
                    end if
#endif
!_______________________________________________

                  ! Compute contribution
                    call eval_singlex2( pstring(plocate(l)+m,2), belec, orbitals, bdets, &
                                   astring(j), pexits1(k), moints1, moints2, moints1len, &
                                   moints2len, int2e2 )
                  ! Find index of determinant C[r(pjv*,q)]
                    vecindx1 = indxk( xindx1, pstring(plocate(l)+m,2), belec, orbitals )
                  ! Find index of determinant V[(p,q)]
                    vecindx2 = indxk( pdets(i), pstring(plocate(l)+m,2), belec, orbitals )
!________________________________________________
#ifdef HVACTION
                    if ( pdets(i) .eq. 27 ) then
                      write( unit=debugunit, fmt=13 ) "        Adding contribution from ", vecindx1, "to", vecindx2
                      write( unit=debugunit, fmt=14 ) "        vector(vecindx1) is ", vector1(vecindx1)
                      write( unit=debugunit, fmt=14 ) "        vector(1) is ", vector2(1)
                    end if
#endif
!________________________________________________
 
                  ! Add the stored integrals and multiply by parity
                    vector2(vecindx2) = vector2(vecindx2) + eps1*( int1e1 + int2e1 + int2e2 )*&
                                         vector1(vecindx1)
!_______________________________________________
#ifdef HVACTION
                    if ( pdets(i) .eq. 27 ) then
                      write( unit=debugunit,fmt=14) "         vector(1) is ", vector2(1)
                    end if
#endif
!_______________________________________________
                    exit test_p2
                  end if
                end do test_p2
              end do

    ! Loop over additional replacements in alpha strings. ( m > j; n > k )
!________________________________________________
#ifdef HVACTION
              if ( pdets(i) .eq. 27 ) then
                write( unit=debugunit, fmt=9 ) " Now for addition replacements..."
              end if
#endif
!_________________________________________________
              do m=j+1, aelec      
                do n=k+1, orbitals-aelec
                  call doublerepinfo( astring, aelec, pexits1(k), j, pexits1(n), m, &
                                      orbitals, eps2, xindx2 )
                  int3elweps = eps2*(moints2(index2e2(astring(j),pexits1(k),astring(m),&
                                     pexits1(n) )) - moints2( index2e2(astring(j),pexits1(n),&
                                     astring(m),pexits1(k) )))

                  ! Test if xindx2 is in expansion
                  test_x2: do o=1, pdetstrunc
                    if ( xindx2 .eq. pdets(o) ) then
                    ! Loop over q string without generating orbital index set
                    ! Again, q must correspond to both p and p**
                      do p=1, pstep(o)
                        do q=1, pstep(i)
                          if ( pstring(plocate(o)+p,2) .eq. pstring(plocate(i)+q,2) ) then
                          ! Do NOT generate orbital index strings. All that is needed is the indoex.
                            vecindx1 = indxk( xindx2, pstring(plocate(o)+p,2), belec, orbitals )
                            vecindx2 = indxk( pdets(i), pstring(plocate(o)+p,2), belec, orbitals )
                            vector2(vecindx2) = vector2(vecindx2) + int3elweps*vector1(vecindx1)
!_______________________________________________
#ifdef HVACTION            
                            if ( pdets(i) .eq. 27 ) then
                              write( unit=debugunit,fmt=13) " Adding the contribution from ", vecindx1, "to", vecindx2
                            end if
#endif
!_______________________________________________
                          end if ! If q corresponds to both p and p**
                        end do ! Testing loop
                      end do ! Testing loop
                      exit test_x2
                    end if ! p** is in expansion
                  end do test_x2! Testing loop
                end do ! Loop over excitations
              end do ! Loop over electrons
      
      ! We now will loop over single excitations in both alpha and beta strings.
      ! Loop over beta strings in expansion associated with alpha string p.
      ! This information is in pstring.  We are looping over the second column.
!__________________________________________________
#ifdef HVACTION
              if ( pdets(i) .eq. 27 ) then
                write( unit=debugunit,fmt=9) " Now entering single excitations in both alpha and beta strings..."
                write( unit=debugunit,fmt=14) "   vector2(1): ", vector2(1)
              end if
#endif
!___________________________________________________
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
                    test_p3: do p=1, pstep(l)  ! l is the indice of p* in our p string list
                      if ( xindx1 .eq. pstring(plocate(l)+p,2) ) then
                        int3e2 = moints2( index2e2( astring(j), pexits1(k), bstring(n), qexits1(o) ))
                        vecindx1 = indxk(srepinfo(2), xindx1, belec, orbitals )
                        vecindx2 = indxk( pdets(i),pstring(plocate(i)+m,2), belec, orbitals )
                        vector2(vecindx2) = vector2(vecindx2) + srepinfo(1)*eps3*int3e2*vector1(vecindx1)
                        exit test_p3
                      end if
                    end do test_p3
                  end do ! Excitations
                end do ! Electrons
              end do ! q strings
              exit test_p
            end if ! Test
          end do test_p! Test loop
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
    pstring, pstep, plocate, pdets, qstring, qstep, qlocate,    &
    qdets, cidim, pdetstrunc, qdetstrunc, adets, bdets, aelec,  &
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
          test_q: do l=1, qdetstrunc
            if ( xindx1 .eq. qdets(l) ) then
              srepinfo(1)=eps1
              srepinfo(2)=xindx1

              call eval_singlex1( bstring, pexits1, belec, moints1, moints2, moints1len, &
                                  moints2len, orbitals, j, k, int1e1, int2e1 )
                ! Loop over corresponding p strings
              do m=1, qstep(l)
                test_q2: do n=1, qstep(i)
                ! Test if corresponding p of qdets(i) corresponds to qdets(l)
                  if ( qstring(qlocate(l)+m,2) .eq. qstring(qlocate(i)+n,2)) then
                    ! Compute contribution
                    call eval_singlex2( qstring(qlocate(l)+m,2), aelec, orbitals, adets, &
                                        bstring(j), pexits1(k), moints1, moints2, moints1len, &
                                        moints2len, int2e2 )
                    ! Find index of determinant C[r(p,qjv*)]
                    vecindx1 = indxk( qstring(qlocate(l)+m,2),xindx1, belec, orbitals )
                    ! Find index of determinant V[(p,q)]
                    vecindx2 = indxk( qstring(qlocate(l)+m,2),qdets(i), belec, orbitals )
                    ! Add the stored integerals and multiply by parity
                    vector2(vecindx2) = vector2(vecindx2) + eps1*(int1e1+int2e1+int2e2)*&
                                        vector1(vecindx1)
                    exit test_q2
                  end if
                end do test_q2
              end do

              ! Loop over additional replacements in beta strings ( m > j; n> k )
              do m=j+1, belec
                do n=k+1, orbitals-belec
                  call doublerepinfo( bstring, belec, pexits1(k), j, pexits1(n), m, &
                                      orbitals, eps2, xindx2 )
                  int3elweps = eps2*( moints2(index2e2(bstring(j),pexits1(k),bstring(m), &
                                        pexits1(n) ) ) - moints2( index2e2( bstring(j),   &
                                        pexits1(n), bstring(m), pexits1(k) )))
                  ! Test if xindx2 is in expansion
                  test_x2: do o=1, qdetstrunc
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
                      exit test_x2
                    end if
                  end do test_x2
                end do
              end do
              exit test_q
            end if
          end do test_q
        end do loopexite
      end do loopelec
    end do
    return
  end subroutine
!====================================================================
!====================================================================
!>eval_singlex1
!
!Subroutine to compute single exitation contributions
!--------------------------------------------------------------------
  subroutine eval_singlex1( spstring, spexits, spelec, moints1,  moints2, &
    moints1len, moints2len, orbitals, bound_elec, excited_elec, cont_1elec,&
    cont_2elec)
  !Input:
  ! spstring    = alpha/beta string
  ! spexits     = alpha/beta exited orbital string
  ! spelec      = # of alpha/beta electrons
  ! moints1     = 1-e integrals
  ! moints2     = 2-e integrals
  ! moints1len  = length of moints1
  ! moints2len  = length of moints2
  ! orbitals    = # of MO's
  ! bound_elec  = index of electron being excited in spstring
  ! exited_elec = index of new orbital in spexits
  ! cont_1elec  = one electron contribution
  ! cont_2elec  = two electron contribution
    use detci5
    implicit none
  ! ...input integer scalars...
    integer, intent(in) :: spelec, moints1len, moints2len, orbitals, bound_elec,&
                         excited_elec
  ! ...input integer arrays...
    integer, dimension(spelec), intent(in)          :: spstring
    integer, dimension(orbitals-spelec), intent(in) :: spexits
  ! ...input real*8 arrays...
    real*8, dimension(moints1len),intent(in)        :: moints1
    real*8, dimension(moints2len),intent(in)        :: moints2
  ! ...output real*8 scalars...
    real*8, intent(out) :: cont_1elec, cont_2elec
  ! ...loop integer scalars...
    integer :: i, j, k
  !--------------------------------------------------------------------
  ! One electron contribution
    cont_1elec = moints1(ind2val(spstring(bound_elec),spexits(excited_elec)))
  ! Two elecron contribution
    cont_2elec = 0d0
    do i=1, spelec
      if ( i .ne. bound_elec ) then
        cont_2elec = cont_2elec + moints2( index2e2( spstring(i), spstring(i), &
                     spstring(bound_elec),spexits(excited_elec))) - moints2(   &
                     index2e2( spstring(i), spstring(bound_elec), spstring(i), &
                     spexits(excited_elec)))
      end if
    end do
  ! Return
    return
  end subroutine eval_singlex1
!====================================================================
!====================================================================
!>eval_singlex2
!
! Subroutine to compute 3 integral of single excitations
!--------------------------------------------------------------------
  subroutine eval_singlex2( sp2index, sp2elec, orbitals, sp2_dets, sp1orbital_bnd,&
    sp1orbital_ex , moints1, moints2, moints1len, moints2len, cont_2elec)
  !Input:
  ! sp2index       = beta/alpha string index
  ! sp2elec        = beta/alpha electron number
  ! orbitals       = # of MO's
  ! sp2_dets       = non-truncated beta/alpha determinants
  ! sp2_string     = beta/alpha string
  ! sp1orbital_bnd = orbital of electron being excited FROM
  ! sp1orbital_ex  = orbital electron is being excited TO
  ! cont_2elec     = 2-e contribution
    use detci2
    use detci5
    implicit none
    integer, intent(in) :: sp2index, sp2elec, orbitals, sp2_dets, sp1orbital_bnd, &
                           sp1orbital_ex, moints1len, moints2len
    real*8, dimension(moints1len), intent(in) :: moints1
    real*8, dimension(moints2len), intent(in) :: moints2
    real*8, intent(out) :: cont_2elec
    integer, dimension(sp2elec) :: sp2_string
    integer :: i
    !----------------------------------------------------------------
    ! Generate orbital index string
      call genorbstring( sp2index, sp2elec, orbitals, sp2_dets, sp2_string )
    ! Compute contribution
    cont_2elec = 0d0
    do i=1, sp2elec
      cont_2elec = cont_2elec + moints2( index2e2( sp2_string(i), sp2_string(i), &
                   sp1orbital_bnd, sp1orbital_ex ))
    end do
    ! Return
    return
  end subroutine eval_singlex2
!====================================================================
!====================================================================
end module