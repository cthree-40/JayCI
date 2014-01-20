! This module contains subroutines necessary for truncating the 
!  determinant expasion
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chem., The Johns Hopkins University
!
! Last edit: 11-18-13
!====================================================================
module truncation
  use detci2
contains
!====================================================================
!====================================================================
! Subroutine genstrings() generates alpha and beta strings
!--------------------------------------------------------------------
  subroutine genstrings( astrings, bstrings, adets, bdets )
    implicit none
! ...input integer scalars...
    integer, intent(in) :: adets, bdets
! ...input integer arrays...
    integer, dimension(adets) :: astrings
    integer, dimension(bdets) :: bstrings
! ...loop integer scalars...    
    integer :: i    
!--------------------------------------------------------------------
    do i=1, adets
      astrings(i) = i
    end do
    do i=1, bdets
      bstrings(i) = i
    end do
    return
  end subroutine
!====================================================================
!====================================================================
! Subrouinte enffrozen() enforces the frozen orbital restrictions on
!  either alpha or beta strings. Makes det indices not satisfying the
!  restriction 0.
!--------------------------------------------------------------------
  subroutine enffrozen( string, nfrozen, elecs, orbitals, &
                        determs )
! Input:
!  string   = alpha or beta string indices      integer array  1-d
!  nfrozen  = number of frozen orbitals         integer scalar
!  elecs    = number of alpha/beta electrons    integer scalar
!  orbitals = number of MO's                    integer scalar
!  determs  = number of alpha/beta determinants integer scalar
!--------------------------------------------------------------------
    implicit none
! ...input integer scalars...
    integer, intent(in) :: nfrozen, elecs, orbitals, determs
! ...input integer arrays...
    integer, dimension(determs),intent(inout) :: string
! ...loop integer scalars...
    integer :: i, j
! ...integer scalars...
    integer :: test
! ...integer arrays...
    integer, dimension(elecs) :: string1, string2
!--------------------------------------------------------------------
! Loop over determinants
    lpa: do i=1, determs
      if ( i .eq. 1 ) then
        do j=1, elecs
          string1(j) = j
        end do
      else
        call nstrfnd( string2, elecs, orbitals, determs, string1 )
      end if
! Test if the first nfrozen orbitals are occupied
      test = 0
      lpb: do j=1, nfrozen
        if ( string1(j) .eq. j ) then
          cycle lpb
        else
          test = test + 1
        end if
        if ( test .ne. 0 ) then
          string(i) = 0
          exit lpb
        end if
      end do lpb
      string2 = string1
      string1 = 0
    end do lpa
  end subroutine enffrozen
!====================================================================
!====================================================================
! Subroutine to enforce DOCC restriction on either alpha or beta
!  strings. Makes indices that do not satisfy DOCC restrictions 0.
!--------------------------------------------------------------------
  subroutine enfdocc( string, determs, elecs, orbitals, nfrozen, &
                      ndocc, nactive, xlevel )
! Input:
!  string   = alpha/beta string indices         integer array  1-d
!  determs  = alpha/beta determinants           integer scalar
!  elecs    = alpha/beta electrons              integer scalar
!  orbitals = number of MO's                    integer scalar
!  nfrozen  = number of frozen orbitals         integer scalar
!  ndocc    = number of docc orbitals           integer scalar
!  xlevel   = excitation level                  integer scalar
!--------------------------------------------------------------------
    implicit none
! ...input integer scalars...
    integer, intent(in) :: determs, elecs, orbitals, nfrozen, ndocc,&
                           xlevel, nactive
! ...input/output integer arrays...
    integer, dimension(determs), intent(inout) :: string
! ...loop integer scalars...
    integer :: i, j
! ...integer scalars...
    integer :: test
! ...integer arrays...
    integer, dimension(elecs) :: string1, string2
!--------------------------------------------------------------------
! Loop through determinants
    lpa: do i=1, determs
! Test if string is in expansion
      if ( string(i) .eq. 0 ) then
        cycle lpa
      end if
! Construct string for indice i
      if ( string(i) .eq. 1 ) then
        do j=1, elecs
          string1(j) = j
        end do
      else
        call nstrfnd( string2, elecs, orbitals, determs, string1 )
      end if
      
      test=0
! Count excitations in string1
      do j=1, elecs
        if ( string1(j) > nfrozen+ndocc ) then
          test=test+1
        end if
      end do

      if ( (test-nactive) > xlevel) then
        string(i) = 0   ! Remove from expansion
      end if
      
      string2 = string1
      string1 = 0
    end do lpa
    return
  end subroutine enfdocc
!====================================================================
!====================================================================
! Subroutine to enforce CAS restrictions on alpha/beta det. strings
! Makes indices that do not satisfy restrictions 0.
!--------------------------------------------------------------------
  subroutine enfactive( string, determs, orbitals, elecs, nfrozen, &
                        ndocc, nactive, xlevel )
! Input:
!  string   = alpha/beta string indices           integer array  1-d
!  determs  = alpha/beta determinants             integer scalar
!  orbitals = number of MO's                      integer scalar
!  elecs    = number of electrons                 integer scalar
!  nfrozen  = number of frozen orbitals           integer scalar
!  ndocc    = number of DOCC orbitals             integer scalar
!  nactive  = number of CAS orbitals              integer scalar
!  xlevel   = excitation level                    integer scalar
!--------------------------------------------------------------------
    implicit none
! ...input integer scalars...
    integer, intent(in) :: determs, orbitals, elecs, nfrozen, &
                           ndocc, nactive, xlevel
! ...input integer arrays...
    integer, dimension(determs), intent(inout) :: string
! ...loop integer scalars...
    integer :: i, j
! ...integer arrays...
    integer, dimension(elecs) :: string1, string2
! ...integer scalars...
    integer :: test
!--------------------------------------------------------------------
! Loop over determinants
    lpa: do i=1, determs
! Test if string is in expansion
      if ( string(i) .eq. 0 ) then
        cycle lpa
      end if
! Construct string
      if ( string(i) .eq. 1 ) then
        do j=1, elecs
          string1(j) = j
        end do
      else
        call nstrfnd( string2, elecs, orbitals, determs, string1 )
      end if
! Count excitations in string
      test = 0
      do j=1, elecs
        if (string1(j) > nfrozen + ndocc + nactive ) then
          test = test + 1
        end if
      end do
      if ( test > xlevel ) then
        string(i) = 0
      end if
! Close loop
      string2=string1
      string2=0
    end do lpa
    return
  end subroutine enfactive
!====================================================================
!====================================================================
! Subroutine to return number of 'active' alpha/beta strings
!--------------------------------------------------------------------
  subroutine modspstrlen( string, length, modlength )
! Input:
!  string    = alpha/beta string idices              integer array  1-d
!  length    = number of alpha/beta determinants     integer scalar
! Output:
!  modlength = length of active strings              integer scalar
    implicit none
! ...input integer scalars...
    integer, intent(in) :: length
! ...input integer arrays...
    integer, dimension(length), intent(in) :: string
! ...output integer scalars...
    integer, intent(out) :: modlength
! ...loop integer scalars...
    integer :: i, j
!--------------------------------------------------------------------
    modlength = 0
    do i=1, length
      if ( string(i) .ne. 0 ) then
        modlength = modlength + 1
      end if
    end do
    return
  end subroutine modspstrlen
!====================================================================
!====================================================================
! Subroutine to generate list of determinants from the 'active'
!  alpha and beta strings
!--------------------------------------------------------------------
  subroutine gendetlist( alphastring, betastring, alphalen, betalen, &
                         belec, orbitals, determlist )
! Input:
!  alphastring  = 'active' alpha strings            integer array  1-d
!  betastring   = 'active' beta strings             integer array  1-d
!  alphalen     = length of alphastrings()          integer scalar
!  betalen      = length of betastrings()           integer scalar
!  belec        = number of beta electrons          integer scalar
!  orbitals     = number of MO's                    integer scalar
! Output:
!  determlist   = list of determinants              integer array 1-d
!--------------------------------------------------------------------
    implicit none
! ...input integer scalars...
    integer, intent(in) :: alphalen, betalen, belec, orbitals
! ...input integer arrays...
    integer, dimension(alphalen), intent(in) :: alphastring
    integer, dimension(betalen ), intent(in) :: betastring
! ...output integer arrays...
    integer, dimension(alphalen*betalen), intent(inout) :: determlist
! ...loop integer scalars...
    integer :: i, j, k
!--------------------------------------------------------------------
! loop over alpha strings
    do i=1, alphalen
! loop over beta strings
      do j=1, betalen
        k = ( i-1 )*betalen + j
        determlist(k) = indxk( alphastring(i), betastring(j), belec, orbitals )
      end do
    end do
    return
  end subroutine gendetlist
!====================================================================
!====================================================================
! Subroutine to enforce docc restrictions on determinants
!--------------------------------------------------------------------
  subroutine enfdoccdet( detlist, detlistlen, aelec, adets, belec, &
    bdets, orbitals, nfrozen, ndocc, nactive, xlevel, remdet )
! Input:
!  detlist     = list of determinants                integer array  1-d
!  detlistlen  = length of detlist()                 integer scalar
!  aelec       = alpha electrons                     integer scalar
!  adets       = alpha strings ( non-truncated )     integer scalar
!  belec       = beta electrons                      integer scalar
!  bdets       = beta strings ( non-truncated )      integer scalar
!  orbitals    = number of MO's                      integer scalar
!  nfrozen     = number of frozen orbitals           integer scalar
!  ndocc       = number of DOCC orbitals             integer scalar
!  xlevel      = excitation level                    integer scalar
!  remdet      = removed determinants                integer scalar
!--------------------------------------------------------------------
    implicit none
! ...input integer scalars...
    integer, intent(in) :: detlistlen, aelec, belec, adets, bdets, &
                           orbitals, nfrozen, ndocc, nactive, xlevel
! ...input/output integer scalars...
    integer, intent(inout) :: remdet
! ...input/output integer arrays...
    integer, dimension(detlistlen), intent(inout) :: detlist
! ...loop integer scalars...
    integer :: i, j
! ...integer arrays...
    integer, dimension(adets, aelec) :: alphamat
    integer, dimension(bdets, belec) :: betamat
! ...integer scalars...
    integer :: p, q, test
!--------------------------------------------------------------------
! Generate string lists
    call strfnd( adets, aelec, orbitals, adets, alphamat )
    call strfnd( bdets, belec, orbitals, bdets, betamat  )
! Loop over determinants
    do i=1, detlistlen
      call k2indc( detlist(i), belec, orbitals, p, q )
      test=0
! Check alpha string
      do j=nfrozen+1, aelec
        if ( alphamat(p,j) > nfrozen + ndocc ) then
          test = test + 1
        end if
      end do
! Check beta string
      do j= 1, belec
        if ( betamat(q, j) > nfrozen + ndocc ) then
          test = test + 1
        end if
      end do
! If test > xlevel, throw determinant away
      if ( (test-nactive) > xlevel ) then
        detlist(i) = 0
        remdet = remdet + 1
      end if
    end do
    return
  end subroutine enfdoccdet
!====================================================================
!====================================================================
! Subroutine to enforce CAS restrictions on determinants
!--------------------------------------------------------------------
  subroutine enfactivedet( detlist, detlistlen, aelec, adets, belec, &
    bdets, orbitals, nfrozen, ndocc, nactive, xlevel, remdet )
! Input:
!  detlist    = determinant list
!  detlistlen = determinant list length
!  aelec      = alpha electrons
!  adets      = alpha strings (non-truncated)
!  belec      = beta electrons
!  bdets      = beta strings  (non-truncated)
!  orbitals   = number of MO's
!  nfrozen    = frozen orbitals
!  ndocc      = number of DOCC orbitals
!  nactive    = active orbitals
!  xlevel     = excitation level
!  remdet     = removed determinants
!--------------------------------------------------------------------
    implicit none
! ...input integer scalars...
    integer, intent(in) :: detlistlen, aelec, adets, belec, bdets, &
                           orbitals, nfrozen, ndocc, nactive,      &
                           xlevel
! ...input/output integer scalars...
    integer, intent(inout) :: remdet
! ...input/output integer arrays...
    integer, dimension(detlistlen), intent(inout) :: detlist
! ...loop integer scalars...
    integer :: i, j
! ...integer arrays...
    integer, dimension(adets,aelec) :: alphamat
    integer, dimension(bdets,belec) :: betamat 
! ...integer scalars...
    integer :: p, q, test
!--------------------------------------------------------------------
! Generate spin string arrays
    call strfnd( adets, aelec, orbitals, adets, alphamat)
    call strfnd( bdets, belec, orbitals, bdets, betamat )
! Loop over determinants
    do i=1, detlistlen
      if ( detlist(i) .eq. 0 ) then
        cycle
      end if
      call k2indc( detlist(i), belec, orbitals, p, q )
      test=0
! Check alpha string
      do j=1, aelec
        if ( alphamat(p,j) > nfrozen + ndocc + nactive ) then
          test = test + 1
        end if
      end do
! Check beta string
      do j=1, belec
        if ( betamat(q,j) > nfrozen + ndocc + nactive ) then
          test = test + 1
        end if
      end do
! Test if test > xlevel
      if ( test > xlevel ) then
        detlist(i) = 0
        remdet=remdet+1
      end if
    end do
    return
  end subroutine enfactivedet
!====================================================================
!====================================================================
! Subroutine to form string pairs
!--------------------------------------------------------------------
  subroutine genstrpairs( leadspinstr, leadspinstrlen, secondspinstr,  &
    secondspinstrlen, determlist, determlistlen, belec, orbitals,      &
    leadspin, stringpairs, stringpairslen, stringsteps, stringstepslen,&
    stringdets, locate )
! Input:
!  leadspinstr     = alpha/beta strings                integer array  1-d
!  leadspinstrlen  = $spinstr length                   integer scalar
!  secondspinstr   = beta/alpha  strings               integer array  1-d
!  secondspinstrlen= $spinstr  length                  integer scalar
!  determlist      = determinants                      integer array  1-d
!  determlistlen   = length of determlist              integer scalar
!  belec           = beta electrons                    integer scalar
!  orbitals        = orbitals                          integer scalar
!  leadspin        = leading sping    
!  locate          = array of locations of new string  integer array  1-d
! Output:
!  stringpairs     = string pairings                   integer array  2-d
!  stringpairslen  = rows of string pairings           integer scalar
!  stringsteps     = step array for stringpairs        integer array  1-d
!  stringstepslen  = length of stringsteps             integer scalar
!  stringdets      = determinant list from pairs
!--------------------------------------------------------------------
    implicit none
    character*1, intent(in) :: leadspin
! ...input integer scalars...
    integer, intent(in) :: leadspinstrlen, secondspinstrlen, determlistlen,&
                           belec, orbitals, stringpairslen, stringstepslen
! ...input integer arrays...
    integer, dimension( leadspinstrlen    ), intent(in) :: leadspinstr
    integer, dimension( secondspinstrlen  ), intent(in) :: secondspinstr
    integer, dimension( determlistlen     ), intent(in) :: determlist
! ...output integer arrays...
    integer, dimension( stringpairslen, 2 ), intent(inout) :: stringpairs
    integer, dimension( stringstepslen), intent(inout)     :: stringsteps
    integer, dimension( stringstepslen), intent(inout)     :: locate
    integer, dimension( determlistlen ), intent(inout)     :: stringdets
! ...loop integer scalars...
    integer :: i, j, k, l
! ...integer scalars...
    integer :: step, index1
!--------------------------------------------------------------------
    l=1
    locate = 0
! Loop over leading string
    do i=1, leadspinstrlen
      step=0
! Loop over second string
      do j=1, secondspinstrlen
! Find the index of the determinant composed of the two strings
        if ( leadspin .eq. 'a' ) then
          index1 = indxk( leadspinstr(i), secondspinstr(j), belec, orbitals )
        else if ( leadspin .eq. 'b' ) then
          index1 = indxk( secondspinstr(j), leadspinstr(i), belec, orbitals )
        else
          STOP " NO SPIN SPECFIED. EXITING."
        end if
! Test if index1 is in expansion
        do k=1, determlistlen
          if ( index1 .eq. determlist(k) ) then
            stringpairs(l,1) = i                 ! START FROM PREVIOUS INDEX
            stringpairs(l,2) = j
            stringdets(l) = index1
            l=l+1
            step=step+1
            exit
          end if
        end do
      end do
      stringsteps(i) = step
      do j=1,i-1
        locate(i) = locate(i) + stringsteps(j)
      end do
    end do
    return
  end subroutine genstrpairs
!====================================================================
!====================================================================
! Subroutine to generate the determinant cross refernce list
!--------------------------------------------------------------------
  subroutine detcrossref( pdeterms, qdeterms, determs, length, &
                          pcrossreflist, qcrossreflist ) 
! Input:
!  pdeterms     = determinants ordered by p           integer array  1-d
!  qdeterms     = determinants ordered by q           integer array  1-d
!  length       = number of determinants              integer scalar
! Output:
!  crossreflist = cross reference list                integer scalar 1-d
!--------------------------------------------------------------------
    implicit none
! ...input integer scalars...
    integer, intent(in) :: length
! ...input integer arrays...
    integer, dimension( length ), intent(in) :: pdeterms, qdeterms, determs
! ...output itneger arrays...
    integer, dimension( length ), intent(inout) :: pcrossreflist, qcrossreflist
! ...loop integer scalars...
    integer :: i, j
! ...integer scalars...
    integer :: l
!--------------------------------------------------------------------
! Loop over p list
    l = 1
    do i=1, length
      do j=1, length
        if ( pdeterms(i) .eq. determs(j) ) then
          pcrossreflist(l) = j
          l=l+1
          exit
        end if
      end do
    end do
 
! Loop over q list
    l = 1
    do i=1, length
      do j=1, length
        if ( qdeterms(i) .eq. determs(j) ) then
          qcrossreflist(l) = j
          l=l+1
          exit
        end if
      end do
    end do
    return
  end subroutine detcrossref
!====================================================================
!====================================================================
!>gen_alphastrpr
!
! Subroutine to generate alpha string paris
!--------------------------------------------------------------------
  subroutine gen_alphastrpr( determlist, determlistlen, belec, orbitals, &
    alpha_strpair, alpha_step, alpha_locate, alpha_detlen )
    implicit none
    integer, intent(in)  :: determlistlen, belec, orbitals, alpha_detlen
    integer, dimension( determlistlen ), intent(in)  :: determlist
    integer, dimension( determlistlen,2 ), intent(out) :: alpha_strpair 
    integer, dimension( alpha_detlen ),  intent(out) :: alpha_step, alpha_locate 
    integer :: l, i, j, step, locate, p, q, tot_step
  !--------------------------------------------------------------------
  ! Loop through determinants generating alpha_strpair
    alpha_step = 1   ! All determinants will have at least one step
    do i=1, determlistlen
      call k2indc( determlist(i), belec, orbitals, p, q )
      alpha_strpair(i,1) = p
      alpha_strpair(i,2) = q
    end do
  ! Loop through alpha_strpair generating step array and location array
    tot_step=0
    step = 1
    j=1
    l=1
    do i=1, determlistlen-1
      if ( alpha_strpair(i,1) .eq. alpha_strpair(i+1,1) ) then
        step = step + 1
      else
        alpha_step(l) = step
        tot_step = tot_step + step
        step = 1
        l = l+1
        if ( l .eq. alpha_detlen ) then
          alpha_step(l) = determlistlen - tot_step
          exit
        end if
      end if
    end do
    do i=1, alpha_detlen
      locate = 0
      do j=1, i-1
        locate = alpha_step(j) + locate
      end do
      alpha_locate(i) = locate
    end do
    return
  end subroutine gen_alphastrpr
!======================================================================
!======================================================================
!>gen_betastrpr
!
! Subroutine to generate beta string pairings
!----------------------------------------------------------------------
  subroutine gen_betastrpr( alpha_strpr, cidimension, alpha_step,  &
    alpha_detlen, alpha_locate, beta_strpr, beta_step, beta_locate,&
    beta_det, beta_detlen )        
    implicit none
    integer, intent(in) :: cidimension, alpha_detlen, beta_detlen
    integer, dimension( cidimension, 2 ), intent(in) :: alpha_strpr
    integer, dimension( alpha_detlen ),   intent(in) :: alpha_step, alpha_locate
    integer, dimension( beta_detlen ),    intent(in) :: beta_det
    integer, dimension( beta_detlen ),   intent(out) :: beta_step, beta_locate
    integer, dimension( cidimension, 2 ),intent(out) :: beta_strpr
    integer :: i, j, step, locate, l, tot_step
  !--------------------------------------------------------------------
  ! Loop over q strings
    beta_step = 1 ! All strings will have at least one step
    l=1
    do i=1, beta_detlen
      ! Loop over p string pairings
      do j=1, cidimension
        if ( alpha_strpr(j,2) .eq. beta_det(i) ) then
          beta_strpr(l,1) = alpha_strpr(j,2)
          beta_strpr(l,2) = alpha_strpr(j,1)
          l=l+1
        end if
      end do
    end do
  ! Generate step array and location array1
    step = 1
    j=1
    l = 1
    tot_step = 0
    do i=1, cidimension-1
      if ( beta_strpr(i,1) .eq. beta_strpr(i+1, 1) ) then
        step = step + 1
      else
        beta_step(l) = step
        tot_step = tot_step + step
        step=1
        l = l+1
        if ( l .eq. beta_detlen ) then
          beta_step(l) = cidimension - tot_step
          exit
        end if
      end if
    end do
    do i=1, beta_detlen
      locate = 0
      do j=1, i-1
        locate = beta_step(j) + locate
      end do
      beta_locate(i) = locate
    end do
    return
  end subroutine gen_betastrpr
!======================================================================
!======================================================================
!>xreflist_gen
!
! Generate cross reference list ! Scales like n!
!----------------------------------------------------------------------
  subroutine xreflist_gen( length, detlist1, qpairings, belec, &
    orbitals, xreflist )
    use detci2, only: indxk
    implicit none
    integer, intent(in) :: length, belec, orbitals
    integer, dimension( length ),    intent(in)  :: detlist1
    integer, dimension( length, 2),  intent(in)  :: qpairings
    integer, dimension( length ),    intent(out) :: xreflist
    integer, dimension(:,:),           allocatable :: temp
    integer, dimension(:),             allocatable :: detlist2
    integer :: k, i, j, trunc, element
!----------------------------------------------------------------------
    ! Use string pairings to generate det list2
    allocate( detlist2( length ) )
    do i=1, length
      detlist2(i) = indxk( qpairings(i,2), qpairings(i,1), belec, orbitals )
    end do
    
    ! Generate temp list
    allocate( temp( length, 2 ) )
    do i=1, length
      temp(i,1) = detlist1(i)
      temp(i,2) = i
    end do
    trunc = 0
    element = 1
   
    ! Loop over detlist2
    loopa: do i=1, length
      ! Loop over detlist1
      loopb: do j=1, length-trunc
        if ( detlist2(i) .eq. temp(j,1) ) then
          xreflist(element) = temp(j,2)
          element = element + 1
          trunc = trunc + 1
          ! Move elements down
          do k=j+1, length
            temp(k-1,1) = temp(k,1)
            temp(k-1,2) = temp(k,2)
          end do
          exit loopb
        end if
      end do loopb
    end do loopa

    deallocate(temp)

    return
  end subroutine xreflist_gen   
          
    
end module
