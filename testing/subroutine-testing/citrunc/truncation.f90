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
                      ndocc, xlevel )
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
                           xlevel
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

      if ( test > xlevel ) then
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
    bdets, orbitals, nfrozen, ndocc, xlevel, remdet )
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
                           orbitals, nfrozen, ndocc, xlevel
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
      do j=1, aelec
        if ( alphamat(p,j) > nfrozen + ndocc ) then
          test = test + 1
        end if
      end do
! Check beta string
      do j=1, belec
        if ( betamat(q, j) > nfrozen + ndocc ) then
          test = test + 1
        end if
      end do
! If test > xlevel, throw determinant away
      if ( test > xlevel ) then
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
    stringdets )
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
    integer, dimension( determlistlen ), intent(inout)     :: stringdets
! ...loop integer scalars...
    integer :: i, j, k, l
! ...integer scalars...
    integer :: step, index1
!--------------------------------------------------------------------
    l=1
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
            stringpairs(l,1) = i
            stringpairs(l,2) = j
            stringdets(l) = index1
            l=l+1
            step=step+1
            exit
          end if
        end do
      end do
      stringsteps(i) = step
    end do
    return
  end subroutine genstrpairs
!====================================================================
!====================================================================
! Subroutine to generate the determinant cross refernce list
!--------------------------------------------------------------------
  subroutine detcrossref( pdeterms, qdeterms, length, crossreflist ) 
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
    integer, dimension( length ), intent(in) :: pdeterms, qdeterms
! ...output itneger arrays...
    integer, dimension( length ), intent(inout) :: crossreflist
! ...loop integer scalars...
    integer :: i, j
! ...integer scalars...
    integer :: l
!--------------------------------------------------------------------
! Loop over q list
    l = 1
    do i=1, length
      do j=1, length
        if ( qdeterms(i) .eq. pdeterms(j) ) then
          crossreflist(l) = j
          l=l+1
          exit
        end if
      end do
    end do
    return
  end subroutine detcrossref
!====================================================================
end module
