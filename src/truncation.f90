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
end module
