module detci1
!===============================================================
! This module computes binomials and factorials
! 
! Christopher L. Malbon
! Yarkony Group
! Dept. of Chemistry, The Johns Hopkins University
! 
! Last edit:
!  09-27-13 - removal of subroutine LWVAL
!===============================================================
  implicit none
  public :: factl, binom

contains
!---------------------------------------------------------------
! fact(n) returns the factorial of n
!---------------------------------------------------------------
  integer function factl(n)
!----------------------------
    implicit none
    integer, intent(in) :: n
    
    integer :: i
!----------------------------
!
    factl = product( (/ (i, i=1, n) /))
  end function factl

!---------------------------------------------------------------
! binom(n,k)
! 
! computes the binomial coefficient of n and k
!---------------------------------------------------------------
  integer function binom(n,k)
    implicit none
    integer, intent(in) :: n, k
    integer :: i, j, w, q
    integer, dimension(n-k) :: a, d
    integer, dimension(k) :: b, c
!----------------------------
!
    if ( (n-k) < k ) then
! Construct numerator and denominator
      do i=1, (n-k)
        a(i) = i      ! Denominator
        d(i) = k + 1  ! Numerator
      end do

! Begin with the largest element of string(s).  If element of
!  denominator divides element of numerator, then the element
!  is replaced by the quotient
      do w=(n-k), 2, -1
        do q=(n-k), 1, -1
          if ( mod( d(q),a(w) ) == 0 ) then
            d(q) = d(q) / a(w)
            a(w) = 1
          end if
        end do
      end do
      
      binom = product(d) / product(a)

    else if ( (n-k) .ge. k ) then
! Construction of numerator and denominator
      do i=1, k
        b(i) = (n-k) + i  !Numerator
        c(i) = i          !Denominator
      end do

! Same as above
      do w=k, 2, -1
        do q=k, 1, -1
          if ( mod( b(q), c(w) ) == 0 ) then
            b(q) = b(q) / c(w)
            c(w) = 1
          end if
        end do
      end do

      binom = product(b) / product(c)

    end if
  end function binom

!---------------------------------------------------------------
! alphbet( elec, alpha, beta )
!
! returns the number of alpha and beta electrons based upon the
!  number of electrons in the system
! alpha > beta
!---------------------------------------------------------------
  subroutine alphbet( elec, alpha, beta )
    implicit none
    integer, intent(in) :: elec
    integer, intent(out) :: beta, alpha
!----------------------------
!
    if ( mod(elec, 2) == 0 ) then
      alpha = elec / 2
      beta  = elec / 2
    else
      alpha = ( elec + 1 ) / 2
      beta  = ( elec - 1 ) / 2
    end if
    return
  end subroutine alphbet

!---------------------------------------------------------------
!===============================================================
end module
