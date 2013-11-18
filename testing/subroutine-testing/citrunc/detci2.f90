module detci2
!===============================================================
! This module contains subroutines for addresssing the
! determinant strings
!
! Subroutines:
!  strfnd
!  nstrfnd
!  adrfnd
!  confadr
!  binadd
!  indxk
!  k2indc
!  k2indcp
!  k2indcq
!---------------------------------------------------------------
! Christopher L. Malbon
! Yarkony Group
! Dept. Chemistry, The Johns Hopkins University
!
! Last edit:
!  09-27-13 - 
!===============================================================
  implicit none

contains
!---------------------------------------------------------------
! strfnd( address, num_elec, num_orb, num_det, A )
! 
! locates the string associated with an address. This uses the 
!  algorithm developed by Ruedenberg and Ivanic:
!   Theo. Chem. Acc. (2001) 106:339-351
!---------------------------------------------------------------
  subroutine strfnd( address, num_elec, num_orb, num_det, A )
!----------------------------
! Input:
!  address = address of string( alpha or beta )
!  num_elec = number of electrons ( alpha or beta )
!  num_orb = number of orbitals
!  num_det = number of total determinants (not truncated)
! Output:
!  A = array A(p,j) - orbital occupied by jth electron in pth string
!----------------------------
    implicit none
    integer,intent(in) :: address, num_elec, num_orb, num_det
    integer,dimension(num_det, num_elec) :: A

    integer :: k, l, i, j, p, q
!----------------------------
!
! Initialize by setting A(1,j) = 1, 2, ..., num_elec
    do j=1, num_elec
      A(1,j) = j
    end do
!
    a_loop: do p=1, (num_det - 1) !Loop over all strings p
      b_loop: do i=0, (num_elec - 1)
 
        if ( A(p,(num_elec - i)) == (num_orb - i)) then
          cycle !Advance i -> i + 1
        
        else !Determine A(p+1,j)
          do j=1, (num_elec - i - 1)
            A((p+1),j) = A(p,j)
          end do
        
          do j=(num_elec - i), num_elec
            A((p+1),j) = A(p, (num_elec - i)) + 1 + j - (num_elec - i)
          end do
        
        end if
        
        exit b_loop
      
      end do b_loop
    
    end do a_loop
! Return desired string
    return
  end subroutine strfnd
!
!---------------------------------------------------------------
! nstrfnd( oldstring, num_elec, num_orbitals, num_det. newstring )
! 
! locate string given index and old string
!---------------------------------------------------------------
  subroutine nstrfnd( oldstring, num_elec, num_orbitals, &
    num_det, newstring )

    implicit none

    integer, intent(in) :: num_elec, num_orbitals, num_det

    integer, dimension( num_elec ) :: oldstring

    integer, dimension( num_elec ) :: newstring

    integer :: i, j, k

!----------------------------
! Loop over electrons i
    do i=0, num_elec - 1
      if ( oldstring( num_elec - i ) == ( num_orbitals - i ) ) then
        cycle ! Advance i -> i+1
      else ! Determine newstring
        do j=1, (num_elec - i - 1)
          newstring(j) = oldstring(j)
        end do
        do j=(num_elec - i), num_elec
          newstring(j) = oldstring( (num_elec - i)) + 1 + j - (num_elec - i)
        end do
      end if
      exit
    end do
    return
  end subroutine
! 
!---------------------------------------------------------------
! adrfnd( detstr, elects, orbits, addrs )
!
! computes the address of a string of orbital indices. This uses
!  the algorithm developed by Ruedenberg and Ivanic:
!  Theo. Chem. Acc. (2001) 106: 339-351
!---------------------------------------------------------------
  subroutine adrfnd( detstr, elects, orbits, addrss )
!----------------------------
! Input:
!  detstr = determinant string
!  elects = number of electrons (alpha / beta)
!  orbits = nubmer of orbitals
! Output:
!  addrss = address
!----------------------------
    implicit none
    integer, intent(in) :: elects, orbits
    integer, dimension(elects), intent(in) :: detstr
    integer :: addrss

    integer :: i, j, sum1
!----------------------------
!
    sum1 = 1

! Sum over all electrons
    do i=1, elects
      j = i - 1

      if ( j == 0 ) then ! For the first electron
        sum1 = sum1 + confadr(i, detstr(i), j, elects, orbits )
      else
        sum1 = sum1 + confadr(i, detstr(i), detstr(j), elects, orbits)
      end if
    
    end do
    addrss = sum1
    return
  end subroutine adrfnd

!----------------------------------------------------------
! confadr( electn, orbitl, preorb, numelc, numorb )
!
! computes numerical value of each electron's location
!----------------------------------------------------------
  integer function confadr( electn, orbitl, preorb, numelc, &
    numorb )
!----------------------------
! Input:
!  electn    = electron we are summing over
!  orbitl        = orbital the electron occupies
!  preorb    = orbital occupied by (electn - 1) electron
!  numelc    = number of electrons (alpha / beta)
!  numorb    = number of orbitals
! Output:
!  cnfiga    = value of the particular configuration
!----------------------------
    implicit none
    integer, intent(in) :: electn, orbitl, preorb, numelc, &
      numorb

    integer :: sum1, sumfrm, sumuto, j
!----------------------------
!
    sumfrm = preorb + 1
    sumuto = orbitl - 1
    sum1 = 0

    do j=sumfrm, sumuto
      sum1 = sum1 + binadd( numelc, numorb, electn, j)
    end do

    confadr = sum1
  end function confadr

!----------------------------------------------------------
! binadd(numels, orbtls, electn, orbind )
!
! used by confadr to evaluate each binomial
!----------------------------------------------------------
  integer function binadd( numels, orbtls, electn, orbind )
!----------------------------
! Input:
!  numels  = number of electrons (alpha / beta)
!  orbtls  = number of orbitals
!  electn  = electron we are summing over
!  orbind  = summation index from config_address
! Output:
!  binadd  = binomial for config
!----------------------------
    use detci1, only: binom
!----------------------------
    implicit none
    integer, intent(in) :: numels, orbtls, electn, orbind

    integer :: j, k
!----------------------------
!
    j = numels - electn
    k = orbtls - orbind
    binadd = binom(k, j)
  end function binadd

!----------------------------------------------------------
! indxk( p, q, elec, orbs )
!
! computes the index of the determinant K=|p,q>
!----------------------------------------------------------
  integer function indxk( p, q, belec, orbs )
!----------------------------
! Input:
!  p = index of alpha string
!  q = index of beta string
!  belec = beta electrons
!  orbs = orbitals
!----------------------------
     use detci1, only: binom
!----------------------------
     implicit none
     integer,intent(in) :: p, q, belec, orbs
!----------------------------
!
     indxk = (p - 1)*binom(orbs, belec) + q
   end function indxk

!----------------------------------------------------------
! k2indc(k, aelec, belec, orbs, p, q )
!
! returns the value of the indices of alpha and beta strings
!  given the index of K = |p, q>
!----------------------------------------------------------
  subroutine k2indc( k, belec, orbs, p, q )
!----------------------------
! Input: 
!  k = index of | p, q >
!  belec = number of beta electrons
!  orbs = number of orbitals
! Output:
!  p = alpha string
!  q = beta string
!----------------------------
    implicit none
    integer,intent(in) :: k, belec, orbs
    integer,intent(out) :: p, q
!----------------------------
!

!#ifdef K2_IND
!  print *, "Debugging k2ind"
!  write(*,fmt=10) "Number of beta electrons =", belec
!10 format(1x, A, I5)
!#endif


    p = k2indcp(k, belec, orbs)
    q = k2indcq(k, belec, orbs, p)

  end subroutine k2indc

!----------------------------------------------------------
! k2indcp( k, belec, orbs )
!
! computes alpha string index. used by k2indc (above)
!----------------------------------------------------------
  integer function k2indcp(k, belec, orbs)
!----------------------------
! Input:
!  k = determinant index
!  belec = beta electrons
!  orbs = orbitals
!----------------------------
    use detci1, only: binom
!----------------------------
    implicit none
    integer, intent(in) :: k, belec, orbs
    integer :: test
    real*8 :: frac
!----------------------------
!
    frac = dble(k) / dble(binom(orbs,belec))
    test = ceiling( frac )
    k2indcp = test

  end function k2indcp

!----------------------------------------------------------
! k2indcq( k, belec, orbs, p )
!
! computes beta string index. used by k2indc (above)
!----------------------------------------------------------
  integer function k2indcq( k, belec, orbs, p)
!----------------------------
! Input:
!  k = determinant index
!  belec = beta electrons
!  orbs = orbitals
!  p = alpha string index
!----------------------------
    use detci1, only: binom
!----------------------------
    implicit none
    integer, intent(in) :: k, belec, orbs, p
!----------------------------
!
    k2indcq = k - (p - 1)*binom(orbs,belec)

  end function k2indcq

!----------------------------------------------------------
!==========================================================
end module  
     











