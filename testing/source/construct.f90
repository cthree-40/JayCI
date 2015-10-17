module construct
  implicit none
contains
!--------------------------------------------------------------------
  subroutine orbdiffs( pstring1, pstring2, qstring1, qstring2, aelec, belec, &
    diffs )
    implicit none
    integer, intent(in)  :: aelec, belec
    integer, dimension( aelec ), intent(in) :: pstring1, pstring2
    integer, dimension( belec ), intent(in) :: qstring1, qstring2
    integer, intent(out) :: diffs

    integer :: pstdiff, qstdiff
!--------------------------------------------------------------------
! Compare pstrings
    call compstrings( pstring1, pstring2, aelec, pstdiff )
    call compstrings( qstring1, qstring2, belec, qstdiff )
! Compute diffs
    diffs = pstdiff+qstdiff
    return
  end subroutine
!====================================================================
!====================================================================
!> compstrings
!
! Subroutine to compute differences in strings
!--------------------------------------------------------------------
  subroutine compstrings( string1, string2, length, diff )
    implicit none
    integer, intent(in) :: length
    integer, dimension( length), intent(in) :: string1, string2
    integer, intent(out):: diff
    integer :: i, j, test
!--------------------------------------------------------------------
! Loop through strings
    diff = 0
    loop1: do i=1, length
      test = 0
      loop2: do j=1, length
        if ( string1(i) .eq. string2(j) ) then
          test = test + 1
        end if
      end do loop2
      if ( test .eq. 0 ) then
        diff = diff + 1
      end if
    end do loop1
    return
  end subroutine
!====================================================================
!====================================================================
!> dblexcitations
!
! real*8 function to compute matrix elements between determinants with
!  2 diff orbitals
!--------------------------------------------------------------------
  real*8 function dblexcitations( pstring1, pstring2, qstring1, qstring2, &
    aelec, belec, moints1, moints1len, moints2, moints2len )
    use integral
    implicit none
    integer, intent(in) :: aelec, belec, moints1len, moints2len
    integer, dimension( aelec ),     intent(in) :: pstring1, pstring2
    integer, dimension( belec ),     intent(in) :: qstring1, qstring2
    real*8, dimension( moints1len ), intent(in) :: moints1
    real*8, dimension( moints2len ), intent(in) :: moints2
    integer, dimension(:,:),ALLOCATABLE :: pdiffs, qdiffs
    integer :: pd, qd, PermInd1, PermInd2
!--------------------------------------------------------------------
    ALLOCATE(pdiffs(2,2))
    ALLOCATE(qdiffs(2,2))
! Find orbitals differing in alpha strings
    call stringdiffs( pstring1, pstring2, aelec, pdiffs, pd, PermInd1 )
! Find orbitals differing in beta strings
    call stringdiffs( qstring1, qstring2, belec, qdiffs, qd, PermInd2 )
    if ( pd .eq. 2 ) then
            dblexcitations = PermInd1*(moints2( Index2E( pstring1(pdiffs(1,1)), pstring2(pdiffs(1,2)),&
                              pstring1(pdiffs(2,1)), pstring2(pdiffs(2,2)))) - &
                       moints2( Index2E( pstring1(pdiffs(1,1)), pstring2(pdiffs(2,2)),&
                              pstring1(pdiffs(2,1)), pstring2(pdiffs(1,2)))))
    else if ( qd .eq. 2 ) then
            dblexcitations =PermInd2*( moints2( Index2E( qstring1(qdiffs(1,1)), qstring2(qdiffs(1,2)),&
                                qstring1(qdiffs(2,1)), qstring2(qdiffs(2,2)))) - &
                       moints2( Index2E( qstring1(qdiffs(1,1)), qstring2(qdiffs(2,2)),&
                              qstring1(qdiffs(2,1)), qstring2(qdiffs(1,2)))))
    else
            dblexcitations = PermInd1*PermInd2*( moints2( Index2E( pstring1(pdiffs(1,1)),pstring2(pdiffs(1,2)), &
                                qstring1(qdiffs(1,1)), qstring2(qdiffs(1,2)))))
    end if
    DEALLOCATE(pdiffs,qdiffs)
  end function
!====================================================================
!====================================================================
!> stringdiffs
!
! Subroutine to find difference between two strings and return the 
!  differing orbitals
!--------------------------------------------------------------------
  subroutine stringdiffs( string1, string2, length, diff_mat, diff_num, PermInd )
    use string_util
    implicit none
    integer, intent(in)                       :: length
    integer, dimension( length ), intent(in)  :: string1, string2
    integer, dimension( length )              :: TempString
    integer, dimension(2,2), intent(out)      :: diff_mat
    integer, intent(out)                      :: diff_num , PermInd
    integer                                   :: i, j, l
    integer                                   :: test, eps1, esp2, DUMMY
!--------------------------------------------------------------------
    diff_num = 0
    diff_mat = 0
    l = 1
    loopa: do i=1, length
      test = 0
      loopb: do j=1, length
        if ( string1(i) .eq. string2(j) ) then
          test = test + 1
        end if
      end do loopb
      if ( test .eq. 0 ) then
        diff_mat(l,1) = i
        l = l + 1
        diff_num = diff_num + 1
      end if
    end do loopa
    l=1
    do i=1, length
      test = 0
      do j=1, length
        if ( string2(i) .eq. string1(j) ) then
          test = test + 1
        end if
      end do
      if ( test .eq. 0 ) then
        diff_mat(l,2) = i
        l=l+1
      end if
    end do
    ! Replace orbital with excitation to get permutation index
    TempString = string1

    if ( diff_num .eq. 1 ) then
            TempString(diff_mat(1,1)) = string2(diff_mat(1,2))
            CALL cannon1swp(TempString,length,diff_mat(1,1),PermInd)
            !call singrepinfo( TempString, length, string2( diff_mat(1,2) ), &
            !        diff_mat(1,1), 0, PermInd, DUMMY ) 
    else if ( diff_num .eq. 2 ) then ! Should not be necessary to test again

            CALL cannon2(TempString,length,diff_mat(1,1),string2(diff_mat(1,2)),&
                  diff_mat(2,1),string2(diff_mat(2,2)),PermInd)
            !call doublerepinfo( TempString, length, string2(diff_mat(1,2) ),&
            !        diff_mat(1,1), string2(diff_mat(2,2)), diff_mat(2,1), &
            !        0, PermInd, DUMMY )
    end if
    
    return
  end subroutine
!====================================================================
!====================================================================
!> singlexcitations
!
! real*8 function to compute single excitations
!--------------------------------------------------------------------
  real*8 function singlexcitations( pstring1, pstring2, qstring1, qstring2, &
    aelec, belec, moints1, moints1len, moints2, moints2len, ind1, ind2 )
    use integral
    implicit none
    integer, intent(in)   :: aelec, belec, moints1len, moints2len, ind1, ind2
    integer, dimension( aelec ),    intent(in) :: pstring1, pstring2     
    integer, dimension( belec ),    intent(in) :: qstring1, qstring2
    real*8,  dimension(moints1len), intent(in) :: moints1
    real*8,  dimension(moints2len), intent(in) :: moints2
    real*8 :: val
    integer :: i
    integer :: pd, qd
    integer :: PermInd1, PermInd2
    integer, dimension( :, : ),ALLOCATABLE :: pdiffs, qdiffs
!--------------------------------------------------------------------
    ALLOCATE(qdiffs(2,2))
    ALLOCATE(pdiffs(2,2))
    qdiffs = 0
    pdiffs = 0
! Find orbitals differing in alpha strings
    call stringdiffs( pstring1, pstring2, aelec, pdiffs, pd , PermInd1)
! Find orbitals differing in beta strings
    call stringdiffs( qstring1, qstring2, belec, qdiffs, qd , PermInd2)
! Test which string has excitation
    val = 0d0
    if ( pd .eq. 1 ) then
      val = moints1( ind2val(pstring1(pdiffs(1,1)), pstring2(pdiffs(1,2))))
      do i=1, aelec
        if ( pstring1(i) .ne. pstring1(pdiffs(1,1)) ) then
          val = val + moints2( Index2E( pstring1(i), pstring1(i),           &
                      pstring1(pdiffs(1,1)), pstring2(pdiffs(1,2)) ) ) -     &
                      moints2( Index2E( pstring1(i), pstring1(pdiffs(1,1)), &
                      pstring1(i), pstring2(pdiffs(1,2))))
        end if
      end do
      do i=1, belec
        val = val + moints2( Index2E( qstring1(i), qstring1(i), &
                    pstring1(pdiffs(1,1)), pstring2(pdiffs(1,2)) ) ) !-       &
      end do

#ifdef DEBUGGING
!      write(*, 11) qstring1, qstring2, pstring1, pstring2, PermInd1, pdiffs(1,1:2)
!11    format(1x,/,20('-'),/,'b1: ',4i3,/,'b2: ',4i3,/,'a1: ',4i3,/,'a2: ',4i3,/,'p: ',i3,/, &
!        'pdiffs:',2(2i2,/))
#endif
      val = PermInd1*val
    else if ( qd .eq. 1 ) then
      val = moints1( ind2val(qstring1(qdiffs(1,1)), qstring2(qdiffs(1,2))))
      do i=1, belec
        if ( qstring1(i) .ne. qstring1(qdiffs(1,1)) ) then
          val = val + moints2( Index2E( qstring1(i), qstring1(i),           &
                      qstring1(qdiffs(1,1)), qstring2(qdiffs(1,2)))) -       &
                      moints2( Index2E( qstring1(i), qstring1(qdiffs(1,1)), &
                      qstring1(i), qstring2(qdiffs(1,2))))
        end if
      end do
      do i=1, aelec
        val = val + moints2( Index2E( pstring1(i), pstring1(i), &
                  qstring1(qdiffs(1,1)), qstring2(qdiffs(1,2)) ) ) !-       &
      end do
      val = PermInd2*val
    end if
    singlexcitations = val
    DEALLOCATE(pdiffs,qdiffs)
  end function
!====================================================================
!====================================================================
!>
!
!
!--------------------------------------------------------------------
  real*8 function Diagonal_Mat( pstring1, qstring1,  aelec, &
                       belec, moints1, moints1len, moints2, moints2len )
    use integral
    implicit none
    integer, intent(in) :: belec, aelec, moints1len, moints2len
    integer, dimension( aelec ),    intent(in) :: pstring1
    integer, dimension( belec ),    intent(in) :: qstring1
    real*8,  dimension(moints1len), intent(in) :: moints1
    real*8,  dimension(moints2len), intent(in) :: moints2
    integer :: i, j
    real*8  :: val
!--------------------------------------------------------------------
! Loop over alpha string
    val = 0d0
! 1-e contribution
    do i=1, aelec
      val = val + moints1(ind2val(pstring1(i),pstring1(i)))
    end do
! 2-e contribution
    do i=1, aelec
      do j=1, i-1
        val = val + moints2( Index2E( pstring1(i), pstring1(i), &
                    pstring1(j), pstring1(j))) -                 &
                    moints2( Index2E( pstring1(i), pstring1(j), &
                    pstring1(i), pstring1(j)))
      end do
    end do
! Loop over beta string
! 1-e contribution
    do i=1, belec
      val = val + moints1(ind2val(qstring1(i), qstring1(i)))
    end do
! 2-e contribution
    do i=1, belec
      do j=1, i-1
        val = val + moints2( Index2E( qstring1(i), qstring1(i), &
                    qstring1(j), qstring1(j))) -                 &
                    moints2( Index2E( qstring1(i), qstring1(j), &
                    qstring1(i), qstring1(j)))
      end do
    end do
! Alpha and beta 2-e contribution
    do i=1, aelec
      do j=1, belec
        val = val + moints2( Index2E( pstring1(i), pstring1(i), &
                    qstring1(j), qstring1(j)))
      end do
    end do
    diagonal_mat = val
  end function
!====================================================================
!====================================================================
!====================================================================
!> exp_construct
!
! Subroutine to explicitly construct H by finding value of each matrix
!  element.
!--------------------------------------------------------------------
  subroutine exp_construct( moints1, moints1len, moints2, moints2len, & 
    cidim, aelec, belec, orbitals, determs, hamiltonian )
    implicit none
    integer, intent(in) :: moints1len, moints2len, cidim, aelec, belec, &
                         orbitals
    real*8,  dimension( moints1len ), intent(in) :: moints1
    real*8,  dimension( moints2len ), intent(in) :: moints2
    integer, dimension( cidim ), intent(in)      :: determs
    real*8,  dimension( cidim, cidim ), intent(out) :: hamiltonian
    integer :: i, j

!--------------------------------------------------------------------
! Construct hamiltonian
    do i=1, cidim
      do j=1, cidim
        hamiltonian(j,i) = ham_element( determs(j), determs(i), moints1,    &
                           moints1len, moints2, moints2len, aelec,   &
                           belec, orbitals)
      end do
    end do
    return
  end subroutine exp_construct
!======================================================================
!======================================================================
!> ham_element
!
! real*8 function to compute hamiltonian element(i,j)
!--------------------------------------------------------------------
  real*8 function ham_element( ind1, ind2, moints1, moints1len, moints2, &
    moints2len, aelec, belec, orbitals)
    use combinatorial,  only: Binom
    use addressing,     only: K2Indc,GenOrbString
    implicit none
    integer, intent(in) :: ind1, ind2, moints1len, moints2len, &
                           aelec, belec, orbitals
    real*8, dimension( moints1len ), intent(in) :: moints1
    real*8, dimension( moints2len ), intent(in) :: moints2
    integer :: p1, q1, p2, q2
    integer, dimension( aelec ) :: pstring1, pstring2
    integer, dimension( belec ) :: qstring1, qstring2

    integer :: diffs, adets, bdets
!--------------------------------------------------------------------
    adets = Binom( orbitals, aelec )
    bdets = Binom( orbitals, belec )
! Find determinant string indices for ind1 and ind2
    call K2Indc( ind1, belec, orbitals, p1, q1 )
    call K2Indc( ind2, belec, orbitals, p2, q2 )
! Find respective strings for p1, q1, p2, q2
    call GenOrbString( p1, aelec, orbitals, adets, pstring1 )
    call GenOrbString( p2, aelec, orbitals, adets, pstring2 )
    call GenOrbString( q1, belec, orbitals, bdets, qstring1 )
    call GenOrbString( q2, belec, orbitals, bdets, qstring2 )
! Test differences in strings. If > 2 orbitals, element is 0
    call orbdiffs( pstring1, pstring2, qstring1, qstring2, aelec, belec, &
                   diffs )
    if ( diffs .gt. 2 ) then
      ham_element = 0d0
    else if ( diffs .eq. 2 ) then
      ham_element = dblexcitations( pstring1, pstring2, qstring1, qstring2, aelec, &
                           belec, moints1, moints1len, moints2, moints2len )
    else if ( diffs .eq. 1 ) then
      ham_element =  singlexcitations( pstring1, pstring2, qstring1, qstring2, aelec, &
                           belec, moints1, moints1len, moints2, moints2len, ind1, ind2 )
    else 
      ham_element =  diagonal_mat( pstring1, qstring1, aelec, &
                           belec, moints1, moints1len, moints2, moints2len )
    end if
  end function
!======================================================================
!======================================================================

!====================================================================
!====================================================================
!>ham_element_diag
!
! Subroutine to compute diagonal elements...see above for input details
!--------------------------------------------------------------------
  real*8 function ham_element_diag( ind1, moints1, moints1len, moints2, &
    moints2len, aelec, belec, orbitals)
    use combinatorial,   only: Binom
    use addressing,     only: GenOrbString,K2Indc
    implicit none
    integer, intent(in) :: ind1, moints1len, moints2len, &
                         aelec, belec, orbitals
    real*8, dimension( moints1len ), intent(in) :: moints1
    real*8, dimension( moints2len ), intent(in) :: moints2
    integer :: p1, q1, p2, q2
    integer, dimension( aelec ) :: pstring1, pstring2
    integer, dimension( belec ) :: qstring1, qstring2

    integer :: diffs, adets, bdets
!--------------------------------------------------------------------
    adets = Binom( orbitals, aelec )
    bdets = Binom( orbitals, belec )
! Find determinant string indices for ind1 and ind2
    call K2Indc( ind1, belec, orbitals, p1, q1 )
    call K2Indc( ind1, belec, orbitals, p2, q2 )
! Find respective strings for p1, q1, p2, q2
    call GenOrbString( p1, aelec, orbitals, adets, pstring1 )
    call GenOrbString( p2, aelec, orbitals, adets, pstring2 )
    call GenOrbString( q1, belec, orbitals, bdets, qstring1 )
    call GenOrbString( q2, belec, orbitals, bdets, qstring2 )
! Test differences in strings. If > 2 orbitals, element is 0
    ham_element_diag = diagonal_mat( pstring1, qstring1, aelec, &
                           belec, moints1, moints1len, moints2, moints2len )
  end function ham_element_diag
!====================================================================
!====================================================================
end module
