subroutine cannonreal( string, length, permind )
!==========================================================
! This subroutine returns a string to cannonical ordering
!  and gives the character of the permutation
!==========================================================
!----------------------------
! Input:
!  length = length of string
!  string = string of dimension length
! Output:
!  string = reordered string
!  permind = index of the permutations
!----------------------------
  implicit none
  integer, intent(in) :: length
  real*8, dimension(length) :: string
  integer, intent(out) :: permind

  integer :: k, i, j
  real*8 :: temp
!----------------------------
! 

! permind set to 1
  permind = 1

! Loop through string(:)
  do i=1, (length - 1)

! Test if string(i) > string(i+1)
    if ( string(i) > string(i+1) ) then

! Loop through previous elements to see where string(i+1)
!  belongs
      do j=1, i

! If the i+1 element is less than the jth element, it 
!  must occupy the jth element's position. Bump all other
!  elements up.
        if ( string(i+1) < string(j) ) then
          temp = string(i+1)
          do k=i, j, -1
            string(k+1) = string(k)
! Multiply the permutational index by -1
            permind = permind + 1
          end do
! Now place i+1 in the jth position
          string(j) = temp 
        else

! Try the next element j (j -> j+1)
          cycle
        end if
      end do
    else
! Try the next element i (i -> i+1)
      cycle
    end if
  end do

  return        

end subroutine
