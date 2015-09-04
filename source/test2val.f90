function test2val(val1, val2, thrsh) return(test)
  ! test2val
  ! --------
  ! Purpose: test 2 real*8 values and see if they are within a set threshold
  !-------------------------------------------------------------------------
  implicit none

  ! .. input arguments ..
  real*8, intent(in) :: val1, val2, thrsh

  ! .. output arguments ..
  integer :: test

  if (abs(val1 - val2) .lt. thrsh) then
          test = 0
  else
          test = 1
  end if

  return
end function test2val
