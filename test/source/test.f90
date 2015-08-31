program test
  ! program for testing jayci subroutines
  use test_module
  
  implicit none

  call test_citrunc()

  call test_expdiag()
  
end program test
