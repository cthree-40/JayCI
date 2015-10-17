program test
  ! program for testing jayci subroutines
  use test_module
  
  implicit none

  call test_citrunc()

  call test_cannon4()
  
  call test_strfind3()
  
end program test
