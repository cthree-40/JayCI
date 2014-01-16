program test
  
  use orthogroutines

  implicit none

  integer, parameter :: n=10, m=10

  real*8, dimension( m, n ) :: matrix
  
  call modgramschmidt( matrix, m, m, n )
end program

