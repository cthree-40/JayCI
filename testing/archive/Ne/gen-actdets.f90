program genactdets
 
  implicit none

  integer :: i


  open(unit=1,file='actdets.in',status='new')
  do i=1, 4008005
    write(unit=1,fmt=10) i
  end do
  close(unit=1)
10 format(1x,I10)
end program
