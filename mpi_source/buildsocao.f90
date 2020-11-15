! Subroutine to build SCAO reading soinfo.dat file.
subroutine buildsocao (natoms, atominfo)
  use fortran_errorlib
  implicit none

  ! ..input variables..
  integer, intent(in) :: natoms

  ! ..output variables..
  character(80), dimension(natoms), intent(out) :: atominfo
  ! ..local variables..
  integer, parameter :: maxbfn   = 1023
  integer, parameter :: maxatoms =  255

  ! Check atom number
  if (natoms .gt. maxatoms) then
          call error_dropout("natoms > maxatoms", "buildsocao")
  end if

  
  
  
end subroutine buildsocao
          
  
  
