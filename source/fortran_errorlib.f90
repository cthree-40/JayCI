! File: fortran_errorlib
!
! Routines for error processing within fortran subroutines.
! Contains:
!  error_dropout
module fortran_errorlib

  implicit none

contains

  ! error_dropout: print error message and leave execution.
  subroutine error_dropout (estr, pname)
    implicit none
    character(255), intent(in) :: estr
    character(255), intent(in) :: pname

    write(*,"(A,A,A)") trim(adjustl(pname)),": ",trim(adjustl(estr))
    stop
    return
  end subroutine error_dropout
end module fortran_errorlib
