! File: dlamch_fcn.f90
!*********************************************************************
! dlamch_fcn: returns output of dlamch call for call within C program
double precision function dlamch_fcn()
  implicit none
  character*1 :: safe_min = "s"
  double precision, external :: dlamch
  
  dlamch_fcn = dlamch(safe_min)

end function dlamch_fcn
  
