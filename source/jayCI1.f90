!********************************************************************
! DRIVER 1 for JayCI
!--------------------------------------------------------------------
PROGRAM jayCI1
  use input_proc, only: alpha_beta
  use combinatorial
  IMPLICIT NONE
  character(len=20) :: in1_file,in2_file,out1_file
  integer           :: electrons,openStat,nFrozen,nDOCC,nActive,xLevel,&
    Orbitals,taDLen,tbDLen,aElec,bElec,aDets,bDets,ciDim
  integer           :: unit_in1,unit_in2,unit_out1,unit_detlist, &
    unit_beta,unit_alpha,unit_pStr,unit_qStr,unit_qStep,unit_pStep,&
    unit_qLoc,unit_pLoc,unit_Xref
  NAMELIST /CIinfo/ electrons,Orbitals,nFrozen,nDOCC,nActive,xLevel
  !--------------------------------------------------------------
  !Assign unit numbers
  unit_in1     = 1        !jayCI.in
  unit_in2     = 2        !input.jayCI
  unit_out1    = 3        !jayCI.out
  unit_detlist = 4        !det.list   (truncated)
  unit_beta    = 7        !beta.dets  (truncated)
  unit_alpha   = 8        !alpha.dets (truncated)
  unit_pStr    = 9        !pstring.list ( p,q pairs searching through p )
  unit_qStr    = 10       !qstring.list ( q,p pairs searching through q )
  unit_qStep   = 11       !qstep.list   ( # of q_j corresponding to p_i )
  unit_pStep   = 12       !pstep.list   ( # of p_j corresponding to q_i )
  unit_qLoc    = 13       !qlocate.list ( Location of q_i )
  unit_pLoc    = 14       !plocate.list ( Location of p_i )
  unit_Xref    = 15       !Map of q,p to p,q
  ! Read input file jayCI.in
  in1_file='jayCI.in'
  OPEN(unit=unit_in1,file=in1_file,status='old',IOSTAT=openStat)
  if(openStat.GT.0) STOP "*** CANNOT OPEN INPUT FILE: jayCI.in ***"
  READ(unit=unit_in1,nml=CIinfo)
  CLOSE(unit_in1)
  ! Open output file and write header information
  CALL write_header1(unit_out1)
  ! Write input information
  out1_file="jayCI.out"
  OPEN(unit=unit_out1,file=out1_file,status='old',IOSTAT=openStat)
  if(openStat.GT.0) STOP "*** CANNOT OPEN OUTPUT FILE: jayCI.out ***"
  WRITE(unit=unit_out1,fmt=10) " ---------------------------- 1 -----------------------------------"
  WRITE(unit=unit_out1,fmt=10) " Namelist input: "
  WRITE(unit=unit_out1,fmt=11) "  Electrons in system:   ", electrons
  WRITE(unit=unit_out1,fmt=11) "  Orbitals in system:    ", Orbitals
  WRITE(unit=unit_out1,fmt=11) "  Frozen orbitals:       ", nFrozen
  WRITE(unit=unit_out1,fmt=11) "  DOCC orbitals:         ", nDOCC
  WRITE(unit=unit_out1,fmt=11) "  Active orbitals:       ", nActive
  WRITE(unit=unit_out1,fmt=11) "  Excitation level:      ", xLevel
  WRITE(unit=unit_out1,fmt=10) " "
10 format(1x,A)
11 format(1x,A,I10)
  ! Call alpha_beta to find correct number of alpha and beta electrons
  CALL alpha_beta(Electrons,aElec,bElec)
  WRITE(unit=unit_out1,fmt=11) " Alpha Electrons: ", aElec
  WRITE(unit=unit_out1,fmt=11) " Beta  Electrons: ", bElec
  WRITE(unit=unit_out1,fmt=10) " "
  ! Compute aDets, and bDets
  aDets = Binom(Orbitals,aElec)
  bDets = Binom(Orbitals,bElec)
  WRITE(unit=unit_out1,fmt=11) " Alpha Determinants: ", aDets
  WRITE(unit=unit_out1,fmt=11) " Beta  Determinants: ", bDets
  WRITE(unit=unit_out1,fmt=10) " "
  WRITE(unit=unit_out1,fmt=11) " Size of Full-CI expansion: ", aDets*bDets
  WRITE(unit=unit_out1,fmt=10) " "
  ! Construct input for our truncated space
  CALL citrunc( aDets,bDets,aElec,bElec,Orbitals,nFrozen,nDOCC,nActive,  &
    xLevel,taDLen,tbDlen,ciDim,unit_detlist,unit_beta,unit_alpha,    &
    unit_pStr,unit_qStr,unit_qStep,unit_qLoc,unit_pStep,unit_pLoc,   &
    unit_Xref )
  ! Write size of final expansion
  WRITE(unit=unit_out1,fmt=11) " Size of Trucated expansion: ", ciDim
  ! Generate input for jayCI2.x
  in2_file="input.jayCI"
  OPEN(unit=unit_in2,file=in2_file,status='new',IOSTAT=openStat)
  if(openStat.GT.0) STOP "*** CANNOT OPEN OUTPUT FILE: input.jayCI ***"
  WRITE(unit=unit_in2,fmt=10)"&ExpInfo"
  WRITE(unit=unit_in2,fmt=10)"MOfile='moints'"
  WRITE(unit=unit_in2,fmt=11)"ciDim=",ciDim
  WRITE(unit=unit_in2,fmt=11)"aElec=",aElec
  WRITE(unit=unit_in2,fmt=11)"bElec=",bElec
  WRITE(unit=unit_in2,fmt=11)"aDets=",aDets
  WRITE(unit=unit_in2,fmt=11)"bDets=",bDets
  WRITE(unit=unit_in2,fmt=11)"pDLen=",taDLen
  WRITE(unit=unit_in2,fmt=11)"qDLen=",tbDLen
  WRITE(unit=unit_in2,fmt=10)"/"
  CLOSE(unit=unit_in2)
  ! Completed jayCI1...
  WRITE(unit=unit_out1,fmt=10) " jayCI1 finished..."
  WRITE(unit=unit_out1,fmt=10) " ---------------------------------------------------------------"
  CLOSE(unit=unit_out1)
!--------------------------------------------------------------------
CONTAINS
  !====================================================================
  !>write_header1
  ! Write header, creating the file "jayCI.out"
  !--------------------------------------------------------------------
  subroutine write_header1(file_unit)
    IMPLICIT NONE
    integer,intent(IN)      :: file_unit
    character(len=20)       :: file_name
    integer                 :: OpenStatus
    
    !Open jayCI.out and write header
    file_name="jayCI.out"
    OPEN(unit=file_unit,file=file_name,status='new',iostat=OpenStatus)
    if(OpenStatus.GT.0)STOP "*** CANNOT OPEN OUTPUT FILE: jayCI.out ***"
    !Write header
    write(unit=file_unit,fmt=9) "========================================================================"
    write(unit=file_unit,fmt=9) "--------------------------DETERMINANT CI--------------------------------"
    write(unit=file_unit,fmt=9) "========================================================================"
    write(unit=file_unit,fmt=9) " "
    write(unit=file_unit,fmt=9) " By Christopher L. Malbon"
    write(unit=file_unit,fmt=9) " Yarkony Group"
    write(unit=file_unit,fmt=9) " Dept. of Chemistry, The Johns Hopkins University"
    write(unit=file_unit,fmt=9) " "
    write(unit=file_unit,fmt=9) " "
    write(unit=file_unit,fmt=9) "            ---------- Created on SLACKWARE LINUX ----------              "
    write(unit=file_unit,fmt=9) "                           www.slackware.com                              "
    write(unit=file_unit,fmt=9) "------------------------------------------------------------------------"
    write(unit=file_unit,fmt=9) "Send bug reports to cmalbon1@jhu.edu ------------"
    write(unit=file_unit,fmt=9) "------------------------------------------------------------------------"
    write(unit=file_unit,fmt=9) " "
9   format(1x,A)
    !Close file
    CLOSE(unit=file_unit)
    RETURN
  end subroutine write_header1
  !====================================================================
END PROGRAM jayCI1

