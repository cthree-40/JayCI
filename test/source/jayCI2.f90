!********************************************************************
! DRIVER 2 for JayCI
!--------------------------------------------------------------------
PROGRAM jayCI2
  use integral
  use action_util,        only: diag_element1
  use orthogroutines,     only: modGramSchmidt
  use InitialGuess,       only: preDiagSubBlock
  IMPLICIT NONE
  character(len=20) :: file_in1,file_in2,file_out1,file_alpha,file_beta, &
    file_dets,file_pStep,file_qStep,file_pStr,file_qStr,MOfile,      &
    file_pLoc,file_qLoc,file_xRef
  integer     :: electrons,Orbitals,nFrozen,nDocc,nActive,xLevel
  integer     :: ciDim,aElec,bElec,aDets,bDets,pDLen,qDlen
  integer     :: maxIter,KryMin,KryMax,Roots,preDiagRoutine,initGuessDim
  real*8      :: res_Tol,Nuc_Rep_E
  real*8      :: ddot
  integer,parameter :: Diagional_Guess=30
  integer     :: M1Len,M2Len
  integer,dimension(:),ALLOCATABLE    :: Det_Exp,pDets,qDets,pStep,qStep,&
    pLocate,qLocate,xRefList,pString,qString
  integer,dimension(:,:),ALLOCATABLE  :: pString_temp,qString_temp
  real*8,dimension(:),ALLOCATABLE     :: Diagonals,Eig_Values,MOints1,MOints2
  real*8,dimension(:,:),ALLOCATABLE   :: Eig_Vectors,Init_Vectors
  integer     :: openStat, i,j
  integer     :: type1, subBlockDim   ! These integer scalars are declared in the program
  integer     :: unit_in1,unit_in2,unit_out1,unit_alpha,unit_beta,unit_dets,&
    unit_pStep,unit_qStep,unit_pLoc,unit_qLoc,unit_pStr,unit_qStr,      &
    unit_moint,unit_xRef
  logical     :: moints_exists
  NAMELIST /CIinfo/ electrons,Orbitals,nFrozen,nDocc,nActive,xLevel
  NAMELIST /ExpInfo/ MOfile,ciDim,aElec,bElec,aDets,bDets,pDLen,qDLen
  NAMELIST /DavidInfo/ maxIter,KryMin,KryMax,res_Tol,Roots,Nuc_Rep_E,   &
    preDiagRoutine,initGuessDim
  ! Assing prediagonalization routine
  preDiagRoutine=1 !1=prediagonalization of a subblock of <k|H|k>
  ! Assign file names and unit numbers
  file_in1   = 'jayCI.in'
  unit_in1   = 1
  file_in2   = 'input.jayCI'
  unit_in2   = 2
  file_out1  = 'jayCI.out'
  unit_out1  = 3
  file_dets  = 'det.list'
  unit_dets  = 4
  file_alpha = 'alpha.dets'
  unit_alpha = 7
  file_beta  = 'beta.dets'
  unit_beta  = 8
  file_pStep = 'pstep.list'
  unit_pStep = 9
  file_qStep = 'qstep.list'
  unit_qStep = 10
  file_pStr  = 'pstring.list'
  unit_pStr  = 11
  file_qStr  = 'qstring.list'
  unit_qStr  = 12
  file_pLoc  = 'plocate.list'
  unit_pLoc  = 13
  file_qLoc  = 'qlocate.list'
  unit_qLoc  = 14
  file_xRef  = 'xref.list'
  unit_xRef  = 15
  !Read input from file 1
  OPEN(unit=unit_in1,file=file_in1,status='old',IOSTAT=openStat)
  if(openStat.GT.0)STOP "*** CANNOT OPEN INPUT FILE: jayCI.in *** "
  READ(unit=unit_in1,NML=CIinfo)
  READ(unit=unit_in1,NML=DavidInfo)
  CLOSE(unit=unit_in1)
  !Read input from file 2
  OPEN(unit=unit_in2,file=file_in2,status='old',IOSTAT=openStat)
  if(openStat.GT.0)STOP "*** CANNOT OPEN INPUT FILE: input.jayCI *** "
  READ(unit=unit_in2,NML=ExpInfo)
  CLOSE(unit=unit_in2)
  !Open output file and begin writing output
  OPEN(unit=unit_out1,file=file_out1,status='old',IOSTAT=openStat, &
    POSITION='append')
  if(openStat.GT.0)STOP "*** CANNOT OPEN OUTPUT FILE: jayCI.out *** "
  write(unit=unit_out1,fmt=9)  " Input from &CIinfo "
  write(unit=unit_out1,fmt=10) "  electrons=", electrons
  write(unit=unit_out1,fmt=10) "  orbitals=", Orbitals
  write(unit=unit_out1,fmt=10) "  nFrozen=", nFrozen
  write(unit=unit_out1,fmt=10) "  nDocc=", nDocc
  write(unit=unit_out1,fmt=10) "  nActive=", nActive
  write(unit=unit_out1,fmt=10) "  xLevel=", xLevel
  write(unit=unit_out1,fmt=9)  " "
  write(unit=unit_out1,fmt=9)  " Input from &ExpInfo "
  write(unit=unit_out1,fmt=11) "  MOfile=", MOfile
  write(unit=unit_out1,fmt=10) "  ciDim=", ciDim
  write(unit=unit_out1,fmt=10) "  aElec=", aElec
  write(unit=unit_out1,fmt=10) "  bElec=", bElec
  write(unit=unit_out1,fmt=10) "  aDets=", aDets
  write(unit=unit_out1,fmt=10) "  bDets=", bDets
  write(unit=unit_out1,fmt=10) "  pDLen=", pDLen
  write(unit=unit_out1,fmt=10) "  qDLen=", qDLen
  write(unit=unit_out1,fmt=9)  " "
  write(unit=unit_out1,fmt=9)  " Input from &Davidson "
  write(unit=unit_out1,fmt=10) "  maxiter=", maxIter
  write(unit=unit_out1,fmt=10) "  krmin=", KryMin
  write(unit=unit_out1,fmt=10) "  krmax=", KryMax
  write(unit=unit_out1,fmt=12) "  res_Tol=", res_Tol
  write(unit=unit_out1,fmt=10) "  Roots=", Roots
  write(unit=unit_out1,fmt=12) "  Nuc_Rep_E=", Nuc_Rep_E
  write(unit=unit_out1,fmt=10) "  initGuessDim=", initGuessDim
  write(unit=unit_out1,fmt=9)  " "
9 format(1x,A)
10 format(1x,A,I10)
11 format(1x,A,A)
12 format(1x,A,F12.8)
13 format(1x,I10)
14 format(1x,I10,I10)
15 format(1x,I3,A)
  ! Read in MO's
  WRITE(unit=unit_out1,fmt=9) " Reading in MO's... "
  M1Len=Ind2Val(Orbitals,Orbitals)
  M2Len=Index2E(Orbitals,Orbitals,Orbitals,Orbitals)
  ! Allocate integral arrays
  ALLOCATE(MOints1(M1Len))
  ALLOCATE(MOints2(M2Len))
  MOints1=0d0
  MOints2=0d0
  !type1=1...<-- this chooses what integrals to read in. 
  ! **Modification of type1 will come when dipole integrals are needed, etc
  type1=1
  INQUIRE(file=MOfile,EXIST=moints_exists)!Check if moints is in directory
  IF(.NOT.moints_exists)STOP " No molecular integral file. "
  CALL iwfmt(MOints1,MOints2,type1,Orbitals,MOfile,M1Len,M2Len)
  !Allocate wavefunction information arrays
  ALLOCATE(pDets(pDLen))
  ALLOCATE(qDets(qDLen))
  ALLOCATE(pString_temp(ciDim,2))
  ALLOCATE(qString_temp(ciDim,2))
  ALLOCATE(pLocate(pDLen))
  ALLOCATE(qLocate(qDLen))
  ALLOCATE(pStep(pDLen))
  ALLOCATE(qStep(qDLen))
  ALLOCATE(Det_Exp(ciDim))
  ALLOCATE(xRefList(ciDim))
  !Read in wavefunction information
  WRITE(unit=unit_out1,fmt=9) " Reading in wavefunction information..."
  !Alpha strings
  OPEN(unit=unit_alpha,file=file_alpha,status='old',position='rewind', &
    IOSTAT=openStat)
  if(openStat.GT.0) STOP "*** CANNOT OPEN ALPHA STRING FILE, alpha.dets *** "
  READ(unit=unit_alpha,fmt=13) (pDets(i),i=1,pDLen)
  CLOSE(unit=unit_alpha)
  !Determinant list, searching through alpha strings
  OPEN(unit=unit_pStr,file=file_pStr,status='old',position='rewind', &
    IOSTAT=openStat)
  if(openStat.GT.0) STOP "*** CANNOT OPEN PSTRING(:,2) FILE, pstring.list *** "
  READ(unit=unit_pStr,fmt=14)((pString_temp(i,j),j=1,2),i=1,ciDim)
  CLOSE(unit=unit_pStr)
  ! Transfer second column of pString_temp to pString
  ALLOCATE(pString(ciDim))
  pString(1:ciDim)=pString_temp(1:ciDim,2)
  DEALLOCATE(pString_temp)      !Deallocate pString_temp
  !alpha string step list
  OPEN(unit=unit_pStep,file=file_pStep,status='old',position='rewind', &
    IOSTAT=openStat)
  if(openStat.GT.0) STOP "*** CANNOT OPEN PSTEP(:) FILE, pstep.list *** "
  READ(unit=unit_pStep,fmt=13) (pStep(i),i=1,pDLen)
  CLOSE(unit=unit_pStep)
  !alpha string locate list
  OPEN(unit=unit_pLoc,file=file_pLoc,status='old',position='rewind', &
    IOSTAT=openStat)
  if(openStat.GT.0) STOP "*** CANNOT OPEN PLOCATE(:) FILE, plocate.list *** "
  READ(unit=unit_pLoc,fmt=13) (pLocate(i),i=1,pDLen)
  CLOSE(unit=unit_pLoc)
  !beta string file
  OPEN(unit=unit_beta,file=file_beta,status='old',position='rewind',&
    IOSTAT=openStat)
  if(openStat.GT.0) STOP "*** CANNOT OPEN BETA STRING FILE, beta.dets *** "
  READ(unit=unit_beta,fmt=13) (qDets(i),i=1,qDLen)
  CLOSE(unit=unit_beta)
  !Determinant list, searching through beta strings
  OPEN(unit=unit_qStr,file=file_qStr,status='old',position='rewind', &
    IOSTAT=openStat)
  if(openStat.GT.0) STOP "*** CANNOT OPEN QSTRING(:,2) FILE, qstring.list *** "
  READ(unit=unit_qStr,fmt=14) ((qString_temp(i,j),j=1,2),i=1,ciDim)
  CLOSE(unit=unit_qStr)
  ! Transfer second column of qString_temp to qString
  ALLOCATE(qString(ciDim))
  qString(1:ciDim)=qString_temp(1:ciDim,2)
  DEALLOCATE(qString_temp)      !Deallocate qString_temp
  !beta string step list
  OPEN(unit=unit_qStep,file=file_qStep,status='old',position='rewind', &
    IOSTAT=openStat )
  if(openStat.GT.0) STOP "*** CANNOT OPEN QSTEP(:) FILE, qstep.list *** "
  READ(unit=unit_qStep,fmt=13) (qStep(i),i=1,qDLen)
  CLOSE(unit=unit_qStep)
  !beta string locate list
  OPEN(unit=unit_qLoc,file=file_qLoc,status='old',position='rewind', &
    IOSTAT=openStat)
  if(openStat.GT.0) STOP "*** CANNOT OPEN QLOCATE(:) FILE, qlocate.list *** "
  READ(unit=unit_qLoc,fmt=13) (qLocate(i),i=1,qDLen)
  CLOSE(unit=unit_qLoc)
  !Determinant list
  OPEN(unit=unit_dets,file=file_dets,status='old',position='rewind', &
    IOSTAT=openStat)
  if(openStat.GT.0) STOP "*** CANNOT OPEN DETERMINANT FILE, dets.list *** "
  READ(unit=unit_dets,fmt=13) ( Det_Exp(i),i=1,ciDim )
  CLOSE(unit=unit_dets)
  !Cross reference list
  OPEN(unit=unit_xRef,file=file_xRef,status='old',position='rewind', &
    IOSTAT=openStat)
  if(openStat.GT.0) STOP "*** CANNOT OPEN CROSS REFERENCE FILE, xref.list *** "
  READ(unit=unit_xRef,fmt=13) ( xRefList(i),i=1,ciDim )
  CLOSE(unit=unit_xRef)
  ! Compute diagonal matrix elements
  WRITE(unit=unit_out1,fmt=9) " Computing diagonal matrix elements... "
  ALLOCATE(Diagonals(ciDim)) ! Allocate diagonal element array
  do i=1,ciDim
          Diagonals(i)=diag_element1(Det_Exp(i),MOints1,M1Len,MOints2,M2Len,&
            aElec,bElec,Orbitals)
  end do
  WRITE(unit=unit_out1,fmt=9) " "
  WRITE(unit=unit_out1,fmt=12) " Hartree Fock Energy:   ", (Diagonals(1)+Nuc_Rep_E)
  WRITE(unit=unit_out1,fmt=9) " "
  !Generate initial guess vectors
  !Allocate arrays
  ALLOCATE(Init_Vectors(ciDim,initGuessDim))
  if(preDiagRoutine.EQ.1)then !Subblock diagonalization

! PP Block for explicit hamiltonian diagonalization
#ifdef HAMILTONIAN
          subBlockDim=ciDim
          print *, "Calling diagaonalization subroutine"
          call preDiagSubBlock(ciDim,MOints1,M1Len,MOints2,M2Len, &
            Det_Exp, subBlockDim, initGuessDim, aElec, bElec,     &
            Orbitals, Init_Vectors)
          STOP "Done."
#endif
          subBlockDim=200
          print *, "Calling prediagonalization subroutine "
          CALL preDiagSubBlock(ciDim,MOints1,M1Len,MOints2,M2Len,&
            Det_Exp,subBlockDim,initGuessDim,aElec,bElec,    &
            Orbitals,Init_Vectors)
          WRITE(unit=unit_out1,fmt=15) initGuessDim, " initial vectors generated..."
          WRITE(unit=unit_out1,fmt=9) " "
  else if(preDiagRoutine.GT.1)then !NOT IMPLEMENTED CURRENTLY
          STOP " Chosen prediagonalization routine is not currently implemented. "
  end if
  !Orthogonalize initial vectors
  CALL modGramSchmidt(Init_Vectors,initGuessDim,initGuessDim,ciDim)
  !Call Davidson algorithm
  ALLOCATE(Eig_Values(Roots))
  ALLOCATE(Eig_Vectors(ciDim,Roots))
  print *, "initGuessDim=",initGuessDim
  CALL Davidson(initGuessDim,Init_Vectors,Diagonals,MOints1,M1Len,MOints2,M2Len,ciDim,&
    pString,pLocate,pStep,qString,qLocate,qStep,xRefList,pDets,pDLen,qDets,qDLen, &
    aDets,bDets,aElec,bElec,Orbitals,KryMin,KryMax,res_Tol,Roots,maxIter,nFrozen, &
    nDocc,nActive,Eig_Values,Eig_Vectors)
  DEALLOCATE(Eig_Values,Eig_Vectors)
END PROGRAM jayCI2



