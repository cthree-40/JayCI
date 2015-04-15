!********************************************************************
! DRIVER 1 for JayCI
!--------------------------------------------------------------------
PROGRAM jayCI1
      use input_proc, only: alpha_beta
      use combinatorial, only: Binom
      IMPLICIT NONE
      character(len=20) :: inputfl, outputfl
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
      READ(unit=unit_in1,nml=WFinfo)
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
10    format(1x,A)
11    format(1x,A,I10)
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
9     format(1x,A)
            !Close file
            CLOSE(unit=file_unit)
            RETURN
      end subroutine write_header1
!====================================================================
END PROGRAM
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
            pLocate,qLocate,xRefList
      integer,dimension(:,:),ALLOCATABLE  :: pString,qString
      real*8,dimension(:),ALLOCATABLE     :: Diagonals,Eig_Values,MOints1,MOints2
      real*8,dimension(:,:),ALLOCATABLE   :: Eig_Vectors,Init_Vectors
      integer     :: openStat
      integer     :: unit_in1,unit_in2,unit_out1,unit_alpha,unit_beta,unit_dets,&
            unit_pStep,unit_qStep,unit_pLoc,unit_qLoc,unit_pStr,unit_qStr,      &
            unit_moint,unit_xRef
      NAMELIST /CIinfo/ electrons,Orbitals,nFrozen,nDocc,nActive,xLevel
      NAMELIST /ExpInfo/ MOfile,ciDim,aElec,bElec,aDets,bDets,pDLen,qDLen
      NAMELIST /Davidson/ maxIter,KryMin,KryMax,res_Tol,Roots,Nuc_Rep_E,   &
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
      READ(unit=unit_in1,NML=Davidson)
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
9     format(1x,A)
10    format(1x,A,I10)
11    format(1x,A,A)
12    format(1x,A,F12.8)
13    format(1x,I10)
14    format(1x,I10,I10)
15    format(1x,I3,A)
      ! Read in MO's
      WRITE(unit=unit_out1,fmt=9) " Reading in MO's... "
      M1Len=Ind2Val(Orbitals,Orbitals)
      M2Len=index2e2(Orbitals,Orbitals,Orbitals,Orbitals)
      ! Allocate integral arrays
      ALLOCATE(MOints1(M1Len))
      ALLOCATE(MOints2(M2Len))
      MOints1=0d0
      MOints2=0d0
      !type1=1...<-- this chooses what integrals to read in. 
      ! **Modification of type1 will come when dipole integrals are needed, etc
      type1=1
      CALL iwfmt(MOints1,MOints2,type1,Orbitals,MOfile,M1Len,M2Len)
      !Allocate wavefunction information arrays
      ALLOCATE(pDets(pDLen))
      ALLOCATE(qDets(qDLen))
      ALLOCATE(pString(ciDim,2))
      ALLOCATE(qString(ciDim,2))
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
      READ(unit=unit_pStr,fmt=14)((pString(i,j),j=1,2),i=1,ciDim)
      CLOSE(unit=unit_pStr)
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
      !Determinant list, searching through beta strings
      OPEN(unit=unit_qStr,file=file_qStr,status='old',position='rewind', &
            IOSTAT=openStat)
      if(openStat.GT.0) STOP "*** CANNOT OPEN QSTRING(:,2) FILE, qstring.list *** "
      READ(unit=unit_qStr,fmt=14) ((qString(i,j),j=1,2),i=1,ciDim)
      CLOSE(unit=unit_qStr)
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
            subBlockDim=200
            CALL preDiagSubBlock(ciDim,MOints1,M1Len,MOints2,M2Len,&
                  Det_Exp,subBlockDim,initGuessDim,aElec,bElec,    &
                  Orbitals,Init_Vectors)
            WRITE(unit=unit_out1,fmt=15) initGuessDim, " initial vectors generated..."
            WRITE(unit=unit_out1,fmt=15) " "
      else if(preDiagRoutine.GT.1)then !NOT IMPLEMENTED CURRENTLY
            STOP " Chosen prediagonalization routine is not currently implemented. "
      end if
      !Orthogonalize initial vectors
      CALL modGramSchmidt(Init_Vectors,initGuessDim,initGuessDim,ciDim)
      !Call Davidson algorithm
      CALL davidson(



