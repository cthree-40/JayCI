! Driver for level2 of Det-CI
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chem., The Johns Hopkins University
!
! Last edit: 11-13-13
!====================================================================
! inputfl1    = main input file name                   character*20
! inputfl2    = input file generated from driver1      character*20
! outputfl    = output file                            character*20
! moflnm      = name of integral file                  character*20
! electrons   = total electrons                        integer scalar
! aelec       = alpha electrons                        integer scalar
! belec       = beta electrons                         integer scalar
! orbitals    = number of MO's                         integer scalar
! type1       = type of MO's to extract                integer scalar
! adets       = alpha determinants                     integer scalar
! bdets       = beta determinants                      integer scalar
! tdets       = truncated determinant list             integer array  1-d
! tadets      = truncated alpha string list            integer array  1-d
! tbdets      = truncated beta  string list            integer array  1-d
! tadetslen   = truncated alpha determinant strings    integer scalar
! tbdetslen   = truncated beta  determinant strings    integer scalar
! cidim       = truncated determinant expansion size   integer scalar
! nfrozen     = number of frozen orbitals              integer scalar
! ndocc       = number of docc orbitals                integer scalar
! nactive     = number of active orbitals              integer scalar
! xlevel      = excitation level                       integer scalar
! rtoler      = tolerance of CI convergence            real*8  scalar
! moints1     = 1-e integerals                         real*8  array  1-d
! moints2     = 2-e integerals                         real*8  array  1-d
! moints1len  = length of moints1                      integer scalar
! moints2len  = length of moints2                      integer scalar
! diagonals   = < K | H | K >                          real*8  array  1-d
! eigenvalues = output eigenvalues for roots desired   real*8  array  1-d
! roots       = number of roots desired                integer scalar
! krmin       = minimum dimension of krylov space      integer scalar
! krmax       = maximum dimension of krylov space      integer scalar
! maxiter     = maximum iteration of Davidson Alg.     integer scalar
!====================================================================

program driver2

  use detci1
  use detci5
  use orthogroutines
  use prediag

  implicit none

! ...file names...
  character*20 :: inputfl1, inputfl2, outputfl, astringfl, bstringfl, &
                  determfl, pstepfl, qstepfl, pstringfl, &
                  qstringfl, moflnm, plocfl, qlocfl, xreflsfl

! .....................
! ...NAMELIST INPUTS...

! ...runlevel1 integer scalars...
  integer :: electrons, orbitals, nfrozen, ndocc, nactive, xlevel

! ...expinfo integer scalars...
  integer :: cidim, aelec, belec, adets, bdets, &
             tadetslen, tbdetslen
  
! ...runlevel2 integer scalars...
  integer :: maxiter, krmin, krmax, roots, prediagr, initgdim

! ...runlevel2 real*8 scalars...
  real*8 :: rtoler, nucrep, ddot

! ......................
! ...DRIVER VARIABLES...

! ...driver integer parameters...
  integer, parameter :: diagsguess=30
! ...driver integer scalars...
  integer :: moints1len, moints2len, max1e, max2e, nrm

! ...driver integer arrays...
  integer, dimension(:),   allocatable :: tdets
  integer, dimension(:),   allocatable :: tadets, tbdets
  integer, dimension(:),   allocatable :: pstep, qstep, plocate, qlocate
  integer, dimension(:,:), allocatable :: pstring, qstring
  integer, dimension(:),   allocatable :: xreflist
  
! ...driver real*8 arrays...
  real*8, dimension(:),    allocatable :: diagonals
  real*8, dimension(:),    allocatable :: eigenvalues
  real*8, dimension(:,:),  allocatable :: eigenvectors
  real*8, dimension(:,:),  allocatable :: initvectors
  real*8, dimension(:),    allocatable :: moints1, moints2

! ...driver integer scalars...
  integer :: i, j, k, l
  integer :: openstat
  integer :: infl1, infl2, outfl, astrfl, bstrfl, detfl, steppfl, &
             stepqfl, locpfl, locqfl, strpfl, strqfl, xreffl
  integer :: subblockdim
  integer :: type1
!................
! ...NAMELISTS...
  namelist /runlevel1/ electrons, orbitals, nfrozen, ndocc, nactive, &
                       xlevel

  namelist /expinfo/ moflnm, cidim, aelec, belec, adets, bdets,&
                     tadetslen, tbdetslen

  namelist /runlevel2/ maxiter, krmin, krmax, rtoler, roots, nucrep, prediagr,&
                       initgdim
#ifdef DEBUG
  real*8, dimension(10,3) :: mat_a
#endif
!--------------------------------------------------------------------

! Assign file names and unit numbers
  inputfl1  = 'detci.in'
  infl1     = 1

  inputfl2  = 'input.2'
  infl2     = 3

  outputfl  = 'detci.out'
  outfl     = 2

  astringfl = 'alpha.dets'
  astrfl    = 8

  bstringfl = 'beta.dets'
  bstrfl    = 7

  determfl  = 'det.list'
  detfl     = 4
  
  pstepfl   = 'pstep.list'
  steppfl   = 13

  qstepfl   = 'qstep.list'
  stepqfl   = 11

  plocfl    = 'plocate.list'
  locpfl    = 14

  qlocfl    = 'qlocate.list'
  locqfl    = 12

  pstringfl = 'pstring.list'
  strpfl    = 9

  qstringfl = 'qstring.list'
  strqfl    = 10

  xreflsfl  = 'xref.list'
  xreffl    = 15
  
! Read input file 1
  openstat=0
  open(unit=infl1,file=inputfl1,status='old',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN INPUT FILE, detci.in. ****"
  read(unit=infl1,nml=runlevel1)
  read(unit=infl1,nml=runlevel2)
  close(unit=infl1)

! Read input file 2
  openstat=0
  open(unit=infl2,file=inputfl2,status='old',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN INPUT FILE, input.2 ****"
  read(unit=infl2,nml=expinfo)
  close(unit=infl2)

! Open outputfile, and begin writing output.
  openstat=0
  open(unit=outfl,file=outputfl,status='old',iostat=openstat,position='append')
  if ( openstat > 0 ) stop "**** CANNOT OPEN OUTPUT FILE, detci.out ****"

  write(unit=outfl,fmt=9)  "------------------------------------------------------------------------"
  write(unit=outfl,fmt=9)  "                            RUN LEVEL 2                                 "
  write(unit=outfl,fmt=9)  " "
  write(unit=outfl,fmt=9)  " Input from &runlevel1 "
  write(unit=outfl,fmt=10) "  electrons=", electrons
  write(unit=outfl,fmt=10) "  orbitals=", orbitals
  write(unit=outfl,fmt=10) "  nfrozen=", nfrozen
  write(unit=outfl,fmt=10) "  ndocc=", ndocc
  write(unit=outfl,fmt=10) "  nactive=", nactive
  write(unit=outfl,fmt=10) "  xlevel=", xlevel
  write(unit=outfl,fmt=9)  " "
  write(unit=outfl,fmt=9)  " Input from &expinfo "
  write(unit=outfl,fmt=11) "  moflnm=", moflnm
  write(unit=outfl,fmt=10) "  cidim=", cidim
  write(unit=outfl,fmt=10) "  aelec=", aelec
  write(unit=outfl,fmt=10) "  belec=", belec
  write(unit=outfl,fmt=10) "  adets=", adets
  write(unit=outfl,fmt=10) "  bdets=", bdets
  write(unit=outfl,fmt=10) "  tadetslen=", tadetslen
  write(unit=outfl,fmt=10) "  tbdetslen=", tbdetslen
  write(unit=outfl,fmt=9)  " "
  write(unit=outfl,fmt=9)  " Input from &runlevel2 "
  write(unit=outfl,fmt=10) "  maxiter=", maxiter
  write(unit=outfl,fmt=10) "  krmin=", krmin
  write(unit=outfl,fmt=10) "  krmax=", krmax
  write(unit=outfl,fmt=12) "  rtoler=", rtoler
  write(unit=outfl,fmt=10) "  roots=", roots
  write(unit=outfl,fmt=12) "  nucrep=", nucrep
  write(unit=outfl,fmt=10) "  initgdim=", initgdim
  write(unit=outfl,fmt=9)  " "

9  format(1x,A)
10 format(1x,A,I10)
11 format(1x,A,A)
12 format(1x,A,F10.7)
13 format(1x,I10)
14 format(1x,I10,I10)
15 format(1x,I3,A)

! Read in MO's
  write(unit=outfl,fmt=9) " Reading in MO's..."
  moints1len = ind2val(orbitals,orbitals)
  moints2len = index2e2(orbitals,orbitals,orbitals,orbitals)
  allocate(moints1(moints1len))
  allocate(moints2(moints2len))
  moints1 = 0d0
  moints2 = 0d0
  type1 = 1
  call iwfmt( moints1, moints2, type1, orbitals, moflnm, moints1len, moints2len )  
  write(unit=outfl,fmt=9) " "

! Allocate arrays for CI wavefunction truncation
  allocate(tadets(tadetslen))
  allocate(pstring(cidim,2))
  allocate(pstep(tadetslen))
  allocate(plocate(tadetslen))
  allocate(tbdets(tbdetslen))
  allocate(qstring(cidim,2))
  allocate(qstep(tbdetslen))
  allocate(qlocate(tbdetslen))
  allocate(tdets(cidim))

! Read in determinants
  write(unit=outfl,fmt=9) " Reading in strings and determinants..."

! Alpha strings
  openstat=0
  open(unit=astrfl,file=astringfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN ALPHA STRING FILE, alpha.dets. ****"
  read(unit=astrfl,fmt=13) ( tadets(i), i=1, tadetslen )
  close(unit=astrfl)

! Alpha string det list
  openstat=0
  open(unit=strpfl,file=pstringfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN PSTRING DETERMINANT LIST, pstring.list. ****"
  read(unit=strpfl,fmt=14) ( (pstring(i,j),j=1,2), i=1, cidim )
  close(unit=strpfl)

! Alpha string step list
  openstat=0
  open(unit=steppfl,file=pstepfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANOT OPEN PSTRING STEP FILE, pstep.list. ****"
  read(unit=steppfl,fmt=13) ( pstep(i), i=1, tadetslen )
  close(unit=steppfl)

! Alpha string locate list
  openstat=0
  open(unit=locpfl,file=plocfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN PLOCATE FILE, plocate.list. *****"
  read(unit=locpfl,fmt=13) ( plocate(i), i=1, tadetslen)
  close(unit=locpfl)

! Beta strings
  openstat=0
  open (unit=bstrfl,file=bstringfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN BETA STRING FILE, beta.dets. ****"
  read(unit=bstrfl,fmt=13) ( tbdets(i), i=1, tbdetslen )
  close(unit=bstrfl)

! Beta string det list
  openstat=0
  open(unit=strqfl,file=qstringfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN QSTRING DETERMINANT LIST, qstring.list. ****"
  read(unit=strqfl,fmt=14) ( (qstring(i,j),j=1,2), i=1, cidim )
  close(unit=strqfl)

! Beta string step list
  openstat=0
  open(unit=stepqfl,file=qstepfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN QSTRING STEP FILE, qstep.list. ****"
  read(unit=stepqfl,fmt=13) ( qstep(i), i=1, tbdetslen )
  close(unit=stepqfl)

! Beta sting locate list
  openstat=0
  open(unit=locqfl,file=qlocfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "*** CANNOT OPEN QLOCATE FILE, qlocate.list. ****"
  read(unit=locqfl,fmt=13) ( qlocate(i), i=1, tbdetslen)
  close(unit=locqfl)

! Determinants
  openstat=0
  open (unit=detfl,file=determfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN DETERMINANT FILE, dets.list. ****"
  read(unit=detfl,fmt=13) ( tdets(i), i=1, cidim )
  close(unit=detfl)

! Crossref list
  open (unit=xreffl,file=xreflsfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN CROSS REFERENCING FILE, xref.list. ****"
  read(unit=xreffl,fmt=13) ( xreflist(i), i=1, cidim )
  close(unit=xreffl)

! Compute diagonal matrix elements
  openstat=0
  write(unit=outfl,fmt=9) " Computing diagonal matrix elements..."
  allocate(diagonals(cidim))
  call diagonal( pstring, pstep, tadetslen, qstring, qstep, tbdetslen,       &
                 tadets, tadetslen, tbdets, tbdetslen, tdets, cidim, electrons,   &
                 orbitals, aelec, belec, adets, bdets, moints1, moints1len,&
                 moints2, moints2len, plocate, qlocate, diagonals )   

  write(unit=outfl,fmt=9) " "
  write(unit=outfl,fmt=9) " Diagonal matrix elements computed. "
  write(unit=outfl,fmt=9) " "
  write(unit=outfl,fmt=12)" -- Hartree Fock energy:    ", diagonals(1)-nucrep
  write(unit=outfl,fmt=9) " "

! Generate initial guess vectors
! This can be done three ways: 1) unit vectors 
!                              2) lowest diagonals
!                              3) prediagonalization of H subblock
  allocate( initvectors(cidim, initgdim ))
  if ( prediagr .eq. 1 ) then
    write(unit=outfl,fmt=9) " Using unit vectors for initial guess..."
    call geninituvec( initgdim, cidim, initvectors )
    write(unit=outfl,fmt=15) initgdim, " initial vectors generated."
    write(unit=outfl,fmt=9) " "
  else if ( prediagr .eq. 2 ) then
    write(unit=outfl,fmt=9) " Using lowest diagonals prediag. subroutine for initial guess..."
    call lowdiagprecond( diagonals, cidim, moints1, moints1len, moints2,      &
                         moints2len, pstring, pstep, plocate, qstring, qstep, &
                         qlocate, tadets, tbdets, tadetslen, tbdetslen,adets, &
                         bdets, aelec, belec, orbitals, xreflist, diagsguess, &
                         initgdim ,initvectors )
    write(unit=outfl,fmt=15) initgdim, " initial vectors generated."
    write(unit=outfl,fmt=9) " "
  else if ( prediagr .eq. 3 ) then 
    write(unit=outfl,fmt=9) " Using prediagonalization of H subblock for initial guess..."
    subblockdim = 200
    call prediagsubblock( cidim, moints1, moints2, moints1len, &
                          moints2len, tadets, tadetslen, tbdets, tbdetslen,   &
                          tdets, subblockdim, initgdim, aelec, belec, orbitals,initvectors )
    write(unit=outfl,fmt=15) initgdim, " initial vectors generated."
    write(unit=outfl,fmt=9) " "
  else if ( prediagr .eq. 4) then
    write(unit=outfl,fmt=9) " Using the subdiag1 subroutine... "
    call ovrlp_subspace( 200, initgdim, moints1, moints1len, moints2, moints2len, pstring,  &
                         pstep, plocate, qstring, qstep, qlocate, xreflist, cidim, tadets,  &
                         tbdets,tbdetslen, tadetslen, adets, bdets, aelec, belec, orbitals, &
                         diagonals, initvectors )
    write(unit=outfl,fmt=15) initgdim, " initial vectors generated."
    write(unit=outfl,fmt=9) " "
  end if

! Orthogonalize initial vectors
!  initvectors(1,1) = 1d0
  call modgramschmidt( initvectors, initgdim, initgdim, cidim )

! Call Davidson algorithm
  write(unit=outfl,fmt=9) " Calling Davidson algorithm..."
  allocate(eigenvalues(roots))
  allocate(eigenvectors(cidim,roots))
  call davidson( initgdim, initvectors, diagonals, moints1, moints1len, moints2, &
                 moints2len, cidim, pstring, pstep, plocate, qstring, qstep, qlocate, xreflist, &
                 tadets, tadetslen, tbdets, tbdetslen, adets, bdets, aelec, belec,   &
                 orbitals, krmin, krmax, rtoler, roots, maxiter, eigenvalues, eigenvectors )
  do i=1, roots
    write(unit=outfl,fmt=16) "Total CI energy for root # ", i," = ", (eigenvalues(i) - nucrep)
  end do
16 format(1x, A,I2,A,F10.7 )
  write(unit=outfl,fmt=9) "----------------------------------------"
  write(unit=outfl,fmt=9) " Finished. "
  write(unit=outfl,fmt=9) "-----------"
end program

 
