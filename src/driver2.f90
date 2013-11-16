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
  implicit none

! ...file names...
  character*20 :: inputfl1, inputfl2, outputfl, astringfl, bstringfl, &
                  determfl

! .....................
! ...NAMELIST INPUTS...

! ...runlevel1 integer scalars...
  integer :: electrons, orbitals, nfrozen, ndocc, nactive, xlevel

! ...expinfo integer scalars...
  integer :: cidim, aelec, belec, adets, bdets, &
             tadetslen, tbdetslen

! ...expinfo character*20...
  character*20 :: moflnm
  
! ...runlevel2 integer scalars...
  integer :: maxiter, krmin, krmax, roots

! ...runlevel2 real*8 scalars...
  real*8 :: rtoler

! ......................
! ...DRIVER VARIABLES...

! ...driver integer scalars...
  integer :: moints1len, moints2len, max1e, max2e

! ...driver integer arrays...
  integer, dimension(:), allocatable :: tdets
  integer, dimension(:), allocatable :: tadets, tbdets
  
! ...driver real*8 arrays...
  real*8, dimension(:), allocatable :: diagonals
  real*8, dimension(:), allocatable :: eigenvalues


! ...driver integer scalars...
  integer :: i, j, k, l
  integer :: openstat

!................
! ...NAMELISTS...
  namelist /runlevel1/ electrons, orbitals, nfrozen, ndocc, nactive, &
                       xlevel

  namelist /expinfo/ moflnm, cidim, aelec, belec, adets, bdets,&
                     tadetslen, bdetslen

  namelist /runlevel2/ maxiter, krmin, krmax, rtoler, roots

!--------------------------------------------------------------------

! Assign file names
  inputfl1 = 'detci.in'
  inputfl2 = 'input.2'
  outputfl = 'detci.out'
  astringfl= 'alpha.dets'
  bstringfl= 'beta.dets'
  determfl = 'dets.list'
  
! Read input file 1
  open(unit=1,file=inputfl1,status='old',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN INPUT FILE, detci.in. ****"
  read(unit=1,nml=runlevel1)
  close(unit=1)

! Read input file 2
  open(unit=2,file=inputfl2,status='old',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN INPUT FILE, input.2 ****"
  read(unit=2,nml=expinfo)
  close(unit=2)

! Open outputfile, and begin writing output.
  open(unit=3,file=outputfl,status='old',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN OUTPUT FILE, detci.out ****"

  write(unit=2,fmt=9)  "------------------------------------------------------------------------"
  write(unit=2,fmt=9)  "                            RUN LEVEL 2                                 "
  write(unit=2,fmt=9)  " "
  write(unit=2,fmt=9)  " Input from &runlevel1 "
  write(unit=2,fmt=10) "  electrons=", electrons
  write(unit=2,fmt=10) "  orbitals=", orbitals
  write(unit=2,fmt=10) "  nfrozen=", nfrozen
  write(unit=2,fmt=10) "  ndocc=", ndocc
  write(unit=2,fmt=10) "  nactive=", nactive
  write(unit=2,fmt=10) "  xlevel=", xlevel
  write(unit=2,fmt=9)  " "
  write(unit=2,fmt=9)  " Input from &expinfo "
  write(unit=2,fmt=11) "  moflnm=", moflnm
  write(unit=2,fmt=10) "  cidim=", cidim
  write(unit=2,fmt=10) "  aelec=", aelec
  write(unit=2,fmt=10) "  belec=", belec
  write(unit=2,fmt=10) "  adets=", adets
  write(unit=2,fmt=10) "  bdets=", bdets
  write(unit=2,fmt=10) "  tadetslen=", tadetslen
  write(unit=2,fmt=10) "  tbdetslen=", tbdetslen
  write(unit=2,fmt=9)  " "
  write(unit=2,fmt=9)  " Input from &runlevel2 "
  write(unit=2,fmt=10) "  maxiter=", maxiter
  write(unit=2,fmt=10) "  krmin=", krmin
  write(unit=2,fmt=10) "  krmax=", krmax
  write(unit=2,fmt=12) "  rtoler=", rtoler
  write(unit=2,fmt=10) "  roots=", roots
  write(unit=2,fmt=9)  " "

9  format(1x,A)
10 format(1x,A,I10)
11 format(1x,A,A)
12 format(1x,A,ES10.1)
13 format(1x,I10)

! Read in MO's
  write(unit=2,fmt=9) " Reading in MO's..."
  moints1len = ind2val(orbitals,orbitals)
  moints2len = index2e2(orbitals,orbitals,orbitals,orbitals)
  allocate(moints1(moints1len))
  allocate(moints2(moints2len))
  moints1 = 0d0
  moints2 = 0d0
  call iwfmt( moints1, moints2, type1, orbtls, moflnm, moints1len, moints2len )  
  write(unit=2,fmt=9) " "

! Read in determinants
  write(unit=2,fmt=9) " Reading in strings and determinants..."
! Alpha strings
  open(unit=4,file=astringfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN ALPHA STRING FILE, alpha.dets. ****"
  read(unit=4,fmt=13) ( tadets(i), i=1, tadetslen )
  close(unit=4)
! Beta strings
  open (unit=5,file=bstringfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN BETA STRING FILE, beta.dets. ****"
  read(unit=5,fmt=13) ( tbdets(i), i=1, tbdetslen )
  close(unit=5)
! Determinants
  open (unit=6,file=determfl,status='old',position='rewind',iostat=openstat)
  if ( openstat > 0 ) stop "**** CANNOT OPEN DETERMINANT FILE, dets.list. ****"
  read(unit=6,fmt=13) ( tdets(i), i=1, cidim )
  close(unit=6)
  write(unit=2,fmt=9) " "


! Compute diagonal matrix elements
  write(unit=2,fmt=9) " Computing diagonal matrix elements..."



















 
