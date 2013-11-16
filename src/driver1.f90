! Driver for level1 of Det-CI
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chem., The Johns Hopkins University
!
! Last edit: 11-13-13
!====================================================================
! inputfl   = input file name                    character*20
! outputfl  = output file name                   character*20
! electrons = total electrons                    integer scalar
! orbitals  = number of MO's                     integer scalar
! openstat  = status of open file                integer scalar
! nfrozen   = number of frozen orbitals          integer scalar
! ndocc     = number of docc orbitals            integer scalar
! nactive   = number of active orbitals          integer scalar
! xlevel    = excitation level                   integer scalar
!==================================================================== 
program driver1 
  
  use detci1
  use detci2
  use citrunc

  implicit none


! ...file names...
  character*20 :: inputfl, outputfl

! ...NAMELIST integer scalars...
  integer :: electrons, openstat, nfrozen, ndocc, nactive, xlevel, &
             orbitals, tadetslen, tbdetslen
 
! ...NAMELIST input...
  namelist /runlevel1/ electrons, orbitals, nfrozen, ndocc, nactive, &
                       xlevel

! ...integer scalars...
  integer :: openstatus
  integer :: aelec, belec
  integer :: adets, bdets
  integer :: cidimension
! -------------------------------------------------------------------

! Read input file
  inputfl = 'detci.in'
  open(unit=1,file=inputfl,status='old',iostat=openstatus)
  if (openstatus > 0) stop "**** CANNOT OPEN INPUT FILE. ****"
  read(unit=1,nml=runlevel1)
  close(unit=1)

! Create and begin writing to output file
  outputfl = 'detci.out'
  open(unit=2,file=outputfl,status='new',iostat=openstatus)
  if (openstatus > 0) stop "**** CANNOT OPEN OUTPUT FILE. ****"
  
! Write header
  write(unit=2,fmt=9) "========================================================================"
  write(unit=2,fmt=9) "--------------------------DETERMINANT CI--------------------------------"
  write(unit=2,fmt=9) "========================================================================"
  write(unit=2,fmt=9) " "
  write(unit=2,fmt=9) " By Christopher L. Malbon"
  write(unit=2,fmt=9) " Yarkony Group"
  write(unit=2,fmt=9) " Dept. of Chemistry, The Johns Hopkins University"
  write(unit=2,fmt=9) " "
  write(unit=2,fmt=9) " "
  write(unit=2,fmt=9) "            ---------- Created on GENTOO LINUX ----------               "
  write(unit=2,fmt=9) "                           www.gentoo.org                               "
  write(unit=2,fmt=9) "------------------------------------------------------------------------"
  write(unit=2,fmt=9) "Send bug reports to cmalbon1@jhu.edu ------------"
  write(unit=2,fmt=9) "------------------------------------------------------------------------"
  write(unit=2,fmt=9) "                             RUN LEVEL 1                                "
  write(unit=2,fmt=9) " Namelist input: "
  write(unit=2,fmt=11)"  Electrons in system: ", electrons
  write(unit=2,fmt=11)"  Orbitals in system:  ", orbitals
  write(unit=2,fmt=11)"  Frozen orbitals:     ", nfrozen
  write(unit=2,fmt=11)"  DOCC   orbitals:     ", ndocc
  write(unit=2,fmt=11)"  Active orbitals:     ", nactive
  write(unit=2,fmt=11)"  Excitation level:    ", xlevel
  write(unit=2,fmt=9) " "
9  format(1x, A)
11 format(1x, A, I10)

! Call alphbet to find the correct numbers of alpha and beta electrons
  call alphbet( electrons, aelec, belec )
  write(unit=2,fmt=11)" Alpha electrons: ", aelec
  write(unit=2,fmt=11)" Beta  electrons: ", belec
  write(uint=2,fmt=9) " "

! Compute adets, bets
  adets = binom( orbitals, aelec )
  bdets = binom( orbitals, belec )
  write(unit=2,fmt=11)" Alpha determinants: ", adets
  write(unit=2,fmt=11)" Beta  determinants: ", bdets
  write(unit=2,fmt=9) " "

  write(unit=2,fmt=11)" Size of Full-CI expansion:      ", adets*bdets

! Truncate the CI space by calling citrunc
  call citrunc( adets, bdets, aelec, belec, orbitals, nfrozen, &
                ndocc, nactive, xlevel, tadetslen, tbdetslen, cidimension )

  write(unit=2,fmt=11)" Size of Truncated-CI expansion: ", cidimension

! Generate input for Run Level 2
  write(unit=2,fmt=9)" "
  write(unit=2,fmt=9)" Generating the necessary input for RUN LEVEL 2... "
  write(unit=2,fmt=9)"  Files: dets.list*, alpha.dets*, beta.dets*, input.2 "
  write(unit=2,fmt=9)"       *These files are created in the space-truncating subroutine."
  
  open(unit=3,file='input.2',status='new',iostat=openstatus)
  if ( openstatus > 0 ) stop "**** CANNOT RUN LEVEL 2 INPUT FILE! ****"
  write(unit=3,fmt=9) "&expinfo"
  write(unit=3,fmt=9) " moflnm='moints'"
  write(unit=3,fmt=11)" cidim =",cidimension
  write(unit=3,fmt=11)" aelec =",aelec
  write(unit=3,fmt=11)" belec =",belec
  write(unit=3,fmt=11)" adets =",adets
  write(unit=3,fmt=11)" bdets =",bdets
  write(unit=3,fmt=11)" tadetslen=",tadetslen
  write(unit=3,fmt=11)" tbdetslen=",tbdetslen
  write(unit=3,fmt=9) "/"
  close(unit=3)
  
! Run Level 1 has completed
  write(unit=2,fmt=9) " "
  write(unit=2,fmt=9) "                        RUN LEVEL 1 COMPLETED                          "
  close(unit=2)

! Exit
end program driver1

