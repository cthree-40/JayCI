program run_level1
!==========================================================
! Computes a determinant CI (Run_level1)
! 
! By Christopher L. Malbon
! Yarkony Group
! The Johns Hopkins University
!==========================================================
! # Modules
  use detci1 !binom,fact,alphbet,lwval
  use detci2 !strfnd,adrfnd
  use detci3new !frzn,docc,active
!==========================================================
! Inputs:
!	inptfl   = input file name
!	outpfl   = output file name
!	moflnm   = MO integral file name
! 	totels   = total electrons
!	orbtls   = total orbitals
!	type1    = 1
!	krmin    = minium dimension of krylov space
!	krmax    = maximimum dimension of krylov space
! 	openstat = status of file
! nfrzn    = number of frozen orbitals
!	ndocc     = number of docc orbitals
! ncas     = number of cas orbitals
!	xlev     = excitation level in CI
! 	roots    = number of roots to compute
! Other variables
!	dgls     = < K | H | K >
!	moints1  = 1e- integrals
!   moints2  = 2e- integrals
!   alphmat  = alpha string matrix
!   betamat  = betastring matrix
!	invecs	 = initial vector guess
!   fnlvecs  = output vectors from davidson
!   actstr   = string of active determinants
!   nactstr  = string of active determinants (final)
!   aelec	 = number of alpha electrons
!   belec    = number of beta electrons
!   adets    = number of alpha determinants
!   bdets    = number of beta determinants
!   totdets  = total number of determinants
!   lactstr  =
!----------------------------------------------------------
  implicit none
  character*20 :: inptfl, outpfl, r2input, actdets
!
  character*20 :: moflnm
  integer :: totels, orbitals, type1, krmin, krmax, &
             openstat, nfrzn, ndocc, ncas, xlevel, roots, maxiter
  real*8 :: rtol
  namelist /general/ totels, orbitals, moflnm, type1, krmin,  &
                     krmax,  rtol,  nfrzn, ndocc,   ncas, xlevel, roots, maxiter
  real*8,dimension(:),allocatable :: dgls
  real*8,dimension(:,:),allocatable :: moints1
  real*8,dimension(:,:,:,:),allocatable :: moints2
  integer,dimension(:,:),allocatable :: alphmat,betamat
  real*8,dimension(:,:),allocatable :: invecs
  real*8,dimension(:,:),allocatable :: fnlvec
  integer,dimension(:),allocatable :: frzstr,dcclist,actlist
  integer :: aelec, belec, adets, bdets, totdets
  integer :: i, j, k, l, lnils, live, npa, npb,dcclen,actlen 
!----------------------------------------------------------
!

! Read input file
  inptfl = 'detci.in'
  open(unit=1,file=inptfl,status='old',iostat=openstat)
  if (openstat>0) stop "*** Cannot open input file. ***"
  read(unit=1,nml=general)

! Create and begin writing to output file
  outpfl = 'detci.out'
  open(unit=2,file=outpfl,status='new',iostat=openstat)
  if (openstat>0) stop "*** Cannot open output file. ***"

  write(unit=2,fmt=9) "------------------------ Determinant CI Program ------------------------ "
  write(unit=2,fmt=9) " By Christopher L. Malbon "
  write(unit=2,fmt=9) " Yarkony Group"
  write(unit=2,fmt=9) " Dept. of Chemistry, The Johns Hopkins University"
  write(unit=2,fmt=9) "------------------------------------------------------------------------- "
  write(unit=2,fmt=9) "************************************************************************* "
  write(unit=2,fmt=9) "----------------------------- Run Level 1 ------------------------------- "
  write(unit=2,fmt=9) "************************************************************************* "
  write(unit=2,fmt=11) " Electrons in system:  ", totels
  write(unit=2,fmt=11) " Orbitals  in system:  ", orbitals
  write(unit=2,fmt=10) " Integral file name:       ", moflnm
  write(unit=2,fmt=11) " Frozen orbitals:      ", nfrzn
  write(unit=2,fmt=11) " DOCC orbitals:        ", ndocc
  write(unit=2,fmt=11) " CAS orbitals:         ", ncas
  write(unit=2,fmt=11) " Excitation level:     ", xlevel
9  format(1x, A)
10 format(1x, A, A)
11 format(1x, A, I10)
  write(unit=2,fmt=9) "------------------------------------------------------------------------- "
  write(unit=2,fmt=9) "Calling alphbet (detci1)"
! # Compute alpha and beta electron numbers
  call alphbet( totels, aelec, belec )
  write(unit=2,fmt=11) "alpha electrons:  ", aelec
  write(unit=2,fmt=11) "beta  electrons:  ", belec
! # Compute total number of determinants
  adets  = binom(orbitals,aelec)
  bdets  = binom(orbitals,belec)
  totdets= adets * bdets
  write(unit=2,fmt=9) "------------------------------------------------------------------------- "
  write(unit=2,fmt=11) "alpha determinants:  ", adets
  write(unit=2,fmt=11) "beta  determinants:  ", bdets
  write(unit=2,fmt=11) "total determinants:  ", totdets
! # Generate alpha and beta string index matrices
  write(unit=2,fmt=9) "------------------------------------------------------------------------- "
  write(unit=2,fmt=9) "Allocating Alpha/Beta string matrices"
  allocate(alphmat(adets,aelec))
  allocate(betamat(bdets,belec))
  write(unit=2,fmt=9) "Calling strfnd"
  call strfnd( adets, aelec, orbitals, adets, alphmat)
  call strfnd( bdets, belec, orbitals, bdets, betamat)
! # Generate strings contained by active space def.
  write(unit=2,fmt=9) "------------------------------------------------------------------------- "
  write(unit=2,fmt=9) "Calling frozen to enforce frozen orbital restriction."
  npa = binom(orbitals-nfrzn,aelec-nfrzn)
  npb = binom(orbitals-nfrzn,belec-nfrzn)
  allocate(frzstr(npa*npb))
  call frozen(nfrzn,aelec,belec,totels,orbitals,adets,bdets,alphmat,&
							betamat,npa*npb,npa,npb,frzstr)
	write(unit=2,fmt=13)"CI truncated by  ", (adets*bdets)-(npa*npb), " terms."
13 format(1x,A,I10,A)
	write(unit=2,fmt=9) "------------------------------------------------------------------------- "
	write(unit=2,fmt=9) "Calling docc to enforce docc orbital restriction."
	call docc(ndocc,xlevel,nfrzn,aelec,belec,totels,orbitals,adets,bdets,alphmat,&
						betamat, frzstr,npa*npb,dcclist,dcclen)
	write(unit=2,fmt=13)"CI further truncated by  ", (npa*npb)-dcclen," terms."
	write(unit=2,fmt=9) "------------------------------------------------------------------------- "
	write(unit=2,fmt=9) "Calling active to enforce active space restriction"
	call active(ncas,ndocc,nfrzn,xlevel,aelec,belec,totels,orbitals,adets,bdets,&
							alphmat,betamat,dcclist,dcclen,actlist,actlen)
	write(unit=2,fmt=13)"CI further truncated by  ", dcclen - actlen," terms."
	write(unit=2,fmt=13)"Size of CI expansion is  ", actlen, " determinants."
	write(unit=2,fmt=9) "------------------------------------------------------------------------- "
! # Remove 0's from actsr
  write(unit=2,fmt=9) "Writing actlen(:) to the file actdets.in for use in run level 2."
  actdets = 'actdets.in'
  open(unit=4,file=actdets,status='new',iostat=openstat)
  if (openstat > 0) stop "*** Cannot open ./actdets.in ***"
  do i=1, actlen
	  write(unit=4,fmt=12) actlist(i)
  end do
  close(unit=4)
  write(unit=2,fmt=9) "------------------------------------------------------------------------- "
  write(unit=2,fmt=9) "Writing file r2.in"
  r2input='r2.in'
  open(unit=5,file=r2input,status='new',iostat=openstat)
  if (openstat > 0) stop "*** Cannot open ./r2.in ***"
  write(unit=5,fmt=9) "&runlevel2"
  write(unit=5,fmt=9) "moflnm='moints'"
  write(unit=5,fmt=11) "live=", actlen
  write(unit=5,fmt=11) "aelec=", aelec
  write(unit=5,fmt=11) "belec=", belec
  write(unit=5,fmt=11) "orbtls=", orbitals
  write(unit=5,fmt=11) "adets=", adets
  write(unit=5,fmt=11) "bdets=", bdets
  write(unit=5,fmt=14) "rtol=", rtol
  write(unit=5,fmt=11) "krmin=", krmin
  write(unit=5,fmt=11) "krmax=", krmax
  write(unit=5,fmt=11) "roots=", roots
  write(unit=5,fmt=11) "nfrzn=", nfrzn
  write(unit=5,fmt=11) "ndocc=", ndocc
  write(unit=5,fmt=11) "ncas=", ncas
  write(unit=5,fmt=11) "xlevel=", xlevel
  write(unit=5,fmt=11) "maxiter=",maxiter
  write(unit=5,fmt=9) "/"
12 format(1x,I10)
14 format(1x,A,1pe20.2)
  write(unit=2,fmt=9) "------------------------------------------------------------------------- "
  write(unit=2,fmt=9) "------------------------ Run Level 1 Complete. -------------------------- "
  write(unit=2,fmt=9) "------------------------------------------------------------------------- "
end program run_level1
