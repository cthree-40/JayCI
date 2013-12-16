program run_level2
!=====================================
! Computes a determinant CI (Run_level 2)
! 
! By Christopher L. Malbon
! Yarkony Group
! The Johns Hopkins University
!=====================================
  use detci5
  use detci2
 ! use detci6
  implicit none
  character*20 :: inputfl, outputfl, r3input, activein
!
  character*20 :: moflnm
  integer :: live, aelec, belec, orbtls, adets, bdets, max1e, max2e, &
    nfrzn, ndocc, ncas, xlevel, krmin, krmax, roots, maxiter
  real*8 :: rtol
  integer,dimension(:,:),allocatable :: alphmat, betamat
  integer,dimension(:),allocatable :: nactstr
  real*8,dimension(:),allocatable :: dgls, moints1, moints2
  real*8,dimension(:,:),allocatable :: outvec, initial
  real*8,dimension(:),allocatable :: outval
  integer,parameter :: type1=1
  namelist /runlevel2/ moflnm, live, aelec, belec, orbtls, adets, &
                   bdets, rtol, krmin, krmax, roots, nfrzn, ndocc, ncas, &
                   xlevel, maxiter
!
  integer :: openstat, i, j, totdets,totels!, index2e2, index2e, ind2val
  integer :: initdim
!----------------------------------------------------------
!

#ifdef RUN2
  print *, "We are debugging runlevel2."
#endif

! # Read input file
  inputfl = 'r2.in'
  open(unit=1,file=inputfl,status='old',iostat=openstat)
  if (openstat > 0) stop "*** Cannot open input file ./r2.in ***"
  read(unit=1,nml=runlevel2)
! # Create and begin writing to output file
  close(unit=1)
  outputfl = 'detci.out'
  open(unit=2,file=outputfl,status='old',position='append',iostat=openstat)
  if (openstat > 0) stop "*** Cannot open output file ./detci.out ***"
  write(unit=2,fmt=9) "------------------------------- "
  write(unit=2,fmt=9) "******************************* "
  write(unit=2,fmt=9) "--------- Run Level 2 --------- "
  write(unit=2,fmt=9) "******************************* "
9  format(1x,A)
10 format(1x, A, A)
11 format(1x, A, I10)
12 format(1x, A, ES10.1)
  write(unit=2,fmt=9)  "Values read in from input file."
  write(unit=2,fmt=10) " Integral file name:              ", moflnm
  write(unit=2,fmt=11) " Number of 'active' determinants: ", live
  write(unit=2,fmt=11) " Alpha electrons:                 ", aelec
  write(unit=2,fmt=11) " Beta electrons:                  ", belec
  write(unit=2,fmt=11) " Number of orbitals               ", orbtls
  write(unit=2,fmt=11) " Number of alpha determinants:    ", adets
  write(unit=2,fmt=11) " Number of beta determinants:     ", bdets
  write(unit=2,fmt=12) " CI convergence tolerance:        ", rtol
  write(unit=2,fmt=11) " Minimum dimension of CI space:   ", krmin
  write(unit=2,fmt=11) " Maximum dimension of CI space:   ", krmax
  write(unit=2,fmt=11) " Maximum number of CI iterations: ", maxiter
! # Generate alpha and beta string index matrices
  write(unit=2,fmt=9) "------------------------------- "
!  write(unit=2,fmt=9) "Allocating Alpha/Beta string matrices"
!  allocate(alphmat(adets,aelec))
!  allocate(betamat(bdets,belec))
!  write(unit=2,fmt=9) "Calling strfnd"
!  call strfnd( adets, aelec, orbtls, adets, alphmat)
!  call strfnd( bdets, belec, orbtls, bdets, betamat)
!  print *, alphmat(1001,1), alphmat(1001,1)
!  print *, alphmat(1002,1), alphmat(1002,2)
!  print *, betamat(1001,1), betamat(1001,1)
!  print *, betamat(1002,1), betamat(1002,1)
! # Read in active determinant array
  activein='actdets.in'
  write(unit=2,fmt=9) "------------------------------- "
  write(unit=2,fmt=9) "Reading in active determinant string and allocating nactstr(i)"
  open(unit=3,file=activein,status='old',iostat=openstat)
  if (openstat > 0) stop "*** Cannot open ./actdets.in' ***"
  allocate(nactstr(live))
  read(unit=3,fmt=13) nactstr
  close(unit=3)
13 format(1x,I10)
! # Read in 1-e  and 2-e integrals
  write(unit=2,fmt=9) "------------------------------- "
! # Allocate moints1 and moints2
  write(unit=2,fmt=9) "Allocating moints1 and moints2 "
  max1e = ind2val(orbtls,orbtls)
  max2e = index2e2(orbtls,orbtls,orbtls,orbtls)
  write(unit=2,fmt=11)"  Length of moints1:         ", max1e
  write(unit=2,fmt=11)"  Length of moints2:         ", max2e
  allocate(moints1(max1e))
  allocate(moints2(max2e))
  do i=1, max1e
    moints1(i) = 0d0
  end do
  do i=1, max2e
    moints2(i) = 0d0
  end do
  write(unit=2,fmt=9) "Calling iwfmt.f"
  call iwfmt( moints1, moints2, type1, orbtls, moflnm, &
       max1e, max2e)    
! # Compute the diagonal matrix elements and hold them
  write(unit=2,fmt=9) "------------------------------- "
  write(unit=2,fmt=9) "Computing diagonal matrix elements"
  totels=aelec+belec
  allocate(dgls(live))

#ifdef RUN2
  print *, "Computing diagonal elements"
#endif
  call diagonal( nactstr, live, totels, aelec, belec, orbtls, adets, &
                 bdets, max1e, max2e, moints1, moints2, dgls )

!  do i = 1, live
!
!#ifdef RUN2_DIAG
!  write(*,fmt=11) "Currently computing diagonal element", i
!!#endif
!   print *, nactstr(i)
!   print *, nactstr(i+1)
!   dgls(i) = diagonal(nactstr(i), totels, aelec, belec, orbtls, alphmat, betamat, &
!     adets, bdets, max1e, max2e, moints1, moints2)
!  end do

  write(unit=2,fmt=9) "------------------------------- "
!  write(unit=2,fmt=9) "Deallocating alphmat and betamat."
!  write(unit=2,fmt=9) "Why? Because they will be generated in each action subroutine"
!  deallocate(alphmat)
!  deallocate(betamat)
! # Generate initial vector guess
  write(unit=2,fmt=9) "------------------------------- "
  write(unit=2,fmt=9) "Generating initial vector guess "
  write(unit=2,fmt=9) "Calling initvec"
!  call initvec(
  write(unit=2,fmt=9) " Currently, we will only use the HF det. for an initial guess - CLM 6/11/13"
#ifdef NOPREDIAG
  allocate(initial(live,krmin))
  do i = 1, krmin
    do j=1, live
      initial(j, i) = 0d0
      if ( j == i ) then
        initial(j,i) = 1d0
      end if
    end do  
  end do
#endif
! # Perform Davidson
  write(unit=2,fmt=9) "------------------------------- "
  write(unit=2,fmt=9) "Performing davidson algorithm. Calling davidson"
  close(unit=2)
  allocate(outval(roots))

#ifdef RUN2
  print *, "Calling davidson subroutine."
#endif
  allocate(outvec(live,roots))
#ifdef NOPREDIAG
  call davidson( initial, krmin, live, moints1, moints2, max1e, max2e, nactstr,  &
                 rtol, aelec, belec, totels, orbtls, adets, bdets, ndocc,       &
                 nfrzn, ncas, xlevel, roots, krmin, krmax, dgls, maxiter,&
                 outval, outvec )
#endif
#ifdef PREDIAG
  initdim = 4
  call davidson( initdim, live, moints1, moints2, max1e, max2e, nactstr, &
                 rtol, aelec, belec, totels, orbtls, adets, bdets, ndocc,&
                 nfrzn, ncas, xlevel, roots, krmin, krmax, dgls, maxiter,&
                 outval, outvec )
#endif
!  call davidson( initial, dgls, totels, aelec, belec, orbtls, adets, bdets, &
!						max1e, max2e, moints1, moints2, ndocc, nfrzn, ncas,     &
!						xlevel, nactstr, live, rtol, roots, krmin, krmax, maxiter,  &
!						outval)
  open(unit=2,file=outputfl,status='old',position='append',iostat=openstat)
  if (openstat > 0) stop "*** Cannot open output file ./detci.out ***"
  write(unit=2,fmt=12) "CI energy", outval(1)
  close(unit=2)
end program
  
