program test
! PROGRAM TO TEST VARIOUS SUBROUTINES
  use detci2
  use detci5
  implicit none

  integer :: adets, bdets, aelec, belec, orbitals, nfrozen, &
             ndocc, nactive, xlevel, cidimension, adetlen, bdetlen, &
             moints1len, moints2len, tadetslen, tbdetslen, cidim, &
             openstat
  character*20 :: inputfl1, inputfl2, outputfl, astringfl, bstringfl, &
                  determfl, pxreffile, pstepfl, qstepfl, pstringfl, qstringfl,&
                  moflnm, qxreffile
#ifdef GENORB
  integer, dimension(:), allocatable :: string1
  integer :: i
#endif
#ifdef HVDIAG
  integer, dimension(:), allocatable :: tadets, pstep, tbdets,plocate,qlocate,&
                                        qstep, tdets, qcrossref, pcrossref
  integer, dimension(:,:),allocatable::pstring,qstring
  real*8, dimension(:), allocatable  :: moints1, moints2
  real*8, dimension(:), allocatable :: diagonals
  integer ::i, j, k, l
#endif
!-------------------------

   print *, "Starting..."
#ifdef CITRUNC
   print *, "Testing citrunc()..."

   adets    = 20
   bdets    = 20
   aelec    = 3
   belec    = 3
   orbitals = 6
   nfrozen  = 1
   ndocc    = 2
   nactive     = 2
   xlevel   = 1

   call citrunc( adets, bdets, aelec, belec, orbitals, nfrozen,&
                 ndocc, nactive, xlevel, adetlen, bdetlen, cidimension ) 
   print *, cidimension
#endif
#ifdef GENORB
   belec = 3
   orbitals = 6
   i = 5
   bdets = 20
   allocate(string1(belec))
   call genorbstring(i,belec,orbitals,bdets,string1)
   print *, string1
#endif
#ifdef HVDIAG
  inputfl1 = 'detci.in'
  inputfl2 = 'input.2'
  outputfl = 'detci.out'
  astringfl= 'alpha.dets'
  bstringfl= 'beta.dets'
  determfl = 'det.list'
  pxreffile = 'pcross.ref'
  qxreffile = 'qcross.ref'
  pstepfl  = 'pstep.list'
  qstepfl  = 'qstep.list'
  pstringfl= 'pstring.list'
  qstringfl= 'qstring.list'
  moflnm   = 'moints'


9  format(1x,A)
10 format(1x,A,I10)
11 format(1x,A,A)
12 format(1x,A,ES10.1)
13 format(1x,I10)
14 format(1x,I10,I10)

   print *, "Testing hvdiag()..."
   print *, "System info:       "
   print *, "  H2"
   print *, "  Orbitals:                           28 "
   print *, "  Electrons:                           2 "
   print *, "  CI-dimension:                      784 "
   print *, " HF electronic energy should be:  -1.133 " 

   adets       = 28
   bdets       = 28
   aelec       = 1
   belec       = 1
   cidimension = 1575
   orbitals    = 28
   nfrozen     = 0
   ndocc       = 0
   nactive     = 10
   xlevel      = 2
   
   print *, " "
   print *, "Calling citrunc()..."
   call citrunc( adets, bdets, aelec, belec, orbitals, nfrozen, &
                 ndocc, nactive, xlevel, tadetslen, tbdetslen, cidim )
   print *, " "
   print *, "Reading in results of citrunc()..."
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
   allocate(pcrossref(cidim))
   allocate(qcrossref(cidim))

! Read in determinants
   write(unit=2,fmt=9) " Reading in strings and determinants..."
! Alpha strings
   open(unit=4,file=astringfl,status='old',position='rewind',iostat=openstat)
   if ( openstat > 0 ) stop "**** CANNOT OPEN ALPHA STRING FILE, alpha.dets. ****"
   read(unit=4,fmt=13) ( tadets(i), i=1, tadetslen )
   close(unit=4)

! Alpha string det list
   open(unit=7,file=pstringfl,status='old',position='rewind',iostat=openstat)
   if ( openstat > 0 ) stop "**** CANNOT OPEN PSTRING DETERMINANT LIST, pstring.list. ****"
   read(unit=7,fmt=14) ( (pstring(i,j), j=1, 2),i=1, cidim )
   close(unit=7)
! Alpha string step list
   open(unit=10,file=pstepfl,status='old',position='rewind',iostat=openstat)
   if ( openstat > 0 ) stop "**** CANOT OPEN PSTRING STEP FILE, pstep.list. ****"
   read(unit=10,fmt=13) ( pstep(i), i=1, tadetslen )
   close(unit=10)

! Beta strings
   open (unit=5,file=bstringfl,status='old',position='rewind',iostat=openstat)
   if ( openstat > 0 ) stop "**** CANNOT OPEN BETA STRING FILE, beta.dets. ****"
   read(unit=5,fmt=13) ( tbdets(i), i=1, tbdetslen )
   close(unit=5)

! Beta string det list
   open(unit=8,file=qstringfl,status='old',position='rewind',iostat=openstat)
   if ( openstat > 0 ) stop "**** CANNOT OPEN QSTRING DETERMINANT LIST, qstring.list. ****"
   read(unit=8,fmt=14) (( qstring(i,j), j=1,2), i=1, cidim )
   close(unit=8)

! Beta string step list
   open(unit=11,file=qstepfl,status='old',position='rewind',iostat=openstat)
   if ( openstat > 0 ) stop "**** CANNOT OPEN QSTRING STEP FILE, qstep.list. ****"
   read(unit=11,fmt=13) ( qstep(i), i=1, tbdetslen )
   close(unit=11)

! Determinants
   open (unit=6,file=determfl,status='old',position='rewind',iostat=openstat)
   if ( openstat > 0 ) stop "**** CANNOT OPEN DETERMINANT FILE, dets.list. ****"
   read(unit=6,fmt=13) ( tdets(i), i=1, cidim )
   close(unit=6)

! Read in cross reference list
   open (unit=9,file=pxreffile,status='old',position='rewind',iostat=openstat)
   if ( openstat > 0 ) stop "**** CANNOT OPEN CROSS REFERENCE FILE, pcross.ref. ****"
   read(unit=9,fmt=13) (pcrossref(i), i=1, cidim )
   close(unit=9)
   open(unit=12,file=qxreffile,status='old',position='rewind',iostat=openstat)
   if ( openstat > 0 ) stop "**** CANNOT OPEN CROSS REFERENCE FILE, qcross.ref. ****"
   read(unit=12,fmt=13) (qcrossref(i), i=1, cidim )
   close(unit=12)
! Read in locate lists
   open (unit=13,file='plocate.list',status='old',position='rewind',iostat=openstat)
   if ( openstat > 0 ) stop "**** CANNOT OPEN LOCATE FILE, plocate.list ****"
   read(unit=13,fmt=13) (plocate(i), i=1, tadetslen )
   close(unit=13)
   open (unit=14,file='qlocate.list',status='old',position='rewind',iostat=openstat)
   if (openstat > 0 ) stop "**** CANNOT OPEN LOCATE FILE, qlocate.list. ****"
   read(unit=14,fmt=13) (qlocate(i), i=1, tbdetslen )
   close(unit=14)

! Allocate MO arrays
   moints1len=ind2val(orbitals,orbitals)
   moints2len=index2e2(orbitals,orbitals,orbitals,orbitals)
   allocate(moints1(moints1len))
   allocate(moints2(moints2len))
   call iwfmt( moints1, moints2, 1, orbitals, moflnm, moints1len, moints2len )
! Call subroutine
   print *, " "
   print *, "Calling diagonal()..."
   allocate(diagonals(cidim))
   call diagonal( pstring, pstep, tadetslen, qstring, qstep, tbdetslen, tadets, &
                  tadetslen, tbdets, tbdetslen, tdets, cidim, 3, orbitals, aelec, belec,&
                  adets, bdets, moints1, moints1len, moints2, moints2len, pcrossref, qcrossref,cidim,&
                  plocate, qlocate, diagonals )
   print *, diagonals(1)
   print *, diagonals(2)
   print *, diagonals(3)
   print *, diagonals(4)
#endif
end program
