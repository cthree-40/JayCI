program test
! PROGRAM TO TEST VARIOUS SUBROUTINES
  use detci2
  use detci5
  use prediag
  implicit none

  integer :: adets, bdets, aelec, belec, orbitals, nfrozen, &
             ndocc, nactive, xlevel, cidimension, adetlen, bdetlen, &
             moints1len, moints2len, tadetslen, tbdetslen, cidim, &
             openstat, i, j, k, l
  character*20 :: inputfl1, inputfl2, outputfl, astringfl, bstringfl, &
                  determfl, pxreffile, pstepfl, qstepfl, pstringfl, qstringfl,&
                  moflnm, qxreffile
#ifdef GENORB
  integer, dimension(:), allocatable :: string1
  integer :: h
#endif
#ifdef HVDIAG
  integer, dimension(:), allocatable :: tadets, pstep, tbdets,plocate,qlocate,&
                                        qstep, tdets, qcrossref, pcrossref
  integer, dimension(:,:),allocatable::pstring,qstring
  real*8, dimension(:), allocatable  :: moints1, moints2
  real*8, dimension(:), allocatable :: diagonals
#endif
#ifdef HVACTION
  real*8, dimension(:), allocatable :: initialguess, resultingvec
  integer :: m
#endif
#ifdef ORBSTRINGS
  integer, dimension(3) :: pstring1, pstring2
  integer, dimension(3) :: qstring1, qstring2
  integer :: diffs
#endif
#ifdef STRINGDIFFS
  integer, dimension(:,:), allocatable :: diff_mat
#endif
#ifdef EXPBUILD
  real*8, dimension(:,:), allocatable :: hamiltonian
  real*8 :: abstol, dlamch
  integer :: lwork, liwork, il, iu, info, o, vu, vl
  integer, dimension(:),   allocatable :: isuppz
  real*8,  dimension(:),   allocatable :: work, iwork
  real*8,  dimension(:,:), allocatable :: eigvec
  real*8,  dimension(:),   allocatable :: eigval
#endif
#ifdef HVBUILDMAT
  real*8, dimension(:,:), allocatable :: unitmat
#endif
#ifdef GENPROBEL
  integer, dimension( 28, 2 ) :: probstr
  integer, dimension(1) :: apstring, bpstring
#endif
#if PREDIAG=2
  integer :: initguessdim1, initdiagsguess, dim1, dim2, dim3
  real*8, dimension(:,:),  allocatable :: outputvectors1, test_mat, test_vecs, &
                                          mat_a, mat_b, mat_out
#endif
#if PREDIAG=3
  real*8, dimension(:,:), allocatable :: outvectors
  integer :: subblockdim, initgdim
#endif
  integer :: testcol    ! Testing column for Hv
  integer :: probdet, probb, probp
!-------------------------
  testcol = 729

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
   h = 5
   bdets = 20
   allocate(string1(belec))
   call genorbstring(h,belec,orbitals,bdets,string1)
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
15 format(1x,F10.6)

   print *, "Testing hvdiag()..."
   print *, "System info:       "
   print *, "  H2"
   print *, "  Orbitals:                              28 "
   print *, "  Electrons:                              2 "
   print *, "  CI-dimension:                         784 "
   print *, " HF electronic energy should be:  -1.839941 "
   print *, " Hartree Fock GS is:              -1.132854 "
   print *, " Nuclear Repulsion energy is:      0.707107 " 

   adets       = 28
   bdets       = 28
   aelec       = 1
   belec       = 1
   cidimension = 784
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
!   allocate(pcrossref(cidim))
!   allocate(qcrossref(cidim))

! Read in determinants
   print *, " Reading in strings and determinants..."
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
!   open (unit=9,file=pxreffile,status='old',position='rewind',iostat=openstat)
!   if ( openstat > 0 ) stop "**** CANNOT OPEN CROSS REFERENCE FILE, pcross.ref. ****"
!   read(unit=9,fmt=13) (pcrossref(i), i=1, cidim )
!   close(unit=9)
!   open(unit=12,file=qxreffile,status='old',position='rewind',iostat=openstat)
!   if ( openstat > 0 ) stop "**** CANNOT OPEN CROSS REFERENCE FILE, qcross.ref. ****"
!   read(unit=12,fmt=13) (qcrossref(i), i=1, cidim )
!   close(unit=12)
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
                  adets, bdets, moints1, moints1len, moints2, moints2len, &
                  plocate, qlocate, diagonals )
!   print *, diagonals(1)
!   print *, diagonals(2)
!   print *, diagonals(3)
!   print *, diagonals(4)
#endif

#ifdef HVACTION
    print *, " Generating unit vector guess..."
    allocate(initialguess(cidim))
    allocate(resultingvec(cidim))
    initialguess=0d0
    initialguess(729)=1d0
    print *, " Calling acthv()..."
    call acthv( initialguess, moints1, moints2, moints1len, moints2len, pstring, pstep, plocate, &
                qstring, qstep, qlocate, cidim, tadets, tbdets, tadetslen, tbdetslen, adets, bdets,&
                aelec, belec, orbitals, diagonals, resultingvec )
    print *, "Resulting vector being printed to file: output_acthv.1"
    open( unit=15, file='output_acthv.1', status='new', position='rewind', iostat=openstat )
    if ( openstat .ne. 0 ) stop "**** COULD NOT OPEN output_acthv.1 ****"
    do m=1, cidim
      write( unit=15, fmt= 15 ) resultingvec(m)
    end do
    close( unit=15 )
    print *, "Done."
#endif
#ifdef ORBSTRINGS
    print *, "Testing orbdiffs..."
    do i=1, 3
      pstring1(i) = i
      pstring2(i) = i
      qstring1(i) = i
      qstring2(i) = i
    end do
    pstring2(3) = 5
    qstring2(3) = 6
    pstring2(2) = 4
    call orbdiffs( pstring1, pstring2, qstring1, qstring2, 3, 3, diffs )
    print *, diffs, " should equal ", 3
#endif
#ifdef STRINGDIFFS
    print *, "Testing stringdiffs..."
    do i=1, 3
      pstring1(i) = i
      pstring2(i) = i
      qstring1(i) = i
      qstring2(i) = i
    end do
    pstring2(2) = 4
    pstring2(3) = 5
    allocate( diff_mat(2,2) )
    call stringdiffs( pstring1, pstring2, 3, diff_mat, diffs )
    print *, diffs, " should equal ", 2
    print *, "and..."
    print *, diff_mat(1,1), diff_mat(1,2)
    print *, "should be equal."
    print *, "And below should be the excitation(s)..."
    print *, pstring1(2), pstring2(diff_mat(1,2))
    print *, pstring1(3), pstring2(diff_mat(2,2))
#endif
#ifdef EXPBUILD
    print *, "Building the Hamiltonian explicitly... "
    allocate( hamiltonian(cidim, cidim) )
    call exp_construct( moints1, moints1len, moints2, moints2len, cidim, aelec, belec, &
                    orbitals, tdets, hamiltonian )
    print *, " Hartree Fock GS: ", hamiltonian(1,1)
    open( unit=30, file='hamiltonian.col1', status='new', position='rewind', iostat = openstat )
    if ( openstat .ne. 0 ) stop "**** COULD NOT OPEN hamiltonian file ****"
    print *, " Printing column 1 of Hamiltonian... "
    do m=1, cidim
      write( unit=30, fmt=20 ) hamiltonian(m,1)
    end do
20 format(1x,F10.6,F10.6)
    close(unit=30)
    print *, "Diagonalizing this matrix..."
    print *, "Calling DSYEVR..."
    abstol = dlamch( 'safe minimum' )
    liwork = 12*cidim
    lwork  = 27*cidim
    allocate(isuppz(2*cidim))
    allocate(work(lwork))
    allocate(iwork(liwork))
    allocate(eigval(cidim))
    allocate(eigvec(cidim,cidim))
    open(unit=98, file='ham_debug.all', status='new')
    do i=1, cidim
      do j=1, cidim
        write(unit=98,fmt=98) j, i, hamiltonian(j,i)
      end do
    end do
98 format(1x,I4,I4,F10.6)
    call dsyevr( 'v','a','u', cidim, hamiltonian, cidim, vl, vu, il, iu, abstol, o, eigval, &
                  eigvec, cidim, isuppz, work, lwork, iwork, liwork, info )
    print *, "Eigenvalue 1: ", eigval(1)
    do i=1, cidim
      print *, eigvec(i,1)
    end do
#endif
#ifdef HVBUILDMAT
    print *, "Building Hamiltonian by performing Hv..."
    allocate(unitmat(cidim,cidim))
    unitmat=0d0
    do i=1, cidim
      unitmat(i,i) = 1d0
    end do
    print *, "Reallocating hamiltonian(:,:)..."
    deallocate(hamiltonian)
    allocate(hamiltonian(cidim,cidim))
    print *, "Peforming Hv_i = M_ij for i=1,2,...,cidim ..."
    do i=1, cidim
      print *, i
      call acthv( unitmat(1,i), moints1, moints2, moints1len, moints2len, pstring, pstep, plocate,&
                qstring, qstep, qlocate, tadets, tbdets, tadetslen, tbdetslen, adets, bdets,&
                aelec, belec, orbitals, diagonals, hamiltonian(1,i) )
    end do
    open( unit=90, file='hv-hamiltonian.col729', status='new', iostat=openstat )
    if ( openstat .ne. 0 ) stop "*** COULD NOT OPEN hv-hamiltonian.col729 ***"
    do i=1, cidim
      write(unit=90,fmt=20) hamiltonian(i,729), hamiltonian(729,i)
    end do
    close (unit=90)
    print *, "Diagonalizing constructed Hamiltonian..."
    deallocate(isuppz)
    deallocate(work)
    deallocate(iwork)
    deallocate(eigval)
    deallocate(eigvec)
    allocate(isuppz(2*cidim))
    allocate(work(lwork))
    allocate(iwork(liwork))
    allocate(eigval(cidim))
    allocate(eigvec(cidim,cidim))
    call dsyevr( 'n', 'a', 'u', cidim, hamiltonian, cidim, vl, vu, il, iu, abstol, o, eigval, &
                  eigvec, cidim, isuppz, work, lwork, iwork, liwork, info )
    print *, "Eigenvalue 1: ", eigval(1)
#endif
#ifdef GENPROBEL
    print *, "Generating strings of problem determinants..."
    do i=729, 756
      call k2indc( i, belec, orbitals, probstr(i-728,1), probstr(i-728,2) )
    end do
    print *, "Writing problem strings to file..."
    open( unit=60, file='problem.strings', status='new', iostat=openstat )
    if ( openstat .ne. 0 ) stop "**** CANNOT OPEN problem.strings ******"
    do i=1,28
      write( unit=60, fmt=60 ) probstr(i,1), probstr(i,2)
    end do
    close(unit=60)
60 format(1x,I5,I5)
    print *, "Writing problem alpha string arrays to file..."
    open( unit=61, file='problem.alpha', status='new', iostat=openstat )
    if ( openstat .ne. 0 ) stop "*** CANNOT OPEN problem.alpha ******"
    do i=1, 28
      call genorbstring( probstr(i,1), aelec, orbitals, adets, apstring )
      write( unit=61, fmt=61 ) apstring
    end do
    close( unit=61 )
61 format(1x,I10)
    print *, "Writing problem beta string arrays to file..."
    open( unit=62, file='problem.beta', status='new', iostat=openstat )
    if ( openstat .ne. 0 ) stop "*** CANNOT OPEN problem.beta ******"
    do i=1, 28
      call genorbstring( probstr(i,2), belec, orbitals, bdets, bpstring )
      write( unit=62, fmt=61 ) bpstring
    end do
    close( unit=62 )
    print *, " What determinant is causing problems? "
    read *, probdet
    call k2indc( probdet, belec, orbitals, probp, probb )
    print *, " This det is composed of strings: ", probp, probb
#endif    
!#if PREDIAG=2
!    print *, "Testing diagonalization subroutine..."
!    allocate( test_mat(4,4))
!   ! if ( allocated( test_mat) ) print *, "Allocated!"
!    allocate( test_vecs(4,4))
!   ! if ( allocated( test_mat) ) print *, "Allocated!"
!    do i=1,4
!      do j=1, 4
!        test_mat(i,j) = 1d0
!      end do
!    end do
!    print *, "Calling diag_hamsub..."
!    call diag_hamsub( test_mat, 4, 4, test_vecs )
!    print *, "Printing test_vecs"
!    do i=1, 4
!      print *, "________________________"
!      do j=1, 4
!  !      print *, test_vecs(j,i)
!      end do
!    end do
!    print *, "Calling ritz_vec..."
!    dim1 = 100
!    dim2 = 20
!    dim3 = 3
!    allocate( mat_b(dim2,dim2))
!    allocate( mat_a(dim1,dim2))
!    allocate( mat_out(dim1,dim3))
!    do i=1, 20
!      do j=1, 100
!        mat_a(j,i) = real(i)
!      end do
!    end do
!    mat_b = 0d0
!    do i=1, dim2
!      mat_b(i,i) = 1d0
!    end do
!    call ritz_vec( mat_b, dim1, dim2, dim3, mat_a, mat_out )
!    do i=1, dim3
    !  print *, mat_out(:,i)
!    end do
!    print *, "Testing lowdiag preconditioner subroutine..."
!    initguessdim1 = 10
!    initdiagsguess= 5
!    allocate(outputvectors1( cidim, initdiagsguess ))
!    call lowdiagprecond( diagonals, cidim, moints1, moints1len, moints2, &
!                         moints2len, pstring, pstep, plocate, qstring, qstep, &
!                         qlocate, tadets, tbdets, tadetslen, tbdetslen, &
!                         adets, bdets, aelec, belec, orbitals, initguessdim1, &
!                         initdiagsguess, outputvectors1 )
!    print *, "Test complete."
!    deallocate( outputvectors1 )
!#endif
!#if PREDIAG=3
!    subblockdim = 20
!    initgdim = 4
!    if (allocated( outvectors)) deallocate(outvectors)
!    allocate( outvectors( cidim, initgdim ) )
!    print *, "Call prediag...."
!    call prediagsubblock( cidim, moints1, moints2, moints1len, &
!     moints2len, tadets, tadetslen, tbdets, tbdetslen, tdets, subblockdim, initgdim,&
!     aelec, belec, orbitals, outvectors )    
!    print *, "PRINTING PREDIAG OUTPUT "
!    do i=1, initgdim
!      print *, "----------------------------------------"
!      do j=1, cidim
!        print *, outvectors(j,i)
!      end do
!    end do
!#endif
end program
