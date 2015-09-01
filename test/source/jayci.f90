program jayci
  !=============================================================================
  ! jayci
  ! -----
  ! Purpose: Determinant CI program.
  !
  ! Must execute jayci_exp.x before executation
  !
  ! Written by: cthree-40
  ! Dept. of Chemistry, The Johns Hopkins University
  !=============================================================================
  use input_proc,     only: alpha_beta
  use combinatorial,  only: binom
  use integral,       only: ind2val, index2e, readmoints
  use construct,      only: ham_element_diag
  use initialguess,   only: diagrefblock
  use orthogroutines, only: modgramschmidt
  implicit none

  ! -- IN/OUT FILES ------------------------------------------------------------
  ! ipfl1_name  = master input file name
  ! ipfl1_unit  = master input file unit
  ! ipfl2_name  = jayci_exp-generated input file name
  ! ipfl2_unit  = jayci_exp-generated input file unit
  ! outfl_name  = output file name
  ! outfl_unit  = output file unit
  character*25  :: ipfl1_name = "jayci.in", ipfl2_name = "input.jayci"
  character*25  :: outfl_name = "jayci.out"
  integer       :: ipfl1_unit = 10, ipfl2_unit = 11, outfl_unit = 12

  ! -- EXPANSION INPUT FILES ---------------------------------------------------
  ! dtlst_name  = determinant list file name
  ! dtlst_unit  = determinant list file unit
  character*25  :: dtlst_name = "det.list"
  integer       :: dtlst_unit

  ! -- general NAMELIST --------------------------------------------------------
  ! electrons = electrons in system
  ! orbitals  = orbitals in system
  ! nfrozen   = number of frozen core orbitals
  ! ndocc     = number of doubly-occupied orbitals
  ! nactive   = number of active orbitals
  ! xlevel    = excitation level
  ! nfrzvirt  = number of frozen virtual orbitals
  ! printlvl  = level of output printing:    0 = minimal
  !                                       1..5 = reserved
  !                                         >5 = debugging
  integer :: electrons, orbitals, nfrozen, ndocc, nactive, xlevel, nfrzvirt
  integer :: printlvl
  namelist /general/ electrons, orbitals, nfrozen, ndocc, nactive, xlevel, &
    printlvl, nfrzvirt

  ! -- expansion NAMELIST ------------------------------------------------------
  ! ci_electrons = electrons subject to excitations in CI expansion
  ! ci_orbitals  = orbitals  subject to excitations in CI expansion
  ! astr_len     = number of alpha strings in expansion
  ! bstr_len     = number of beta  strings in expansion
  ! dtrm_len     = number of determinants  in expansion
  integer :: ci_electrons, ci_orbitals, astr_len, bstr_len, dtrm_len
  namelist /expansion/ ci_electrons, ci_orbitals, astr_len, bstr_len, dtrm_len

  ! -- dalginfo NAMELIST -------------------------------------------------------
  ! maxiter   = maximum iterations
  ! krymin    = minimum dimension of krylov space
  ! krymax    = maximum dimension of krylov space
  ! restol    = convergence tolerance of residual
  ! nroots    = number of roots to converge
  ! prediagr  = prediagonalization routine
  ! iguessdim = dimension of initial guess
  ! refdim    = reference space dimension
  integer :: maxiter, krymin, krymax, nroots, prediagr, iguessdim, refdim
  real*8  :: restol
  namelist /dalginfo/ maxiter, krymin, krymax, nroots, prediagr, iguessdim, &
    restol, refdim

  ! .. LOCAL scalars ..
  ! m1len    = number of 1-e integrals
  ! m2len    = number of 2-e integrals
  ! ierr     = error handling
  ! ci_aelec = ci alpha electrons
  ! ci_belec = ci beta  electrons
  ! nuc_rep  = nuclear repulsion energy
  ! frz_core = frozen core energy
  ! astrings   = binom(ci_orbitals, ci_aelec)
  ! bstrings   = binom(ci_orbitals, ci_belec)
  integer :: m1len, m2len
  integer :: ierr
  integer :: ci_aelec, ci_belec
  integer :: i
  real*8  :: nuc_rep, frz_core, tot_frz
  integer :: astrings, bstrings
  
  ! .. LOCAL arrays ..
  ! moints1   = 1-e integrals
  ! moints2   = 2-e integrals
  ! eigvecs   = eigenvectors of hamiltonian
  ! eigvals   = eigenvalues  of hamiltonian
  ! diagvals  = diagonal values of hamiltonian
  ! dtrm_list = determinants in expansion
  ! astr_list = alpha string list
  ! bstr_list = beta  string list
  ! dtrm_alpha_listq = beta  string values when searching determinant list by
  !                    alpha string values
  ! dtrm_beta_listp  = alpha string values when searching determinant list by
  !                    beta  string values
  ! pstep      = number of beta  strings for each alpha string
  ! qstep      = number of alpha strings for each beta  string
  ! plocate    = location of ith alpha string
  ! qlocate    = locateion of ith beta string
  ! init_vecs  = initial vectors
  ! mofile     = molecular orbital file
  ! xreflist   = cross reference list for determinants
  real*8,  dimension(:),   allocatable :: moints1, moints2, eigvals, diagvals
  real*8,  dimension(:,:), allocatable :: eigvecs, init_vecs
  integer, dimension(:),   allocatable :: dtrm_list, astr_list, bstr_list
  integer, dimension(:),   allocatable :: dtrm_alpha_listq, dtrm_beta_listp
  integer, dimension(:),   allocatable :: pstep, qstep, plocate, qlocate
  integer, dimension(:),   allocatable :: xreflist
  character*25 :: mofile = "moints"
  
  ! .. LOCAL logicals ..
  ! file_exists = logical test for existance of necessary files
  logical :: file_exists

  ! default &general namelist values
  electrons = 0
  orbitals  = 0
  nfrozen   = 0
  ndocc     = 0
  nactive   = 0
  nfrzvirt  = 0
  xlevel    = 2 ! This won't return error. Default CI program is SDCI.
  printlvl  = 0 ! Minimal printing by default.

  ! default &expansion namelist values
  ci_electrons = 0
  ci_orbitals  = 0
  astr_len     = 0
  bstr_len     = 0
  dtrm_len     = 0

  ! default &dalginfo namelist values
  maxiter   =  25
  krymin    =   3
  krymax    =   8
  nroots    =   1
  prediagr  =   1
  iguessdim =   4
  refdim    = 500
  restol    = 1.0d-5

  ! open input file. read &general and &dalginfo namelists
  open(file = ipfl1_name, unit = ipfl1_unit, action = "read", status = "old", &
    iostat = ierr)
  if (ierr .ne. 0) stop "*** Could not open input file! ***"
  read(unit = ipfl1_unit, nml = general)
  read(unit = ipfl1_unit, nml = dalginfo)
  ! close input file.
  close(ipfl1_unit)

  ! open input file from jaci_exp.x. read &expansion namelist
  open(file = ipfl2_name, unit = ipfl2_unit, action = "read", status = "old", &
    iostat = ierr)
  if (ierr .ne. 0) stop "*** Could not open input file 2! ***"
  read(unit = ipfl2_unit, nml = expansion)
  ! close input file
  close(ipfl2_unit)
  
  ! open output file and echo namelist input
  open(file = outfl_name, unit = outfl_unit, action = "write", &
    status = "unknown", position = "rewind", iostat = ierr)
  if (ierr .ne. 0) stop "*** Could not open output file! ***"
  
  call write_header(outfl_unit)
  call write_gennml(electrons, orbitals, nfrozen, ndocc, nactive, xlevel, &
    outfl_unit)
  call write_expinml(ci_electrons, ci_orbitals, astr_len, bstr_len, dtrm_len, &
    outfl_unit)
  call write_dainml(maxiter, krymin, krymax, nroots, prediagr, iguessdim, &
    restol, refdim, outfl_unit)

  ! compute ci_aelec and ci_belec
  call alpha_beta(ci_electrons, ci_aelec, ci_belec)
  astrings = binom(ci_orbitals, ci_aelec)
  bstrings = binom(ci_orbitals, ci_belec)
  
  ! read in molecular orbitals
  m1len = ind2val(orbitals, orbitals)
  m2len = index2e(orbitals, orbitals, orbitals, orbitals)
  allocate(moints1(m1len))
  allocate(moints2(m2len))
  inquire(file = mofile, exist = file_exists)
  if (.not. file_exists) stop "*** No molecular integral file! ***"
  call readmoints(moints1, moints2, 1, ci_orbitals, mofile, m1len, m2len, &
    nuc_rep, frz_core)
  tot_frz = nuc_rep + frz_core
  write(outfl_unit, "(A)") ""
  write(outfl_unit, "(1x,A,f15.8)") "Nuclear repulsion energy:", nuc_rep
  write(outfl_unit, "(1x,A,f15.8)") "Frozen-core contribution:", frz_core
  write(outfl_unit, "(1x,A,f15.8)") "Total contribution:      ", tot_frz
  write(outfl_unit, "(A)") ""

  ! read in expansion information: dtrm_list, dtrm_alpha_listq, pstep
  allocate(dtrm_list(dtrm_len))
  allocate(pstep(astr_len))
  allocate(dtrm_alpha_listq(dtrm_len))
  call read_exp_info1(dtrm_list, dtrm_alpha_listq, pstep, dtrm_len, astr_len, &
    ierr)
  if (ierr .ne. 0) stop "*** Error reading expansion information! ***"

  ! read in expansion information: astr_list, bstr_list
  allocate(astr_list(astr_len))
  allocate(bstr_list(bstr_len))
  call read_exp_info2(astr_list, bstr_list, astr_len, bstr_len, ierr)

  ! generate dtrm_beta_listp, qstep and xreflist
  allocate(dtrm_beta_listp(dtrm_len))
  allocate(qstep(bstr_len))
  allocate(xreflist(dtrm_len))
  call gen_exp_info1(dtrm_beta_listp, qstep, bstr_list, astr_list, &
    dtrm_len, bstr_len, astr_len, xreflist, dtrm_list, ierr)

  ! generate qlocate and plocate
  allocate(plocate(astr_len))
  allocate(qlocate(bstr_len))
  call gen_locate(pstep, astr_len, plocate)
  call gen_locate(qstep, bstr_len, qlocate)
  
  ! compute diagonal matrix elements
  allocate(diagvals(dtrm_len))
  do i = 1, dtrm_len
          diagvals(i) = ham_element_diag(dtrm_list(i), moints1, m1len, moints2, &
            m2len, ci_aelec, ci_belec, ci_orbitals)
  end do
  write(outfl_unit, "(1x,A,2f15.8)") "Hartree Fock Energy: ", &
    (diagvals(1) + tot_frz)

  ! allocate initial vectors
  allocate(init_vecs(dtrm_len, iguessdim))

  ! perform prediagonalization routine
  if (prediagr .eq. 1) then
          ! diagonalization of a reference block of hamiltonian
          call diagrefblock(dtrm_len, refdim, moints1, m1len, moints2, m2len, &
            dtrm_list, ci_aelec, ci_belec, ci_orbitals, iguessdim, init_vecs, &
            tot_frz)

  else if (prediagr .eq. 3) then
          ! diagonalize ENTIRE hamiltonian. This is done to check Hv.
          call diagrefblock(dtrm_len, dtrm_len, moints1, m1len, moints2, m2len,&
            dtrm_list, ci_aelec, ci_belec, ci_orbitals, iguessdim, init_vecs,   &
            tot_frz)
  else
          ! unknown prediagonalization routine
          stop " *** Unknown prediagonalization routine! ***"
          
  end if

  ! orthogonalize initial guess vectors
  call modgramschmidt(init_vecs, iguessdim, 0, dtrm_len)

  ! call davidson algorithm
  allocate(eigvals(nroots))
  allocate(eigvecs(dtrm_len, nroots))

  call davidson(iguessdim, init_vecs, diagvals, moints1, m1len, moints2, m2len, &
    dtrm_len, dtrm_alpha_listq, plocate, pstep, dtrm_beta_listp, qlocate, qstep,&
    xreflist, astr_list, astr_len, bstr_list, bstr_len, astrings, bstrings,     &
    ci_aelec, ci_belec, ci_orbitals, krymin, krymax, restol, nroots, maxiter,   &
    nfrozen, ndocc, nactive, tot_frz, eigvals, eigvecs)

  ! write out final energy
  write(*,"(1x,'Final CI energy = ',f15.8)") (eigvals(1) + tot_frz)

  deallocate(eigvals, eigvecs, init_vecs)
contains

  subroutine cmp_arrays(array1, array2, length, xref)
    !===========================================================================
    ! cmp_arrays
    ! ----------
    ! Purpose: compare two integer arrays and return mapping of 2 -> 1.
    !
    ! Input:
    !  array1 = integer array
    !  array2 = integer array
    !  length = length of array1 and array2
    !
    ! Output:
    !  xref   = integer array mapping of 2 -> 1
    !---------------------------------------------------------------------------
    use search_fcns, only: int_search2
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: length
    integer, dimension(length), intent(in) :: array1, array2

    ! .. OUTPUT arguments ..
    integer, dimension(length), intent(out) :: xref

    ! .. LOCAL scalars ..
    integer :: i, loc

    ! loop over elements of array2
    do i = 1, length

            ! find location of array2(i) in array1
            call int_search2(array2(i), array1, length, loc)

            xref(i) = loc

    end do

    return
  end subroutine cmp_arrays

  subroutine gen_exp_info1(dblistp, qstep,  blist, alist, ndets, nbstr, &
    nastr, xrlist, dtlist, ierr)
    !===========================================================================
    ! gen_exp_info1
    ! -------------
    ! Purpose: generate list of p string pairings for searching determinant list
    !          by q strings, qstep list 
    !
    ! Input:
    !  ndets  = number of determinants
    !  nbstr  = number of beta strings
    !  alist  = alpha strings list
    !  blist  = beta strings list
    !  nastr  = number of alpha strings
    !  dtlist = determinant list
    ! Output:
    !  dblistp = p pairings for q strings
    !  qstep   = number of p pairings for each q
    !  xrlist  = cross-reference list
    !  ierr    = error handling
    !---------------------------------------------------------------------------
    use search_fcns, only: int_search2
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: ndets, nbstr, nastr
    integer, dimension(nbstr), intent(in) :: blist, alist
    integer, dimension(ndets), intent(in) :: dtlist
    
    ! .. OUTPUT arguments ..
    integer, dimension(ndets), intent(out) :: dblistp
    integer, dimension(nbstr), intent(out) :: qstep
    integer, dimension(ndets), intent(out) :: xrlist
    integer :: ierr
    
    ! .. LOCAL scalars ..
    integer :: p, q, pi, qi, deti
    integer :: i, j, iter, cntr
    
    ! .. LOCAL arrays ..
    integer, dimension(:,:,:), allocatable :: scr
    integer, dimension(:),     allocatable :: qdetlist
    
    ! allocate scratch array and qdetlist
    allocate(scr(nastr, nbstr, 2))
    allocate(qdetlist(ndets))
    
    scr = 0
    
    ! open det.list file
    open(file = "det.list", unit = 15, status = "old", action = "read", &
      iostat = ierr)
    if (ierr .ne. 0) return

    ! loop over determinants
    do i = 1, ndets
            read(15, 10) p, q, deti

            ! find p and q in determinant lists
            call int_search2(p, alist, nastr, pi)
            call int_search2(q, blist, nbstr, qi)

            scr(pi,qi,1) = p
            scr(pi,qi,2) = deti
    end do

    ! add p values in columns to dblistp and determinant indexes to qdetlist
    iter = 1
    do i = 1, nbstr
    
            cntr = 0
            do j = 1, nastr

                    if (scr(j, i, 1) .ne. 0) then
                            dblistp(iter) = scr(j, i, 1)
                            qdetlist(iter)= scr(j, i, 2)
                            iter = iter + 1
                            cntr = cntr + 1
                    end if
            end do
            qstep(i) = cntr ! add step quantity
    end do

    ! build cross reference list
    call cmp_arrays(dtlist, qdetlist, ndets, xrlist)

    deallocate(scr, qdetlist)
    return
10  format(i15,i15,i15)
  end subroutine gen_exp_info1

  subroutine gen_locate(step, nstr, locate)
    !===========================================================================
    ! gen_locate
    ! ----------
    ! Purpose: generate location array using step array.
    !
    ! Input:
    !  step  = array showing number of string pairings for each string
    !  nstr  = total number of strings
    !
    ! Output:
    !  locate = location of each string
    !---------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: nstr
    integer, dimension(nstr), intent(in) :: step

    ! .. OUTPUT arguments ..
    integer, dimension(nstr), intent(out) :: locate

    ! .. LOCAL scalars ..
    integer :: i

    locate(1) = 0

    do i = 2, nstr

            locate(i) = step(i-1) + locate(i - 1)

    end do

    return
  end subroutine gen_locate

  subroutine read_exp_info1(dlist, dalistq, pstep, dlist_len, pstp_len, ierr)
    !===========================================================================
    ! read_exp_info1
    ! --------------
    ! Purpose: Read in determinant list, bulid alpha string step list and
    !          read in determinant list q values.
    !
    ! Input:
    !  dlist_len, pstep_len
    !
    ! Output:
    !  dlist, dalistq, pstep, ierr
    !---------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: dlist_len, pstp_len

    ! .. OUTPUT arguments ..
    integer, intent(out) :: ierr
    integer, dimension(dlist_len), intent(out) :: dlist, dalistq
    integer, dimension(pstp_len),  intent(out) :: pstep

    ! .. LOCAL scalars ..
    integer :: dtlst_flun = 15
    integer :: i
    integer :: pi, ai1, ai2, cntr
    
    ! .. LOCAL arrays ..
    character*25 :: dtlst_flnm = "det.list"
    
    ! open det.list file
    open(file = dtlst_flnm, unit = dtlst_flun, status = "old", action = "read",&
      iostat = ierr)
    if (ierr .ne. 0) return

    ! read in arrays
    cntr = 0
    pi   = 1
    ai2  = 1
    ai1  = 1
    ierr = 0
    do i = 1, dlist_len
            read(dtlst_flun, 10, iostat = ierr) ai2, dalistq(i), dlist(i)
            cntr = cntr + 1

            ! test if ai has changed
            if (ai2 .ne. ai1) then
                    pstep(pi) = cntr - 1
                    ! iterate alpha string index
                    pi = pi + 1
                    ! reset counter to 1 because we just read in a new p
                    cntr = 1
            end if
            ai1 = ai2
    end do
    ! last element
    pstep(pi) = cntr

    return
10  format(i15,i15,i15)
  end subroutine read_exp_info1

  subroutine read_exp_info2(alist, blist, alen, blen, ierr)
    !==========================================================================
    ! read_exp_info2
    ! --------------
    ! Purpose: Read alist and blist from str.list file.
    !--------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: alen, blen

    ! .. OUTPUT arguments ..
    integer, dimension(alen), intent(out) :: alist
    integer, dimension(blen), intent(out) :: blist
    integer, intent(out) :: ierr

    ! .. LOCAL scalars ..
    integer :: i
    integer :: str_unit = 15

    ! .. LOCAL arrays ..
    character*25 :: str_name = "str.list"
    character*5  :: spin, astr = "ALPHA", bstr = "BETA"
    ! open str.list file
    open(file = str_name, unit = str_unit, status = "old", action = "read", &
      iostat = ierr)
    if (ierr .ne. 0) return

    ! read first line
    read(str_unit, 10)
    ! read in beta strings
    do i = 1, blen
            read(str_unit, 11) blist(i)
    end do
    ! read "ALPHA" line
    read(str_unit, 10)
    do i = 1, alen
            read(str_unit, 11) alist(i)
    end do
    
    return
10  format(a5)
11  format(i15)
  end subroutine read_exp_info2

  subroutine write_dainml(miter, kmin, kmax, nrts, pdr, igd, rstl, refdim, &
    outfl_unit)
    !===========================================================================
    ! write_dainml
    ! ------------
    ! Purpose: Write &dalginfo namelist input values
    !---------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: miter, kmin, kmax, nrts, pdr, igd, refdim, outfl_unit
    real*8,  intent(in) :: rstl
    
    ! write input values
    write(outfl_unit, 10) ""
    write(outfl_unit, 10) "&dalginfo input values:"
    write(outfl_unit, 11) "maxiter      = ", miter
    write(outfl_unit, 11) "krymin       = ", kmin
    write(outfl_unit, 11) "krymax       = ", kmax
    write(outfl_unit, 11) "nroots       = ", nrts
    write(outfl_unit, 11) "prediagr     = ", pdr
    write(outfl_unit, 11) "iguessdim    = ", igd
    write(outfl_unit, 12) "restol       = ", rstl
    write(outfl_unit, 11) "refdim       = ", refdim
    write(outfl_unit, 10) ""
    
    return
10  format(1x,a)
11  format(2x,a,i15)
12  format(2x,a,f15.8)
  end subroutine write_dainml

  subroutine write_expinml(elec, orb, alen, blen, dlen, outfl_unit)
    !===========================================================================
    ! write_expinml
    ! -------------
    ! Purpose: Write &expansion namelist input values
    !---------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: elec, orb, alen, blen, dlen, outfl_unit

    ! write input values
    write(outfl_unit, 10) ""
    write(outfl_unit, 10) "&expansion input values:"
    write(outfl_unit, 11) "ci_electrons = ", elec
    write(outfl_unit, 11) "ci_orbitals  = ", orb
    write(outfl_unit, 11) "alpha strings= ", alen
    write(outfl_unit, 11) "beta  strings= ", blen
    write(outfl_unit, 11) "determinants = ", dlen
    write(outfl_unit, 10) ""
    
    return
10  format(1x,a)
11  format(2x,a,i15)
  end subroutine write_expinml

  subroutine write_gennml(elec, orb, nfrzn, ndocc, nactv, xlvl, outfl_unit)
    !===========================================================================
    ! write_gennml
    ! ------------
    ! Purpose: Write &general namelist input values.
    !
    ! elec  = number of electrons
    ! orb   = number of orbitals
    ! nfrzn = number of frozen core orbitals
    ! ndocc = number of docc orbitals
    ! nactv = number of active orbitals
    ! xlvl  = excitation level
    !---------------------------------------------------------------------------
    implicit none
    
    ! .. INPUT arguments ..
    integer, intent(in) :: outfl_unit, elec, orb, nfrzn, ndocc, nactv, xlvl

    ! write input values
    write(outfl_unit, 10) ""
    write(outfl_unit, 10) "&general input values:"
    write(outfl_unit, 11) "electrons =", elec
    write(outfl_unit, 11) "orbitals  =", orb
    write(outfl_unit, 11) "nfrozen   =", nfrzn
    write(outfl_unit, 11) "ndocc     =", ndocc
    write(outfl_unit, 11) "nactive   =", nactv
    write(outfl_unit, 11) "xlevel    =", xlvl
    write(outfl_unit, 10) ""

    return
10  format(1x,A)
11  format(2x,A,i5)
  end subroutine write_gennml
  
  subroutine write_header(outfl_unit)
    !===========================================================================
    ! write_header
    ! ------------
    ! Purpose: Write header for jayci_exp.
    !---------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: outfl_unit

    write(unit = outfl_unit, fmt = "(A)") " ------------------------------------"
    write(unit = outfl_unit, fmt = "(A)") " jayci.x "
    write(unit = outfl_unit, fmt = "(A)") " ------------------------------------"

    return
  end subroutine write_header
end program jayci
